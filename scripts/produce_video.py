# render.py
# this python scripts renders a video starting from input matrices, expecting an input directory path and wether to produce
# an animation regarding velocities or densities.
#
# The basic behavior of the script is the following:
#   1. load the file paths in the given input directory
#   2. for each file, translate the matrix (or matrices if using velocities) into an image.
#       2.1 if requested, save the image for direct qualitative analysis
#       2.2 build a frame of the animation with an image
#   3. save the final animation
#
# Usage (from root directory of the project):
#   python ./scripts/render.py path_to_input_dir (rho|u) [save path_to_output_dir]
#   
#   path_to_input_dir -> input directory. if the directory contains images then they are used to build the animation.
#   rho               -> produces an animation regarding the density
#   u                 -> produces animations with velocity components and 2-norm of the velocity
#   save              -> saves images in path_to_output_dir/images/
#
#   The program performs a first analysis of input data to determine the maximum and minimum value present from the result files.
#
#


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import seaborn as sns
import os
import sys
import re

def extract_number(filename):
    # Use a regular expression to extract the number from the filename
    match = re.search(r'\d+', filename)
    return int(match.group()) if match else float('inf')


class Renderer:
    def __init__(self) -> None:

        # directory related stuff
        self.arg_num        :int            = len(sys.argv)
        self.usage_message  :str            = """
USAGE: python ./scripts/render.py path_to_input_dir (rho|u) path_to_output_dir [save]
    path_to_input_dir -> input directory. if the directory contains images then they are used to build the animation. If directory contains 
                         both images and result txt files, images are preferred.
    rho               -> produces an animation regarding the density
    u                 -> produces animations with velocity components and 2-norm of the velocity
    save              -> saves images in path_to_output_dir/images/"""    
         
        if (self.arg_num < 3):
            print(self.usage_message)
            exit(1)
        
        
        self.save_imgs      :bool           = False

        if self.arg_num == 5 and sys.argv[3] == "save":
            self.save_imgs = True
            self.output_dir     :str            = sys.argv[4]
        elif self.arg_num == 5 and sys.argv[3] != "save":
            print(f"{sys.argv[3]} was not recognized")
            print("\n\n"+self.usage_message)
        

        self.input_dir      :str            = sys.argv[1]
        self.quantity       :str            = sys.argv[2]
        self.file_paths     :list           = []
        self.animation_name :str            = "animation.mp4"
        self.img_extension  :str            = ".png"
        self.load_images    :bool           = False
        self.images         :list           = []
        self.interval_ms    :int            = 50
        self.animator       :FuncAnimation  = None
        self.quantity_rho   :str            = "rho"
        self.quantity_u     :str            = "u"
        self.dimensions     :int            = 2
        self.dpi            :int            = 200

        self.vmin           :float          = 0.0
        self.vmax           :float          = 0.2
        self.interpolation  :str            = "spline16"
        self.frame_index    :int            = 0
        self.frame_number   :int            = 0

    def check_input_data(self) -> bool:
        """
        Checks if the directory contains png images.
        """
        self.all_files = os.listdir(self.input_dir)
        if any(file.endswith(self.img_extension) for file in self.all_files):
            print("Images detected.")
            return True
        else:
            print("Results file detected.")
            return False
    
    def load_input_images(self) -> list:
        """
        Loads all the image files from the input directory
        """
        print(f"Loading images from {self.input_dir}")
        
        images = [file for file in self.all_files if file.endswith(".png")]
        images = sorted(images, key=extract_number)
        return images

    def update_frame(self, frame) -> None:
        plt.clf()
        plt.axis("off")
        plt.imshow(plt.imread(self.images[frame]))

    def animate_input_images(self) -> None:
        """
        Produces an animation starting from images loaded from disk
        """
        print("Building animation")
        figure = plt.figure()
        self.animator = FuncAnimation(fig=figure, func=self.update_frame, frames=len(self.images), interval=self.interval_ms)
        self.animator.save(self.animation_name)
        print("Finished")
    
    def update__rho_frame(self, frame):
        
        print(f"Frame {self.frame_index} out of {self.frame_number}")
        if self.save_imgs:
            plt.savefig(f"{self.output_dir}/image{self.frame_index}.png")
        self.frame_index += 1
        self.matrix = np.loadtxt(f"{self.input_dir}output-{self.frame_index}{self.quantity_rho}.txt")
        self.im.set_array(self.matrix)
        return self.im

    def render_and_animate_rho(self):

        figure, ax = plt.subplots(1,1)
        
        self.matrix = np.loadtxt(f"{self.input_dir}output-{self.frame_index}{self.quantity_rho}.txt")
        self.im = ax.imshow(self.matrix, 
                            cmap='coolwarm',
                            vmin=self.vmin, 
                            vmax=self.vmax, 
                            interpolation=self.interpolation)
        
        self.animator = FuncAnimation(fig=figure, func=self.update__rho_frame, frames=self.frame_number, interval=self.interval_ms)
        self.animator.save(self.animation_name)
        print("Finished")


    def update_u_frame_2D(self, frame):
        print(f"Frame {self.frame_index} out of {self.frame_number}")

        if self.save_imgs:
            plt.savefig(f"{self.output_dir}/image{self.frame_index}.png")
        if self.frame_index < self.frame_number:
            self.frame_index += 1
        
        self.matrix_x = np.loadtxt(f"{self.input_dir}output-{self.frame_index}{self.quantity_u}_x.txt")
        self.matrix_y = np.loadtxt(f"{self.input_dir}output-{self.frame_index}{self.quantity_u}_y.txt")
        self.matrix_norm = np.sqrt(np.multiply(self.matrix_x, self.matrix_x) + np.multiply(self.matrix_y, self.matrix_y))

        self.im_x.set_array(self.matrix_x)        
        self.im_y.set_array(self.matrix_y)
        self.im_n.set_array(self.matrix_norm)
        return self.im_x, self.im_y, self.im_n

    def render_and_animate_u_2D(self):
        self.figure, (self.ax_x, self.ax_y, self.ax_n) = plt.subplots(1, 3)
        self.matrix_x = np.loadtxt(f"{self.input_dir}output-{self.frame_index}{self.quantity_u}_x.txt")
        self.matrix_y = np.loadtxt(f"{self.input_dir}output-{self.frame_index}{self.quantity_u}_y.txt")
        self.matrix_norm = np.sqrt(np.multiply(self.matrix_x, self.matrix_x) + np.multiply(self.matrix_y, self.matrix_y))


        self.ax_x.set_title("x velocity", fontsize=10)
        self.ax_y.set_title("y velocity", fontsize=10)
        self.ax_n.set_title("velocity norm", fontsize=10)        

        self.im_x = self.ax_x.imshow(self.matrix_x,
            cmap='RdBu_r',
            vmin=-self.vmax,
            vmax=self.vmax,
            interpolation=self.interpolation)  
        
        self.c_x = plt.colorbar(self.im_x, ax=self.ax_x, shrink=0.3)
        self.c_x.ax.tick_params(labelsize=8)

        self.im_y = self.ax_y.imshow(self.matrix_y,
            cmap='RdBu_r',
            vmin=-self.vmax,
            vmax=self.vmax,
            interpolation=self.interpolation)
        
        self.c_y = plt.colorbar(self.im_y, ax=self.ax_y, shrink=0.3)
        self.c_y.ax.tick_params(labelsize=8)

        self.im_n = self.ax_n.imshow(self.matrix_norm,
            cmap='RdBu_r',
            vmin=self.vmin,
            vmax=self.vmax,
            interpolation=self.interpolation)
        
        self.c_n = plt.colorbar(self.im_n, ax=self.ax_n, shrink=0.3)
        self.c_n.ax.tick_params(labelsize=8)
        plt.tight_layout()

        self.animator = FuncAnimation(fig=self.figure, func=self.update_u_frame_2D, frames=self.frame_number, interval=self.interval_ms)
        self.animator.save(self.animation_name, dpi=self.dpi, fps=15)
        print("Finished")

    def update_u_frame_3D(self, frame):
        pass

    def render_and_animate_u_3D():
        print("ANIMATION OF 3D VELOCITIES NOT YET IMPLEMENTED.")
        pass

    def run(self) -> None:
        """
        The main function: controls the flow and based on input data manages its loading and rendering
        """
        self.load_images = self.check_input_data()
        if self.load_images:
            # images are already built
            self.images = self.load_input_images()
            self.animate_input_images()
        else:
            if self.save_imgs and not os.path.exists(self.output_dir):
                os.mkdir(self.output_dir)

            print(f"Loading data describing {self.quantity}")
            
            # based on whether the files describe a three dimensional velocity field or a two dimensional velocity field
            # the program must know which files to load. 
            if self.quantity == self.quantity_u and any(file.endswith("u_z.txt") for file in self.all_files):
                self.dimensions = 3
                print("Velocity field is three dimensional")
            elif self.quantity == self.quantity_u:
                self.dimensions = 2
                print("Velocity field is two dimensional")
            self.frame_number = len([file.endswith(".txt") for file in self.all_files])//3 - 1
            print(f"Animation will have {self.frame_number + 1} frames")
            # now onto the rendering and animating of inputs
            if self.quantity == self.quantity_rho:
                self.render_and_animate_rho()
            elif self.dimensions == 2:
                self.render_and_animate_u_2D()
            elif self.dimensions == 3:
                self.render_and_animate_u_3D()

if __name__ == "__main__":
    renderer = Renderer()
    renderer.run()
    print("Done!")