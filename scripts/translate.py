# translate.py
#   python script used to convert input images (.png extension) that either describe
#       1. a computational domain that represents the structure of the lattice
#       2. a density field over that domain
#   into a .mtx file with values based on colors.
#   
#   More specifically the following scheme is adopted
#   + For case number 1 
#       - black => boundary
#       - white => fluid
#   + For case number 2
#       - white => zero density
#       - tonalities of red => increasing density until 1
#   
#   As far as usage is concerned, the python script is called by running the following command
#   python translate.py path_to_image [lattice|rho] path_to_output_dir
#

import sys
import os
import random
from PIL import Image

OPEN_BOUNDARY_NODE = 1
BOUNDARY_NODE = 2

def is_boundary(pixel: tuple[int, int, int]) -> bool: 
    """
    A pixel is considered a boundary if its color approaches black
    """
    return pixel[0] < 100 and pixel[1] < 100 and pixel[2] < 100

def is_fluid(pixel: tuple[int, int, int]) -> bool: 
    """
    A pixel is considered a fluid node if its color approaches white
    """
    return pixel[0] > 200 and pixel[1] > 200 and pixel[2] > 200

# def is_solid(pixel):
#     if (pixel[0]<100 and pixel[1]<100 and pixel[2]>100):
#         return True
#     return False


def produce_lattice(img: Image) -> tuple[int, int, int, list]:
    """
    Reads the image passed as input and returns the width, the height, the number of
    non zero pixels as given by our encoding and the list of said pixels
    """
    width, height = img.size
    pixels = list(img.getdata())

    non_zero = 0
    i = 0
    j = 0
    non_zeroes = []
    # number of boundaries, either solid or open.
    # the linearized index is idx = (i * width) + j
    # which means that j = idx % width and i = (idx - j)/width
    for idx, pixel in enumerate(pixels):
        j = idx % (width)
        i = int((idx - j) / width)
        if is_fluid(pixel):
            if i == 0 or j == 0 or i == height - 1 or j == height - 1:
                non_zero += 1 
                non_zeroes.append((OPEN_BOUNDARY_NODE, i, j))
        elif is_boundary(pixel):
                non_zero += 1
                non_zeroes.append((BOUNDARY_NODE, i, j))

    non_zero = len(non_zeroes)
    return (width, height, non_zero, non_zeroes)

def produce_rho(img: Image) -> tuple[int, int, int, list]:
    """
    Reads the image passed as input and returns the width, the height, the number of
    non zero pixels as given by our encoding and the list of said pixels
    """
    width, height = img.size
    pixels = list(img.getdata())

    non_zero = 0
    i = 0
    j = 0
    non_zeroes = []
    # number of boundaries, either solid or open.
    # the linearized index is idx = (i * width) + j
    # which means that j = idx % width and i = (idx - j)/width
    for idx, pixel in enumerate(pixels):
        j = idx % (width)
        i = int((idx - j) / width)
        
        if not is_boundary(pixel):
            rho=((1 - ((pixel[1] + pixel[2])/ 510)))
            if rho != 0:
                non_zero += 1 
                non_zeroes.append((rho, i, j))
                    
            

    non_zero = len(non_zeroes)
    return (width, height, non_zero, non_zeroes)

def produce_velocity(img: Image) -> tuple[int,int,int,list]:
    width, height = img.size
    pixels = list(img.getdata())

    non_zero = 0
    i = 0
    j = 0
    non_zeroes = []
    
    const = (width/1.8)
    print(f"{const}")
    for idx, pixel in enumerate(pixels):
        j = idx % (width)
        i = int((idx - j) / width)
        if is_fluid(pixel):
            if i == 0 or j == 0 or i == height - 1 or j == height - 1:
                non_zero += 1 
                non_zeroes.append((random.uniform((-1+(j/const)),(-0.8+(j/const))), i, j))
    non_zero = len(non_zeroes)
    return (width, height, non_zero, non_zeroes)

   



def execute_translate():
    """
    Performs the translation by reading image data
    """
    args = sys.argv
    img = Image.open(f"{args[1]}")

    width = 0
    height = 0
    non_zero = 0
    non_zeroes = []


    if "lattice" == args[2]:
        width, height, non_zero, non_zeroes = produce_lattice(img)
    elif "rho" == args[2]:
        width, height, non_zero, non_zeroes = produce_rho(img)
    elif "ux" == args[2] or "uy" == args[2]:
        width, height, non_zero, non_zeroes = produce_velocity(img)
    else:
        print("Command not recognized. Aborting.")
    
    # make sure that the output directory path exists
    if not os.path.exists(args[3]):
        os.mkdir(args[3])

    with open(f"{args[3]}" + f"{args[2]}.mtx", "w") as file:
        file.write("%%MatrixMarket matrix coordinate real general\n")
        file.write(f"{height} {width} {non_zero}\n")
        [file.write(f"{nz[1] + 1} {nz[2] + 1} {nz[0]}\n") for nz in non_zeroes]
    file.close()
    #             mtx.write("\n")

    #             mtx_rho.write("%%MatrixMarket matrix coordinate real general\n")
    #             mtx_rho.write("\n")
    #             for x in range(height):
    #                 for y in range(length):
    #                     pixel=img.getpixel((y,x))
    #                     x_=x+1
    #                     y_=y+1
    #                     #print(f"Pixel alla posizione ({x}, {y}): {pixel}")
    #                     if(is_boundary(pixel)) :
    #                         file.write("B")
    #                         nonzeros=nonzeros+1
    #                         mtx.write(f"{x_} {y_} 2\n")
    #                         #file.write("\033[30ml\033[0m")
    #                     elif(is_fluid(pixel)):
    #                         if(x == 0 or y == 0 or x == (height-1) or y == (length-1)):
    #                             file.write("O")
    #                             mtx.write(f"{x_} {y_} 1\n")
    #                             nonzeros=nonzeros+1
    #                         else:
    #                             file.write("W")
                            
    #                         pixel_rho =density_img.getpixel((y,x))
    #                         rho = 1 - ((pixel_rho[1] + pixel_rho[2])/ 510)
    #                         print(f"{rho}\n")
    #                         if(rho != 0):
    #                             nonzeros_rho = nonzeros_rho +1
    #                             mtx_rho.write(f"{x_} {y_} {rho}\n")
    #                         #file.write("\033[37ml\033[0m")
    #                     elif(is_solid(pixel)):
    #                         file.write("S")
    #                         nonzeros=nonzeros+1
    #                         mtx.write(f"{x_} {y_} 3\n")
    #                         #file.write("\033[34ml\033[0m")
    #                     else: 
    #                         file.write("?")
    #                 file.write("\n")        
    #                 #file.write("                                                                                                                                                                           ")
    #                 # Blu
    #                 # Nero

    # with open(f"resources/lattices/{args[3]}/{args[2]}.mtx","r") as mtx:  
    #     lines=mtx.readlines()
    # lines[1]=f"{length} {height} {nonzeros}\n"

    # with open(f"resources/lattices/{args[3]}/{args[2]}.mtx","w") as mtx:
    #     mtx.writelines(lines)

    # with open(f"resources/lattices/{args[3]}/rho.mtx","r") as mtx_rho:  
    #     lines=mtx_rho.readlines()
    # lines[1]=f"{length} {height} {nonzeros_rho}\n"

    # with open(f"resources/lattices/{args[3]}/rho.mtx","w") as mtx_rho:
    #     mtx_rho.writelines(lines)

    # img.close()
    # density_img.close()

if __name__ == "__main__":
    execute_translate()