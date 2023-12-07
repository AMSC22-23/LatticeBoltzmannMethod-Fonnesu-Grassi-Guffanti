#
#
#
#
#
#
#
#
#

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import sys

args = sys.argv

MINIMUM_VALUE = 0.0
MAXIMUM_VALUE = 1.0
INTERPOLATION = 'gaussian'

def load_directory():
    list = []
    print(f"    loading matrices in: {args[1]}")
    for file in os.listdir(f"{args[1]}"):
        if file.endswith("rho.txt"):
            list.append(f"{args[1]}{file}")
    list = sorted(list)
    list.sort(key=len)
    return list


def render_image(matrix,i): 
    if not os.path.exists(f"{args[3]}"):
        os.mkdir(f"{args[3]}")

    plt.imshow(matrix, cmap ='coolwarm', 
               vmin=MINIMUM_VALUE, 
               vmax=MAXIMUM_VALUE, 
               interpolation=INTERPOLATION
               ) 
    plt.title(f'{args[4]}')
    plt.savefig(f"{args[3]}image{i+1}.png")


def render_images():
    print(f"rendering: {args[2]}")
    file_paths = load_directory()
    print("     matrices loaded")
    for i, file in enumerate(file_paths):
        print(f"rendering {file}: {i}/{len(file_paths)}")
        matrix = np.loadtxt(file)
        render_image(matrix, i)


if __name__ == "__main__":
    if len(args) != 5:
        print("usage: python render.py path_to_input_dir rho|u path_to_output_dir title")
    else:
        render_images()

