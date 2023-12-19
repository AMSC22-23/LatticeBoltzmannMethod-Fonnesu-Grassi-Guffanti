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
#   As far as usage is concerned, the python script is called by running the following command, for 
#   density and structure:
#   python translate.py path_to_image lattice|rho path_to_output_dir
#   
#   To produce the velocity field the command to call is the following
#   python translate.py path_to_image ux|uy path_to_output_dir [open]
#      
#
#

import sys
import os
import random
import numpy
import math
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

def produce_velocity(img: Image, open) -> tuple[int,int,int,list]:
    
    width, height = img.size
    pixels = list(img.getdata())

    non_zero = 0
    i = 0
    j = 0
    non_zeroes = []
    
    const = (width/1.8)
    if open:
        for idx, pixel in enumerate(pixels):
            j = idx % (width)
            i = int((idx - j) / width)
            if is_fluid(pixel):
                if i == 0 or j == 0 or i == height - 1 or j == height - 1:
                    non_zero += 1 
                    non_zeroes.append((random.uniform((-1+(j/const)),(-0.8+(j/const))), i, j))
    else:
        for idx, pixel in enumerate(pixels):
            j = idx % (width)
            i = int((idx - j) / width)
            if is_fluid(pixel):
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
        if len(args) == 5:
            if "open" == args[4]:
                width, height, non_zero, non_zeroes = produce_velocity(img,True)
            else:
                print(f"{args[4]} is not an accepted parameter")
        else:
            width, height, non_zero, non_zeroes = produce_velocity(img,False)
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

if __name__ == "__main__":
    execute_translate()