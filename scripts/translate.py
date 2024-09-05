# translate.py
#   python script used to convert input images (.png extension) that either describe
#       1. a computational domain that represents the structure of the lattice
#       2. a density field over that domain
#   into a .txt file with values based on colors.
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
#   python translate.py path_to_image lattice|rho path_to_output_dir outputfilename
#   
#   To produce the velocity field the command to call is the following
#   python translate.py path_to_image ux|uy path_to_output_dir [open]
#      
#   1 S 0 Fluid 2 Boundary 3 Inlet 4 Outlet 5 Obstacle
#

import sys
import os
import random
import numpy
import math
from PIL import Image

SOLID = 1 #GREY I THINK
BOUNDARY = 2 #BLACK
INLET = 3 #BLUE
OUTLET = 4 #RED
OBSTACLE = 5 #GREEN

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
def is_inlet(pixel: tuple[int, int, int]) -> bool:
    """
    A pixel is considered a inlet node if its color approaches blue
    """
    return pixel[0] < 100 and pixel[1] < 100 and pixel[2] > 200
def is_outlet(pixel: tuple[int, int, int]) -> bool:
    """
    A pixel is considered a outlet node if its color approaches red
    """
    return pixel[0] > 200 and pixel[1] < 100 and pixel[2] < 100
def is_obstacle(pixel: tuple[int, int, int]) -> bool:
    """
    A pixel is considered a outlet node if its color approaches green
    """
    return pixel[0] < 100 and pixel[1] > 200 and pixel[2] < 100

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
    solids = 0
    boundaries = 0
    inlets = 0
    outlets = 0
    obstacles = 0
    i = 0
    j = 0
    non_zeroes = []
    
    # number of boundaries, either solid or open.
    # the linearized index is idx = (i * width) + j
    # which means that j = idx % width and i = (idx - j)/width
    for idx, pixel in enumerate(pixels):
        j = idx % (width)
        i = int((idx - j) / width)

    
        if is_boundary(pixel):
                boundaries += 1
                non_zeroes.append((BOUNDARY, i, j))
        elif is_inlet(pixel):
                inlets +=1
                non_zeroes.append((INLET, i, j))
        elif is_outlet(pixel):
                outlets +=1
                non_zeroes.append((OUTLET, i, j))
        elif is_obstacle(pixel):
                obstacles +=1
                non_zeroes.append((OBSTACLE, i, j))
        elif not(is_fluid(pixel)):
                solids += 1 
                non_zeroes.append((SOLID, i, j))
            

    non_zero = len(non_zeroes)
    fluids = (width * height) - non_zero
    return (width, height, fluids, solids, boundaries, inlets, outlets, obstacles, non_zeroes)

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
        width, height, fluids, solids, boundaries, inlets, outlets, obstacles, non_zeroes = produce_lattice(img)
    elif "rho" == args[2]:
        width, height, non_zero, non_zeroes = produce_rho(img)
    elif "ux" == args[2] or "uy" == args[2]:
        if len(args) == 6:
            if "open" == args[5]:
                width, height, non_zero, non_zeroes = produce_velocity(img,True)
            else:
                print(f"{args[5]} is not an accepted parameter")
        else:
            width, height, non_zero, non_zeroes = produce_velocity(img,False)
    else:
        print("Command not recognized. Aborting.")
    
    # make sure that the output directory path exists
    if not os.path.exists(args[3]):
        os.mkdir(args[3])

    with open(f"{args[3]}" + f"{args[4]}-{args[2]}.txt", "w") as file:
        file.write("2\n")
        file.write(f"{height} {width}\n{fluids} {solids} {boundaries} {inlets} {outlets} {obstacles}\n")
        [file.write(f"{nz[1]} {nz[2]} {nz[0]}\n") for nz in non_zeroes]
    file.close()

if __name__ == "__main__":
    execute_translate()