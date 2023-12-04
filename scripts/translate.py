#comando per runnare
#python scripts/translate.py immagine_input.jpg nome_output sottocartella

import sys
import os
from PIL import Image

def is_boundary(pixel): #tendente al bianco
    if (pixel[0]<100 and pixel[1]<100 and pixel[2]<100):
        return True
    return False

def is_fluid(pixel): #tendente al nero
    if (pixel[0]>200 and pixel[1]>200 and pixel[2]>200):
        return True
    return False

def is_solid(pixel): #condizioni da rivedere per ora pensavo qualsiasi sfumatura di blu come solido per semplicità
    if (pixel[0]<100 and pixel[1]<100 and pixel[2]>100):
        return True
    return False

args=sys.argv

img = Image.open(f"resources/img/{args[1]}")
density_img= Image.open(f"resources/img/{args[4]}")

length, height =img.size
nonzeros = 0
nonzeros_rho = 0

if not os.path.exists("resources/patterns"):
    os.mkdir("resources/patterns")

with open(f"resources/patterns/{args[2]}.txt", "w") as file:
    with open(f"resources/lattices/{args[3]}/{args[2]}.mtx","w") as mtx: # 0 fluid 1 open boundary 2 boundary 3 solid
        with open(f"resources/lattices/{args[3]}/rho.mtx","w") as mtx_rho:
            mtx.write("%%MatrixMarket matrix coordinate real general\n")
            mtx.write("\n")

            mtx_rho.write("%%MatrixMarket matrix coordinate real general\n")
            mtx_rho.write("\n")
            for x in range(height):
                for y in range(length):
                    pixel=img.getpixel((y,x))
                    x_=x+1
                    y_=y+1
                    #print(f"Pixel alla posizione ({x}, {y}): {pixel}")
                    if(is_boundary(pixel)) :
                        file.write("B")
                        nonzeros=nonzeros+1
                        mtx.write(f"{x_} {y_} 2\n")
                        #file.write("\033[30ml\033[0m")
                    elif(is_fluid(pixel)):
                        if(x == 0 or y == 0 or x == (height-1) or y == (length-1)):
                            file.write("O")
                            mtx.write(f"{x_} {y_} 1\n")
                            nonzeros=nonzeros+1
                        else:
                            file.write("W")
                        
                        pixel_rho =density_img.getpixel((y,x))
                        rho = 1 - ((pixel_rho[1] + pixel_rho[2])/ 510)
                        print(f"{rho}\n")
                        if(rho != 0):
                            nonzeros_rho = nonzeros_rho +1
                            mtx_rho.write(f"{x_} {y_} {rho}\n")
                        #file.write("\033[37ml\033[0m")
                    elif(is_solid(pixel)):
                        file.write("S")
                        nonzeros=nonzeros+1
                        mtx.write(f"{x_} {y_} 3\n")
                        #file.write("\033[34ml\033[0m")
                    else: 
                        file.write("?")
                file.write("\n")        
                #file.write("                                                                                                                                                                           ")
                # Blu
                # Nero

with open(f"resources/lattices/{args[3]}/{args[2]}.mtx","r") as mtx:  
    lines=mtx.readlines()
lines[1]=f"{length} {height} {nonzeros}\n"

with open(f"resources/lattices/{args[3]}/{args[2]}.mtx","w") as mtx:
    mtx.writelines(lines)

with open(f"resources/lattices/{args[3]}/rho.mtx","r") as mtx_rho:  
    lines=mtx_rho.readlines()
lines[1]=f"{length} {height} {nonzeros_rho}\n"

with open(f"resources/lattices/{args[3]}/rho.mtx","w") as mtx_rho:
    mtx_rho.writelines(lines)

img.close()
density_img.close()
