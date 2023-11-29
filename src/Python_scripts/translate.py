#comando per runnare
#python src/Python_scripts/translate.py immagine_input.jpg nome_output

import sys
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

length, height =img.size
nonzeros = 0

with open(f"resources/patterns/{args[2]}.txt", "w") as file:
    with open(f"resources/patterns/{args[2]}.mtx","w") as mtx:
        mtx.write("%%MatrixMarket matrix coordinate real general\n")
        mtx.write("\n")
        for x in range(height):
            for y in range(length):
                pixel=img.getpixel((y,x))
                #print(f"Pixel alla posizione ({x}, {y}): {pixel}")
                if(is_boundary(pixel)) :
                    file.write("B")
                    nonzeros=nonzeros+1
                    mtx.write(f"{x} {y} 1\n")
                    #file.write("\033[30ml\033[0m")
                elif(is_fluid(pixel)):
                    file.write("W")
                    #file.write("\033[37ml\033[0m")
                elif(is_solid(pixel)):
                    file.write("S")
                    nonzeros=nonzeros+1
                    mtx.write(f"{x} {y} 3\n")
                    #file.write("\033[34ml\033[0m")
                else: 
                    file.write("?")
            file.write("\n")        
            #file.write("                                                                                                                                                                           ")
            # Blu
            # Nero

with open(f"resources/patterns/{args[2]}.mtx","r") as mtx:  
    lines=mtx.readlines()
lines[1]=f"{length} {height} {nonzeros}\n"

with open(f"resources/patterns/{args[2]}.mtx","w") as mtx:
    mtx.writelines(lines)


img.close

