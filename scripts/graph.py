import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

def visualizza_griglia(matrice,i):
    plt.imshow(matrice, cmap='coolwarm') 
    plt.title('Lid-Driven Cavity')
    plt.gca().invert_yaxis()  
    plt.savefig(f'scripts/image{i+1}.png')

list =[]
for file in os.listdir("resources/lattices/lid_driven_cavity/results"):
    if file.endswith("rho.txt"):
        list.append(f"resources/lattices/lid_driven_cavity/results/{file}")

list=sorted(list)
list.sort(key=len)

for i , file in enumerate(list):
    matrice_input = np.loadtxt(file, usecols= range(100))
    visualizza_griglia(matrice_input,i)


