import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os

def visualizza_griglia(matrice,i):
    plt.imshow(matrice, cmap='coolwarm') 
    plt.title('Lid-Driven Cavity')
    plt.gca().invert_yaxis()  
    plt.savefig(f'scripts/image{i+1}.png')

for i in range(50):
    matrice_input = np.random.rand(100, 100)
    visualizza_griglia(matrice_input,i)


