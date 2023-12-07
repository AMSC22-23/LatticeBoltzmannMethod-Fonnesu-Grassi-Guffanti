import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os
import sys


args = sys.argv
images= []

print(f"animating images in {args[1]}")

for file in os.listdir(f"{args[1]}"):
    if file.endswith(".png"):
        images.append(f"{args[1]}{file}")

images = sorted(images)
images.sort(key=len)

def update(frame):
    plt.clf()
    plt.axis("off")
    plt.imshow(plt.imread(images[frame]))

fig = plt.figure()


# Crea l'animazione utilizzando FuncAnimation di Matplotlib
ani = FuncAnimation(fig, update, frames=len(images), interval=50)  # Intervallo in millisecondi tra le immagini

ani.save(f'{args[2]}.mp4')
