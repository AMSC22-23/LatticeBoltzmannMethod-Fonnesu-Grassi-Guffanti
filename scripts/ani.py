import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import os


images= []

for file in os.listdir("scripts/"):
    if file.endswith(".png"):
        images.append(f"scripts/{file}")

images= sorted(images)
images.sort(key=len)

def update(frame):
    plt.clf()
    plt.axis("off")
    plt.imshow(plt.imread(images[frame]))

fig = plt.figure()


# Crea l'animazione utilizzando FuncAnimation di Matplotlib
ani = FuncAnimation(fig, update, frames=len(images), interval=100)  # Intervallo in millisecondi tra le immagini

ani.save('scripts/animazione.mp4', writer='ffmpeg')
