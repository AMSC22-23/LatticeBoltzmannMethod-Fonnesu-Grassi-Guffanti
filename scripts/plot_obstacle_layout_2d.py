# plot_obstacle_layout_2d.py - A script to plot the obstacle layout in 2D in a scatter plot

import os
import sys
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

# Read the first argument as the data file
data_file = sys.argv[1]

# Read the file as a CSV
data = pd.read_csv(data_file)

# Extract X, Y from the data
X = data.iloc[:, 1]
Y = data.iloc[:, 0]

# Print the read points
print(X,Y)

# Create a new figure for the obstacle layout
fig, ax = plt.subplots(figsize=(7, 6))

# Plot the obstacle layout as a scatter plot
sc_obstacle = ax.scatter(X, Y, c='black', marker='s')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('Obstacle Layout')
plt.savefig('plot_obstacle_layout.png')
plt.close(fig) 
