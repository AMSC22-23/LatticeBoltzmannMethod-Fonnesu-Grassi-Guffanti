from matplotlib import pyplot as plt
import sys
import numpy as np


vec = np.loadtxt(f"{sys.argv[1]}")
print(vec)
x_data, y_data = vec[:, 0], vec[:, 1]


plt.scatter(x_data, y_data, marker='o')
plt.xticks(np.arange(min(x_data), max(x_data)+2, 2))
plt.title("STRONG SCALABILITY ANALYSIS")
plt.xlabel('Threads', fontsize=20)
plt.ylabel('Time (s)', fontsize=20)
plt.legend()
plt.grid(True)

plt.savefig("strong_scalability_plot.png")
print("Plot saved as 'strong_scalability_plot.png'")
