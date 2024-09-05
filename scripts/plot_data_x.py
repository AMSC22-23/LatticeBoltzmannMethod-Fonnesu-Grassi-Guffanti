# Plot the selected velocity of a specified horizontal line

# How to use
# python plot_data_x.py <file_path_x> <file_path_y> <line>

import sys
import numpy as np
import matplotlib.pyplot as plt

def plot_selected_row(file_path_x, file_path_y, row_number):
    data_x = np.loadtxt(file_path_x)
    data_y = np.loadtxt(file_path_y)
     
    if row_number < 1 or row_number > data_x.shape[0]:
        print(f"Error: file has only {data_x.shape[0]} rows.")
        sys.exit(1)
    
    if row_number < 1 or row_number > data_y.shape[0]:
        print(f"Error: file has only {data_y.shape[0]} rows.")
        sys.exit(1)
    
    selected_row = data_x[row_number - 1]  
    x_data_x = np.arange(1, len(selected_row) + 1)
    y_data_x = selected_row

    selected_row = data_y[row_number - 1]
    x_data_y = np.arange(1, len(selected_row) + 1)
    y_data_y = selected_row

    plt.scatter(x_data_x, y_data_x, marker='o', label='$u_x$')
    plt.scatter(x_data_y, y_data_y, marker='o', label='$u_y$')

    plt.xticks(np.arange(10, max(x_data_x) + 1, 10))
    plt.title("velocities along Horizontal Line through Geometric Center of Cavity")
    plt.xlabel('x', fontsize=20)
    plt.ylabel('Value', fontsize=20)
    plt.grid(True)

    plt.legend()

    plt.savefig("h-velocities.pdf")
    plt.show()
    print("Plot saved as 'h-velocities.pdf'")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Use: python plot_data.py <file_path_x> <file_path_y> <line>")
        sys.exit(1)

    file_path_x = sys.argv[1]
    file_path_y = sys.argv[2]
    row_number = int(sys.argv[3])

    plot_selected_row(file_path_x, file_path_y, row_number)