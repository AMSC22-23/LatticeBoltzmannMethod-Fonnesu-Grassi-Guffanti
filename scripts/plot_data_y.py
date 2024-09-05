# Plot the selected velocity of a specified horizontal line

# How to use
# python plot_data_x.py <file_path_x> <file_path_y> <line>

import sys
import numpy as np
import matplotlib.pyplot as plt

def plot_selected_row(file_path_x, file_path_y, column_number):
    data_x = np.loadtxt(file_path_x)
    data_y = np.loadtxt(file_path_y)
    
    if column_number < 1 or column_number > data_x.shape[1]:
        print(f"Error: file X has only {data_x.shape[1]} columns.")
        sys.exit(1)
    
    if column_number < 1 or column_number > data_y.shape[1]:
        print(f"Error: file Y has only {data_y.shape[1]} columns.")
        sys.exit(1)
    
    selected_column_x = data_x[:, column_number - 1]
    x_data_x = np.arange(1, len(selected_column_x) + 1)
    y_data_x = selected_column_x

    selected_column_y = data_y[:, column_number - 1]
    x_data_y = np.arange(1, len(selected_column_y) + 1)
    y_data_y = selected_column_y

    plt.scatter(x_data_x, y_data_x, marker='o', label='$u_x$')
    plt.scatter(x_data_y, y_data_y, marker='o', label='$u_y$')

    plt.xticks(np.arange(10, max(x_data_x) + 1, 10))
    plt.title("velocities along Vertical Line through Geometric Center of Cavity")
    plt.xlabel('y', fontsize=20)
    plt.ylabel('Value', fontsize=20)
    plt.grid(True)

    plt.legend()

    plt.savefig("v-velocities.pdf")
    plt.show()
    print("Plot saved as 'v-velocities.pdf'")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Use: python plot_data.py <file_path_x> <file_path_y> <line>")
        sys.exit(1)

    file_path_x = sys.argv[1]
    file_path_y = sys.argv[2]
    row_number = int(sys.argv[3])

    plot_selected_row(file_path_x, file_path_y, row_number)