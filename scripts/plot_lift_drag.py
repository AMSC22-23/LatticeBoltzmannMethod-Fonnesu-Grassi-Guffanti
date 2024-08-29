import os
import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

# Read the first argument as the directory
DIRECTORY = sys.argv[1]

# Read the content of the directory
files = os.listdir(DIRECTORY)
# Pop from the files list the ones that are not CSV files
for file in files:
    if not file.endswith('.csv'):
        files.remove(file)
# And the global results one
if 'total_lift_drag.csv' in files:
    files.remove('total_lift_drag.csv')

# 1. Start by creating a directory to contain the plots
if not os.path.exists('flow_analysis_plots'):
    os.makedirs('flow_analysis_plots')
else:
    # If the directory already exists, remove the files in it
    for file in os.listdir('flow_analysis_plots'):
        os.remove(os.path.join('flow_analysis_plots', file))

# 2 Read the total_lift_drag csv file and plot the values
total_data = pd.read_csv(os.path.join(DIRECTORY, 'total_lift_drag.csv'))
# Extract the Lift and Drag values
Iteration = total_data.iloc[:,0]
Lift = total_data.iloc[:,1]
Drag = total_data.iloc[:,2]
# Create a new figure and plot the Lift and Drag values
fig, ax = plt.subplots(figsize=(7, 6))
ax.plot(Iteration, Lift, label='Lift', color='blue')
ax.plot(Iteration, Drag, label='Drag', color='red')
ax.set_xlabel('Iteration')
ax.set_ylabel('Force')
ax.set_title('Total Lift and Drag per Iteration')
ax.legend()
plt.savefig(os.path.join('flow_analysis_plots', 'total_lift_drag.pdf'), format="pdf")
plt.close()

# Order the contributions files by iteration
files.sort(key=lambda x: int(x.split('_')[1].split('.')[0]))
print(files)
# 3. Loop over the files in the directory and for each file render a frame of the animation
for file in files:
    # Read the file as a CSV
    data = pd.read_csv(os.path.join(DIRECTORY, file))
    # Merge all that in a big dataframe
    if file == files[0]:
        all_data = data
    else:
        all_data = pd.concat([all_data, data])

# Extract the maximum and mimnum values of Lift and Drag
max_Lift = all_data.iloc[:,2].max()
min_Lift = all_data.iloc[:,2].min()
max_Drag = all_data.iloc[:,3].max()
min_Drag = all_data.iloc[:,3].min()

# Now, for each file in the directory crete a plot with the Lift and Drag values and an animation
for file in files:
    # Read the file as a CSV
    data = pd.read_csv(os.path.join(DIRECTORY, file))
    # Extract X, Y, Lift, and Drag from the data
    X = data.iloc[:, 1]
    Y = data.iloc[:, 0]
    Lift = data.iloc[:, 2]
    Drag = data.iloc[:, 3]
    # Create subplots for Lift and |Lift|
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    # Plot Lift with the colorbar
    sc_lift = axes[0].scatter(X, Y, c=Lift, cmap='viridis', marker='o', alpha=0.6, vmin=min_Lift, vmax=max_Lift)
    axes[0].set_xlabel('X')
    axes[0].set_ylabel('Y')
    axes[0].set_title('Lift')
    axes[0].invert_yaxis()
    fig.colorbar(sc_lift, ax=axes[0])
    # Plot |Lift| with the colorbar
    sc_abs_lift = axes[1].scatter(X, Y, c=np.abs(Lift), cmap='viridis', marker='o', vmin=0, vmax=max_Lift, s=np.abs(Lift)/np.max(np.abs(Lift))*100, alpha=0.6)
    axes[1].set_xlabel('X')
    axes[1].set_ylabel('Y')
    axes[1].set_title('|Lift|')
    axes[1].invert_yaxis()
    fig.colorbar(sc_abs_lift, ax=axes[1])
    # Save the Lift and |Lift| plot to a file
    plt.savefig(os.path.join('flow_analysis_plots', file.split('.')[0] + '_lift.png'), format="png")
    plt.close(fig)
    # Create subplots for Drag and |Drag|
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    # Plot Drag with the colorbar  
    sc_Drag = axes[0].scatter(X, Y, c=Drag, cmap='viridis', marker='o', alpha=0.6, vmin=min_Drag, vmax=max_Drag)
    axes[0].set_xlabel('X')
    axes[0].set_ylabel('Y')
    axes[0].set_title('Drag')
    axes[0].invert_yaxis()
    fig.colorbar(sc_Drag, ax=axes[0])
    # Plot |Drag| with the colorbar
    sc_abs_Drag = axes[1].scatter(X, Y, c=np.abs(Drag), cmap='viridis', marker='o', vmin=0, vmax=max_Drag, s=np.abs(Drag)/np.max(np.abs(Drag))*100, alpha=0.6)
    axes[1].set_xlabel('X')
    axes[1].set_ylabel('Y')
    axes[1].set_title('|Drag|')
    axes[1].invert_yaxis()
    fig.colorbar(sc_abs_Drag, ax=axes[1])
    # Save the Drag and |Drag| plot to a file
    plt.savefig(os.path.join('flow_analysis_plots', file.split('.')[0] + '_drag.png'), format="png")
    plt.close(fig)