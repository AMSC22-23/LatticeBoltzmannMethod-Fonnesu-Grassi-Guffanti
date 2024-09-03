import numpy as np
from PIL import Image
import sys

# Load the image
image_path = sys.argv[1]
output_path = sys.argv[2]
image = Image.open(image_path)

# Convert the image to grayscale
gray_image = image.convert('L')

# Convert the grayscale image to a numpy array
image_array = np.array(gray_image)
# Extract the coordinates of black pixels
black_pixel_coords = np.column_stack(np.where(image_array == 0))

# Get the number of black pixels
num_black_pixels = black_pixel_coords.shape[0]

# Prepare the header with the number of black pixels
header = f"{num_black_pixels}\n"

# Convert the coordinates to the desired format (each coordinate pair on a new line)
formatted_coords = '\n'.join([f"{x} {y}" for x, y in black_pixel_coords])

# Save the header and coordinates to a text file
with open(output_path, 'w') as f:
    f.write(header)
    f.write(formatted_coords)

# Print the first few coordinates as a sample
print(header)
print(formatted_coords[:100])  # Print the first few coordinates (first 100 characters as a preview)
