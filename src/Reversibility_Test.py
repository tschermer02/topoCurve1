import numpy as np
import matplotlib.pyplot as plt
from SpectralFiltering import SpectralFiltering  # Import the spectral filtering class

# Generate a known surface using sine and cosine functions
def generate_surface(size=100, scale=10):
    x = np.linspace(-np.pi, np.pi, size)
    y = np.linspace(-np.pi, np.pi, size)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X) * np.cos(Y) * scale
    return X, Y, Z

# Function to apply spectral filtering
def apply_spectral_filtering(Z, dx, dy):
    """Applies spectral filtering to the surface."""
    spectral_filter = SpectralFiltering.__new__(SpectralFiltering)  # Instantiate without calling __init__
    
    # Manually assign necessary attributes
    spectral_filter.z_array = Z 
    spectral_filter.dx = [dx, dy]

    # Ensure detrending happens before FFT
    spectral_filter.detrend()

    # Apply filtering
    dx_filt, dy_filt, Z_filtered = spectral_filter.FFT([0, 200], 'lowpass', 0)
    
    return Z_filtered, dx_filt, dy_filt

# **Step 1: Generate the original surface**
X, Y, Z_original = generate_surface()
dx = dy = np.abs(X[0, 1] - X[0, 0])  # Calculate grid spacing

# **Step 2: Apply Spectral Filtering in a Loop**
num_iterations = 50  # Run the filtering 50 times

for i in range(1, num_iterations + 1):
    # Regenerate surface for each iteration
    X, Y, Z_new = generate_surface()
    
    # Apply Spectral Filtering
    Z_filtered, dx_filtered, dy_filtered = apply_spectral_filtering(Z_new, dx, dy)

    # Calculate error as absolute difference from the original surface
    error = np.abs(Z_original - Z_filtered)
    mean_error = np.mean(error)  # Compute mean error for summary

    # Display results every 10 iterations
    if i % 10 == 0:
        print(f"Iteration {i}: Mean Absolute Error = {mean_error:.4f}")

        fig, ax = plt.subplots(1, 3, figsize=(15, 5))
        
        ax[0].imshow(Z_new, cmap='viridis')
        ax[0].set_title(f'Iteration {i}: Regenerated Surface')
        
        ax[1].imshow(Z_filtered, cmap='viridis')
        ax[1].set_title(f'Iteration {i}: Spectrally Filtered Surface')

        ax[2].imshow(error, cmap='inferno')
        ax[2].set_title(f'Iteration {i}: Error (|Original - Filtered|)')

        plt.show()
