import numpy as np

# Given array
array = np.array([[4,3,3,4,4,3],
                  [2,1,1,2,2,1],
                  [2,1,1,2,2,1],
                  [4,3,3,4,4,3],
                  [4,3,3,4,4,3],
                  [2,1,1,2,2,1]])

# Define Tukey window parameters
alpha = 0.5  # Adjust this value for different window shapes

# Create Tukey window
M, N = array.shape
tukey_window = np.zeros((M, N))
for i in range(M):
    for j in range(N):
        tukey_window[i, j] = 0.5 * (1 + np.cos(np.pi * (-1 + 2 * (j + 1) / (N - 1)))) if ((j + 1) >= (alpha * (N - 1) / 2)) and ((j + 1) <= (N * (1 - alpha / 2))) else 0

# Element-wise multiplication
array_with_tukey = array * tukey_window

# Print result
print(array_with_tukey)
