import numpy as np
import matplotlib.pyplot as plt

# Load data from text files
x_values = np.loadtxt("build/adjoint_eps_vec.txt") # or ../build/adjoint_eps_vec.txt
y_values_mean = np.loadtxt("build/adjoint_mean_error.txt")
y_values_mse = np.loadtxt("build/adjoint_MSE.txt")

# Create the plot
plt.figure(figsize=(8, 6))
plt.plot(x_values, y_values_mean, label="Mean Error", color="blue", marker="o")
plt.plot(x_values, y_values_mse, label="MSE", color="red", marker="s")

# Labels and title
plt.xlabel("Epsilon Values")
plt.ylabel("Error Values")
plt.title("Adjoint Sensitivity Analysis")
plt.legend()
plt.grid(True)

# Show the plot
plt.show()