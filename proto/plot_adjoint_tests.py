import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

# Load data from text files
x_values = np.loadtxt("build/adjoint_eps_vec.txt")
y_values_mean = np.loadtxt("build/adjoint_mean_error.txt")
y_values_mse = np.loadtxt("build/adjoint_MSE.txt")

# Create the log-log plot
plt.figure(figsize=(8, 6))
#plt.loglog(x_values, y_values_mean, label="Mean Error", color="blue", marker="o")
plt.loglog(x_values, y_values_mse, label="MSE", color="red", marker="s")

# Labels and title
plt.xlabel("Epsilon Values (log scale)")
plt.ylabel("Error Values (log scale)")
plt.title("Adjoint Sensitivity Analysis (Log-Log Plot)")
plt.legend()
plt.grid(True, which="both", ls="--")

# Show the plot
plt.show()
