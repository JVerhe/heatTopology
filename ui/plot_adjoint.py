import matplotlib.pyplot as plt
import numpy as np

# import matplotlib
# matplotlib.use('TkAgg')

# Penalty values
p = [1, 2, 3, 4, 5]

# Maximum deviations for each filtering type
max_dev_ft0 = [0.499576, 0.493083, 0.499551, 0.497682, 0.458566]
max_dev_ft1 = [0.472794, 0.49984, 0.492453, 0.372409, 0.222479]
max_dev_ft2 = [0.49601, 0.496399, 0.495673, 0.488637, 0.298378]

# Mean deviations for each filtering type
mean_dev_ft0 = [0.19392, 0.0617476, 0.0522735, 0.0288398, 0.00719593]
mean_dev_ft1 = [0.0141053, 0.0205385, 0.0123052, 0.00441697, 0.00166143]
mean_dev_ft2 = [0.193549, 0.0589912, 0.0402547, 0.0249582, 0.00720434]
tolerance = [0.2,0.2/2,0.2/3,0.2/4,0.2/5]

# Create subplots
plt.figure(figsize=(7, 5))

# Plot for maximum deviations
# ax1.plot(p, max_dev_ft0, marker='o', label='No Filtering (ft=0)')
# ax1.plot(p, max_dev_ft1, marker='o', label='Sensitivity Filtering (ft=1)')
# ax1.plot(p, max_dev_ft2, marker='o', label='Density Filtering (ft=2)')
# ax1.set_xlabel('Penalty (p)')
# ax1.set_ylabel('Maximum Deviation from 0 or 1')
# ax1.set_title('Maximum Deviation vs Penalty')
# ax1.legend()
# ax1.grid(True)

# Plot for mean deviations
plt.plot(p, mean_dev_ft0, marker='o', label='No Filtering')
plt.plot(p, mean_dev_ft1, marker='o', label='Sensitivity Filtering')
plt.plot(p, mean_dev_ft2, marker='o', label='Density Filtering')
plt.plot(p,tolerance,marker='o',label='Test Tolerance')
plt.xlabel('Penalty (p)')
plt.ylabel('Mean Deviation from 0 or 1')
plt.title('Mean Deviation vs Penalty')
plt.legend()
plt.grid(True)

plt.show()


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
