import matplotlib.pyplot as plt 
import numpy as np 



log_h = np.array([-3.38, -3.69, -3.86923, -3.9956, -4.09321, -4.240549, -4.298853, -4.3961])
log_E = [-8.022, -9.26, -9.97, -10.484044, -10.877444, -11.4859, -11.7461, -12.3464]

coeffs = np.polyfit(log_h, log_E, 1)
fit_line = np.poly1d(coeffs)

plt.figure(figsize=(8, 6))
plt.scatter(log_h, log_E, label="Data", color='dodgerblue', edgecolors='black', s=80, alpha=0.7, marker='o')
plt.plot(log_h, fit_line(log_h), label=f"Linear fit: y={coeffs[0]:.2f}x + {coeffs[1]:.2f}", color='tomato', linewidth=2, linestyle='--')

plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.xlabel("log(h)", fontsize=14, fontweight='bold')
plt.ylabel("log(E)", fontsize=14, fontweight='bold')
plt.title("Linear Fit of log(E) vs log(h)", fontsize=16, fontweight='bold')
plt.legend(fontsize=12)

plt.tight_layout()
plt.show()



################ K_VARIABLE #########


log_h = np.array([-3.38, -3.59, -3.77, -3.89, -3.99, -4.17, -4.29, -4.39, -4.47, -4.54, -4.6, -4.65, -4.77])
log_E = [-5.32, -5.7, -6.04, -6.28, -6.47, -6.8, -7.0533, -7.24, -7.3955, -7.5256, -7.638, -7.7587, -7.9827]

coeffs = np.polyfit(log_h, log_E, 1)
fit_line = np.poly1d(coeffs)

plt.figure(figsize=(8, 6))
plt.scatter(log_h, log_E, label="Data", color='royalblue', edgecolors='black', s=80, alpha=0.7, marker='o')
plt.plot(log_h, fit_line(log_h), label=f"Linear fit: y={coeffs[0]:.2f}x + {coeffs[1]:.2f}", color='crimson', linewidth=2, linestyle='--')

plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.xlabel("log(h)", fontsize=14, fontweight='bold')
plt.ylabel("log(E)", fontsize=14, fontweight='bold')
plt.title("Linear Fit of log(E) vs log(h) with Variable k", fontsize=16, fontweight='bold')
plt.legend(fontsize=12)

plt.tight_layout()
plt.show()
