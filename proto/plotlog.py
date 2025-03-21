import matplotlib.pyplot as plt 
import numpy as np 


#log_h =np.array([-7.78322,-8.49699,-9.20029,-9.598475,-10.305613,-10.12262,-9.764225485202621])
log_h =np.array([-3.38,-3.69,-3.86923,-3.9956,-4.09321,-4.240549,-4.298853,-4.3961])
slope_2 = log_h*2

#log_E =[-17.7752,-20.631928,-23.447089,-26.354055,-28.72415,-27.735498,-25.75414831693885]
log_E =[-8.022,-9.26,-9.97,-10.484044,-10.877444,-11.4859,-11.7461,-12.3464]

coeffs = np.polyfit(log_h, log_E, 1)
fit_line = np.poly1d(coeffs)


plt.scatter(log_h, log_E, label="Données", color='blue')
plt.plot(log_h, fit_line(log_h), label=f"Fit linéaire: y={coeffs[0]:.2f}x + {coeffs[1]:.2f}", color='red')
plt.xlabel("log(h)")
plt.ylabel("log(E)")
plt.legend()
plt.title("Ajustement linéaire de log(E) en fonction de log(h)")
plt.show()