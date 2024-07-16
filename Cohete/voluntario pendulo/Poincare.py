import numpy as np
import matplotlib.pyplot as plt

# Leer los datos del fichero
data = np.loadtxt('Poincare3.txt', delimiter=',')

# Asumimos que el fichero tiene dos columnas:
# Primera columna: phi
# Segunda columna: psi
phi = data[:, 0]
psi = data[:, 1]

# Crear el mapa de Poincaré
plt.figure()
plt.scatter(phi, psi, s=1, color='blue')  # Usamos scatter para un mapa de Poincaré
plt.xlabel('phi')
plt.ylabel('phi_punto')
plt.title('')
plt.grid(True)

# Mostrar el gráfico
plt.show()
