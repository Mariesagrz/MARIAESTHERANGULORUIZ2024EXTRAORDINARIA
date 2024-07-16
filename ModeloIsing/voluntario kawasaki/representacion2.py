import numpy as np
import matplotlib.pyplot as plt

# Leer los datos del fichero
data = np.loadtxt('magnetizacion.txt')

# Asumimos que el fichero tiene tres columnas:
# Primera columna: variable independiente (x)
# Segunda columna: primera variable dependiente (y1)
# Tercera columna: segunda variable dependiente (y2)
x = data[:, 0]
y1 = data[:, 1]
y2 = data[:, 2]

# Crear el gráfico con ambas variables dependientes
plt.figure()
plt.plot(x, y1, label='Magnetizacion Superior', color='blue', marker='o')
plt.plot(x, y2, label='Magnetizacion Inferior', color='red', marker='x')
plt.xlabel('Tiempo')
plt.ylabel('Magnetización')
plt.title('')
plt.legend()
plt.grid(True)

# Mostrar el gráfico
plt.show()
