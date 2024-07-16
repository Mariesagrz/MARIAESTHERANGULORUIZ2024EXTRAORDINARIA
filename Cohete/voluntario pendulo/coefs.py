import numpy as np
import matplotlib.pyplot as plt

# Leer los datos del fichero
data = np.loadtxt('lyapunov.txt', delimiter=',')
t = data[:, 0]
dist = data[:, 1]

# Calcular el exponente de Lyapunov
d0 = 1e-8
lyapunov = np.log(dist / d0) / t

# Graficar el exponente de Lyapunov
plt.figure()
plt.plot(t, lyapunov, label='Coeficiente de Lyapunov')
plt.xlabel('Tiempo')
plt.ylabel('Coeficiente de Lyapunov')
plt.title('Coeficiente de Lyapunov vs Tiempo')
plt.legend()
plt.grid(True)

# Mostrar el gráfico
plt.show()

# Calcular el coeficiente de Lyapunov promedio
lyapunov_avg = np.mean(lyapunov[int(len(lyapunov)/2):])  # Usar la segunda mitad de los datos para el promedio
print(f'Coeficiente de Lyapunov promedio: {lyapunov_avg}')
