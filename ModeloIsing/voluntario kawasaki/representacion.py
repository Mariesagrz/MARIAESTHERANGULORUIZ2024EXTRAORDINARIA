import matplotlib.pyplot as plt

# Leer el fichero
datos = []
with open('tempertura.txt', 'r') as file:
    for line in file:
        columnas = line.strip().split('\t')
        datos.append((float(columnas[0]), float(columnas[1])))

# Separar las columnas en dos listas
x, y = zip(*datos)

# Crear la gráfica
plt.figure(figsize=(10, 6))
plt.plot(x, y, marker='o', linestyle='-', color='b', label='Datos')
plt.xlabel('Temperatura')
plt.ylabel('Magnetizacion Total N=16')
plt.title('')
plt.legend()
plt.grid(True)
plt.show()
