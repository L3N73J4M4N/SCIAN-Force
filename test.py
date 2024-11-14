import numpy as np
import pickle
import matplotlib.pyplot as plt

with open('rusienkiewicz.pkl', 'rb') as archivo:
    rusinki = pickle.load(archivo)

# Datos de ejemplo (asegúrate de reemplazar con tus valores reales)
valores_calculados = rusinki['calc_curvatures'][::-1]
valores_esperados = rusinki['real_curvatures'][::-1]
radios = rusinki['teoric_radius'][::-1]

# Cálculo del error cuadrático medio (MSE)
mse = np.sqrt(np.mean((valores_esperados - valores_calculados) ** 2))

# Graficar
plt.plot(1 / valores_esperados, 1 / valores_esperados, color='#0072B2', label='Expected Radius')
plt.plot(1 / valores_esperados, 1 / valores_calculados, '--', color='#D55E00', label='Calculated Radius')
plt.legend(loc='upper left', fontsize=10)

# Títulos y etiquetas
plt.title(f'Calculated vs Expected Radius', fontsize=14)
plt.xlabel(r'Expected Radius [$\mu$m]', fontsize=12)
plt.ylabel(r'Calculated Radius [$\mu$m]', fontsize=12)

# Estilo y cuadrícula
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()

# Mostrar el gráfico
plt.show()

# Elipses
J = np.linspace(2, 10, 20)
mse_x = np.load('MSE_ellipsoides_x.npy')
mse_y = np.load('MSE_ellipsoides_y.npy')
mse_z = np.load('MSE_ellipsoides_z.npy')
plt.plot(J, mse_x, label='eje x')
plt.plot(J, mse_y, '--', label='eje y')
plt.plot(J, mse_z, label='eje z')
plt.legend()
plt.show()

# Elipses (R)
# Cargar los datos (asegúrate de reemplazar con tus datos reales)
H_x = np.load('H_ellipsoides_x.npy')
H_y = np.load('H_ellipsoides_y.npy')
H_z = np.load('H_ellipsoides_z.npy')
R_x = np.load('R_ellipsoides_x.npy')
R_y = np.load('R_ellipsoides_y.npy')
R_z = np.load('R_ellipsoides_z.npy')

# Graficar
plt.plot(1 / R_x, 1 / np.abs(H_x), label='Radius incrementing x-axis', color='#0072B2')
plt.plot(1 / R_y, 1 / np.abs(H_y), '--', label='Radius incrementing y-axis', color='#D55E00')
plt.plot(1 / R_z, 1 / np.abs(H_z), label='Radius incrementing z-axis', color='#009E73')
plt.plot(1 / R_x, 1 / R_x, label='Expected radius', linestyle=':', color='gray')

# Leyenda
plt.legend(loc='upper left', fontsize=10)

# Títulos y etiquetas
plt.title('Effect of Directional Radius Increases on Ellipsoids Curvatures')
plt.xlabel(r'Isotropic Radius from Volume [$\mu$m]')
plt.ylabel(r'Isotropic Radius from Curvatures [$\mu$m]')

# Estilo y cuadrícula
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()

# Mostrar el gráfico
plt.show()