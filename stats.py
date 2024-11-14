import numpy as np
from forces import Forces3D
import matplotlib.pyplot as plt
import utility_forces as uf

args_x = []
args_y, args_z = args_x.copy(), args_x.copy()
for i in range(1, 21):
    fx = Forces3D(r"C:\Users\matia\Desktop\ellipsoids\x_axis\IDL_Obj " + f"({i}).off",
                  interval='0.163x0.163x0.5')
    fy = Forces3D(r"C:\Users\matia\Desktop\ellipsoids\y_axis\IDL_Obj " + f"({i}).off",
                  interval='0.163x0.163x0.5')
    fz = Forces3D(r"C:\Users\matia\Desktop\ellipsoids\z_axis\IDL_Obj " + f"({i}).off",
                  interval='0.163x0.163x0.5')
    vx, vy, vz = fx.mesh[0], fy.mesh[0], fz.mesh[0]
    args_x.append(np.max(vx[:, 0]) - np.min(vx[:, 0]))
    args_y.append(np.max(vy[:, 1]) - np.min(vy[:, 1]))
    args_z.append(np.max(vz[:, 2]) - np.min(vz[:, 2]))

args_x = (np.argsort(args_x) + np.ones(20)).astype(int).tolist()
args_y = (np.argsort(args_y) + np.ones(20)).astype(int).tolist()
args_z = (np.argsort(args_z) + np.ones(20)).astype(int).tolist()
mse_x, mse_y, mse_z = np.zeros(20), np.zeros(20), np.zeros(20)
J = np.linspace(2, 10, 20)

print('eje x')
k = 1
for i in args_x:
    fx = Forces3D(r"C:\Users\matia\Desktop\ellipsoids\x_axis\IDL_Obj " + f"({i}).off",
                  interval='0.163x0.163x0.5')
    curves = np.abs(fx.rusinkiewicz_curvature())
    verts = fx.mesh[0]
    n = len(verts)
    a, b, c = J[k - 1], 2, 2
    error = np.zeros(n)
    for j in range(0, n):
        x, y, z = verts[j]
        ht = uf.theoretic_curvature_ellipsoid(x, y, z, a, b, c)
        error[j] = ht - curves[j]
    mse_x[k - 1] = 1 / n * np.sum(error ** 2)
    print(k, '/', 20)
    k += 1

print('\neje y')
k = 1
for i in args_y:
    fy = Forces3D(r"C:\Users\matia\Desktop\ellipsoids\y_axis\IDL_Obj " + f"({i}).off",
                  interval='0.163x0.163x0.5')
    curves = np.abs(fy.rusinkiewicz_curvature())
    verts = fy.mesh[0]
    n = len(verts)
    a, b, c = 2, J[k - 1], 2
    error = np.zeros(n)
    for j in range(0, n):
        x, y, z = verts[j]
        ht = uf.theoretic_curvature_ellipsoid(x, y, z, a, b, c)
        error[j] = ht - curves[j]
    mse_y[k - 1] = 1 / n * np.sum(error ** 2)
    print(k, '/', 20)
    k += 1

print('\neje z')
k = 1
for i in args_z:
    fz = Forces3D(r"C:\Users\matia\Desktop\ellipsoids\z_axis\IDL_Obj " + f"({i}).off",
                  interval='0.163x0.163x0.5')
    curves = np.abs(fz.rusinkiewicz_curvature())
    verts = fz.mesh[0]
    n = len(verts)
    a, b, c = 2, 2, J[k - 1]
    error = np.zeros(n)
    for j in range(0, n):
        x, y, z = verts[j]
        ht = uf.theoretic_curvature_ellipsoid(x, y, z, a, b, c)
        error[j] = ht - curves[j]
    mse_z[k - 1] = 1 / n * np.sum(error ** 2)
    print(k, '/', 20)
    k += 1
np.save('MSE_ellipsoides_x.npy', mse_x)
np.save('MSE_ellipsoides_y.npy', mse_y)
np.save('MSE_ellipsoides_z.npy', mse_z)
plt.plot(J, mse_x, label='eje x')
plt.plot(J, mse_y, label='eje y')
plt.plot(J, mse_z, label='eje z')
plt.legend()
plt.show()
