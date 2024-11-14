import numpy as np
from tifffile import imwrite
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize, to_hex


# módulo con funciones de utilidad para poder definir la clase de fuerzas
####################################################################################
# STATS
def stats(obj):
    vertices, faces = obj.mesh
    print(len(vertices))
    print(len(faces))
    print(obj.get_volume())
    print(obj.get_radius())
    curves = obj.rusinkiewicz_curvature()
    print(np.mean(curves))
    print(np.max(curves))
    print(np.min(curves))
    print(np.median(curves))
    print(np.std(curves))
    curves = obj.laplacian_curvature()
    print(np.mean(curves))
    print(np.max(curves))
    print(np.min(curves))
    print(np.median(curves))
    print(np.std(curves))


def histogram(obj, radius, title, value=1):
    # laplaciano
    lap = obj.laplacian_curvature()
    rad = obj.get_radius() * value
    plt.hist(lap, label='Laplacian', bins='auto')
    plt.plot(np.ones(50) * 1 / rad, np.linspace(0, np.max(np.histogram(lap, bins='auto')[0])),
             label=r'$\frac{1}{r}=$'+str(round(1/rad, 5)), color='red')
    plt.title(f'$R={radius}$ y $r={round(rad, 5)}$ ' + obj.method + ' ' + title)
    plt.legend()
    plt.show()
    # rusinkiewicz
    rus = obj.rusinkiewicz_curvature()
    plt.hist(rus, color='green', label='Rusinkiewicz', bins='auto')
    plt.plot(1 / rad * np.ones(50), np.linspace(0, np.max(np.histogram(rus, bins='auto')[0])),
             color='red', label=r'$\frac{1}{r}=$' + str(round(1 / rad, 5)))
    plt.title(f'$R={radius}$ y $r={round(rad, 5)}$ ' + obj.method + ' ' + title)
    plt.legend()
    plt.show()


def histogram_unit(obj, radius, title, value=1, value2=None, original=True):
    value2 = value if value2 is None else value2
    rad = obj.get_radius() * value
    curves = obj.meanCurvature
    plt.hist(curves, color='green', label=obj.method, bins='auto')
    plt.plot(1 / radius * value2 * np.ones(50), np.linspace(0, np.max(np.histogram(obj.meanCurvature, bins='auto')[0])),
             color='black', label='original =' + str(1 / radius * abs(value2)))
    plt.plot(1 / rad * np.ones(50), np.linspace(0, np.max(np.histogram(obj.meanCurvature, bins='auto')[0])),
             color='red', label=r'$\frac{1}{r}=$' + str(round(1 / rad, 5)))
    plt.title(f'$R={radius}$ y $r={round(rad, 5)}$ ' + obj.method + ' ' + title)
    plt.legend()
    plt.show()


#######################################################################################
############################################################################
# convierte un string con números (sin espacios al final y al inicio),
# separados por espacios, en números.
# ejemplo: str2num('1 2 3') -> [1.0, 2.0, 3.0]
def str2num(string, type_num='float'):
    string_list = []
    num_list = []
    string += ' '
    k = 0
    for i in range(0, len(string) - 1):
        if string[i + 1] == ' ' and string[i] != ' ':
            string_list.append(string[k: i + 1])
        if string[i] == ' ' and string[i + 1] != ' ':
            k = i + 1
    for s in string_list:
        if type_num == 'float':
            num_list.append(float(s))
        if type_num == 'int':
            num_list.append(int(s))
    return num_list


# en caso de tener un vector con elementos no deseados, los elimina.
# ejemplo: anti_nan([1, 2, np.nan, 3]) -> [1, 2, 3]
def anti_nan(vector):
    new_vector = []
    for p in vector:
        if p == np.nan or p == np.inf or p == -np.inf:
            continue
        new_vector.append(p)
    return new_vector


# cotangente, la inversa de la tangente
def cot(angle):
    return 1 / np.tan(angle)


############################################################################
# retorna la suma de las áreas baricéntricas de los triángulos
# formados por el punto i. Obtienes los lados de cada
# triángulo al cual pertenece i y se aplica la fórmula de Heron
def cell_area(i, vertices, triangles, three=True):
    area = 0
    a = vertices[i]
    tri = triangles[i]
    for j in range(0, len(tri)):
        b, c = tri[j]
        b, c = vertices[b], vertices[c]
        ab, bc, ca = np.linalg.norm(a - b), np.linalg.norm(b - c), np.linalg.norm(c - a)
        s = (ab + bc + ca) / 2
        area += (s * (s - ab) * (s - bc) * (s - ca)) ** (1 / 2)  # fórmula de Heron
    return (1 / 3) * area if three else area


def mixed_area(i, vertices, faces, triangles, angles):
    area = 0
    tri = triangles[i]
    for t in range(0, len(tri)):
        face = faces[tri[t]]
        k = np.where(face == i)[0]
        angles_tri = angles[tri[t]]
        if any(np.abs(angles_tri)) > np.pi / 2:  # no es obtuso
            p = vertices[i]
            k1, k2 = np.where(face != i)[0]
            print(k1, k2)
            q, r = vertices[face[k1]], vertices[face[k2]]
            # area de Voronoi
            area += (1 / 8) * (np.linalg.norm(p - r) ** 2 * cot(angles_tri[k1])
                               + np.linalg.norm(p - q) ** 2 * cot(angles_tri[k2]))
        else:  # es obtuso
            a, b, c = face
            a, b, c = vertices[a], vertices[b], vertices[c]
            ab, bc, ca = np.linalg.norm(a - b), np.linalg.norm(b - c), np.linalg.norm(c - a)
            s = (ab + bc + ca) / 2
            area_tri = (s * (s - ab) * (s - bc) * (s - ca)) ** (1 / 2)  # fórmula de Heron
            if angles_tri[k] > np.pi / 2:
                area += area_tri / 2
            else:
                area += area_tri / 4
    return area


# dados dos puntos, i y j, los cuales son vecinos y forman
# un borde. Calcula los ángulos opuestos en cada triángulo
def angles_triangles(i, j, neighbors, angles, faces):
    try:
        tri1, tri2 = neighbors[i, j]
    except KeyError:  # el orden en que se definen los bordes es desconocido
        tri1, tri2 = neighbors[j, i]
    angles1, angles2 = angles[tri1], angles[tri2]
    k1 = np.where(np.logical_and(faces[tri1] != i, faces[tri1] != j))
    k2 = np.where(np.logical_and(faces[tri2] != i, faces[tri2] != j))
    return angles1[k1], angles2[k2]


# retorna un arreglo donde cada índice indica una cara
# junto a sus 3 ángulos, en orden según el orden
# en que están sus vertices ordenados de menor a mayor
def get_angles(vertices, faces):
    angles = np.zeros((len(faces), 3))
    for i in range(0, len(faces)):
        a, b, c = faces[i]
        p, q, r = vertices[a], vertices[b], vertices[c]
        pq, qp = q - p, p - q
        pr, rp = r - p, p - r
        qr, rq = r - q, q - r
        angle1 = np.dot(pq, pr) / (np.linalg.norm(pq) * np.linalg.norm(pr))
        angle2 = np.dot(qp, qr) / (np.linalg.norm(qp) * np.linalg.norm(qr))
        angle3 = np.dot(rq, rp) / (np.linalg.norm(rq) * np.linalg.norm(rp))
        angle1, angle2, angle3 = np.arccos(angle1), np.arccos(angle2), np.arccos(angle3)
        #             a -> angle1
        #            / \
        # angle3 <- c _ b -> angle2
        angles[i] = [angle1, angle2, angle3]
    return angles


############################################################################
# se calculará el área y el ángulo a la vez para el caso gaussiano.
def angle_area(i, vertices, triangles):
    point = vertices[i]
    triangle = triangles[i]
    angle, area = 0, 0
    for j in range(0, len(triangle)):
        a, b = triangle[j]
        a, b = vertices[a], vertices[b]
        v1, v2 = a - point, b - point
        angle += np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
        e1, e2, e3 = np.linalg.norm(v1), np.linalg.norm(v2), np.linalg.norm(a - b)
        s = (e1 + e2 + e3) / 2
        area += (s * (s - e1) * (s - e2) * (s - e3)) ** (1 / 2)  # fórmula de Heron
    return angle, area


"""
3D Curvature Computation Module for Rusinkiewicz Algorithm

This module provides functions to compute 3D curvatures on a surface mesh based on the Rusinkiewicz algorithm.
The functions included allow for calculating the principal curvatures, diagonalizing curvature tensors, and 
projecting curvature components to new coordinate systems. Additional utilities for matrix decompositions and 
vertex area calculations are provided to support the curvature calculations.

Functions:
- ldltdc: Performs the LDL^T (Cholesky) decomposition.
- ldltsl: Solves a system of equations using LDL^T decomposition.
- rot_coord: Rotates coordinate directions to align with a new normal vector.
- proj_curves: Projects curvature coefficients to new directions.
- diagonalize_curves: Diagonalizes the curvature tensor to find principal curvatures and directions.
- proj_dcurves: Projects higher-order curvature terms to new directions.
- ldl_solve: Solves a system of linear equations using LDL decomposition.
- get_point_area: Calculates vertex area contributions in a mesh.

This module is useful for applications in 3D surface analysis, computer graphics, and geometric modeling 
that require accurate curvature computations on discrete surface meshes.
"""


def ldltdc(A, diag):
    """
    Performs the LDL^T decomposition (Cholesky decomposition) on matrix A.

    Parameters:
    - A (numpy.ndarray): The matrix to decompose.
    - diag (numpy.ndarray): An array to store the diagonal elements D of the decomposition.

    Returns:
    - bool: True if the decomposition is possible, False otherwise.
    """
    N = len(A)
    # si no es matriz, nada que hacer
    if N < 1:
        return False
    # caso especial, qué sucede si N <= 3
    elif N <= 3:
        d0 = A[0, 0]
        diag[0] = 1 / d0
        if N == 1:
            return d0 != 0
        A[1, 0] = A[0, 1]
        l10 = diag[0] * A[1, 0]
        d1 = A[1, 1] - l10 * A[1, 0]
        diag[1] = 1 / d1
        if N == 2:
            return d0 != 0 and d1 != 0
        d2 = A[2, 2] - diag[0] * A[2, 0] ** 2 - diag[1] * A[2, 1] ** 2
        diag[2] = 1 / d2
        A[2, 0] = A[0, 2]
        A[2, 1] = A[1, 2] - l10 * A[2, 0]
        return d0 != 0 and d1 != 0 and d2 != 0
    v = np.zeros(N - 1)
    for i in range(0, N):
        for k in range(0, i):
            v[k] = A[i, k] * diag[k]
        for j in range(i, N):
            sum = A[i, j]
            for k in range(0, i):
                sum -= v[k] * A[j, k]
            if i == j:
                if sum == 0:
                    return False
                diag[i] = 1 / sum
            else:
                A[j, i] = sum
    return True


def ldltsl(A, rdiag, b, same=False):
    """
    Solves the system of equations using forward and backward substitution on an LDL^T decomposed matrix.

    Parameters:
    - A (numpy.ndarray): The matrix after LDL^T decomposition.
    - rdiag (numpy.ndarray): Reciprocal of diagonal elements from the decomposition.
    - b (numpy.ndarray): The right-hand side vector.
    - same (bool): If True, solves Ly = x -> L^T x = y. Otherwise, solves Ly = b -> L^T x = y.

    Returns:
    - numpy.ndarray: Solution vector x.
    """
    x = b if same is True else np.zeros(len(A))
    for i in range(0, len(A)):
        sum = b[i]
        for k in range(0, i):
            sum -= A[i, k] * x[k]
        x[i] = sum * rdiag[i]
    for i in range(len(A) - 1, -1, -1):
        sum = 0
        for k in range(i + 1, len(A)):
            sum += A[k, i] * x[k]
        x[i] -= sum * rdiag[i]
    return x


def rot_coord(u_old, v_old, n_new):
    """
    Rotates coordinates u and v to align with a new normal vector.

    Parameters:
    - u_old (numpy.ndarray): Original u direction.
    - v_old (numpy.ndarray): Original v direction.
    - n_new (numpy.ndarray): New normal vector to align with.

    Returns:
    - tuple: The rotated u and v directions (u_new, v_new).
    """
    n_old = np.cross(u_old, v_old)
    n_dot = np.dot(n_old, n_new)
    u_new, v_new = u_old, v_old
    if n_dot <= -1:
        return -u_new, -v_new
    perp_old = n_new - n_dot * n_old
    d_perp = 1 / (1 + n_dot) * (n_old + n_new)
    u_new -= d_perp * np.dot(u_new, perp_old)
    v_new -= d_perp * np.dot(v_new, perp_old)
    return u_new, v_new


def proj_curves(u_old, v_old, ku_old, kuv_old, kv_old, u_new, v_new):
    """
    Projects curvatures from old directions to new directions.

    Parameters:
    - u_old, v_old (numpy.ndarray): Original directions.
    - ku_old, kuv_old, kv_old (float): Curvature coefficients in original directions.
    - u_new, v_new (numpy.ndarray): New directions for projection.

    Returns:
    - tuple: Curvature coefficients (ku_new, kuv_new, kv_new) in the new directions.
    """
    ru_new, rv_new = rot_coord(u_new, v_new, np.cross(u_old, v_old))
    u1, v1 = np.dot(ru_new, u_old), np.dot(ru_new, v_old)
    u2, v2 = np.dot(rv_new, u_old), np.dot(rv_new, v_old)
    # coeficientes de curvatura
    ku_new = ku_old * u1 * u1 + kuv_old * (2 * u1 * v1) + kv_old * v1 * v1
    kuv_new = ku_old * u1 * u2 + kuv_old * (u1 * v2 + u2 * v1) + kv_old * v1 * v2
    kv_new = ku_old * u2 * u2 + kuv_old * (2 * u2 * v2) + kv_old * v2 * v2
    return ku_new, kuv_new, kv_new


def diagonalize_curves(u_old, v_old, ku, kuv, kv, n_new):
    """
    Diagonalizes the curvature tensor to find principal curvatures and directions.

    Parameters:
    - u_old, v_old (numpy.ndarray): Original directions.
    - ku, kuv, kv (float): Curvature coefficients in original directions.
    - n_new (numpy.ndarray): New normal vector.

    Returns:
    - tuple: Principal curvatures (k1, k2) and principal directions (p1, p2).
    """
    ru_old, rv_old = rot_coord(u_old, v_old, n_new)
    c, s, tt = 1, 0, 0
    if kuv != 0:
        h = 0.5 * (kv - ku) / kuv
        tt = 1 / (h - np.sqrt(1 + h ** 2)) if h < 0 else 1 / (h + np.sqrt(1 + h ** 2))
        c = 1 / np.sqrt(1 + tt ** 2)
        s = tt * c
    k1 = ku - tt * kuv
    k2 = kv + tt * kuv
    if abs(k1) >= abs(k2):
        p1 = c * ru_old - s * rv_old
    else:
        k1, k2 = k2, k1
        p1 = s * ru_old + c * rv_old
    p2 = np.cross(n_new, p1)
    return k1, k2, p1, p2


def proj_dcurves(old_u, old_v, old_dcurv, new_u, new_v):
    """
    Projects higher-order curvature terms from old directions to new directions.

    Parameters:
    - old_u, old_v (numpy.ndarray): Original directions.
    - old_dcurv (numpy.ndarray): Higher-order curvature terms in original directions.
    - new_u, new_v (numpy.ndarray): New directions for projection.

    Returns:
    - numpy.ndarray: Projected higher-order curvature terms.
    """
    r_new_u, r_new_v = rot_coord(new_u, new_v, np.cross(old_u, old_v))
    u1 = np.dot(r_new_u, old_u)
    v1 = np.dot(r_new_u, old_v)
    u2 = np.dot(r_new_v, old_u)
    v2 = np.dot(r_new_v, old_v)
    new_dcurv = np.zeros(4)
    new_dcurv[0] = old_dcurv[0] * u1 ** 3 + old_dcurv[1] * 3 * u1 ** 2 * v1 + \
        old_dcurv[2] * 3 * v1 ** 2 * u1 + old_dcurv[3] * v1 ** 3
    new_dcurv[1] = old_dcurv[0] * u1 ** 2 * u2 + old_dcurv[1] * (u1 ** 2 * v2 + 2.0 * u2 * u1 * v1) + \
        old_dcurv[2] * (u2 * v1 ** 2 + 2 * u1 * v1 * v2) + old_dcurv[3] * v1 ** 2 * v2
    new_dcurv[2] = old_dcurv[0] * u1 * u2 ** 2 + old_dcurv[1] * (u2 ** 2 * v1 + 2.0 * u1 * u2 * v2) + \
        old_dcurv[2] * (u1 * v2 ** 2 + 2.0 * u2 * v2 * v1) + old_dcurv[3] * v1 * v2 ** 2
    new_dcurv[3] = old_dcurv[0] * u2 ** 3 + old_dcurv[1] * 3 * u2 ** 2 * v2 + \
        old_dcurv[2] * 3 * u2 * v2 ** 2 + old_dcurv[3] * v2 ** 3
    return new_dcurv


def ldl_solve(a, b):
    """
    Solves a system of linear equations using the LDL decomposition.

    Parameters:
    - a (numpy.ndarray): The decomposed matrix.
    - b (numpy.ndarray): The right-hand side vector.

    Returns:
    - numpy.ndarray: Solution vector x.
    """
    l = np.linalg.cholesky(a)
    # a = l * d * lT
    y = np.linalg.lstsq(l, b)[0]
    x = np.linalg.lstsq(l.T, y)[0]
    return x


def get_point_area(vertices, faces):
    """
    Calculates the area associated with each vertex in a mesh, using barycentric weights.

    Parameters:
    - vertices (numpy.ndarray): Array of vertex coordinates.
    - faces (numpy.ndarray): Array of face indices, each row representing a triangle.

    Returns:
    - tuple: (point_areas, corner_areas) where point_areas is an array of areas for each vertex,
      and corner_areas is an array of areas for each corner of each face.
    """
    point_areas = np.zeros(len(vertices))
    corner_areas = np.zeros((len(faces), 3))
    for i in range(0, len(faces)):
        edges = np.array([vertices[faces[i, 2]] - vertices[faces[i, 1]],
                          vertices[faces[i, 0]] - vertices[faces[i, 2]],
                          vertices[faces[i, 1]] - vertices[faces[i, 0]]])
        # cálculo de áreas por punto a través de pesos
        area = 0.5 * np.linalg.norm(np.cross(edges[0], edges[1]))
        l2 = np.array([np.linalg.norm(edges[0]) ** 2,
                       np.linalg.norm(edges[1]) ** 2,
                       np.linalg.norm(edges[2]) ** 2])
        bcw = np.array([l2[0] * (l2[1] + l2[2] - l2[0]),
                        l2[1] * (l2[2] + l2[0] - l2[1]),
                        l2[2] * (l2[0] + l2[1] - l2[2])])
        if bcw[0] <= 0:
            corner_areas[i, 1] = -0.25 * l2[2] * area / np.dot(edges[0], edges[2])
            corner_areas[i, 2] = -0.25 * l2[1] * area / np.dot(edges[0], edges[1])
            corner_areas[i, 0] = area - corner_areas[i, 1] - corner_areas[i, 2]
        elif bcw[1] <= 0:
            corner_areas[i, 2] = -0.25 * l2[0] * area / np.dot(edges[1], edges[0])
            corner_areas[i, 0] = -0.25 * l2[2] * area / np.dot(edges[1], edges[2])
            corner_areas[i, 1] = area - corner_areas[i, 2] - corner_areas[i, 0]
        elif bcw[2] <= 0:
            corner_areas[i, 0] = -0.25 * l2[1] * area / np.dot(edges[2], edges[1])
            corner_areas[i, 1] = -0.25 * l2[0] * area / np.dot(edges[2], edges[0])
            corner_areas[i, 2] = area - corner_areas[i, 0] - corner_areas[i, 1]
        else:
            scale = 0.5 * area / (bcw[0] + bcw[1] + bcw[2])
            for j in range(0, 3):
                corner_areas[i, j] = scale * (bcw[(j + 1) % 3] + bcw[(j - 1) % 3])
        point_areas[faces[i, 0]] += corner_areas[i, 0]
        point_areas[faces[i, 1]] += corner_areas[i, 1]
        point_areas[faces[i, 2]] += corner_areas[i, 2]
    return point_areas, corner_areas


"""
Custom Colormap Generation Module

This module provides functions to create custom colormaps and color scales based on input data.
The functions included allow users to generate color mappings based on data distribution, such as
positive/negative dominance or specific function values, and produce color scales using predefined
colormaps like RdBu. Each function allows for normalization and customization of the generated color scale.

Functions:
- colormap (DEPRECATED): Generates a colormap based on the proportion of positive, negative, and zero values in a function array.
- cmap_go (DEPRECATED): Normalizes curves and extracts key values for color scale mapping.
- create_standard_colorscale: Creates an RdBu color scale normalized between -a and a.
- get_colorscale: Extracts a subset of an RdBu color scale based on the specified value range.

The functions are useful for visualizing data distributions, highlighting key values in a dataset, and
creating visually distinct color gradients for plots and heatmaps.
"""


def colormap(function):
    """
    Generates a colormap based on the distribution of positive, negative, and zero values in a function array.

    Parameters:
    - function (numpy.ndarray): 1D array containing function values to analyze.

    Returns:
    - list: A list of RGB tuples representing the generated colormap.
    """
    fun_color = function / abs(np.max(function))
    more0, less0, cero = 0, 0, 0
    for i in range(0, len(fun_color)):
        if fun_color[i] > 0:
            more0 += 1
        if fun_color[i] < 0:
            less0 += 1
        if fun_color[i] == 0:
            cero += 1
    if less0 + cero == 0:  # h > 0
        cmap = [(1.00, 1.00, 0.00),  # amarillo
                (1.00, 0.00, 0.00)]  # rojo
    elif more0 + cero == 0:  # h < 0
        cmap = [(0.10, 0.43, 0.86),  # celeste
                (0.63, 0.78, 1.00)]  # azul
    elif less0 == 0 and cero != 0:  # h >= 0
        cmap = [(1.00, 1.00, 1.00),  # blanco
                (1.00, 1.00, 0.00),  # amarillo
                (1.00, 0.00, 0.00)]  # rojo
    elif more0 == 0 and cero != 0:  # h <= 0
        cmap = [(0.10, 0.43, 0.86),  # azul
                (0.63, 0.78, 1.00),  # celeste
                (1.00, 1.00, 1.00)]  # blanco
    elif less0 < more0 / 2:  # hay muchos más h > 0 que h < 0
        cmap = [(0.63, 0.78, 1.00),  # celeste
                (1.00, 1.00, 1.00),  # blanco
                (1.00, 1.00, 0.00),  # amarillo
                (1.00, 0.00, 0.00)]  # rojo
    elif more0 < less0 / 2:  # hay muchos más h < 0 que h > 0
        cmap = [(0.10, 0.43, 0.86),  # azul
                (0.63, 0.78, 1.00),  # celeste
                (1.00, 1.00, 1.00),  # blanco
                (1.00, 1.00, 0.00)]  # amarillo
    else:
        cmap = [(0.10, 0.43, 0.86),  # azul
                (0.63, 0.78, 1.00),  # celeste
                (1.00, 1.00, 1.00),  # blanco
                (1.00, 1.00, 0.00),  # amarillo
                (1.00, 0.00, 0.00)]  # rojo
    return cmap


def cmap_go(curves, max_value=None, min_value=None):
    """
    Normalizes a set of curves and determines key function values based on their distribution.

    Parameters:
    - curves (numpy.ndarray): Array of function values for which the color scale will be generated.
    - max_value (float, optional): Maximum value to set for normalization.
    - min_value (float, optional): Minimum value to set for normalization.

    Returns:
    - list: A list of normalized values to map to a color scale.
    """
    curves_norm = curves / np.max(np.abs(curves))
    index = np.zeros((5, len(curves)))
    for i in range(0, len(curves)):
        index[0, i] = abs(curves_norm[i] + 1)
        index[1, i] = abs(curves_norm[i] + 0.5)
        index[2, i] = abs(curves_norm[i])
        index[3, i] = abs(curves_norm[i] - 0.5)
        index[4, i] = abs(curves_norm[i] - 1)
    func = [curves[np.argmin(index[j])] for j in range(0, 5)]
    func_sort = sorted(func)
    if max_value is not None:
        func_sort[-1] = max_value
    if min_value is not None:
        func_sort[0] = min_value
    values = [(v - func_sort[0]) / (func_sort[-1] - func_sort[0]) for v in func_sort]
    for i in range(0, len(values)):
        if values[i] > 1.0:
            values[i] = 1.0
    return values


def create_standard_colorscale(a, num_colors=500):
    """
    Creates a standard color scale based on the RdBu colormap, normalized between -a and a.

    Parameters:
    - a (float): Maximum absolute value for normalization.
    - num_colors (int): Number of color entries in the generated colorscale.

    Returns:
    - list: A list of color hex codes representing the colorscale from -a to a.
    """
    # Normalizamos la escala entre -a y a
    norm = Normalize(vmin=-a, vmax=a)
    # Usamos una escala de matplotlib RdBu, convertida a hexadecimal
    colors = [to_hex(cm.RdBu(norm(x))) for x in np.linspace(-a, a, num_colors)]
    return colors[::-1]


def get_colorscale(values, a=8, num_colors=256):
    """
    Extracts a subset of the standard RdBu color scale based on a specified value range.

    Parameters:
    - values (numpy.ndarray): Array of values used to determine the color range.
    - a (float): Maximum absolute value for normalization.
    - num_colors (int): Number of color entries in the extracted colorscale subset.

    Returns:
    - list: A list of [position, color] pairs representing the subset of the colorscale.
    """
    b, c = np.min(values), np.max(values)
    norm = Normalize(vmin=-a, vmax=a)
    # Usamos una escala de matplotlib RdBu, convertida a hexadecimal
    colorscale = [to_hex(cm.RdBu(norm(x))) for x in np.linspace(-a, a, num_colors)][::-1]
    # Mapear el rango b a c en el espacio -a a a
    start_idx = int((b + a) / (2 * a) * (len(colorscale) - 1))
    end_idx = int((c + a) / (2 * a) * (len(colorscale) - 1))

    # Extraer el subconjunto relevante
    sub_colorscale = colorscale[start_idx:end_idx + 1]

    # Interpolar si el número de colores no coincide con num_colors
    if len(sub_colorscale) != num_colors:
        indices = np.linspace(0, len(sub_colorscale) - 1, num_colors).astype(int)
        sub_colorscale = [sub_colorscale[i] for i in indices]

    return [[i / (num_colors - 1), color] for i, color in enumerate(sub_colorscale)]


"""
3D Geometric Shape Generation and Saving Module

This module provides functions to generate various 3D geometric shapes (tube, sphere, ellipsoid, cylinder, and torus)
within a voxel grid represented as a 3D NumPy array. Each function allows for customization of shape dimensions, 
voxel sizes, and placement within the 3D grid. Additionally, the shapes can be saved as TIFF images, either as a 
single 3D file or as individual slices for each z-plane.

Functions:
- make_tube (DEPRECATED): Generates a cylindrical tube in a 3D grid and optionally saves it as a TIFF file.
- make_sphere: Generates a sphere in a 3D grid and optionally saves it as a TIFF file.
- save_a_lot_of_spheres: Creates multiple spheres of varying radii in a single volume and saves it as a TIFF file.
- make_ellipsoid: Generates an ellipsoid in a 3D grid and optionally saves it as a TIFF file.
- make_cylinder: Generates a cylinder in a 3D grid and optionally saves it as a TIFF file.
- make_torus: Generates a torus in a 3D grid and optionally saves it as a TIFF file.

Each function includes options to:
- Customize the size and placement of the shape within the 3D grid.
- Specify the size of each voxel to control the scaling of the shapes.
- Save the generated 3D grid as either a single TIFF file or individual slice TIFF images.

All functions return the generated 3D volume as a NumPy array. This module is suitable for applications in 3D imaging,
modeling, or simulation that require voxel-based representations of geometric shapes.
"""


def make_tube(path: str,
              r: int,
              h: int,
              dim=(128, 128, 64),
              one_tiff=True,
              name=False):
    """
    Creates a 3D cylindrical tube within a stack of images.

    Parameters:
    - path (str): Directory path for saving the output.
    - r (int): Radius of the tube.
    - h (int): Height of the tube.
    - dim (tuple): Dimensions of the volume as (x, y, z).
    - one_tiff (bool): If True, saves the volume as a single TIFF file.
    - name (bool): If True, uses the specified path for the filename.

    Returns:
    - numpy.ndarray: A 3D array representing the tube volume.
    """
    x_dim, y_dim, z_dim = dim
    if r > x_dim:
        print('radio mayor a las dimensiones')
        x_dim = r
    if r > y_dim:
        y_dim = r
        print('radio mayor a las dimensiones')
    if h > z_dim:
        z_dim = h
        print('altura mayor a las dimensiones')
    zl, zh = z_dim // 2 - h / 2, z_dim // 2 + h / 2 - 1
    tiff = np.zeros((z_dim, y_dim, x_dim))
    for k in range(0, z_dim):
        if k > zh or k < zl:
            continue
        for i in range(0, x_dim):
            for j in range(0, y_dim):
                if (i - x_dim // 2) ** 2 + (j - y_dim // 2) ** 2 <= r ** 2:
                    tiff[k, j, i] = 255
    if not os.path.isdir(path):
        os.makedirs(path)
    path += r'\tube_' + f'r{r}_h{h}_x{x_dim}_y{y_dim}_z{z_dim}' if name is False else path
    if one_tiff is True:
        path += '.tif'
        imwrite(path, tiff.astype(np.uint8))
    else:
        for k in range(0, z_dim):
            sub_path = path + f'_t000_ch_000_z00{k}.tif'
            imwrite(sub_path, tiff[k].astype(np.uint8))
    return tiff


def make_sphere(r: int,
                dim=(-1, -1, -1),
                voxel_size=(1.0, 1.0, 1.0),
                c=None,
                one_tiff=True,
                path=r'C:\\',
                save=False):
    """
    Creates a 3D sphere within a stack of images.

    Parameters:
    - r (int): Radius of the sphere.
    - dim (tuple): Dimensions of the volume as (x, y, z).
    - voxel_size (tuple): Size of each voxel in (x, y, z).
    - c (tuple): Center coordinates of the sphere (if None, defaults to center).
    - one_tiff (bool): If True, saves the volume as a single TIFF file.
    - path (str): Directory path for saving the output.
    - save (bool): If True, saves the generated volume as a TIFF file.

    Returns:
    - numpy.ndarray: A 3D array representing the sphere volume.
    """
    if 2 * r > dim[0] or 2 * r > dim[1] or 2 * r > dim[2]:
        if (2 * r + 1) % 2 == 1:
            dim = ((2 * r + 2) // voxel_size[0],
                   (2 * r + 2) // voxel_size[1],
                   (2 * r + 2) // voxel_size[2])
        else:
            dim = ((2 * r + 1) // voxel_size[0],
                   (2 * r + 1) // voxel_size[1],
                   (2 * r + 1) // voxel_size[2])
        dim = (int(dim[0]), int(dim[1]), int(dim[2]))
    tiff = np.zeros(dim)
    if c is None:
        c = (voxel_size[0] * dim[0] / 2,
             voxel_size[1] * dim[1] / 2,
             voxel_size[2] * dim[2] / 2)
    for z in range(dim[0]):
        for y in range(dim[1]):
            for x in range(dim[2]):
                if (x * voxel_size[2] - c[2]) ** 2 + \
                        (y * voxel_size[1] - c[1]) ** 2 + \
                        (z * voxel_size[0] - c[0]) ** 2 <= r ** 2:
                    tiff[z, y, x] = 255
    if save is True:
        if not os.path.isdir(path):
            os.makedirs(path)
        path += r'\sphere_' + f'r{r}'
        if one_tiff is True:
            path += '.tif'
            imwrite(path, tiff.astype(np.uint8))
        else:
            for k in range(0, dim[0]):
                sub_path = path + f'_t000_ch_000_z00{k}.tif'
                imwrite(sub_path, tiff[k].astype(np.uint8))
    return tiff


def save_a_lot_of_spheres(path, name, x, one_tiff=False):
    """
    Creates and saves multiple spheres in a single image stack, aligned in sequence.

    Parameters:
    - path (str): Directory path for saving the output.
    - name (str): Base name for the saved file.
    - x (list): List of sphere radii.
    - one_tiff (bool): If True, saves the volume as a single TIFF file.
    """
    max_sphere = make_sphere(x[-1], voxel_size=(0.5, 0.163, 0.163))
    z, m, n = np.shape(max_sphere)
    spheres = np.zeros((z, m, n * len(x)))
    v = 0
    for rad in x:
        if rad == x[-1]:
            sphere = max_sphere
        else:
            sphere = make_sphere(rad, dim=(z, m, n), voxel_size=(0.5, 0.163, 0.163))
        for k in range(z):
            for j in range(m):
                for i in range(n):
                    spheres[k, j, i + v] = sphere[k, j, i]
        v += n
    if not os.path.isdir(path):
        os.makedirs(path)
    path += r'\\' + name
    if one_tiff is True:
        path += '.tif'
        imwrite(path, spheres.astype(np.uint8))
    else:
        for k in range(0, z):
            sub_path = path + f'_ch_000_z00{k}.tif'
            imwrite(sub_path, spheres[k].astype(np.uint8))


def make_ellipsoid(radii: tuple,
                   dim=(-1, -1, -1),
                   voxel_size=(1.0, 1.0, 1.0),
                   c=None,
                   one_tiff=True,
                   path=r'C:\\',
                   save=False,
                   debug_dim=False):
    """
    Creates a 3D ellipsoid within a stack of images.

    Parameters:
    - radii (tuple): Radii of the ellipsoid along (z, y, x).
    - dim (tuple): Dimensions of the volume as (z, y, x).
    - voxel_size (tuple): Size of each voxel in (z, y, x).
    - c (tuple): Center coordinates of the ellipsoid (if None, defaults to center).
    - one_tiff (bool): If True, saves the volume as a single TIFF file.
    - path (str): Directory path for saving the output.
    - save (bool): If True, saves the generated volume as a TIFF file.

    Returns:
    - numpy.ndarray: A 3D array representing the ellipsoid volume.
    """
    dim = list(dim)
    r_z, r_y, r_x = radii

    if 2 * r_z > dim[0] or 2 * r_y > dim[1] or 2 * r_x > dim[2]:
        if (2 * r_z + 1) % 2 == 1:
            dim[0] = (2 * r_z + 2) // voxel_size[0]
        else:
            dim[0] = (2 * r_z + 1) // voxel_size[0]
        if (2 * r_y + 1) % 2 == 1:
            dim[1] = (2 * r_y + 2) // voxel_size[1]
        else:
            dim[1] = (2 * r_y + 1) // voxel_size[1]
        if (2 * r_x + 1) % 2 == 1:
            dim[2] = (2 * r_x + 2) // voxel_size[2]
        else:
            dim[2] = (2 * r_x + 1) // voxel_size[2]

        dim = (int(dim[0]), int(dim[1]), int(dim[2]))

    # Crear el volumen en 3D (vacío) con las dimensiones ajustadas
    tiff = np.zeros(dim)

    # Ajustar el centro si no se especifica
    if c is None:
        c = [voxel_size[0] * dim[0] / 2,
             voxel_size[1] * dim[1] / 2,
             voxel_size[2] * dim[2] / 2]

    if debug_dim:
        return dim, c

    # Rellenar el volumen con el elipsoide
    for z in range(dim[0]):
        for y in range(dim[1]):
            for x in range(dim[2]):
                if (x * voxel_size[2] - c[2]) ** 2 / r_x ** 2 + \
                        (y * voxel_size[1] - c[1]) ** 2 / r_y ** 2 + \
                        (z * voxel_size[0] - c[0]) ** 2 / r_z ** 2 <= 1:
                    tiff[z, y, x] = 255

    # Guardar el volumen si se solicita
    if save is True:
        if not os.path.isdir(path):
            os.makedirs(path)
        path += r'\ellipsoid_' + f'rx{r_x}_ry{r_y}_rz{r_z}'
        if one_tiff is True:
            path += '.tif'
            imwrite(path, tiff.astype(np.uint8))
        else:
            for k in range(0, dim[0]):
                sub_path = path + f'_t000_ch_000_z00{k}.tif'
                imwrite(sub_path, tiff[k].astype(np.uint8))

    return tiff


def make_cylinder(radius: int,
                  height: int,
                  dim=(-1, -1, -1),
                  voxel_size=(1.0, 1.0, 1.0),
                  c=None,
                  one_tiff=True,
                  path=r'C:\\',
                  save=False):
    """
    Creates a 3D cylinder within a stack of images.

    Parameters:
    - radius (int): Radius of the cylinder.
    - height (int): Height of the cylinder.
    - dim (tuple): Dimensions of the volume as (z, y, x).
    - voxel_size (tuple): Size of each voxel in (z, y, x).
    - c (tuple): Center coordinates of the cylinder (if None, defaults to center).
    - one_tiff (bool): If True, saves the volume as a single TIFF file.
    - path (str): Directory path for saving the output.
    - save (bool): If True, saves the generated volume as a TIFF file.

    Returns:
    - numpy.ndarray: A 3D array representing the cylinder volume.
    """
    dim = list(dim)

    # Ajustar las dimensiones si es necesario
    if height > dim[0] or 2 * radius > dim[1] or 2 * radius > dim[2]:
        if (2 * radius + 1) % 2 == 1:
            dim[1] = (2 * radius + 2) // voxel_size[1]
            dim[2] = (2 * radius + 2) // voxel_size[2]
        else:
            dim[1] = (2 * radius + 1) // voxel_size[1]
            dim[2] = (2 * radius + 1) // voxel_size[2]
        dim[0] = height // voxel_size[0]

    dim = (int(dim[0]), int(dim[1]), int(dim[2]))
    # Inicializar el volumen en 3D
    tiff = np.zeros(dim, dtype=np.uint8)

    # Ajustar el centro si no se especifica
    if c is None:
        c = [dim[0] * voxel_size[0] / 2,
             dim[1] * voxel_size[1] / 2,
             dim[2] * voxel_size[2] / 2]

    # Llenar el volumen con el cilindro
    for z in range(dim[0]):
        for y in range(dim[1]):
            for x in range(dim[2]):
                if (x * voxel_size[2] - c[2]) ** 2 + \
                        (y * voxel_size[1] - c[1]) ** 2 <= radius ** 2 and \
                        abs(z * voxel_size[0] - c[0]) <= height / 2:
                    tiff[z, y, x] = 255

    # Guardar el volumen si se solicita
    if save:
        if not os.path.isdir(path):
            os.makedirs(path)
        path += f'\\cylinder_r{radius}_h{height}.tif'
        if one_tiff:
            imwrite(path, tiff)
        else:
            for k in range(dim[0]):
                slice_path = f"{path.rstrip('.tif')}_z{k:03d}.tif"
                imwrite(slice_path, tiff[k])

    return tiff


def make_torus(inner_radius: float,
               tube_radius: float,
               dim=(-1, -1, -1),
               voxel_size=(1.0, 1.0, 1.0),
               c=None,
               one_tiff=True,
               path=r'C:\\',
               save=False,
               return_center=False,
               debug_dim=False):
    """
    Creates a 3D torus within a stack of images.

    Parameters:
    - inner_radius (float): Radius from the torus center to the center of the tube.
    - tube_radius (float): Radius of the tube of the torus.
    - dim (tuple): Dimensions of the volume as (z, y, x).
    - voxel_size (tuple): Size of each voxel in (z, y, x).
    - c (tuple): Center coordinates of the torus (if None, defaults to center).
    - one_tiff (bool): If True, saves the volume as a single TIFF file.
    - path (str): Directory path for saving the output.
    - save (bool): If True, saves the generated volume as a TIFF file.
    - return_center (bool): If True, returns the center coordinates along with the volume.
    - debug_dim (bool): If True, returns only the dimensions and center coordinates.

    Returns:
    - numpy.ndarray: A 3D array representing the torus volume, or a tuple with the volume and center if return_center is True.
    """
    dim = list(dim)

    # Ajustar las dimensiones si es necesario
    if 2 * tube_radius > dim[0] or \
            2 * (inner_radius + tube_radius) > dim[1] or \
            2 * (inner_radius + tube_radius) > dim[2]:
        if (2 * tube_radius + 1) % 2 == 1:
            dim[0] = (2 * tube_radius + 2) // voxel_size[0]
        else:
            dim[0] = (2 * tube_radius + 1) // voxel_size[0]
        if (2 * (inner_radius + tube_radius) + 1) % 2 == 1:
            dim[1] = (2 * (inner_radius + tube_radius) + 2) // voxel_size[1]
            dim[2] = (2 * (inner_radius + tube_radius) + 2) // voxel_size[2]
        else:
            dim[1] = (2 * (inner_radius + tube_radius) + 1) // voxel_size[1]
            dim[2] = (2 * (inner_radius + tube_radius) + 1) // voxel_size[2]

        dim = (int(dim[0]) + 1, int(dim[1]) + 1, int(dim[2]) + 1)

    # Crear el volumen en 3D
    tiff = np.zeros(dim, dtype=np.uint8)

    # Ajustar el centro si no se especifica
    if c is None:
        c = [dim[0] * voxel_size[0] / 2,  # Centro en Z
             dim[1] * voxel_size[1] / 2,  # Centro en Y
             dim[2] * voxel_size[2] / 2]  # Centro en X

    if debug_dim:
        return dim, c

    # Llenar el volumen con el toroide
    for z in range(dim[0]):
        for y in range(dim[1]):
            for x in range(dim[2]):
                # Calcular las coordenadas en el plano radial
                y_dist = (y * voxel_size[1] - c[1])
                x_dist = (x * voxel_size[2] - c[2])
                radial_dist = np.sqrt(y_dist ** 2 + x_dist ** 2)

                # Calcular la distancia desde el círculo central del toroide
                if (radial_dist - inner_radius) ** 2 + (z * voxel_size[0] - c[0]) ** 2 <= tube_radius ** 2:
                    tiff[z, y, x] = 255

    # Guardar el volumen si se solicita
    if save:
        if not os.path.isdir(path):
            os.makedirs(path)
        path += f'\\torus_ir{inner_radius}_tr{tube_radius}.tif'
        if one_tiff:
            imwrite(path, tiff)
        else:
            for k in range(dim[0]):
                slice_path = f"{path.rstrip('.tif')}_z{k:03d}.tif"
                imwrite(slice_path, tiff[k])

    return (tiff, c) if return_center else tiff


def torus_hyperstack(I, J, path=None):
    dims, centers, dx = [0, 0, 0], [], 0
    for (inner_radius, tube_radius) in zip(I, J):
        dim, c = make_torus(inner_radius, tube_radius, voxel_size=(0.5, 0.163, 0.163), debug_dim=True)
        dims[2] += dim[2] + 10
        if dims[0] == 0:
            dims[1] += dim[1]
            dims[0] += dim[0]
        centers.append([dims[0] / 2 * 0.5, c[1], dx + c[2]])
        dx += dim[2] * 0.163
    torus = np.zeros(dims)
    k = 0
    for (inner_radius, tube_radius) in zip(I, J):
        torus += make_torus(inner_radius, tube_radius, voxel_size=(0.5, 0.163, 0.163), dim=dims, c=centers[k])
        k += 1
        print(k)
    if path is not None:
        for k in range(len(torus)):
            slice_path = f"{path.rstrip('.tif')}_z{k:03d}.tif"
            imwrite(slice_path, torus[k])
    return torus


def ellipsoids_hyperstack(I, J, K, path=None):
    dims, centers, dx, dy = [0, 0, 0], [], 0, 0
    for (a, b, c) in zip(I, J, K):
        dim, c = make_ellipsoid((c, b, a), voxel_size=(0.5, 0.163, 0.163), debug_dim=True)
        dims[2] += dim[2] + 1
        dims[1] += dim[1] + 1
        if dim[0] > dims[0]:
            dims[0] = dim[0]
        centers.append([dims[0] / 2 * 0.5, dy + c[1], dx + c[2]])
        dx += dim[2] * 0.163
        dy += dim[1] * 0.163
    ellipsoid = np.zeros(dims)
    k = 0
    for (a, b, c) in zip(I, J, K):
        ellipsoid += make_ellipsoid((c, b, a), voxel_size=(0.5, 0.163, 0.163), dim=dims, c=centers[k])
        k += 1
        print(k)
    if path is not None:
        for k in range(len(ellipsoid)):
            slice_path = f"{path.rstrip('.tif')}_z{k:03d}.tif"
            imwrite(slice_path, ellipsoid[k])
    return ellipsoid


def theoretic_curvature_ellipsoid(x, y, z, a, b, c):
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    th = np.arctan2(y, x)
    phi = np.arccos(z / r)
    return 1 / 2 * (1 / (a ** 2 * np.sin(th) ** 2 + b ** 2 * np.cos(th) ** 2) +
                    1 / (c ** 2 * np.sin(phi) ** 2 + a ** 2 * np.cos(phi) ** 2))


def theoretic_curvature_torus(x, y, z, c, ir, tr):
    dist_rad = np.sqrt((x - c[2]) ** 2 + (y - c[1]) ** 2)
    th = np.arctan2(z - c[0], dist_rad - ir)
    return 1 / 2 * (1 / tr + np.cos(th) / ir)


if __name__ == '__main__':
    I = np.ones(20) * 2
    K = np.linspace(2, 10, 20)
    ellipsoids_hyperstack(I, I, K, path=r'C:\Users\matia\Desktop\ellipsoids\ellipsoid_z')
    ellipsoids_hyperstack(I, K, I, path=r'C:\Users\matia\Desktop\ellipsoids\ellipsoid_y')
    ellipsoids_hyperstack(K, I, I, path=r'C:\Users\matia\Desktop\ellipsoids\ellipsoid_x')

