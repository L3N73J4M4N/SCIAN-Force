import plotly.graph_objects as go
from skimage import io
from skimage import measure
import trimesh
import pymeshfix
import os
import numpy as np
from forces import Forces3D
from time import time


class Mesh3D:
    def __init__(self, path=None, stack=None, mesh=None, objects=None, messages=True,
                 spacing=(0.5, 0.163, 0.163), n=20, fix=True, split=True):
        self.pm = messages
        t0 = time()
        if stack is None and path is not None:
            stack = io.imread(path)
        if mesh is None and stack is not None:
            vertices, faces, _, _ = measure.marching_cubes(stack, 0, allow_degenerate=0, spacing=spacing)
        else:
            vertices, faces = mesh
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
        if split is True:
            meshes = mesh.split()
        else:
            meshes = [mesh]
        meshes_split = []
        self.name = os.path.split(path)[1] if path is not None else 'mesh'
        if self.pm:
            print('\n--------------------------')
            print(self.name + ' - loaded', "--- %s seconds ---" % abs(time() - t0))
            print('\n--------------------------')
        for i in range(0, len(meshes)):
            v, f = meshes[i].vertices, meshes[i].faces
            if len(v) < n:
                continue
            if fix is True:
                mesh_fix = pymeshfix.MeshFix(v, f)
                mesh_fix.repair()
                v, f = mesh_fix.v, mesh_fix.f
            meshes_split.append([v, f])
            if self.pm:
                print('obj_' + str(i) + ': mesh generated', "--- %s seconds ---" % abs(time() - t0))
        self.meshes = []
        if objects is None:
            self.meshes = meshes_split
            self.n = len(self.meshes)
            string = 'all'
        else:
            for i in range(0, len(meshes_split)):
                if i in objects:
                    self.meshes.append(meshes_split[i])
            self.n = len(self.meshes)
            string = ''
            for i in range(0, self.n):
                if i > 0:
                    string += ', ' + str(objects[i])
                else:
                    string += str(objects[i])
        if self.pm:
            print('\n--------------------------')
            print('total objects =', len(meshes))
            print('objects =', len(meshes_split))
            print('deleted objects =', len(meshes) - len(meshes_split))
            print('work with ' + str(self.n) + ' objects = ' + string)

    def smooth(self, lamb=0.6, nu=0.5, iterations=10, laplacian=True):
        t0 = time()
        if self.pm:
            print('\n--------------------------')
        for i in range(0, self.n):
            v, f = self.meshes[i]
            mesh_smoothed = trimesh.Trimesh(v, f)
            if laplacian is True:
                mesh_smoothed = trimesh.smoothing.filter_laplacian(mesh_smoothed,
                                                                   lamb=lamb,
                                                                   iterations=iterations)
            else:
                mesh_smoothed = trimesh.smoothing.filter_taubin(mesh_smoothed,
                                                                lamb=lamb,
                                                                nu=nu,
                                                                iterations=iterations)
            self.meshes[i] = [mesh_smoothed.vertices, mesh_smoothed.faces]
            if self.pm:
                print('obj_' + str(i) + ': smoothed', "--- %s seconds ---" % abs(time() - t0))
        return self.meshes

    def plot(self, curvatures=None, title='default', edges=False, rotate=True, show=True, fig=None):
        title = self.name if title == 'default' else title
        if fig is None:
            fig = go.Figure(layout={'title': title})
        x_max, y_max, z_max = -np.inf, -np.inf, -np.inf
        x_min, y_min, z_min = np.inf, np.inf, np.inf
        for i in range(0, self.n):
            name = 'obj_' + str(i)
            if curvatures is not None:
                curve = curvatures[i]
                text = np.ndarray.astype(curve, str)
                for j in range(0, len(text)):
                    text[j] += '<br>' + name
            v, f = self.meshes[i]
            if rotate is True:
                x, y, z = v[:, 2], v[:, 1], v[:, 0]
            else:
                x, y, z = v[:, 0], v[:, 1], v[:, 2]
            xm0, ym0, zm0 = np.max(x), np.max(y), np.max(z)
            xm1, ym1, zm1 = np.min(x), np.min(y), np.min(z)
            x_max = xm0 if x_max < xm0 else x_max
            y_max = ym0 if y_max < ym0 else y_max
            z_max = zm0 if z_max < zm0 else z_max
            x_min = xm1 if x_min > xm1 else x_min
            y_min = ym1 if y_min > ym1 else y_min
            z_min = zm1 if z_min > zm1 else z_min
            i, j, k = f[:, 0], f[:, 1], f[:, 2]
            if curvatures is None:
                fig.add_mesh3d(x=x, y=y, z=z, i=i, j=j, k=k, color='red', name=name)
            else:
                fig.add_mesh3d(x=x, y=y, z=z, i=i, j=j, k=k,
                               intensity=curve,
                               coloraxis='coloraxis',
                               text=text)
            if edges is True:
                tri_points = v[f]
                xe, ye, ze = [], [], []
                for t in tri_points:
                    xe.extend([t[k % 3][0] for k in range(4)] + [None])
                    ye.extend([t[k % 3][1] for k in range(4)] + [None])
                    ze.extend([t[k % 3][2] for k in range(4)] + [None])
                fig.add_scatter3d(x=ze, y=ye, z=xe, mode='lines', name='',
                                  line=dict(color='rgb(70,70,70)', width=1))
        xyz = np.array([x_max - x_min, y_max - y_min, z_max - z_min])
        xyz /= np.max(xyz)
        fig.layout.scene.aspectratio = {'x': xyz[0], 'y': xyz[1], 'z': xyz[2]}
        if curvatures is not None:
            fig.update_layout(coloraxis={'colorscale': [[0, 'rgb(30,100,170)'],
                                                        [0.5, 'rgb(255,255,255)'],
                                                        [1, 'rgb(200,35,35)']]})
        if show:
            fig.show()
        fig.update_layout(
            scene=dict(
                xaxis=dict(visible=False),  # Oculta el eje x
                yaxis=dict(visible=False),  # Oculta el eje y
                zaxis=dict(visible=False)  # Oculta el eje z
            )
        )
        return fig

    def save_off(self, i_plot=-1, filename='default'):
        t0 = time()
        filename = self.name if filename == 'default' else str(filename)
        if i_plot == -1:
            i_plot = self.n
        for i in range(0, i_plot):
            v, f = self.meshes[i]
            with open(filename + '_obj_' + str(i) + '.off', 'w') as file:
                file.write('OFF\n')
                file.write(' ' + str(len(v)) +
                           ' ' + str(len(f)) +
                           ' 0\n')
                for vertex in v:
                    file.write('      ' +
                               str(vertex[2]) + '       ' +
                               str(vertex[1]) + '       ' +
                               str(vertex[0]) + '\n')
                for face in f:
                    file.write('           3          ' +
                               str(face[0]) + '          ' +
                               str(face[1]) + '          ' +
                               str(face[2]) + '\n')
            if self.pm:
                print('obj_' + str(i) + ': saved', "--- %s seconds ---" % abs(time() - t0))


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


def off2mesh(namepath):
    with open(namepath, 'r') as path:
        lines = path.readlines()
        if lines[0].strip() != 'OFF':
            raise ValueError('The file have not OFF format.')
        n, m, _ = list(map(int, lines[1].split()))
        vertices = np.zeros((n, 3))
        faces = np.zeros((m, 3), int)
        for i in range(2, n + 2):
            line = str2num(lines[i].strip())
            vertices[i - 2] = [float(line[0]), float(line[1]), float(line[2])]
        for i in range(n + 2, n + m + 2):
            line = str2num(lines[i].strip(), 'int')
            faces[i - n - 2] = [int(line[1]), int(line[2]), int(line[3])]
        return vertices, faces


def plot(mesh, curves=None, title='title', edges=False, rotate=False, show=True, fig=None):
    v, f = mesh
    if rotate is True:
        x, y, z = v[:, 2], v[:, 1], v[:, 0]
    else:
        x, y, z = v[:, 0], v[:, 1], v[:, 2]
    i, j, k = f[:, 0], f[:, 1], f[:, 2]
    if fig is None:
        fig = go.Figure(layout={'title': title})
    if curves is None:
        fig.add_mesh3d(x=x, y=y, z=z, i=i, j=j, k=k, color='red')
    else:
        fig.add_mesh3d(x=x, y=y, z=z, i=i, j=j, k=k,
                       intensity=curves,
                       coloraxis='coloraxis',
                       text=np.ndarray.astype(curves, str))
    if edges is True:
        tri_points = v[f]
        xe, ye, ze = [], [], []
        for t in tri_points:
            xe.extend([t[k % 3][0] for k in range(4)] + [None])
            ye.extend([t[k % 3][1] for k in range(4)] + [None])
            ze.extend([t[k % 3][2] for k in range(4)] + [None])
        if rotate is True:
            xe, ze = ze, xe
        fig.add_scatter3d(x=xe, y=ye, z=ze, mode='lines', name='',
                          line=dict(color='rgb(70,70,70)', width=1))
    fig.update_layout(coloraxis={'colorscale': [[0, 'rgb(30,100,170)'],
                                                [0.5, 'rgb(255,255,255)'],
                                                [1, 'rgb(200,35,35)']]})
    if show:
        fig.show()
    return fig


def obj2mesh(path):
    vertices, faces = [], []
    with open(path, 'r') as file:
        for lines in file:
            line = lines.split()
            if not line:
                continue
            elif line[0] == 'v':
                vertices.append([float(line[1]), float(line[2]), float(line[3])])
            elif line[0] == 'f':
                faces.append([int(line[1]) - 1, int(line[2]) - 1, int(line[3]) - 1])
    return np.array(vertices), np.array(faces)


def nCurvatures(meshes, method='r', umbral=2.5):
    t0 = time()
    n = len(meshes)
    curves = []
    print('\n--------------------------')
    for i in range(0, n):
        mesh = meshes[i]
        curv_mesh = Forces3D(mesh=mesh)
        if method == 'r':
            curv_mesh.rusinkiewicz_curvature()
        else:
            curv_mesh.laplacian_curvature()
        mean_curvature = curv_mesh.meanCurvature
        for j in range(0, len(mean_curvature)):
            if mean_curvature[j] > umbral:
                mean_curvature[j] = umbral
            if mean_curvature[j] < -umbral:
                mean_curvature[j] = -umbral
        curves.append(mean_curvature)
        print('obj_' + str(i) + ': curvatures calculated', "--- %s seconds ---" % abs(time() - t0))
    return curves


if __name__ == '__main__':
    m = Mesh3D(path=r"C:\Users\matia\Desktop\SCIAN\stacks\droplets_ROI_02.tif")
    meshes = m.smooth()
    m.plot(curvatures=nCurvatures(meshes=meshes))
