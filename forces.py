import numpy as np
import plotly.figure_factory as ff
import plotly.graph_objects as go
import os
import pymeshfix
import trimesh
from time import time
import utility_forces as uf


class Forces3D:
    """
    Class for calculating forces and geometric properties in 3D for meshes.
    This class provides methods to obtain normals, neighbors, areas, curvatures,
    and other relevant parameters for the 3D geometry of a mesh. It also enables
    reading mesh files, plotting, and coloring objects based on geometric and
    physical properties.
    """

    def __init__(self, path=None, mesh=None, interval='1x1x1', units='px', rotate=False):
        """
        Initializes an instance of Forces3D.

        Parameters
        ----------
        path : str, optional
            Path to the mesh file.
        mesh : tuple, optional
            Tuple containing the vertices and faces of the mesh.
        interval : str, optional
            Scaling factor for the vertices.
        units : str, optional
            Unit of measure for the mesh (default is 'px').
        rotate : bool, optional
            Whether to rotate the mesh.
        """

        self.path = path  # ruta del archivo
        self.colorFun = None  # array que le da color a cada cara
        self.normalsFaces = None  # array 3D con las normales unitarias de cada cara
        self.normalsPoints = None  # array 3D con las normales unitarias promedio en cada punto
        self.neighborsFaces = None  # dict - neighborsFaces[edge = (p, q)] -> [caras adyacentes a edge]
        self.neighborsPoints = None  # dict - neighborsFaces[point] -> [puntos conectados con point]
        self.meanCurvature = None  # array con la curvatura promedio para cada punto
        self.method = None  # str con el nombre del método utilizado para la curvatura
        self.trianglePoints = None  # dict - trianglePoints[v] -> [[otros v en face_i], ...]
        self.areas = None  # áreas de los triángulos / 3
        self.gaussCurvature = None  # array con la curvatura gaussiana para cada punto
        self.principals = None  # curvaturas principales, k1 y k2
        # array con malla de la superficie, vertices y caras
        if mesh is not None:
            self.name = self.path  # nombre del archivo
            if self.name is None:
                self.name = 'unnamed'
            vertices, faces = mesh
        else:
            self.name = os.path.split(path)[1]  # nombre del archivo
            self.type = 'obj' if path[-1] == 'j' else 'off'  # tipo de archivo
            vertices, faces = self.off2mesh() if self.type == 'off' else self.obj2mesh()
        interval = interval.split('x')
        for i in range(0, 3):
            vertices[:, i] *= float(interval[i])
        if rotate:
            x, y, z = vertices.T
            vertices = np.array((z, y, x)).T
        self.mesh = vertices, faces
        self.units = units
        self.trimesh = None  # malla del método Trimesh
        self.volume = None  # volumen de la malla
        self.radius = None  # radio original de la gota
        self.stress_bool = False  # estrés normal, verifica si se ha calculado o no

    # convierte el archivo .off de la ruta en una malla triangular
    # de la superficie que definen los puntos con las caras
    def off2mesh(self):
        """
        Converts an .off file into a triangular surface mesh.

        Returns
        -------
        tuple
            Vertices and faces of the converted mesh.
        """

        with open(self.path, 'r') as path:
            lines = path.readlines()
            if lines[0].strip() != 'OFF':
                raise ValueError('The file have not OFF format.')
            n, m, _ = list(map(int, lines[1].split()))
            vertices = np.zeros((n, 3))
            faces = np.zeros((m, 3), int)
            for i in range(2, n + 2):
                line = uf.str2num(lines[i].strip())
                vertices[i - 2] = [float(line[0]), float(line[1]), float(line[2])]
            for i in range(n + 2, n + m + 2):
                line = uf.str2num(lines[i].strip(), 'int')
                if len(line) == 4:
                    faces[i - n - 2] = [int(line[1]), int(line[2]), int(line[3])]
                else:
                    if all(line) != -1:
                        faces[i - n - 2] = [int(line[0]), int(line[1]), int(line[2])]
            return vertices, faces

    # convierte el archivo .obj de la ruta en una malla triangular
    # de la superficie que definen los puntos con las caras
    def obj2mesh(self):
        """
        Converts an .obj file into a triangular surface mesh.

        Returns
        -------
        tuple
            Vertices and faces of the converted mesh.
        """

        vertices, faces = [], []
        with open(self.path, 'r') as file:
            for lines in file:
                line = lines.split()
                if not line:
                    continue
                elif line[0] == 'v':
                    vertices.append([float(line[1]), float(line[2]), float(line[3])])
                elif line[0] == 'f':
                    faces.append([int(line[1]) - 1, int(line[2]) - 1, int(line[3]) - 1])
        return np.array(vertices), np.array(faces)

    # grafica la malla superficial obtenida
    # - use_fun (bool) -> utilizar una función para los colores de la cara
    # - fun (array) -> arreglo del mismo largo que las caras, indica el color de cada una de ellas
    # - cmap (str) -> mapa de colores con lo que se va a pintar, 'paper' utiliza la escala
    #   del paper "dropletsOriginal"; "bw" escala blanco y negro; otro valor utiliza el mismo cmap
    #   los nombres posibles son similares o iguales a Matplotlib
    # - title (str) -> título de la gráfica, por predeterminado es la ruta
    # - edges (bool) -> mostrar los bordes de la malla
    # - view_method (bool) -> mostrar método en el título
    def plot(self, fun=None,
             use_fun=False,
             cmap='default',
             title='path',
             view_method=False,
             edges=False):
        """
        (DEPRECATED) Plots the mesh with customizable properties for color and edges.

        Parameters
        ----------
        fun : array, optional
            Array defining the color of each face.
        use_fun : bool, optional
            Whether to use the color function.
        cmap : str, optional
            Color map.
        title : str, optional
            Plot title.
        view_method : bool, optional
            Whether to show the method in the title.
        edges : bool, optional
            Whether to show mesh edges.
        """

        fun = self.colorFun if fun is None else fun
        title = self.name if title == 'path' else title
        title += ' - ' + self.method if view_method is True else ''
        vertices, faces = self.mesh
        x, y, z = vertices[:, 0], vertices[:, 1], vertices[:, 2]
        if type(cmap) is str:
            if cmap == 'default':
                cmap = uf.colormap(self.colorFun) if self.colorFun is not None else 'paper'
            if cmap == 'paper':
                cmap = [(0.10, 0.43, 0.86),
                        (0.63, 0.78, 1.00),
                        (1.00, 1.00, 1.00),
                        (1.00, 1.00, 0.00),
                        (1.00, 0.00, 0.00)]
            if cmap == 'bw':
                cmap = [(1, 1, 1),
                        (0, 0, 0)]
            if cmap == 'red':
                cmap = [(1, 0, 0),
                        (1, 0, 0)]
        else:
            cmap = cmap
        if use_fun is False:
            fig = ff.create_trisurf(x, y, z,
                                    simplices=faces,
                                    colormap=cmap,
                                    title=title,
                                    plot_edges=edges
                                    )
        else:
            fig = ff.create_trisurf(x, y, z,
                                    simplices=faces,
                                    colormap=cmap,
                                    color_func=fun,
                                    title=title,
                                    plot_edges=edges
                                    )
        fig.show()

    def plot_go(self,
                curvatures=True,
                edges=False,
                max_value=None,
                debug_edges=False,
                color='red',
                color_edges='rgb(70,70,70)',
                axis=True,
                figure=None,
                show=True,
                opacity=1,
                layout=False):
        """
        Plots the mesh with Plotly and displays curvatures or color properties.

        Parameters
        ----------
        curvatures : bool, optional
            Whether to display curvatures.
        edges : bool, optional
            Whether to show mesh edges.
        max_value : float, optional
            Maximum value for the color scale.
        debug_edges: bool, optional
            Debugging the edges
        color: string, optional
            Color of plotted object
        color_edges: string, optional
            Color of edges
        axis: bool, optional
            Whether to show axis and grid.
        figure: plotly.graph_objects.Figure, optional
            The Plotly figure object to be displayed
        show: bool, optional
            Whether to display the figure
        opacity: float, optional
            Value of opacity
        layout: bool, optional
            Whether to display the layout

        Returns
        -------
        plotly.graph_objects.Figure
            Plotly figure with object
        """
        vertices, faces = self.mesh
        x, y, z = vertices[:, 2], vertices[:, 1], vertices[:, 0]
        i, j, k = faces[:, 0], faces[:, 1], faces[:, 2]
        if curvatures is True:
            if max_value is None:
                max_value = np.max(np.abs(self.meanCurvature))
            title_colorbar = 'Mean<br>Curvatures' if self.stress_bool is False else 'Normal<br>Stress<br>[nN/μm²]'
            title = self.name + ' - ' + self.method
            if figure is None:
                fig = go.Figure(layout={'title': title})
            else:
                fig = figure
            colors = uf.get_colorscale(self.meanCurvature, a=max_value)
            fig.add_mesh3d(x=x, y=y, z=z, i=i, j=j, k=k,
                           intensity=self.meanCurvature,
                           colorbar={'title': title_colorbar},
                           text=np.ndarray.astype(self.meanCurvature, str),
                           colorscale=colors,
                           opacity=opacity)
            if layout:
                def change_colors(colors, value):
                    new_colors = colors.copy()
                    for i in range(0, len(colors)):
                        if self.meanCurvature[i] <= value:
                            new_colors[i][1] = '#808080'
                        else:
                            new_colors[i] = colors[i]
                    return new_colors

                # Agregar sliders al layout de la figura
                k = len(self.meanCurvature) // 20
                sorted_curves = np.sort(self.meanCurvature)
                fig.update_layout(
                    sliders=[{'active': 0,
                              'currentvalue': {'prefix': 'Value: '},
                              'pad': {'t': 50},
                              'steps': [{'method': 'update', 'label': f'{sorted_curves[i]:.4f}',
                                         'args': [{'colorscale': [change_colors(colors, sorted_curves[i])]}]}
                                        for i in range(0, len(self.meanCurvature), k)]
                              }]
                )
        else:
            if figure is None:
                fig = go.Figure(layout={'title': self.name})
            else:
                fig = figure
            fig.add_mesh3d(x=x, y=y, z=z, i=i, j=j, k=k, color=color, opacity=opacity)
        if edges:
            tri_points = vertices[faces]
            xe, ye, ze = [], [], []
            for t in tri_points:
                xe.extend([t[k % 3][0] for k in range(4)] + [None])
                ye.extend([t[k % 3][1] for k in range(4)] + [None])
                ze.extend([t[k % 3][2] for k in range(4)] + [None])
            if debug_edges:
                (xe, ze) = (ze, xe)
            fig.add_scatter3d(x=ze, y=ye, z=xe, mode='lines', name='',
                              line=dict(color=color_edges, width=1))
        if not axis:
            fig.update_layout(
                scene=dict(
                    xaxis=dict(visible=False),  # Oculta el eje x
                    yaxis=dict(visible=False),  # Oculta el eje y
                    zaxis=dict(visible=False)  # Oculta el eje z
                )
            )

        if show:
            fig.show()
        return fig

    def normals_per_face(self):
        """
        Calculates normals for each face of the mesh.

        Returns
        -------
        array
            Normals for each face of the mesh.
        """

        vertices, faces = self.mesh
        normals = np.zeros((len(faces), 3))
        for i in range(0, len(faces)):
            a, b, c = faces[i]
            va, vb, vc = vertices[a], vertices[b], vertices[c]
            e1 = va - vb
            e2 = va - vc
            n = np.cross(e1, e2)
            n /= np.linalg.norm(n)
            normals[i] = n
        self.normalsFaces = normals
        return normals

    # busca los vecinos de cada cara. Para ello,
    # se agruparán los bordes, asumiendo que cada borde tiene
    # solamente dos caras adyacentes, vale decir
    # es una malla "regular". Retorna un dict con lo anterior
    def search_neighbors(self):
        """
        Finds neighbors for each face of the mesh, assuming each edge has two adjacent faces.

        Returns
        -------
        dict
            Dictionary of neighbors for each face of the mesh.
        """

        neighbors = {}
        _, faces = self.mesh
        for i in range(0, len(faces)):
            a, b, c = faces[i]
            # primer borde
            e1 = [(a, b), (b, a)]
            r1, r2 = e1
            if r1 not in neighbors and r2 not in neighbors:
                neighbors[r1] = [i]
            elif r2 in neighbors and r1 not in neighbors:
                neighbors[r2] += [i]
            elif r2 not in neighbors and r1 in neighbors:
                neighbors[r1] += [i]

            # segundo borde
            e2 = [(a, c), (c, a)]
            r1, r2 = e2
            if r1 not in neighbors and r2 not in neighbors:
                neighbors[r1] = [i]
            elif r2 in neighbors and r1 not in neighbors:
                neighbors[r2] += [i]
            elif r2 not in neighbors and r1 in neighbors:
                neighbors[r1] += [i]

            # tercer borde
            e3 = [(b, c), (c, b)]
            r1, r2 = e3
            if r1 not in neighbors and r2 not in neighbors:
                neighbors[r1] = [i]
            elif r2 in neighbors and r1 not in neighbors:
                neighbors[r2] += [i]
            elif r2 not in neighbors and r1 in neighbors:
                neighbors[r1] += [i]
        self.neighborsFaces = neighbors
        return neighbors

    # retorna los puntos conectados al vértice i
    # es un diccionario, pues no todos tienen la misma
    # cantidad de puntos conectados. También, retorna los
    # vértices conectados agrupados en triángulos.
    def neighbors_points(self):
        """
        Finds connected points for each vertex and groups vertices into triangles.

        Returns
        -------
        tuple
            Dictionary of connected points and dictionary of connected triangles.
        """

        _, faces = self.mesh
        connected = {}
        triangles = {}
        for i in range(0, len(faces)):
            a, b, c = faces[i]
            # primero, obtención de los puntos conectados
            # análisis del primer vértice
            if a not in connected:
                connected[a] = [b, c]
            else:
                if b not in connected[a]:
                    connected[a] += [b]
                if c not in connected[a]:
                    connected[a] += [c]
            # análisis del segundo vértice
            if b not in connected:
                connected[b] = [a, c]
            else:
                if a not in connected[b]:
                    connected[b] += [a]
                if c not in connected[b]:
                    connected[b] += [c]
            # análisis del tercer vértice
            if c not in connected:
                connected[c] = [a, b]
            else:
                if b not in connected[c]:
                    connected[c] += [b]
                if a not in connected[c]:
                    connected[c] += [a]
            # luego, obtención de los puntos asociados en cada triángulo
            # primero, para el vértice a
            if a not in triangles:
                triangles[a] = [i]
            else:
                triangles[a] += [i]
            # luego, para el vértice b
            if b not in triangles:
                triangles[b] = [i]
            else:
                triangles[b] += [i]
            # y por último, el vértice c
            if c not in triangles:
                triangles[c] = [i]
            else:
                triangles[c] += [i]
        self.neighborsPoints = connected
        self.trianglePoints = triangles
        return connected, triangles

    # le otorga un color a cada cara
    # dependiendo de los valores de fun. También,
    # permite colorear según la curvatura gaussiana
    # dependiendo del valor booleano de gauss
    def colormap(self, fun=None, gauss=False, normalize=False, max_value=None, type_plot='mean'):
        """
        (DEPRECATED) Colors each face of the mesh according to curvature values or other properties.

        Parameters
        ----------
        fun : array, optional
            Values to use for coloring the mesh.
        gauss : bool, optional
            Whether to use Gaussian curvature for coloring.
        normalize : bool, optional
            Whether to normalize the colors.
        max_value : float, optional
            Maximum value for normalization.
        type_plot : str, optional
            Type of coloring ('mean', 'max', 'min', etc.).
        """

        fun = self.meanCurvature if fun is None else fun
        fun = self.gaussCurvature if gauss is True else fun
        vertices, faces = self.mesh
        colors = np.zeros(len(faces))
        for i in range(0, len(faces)):
            a, b, c = faces[i]
            vect = uf.anti_nan([fun[a], fun[b], fun[c]])
            # el máximo entre cada vértice
            if type_plot == 'max':
                colors[i] = max(vect)
            # el mínimo entre cada vértice
            if type_plot == 'min':
                colors[i] = min(vect)
            # el promedio de los vértices
            if type_plot == 'mean':
                colors[i] = np.mean(vect)
            # solo las caras con todos sus vertices no negativos
            if type_plot == 'one':
                if all(vect) >= 0:
                    colors[i] = 1
                if all(vect) <= 0:
                    colors[i] = -1
                if all(vect) == 0:
                    colors[i] = 0
        if normalize is True:
            max_value = np.max(colors) if max_value is None else max_value
            colors /= max_value
        self.colorFun = colors
        return colors

    # arregla la malla proporcionada. Es decir, repara
    # las caras que tienen más de 3 vecinos (y otros defectos) utilizando
    # el software de Marco Attene. Para más información ver:
    # https://pymeshfix.pyvista.org/
    def fix_mesh(self):
        """
        Repairs the mesh using pymeshfix software.

        Returns
        -------
        tuple
            Vertices and faces of the repaired mesh.
        """

        vertices, faces = self.mesh
        meshfix = pymeshfix.MeshFix(vertices, faces)
        meshfix.repair()
        self.mesh = meshfix.v, meshfix.f
        return self.mesh

    # forma discreta de calcular la curvatura basada
    # en el algoritmo visto en "Keenan Crane’s lecture" disponible
    # en https://youtu.be/sokeN5VxBB8 y https://youtu.be/NlU1m-OfumE. El primer
    # link es para la fórmula de la curvatura. El segundo para el ángulo.
    # Desde el minuto 3:16 y 4:20 respectivamente.
    def discrete_mean_curvature(self):
        """
        Calculates discrete mean curvature based on Keenan Crane's method.

        Returns
        -------
        array
            Mean curvature for each point of the mesh.
        """

        time0 = time()
        if self.method == 'Discrete Method':
            return self.meanCurvature
        if self.neighborsFaces is None:
            self.search_neighbors()
        if self.normalsFaces is None:
            self.normals_per_face()
        if self.neighborsPoints is None:
            self.neighbors_points()
        neighbors = self.neighborsFaces
        normals = self.normalsFaces
        connected = self.neighborsPoints
        vertices, faces = self.mesh
        mean_curvature = np.zeros(len(vertices))
        for i in range(0, len(vertices)):
            points = connected[i]
            v1 = vertices[i]
            hi = 0
            for p in points:
                v2 = vertices[p]
                length = np.linalg.norm(v1 - v2)
                edge = (v1 - v2) / length
                try:
                    f1, f2 = neighbors[i, p]
                except KeyError:
                    f1, f2 = neighbors[p, i]
                n1, n2 = normals[f1], normals[f2]
                phi = np.arctan2(np.dot(edge, np.cross(n1, n2)),
                                 np.dot(n1, n2))
                hi += phi * length
            mean_curvature[i] = hi / 4
        self.meanCurvature = mean_curvature
        self.method = 'Discrete Method'
        print(f'execution time ({self.method}) = {(time() - time0):.4f} [s]')
        return mean_curvature


    def laplacian_curvature(self):
        """
        Calculates mean curvature using the Laplacian method.
        http://rodolphe-vaillant.fr/entry/33/curvature-of-a-triangle-mesh-definition-and-computation

        Returns
        -------
        array
            Mean curvature for each point of the mesh.
        """

        time0 = time()
        if self.method == 'Laplacian Method':
            return self.meanCurvature
        if self.neighborsFaces is None:
            self.search_neighbors()
        if self.normalsPoints is None:
            self.normals_per_point()
        if self.neighborsPoints is None:
            self.neighbors_points()
        neighbors = self.neighborsFaces
        points = self.neighborsPoints
        normals = self.normalsPoints
        triangles = self.trianglePoints
        vertices, faces = self.mesh
        angles = uf.get_angles(vertices, faces)  # ángulos para cada cara
        # obtengamos el laplaciano para cada punto, a través
        # del método de las cotangentes, o su expresión discreta.
        areas = np.zeros(len(vertices))  # áreas
        alpha = np.zeros((len(vertices), len(vertices)))  # ángulo 1
        beta = np.zeros((len(vertices), len(vertices)))  # ángulo 2
        for i in range(0, len(vertices)):
            areas[i] = uf.mixed_area(i, vertices, faces, triangles, angles)  # función para áreas
            connected = points[i]
            for c in connected:
                alpha[i, c], beta[i, c] = uf.angles_triangles(i, c, neighbors, angles, faces)  # función para ángulos
        self.areas = areas
        laplace = np.zeros((len(vertices), 3))
        for i in range(0, len(vertices)):
            laplace[i] = (1 / (2 * areas[i])) \
                         * sum((uf.cot(alpha[i, j]) + uf.cot(beta[i, j]))
                               * (vertices[j] - vertices[i]) for j in points[i])
        # sigue que, necesitamos determinar el signo de H, para esto
        # se realizará el producto punto entre la normal en el punto
        # y su laplaciano respectivo negativo
        sign = np.zeros(len(vertices))
        for i in range(0, len(vertices)):
            dot = np.dot(normals[i], -laplace[i])
            if dot >= 0:
                sign[i] = 1
            else:
                sign[i] = -1

        # por último, la curvatura promedio
        mean_curvature = np.zeros(len(vertices))
        for i in range(0, len(vertices)):
            mean_curvature[i] = sign[i] * np.linalg.norm(laplace[i]) / 2
        self.meanCurvature = mean_curvature
        self.method = 'Laplacian Method'
        print(f'execution time ({self.method}) = {(time() - time0):.4f} [s]')
        return mean_curvature

    # otra forma de calcular la curvatura promedio, basado en el algoritmo propuesto en
    # http://rodolphe-vaillant.fr/entry/33/curvature-of-a-triangle-mesh-definition-and-computation
    def optimized_laplacian(self):
        """
        Calculates mean curvature using an optimized Laplacian method.

        Returns
        -------
        array
            Mean curvature for each point of the mesh.
        """

        time0 = time()
        if self.method == 'Optimized Laplacian Method':
            return self.meanCurvature
        neighbors = self.search_neighbors() if self.neighborsFaces is None else self.neighborsFaces
        points = self.neighbors_points()[0] if self.neighborsPoints is None else self.neighborsPoints
        normals = self.normals_per_point() if self.normalsPoints is None else self.normalsPoints
        triangles = self.trianglePoints
        vertices, faces = self.mesh
        angles = uf.get_angles(vertices, faces)  # ángulos para cada cara
        mean_curvature = np.zeros(len(vertices))  # curvaturas promedio
        # obtengamos el laplaciano para cada punto, a través
        # del método de las cotangentes, o su expresión discreta.
        for i in range(0, len(vertices)):
            area = uf.mixed_area(i, vertices, faces, triangles, angles)
            sum_laplacian = 0
            for j in points[i]:
                alph, bet = uf.angles_triangles(i, j, neighbors, angles, faces)
                sum_laplacian += (uf.cot(alph) + uf.cot(bet)) * (vertices[j] - vertices[i])
            laplace = (1 / (2 * area)) * sum_laplacian
            dot = np.dot(normals[i], -laplace)
            # sigue que, necesitamos determinar el signo de H, para esto
            # se realizará el producto punto entre la normal en el punto
            # y su laplaciano respectivo negativo
            if dot < 0:
                mean_curvature[i] = np.linalg.norm(laplace) / 2
            else:
                mean_curvature[i] = -np.linalg.norm(laplace) / 2
        self.meanCurvature = mean_curvature
        self.method = 'Optimized Laplacian Method'
        print(f'execution time ({self.method}) = {(time() - time0):.4f} [s]')
        return mean_curvature

    # retorna la curvatura gaussiana según lo visto en
    # http://rodolphe-vaillant.fr/entry/33/curvature-of-a-triangle-mesh-definition-and-computation
    def gaussian_curvature(self, method=False):
        """
        Calculates Gaussian curvature.
        http://rodolphe-vaillant.fr/entry/33/curvature-of-a-triangle-mesh-definition-and-computation

        Parameters
        ----------
        method : bool, optional
            Whether to update the calculation method.

        Returns
        -------
        array
            Gaussian curvature for each point of the mesh.
        """

        time0 = time()
        if self.gaussCurvature is not None:
            if method is True:
                self.method = 'Gaussian Method'
            return self.gaussCurvature
        if self.neighborsPoints is None:
            self.neighbors_points()
        triangles = self.trianglePoints
        vertices, faces = self.mesh
        gauss = np.zeros(len(vertices))
        for i in range(0, len(vertices)):
            th, ar = uf.angle_area(i, vertices, triangles)
            gauss[i] = (2 * np.pi - th) / ar
        self.gaussCurvature = gauss
        if method is True:
            self.method = 'Gaussian Method'
        print(f'execution time ({self.method}) = {(time() - time0):.4f} [s]')
        return gauss

    # obtención de las curvaturas principales
    # a través del método descrito por Taubin, 1995
    def taubin_method(self):
        """
        Calculates principal curvatures using Taubin's method.

        Returns
        -------
        array
            Mean curvature for each point of the mesh.
        """

        time0 = time()
        if self.normalsPoints is None:
            self.normals_per_point()
        if self.neighborsPoints is None:
            self.neighbors_points()
        if self.method == 'Taubin Method':
            return self.meanCurvature
        vertices, faces = self.mesh
        normals = self.normalsPoints
        neighbors = self.neighborsPoints
        k1, k2 = np.zeros(len(vertices)), np.zeros(len(vertices))
        for i in range(0, len(vertices)):
            m = np.zeros((3, 3))
            for j in range(0, len(neighbors[i])):
                p = vertices[i] - vertices[neighbors[i][j]]
                k = 2 * np.dot(normals[i], p)
                k /= np.linalg.norm(p) ** 2
                m += k * np.dot(np.transpose(p), p)
            eigen = np.real(np.linalg.eigvals(m))
            eigen[np.argmin(np.abs(eigen))] = 0
            k1[i], k2[i] = np.min(eigen), np.max(eigen)
        self.meanCurvature = (k1 + k2) / 2
        self.principals = k1, k2
        self.method = 'Taubin Method'
        print(f'execution time ({self.method}) = {(time() - time0):.4f} [s]')
        return self.meanCurvature

    def normals_per_point(self):
        """
        Calculates average normals at each point of the mesh.

        Returns
        -------
        array
            Average normals for each point of the mesh.
        """

        if self.normalsFaces is None:
            self.normals_per_face()
        normals = self.normalsFaces
        vertices, faces = self.mesh
        self.normalsPoints = trimesh.geometry.mean_vertex_normals(len(vertices), faces, normals)
        return self.normalsPoints

    def normals_rusinkiewicz(self):
        """
        Calculates normals at each point using Rusinkiewicz's method.

        Returns
        -------
        array
            Normals for each point of the mesh.
        """

        vertices, faces = self.mesh
        normals = np.zeros((len(vertices), 3))
        for i in range(0, len(faces)):
            p0 = vertices[faces[i, 0]]
            p1 = vertices[faces[i, 1]]
            p2 = vertices[faces[i, 2]]
            a, b, c = p0 - p1, p1 - p2, p2 - p0

            def len2(x):
                return x[0] ** 2 + x[1] ** 2 + x[2] ** 2

            l2a, l2b, l2c = len2(a), len2(b), len2(c)
            if any([l2a, l2b, l2c]) == 0:
                continue
            facenormal = np.cross(a, b)
            normals[faces[i, 0]] += facenormal * (1 / (l2a * l2c))
            normals[faces[i, 1]] += facenormal * (1 / (l2b * l2a))
            normals[faces[i, 2]] += facenormal * (1 / (l2c * l2b))
        for i in range(0, len(vertices)):
            normals[i] /= np.linalg.norm(normals[i])
        self.normalsPoints = normals
        return normals

    def rusinkiewicz_curvature(self, force_normal=True, ldl_numpy=False, invert=True):
        """
        Calculates mean curvature using Rusinkiewicz's method.
        https://gfx.cs.princeton.edu/pubs/Rusinkiewicz_2004_ECA/index.php

        Parameters
        ----------
        force_normal : bool, optional
            Whether to force normals.
        ldl_numpy : bool, optional
            Whether to use LDL decomposition in numpy.
        invert : bool, optional
            Whether to invert curvatures sign.

        Returns
        -------
        array
            Mean curvature for each point of the mesh.
        """

        time0 = time()
        if self.method == 'Rusinkiewicz Method':
            return self.meanCurvature
        if self.normalsPoints is None or force_normal is True:
            self.normals_rusinkiewicz()
        vertices, faces = self.mesh
        normals = self.normalsPoints
        pdir1 = np.zeros((len(vertices), 3))
        pdir2 = np.zeros((len(vertices), 3))
        curve1, curve2, curve12 = np.zeros(len(vertices)), np.zeros(len(vertices)), np.zeros(len(vertices))

        # coordenadas iniciales:
        for i in range(0, len(faces)):
            face = faces[i]
            pdir1[face[0]] = vertices[face[1]] - vertices[face[0]]  # borde del triángulo
            pdir1[face[1]] = vertices[face[2]] - vertices[face[1]]  # ''
            pdir1[face[2]] = vertices[face[0]] - vertices[face[2]]  # ''
        for i in range(0, len(vertices)):
            pdir1[i] = np.cross(pdir1[i], normals[i])  # p x n
            pdir1[i] /= np.linalg.norm(pdir1[i])  # normalización
            pdir2[i] = np.cross(normals[i], pdir1[i])  # n x p

        # áreas
        point_areas, corner_areas = uf.get_point_area(vertices, faces)

        # cálculo de las curvaturas en cada cara:
        for i in range(0, len(faces)):
            edges = np.array([vertices[faces[i, 2]] - vertices[faces[i, 1]],
                              vertices[faces[i, 0]] - vertices[faces[i, 2]],
                              vertices[faces[i, 1]] - vertices[faces[i, 0]]])
            # N-T-B por cara junto a la matriz de Weingarten w
            t = edges[0] / np.linalg.norm(edges[0])
            n = np.cross(edges[0], edges[1])
            b = np.cross(n, t)
            b /= np.linalg.norm(b)  # normalización
            m, w = np.zeros(3), np.zeros((3, 3))
            for j in range(0, 3):
                u, v = np.dot(edges[j], t), np.dot(edges[j], b)
                w[0, 0] += u * u
                w[0, 1] += u * v
                w[2, 2] += v * v
                dn = normals[faces[i, (j - 1) % 3]] - normals[faces[i, (j + 1) % 3]]
                dnu = np.dot(dn, t)
                dnv = np.dot(dn, b)
                m[0] += dnu * u
                m[1] += dnu * v + dnv * u
                m[2] += dnv * v
            w[1, 1], w[1, 2] = w[0, 0] + w[2, 2], w[0, 1]
            if ldl_numpy is False:
                diag = np.zeros(3)
                value = uf.ldltdc(w, diag)
                if value is False:
                    continue
                m = uf.ldltsl(w, diag, m, same=True)
            else:
                m = uf.ldl_solve(w, m)
            for j in range(0, 3):
                vj = faces[i, j]
                c1, c12, c2 = uf.proj_curves(t, b, m[0], m[1], m[2], pdir1[vj], pdir2[vj])
                wt = corner_areas[i, j] / point_areas[vj]
                curve1[vj] += wt * c1
                curve12[vj] += wt * c12
                curve2[vj] += wt * c2

        # diagonalización para obtener las direcciones y curvaturas,
        # aunque solo nos interesan las curvaturas
        for i in range(0, len(vertices)):
            curve1[i], curve2[i], pdir1[i], pdir2[i] = uf.diagonalize_curves(pdir1[i], pdir2[i], curve1[i],
                                                                             curve12[i], curve2[i], normals[i])
        self.meanCurvature = (curve1 + curve2) / 2
        if invert:
            self.meanCurvature *= -1
        self.method = 'Rusinkiewicz Method'
        self.principals = curve1, curve2
        print(f'execution time ({self.method}) = {(time() - time0):.4f} [s]')
        return self.meanCurvature

    def trimesh_method(self, radius=1):
        """
        Calculates mean curvature using the Trimesh method.

        Parameters
        ----------
        radius : float, optional
            Radius to consider in the calculation.

        Returns
        -------
        array
            Mean curvature for each point of the mesh.
        """

        time0 = time()
        vertices, faces = self.mesh
        if self.method == 'Trimesh Method':
            return self.meanCurvature
        if self.trimesh is None:
            self.trimesh = trimesh.Trimesh(vertices=vertices, faces=faces)
        self.meanCurvature = trimesh.curvature.discrete_mean_curvature_measure(self.trimesh, vertices, radius)
        self.method = 'Trimesh Method'
        print(f'execution time ({self.method}) = {(time() - time0):.4f} [s]')
        return self.meanCurvature

    def get_radius(self):
        """
        Calculates the radius of the mesh based on its volume, assuming a spherical shape.

        Returns
        -------
        float
            Radius of the mesh.
        """

        if self.volume is None:
            self.get_volume()
        self.radius = (3 * self.volume / (np.pi * 4)) ** (1 / 3)
        return self.radius

    def get_volume(self):
        """
        Calculates the volume of the mesh by summing tetrahedral volumes.

        Returns
        -------
        float
            Volume of the mesh.
        """

        vertices, faces = self.mesh
        vol = 0
        for i in range(0, len(faces)):
            a, b, c = faces[i]
            p, q, r = vertices[a], vertices[b], vertices[c]
            vol += (1 / 6) * np.dot(p, np.cross(q, r))
        self.volume = abs(vol)
        return self.volume

    def stress(self, gamma=1):
        """
        Calculates normal forces in the mesh using curvature and a gamma parameter.

        Parameters
        ----------
        gamma : float, optional
            Scaling factor for stress calculation.

        Returns
        -------
        array
            Normal stresses for each point of the mesh.
        """

        if self.radius is None:
            self.get_radius()
        name_curvatures = self.method
        rad = np.ones(len(self.meanCurvature)) * (1 / self.radius)
        self.meanCurvature = 2 * gamma * (self.meanCurvature - rad)
        self.stress_bool = True
        self.method = 'Normal Stress with ' + name_curvatures
        return self.meanCurvature

    def cut(self, z_min=None, z_max=None):
        """
        Cuts the mesh within a z interval defined by z_min and z_max.

        Parameters
        ----------
        z_min : float, optional
            Lower limit in z for the cut.
        z_max : float, optional
            Upper limit in z for the cut.

        Returns
        -------
        tuple
            Vertices and faces of the cut mesh.
        """

        vertices, faces = self.mesh
        if z_min is None:
            z_min = np.min(vertices[:, 2])
        if z_max is None:
            z_max = np.max(vertices[:, 2])
        new_faces = []
        for i in range(0, len(faces)):
            v1, v2, v3 = vertices[faces[i]]
            if z_max >= v1[2] >= z_min and z_max >= v2[2] >= z_min and z_max >= v3[2] >= z_min:
                new_faces.append(faces[i])
        self.mesh = vertices, np.array(new_faces, np.uint16)
        return self.mesh

    def limit_curvatures(self, ho=None, hf=None):
        assert self.meanCurvature is not None, 'meanCurvature is None'
        mini, maxi = np.min(self.meanCurvature), np.max(self.meanCurvature)
        if ho is None:
            ho = mini
        if hf is None:
            hf = maxi
        if ho is None and hf is None:
            return self.meanCurvature
        if ho < mini:
            ho = mini
        if hf > maxi:
            hf = maxi
        self.meanCurvature = np.clip(self.meanCurvature, ho, hf)
        return self.meanCurvature


if __name__ == '__main__':
    for i in range(20, 25):
        f = Forces3D(path=r"C:\RSI\stack_ROI_02_AC3D_Obj_\obj (" + str(i) + ").off",
                     interval='0.163x0.163x0.5', rotate=True)
        f.fix_mesh()
        f.rusinkiewicz_curvature()
        f.stress(gamma=1)
        f.limit_curvatures(-2, 2)
        f.plot_go(axis=False, max_value=2)

