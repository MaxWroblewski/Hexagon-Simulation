import numpy as np
import coxeter

class Hexagonal_Prism:
    
    def __init__(self, a:float, rho ,h:float):
        # Def hexagon parameters
        self.a = a
        self.h = h
        self.rho = rho
    
    def form_uniform_matrix(self):

        # Call the hexagon dimensions
        h = self.h
        a = self.a

        # Define the angles of a regular hexagon
        angles = np.linspace(0, 2*np.pi, 7)[:-1]

        # Construct vertices
        # For the bottom face at z = -l/2
        bottom_vertices = [(a * np.cos(theta), a * np.sin(theta), -h/2) for theta in angles]
        # For the top face at z = +l/2
        top_vertices = [(a * np.cos(theta), a * np.sin(theta), h/2) for theta in angles]

        # Combine the vertices
        vertices = bottom_vertices + top_vertices

        # Create a Polyhedron from the vertices and faces
        hexagon = coxeter.shapes.ConvexPolyhedron(vertices)

        # Compute the inertia tensor about the centroid (assuming uniform density)
        inertia_tensor = hexagon.inertia_tensor*self.rho

        self.I_xx = inertia_tensor[0][0]
        self.I_yy = inertia_tensor[1][1]
        self.I_zz = inertia_tensor[2][2]

        return inertia_tensor, vertices
    

    def form_matrix_triangles(self, weights = (1, 1, 1, 1, 1, 1)):

        # Call dimensions
        h = self.h
        a = self.a

        # We want to define the inertia tensor of a hexagonal prism as a composite built of six triangular prisms
        # First we need to define all of the vertices

        def triangle_vertices(a, h, theta):

            z_bottom = -h/2
            z_top = +h/2

            cos_t = np.cos(theta)
            sin_t = np.sin(theta)

            v0 = (0, 0, z_bottom)
            v1 = (a * cos_t, a * sin_t, z_bottom)
            v2 = ((a/2) * cos_t - ((a * np.sqrt(3)) / 2) * sin_t, ((a * np.sqrt(3)) / 2) * cos_t + (a/2) * sin_t, z_bottom)

            v3 = (0, 0, z_top)
            v4 = (v1[0], v1[1], z_top)
            v5 = (v2[0], v2[1], z_top)

            return [v0, v1, v2, v3, v4, v5]
        
        inertia_tensor = np.array([
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ], dtype = 'float64')
        vertices = []

        for triangle in range(6):

            v_triangle = triangle_vertices(a, h, triangle * np.pi / 3)

            triangle_i = coxeter.shapes.ConvexPolyhedron(v_triangle)

            inertia_tensor += triangle_i.inertia_tensor * self.rho * weights[triangle]
            vertices += v_triangle

        return inertia_tensor, vertices
        
    
    def run_eval(self, f0):

        T = (1/f0)*np.sqrt(self.I_xx * self.I_zz/(np.abs((self.I_yy-self.I_xx)*(self.I_zz-self.I_yy))))
        
        return T
        
                                