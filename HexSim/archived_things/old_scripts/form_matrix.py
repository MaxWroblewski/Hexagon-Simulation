import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
import coxeter

class Hexagonal_Prism:
    
    def __init__(self, a:float, rho ,l:float, rho_0:float, n_theta:int, n_r:int, n_z:int):
        #Def hexagon parameters
        self.a = a
        self.l = l
        self.rho_0 = rho_0
        self.rho = rho

        #Def resolutions
        self.n_theta = n_theta
        self.n_r = n_r
        self.n_z = n_z

        self.I_xx = 0.0
        self.I_yy = 0.0
        self.I_zz = 0.0
        self.I_xy = 0.0
        self.I_xz = 0.0
        self.I_yz = 0.0

        # Discretize angular, radial, and z ranges
        self.theta_vals = np.linspace(0, 2 * np.pi, self.n_theta, endpoint=False)
        self.z_vals = np.linspace(-self.l / 2, self.l / 2, self.n_z)

    def hexagonal_radius(self, theta, a):
        sector_angle = np.pi / 3
        theta_mod = theta % sector_angle
        return a / np.cos(theta_mod - sector_angle / 2)




    def form_matrix(self):

        self.points = []

        for z in self.z_vals:
            for theta in self.theta_vals:
                r_max = self.hexagonal_radius(theta, self.a)
                r_vals = np.linspace(0, r_max, self.n_r)
                for r in r_vals:
                    # Convert to Cartesian coordinates
                    x = r * np.cos(theta)
                    y = r * np.sin(theta)
                    
                    row = {'x':x, 'y':y, 'z':z}
                    self.points.append(row)

                    # Evaluate density at (x, y)
                    density = self.rho(x, y)
                    # Contribution to diagonal components
                    self.I_xx += density * (r**2 * np.sin(theta)**2 + z**2) * r
                    self.I_yy += density * (r**2 * np.cos(theta)**2 + z**2) * r
                    self.I_zz += density * (r**2) * r
                    # Contribution to off-diagonal components
                    self.I_xy += -density * (r**2 * np.cos(theta) * np.sin(theta)) * r
                    self.I_xz += -density * (r * np.cos(theta) * z) * r
                    self.I_yz += -density * (r * np.sin(theta) * z) * r
        
        self.points = pd.DataFrame(self.points)

        # Scale by integration step sizes
        delta_theta = 2 * np.pi / self.n_theta
        delta_r = self.a / self.n_r
        delta_z = self.l / self.n_z

        self.I_xx *= delta_r * delta_theta * delta_z
        self.I_yy *= delta_r * delta_theta * delta_z
        self.I_zz *= delta_r * delta_theta * delta_z
        self.I_xy *= delta_r * delta_theta * delta_z
        self.I_xz *= delta_r * delta_theta * delta_z
        self.I_yz *= delta_r * delta_theta * delta_z

        # Moment of inertia tensor for non-uniform density
        inertia_tensor_non_uniform = np.array([
            [self.I_xx, self.I_xy, self.I_xz],
            [self.I_xy, self.I_yy, self.I_yz],
            [self.I_xz, self.I_yz, self.I_zz]
        ])

        eigenvalues, eigenvectors = np.linalg.eigh(inertia_tensor_non_uniform)
        self.I = np.diag(eigenvalues)

        # Verify diagonalization: A = PDP^-1
        P = eigenvectors
        P_inv = np.linalg.inv(P)

        # Check if A ≈ P @ D @ P_inv
        I_reconstructed = P @ self.I @ P_inv

        return self.I, P, inertia_tensor_non_uniform, I_reconstructed
    

    def form_uniform_matrix(self):

        # Call the hexagon dimensions
        l = self.l
        a = self.a

        #Define the angles of a regular hexagon
        angles = np.linspace(0, 2*np.pi, 7)[:-1]

        #Construct vertices
        #For the bottom face at z = -l/2
        bottom_vertices = [(a * np.cos(theta), a * np.sin(theta), -l/2) for theta in angles]
        #For the top face at z = +l/2
        top_vertices = [(a * np.cos(theta), a * np.sin(theta), l/2) for theta in angles]

        #Combine the vertices
        vertices = bottom_vertices + top_vertices

        # # Define faces:
        # # Bottom face (using vertices 0 through 5)
        # bottom_face = list(range(0, 6))
        # # Top face (using vertices 6 through 11)
        # top_face = list(range(6, 12))
        # # Side faces (six quadrilaterals connecting bottom and top edges)
        # side_faces = []
        # for i in range(6):
        #     side_face = [i, (i + 1) % 6, ((i + 1) % 6) + 6, i + 6]
        #     side_faces.append(side_face)

        # # Combine all faces
        # faces = [bottom_face, top_face] + side_faces

        # Create a Polyhedron from the vertices and faces
        hexagon = coxeter.shapes.ConvexPolyhedron(vertices)

        # Compute the inertia tensor about the centroid (assuming uniform density)
        inertia_tensor = hexagon.inertia_tensor*self.rho_0

        return inertia_tensor

    def density_colormap(self):
        # Compute density for the points in the DataFrame
        self.points['density'] = self.points.apply(lambda row: self.rho(row['x'], row['y']), axis=1)

        # Create a grid for visualization
        grid_x, grid_y = np.linspace(self.points['x'].min(), self.points['x'].max(), 300), np.linspace(self.points['y'].min(), self.points['y'].max(), 300)
        grid_x, grid_y = np.meshgrid(grid_x, grid_y)

        # Interpolate density values to the grid
        grid_z = griddata((self.points['x'], self.points['y']), self.points['density'], (grid_x, grid_y), method='cubic')

        # Plot the density colormap
        plt.figure(figsize=(8, 6))
        plt.contourf(grid_x, grid_y, grid_z, levels=100, cmap='viridis')
        plt.colorbar(label='Density')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Density Colormap')
        plt.show()

    def hexagon_scatter(self, s =.2, alpha = 1, projection='3d'):
        # Create a 3D scatterplot
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection=projection)

        # Scatter plot
        self.ax.scatter(self.points['x'], self.points['y'], self.points['z'], c='b', marker='o', s = s, alpha = alpha)  # You can customize the color and marker

        # Add labels
        self.ax.set_xlabel('X Axis')
        self.ax.set_ylabel('Y Axis')
        self.ax.set_zlabel('Z Axis')

        self.ax.set_box_aspect([1, 1, 1])
        self.ax.set_xlim([-2*self.a, 2*self.a])
        self.ax.set_ylim([-2*self.a, 2*self.a])
        self.ax.set_zlim([-2*self.a, 2*self.a])
        # Show the plot
        plt.show()
    
    def run_eval(self, w0):
        T = (2*np.pi/w0)*np.sqrt(self.I_xx * self.I_zz/(np.abs((self.I_yy-self.I_xx)*(self.I_zz-self.I_yy))))
        return T, 1/T


    def plot_rho(self):
        
        rho_plot = []
        for i in range(len(self.theta_vals)):
            x = self.hexagonal_radius(self.theta_vals[i],self.a)*np.cos(self.theta_vals[i])
            y = self.hexagonal_radius(self.theta_vals[i],self.a)*np.sin(self.theta_vals[i])
            rho_plot.append(self.rho(x,y))
        
        plt.plot(rho_plot)

                                