import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D

class Hexagonal_Prism:
    
    def __init__(self, a:float, rho:object ,l:float, rho_0:float, n_x:int, n_y:int, n_z:int):
        #Def hexagon parameters
        self.a = a
        self.l = l
        self.rho_0 = rho_0
        self.rho = rho

        #Def resolutions
        self.n_x = n_x
        self.n_y = n_y
        self.n_z = n_z

        self.I_xx = 0.0
        self.I_yy = 0.0
        self.I_zz = 0.0
        self.I_xy = 0.0
        self.I_xz = 0.0
        self.I_yz = 0.0

        # Discretize z range
        self.z_vals = np.linspace(-self.l / 2, self.l / 2, self.n_z)

        self.max_x = self.a
        self.min_x = -self.a
        self.max_y = np.sqrt(3)*self.a/2
        self.min_y = -np.sqrt(3)*self.a/2

        self.delta_x = (self.max_x - self.min_x) / self.n_x
        self.delta_y = (self.max_y - self.min_y) / self.n_y
        self.delta_z = self.l / self.n_z

        self.volume_element = self.delta_x * self.delta_y * self.delta_z

        self.x_vals = np.linspace(self.min_x, self.max_x, self.n_x)
        self.y_vals = np.linspace(self.min_y, self.max_y, self.n_y)

    def is_inside_hexagon(self, x, y):
        # Check if point is within the horizontal limits
        if abs(x) > self.a:
            return False
        # Check if point is within the vertical limits
        if abs(y) > (np.sqrt(3) / 2) * self.a:
            return False
        # Check if point is within the slanted sides
        if abs(y) > np.sqrt(3) * (self.a - abs(x)):
            return False
        return True

    def total_mass(self):
        total_mass = self.points['density'].sum() * self.delta_x * self.delta_y * self.delta_z
        return total_mass
    
    #def form_matrix(self):
        
    #    self.points = []
    #   
    #    for z in self.z_vals:
    #            for x in self.x_vals:
    #                for y in self.y_vals:
    #                    if self.is_inside_hexagon(x, y):
    #                        
    #                        # Evaluate density at (x, y)
    #                        density = self.rho(x, y)

    #                        # Store Point
    #                        row = {'x':x, 'y':y, 'z':z, 'density': density}
    #                        self.points.append(row)


    #                        # Contribution to diagonal components
    #                        self.I_xx += density * (y**2 + z**2) * self.volume_element
    #                        self.I_yy += density * (x**2 + z**2) * self.volume_element
    #                        self.I_zz += density * (x**2 + y**2) * self.volume_element
    #                        # Contribution to off-diagonal components
    #                        self.I_xy += -density * (x * y) * self.volume_element
    #                        self.I_xz += -density * (x * z) * self.volume_element
    #                        self.I_yz += -density * (y * z) * self.volume_element




    def form_matrix(self):
        self.points = []
        
        # Ensure n_x, n_y, n_z are even for Simpson's Rule
        if self.n_x % 2 != 0:
            self.n_x += 1
        if self.n_y % 2 != 0:
            self.n_y += 1
        if self.n_z % 2 != 0:
            self.n_z += 1

        # Update grids and deltas
        self.x_vals = np.linspace(self.min_x, self.max_x, self.n_x)
        self.y_vals = np.linspace(self.min_y, self.max_y, self.n_y)
        self.z_vals = np.linspace(-self.l / 2, self.l / 2, self.n_z)
        self.delta_x = (self.max_x - self.min_x) / (self.n_x - 1)
        self.delta_y = (self.max_y - self.min_y) / (self.n_y - 1)
        self.delta_z = self.l / (self.n_z - 1)
        
        # Simpson's coefficients
        coefficients_x = np.ones(self.n_x)
        coefficients_x[1:-1:2] = 4
        coefficients_x[2:-1:2] = 2

        coefficients_y = np.ones(self.n_y)
        coefficients_y[1:-1:2] = 4
        coefficients_y[2:-1:2] = 2

        coefficients_z = np.ones(self.n_z)
        coefficients_z[1:-1:2] = 4
        coefficients_z[2:-1:2] = 2

        for idx_z, z in enumerate(self.z_vals):
            simpson_coeff_z = coefficients_z[idx_z]
            for idx_x, x in enumerate(self.x_vals):
                simpson_coeff_x = coefficients_x[idx_x]
                for idx_y, y in enumerate(self.y_vals):
                    simpson_coeff_y = coefficients_y[idx_y]
                    if self.is_inside_hexagon(x, y):
                        density = self.rho(x, y)
                        simpson_factor = (simpson_coeff_x * simpson_coeff_y * simpson_coeff_z) / 27
                        volume_element = self.delta_x * self.delta_y * self.delta_z * simpson_factor

                        row = {'x':x, 'y':y, 'z':z, 'density': density}
                        self.points.append(row)
                        # Update moments of inertia
                        self.I_xx += density * (y**2 + z**2) * volume_element
                        self.I_yy += density * (x**2 + z**2) * volume_element
                        self.I_zz += density * (x**2 + y**2) * volume_element
                        self.I_xy += -density * x * y * volume_element
                        self.I_xz += -density * x * z * volume_element
                        self.I_yz += -density * y * z * volume_element


        self.points = pd.DataFrame(self.points)

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
        I_reconstructed = P @ self.I @ P_inv

        return self.I, P, inertia_tensor_non_uniform, I_reconstructed

    
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

    def hexagon_scatter(self, s=0.2, alpha=1, projection='3d'):
        import matplotlib.pyplot as plt
        from matplotlib import cm

        # Create a 3D scatterplot
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection=projection)

        # Normalize density values for the colormap
        norm = plt.Normalize(self.points['density'].min(), self.points['density'].max())
        colors = cm.viridis(norm(self.points['density']))

        # Scatter plot
        scatter = self.ax.scatter(
            self.points['x'], 
            self.points['y'], 
            self.points['z'], 
            c=self.points['density'], 
            cmap='viridis',  # Choose preferred colormap
            marker='o', 
            s=s, 
            alpha=alpha
        )

        # Add a colorbar to represent density
        colorbar = self.fig.colorbar(scatter, ax=self.ax, pad=0.1)
        colorbar.set_label('Density')

        # Add labels
        self.ax.set_xlabel('X Axis')
        self.ax.set_ylabel('Y Axis')
        self.ax.set_zlabel('Z Axis')

        self.ax.set_box_aspect([1, 1, 1])
        self.ax.set_xlim([-2 * self.a, 2 * self.a])
        self.ax.set_ylim([-2 * self.a, 2 * self.a])
        self.ax.set_zlim([-2 * self.a, 2 * self.a])

        # Show the plot
        plt.show()

    def run_eval(self, w0):
        T = (2*np.pi/w0)*np.sqrt(self.I_xx * self.I_zz/(np.abs((self.I_yy-self.I_xx)*(self.I_zz-self.I_yy))))
        return T, 1/T




                                