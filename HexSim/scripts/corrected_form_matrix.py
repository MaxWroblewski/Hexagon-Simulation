import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import griddata

class Hexagonal_Prism:
    def __init__(self, a: float, rho: object, l: float, rho_0: float, n_theta: int, n_r: int, n_z: int):
        # Hexagon parameters
        self.a = a
        self.l = l
        self.rho_0 = rho_0
        self.rho = rho

        # Resolutions
        self.n_theta = n_theta
        self.n_r = n_r
        self.n_z = n_z

        self.I_xx = np.float64(0.0)
        self.I_yy = np.float64(0.0)
        self.I_zz = np.float64(0.0)
        self.I_xy = np.float64(0.0)
        self.I_xz = np.float64(0.0)
        self.I_yz = np.float64(0.0)

        # Discretize angular and z ranges
        self.theta_vals = np.linspace(0, 2 * np.pi, self.n_theta, endpoint=False)
        self.z_vals = np.linspace(-self.l / 2, self.l / 2, self.n_z)
        
        # Ensure n_z is odd for Simpson's Rule
        if self.n_z % 2 == 0:
            self.n_z += 1
            self.z_vals = np.linspace(-self.l / 2, self.l / 2, self.n_z)
        self.delta_z = self.l / (self.n_z - 1)

    def hexagonal_radius(self, theta, a):
        # Corrected hexagonal radius calculation
        theta = theta % (2 * np.pi)
        sector = int(theta / (np.pi / 3))
        angle_in_sector = theta - sector * (np.pi / 3)
        if angle_in_sector > (np.pi / 6):
            angle = np.pi / 3 - angle_in_sector
        else:
            angle = angle_in_sector
        return (a * np.cos(np.pi / 6)) / np.cos(angle)
        
    def form_matrix(self):
        self.points = []

        # Simpson's coefficients for z
        coefficients_z = np.ones(self.n_z)
        coefficients_z[0] = coefficients_z[-1] = 1
        coefficients_z[1:-1:2] = 4
        coefficients_z[2:-1:2] = 2

        for idx_z, z in enumerate(self.z_vals):
            simpson_coeff_z = coefficients_z[idx_z]

            for theta_idx, theta in enumerate(self.theta_vals):
                r_max = self.hexagonal_radius(theta, self.a)
                # Ensure n_r is odd for Simpson's Rule
                n_r = self.n_r
                if n_r % 2 == 0:
                    n_r += 1
                r_vals = np.linspace(0, r_max, n_r)
                delta_r = r_max / (n_r - 1)

                # Simpson's coefficients for r
                coefficients_r = np.ones(n_r)
                coefficients_r[0] = coefficients_r[-1] = 1
                coefficients_r[1:-1:2] = 4
                coefficients_r[2:-1:2] = 2

                for idx_r, r in enumerate(r_vals):
                    simpson_coeff_r = coefficients_r[idx_r]

                    x = r * np.cos(theta)
                    y = r * np.sin(theta)

                    # Evaluate density at (x, y)
                    density = self.rho(x, y)

                    # Volume element in cylindrical coordinates with Simpson's rule
                    volume_element = r * delta_r * self.delta_z * (simpson_coeff_r * simpson_coeff_z) / 9  # 3*3=9

                    # Store Point
                    row = {'x': x, 'y': y, 'z': z, 'density': density}
                    self.points.append(row)

                    # Contribution to diagonal components
                    self.I_xx += density * (y**2 + z**2) * volume_element
                    self.I_yy += density * (x**2 + z**2) * volume_element
                    self.I_zz += density * (x**2 + y**2) * volume_element
                    # Contribution to off-diagonal components
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

        # Check if A ≈ P @ D @ P_inv
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
        # Create a 3D scatterplot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection=projection)

        # Normalize density values for the colormap
        norm = plt.Normalize(self.points['density'].min(), self.points['density'].max())
        colors = plt.cm.viridis(norm(self.points['density']))

        # Scatter plot
        scatter = ax.scatter(
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
        colorbar = fig.colorbar(scatter, ax=ax, pad=0.1)
        colorbar.set_label('Density')

        # Add labels
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Axis')
        ax.set_zlabel('Z Axis')

        ax.set_box_aspect([1, 1, 1])
        ax.set_xlim([-2 * self.a, 2 * self.a])
        ax.set_ylim([-2 * self.a, 2 * self.a])
        ax.set_zlim([-self.l, self.l])

        # Show the plot
        plt.show()

    def run_eval(self, w0):
        T = (2 * np.pi / w0) * np.sqrt(self.I_xx * self.I_zz / (np.abs((self.I_yy - self.I_xx) * (self.I_zz - self.I_yy))))
        return T, 1 / T

