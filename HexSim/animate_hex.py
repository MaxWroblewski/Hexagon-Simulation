import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pandas as pd
import os

class Animate_Hex:
    def __init__(self, timestep: float, frames: int, a: float, h: float, I, df, path:str = None):

        self.timestep = timestep  # Time step for each frame
        self.frames = frames  # Number of frames in the animation
        self.a = a
        self.h = h
        self.I = I
        self.df = df
        
        # Add basis vectors (principal axes)
        _, self.principal_axes = np.linalg.eig(self.I)
        self.basis_vectors = np.array([[0, 0, 0], self.principal_axes[:,0]*a, self.principal_axes[:,1]*a, self.principal_axes[:,2]*a])
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')



        # Set up the hexagonal prism
        self.vertices = np.array(self.hexagon_vertices())
        self.faces = self.hexagon_faces()
        self.prism = Poly3DCollection(self.faces, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.25)
        self.ax.add_collection3d(self.prism)

        self.path = path

        if self.path != None:
            self.out = self.path+'_Ix_'+str(round(self.I[0][0],4))+'Iy_'+str(round(self.I[1][1],4))+'Iz_'+str(round(self.I[2][2],4))+'wx0_'+str(round(self.df['wx'][0],4))+'wy0_'+str(round(self.df['wy'][0],4))+'wz0_'+str(round(self.df['wz'][0],4))+'dt_'+str(self.df['Time'][1])
            if not os.path.exists(self.out):
                os.makedirs(self.out)
        else:
            self.out = None



    # Hexagonal prism vertices and faces definition
    def hexagon_vertices(self):
        # Bottom hexagon vertices (z = -height/2)
        angle = np.pi / 3
        bottom_hex = [
            [self.a * np.cos(i * angle), self.a * np.sin(i * angle), -self.h / 2]
            for i in range(6)
        ]
        # Top hexagon vertices (z = height/2)
        top_hex = [[x, y, self.h / 2] for x, y, _ in bottom_hex]
        return bottom_hex + top_hex

    # Hexagonal prism faces based on vertices
    def hexagon_faces(self):
        # Connect vertices to form the faces of the hexagonal prism
        faces = [
            [self.vertices[i], self.vertices[(i + 1) % 6], self.vertices[(i + 1) % 6 + 6], self.vertices[i + 6]]
            for i in range(6)
        ]
        # Add top and bottom faces
        faces.append(self.vertices[:6])  # Bottom face
        faces.append(self.vertices[6:])  # Top face
        return faces

    # Function to rotate points
    def rotate(self, points, rotation_matrix):
        return np.dot(points, rotation_matrix.T)

    # Function to get the rotation matrix from rotational velocity
    def rotation_matrix_from_velocity(self, rot_vel, dt):
        theta = np.linalg.norm(rot_vel) * dt
        if theta == 0:
            return np.eye(3)
        
        axis = rot_vel / np.linalg.norm(rot_vel)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        ux, uy, uz = axis
        R = np.array([
            [cos_theta + ux**2 * (1 - cos_theta), ux * uy * (1 - cos_theta) - uz * sin_theta, ux * uz * (1 - cos_theta) + uy * sin_theta],
            [uy * ux * (1 - cos_theta) + uz * sin_theta, cos_theta + uy**2 * (1 - cos_theta), uy * uz * (1 - cos_theta) - ux * sin_theta],
            [uz * ux * (1 - cos_theta) - uy * sin_theta, uz * uy * (1 - cos_theta) + ux * sin_theta, cos_theta + uz**2 * (1 - cos_theta)]
        ])
        return R


    
    # Animation function
    def animate(self, i):
        self.w = [self.df['wx'].iloc[i], self.df['wy'].iloc[i], self.df['wz'].iloc[i]]
        R = self.rotation_matrix_from_velocity(self.w, self.timestep)
        
        # Rotate vertices and basis vectors
        self.vertices = self.rotate(self.vertices, R)
        self.basis_vectors = self.rotate(self.basis_vectors, R)
        
        # Clear and redraw
        self.ax.cla()
        self.ax.set_box_aspect([1, 1, 1])
        self.ax.set_xlim([-2*self.a, 2*self.a])
        self.ax.set_ylim([-2*self.a, 2*self.a])
        self.ax.set_zlim([-2*self.a, 2*self.a])

        # Draw hexagonal prism
        faces = self.hexagon_faces()
        prism = Poly3DCollection(faces, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.25)
        self.ax.add_collection3d(prism)

        # Draw principal axes as arrows
        origin = np.array([0, 0, 0])
        for j in range(3):
            self.ax.quiver(*origin, *self.basis_vectors[j+1], color=['r', 'g', 'b'][j], length=self.a, normalize=True)
        
        return prism,
        
    def animation(self, name = 'hexagonal_prism_rotation_aspect.gif', save = False):
        ani = FuncAnimation(self.fig, self.animate, frames=self.frames, interval=50, blit=False)
        # Run the animation
        if save == True and self.out != None: 
            ani.save(self.out+'\\'+name, writer='Pillow', fps=30)
    plt.show()
