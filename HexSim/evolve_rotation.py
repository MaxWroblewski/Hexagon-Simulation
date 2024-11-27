import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import os

class Run_Rotation:

    def __init__(self, I:np.ndarray, w:np.ndarray, dt: float, steps: int, path:str=None):

        #Def simulation parameters
        self.I = I
        self.w = w
        self.dt = dt
        self.I1 = I[0][0]
        self.I2 = I[1][1]
        self.I3 = I[2][2]
        self.steps = steps
        self.path = path
        
        if self.path != None:
            self.out = self.path+'_Ix_'+str(np.round(self.I1,4))+'Iy_'+str(np.round(self.I2,4))+'Iz_'+str(np.round(self.I3,4))+'wx0_'+str(np.round(self.w[0],4))+'wy0_'+str(np.round(self.w[1],4))+'wz0_'+str(np.round(self.w[2],4))+'dt_'+str(self.dt)
            if not os.path.exists(self.out):
                os.makedirs(self.out)
        else:
            self.out = None




    #Simulation code

    #Find w_dot given w (from previous step) and I (constant)
    def find_w_dot(self):
        wx_dot = (self.I2-self.I3)*self.w[1]*self.w[2]/self.I1
        wy_dot = (self.I3-self.I1)*self.w[2]*self.w[0]/self.I2
        wz_dot = (self.I1-self.I2)*self.w[0]*self.w[1]/self.I3
        self.w_dot = np.array([wx_dot, wy_dot, wz_dot])
        return self.w_dot

    #Update w given new w_dot
    #def update_omega(self):
    #    self.w = self.w + self.dt*self.w_dot
    #    #for i in range(len(self.w)):
    #    #    self.w[i] = self.w[i] + self.w_dot[i]*dt
    #    #w[0] = w[0] + w_dot[0]*dt
    #    #w[1] = w[1] + w_dot[1]*dt
    #    #w[2] = w[2] + w_dot[2]*dt
    #    return self.w

    def update_omega(self):
        w0 = self.w
        # Runge-Kutta 4th Order Integration
        k1 = self.dt * self.find_w_dot()
        self.w = self.w + 0.5 * k1
        k2 = self.dt * self.find_w_dot()
        self.w = self.w + 0.5 * k2
        k3 = self.dt * self.find_w_dot()
        self.w = self.w + k3
        k4 = self.dt * self.find_w_dot()

        # Update angular velocity
        self.w = w0 + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)

        return self.w

    #Def magnitude for use in calculation
    def mag(self, v):
        return np.sqrt(v.dot(v))

    #Calculate angular momentum for analysis
    def calculate_L(self):
        Lx = self.I1*self.w[0]
        Ly = self.I2*self.w[1]
        Lz = self.I3*self.w[2]
        L = np.array([Lx, Ly, Lz])
        L_mag = self.mag(L)
        return L, L_mag

    def calculate_energy(self):
        E = 0.5 * (self.I1 * self.w[0]**2 + 
                   self.I2 * self.w[1]**2 + 
                   self.I3 * self.w[2]**2)
        return E

    #This runs a simulation
    def run(self, save = False):
        #Store data in here temporarily
        data = []

        #Do run over given number of steps
        for i in range(self.steps):
            #Find w_dot[i]
            self.w_dot = self.find_w_dot()
            #Calculate L[i] values for analysis
            self.L, self.L_mag = self.calculate_L()
            self.E = self.calculate_energy()
            #Put all [i] data into the i-th row of the df
            row = {'Time': i*self.dt, 
            'wx': self.w[0], 'w_dotx': self.w_dot[0], 'Lx': self.L[0],
            'wy': self.w[1], 'w_doty': self.w_dot[1], 'Ly': self.L[1],
            'wz': self.w[2], 'w_dotz': self.w_dot[2], 'Lz': self.L[2],
            'L_mag':self.L_mag, 'T': self.E
            }
            #Calculate w[i+1]
            data.append(row)
            self.w = self.update_omega()

        df = pd.DataFrame(data)
        self.df = df
        data = []
        if save == True and self.path != None and save != False:
            df.to_csv(self.out+'\\data.csv')
        return df
    
    def plot_L(self, save = False, total = False):
        
        plt.figure()

        plt.rcParams['font.family'] = 'Corbel'
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.titlesize'] = 16
        plt.rcParams['axes.labelsize'] = 14
        plt.rcParams['legend.fontsize'] = 12
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 12
        plt.rcParams['axes.grid'] = True  # Default grid for all plots
        plt.rcParams['grid.alpha'] = 0.6  # Grid transparency
        plt.rcParams['grid.linestyle'] = '--'
        
        plt.rcParams.update({
            'font.family': 'DejaVu Sans',  # Use 'Helvetica' or 'Nimbus Sans L' if preferred
            'font.size': 12,         # Adjust font size for labels and annotations
            'axes.labelsize': 14,    # Axis labels
            'axes.titlesize': 16,    # Title size
            'legend.fontsize': 12,   # Legend
            'xtick.labelsize': 12,   # X-axis tick labels
            'ytick.labelsize': 12,   # Y-axis tick labels
        })

        title = 'Angular Momentum over '+str(round(self.df['Time'][self.df.index[-1]]))+' Seconds for $\omega_{0y}=$'+str(round(self.df['wy'][0], 4))
        plt.title(title)
        plt.ylabel('Angular Momentum')
        plt.xlabel('Time (s)')
        scale = 'linear'
        plt.yscale(scale)
        plt.scatter(self.df['Time'], self.df['Lx'], s = 0.5, label = '$L_x$')
        plt.scatter(self.df['Time'], self.df['Ly'], s = 0.5, label = '$L_y$')
        plt.scatter(self.df['Time'], self.df['Lz'], s = 0.5, label = '$L_z$')
        plt.legend()

        if save == True and self.out != None:
            plt.savefig(self.out+'\\len_'+str(self.steps)+'.png')
            plt.show()
        else:
            plt.show()

        if total == True:
            plt.figure()
            plt.title('Total Angular Momentum over '+str(round(self.df['Time'][self.df.index[-1]]))+' Seconds for $\omega_{0y}=$'+str(self.df['wy'][0]))
            plt.scatter(self.df['Time'], self.df['L_mag'], s = 0.5)
            if save == True and self.out != None:
                plt.savefig(self.out+'\\len_'+str(self.steps)+'totalL.png')
                plt.show()
            else:
                plt.show()
    


    def plot_energy(self, save=False):
        plt.figure()
        plt.plot(self.df['Time'], self.df['T'], label='Kinetic Energy', color='purple')
        plt.xlabel('Time (s)')
        plt.ylabel('Energy')
        plt.title('Conservation of Kinetic Energy')
        plt.legend()
        plt.grid(True)

        if save == True and self.out != None:
            plt.savefig(os.path.join(self.out, 'energy_conservation.png'))
            plt.show()
        else:
            plt.show()

    def plot_poinsot(self, save=False):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        # Calculate the initial kinetic energy E
        E = self.df['T'][0]
        # Generate the inertia ellipsoid surface
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        u, v = np.meshgrid(u, v)

        # Compute the scaling factors for the ellipsoid axes
        a = np.sqrt(2 * E / self.I1)
        b = np.sqrt(2 * E / self.I2)
        c = np.sqrt(2 * E / self.I3)

        # Parametric equations of the ellipsoid
        x = a * np.sin(v) * np.cos(u)
        y = b * np.sin(v) * np.sin(u)
        z = c * np.cos(v)

        # Get the omega components from the dataframe
        omega1 = self.df['wx']
        omega2 = self.df['wy']
        omega3 = self.df['wz']

        # Plot the inertia ellipsoid and polhode curve
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        # Plot the inertia ellipsoid
        ax.plot_surface(x, y, z, color='orange', alpha=0.5, rstride=5, cstride=5, edgecolor='none')

        # Plot the polhode curve
        ax.plot(omega1, omega2, omega3, color='blue', label='Polhode Curve', linewidth=2)

        # Set labels and title
        ax.set_xlabel('$\omega_1$')
        ax.set_ylabel('$\omega_2$')
        ax.set_zlabel('$\omega_3$')
        title = f"Poinsot's Construction: Inertia Ellipsoid and Polhode Curve"
        ax.set_title(title)
        ax.legend()

        if save == True and self.out != None:
            plt.savefig(os.path.join(self.out, 'poinsot_plot.png'))
            plt.show()
        else:
            plt.show()
                
