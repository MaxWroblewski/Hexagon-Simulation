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
            self.out = self.path+'_Ix_'+str(round(self.I1,4))+'Iy_'+str(round(self.I2,4))+'Iz_'+str(round(self.I3,4))+'wx0_'+str(round(self.w[0],4))+'wy0_'+str(round(self.w[1],4))+'wz0_'+str(round(self.w[2],4))+'dt_'+str(self.dt)
            if not os.path.exists(self.out):
                os.makedirs(self.out)
        else:
            self.out = None




    #Simulation code

    #Find w_dot given w (from previous step) and I (constant)
    def find_w_dot(self, w: np.ndarray):
        wx_dot = (self.I2-self.I3)*w[1]*w[2]/self.I1
        wy_dot = (self.I3-self.I1)*w[2]*w[0]/self.I2
        wz_dot = (self.I1-self.I2)*w[0]*w[1]/self.I3
        w_dot = [wx_dot, wy_dot, wz_dot]
        return w_dot

    #Update w given new w_dot
    def update_omega(self, w: np.ndarray, w_dot:np.ndarray, dt: float):
        w[0] = w[0] + w_dot[0]*dt
        w[1] = w[1] + w_dot[1]*dt
        w[2] = w[2] + w_dot[2]*dt
        return w

    #Def magnitude for use in calculation
    def mag(self, v):
        return np.sqrt(v.dot(v))

    #Calculate angular momentum for analysis
    def calculate_L(self, w, I):
        Lx = I[0][0]*w[0]
        Ly = I[1][1]*w[1]
        Lz = I[2][2]*w[2]
        L = np.array([Lx, Ly, Lz])
        L_mag = self.mag(L)
        return L, L_mag

    #This runs a simulation
    def run(self, save = False):
        #Store data in here temporarily
        data = []

        #Do run over given number of steps
        for i in range(self.steps):
            #Find w_dot[i]
            self.w_dot = self.find_w_dot(self.w)
            #Calculate L[i] values for analysis
            self.L, self.L_mag = self.calculate_L(self.w, self.I)
            #Put all [i] data into the i-th row of the df
            row = {'Time': i*self.dt, 
            'wx': self.w[0], 'w_dotx': self.w_dot[0], 'Lx': self.L[0],
            'wy': self.w[1], 'w_doty': self.w_dot[1], 'Ly': self.L[1],
            'wz': self.w[2], 'w_dotz': self.w_dot[2], 'Lz': self.L[2],
            'L_mag':self.L_mag
            }
            #Calculate w[i+1]
            data.append(row)
            self.w = self.update_omega(self.w, self.w_dot, self.dt)

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
                
