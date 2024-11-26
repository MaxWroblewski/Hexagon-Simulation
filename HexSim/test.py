from form_matrix import Hexagonal_Prism
from evolve_rotation import Run_Rotation
from animate_hex import Animate_Hex
import numpy as np

# Parameters
a = 2.5e-6  # Side length of the hexagon
l = 0.27e-6  # Thickness of the prism
rho_0 = 3800  # Uniform density
n_theta = 100  # Angular resolution
n_r = 100  # Radial resolution
n_z = 50  # Z resolution

#Define a function for the mass distribution
def rho(x, y):
    return rho_0 * (1 + 0.05 * max(0,x) + 0.005 * max(0,y))

prism = Hexagonal_Prism(a, rho, l, rho_0, n_theta, n_r, n_z) #Initialize moment of inertia calculation
matrix, _, __, ___ = prism.form_matrix() #Output moment of inertia matrix

print(matrix) #Print moment of inertia matrix


#T, f = prism.run_eval(10)
#print(T)
#print(f)

w = np.array([0.001, 10, 0.001])
dt = 0.0005
steps = 500000



run = Run_Rotation(I = matrix, w = w, dt = dt, steps = steps, path = 'test')
df = run.run(save = True)
print(run.find_w_dot())
run.plot_poinsot()


#print(df)
#run.plot_L(save = True, total = True)


#anim = Animate_Hex(timestep=1, frames = 3000, a=a, h=l, I=matrix, df = df, path = 'test')
#animate = anim.animation(save = True)

