# Hexagon-Simulation
## Overview
In many applications, studying the rotational dynamics of irregularly shaped objects can be approximated using more familiar and symmetric shapes. However, in systems exhibiting the tennis racket effect [link to tennis racket effect], the accuracy of moment of inertia tensor calculations plays a crucial role. This dependency becomes even more pronounced when an object is highly symmetric about two or three axes.

<details> <summary>Sensitivity of the Moment of Inertia Tensor in the Tennis Racket Effect</summary>

___

Consider an object with moment of inertia tensor $I$ that has eigenvalues 
$\lambda_1 < \lambda_2 < \lambda_3$ corresponding to the principal axes.
For a small perturbation $\delta I$, first-order perturbation theory gives
the change in the eigenvalues as:

![Eigenvalue Variation](https://latex.codecogs.com/png.latex?\delta\lambda_i%20=%20\langle\psi_i%20\mid%20\delta%20I%20\mid%20\psi_i\rangle)

where $\psi_i\$ is the eigenvector associated with $\lambda_i$.

Now, consider the relative error for the intermediate moment $\lambda_2$.
Its sensitivity is characterized by the ratio:

![Ratio of Eigenvalue Variations](https://latex.codecogs.com/png.latex?\frac{\delta\lambda_2}{\Delta\lambda}%20\sim%20\frac{\langle\psi_2%20\mid%20\delta%20I%20\mid%20\psi_2\rangle}{\lambda_3%20-%20\lambda_2})

assuming that $\lambda_3 - \lambda_2$ is small.

In the context of the tennis racket effect, Euler’s equations for rotation are:


![Euler's Equations](https://latex.codecogs.com/png.latex?\dot{\omega}_1%20=%20\frac{(\lambda_2%20-%20\lambda_3)}{\lambda_1}%20\omega_2%20\omega_3,%20\quad%20\dot{\omega}_2%20=%20\frac{(\lambda_3%20-%20\lambda_1)}{\lambda_2}%20\omega_3%20\omega_1,%20\quad%20\dot{\omega}_3%20=%20\frac{(\lambda_1%20-%20\lambda_2)}{\lambda_3}%20\omega_1%20\omega_2.)


A small miscalculation in $\lambda_2$ can therefore lead to a large error
in predicting the stability of rotation about the intermediate axis.
This demonstrates the high sensitivity of rotationally asymmetric objects,
where the eigenvalues do not exhibit simple ratios.

___

</details>


We investigate the rotational dynamics of hexagonal prisms and configurations of multiple hexagonal prisms to understand the evolution of a system with these inertial properties exhibiting the tennis racket effect.

## Features
- [X] Simulation of hexagonal prism rotations

- [X] Visualization of the tennis racket effect in 3D

- [X] Customizable parameters for moment of inertia

- [ ] Analysis of stability and eigenvalue sensitivity

## Usage
[Add details on how to use here]

## Results
The project generates visual outputs showing the rotation dynamics and key numerical results related to the moment of inertia sensitivity.

## License
I dont know what this even means
