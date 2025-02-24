# Hexagon-Simulation
## Overview
For many applications, studying the rotational dynamics of irregularly shaped objects can be approximated to a more familiar and symmetric shape. However, in systems which exhibit the tennis racket effect [link to tennis racket effect], we observe that there is often a much more significant dependence on the accuracy of moment of inertia tensor calculations. This heavier dependence on more accurate calculations is further exagurated when an object is highly symmetric about two or three axis [2]. 

We study the rotational dynamics of hexagonal prisms and configurations of multiple hexagonal prisms to understand the evolution of a system with these inertial properties exhibiting the tennis racket effect. 


2*

```latex
% !TEX TS-program = pdflatex
\documentclass{article}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage[margin=1in]{geometry}

\title{Sensitivity of the Moment of Inertia Tensor in the Tennis Racket Effect}
\author{}
\date{}

\begin{document}

\maketitle

Consider an object with moment of inertia tensor \(I\) that has eigenvalues 
\(\lambda_1 < \lambda_2 < \lambda_3\) corresponding to the principal axes. For a 
small perturbation \(\delta I\), first-order perturbation theory gives the change in 
the eigenvalues as
\[
\delta \lambda_i = \langle \psi_i \mid \delta I \mid \psi_i \rangle,
\]
where \(\psi_i\) is the eigenvector associated with \(\lambda_i\).

Now, consider the relative error for the intermediate moment \(\lambda_2\). Its 
sensitivity is characterized by the ratio
\[
\frac{\delta \lambda_2}{\Delta \lambda} \sim \frac{\langle \psi_2 \mid \delta I \mid \psi_2 \rangle}{\lambda_3 - \lambda_2},
\]
assuming that \(\lambda_3 - \lambda_2\) is small.

In the context of the tennis racket effect, Euler’s equations for rotation are given 
by:
\[
\begin{aligned}
\dot{\omega}_1 &= \frac{(\lambda_2 - \lambda_3)}{\lambda_1}\, \omega_2 \omega_3, \\
\dot{\omega}_2 &= \frac{(\lambda_3 - \lambda_1)}{\lambda_2}\, \omega_3 \omega_1, \\
\dot{\omega}_3 &= \frac{(\lambda_1 - \lambda_2)}{\lambda_3}\, \omega_1 \omega_2.
\end{aligned}
\]
A small miscalculation in \(\lambda_2\) can therefore lead to a large error in predicting 
the stability of rotation about the intermediate axis. This demonstrates the high 
sensitivity of irrationally symmetric objects, where the eigenvalues do not exhibit 
simple ratios.

\end{document}
```
