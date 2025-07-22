# Solver to find stationary Schrödinger solution

### Time-dependent Schrödinger equation $\to$ Stationary Schrödinger equation

The Schrödinger equation is a partial differential equation that governs the wave function of a non-relativistic quantum-mechanical system.
\begin{equation}
    i \hbar \frac{d}{dt} \psi = \hat{H} \psi 
\end{equation}
with $\psi$ is the time-dependent wave function, $\hat{H}$ is the Hamiltonien operator and $\hbar$ is the reduced Planck's constant.

We assume that the wave function can be separated into a spatial part $\psi(\mathbf{r})$ and a temporal part $T(t)$:
\begin{equation*}
\Psi(\mathbf{r}, t) = \psi(\mathbf{r}) T(t)
\end{equation*}

Substituting this expression into the time-dependent Schrödinger equation gives us:
\begin{equation*}
i\hbar \psi(\mathbf{r}) \frac{\partial T(t)}{\partial t} = T(t) \hat{H} \psi(\mathbf{r})
\end{equation*}
We rearrange and we find:
\begin{equation*}
i\hbar \frac{1}{T(t)} \frac{\partial T(t)}{\partial t} = \frac{1}{\psi(\mathbf{r})} \hat{H} \psi(\mathbf{r})
\end{equation*}
The left-hand side of the equation depends only on time, while the right-hand side depends only on position. For this equality to be true for all \mathbf{r}\s and \(t\), both sides must be equal to a constant, which we call \(E\) (the energy of the system):

\begin{equation*}
i\hbar \frac{1}{T(t)} \frac{\partial T(t)}{\partial t} = E \quad \text{et} \quad \frac{1}{\psi(\mathbf{r})} \hat{H} \psi(\mathbf{r}) = E
\end{equation*}

The first equation gives a solution for the time part:
\begin{equation*}
T(t) = e^{-iEt/\hbar}
\end{equation*}


The second equation is the time-independent Schrödinger equation:
\begin{equation*}
\hat{H} \psi(\mathbf{r}) = E \psi(\mathbf{r})
\end{equation*}

This equation describes the stationary states of the quantum system, where \(E\) is the energy of the stationary state corresponding to the wave function \(\psi(\mathbf{r})\).

### Radial differential equation

We start from the time-independent Schrödinger equation in three dimensions:

\begin{equation*}
\hat{H} \psi(\mathbf{r}) = E \psi(\mathbf{r})
\end{equation*}

where H^H^ is the Hamiltonian operator given by:

\begin{equation*}
\hat{H} = -\frac{\hbar^2}{2m} \nabla^2 + V(r)
\end{equation*}

In spherical coordinates, the Laplacian ∇2∇2 is given by:

\begin{equation*}
\nabla^2 = \frac{1}{r^2} \frac{\partial}{\partial r} \left( r^2 \frac{\partial}{\partial r} \right) + \frac{1}{r^2 \sin \theta} \frac{\partial}{\partial \theta} \left( \sin \theta \frac{\partial}{\partial \theta} \right) + \frac{1}{r^2 \sin^2 \theta} \frac{\partial^2}{\partial \phi^2}
\end{equation*}

We assume that the wavefunction ψ(r)ψ(r) can be separated into a radial part R(r)R(r) and an angular part Y(θ,ϕ)Y(θ,ϕ):

\begin{equation*}
\psi(\mathbf{r}) = R(r) Y(\theta, \phi)
\end{equation*}

Substituting this expression into the Schrödinger equation, we obtain:

\begin{equation*}
-\frac{\hbar^2}{2m} \left[ Y(\theta, \phi) \frac{1}{r^2} \frac{d}{dr} \left( r^2 \frac{dR(r)}{dr} \right) + \frac{R(r)}{r^2} \nabla^2_{\theta, \phi} Y(\theta, \phi) \right] + V(r) R(r) Y(\theta, \phi) = E R(r) Y(\theta, \phi)
\end{equation*}

where ∇θ,ϕ2∇θ,ϕ2​ is the angular part of the Laplacian.

Dividing by R(r)Y(θ,ϕ)R(r)Y(θ,ϕ) and multiplying by −2mr2ℏ2−ℏ22mr2​, we obtain:

\begin{equation*}
\frac{1}{R(r)} \frac{d}{dr} \left( r^2 \frac{dR(r)}{dr} \right) - \frac{2m r^2}{\hbar^2} (V(r) - E) = \frac{1}{Y(\theta, \phi)} \nabla^2_{\theta, \phi} Y(\theta, \phi)
\end{equation*}

The right-hand side of this equation depends only on the angular variables θθ and ϕϕ, while the left-hand side depends only on the radial variable rr. Therefore, both sides must be equal to a constant, which we denote by −l(l+1)−l(l+1):

\begin{equation*}
\frac{1}{Y(\theta, \phi)} \nabla^2_{\theta, \phi} Y(\theta, \phi) = -l(l+1)
\end{equation*}

By performing the change of function u(r)=rR(r)u(r)=rR(r), and considering atomic units where ℏ=1ℏ=1 and m=1m=1, we obtain the following equation:

\begin{equation*}
-\frac{1}{2} \frac{d^2 u(r)}{dr^2} + \left[ V(r) + \frac{l(l+1)}{2r^2} \right] u(r) = E u(r)
\end{equation*}

### Tridiagonal matrix

In this code, we will numerically solve the differential equation above on a regular grid ranging from rminrmin​ to rmaxrmax​, with a step size hh. We assume the grid consists of NN points, with h=(rmax−rmin)/Nh=(rmax​−rmin​)/N. For this purpose, we define ri=rmin+ihri​=rmin​+ih and ui=u(ri)ui​=u(ri​), for 0≤i≤N0≤i≤N. In this case, using a central difference scheme, the first and second derivatives ui′=u′(ri)ui′​=u′(ri​) and ui′′=u′′(ri)ui′′​=u′′(ri​) are given by:

\begin{equation}
\begin{gathered}
u_i' = \frac{u_{i+1} - u_{i-1}}{2h} \
u_i'' = \frac{u_{i+1} - 2u_i + u_{i-1}}{h^2}
\end{gathered}
\end{equation}

The discretized equation then becomes, for 1≤i≤N−11≤i≤N−1 (with zero boundary conditions):

\begin{equation}
-\frac{u_{i+1} - 2u_i + u_{i-1}}{2h^2} + \left[\frac{\ell(\ell + 1)}{2r_i^2} + V(r_i)\right] u_i = E u_i.
\end{equation}

Let us define a=−1/2h2a=−1/2h2 and di=1/h2+ℓ(ℓ+1)/2ri2+V(ri)di​=1/h2+ℓ(ℓ+1)/2ri2​+V(ri​). The discretized form of the Hamiltonian HH is then written as a tridiagonal matrix of the form:

\begin{equation}
H = \begin{bmatrix}
d_1 & a & 0 & 0 & \cdots & 0 \
a & d_2 & a & 0 & \cdots & 0 \
0 & a & d_3 & a & \cdots & 0 \
\vdots & \ddots & \ddots & \ddots & \ddots & \vdots \
0 & \cdots & 0 & a & d_{N-1} & a \
0 & \cdots & 0 & 0 & a & d_N \
\end{bmatrix}
\end{equation}

This matrix is symmetric. We have used the boundary conditions u(0)=u(rmax)=0u(0)=u(rmax​)=0. The task now is to solve the eigenvalue problem, i.e., to diagonalize this tridiagonal matrix.