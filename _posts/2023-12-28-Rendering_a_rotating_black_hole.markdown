---
title: Rendering a rotating black hole
date: 2023-12-28
categories: [physics, general relativity, black holes]
tags: [black holes]
math: true
---

![Evolution of the renders](../images/miscellaneous/evolution.gif)
![Fly-by](../images/black_hole_renders/fly-by.gif)

## 1 Introduction

### 1.1 Goal

### 1.2 Topics we cover

### 1.3 General disclaimer

___

## 2 Fundamental black hole physics

### 2.1 Special relativity

### 2.2 General relativity and the metric 

#### 2.2.1 The Schwarzschild solution

#### 2.2.2 The Kerr solution

The spacetime around a rotating (charge-neutral) black hole is described by the *Kerr metric*. This metric is commonly given in Boyer-Lindquist (BL) coordinates, where the line element takes the following form

$$
\begin{equation}
\tag{1}
    ds^2 = g_{tt}dt^2 + g_{rr}dr^2 + g_{\theta\theta}d\theta^2 + g_{\phi\phi}d\phi^2 + g_{t\phi}\left(dtd\phi + d\phi dt\right),
\end{equation}
$$

with

$$
\begin{align*}
    g_{tt} &= -\left(1 - \frac{2Mr}{\Sigma}\right), \\ \\
    g_{rr} &= \frac{\Sigma}{\Delta}, \\ \\
    g_{\theta\theta} &= \Sigma, \\ \\
    g_{\phi\phi} &= \frac{\Lambda}{\Sigma}\sin^2 \theta,\\ \\
    g_{t\phi} &= -\frac{2M r a \sin^2 \theta}{\Sigma}
\end{align*}
$$

and

$$
\begin{align*}
    \Sigma &= r^2 + a^2\cos^2\theta, \\
    \Delta &= r^2 - 2Mr + a^2, \\
    \Lambda &= \left(r^2 + a^2\right)^2 - a^2 \Delta \sin^2\theta,
\end{align*}
$$


and where $a = \frac{J}{M}$ describes the angular momentum of the black hole, where $J$ is the angular momentum of the black hole and $M$ is its mass. The parameter $a$ can go from $a = 0$ to $\vert a \vert = M$. The transformation from BL to Cartesian coordinates is given by (see for example [1])

$$
\begin{align*}
    x &= \sqrt{r^2 + a^2}\sin\theta\cos\phi, \tag{2} \\
    y &= \sqrt{r^2 + a^2}\sin\theta\sin\phi, \tag{3} \\
    z &= r\cos\theta. \tag{4}
\end{align*}
$$

In order for the math to be cleaner in the code we'll introduce the dimensionless quantities

$$
\begin{equation}
  r' = \frac{r}{M}, \quad a' = \frac{a}{M}, \quad t' = \frac{t}{M}
\end{equation}
$$

Then

$$
\begin{align*}
  \Sigma &= M^2 \left(\left(r'\right)^2 + \left(a'\right)^2\cos^2\theta\right) \equiv M^2 \Sigma' \\ \\
  \Delta &= M^2 \left(\left(r'\right)^2 - 2r' + \left(a'\right)^2\right) \equiv M^2 \Delta', \\
  \Lambda &= M^4 \left[\left(\left(r'\right)^2 + \left(a'\right)^2\right)^2 - \left(a'\right)^2 \Delta' \sin^2\theta\right] \equiv M^4 \Lambda'
\end{align*}
$$

meaning that 

$$
\begin{align*}
  g_{tt}           &= -\left(1 - \frac{2M^2 r'}{M^2 \Sigma'}\right) = -\left(1 - \frac{2r'}{\Sigma'}\right) \\ \\
  g_{rr}           &= \frac{M^2 \Sigma'}{M^2 \Delta'} = \frac{\Sigma'}{\Delta'} \\ \\
  g_{\theta\theta} &= M^2 \Sigma' \\ \\
  g_{\phi \phi}    &= \frac{M^4 \Lambda'}{M^2 \Sigma'}\sin^2\theta = M^2 \frac{\Lambda'}{\Sigma'}\sin^2\theta \\ \\
  g_{t\phi}        &= -\frac{2M r' a' \sin^2\theta}{\Sigma'},
\end{align*}
$$

in which case the line element can be written as 

$$
\begin{align*}
  ds^2 &= -\left(1 - \frac{2r'}{\Sigma'}\right) M^2 \left(dt'\right)^2 + \frac{\Sigma'}{\Delta'}M^2 \left(dr'\right)^2 + M^2 \Sigma' d\theta^2  \\ \\
  &+ M^2 \frac{\Lambda'}{\Sigma'}\sin^2\theta d\phi^2 - M^2 \frac{4r' a'\sin^2\theta}{\Sigma'} d\phi dt'
\end{align*}
$$

If we now divide everything by $M^2$ we can completely eliminate it from the equations by also defining $ds' = \frac{ds}{M}$. For reasons of brevity we will omit all the primes with the implicit understanding that all quantities are dimensionless. And if we ever want to convert to normal physics units we need to multiply by $M$ and add appropriate factors of $G$ and $c$.

#### 2.2.3 Coordinate transformation between Cartesian and Boyer-Lindquist coordinates

We will need to convert vectors between Cartesian and BL coordinates multiple times, so a coordinate transformation between these systems is necessary. Deriving this is most easily done by studying an infinitesimal displacement vector $d\vec{r}$ expressed in the two coordinate systems. This can be written as 

$$
\begin{equation}
  d\vec{r} = dx \; \vec{\hat{x}} + dy \; \vec{\hat{y}} + dz \; \vec{\hat{z}}
\end{equation}
$$

in Cartesian, and 

$$
\begin{equation}
  d\vec{r} = Rdr \; \vec{\hat{r}} + \Theta d\theta \; \vec{\hat{\theta}} + \Phi d\phi \; \vec{\hat{\phi}}
\end{equation}
$$

in BL coordinates. Here $\vec{\hat{r}}, \vec{\hat{\theta}}$ and $\vec{\hat{\phi}}$ are the unit vectors in the BL coordinate system and $R, \Theta$ and $\Phi$ are scaling parameters defined such that $\vec{\hat{r}}, \vec{\hat{\theta}}$ and $\vec{\hat{\phi}}$ are unit vectors. These two expressions describe the same displacement vector, so they are necessarily equal. The easiest way of finding the expressions for the unit vectors is to write $dx, dy$ and $dz$ as total differentials with respect to $r, \theta$ and $\phi$. Then we set the two formulations of $d\vec{r}$ equal to each other and match the coefficients of the $dr, d\theta$ and $d\phi$ terms on both sides. But first, for the total differentials we get

$$
\begin{align*}
  dx &= \frac{\partial x}{\partial r} dr + \frac{\partial x}{\partial \theta}d\theta + \frac{\partial x}{\partial \phi}d\phi \\ \\
  &= \frac{r}{\sqrt{r^2 + a^2}}\sin\theta\cos\phi \; dr + \sqrt{r^2 + a^2}\cos\theta \cos\phi \; d\theta - \sqrt{r^2 + a^2}\sin\theta \sin\phi \; d\phi \\ \\
  dy &= \frac{\partial y}{\partial r} dr + \frac{\partial y}{\partial \theta} d\theta + \frac{\partial y}{\partial \phi}d\phi \\ \\
  &= \frac{r}{\sqrt{r^2 + a^2}}\sin\theta\sin\phi \; dr + \sqrt{r^2 + a^2}\cos\theta \sin\phi \: d\theta + \sqrt{r^2 + a^2}\sin\theta \cos\phi \; d\phi \\ \\
  dz &= \frac{\partial z}{\partial r}dr + \frac{\partial z}{\partial \theta}d\theta + \frac{\partial z}{\partial \phi}d\phi \\ \\
  &= \cos\theta \; dr - r\sin\theta \; d\theta
\end{align*}
$$

Matching the terms containing $dr$ in $\eqref{eq: infinitesimal_displacement_Cartesian}$ and $\eqref{eq: infinitesimal_displacement_spheroidal}$, by inserting the terms in $dx, dy$ and $dz$ containing $dr$ into $\eqref{eq: infinitesimal_displacement_Cartesian}$ we get

$$
\begin{equation}
  R\vec{\hat{r}} = \frac{r}{\sqrt{r^2 + a^2}}\sin\theta \cos\phi \; \vec{\hat{x}} + \frac{r}{\sqrt{r^2 + a^2}}\sin\theta \sin\phi \; \vec{\hat{y}} + \cos\theta \;\vec{\hat{z}}
\end{equation}
$$

We can then determine the scaling parameter $R$ by finding the norm of $R\vec{\hat{r}}$ and using the defining property $\left|\vec{\hat{r}}\right| = 1$.

$$
\begin{align*}
  R^2 &= \frac{r^2}{r^2 + a^2}\sin^2\theta + \cos^2\theta \\ \\
  &= \frac{r^2 + a^2 \cos^2\theta}{r^2 + a^2} \\ \\
  \implies R &= \sqrt{\frac{r^2 + a^2\cos^2\theta}{r^2 + a^2}}
\end{align*}
$$

And thus the unit vector $\vec{\hat{r}}$ is given by

$$
\begin{equation}
  \vec{\hat{r}} = \frac{1}{\sqrt{r^2 + a^2\cos^2\theta}}\left[r\sin\theta \cos\phi \; \vec{\hat{x}} + r\sin\theta \sin\phi \; \vec{\hat{y}} + \sqrt{r^2 + a^2}\cos\theta \; \vec{\hat{z}}\right].
\end{equation}
$$

Completely analogous procedures can be followed for $\theta$ and $\phi$ by matching the $d\theta$ terms on both sides and then the $d\phi$ terms on both sides. We will skip the steps here, because they are exactly the same as for $\vec{\hat{r}}$, but we end up with the following

$$
\begin{align*}
  \vec{\hat{\theta}} &= \frac{1}{\sqrt{r^2 + a^2 \cos^2\theta}}\left[\sqrt{r^2 + a^2}\cos\theta \cos\phi \; \vec{\hat{x}} + \sqrt{r^2 + a^2}\cos\theta \sin\phi \; \vec{\hat{y}} - r\sin\theta \; \vec{\hat{z}}\right] \\ \\
  \vec{\hat{\phi}} &= -\sin\phi \;\vec{\hat{x}} + \cos\phi \;\vec{\hat{y}}.
\end{align*}
$$

We can therefore write the relation between the BL unit vectors and the Cartesian ones in a matrix form as follows

$$
\begin{equation}  
  \left[\begin{matrix}
          \vec{\hat{r}} \\ \\
          \vec{\hat{\theta}} \\ \\
          \vec{\hat{\phi}}
        \end{matrix}\right] = \left[\begin{matrix}
                                \frac{r\sin\theta \cos\phi}{\sqrt{r^2 + a^2\cos^2\theta}} & \frac{r\sin\theta \sin\phi}{\sqrt{r^2 + a^2\cos^2\theta}} & \frac{\sqrt{r^2 + a^2}\cos\theta}{\sqrt{r^2 + a^2\cos^2\theta}} \\ \\
                                \frac{\sqrt{r^2 + a^2}\cos\theta \cos\phi}{\sqrt{r^2 + a^2\cos^2\theta}} & \frac{\sqrt{r^2 + a^2}\cos\theta \sin\phi}{\sqrt{r^2 + a^2\cos^2\theta}} & -\frac{r\sin\theta}{\sqrt{r^2 + a^2\cos^2\theta}} \\ \\
                                -\sin\phi & \cos\phi & 0
                                \end{matrix}\right] \left[\begin{matrix}
                                                \vec{\hat{x}} \\ \\ 
                                                \vec{\hat{y}} \\ \\
                                                \vec{\hat{z}}
                                                          \end{matrix}\right]
\end{equation}
$$

Let us call this coordinate transformation matrix $M$. It is fairly easy to check that this is an orthogonal matrix, meaning that $M^{-1} = M^T$ (just compute $M^T M$, which is equal to $\mathbb{I}$ for orthogonal matrices $M$). And in that case the inverse transformation, from BL to Cartesian coordinates, is given by $M^T$. Now, of course this is the coordinate transformation between the Cartesian and BL *basis vectors*. And we want the transformation between *vector components*. But it turns out the coordinate transformation for the vector components is the same as for the basis vectors, as is shown by Lutz Lehmann [here](https://math.stackexchange.com/questions/3493647/do-coordinate-components-transform-in-the-same-or-opposite-way-as-their-bases). 
___

## 3 General Relativistic Raymarching

#### 3.1 The geodesic equation 

In order to visualize the black hole we need to trace *geodesics* in the Kerr spacetime, which we do by solving the *geodesic equation*

$$
\begin{equation}
  \frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\rho \sigma}\frac{dx^\rho}{d\lambda}\frac{dx^\sigma}{d\lambda} = 0,
\end{equation}
$$

where $x^\mu$ is the four-position of a photon, and $\lambda$ is an affine parameter for the geodesic, chosen such that $p^\mu \equiv \frac{dx^\mu}{d\lambda}$ is the four-momentum of the photon. $\Gamma^\mu_{\rho \sigma}$ are the *Christoffel symbols* of the metric. Equation $\eqref{eq: geodesic_equation}$ is completely general, and we of course need explicit expressions for the geodesic equation for each value of $\mu$. But computing the Christoffel symbols tends to be very tedious work, and certainly so for the Kerr metric. So instead of deriving the Christoffel symbols by hand we use the Sympy package in Python to derive them, and thus also explicit expressions for each component of the geodesic equation. Code for doing this can be found [here](https://github.com/chrvill/Black-hole-GR).

We have the normalization requirement

$$
\begin{equation}
    p_\nu p^\nu = \mu
\end{equation}
$$

with $\mu = 0$ for massless particles like photons and $\mu = -1$ for massive particles (This is really the normalization of the four-*velocity* of a massive particle, not the four-momentum. But the difference is only a factor of $m^2$, with $m$ being the mass of the particle in question.). This can be written explicitly as 

$$
\begin{equation}
    g_{tt} \left(p^t\right)^2 + g_{rr}\left(p^r\right)^2 + g_{\theta\theta} \left(p^\theta\right)^2 + g_{\phi \phi}\left(p^\phi\right)^2 + 2g_{t\phi} p^t p^\phi = \mu
\end{equation}
$$

When we choose initial directions for our rays we will essentially provide the spatial components of the four-momentum, so we need to compute $p^t$ from $p^r, p^\theta$ and $p^\phi$. Solving \eqref{eq: normalization_explicit} for $p^t$ gives

$$
\begin{equation}
    p^t = -\frac{g_{t\phi}}{g_{tt}} p^\phi \pm \sqrt{\left(\frac{g_{t\phi}}{g_{tt}} p^\phi\right)^2 - \frac{1}{g_{tt}}\left(g_{rr} \left(p^r\right)^2 + g_{\theta\theta} \left(p^\theta\right)^2 + g_{\phi \phi} \left(p^\phi\right)^2 - \mu\right)}
\end{equation}
$$

But we have to find out which sign is appropriate. Consider a case where $p^\phi = 0$, which is only possible outside the ergosphere. Then, since $p^t$ needs to be positive we need to use the $+$ sign outside the ergosphere. Arguing for which sign is appropriate inside the ergosphere proved trickier, but we found that we had to use the $-$ sign inside the ergosphere in order for $p^t$ to increase continuously as you pass the boundary of the ergosphere. 

### 3.2 Initial conditions

### 3.3 Numerically integrating the equations of motion

Numerically \eqref{eq: geodesic_equation} is just a completely standard second order differential equation. So in that sense solving it is just like solving any other second order differential equation numerically. However it should be noted that some of the Christoffel symbols will diverge as the photon approaches the event horizon, which needs to be taken into consideration when solving the geodesic equation. Assume for the purposes of this discussion that we have a first order differential equation of the form

$$
\begin{equation}
  \frac{dy}{dt} = f(t, y).
\end{equation}
$$

Our second order equations of motion can of course be expressed in terms of first order DEs by just by writing the corresponding equations for $\frac{dx^\mu}{d\lambda}$ and $\frac{d^2 x^\mu}{d\lambda^2}$. With an appropriate choice of variable timestep the geodesic equation could most likely be solved using for example the 4th order Runge-Kutta integration scheme. But here, due to an unrelated error that will be discussed later, we chose to use the \textit{Runge-Kutta-Fehlberg} scheme, abbreviated RKF45. This is, as the name implies, also in the family of Runge-Kutta methods. The main benefit of using this scheme is that it allows us to calculate an adaptive timestep $h$ very easily - the scheme itself can choose a small timestep when needed and revert to a bigger timestep when things evolve slowly. The scheme itself is similar to 4th order Runge-Kutta, and is described in detail [here]https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method), where we use the table listed under "FORMULA 2". It is a 4th order accurate method where we perform a few more function evaluations in order to obtain a measure of the error associated with the scheme. When the right hand side of the differential equation is not explicitly a function of time, like here, the algorithm can be written as

$$
\begin{align*}
  k_1 &= f\left(y\right) \\ \\
  k_2 &= f\left(y + B_{21} k_1 h\right) \\ \\
  k_3 &= f\left(y + B_{31} k_1 h + B_{32} k_2 h\right) \\ \\
  k_4 &= f\left(y + B_{41} k_1 h + B_{42} k_2 h + B_{43} k_3 h\right) \\ \\
  k_5 &= f\left(y + B_{51} k_1 h + B_{52} k_3 h + B_{53} k_3 h + B_{54} k_4 h\right) \\ \\
  k_6 &= f\left(y + B_{61} k_1 h + B_{62} k_2 h + B_{63} k_3 h + B_{64} k_4 h + B_{65} k_5 h\right),
\end{align*}
$$

which can be written more succinctly as 

$$
\begin{equation}
    k_n = f\left(y + \sum_{i = 1}^{n - 1} B_{ni}k_i h\right).
\end{equation}
$$

The approximation to the function value at the next step is then

$$
\begin{equation}
  y(t + h) = y(t) + CH_1 k_1 + CH_2 k_2 + CH_3 k_3 + CH_4 k_4 + CH_5 k_5 + CH_6 k_6,
\end{equation}
$$

and the estimate of the error is

$$
\begin{equation}
  TE = \lvert CT_1 k_1 + CT_2 k_2 + CT_3 k_3 + CT_4 k_4 + CT_5 k_5 + CT_6 k_6 \rvert.
\end{equation}
$$

The matrix $B$ and vectors $CH$ and $CT$ can be found in the article linked to. Then after having computed the step we can calculate a new step size

$$
\begin{equation}
  h_\text{new} = 0.9 h \left(\frac{\epsilon}{TE}\right)^{1/5}
\end{equation}
$$

where $\epsilon$ is a tolerance value we can choose in order to specify the level of accuracy we want to achieve. Then, if $TE > \epsilon$ the error is too big and so we replace $h$ with $h_\text{new}$ and repeat the step. We perform this iteration until $TE < \epsilon$. And then in the next step we use this value of $h_\mathrm{new}$ as the new $h$.

### 3.4 Example geodesics 

The following animations show sets of 10 geodesics for massless particles around a Kerr black hole with $a = 1$. These are prograde orbits that start with an initial angular velocity in the same direction as the rotation of the black hole. 

![prograde-geodesics](../images/miscellaneous/geodesics_prograde.gif)

While these are retrograde orbits with initial angular velocity in the opposite direction. 

![retrograde-geodesics](../images/miscellaneous/geodesics_retrograde.gif)
___

## 4 A simple pin-hole camera model

### 4.1 Raymarching & Schwarzschild example

#### 4.1.1 Reference frames

### 4.2 Extension to Kerr and Boyer-Lindquist coordinates

### 4.3 Examples

### 4.4 2D thin disk

### 4.5 Problems and limitations

___

## 5 The non-stationary camera

### 5.1 The motivation

Up until this point the camera has been stationary (by which we mean stationary in the BL coordinates). But it's of course more general to let the camera move in any arbitrary way. This also introduces *relativistic_aberration* and *redshift*/*blueshift*, which we'll go into detail on in a bit. So in addition to just being more general, having a moving camera also induces very interesting visual effects. This more general setup will also very easily allow us to put the camera on a geodesic by just solving the geodesic equation for the camera (keeping in mind that the camera is of course massive, so $\mu = -1$). But before we get into the details on how we do this we should discuss the concept of a *ZAMO* and the *ergosphere*. 

### 5.2 The ZAMO and ergosphere

In GR we can find quantities which are conserved along geodesics by using the concept of \textit{Killing vectors} (note that this section is inspired by \cite{carroll}). These are vectors $K^\mu$ such that $K_\mu p^\mu$ is conserved along a geodesic. One of the Killing vectors of the Kerr metric is the following

$$
\begin{equation}
    R^\mu = \left(\partial_\phi\right)^\mu = \left(0, 0, 0, 1\right)
\end{equation}
$$

This Killing vector is related to rotational symmetry, and we can therefore identify the following conserved quantity analogous to an angular momentum

$$
\begin{align*}
    L &\equiv g_{\mu \nu} R^\mu p^\nu \\ \\
    &= g_{\phi \nu} p^\nu \\ \\
    &= g_{\phi \phi} p^\phi + g_{t\phi} p^t
\end{align*}
$$

along a geodesic. Consider now the expression for $L$ at $\theta = \pi/2$ as $r \to \infty$. Then $g_{t\phi} \to 0$ and $g_{\phi \phi} \to r^2$. So $L \to r^2 p^\phi$. And $p^\phi$ reduces to an angular velocity in flat space. So this is the standard definition of angular momentum in flat space. Meaning that the conserved quantity $L$ corresponds to the angular momentum of an object infinitely far away, in the asymptotically flat region of the spacetime. We can now consider releasing a particle from infinitely far away with $L = 0$ and letting it fall along a geodesic. Then we have $g_{t\phi} p^t + g_{\phi \phi} p^\phi = 0$, meaning that 

$$
\begin{equation}
  p^\phi = -\frac{g_{t\phi}}{g_{\phi \phi}}p^t,
\end{equation}
$$

which is to say that the particle gains an angular velocity even though it starts off falling with zero angular momentum. Observers following geodesics with $L = 0$ are called *Zero Angular Momentum Observers* (ZAMO). The phenomenon causing ZAMO frames to move along with the rotation of the black hole is called *frame dragging*. It turns out that frame dragging has quite dramatic consequences for the motion of particles around a rotating black hole. Consider a photon in the equatorial plane moving only in the $\phi$ direction. Then

$$
\begin{equation*}
  ds^2 = 0 = g_{tt} dt^2 + 2g_{t\phi} dtd\phi + g_{\phi \phi} d\phi^2 \quad \implies \quad g_{\phi \phi} \left(\frac{d\phi}{dt}\right)^2 + 2g_{t\phi} \frac{d\phi}{dt} + g_{tt} = 0.
\end{equation*}
$$

This means that 

$$
\begin{equation*}
  \frac{d\phi}{dt} = -\frac{g_{t\phi}}{g_{\phi \phi}} \pm \sqrt{\frac{g_{t\phi}^2}{g_{\phi \phi}^2} - \frac{g_{tt}}{g_{\phi \phi}}}
\end{equation*}
$$

Consider this at the surface where $\Sigma = 2r$, meaning $g_{tt} = 0$. Then the two solutions are

$$
\begin{equation}
  \frac{d\phi}{dt} = 0, \quad \frac{d\phi}{dt} = - \frac{2g_{t\phi}}{g_{\phi \phi}} = \frac{a}{2 + a^2}
\end{equation}
$$

The second solution is just a photon which moves along with the rotation of the black hole. The other solution is a photon which is instantaneously at rest in the BL coordinates. This is to say that at this boundary where $g_{tt} = 0$ no particles can move against the rotation of the black hole. This is also true for the region below this boundary, which is called the *ergosphere*. 

So this surface where $g_{tt} = 0$ seems special, in that here photons can be instantaneously at rest in the BL coordinates. It is therefore interesting to find out how this surface looks. For $g_{tt}$ to be zero, we need $\Sigma = 2r$, meaning 

$$
\begin{align*}
    &r^2 + a^2 \cos^2\theta - 2r = 0 \\ \\
    \implies &r = 1 \pm \sqrt{1 - a^2 \cos^2\theta}.
\end{align*}
$$

Clearly there are two solutions, so there is both an inner and outer ergosphere boundary just like with the inner and outer event horizons. The following figure illustrates how the event horizons and ergosphere boundaries look with a slice of these surfaces in the $xz$-plane as we vary $a$ from 0 to 1. 

![ergosphere-and-EH](../images/miscellaneous/ergosphere_and_EH.gif)

We can also show the effects of frame dragging, which is done in Figure \ref{fig: gridlines}. Here the inward pointing curves represent photon geodesics with $L = 0$. And the dashed red circle represents the boundary of the ergosphere, while the filled black circle represents the event horizon of the black hole. Clearly the photon geodesics are forced to move with the rotation of the black hole. 

![Gridlines](../images/miscellaneous/gridlines.png)

Frame dragging will be important because, while observers stationary in the BL coordinates might feel quite natural to work with, they no longer exist at all points outside the event horizon in the Kerr spacetime, in particular it is impossible inside the ergosphere. The simple generalization we instead consider is to rather work in ZAMO frames, because there exist local ZAMO frames even inside the ergosphere. 

### 5.4 Rest frame to ZAMO to global transformation

When we emit rays from our camera the components of the four-vectors describing these rays are initialized in the rest frame of the camera. But when solving the geodesic equation we are working in the global BL frame. So we need a way to transform four-vectors from the camera's rest frame to the global frame. This is where the concept of *tetrads* becomes relevant. 

#### 5.4.1 Tetrads

Tetrads, also called *frame fields*, is a set of orthonormal vectors which form a local inertial frame. So the tetrads take us from the global coordinate system in which the metric is $g_{\mu \nu}$ to the local one in which the metric is that of Minkowski $\eta_{\mu \nu}$. The tetrads $e_m^{\;\mu}$ are therefore defined by

$$
    g_{\mu \nu} e_m^{\; \mu} e_n^{\; \nu} = \eta_{mn}.
$$

Keep in mind that Latin indices are here used to refer to the indices in the local inertial frame. We also have the dual tetrads $e^{*m}_\mu$, also called the co-frame field, defined through

$$
\begin{equation}
  e_m^{\mu} e^{*m}_{\nu} = \delta_\nu^\mu, \quad e_m^{\mu} e^{*n}_{\mu} = \delta_m^n.
\end{equation}
$$

From (gotta find that paper again) and [here](https://arxiv.org/abs/1506.01473v2) we have that the frame fields for the local ZAMO frame are given by

$$
\begin{align*}
  e_t^{\mu} &= \delta_t^\mu\sqrt{\frac{\Lambda}{\Delta \Sigma}} + \delta_\phi^\mu \frac{2ar}{\sqrt{\Lambda \Delta \Sigma}} \label{eq: vierbein_tmu} \hspace{1.1cm}\\ \\
  e_r^{\mu} &= \delta_r^\mu \sqrt{\frac{\Delta}{\Sigma}} \label{eq: vierbein_rmu}\\ \\
  e_\theta^{\mu} &= \delta_\theta^\mu \frac{1}{\sqrt{\Sigma}} \label{eq: vierbein_thetamu}\\ \\
  e_\phi^{\mu} &= \delta_\phi^\mu \sqrt{\frac{\Sigma}{\Lambda}}\frac{1}{\sin \theta}. \label{eq: vierbein_phimu}
\end{align*}
$$

while the co-frame fields are given by 

$$
\begin{align*}
  e^{*t}_{\mu} &= \delta_\mu^t \sqrt{\frac{\Delta \Sigma}{\Lambda}} \label{eq: vierbein^tmu}\\ \\
  e^{*r}_{\mu} &= \delta_\mu^r \sqrt{\frac{\Sigma}{\Delta}} \label{eq: vierbein^rmu} \\ \\
  e^{*\theta}_{\mu} &= \delta_\mu^\theta \sqrt{\Sigma} \label{eq: vierbein^thetamu}\\ \\
  e^{*\phi}_{\mu} &= -\delta_\mu^t \frac{2ar\sin\theta}{\sqrt{\Lambda \Sigma}} + \delta_\mu^\phi \sin\theta \sqrt{\frac{\Lambda}{\Sigma}}. \label{eq: vierbein^phimu}
\end{align*}
$$

To transform a vector $x^m$ from the local ZAMO frame to the global frame we then do the following

$$
\begin{equation}
    x'^\mu = e_m^{\mu} x^m.
\end{equation}
$$

And the inverse naturally follows

$$
\begin{equation}
    x^m = e^{*m}_{\mu} x'^\mu.
\end{equation}
$$

#### 5.4.2 Relativistic aberration

The steps 

1. Compute the local velocity of the camera in the local ZAMO frame.
2. Lorentz transform the initial photon four-momentum from the instantaneous rest frame of the camera to the local ZAMO frame.
3. Transform the four-momentum of the photon from hte local ZAMO frame to the global frame using the ZAMO tetrads.

___
##### Step 1

The camera's four-velocity is will primarily be expressed in the global frame. Let's call the four-velocity in that frame $u^\mu$. Imagine then a ZAMO frame at the same position as the camera. The four-velocity of the camera in the ZAMO tetrad frame is given by

$$
\begin{equation}
    \tilde{u}^m = e^m_{\mu} u^\mu
\end{equation}
$$

Since this four-velocity is measured in a local inertial frame it can also be expressed as 

$$
\begin{equation}
    \tilde{u}^m = \gamma \left(1, \vec{v}\right),
\end{equation}
$$

where $\gamma = \frac{1}{\sqrt{1 - \vec{v}^2}}$ and $\vec{v}$ is the local velocity of the camera in the ZAMO frame. Thne $\gamma = \tilde{u}^t$, and so the components $v^i$ of the local velocity $\vec{v}$ of the camera in the ZAMO frame are given by 

$$
\begin{equation}
    v^i = \frac{\tilde{u}^i}{\gamma} = \frac{\tilde{u}^i}{\tilde{u}^t}
\end{equation}
$$

That concludes step 1. 

___ 
##### Step 2

Consider now the rays that we emit from our camera. We now need to transform these to the local ZAMO frame. This is just a Lorentz transformation, whose general form is the following (see for example [2])

$$
\begin{equation}
  \Lambda^m_{\;\;n} = \left[\begin{matrix}\gamma && -\gamma v^x && -\gamma v^y && -\gamma v^z\ \\
                          -\gamma v^x && 1 + \left(\gamma - 1\right)\frac{\left(v^x\right)^2}{\vec{v}^2} && \left(\gamma - 1\right) \frac{v^x v^y}{\vec{v}^2} && \left(\gamma - 1\right)\frac{v^x v^z}{\vec{v}^2} \\
                          -\gamma v^y && \left(\gamma - 1\right)\frac{v^x v^y}{\vec{v}^2} && 1 + \left(\gamma - 1\right)\frac{\left(v^y\right)^2}{\vec{v}^2} && \left(\gamma - 1\right)\frac{v^y v^z}{\vec{v}^2} \\
                          -\gamma v^z && \left(\gamma - 1\right)\frac{v^x v^z}{\vec{v}^2} && \left(\gamma - 1\right)\frac{v^y v^z}{\vec{v}^2} && 1 + \left(\gamma - 1\right)\frac{\left(v^z\right)^2}{\vec{v}^2}
                         \end{matrix}\right]
\end{equation}
$$

in Cartesian coordinates. But the velocit components $v^i$ we have are in BL coordinates. So we first have to convert the velocity vector $\vec{v}$ from BL to Cartesian coordinates. We know that the components of the velocity vector in Cartesian coordinates are 

$$
\begin{align*}
  v^x &= \frac{1}{\sqrt{r^2 + a^2\cos^2\theta}}\left[r\sin\theta \cos\phi \,v^r + \sqrt{r^2 + a^2}\cos\theta \cos\phi \,v^\theta\right] - \sin\phi \,v^\phi \\ \\
  v^y &= \frac{1}{\sqrt{r^2 + a^2 \cos^2\theta}}\left[r\sin\theta \sin\phi \,v^r + \sqrt{r^2 + a^2}\cos\theta \sin\phi \,v^\theta\right] + \cos\phi \,v^\phi \\ \\
  v^z &= \frac{1}{\sqrt{r^2 + a^2 \cos^2\theta}}\left[\sqrt{r^2 + a^2}\cos\theta \,v^r - r\sin\theta \,v^\theta\right].
\end{align*}
$$

when the components are $v^r$, $v^\theta$ and $v^\phi$ in BL coordinates. Let's call the four-momentum of the initial photon expressed in Cartesian coordinates in the camera rest frame $\tilde{p}'^m$. In the ZAMO frame the four-momentum is then given by the unprimed 

$$
\begin{equation}
    \label{eq: camera_four_mom_ZAMO}
    \tilde{p}^m = \Lambda^m_n\left(\vec{v}\right) \tilde{p}'^n
\end{equation}
$$

Then we move back into BL coordinates. And the coordinate transformation is then simply given by the matrix in $\eqref{eq: unit_vector_relation}$. 

___

##### Step 3

The final step is then rather simple. The four-momentum of the ray in the global frame is simply 

$$
\begin{equation}
    p^\mu = e_m^{\; \mu} \tilde{p}^m.
\end{equation}
$$

___

### 5.5 Free-falling camera

#### 5.5.1 Examples

___

## 6 Redshift 

Consider an observer at some point $P$ in spacetime with four-velocity $u^\mu$ in some coordinate system. Assume we also have a photon with four-momentum $p^\mu = \frac{dx^\mu}{d\lambda}$ at the same point $P$. The frequency that the observer measures for the photon at $P$ is 

$$
\begin{equation}
  \label{eq: frequency_definition}
  \nu = -p^\mu u_\mu = -g_{\mu \rho} p^\mu u^\rho,
\end{equation}
$$

which is of course coordinate invariant. Here redshift is of interest because we have an "observer" at the accretion disk around the black hole that emits light with certain frequencies, and that light is then observed by the camera at another point in spacetime (Of course in our simulation it's the other way around, but in reality the photon would be emitted by the accretion disk and be observed by the camera). And since we have solved the geodesic equation for the rays we already know the four-momentum of the photons at all points along their geodesics. In particular we know their four-momenta at the disk and at the camera. Consider a point on the disk with coordinates $(r, \theta, \phi)$ and let the camera be at $(R, \Theta, \Phi)$. Let the parameter $\lambda$ equal 0 at the camera and $\lambda_1 > 0$ at the point on the disk. We will call the four-velocity of the relevant point on the disk $u^\mu_\mathrm{disk}$ and the four-velocity of teh camera $u^\mu_\mathrm{camera}$. Then the redshift $(1 + z)$ for the light emitted at the disk and observed by our camera is given by 

$$
\begin{equation}
  \label{eq: redshift_factor} 
  1 + z = \frac{\lambda_\mathrm{camera}}{\lambda_\mathrm{disk}} = \frac{\nu_\mathrm{disk}}{\nu_\mathrm{camera}} = \frac{g_{\mu \nu}\left(r, \theta, \phi\right) p^\mu\left(\lambda = \lambda_1\right) u^\nu_\mathrm{disk}}{g_{\mu \nu}\left(R, \Theta, \Phi\right) p^\mu \left(\lambda = 0\right) u^\nu_\mathrm{camera}}
\end{equation}
$$

where $\nu_\mathrm{disk}$ and $\nu_\mathrm{camera}$ are the frequencies measured at the disk and camera respectively. And $\lambda_\mathrm{disk}$ and $\lambda_\mathrm{camera}$ are the measured wavelengths at the disk and camera. This general expression has the advantage that it doesn't assume anything about the motion of neither the disk nor the camera. So this places no restrictions on how we can move our camera around (apart from the local velocity not exceeding $c$ of course) or what kind of velocity distribution the disk can have. Also note that by "redshift" we really mean both redshift and blueshift. We are just referring to both decrease and increase in wavelength under the same umbrella term for simplicity (this is really just a habit from cosmology, where we only deal with redshift). It is also important to note that $\eqref{eq: redshift_factor}$ accounts for both the Doppler-like redshift caused by motion and the gravitational redshift. Both of these two kinds of redshift are baked into the metric itself. 

### 6.1 Why does it happen

### 6.2 Gravitational & Doppler redshift

### 6.3 Velocity distribution for the accretion disk

The accretion disk is in reality moving in a spiral towards the black hole. But when dealing with the physics here we will simplify and assume all points on the disk move in circular orbits. From [3] we know that the four-velocity of a massive particle in a prograde, circular orbit in the equatorial plane around a Kerr black hole is given by $u^\mu = \left(u^t, 0, 0, u^\phi\right)$, where 

$$
\begin{equation}
  u^\phi = \frac{\sqrt{Mr}}{r\sqrt{r^2 - 3Mr + 2\lvert a\rvert \sqrt{Mr}}}.\tag{53}
\end{equation}
$$

And $u^t$ is of course determined by the condition $u^\mu u_\mu = -1$. 
___

## 7 Colors

Determining the color of the pixels in the image is of course crucial for the visualizations, and there are multiple ways of doing this. We could just say that each point on the disk and in the jet emits radiation at a certain wavelength. And then we could redshift this according to the description above to get an idea of how light from different parts of the scene is redshifted. That can look something like this 

[Render of single-color emission]

But this is not very realistic. The disk and jet will emit light at all wavelengths, but the particular distribution of light over the wavelengths can be very complicated. Therefore, in order to perform a simplification here, we will assume that the disk and jet act as *blackbodies*. For our purposes all this means is that the distribution of light over the wavelengths is given by the *Planck distribution*, also called the *blackbody distribution*

$$
\begin{equation}
  I_{\lambda}\left(\lambda; T\right) = \frac{2h c^2}{\lambda^5} \frac{1}{e^{hc/\lambda k_B T} - 1},\tag{60}
\end{equation} 
$$

where $h$ and $k_B$ are Planck's constant and Boltzmann's constant respectively. $T$ is the temperature of a particular point on the disk/jet and $\lambda$ is wavelength. This is an *intensity* distribution, and describes how much light is emitted at each wavelength. Now we need a way to convert an intensity distribution to an RGB color that the screen can output. 

The relative amounts of the visible wavelengths that are emitted determines the color of the object. Unfortunately our eyes do not perceive and respond to all visible wavelengths of light equally. So converting an intensity distribution $I_\lambda$, which represents the light that is *actually* observed, to the color our eyes *see* is far from trivial. The procedure for calculating an RGB color comes in two main steps: converting a spectrum to so-called *tristimulus values* $X, Y$ and $Z$, and then converting the tristimulus values to RGB values. We found [5] [this](https://color.org/chardata/rgb/sRGB.pdf) and [this](https://en.wikipedia.org/wiki/SRGB) to be particularly helpful here. 

Assume now that we have a blackbody at some temperature $T$. We can find the tristimulus values $X, Y$ and $Z$ using the so-called *color matching functions* $\bar{x}, \bar{y}$ and $\bar{z}$, which describe how our eyes respond to different wavelengths. The $X, Y$ and $Z$ values are calculated in the following way

$$
\begin{align*}
  X &= \int_{\mathrm{380 nm}}^{\mathrm{780 nm}} I_\lambda \left(\lambda; T\right) \bar{x}\left(\lambda\right)d\lambda \tag{61} \\ \\
  Y &= \int_{\mathrm{380 nm}}^{\mathrm{780 nm}} I_\lambda \left(\lambda; T\right) \bar{y}\left(\lambda\right)d\lambda \tag{62}\\ \\
  Z &= \int_{\mathrm{380 nm}}^{\mathrm{780 nm}} I_\lambda \left(\lambda; T\right) \bar{z}\left(\lambda\right)d\lambda. \tag{63}
\end{align*}
$$

Approximations to these color matching functions can be found in [6]. But it turns out that at low temperatures these approximations are not accurate enough anymore. So we will instead be using the lookup table given in [7]. In the following figure we have plotted the color matching functions as functions of wavelength. The $x$-axis represents wavelength $\lambda$ in nm, while the $y$-axis represents the output of the color matching functions, which can be taken to have units such that $X, Y$ and $Z$ are dimensionless.

![color-matching-functions](../images/miscellaneous/color_matching_functions.png)

While we have color-coded the different color matching functions, the three values $X, Y$ and $Z$ do not directly correspond to red, green and blue. Instead $X, Y$ and $Z$ live in an abstract color space. 

### 7.1 The chromaticity diagram

When we talk about colors we tend to group them into "classes" of colors that are similar, like greens, reds blues etc. And we have an intuition that some of those colors are really the same basic color, just at different brightnesses, or more precisely *luminances*. And the quantity that is the same between different *luminances* is called the *chromaticity*. 

To explain this in more detail we can study the *chromaticity diagram*. This is a 2d plot which shows the chromaticities corresponding to different points in $XYZ$ space. And the color matching functions give us a map between intensity distributions/spectra and the $XYZ$ space. Consider now a spectrum consisting entirely of a single wavelength $\lambda_0$ - which is to say monochromatic light. We can represent this mathematically by $I_\lambda\left(\lambda\right) = \delta\left(\lambda - \lambda_0\right)$ (nevermind the units or the scaling), where $\delta\left(x - y\right)$ is the Dirac-delta "function". Plugging this into (61) - (63) gives 

$$
\begin{align*}
  X &= \bar{x}\left(\lambda_0\right) \tag{64}\\ \\
  Y &= \bar{y}\left(\lambda_0\right) \tag{65}\\ \\
  Z &= \bar{z}\left(\lambda_0\right).\tag{66}
\end{align*}
$$

If we let $\lambda_0$ vary over the visible range we can trace out the curve in $XYZ$ space made by the wavelengths in the visible spectrum. But in order to visualize this space in 2d we define the new coordinates $x$ and $y$ through

$$
\begin{align*}
  x &= \frac{X}{X + Y + Z} \tag{67}\\ \\
  y &= \frac{Y}{X + Y + Z} \tag{68}.
\end{align*}
$$

It also naturally follows that we can define $z = \frac{Z}{X + Y + Z} = 1 - x - y$. We can now plot the $(x, y)$ coordinates of the monochromatic spectra, and the resulting curve is called the *spectral locus*. This is shown in the following figure. The marked points show where different wavelengths fall on the spectral locus. 

![chromaticity-diagram-without-planckian-locus](../images/miscellaneous/chromaticity_diagram_without_planckian_locus.png)

This spectral locus is the boundary of the chromaticity diagram mentioned earlier, and which is shown in [this](https://en.wikipedia.org/wiki/CIE_1931_color_space) Wikipedia article. Here we have colored in the different points along the spectral locus according to which RGB colors they correspond to. Of course we have not yet explained how we can compute these RGB colors from the $XYZ$ values, but we get to that later. 

We can also calculate the curve in the chromaticity diagram which represents the color of blackbodies at different temperatures. This just involves calculating $X, Y$and $Z$ from the spectrum $I_\lambda\left(\lambda; T\right)$ for varying $T$. The resulting curve in $(x, y)$ coordinates is called the *Planckian locus* (or equivalently the *blackbody locus*). Plotted together with the spectral locus the result we get is shown in the following figure. Here the thick curve represents the Planckian locus. We have again calculated the RGB color for each point along the Planckian locus. 

### 7.2 Transformation from XYZ to RGB 

Now we can discuss how you actually get from $X, Y$ and $Z$ to RGB values. There are many different RGB color spaces designed for outputting colors to a screen, and these different RGB color spaces are defined by their so-called *primaries*. When choosing an RGB color space we choose which primaries, red, green and blue to use. These primaries correspond to the vertices in the colored triangle in the following figure.  

![sRGB-gamut](../images/miscellaneous/SRGB_chromaticity_CIE1931.svg){:width="500px"}

The set of points inside the triangle between the primaries is called the *gamut* of the RGB space, and makes up the colors that the given RGB space can represent. And since different RGB spaces have different primaries it follows that the sets of colors that they can represent are also different. 

The points inside the gamut are linear combinations of the primaries, so the primaries form a basis for the RGB space. Red $R$, green $G$ and blue $B$ being the primaries are of course written as 

$$
\begin{equation}
  R = \left[\begin{matrix}
              1 \\
              0 \\
              0
             \end{matrix}\right], \quad
  G = \left[\begin{matrix}
                0 \\
                1 \\
                0
             \end{matrix}\right], \quad
  B = \left[\begin{matrix}
              0 \\
              0 \\
              1
             \end{matrix}\right].\tag{69}
\end{equation}
$$

in the RGB basis. We also have an $XYZ$ basis, in which the primaries can be written as 

$$
\begin{equation}
  R = \left[\begin{matrix}
              X_\mathrm{R} \\
              Y_\mathrm{R} \\
              Z_\mathrm{R}
             \end{matrix}\right], \quad
  G = \left[\begin{matrix}
                X_\mathrm{G} \\
                Y_\mathrm{G} \\
                Z_\mathrm{G}
             \end{matrix}\right], \quad
  B = \left[\begin{matrix}
              X_\mathrm{B} \\
              Y_\mathrm{B} \\
              Z_\mathrm{B}
             \end{matrix}\right].\tag{70}
\end{equation}
$$

The problem of finding the RGB values fro some given set of $XYZ$ values then translates into finding the coordinate transformation between the RGB coordinates and the $XYZ$ coordinates. We can write the transformation as 

$$
\begin{equation}
  \left[\begin{matrix}
          R \\
          G \\
          B
        \end{matrix}\right] = M\left[\begin{matrix}
                                      X \\
                                      Y \\
                                      Z
                                     \end{matrix}\right]\tag{71}
\end{equation}
$$

where $M$ is the transformation matrix we want to find. The standard way of specifying the primaries of an RGB color space is by giving the $x, y$ and $Y$ values of the primaries. So in order to work with the $X, Y$ and $Z$ values we first have to calculate the unknown $X$ and $Z$ values. From (67) we see that 

$$
\begin{equation}
    X = \left(X + Y + Z\right)x = \frac{x}{y}Y.\tag{72}
\end{equation}
$$

And likewise $Z = \frac{z}{y}Y$. Now we want to find the linear transformation $M$ relating the primaries expressed in the RGB basis to the primaries expressed in the XYZ basis, that is

$$
\begin{equation}
  \left[\begin{matrix}
          1 & 0 & 0 \\
          0 & 1 & 0 \\
          0 & 0 & 1
        \end{matrix}\right] = M\left[\begin{matrix}
                                        X_\mathrm{R} & X_\mathrm{G} & X_\mathrm{B} \\
                                        Y_\mathrm{R} & Y_\mathrm{G} & Y_\mathrm{B} \\
                                        Z_\mathrm{R} & Z_\mathrm{G} & Z_\mathrm{B}
                                     \end{matrix}\right] \equiv MA\tag{73}
\end{equation}
$$

where we have defined the matrix $A$ as the matrix containing the $XYZ$ coordinates of the primaries. Then the coordinate transformation matrix $M$ from earlier is just $A^{-1}$. We have of course just assumed that $A$ has an inverse, but we are physicists here, not mathematicians. 

In order to actually compute the values of the entries in $M$ we need to know the $x, y$ and $Y$ values of the primaries. sRGB, which is the color space we will be working with, has the following primaries

$$
\begin{equation}
  \left[\begin{matrix}
          x_\mathrm{R} \\
          y_\mathrm{R} \\
          Y_\mathrm{R}
        \end{matrix}\right] = \left[\begin{matrix}
                                      0.64 \\
                                      0.33 \\
                                      0.2126
                                    \end{matrix}\right], \quad\quad
  \left[\begin{matrix}
          x_\mathrm{G} \\
          y_\mathrm{G} \\
          Y_\mathrm{G}
        \end{matrix}\right] = \left[\begin{matrix}
                                      0.3 \\
                                      0.6 \\
                                      0.7152
                                    \end{matrix}\right], \quad\quad
  \left[\begin{matrix}
          x_\mathrm{B} \\
          y_\mathrm{B} \\
          Y_\mathrm{B}
        \end{matrix}\right] = \left[\begin{matrix}
                                      0.15 \\
                                      0.06 \\
                                      0.0722
                                    \end{matrix}\right].\tag{74}
\end{equation}
$$

With these values we get

$$
\begin{equation}
  M = \left[\begin{matrix}
              3.24156456 & -1.53766524 & -0.49870224 \\
              -0.96920119 & 1.87588535 & 0.04155324 \\
              0.05562416 & -0.20395525 & 1.05685902
            \end{matrix}\right]\tag{75}
\end{equation}
$$

It is also typically advised to perform a final non-linear correction after this linear transformation. But our results looked better without this non-linear correction, so we chose to omit it. 

### 7.3 Redshifted colors

Remember that the reason we wanted to find the transformation from $XYZ$ to RGB was that we know how to calculate the $XYZ$ values corresponding to a given spectrum. So with this coordinate transformation in hand we can convert a spectrum into an RGB color. This is of interest to us here because we know the temperature $T$ and redshift $(1 + z)$ of points on the accretion disk and jet. So we can calculate the blackbody spectrum of each point and convert that to an RGB color to display on screen. Incorporating redshift into this is actually very easy. Consider a blackbody spectrum redshifted such that $\lambda_\text{shifted} = \left(1 + z\right)\lambda$. Then 

$$
\begin{align*}
    I_{\lambda}\left(\lambda_\text{shifted}; T\right) &= \frac{2h c^2 \left(1 + z\right)^5}{\lambda_\text{shifted}^5} \frac{1}{e^{hc \left(1 + z\right)/\lambda_\text{shifted} k_B T} - 1} \\ \\
    &= \frac{2hc^2 \left(1 + z\right)^5}{\lambda_\text{shifted}^5} \frac{1}{e^{hc/\lambda_\text{shifted} k_B T_\text{shifted}} - 1}, \tag{76}
\end{align*}
$$

with $T_\text{shifted} \equiv \frac{T}{1 + z}$. Notice now that this is also a blackbody spectrum, evaluated at a "redshifted temperature". And we have an additional vertical scaling by $(1 + z)^5$. From earlier we know that relativistic beaming modifies the observed intensity by a factor $(1 + z)^{-5}$, and when applying this to (76) these factors of $(1 + z)$ cancel out. This means that we can actually just use the unmodified blackbody spectrum $I_\lambda$, just with a temperature $T_\mathrm{shifted} = \frac{T}{1 + z}$. So this accounts for relativistic beaming. 

We can plot the map showing the blackbody color for different combinations of $T$ and $1 + z$. This is shown in the following figure for temperatures $T \in \left[200, 10000\right]$ K and for $(1 + z) \in \left[0.1, 2\right]$.

![temp-redshift-map](../images/miscellaneous/temp_redshift_map.png)

### 7.6 Example 1: Redshift

[Low temperature disk and jet with clear redshift and beaming turned off, no camera motion]

### 7.7 Example 2: Relativistic beaming

[Same as above]
___

## 8 Review

### 8.1 What we have so far

### 8.2 Why it's about to get a lot worse

___

## 9 The volumetric accretion disk

### 9.1 Disk models

### 9.2 Basic volume rendering 

### 9.3 Challenges

### 9.4 Procedural voronoi based disk noise

### 9.5 Examples

### 9.6 Temperature distributions

### 9.7 Examples

___

## 10 Light travel delay

One aspect of the visualizations which might not be obvious is that we have to take into account how long it takes light to travel from each part of the disk to the camera. Since it takes different amounts of time for light from different parts of the disk to reach us the photons that we receive at some time $t$ must have been emitted from the different parts of the disk at different times. So if we for example take our image at global time $t$ and it takes a time $t_\mathrm{A}$ for light from a point A on the disk to reach us, then the light that reaches us at time $t$ must have left A at time $t - t_\mathrm{A}$. So when the travel time $t_\mathrm{A}$ differs between the points on the disk the parts of the disk we see will correspond to different emission times. Let us assume we have a function $f(\vec{x}, t)$ which calculates the density at the position $\vec{x}$ on the disk evaluated at global time $t$. Then the density that the camera actually observes at point $\vec{x}_i$ when it takes an image at time $t$ is $f(\vec{x}_i, t - t_i)$, where $t$ is the time when the image is taken and $t_i$ is the time it takes light to travel from point $\vec{x}_i$ on the disk to the camera. 


### 10.2 Simple non-relativistic example

### 10.3 Relativistic example

### 10.4 Disk delay angle

A convenient aspect of the geodesic tracing here is that we already calculate this travel time in the global coordinate system when we solve the geodesic equation. And since we assume that each point on the disk moves along a circular orbit this correction due to the travel time only requires a correction for the angle that a given point on the disk moves during the travel time. Furthermore, since the four-velocity of the points on the disk is constant this angular correction is particularly simple:

$$
\begin{equation}
    \Delta \phi = \int_0^{t_\mathrm{travel}} \frac{d\phi}{dt'} dt' = \int_0^{t_\mathrm{travel}} \frac{d\phi}{d\tau} \frac{d\tau}{dt'} dt' = \frac{u^\phi}{u^t}t_\mathrm{travel},\tag{59}
\end{equation}
$$

### 10.5 Example

[Render with clear effects of travel time delay]
___

## 11 The astrophysical jet

### 11.1 Jet models

### 11.2 Procedural jet model with special coordinates

### 11.3 Velocity distribution

### 11.4 Temperature distribution

### 11.5 Examples

### 11.6 Why does it look so weird

### 11.7 Delaying the jet

### 11.8 Examples

### 11.9 Why does it look even weirder

___

## 12 Improving the renders

### 12.1 Super sampling

### 12.2 Anti aliasing

### 12.3 Physical motion blur

### 12.4 FFT Bloom

___

## 13 Test cases

### 13.1 Example 1

### 13.2 Example 2 

### 13.3 Example 3

### 13.4 Contemporary example comparisons

___

## 14 Conclusion

### 14.1 Limitations

### 14.2 Improvements

### 14.3 Special thanks

___

## 15 Gallery

Sample gallery:

![Absolutely godlike](../images/black_hole_renders/below_pretty.png)
![Close up](../images/black_hole_renders/close_up.png)
![Hot disk](../images/black_hole_renders/hot_disk_amazing.png)
![Hot jet](../images/black_hole_renders/hot_jet.png)
![Clear redshift](../images/black_hole_renders/Render3.png)
![Aberration](../images/black_hole_renders/aberration.png)
![R1](../images/black_hole_renders/R1.png)
![Far away](../images/black_hole_renders/MBRender2.png)
![Shadow](../images/black_hole_renders/flat_edge.png)
![Epic jet](../images/black_hole_renders/epic_jet.png)

## 16 References 

[1] Sean M. Carroll. Spacetime and Geometry: An Introduction to General Relativity. Cambridge University Press, 2019

[2] Masud Mansuripur. An exact derivation of the thomas precession rate using the lorentz transformation. In Henri-Jean M. Drouhin, Jean-Eric Wegrowe, and Manijeh Razeghi, editors, Spintronics XIII. SPIE, August 2020

[3] Z. Kh. Kurmakaev. Circular orbits in the Kerr metric. , 18:110, August 1974

[4] Eric Bruneton. Real-time high-quality rendering of non-rotating black holes, 2020

[5] Philippe Colantoni. Color space transformations. 2006

[6] Chris Wyman et al. Simple analytic approximations to the cie xyz color matching  functions. 2013

[7] Kasajima Ichiro. Plotting colors on color circle: Interconversion between xyz values and rgb color system. Current Trends in Analytical and Bioanalytical Chemistry, 1, 05 2017