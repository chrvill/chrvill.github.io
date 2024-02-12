---
title: Rendering a rotating black hole
date: 2023-12-28
categories: [physics, general relativity, black holes]
tags: [black holes]
math: true
---

![alt text](../images/miscellaneous/geodesics.gif)

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
    \label{eq: Kerr_metric} \tag{1}
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
    x &= \sqrt{r^2 + a^2}\sin\theta\cos\phi, \label{eq: x_BL} \tag{2} \\
    y &= \sqrt{r^2 + a^2}\sin\theta\sin\phi, \label{eq: y_BL} \tag{3} \\
    z &= r\cos\theta. \label{eq: z_BL} \tag{4}
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

#### Coordinate transformation between Cartesian and Boyer-Lindquist coordinates

We will need to convert vectors between Cartesian and BL coordinates multiple times, so a coordinate transformation between these systems is necessary. Deriving this is most easily done by studying an infinitesimal displacement vector $d\vec{r}$ expressed in the two coordinate systems. This can be written as 

$$
\begin{equation}
  \label{eq: infinitesimal_displacement_Cartesian}
  d\vec{r} = dx \; \vec{\hat{x}} + dy \; \vec{\hat{y}} + dz \; \vec{\hat{z}}
\end{equation}
$$

in Cartesian, and 

$$
\begin{equation}
  \label{eq: infinitesimal_displacement_spheroidal}
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
  \label{eq: r_unit_vector}
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
  \label{eq: unit_vector_relation}
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

### 2.3 Geodesics
___

## 3 Tracing geodesics in curved spacetime

#### 3.1.2 The geodesic equation 

In order to visualize the black hole we need to trace *geodesics* in the Kerr spacetime, which we do by solving the *geodesic equation*

$$
\begin{equation}
  \label{eq: geodesic_equation} \tag{5}
  \frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\rho \sigma}\frac{dx^\rho}{d\lambda}\frac{dx^\sigma}{d\lambda} = 0,
\end{equation}
$$

where $x^\mu$ is the four-position of a photon, and $\lambda$ is an affine parameter for the geodesic, chosen such that $p^\mu \equiv \frac{dx^\mu}{d\lambda}$ is the four-momentum of the photon. $\Gamma^\mu_{\rho \sigma}$ are the *Christoffel symbols* of the metric. Equation $\eqref{eq: geodesic_equation}$ is completely general, and we of course need explicit expressions for the geodesic equation for each value of $\mu$. But computing the Christoffel symbols tends to be very tedious work, and certainly so for the Kerr metric. So instead of deriving the Christoffel symbols by hand we use the Sympy package in Python to derive them, and thus also explicit expressions for each component of the geodesic equation. Code for doing this can be found [here](https://github.com/chrvill/Geodesic_EOM_deriver).

We have the normalization requirement

$$
\begin{equation}
    \label{eq: normalization} \tag{6}
    p_\nu p^\nu = \mu
\end{equation}
$$

with $\mu = 0$ for massless particles like photons and $\mu = -1$ for massive particles (This is really the normalization of the four-*velocity* of a massive particle, not the four-momentum. But the difference is only a factor of $m^2$, with $m$ being the mass of the particle in question.). This can be written explicitly as 

$$
\begin{equation}
    \label{eq: normalization_explicit} \tag{7}
    g_{tt} \left(p^t\right)^2 + g_{rr}\left(p^r\right)^2 + g_{\theta\theta} \left(p^\theta\right)^2 + g_{\phi \phi}\left(p^\phi\right)^2 + 2g_{t\phi} p^t p^\phi = \mu
\end{equation}
$$

When we choose initial directions for our rays we will essentially provide the spatial components of the four-momentum, so we need to compute $p^t$ from $p^r, p^\theta$ and $p^\phi$. Solving \eqref{eq: normalization_explicit} for $p^t$ gives

$$
\begin{equation}
    \label{eq: p^t expression} \tag{8}
    p^t = -\frac{g_{t\phi}}{g_{tt}} p^\phi \pm \sqrt{\left(\frac{g_{t\phi}}{g_{tt}} p^\phi\right)^2 - \frac{1}{g_{tt}}\left(g_{rr} \left(p^r\right)^2 + g_{\theta\theta} \left(p^\theta\right)^2 + g_{\phi \phi} \left(p^\phi\right)^2 - \mu\right)}
\end{equation}
$$

But we have to find out which sign is appropriate. Consider a case where $p^\phi = 0$, which is only possible outside the ergosphere. Then, since $p^t$ needs to be positive we need to use the $+$ sign outside the ergosphere. Arguing for which sign is appropriate inside the ergosphere proved trickier, but we found that we had to use the $-$ sign inside the ergosphere in order for $p^t$ to increase continuously as you pass the boundary of the ergosphere. 

### 3.2 Initial conditions

### 3.3 Numerically integrating the equations of motion

Numerically \eqref{eq: geodesic_equation} is just a completely standard second order differential equation. So in that sense solving it is just like solving any other second order differential equation numerically. However it should be noted that some of the Christoffel symbols will diverge as the photon approaches the event horizon, which needs to be taken into consideration when solving the geodesic equation. Assume for the purposes of this discussion that we have a first order differential equation of the form

$$
\begin{equation}
    \label{eq: general_diff_eq} \tag{9}
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

### 5.1 Basic frame of reference transformation

### 5.3 The ZAMO and ergosphere

One of the Killing vectors of the Kerr metric is the following

$$
\begin{equation}
    R^\mu = \left(\partial_\phi\right)^\mu = \left(0, 0, 0, 1\right)
\end{equation}
$$

This is related to rotational symmetry, and we can therefore identify the following conserved quantity

$$
\begin{equation}
    L \equiv g_{\mu \nu} R^\mu p^\nu = g_{\phi \nu} p^\nu = g_{\phi \phi} p^\phi + g_{t\phi} p^t
\end{equation}
$$

along a geodesic. Consider now the expression for $L$ at $\theta = \pi/2$ as $r \to \infty$. Then $g_{t\phi} \to 0$ and $g_{\phi \phi} \to r^2$. So $L \to r^2 p^\phi$. And $p^\phi$ reduces to an angular velocity in flat space. So this is the standard definition of angular momentum (per unit mass) in flat space. Meaning that the conserved quantity $L$ corresponds to the angular momentum of an object infinitely far away, in the asymptotically flat region of the spacetime. We can now consider releasing a particle from infinitely far away with $L = 0$ and letting it fall along a geodesic. Then we have $g_{t\phi} p^t + g_{\phi \phi} p^\phi = 0$, meaning that 

$$
\begin{equation}
  \label{eq: frame_dragging_p_phi}
  p^\phi = -\frac{g_{t\phi}}{g_{\phi \phi}}p^t,
\end{equation}
$$

which is to say that the particle gains an angular velocity even though it starts off falling with zero angular momentum. Observers following geodesics with $L = 0$ are called \textit{Zero Angular Momentum Observers} (ZAMO). The phenomenon causing ZAMO frames to move along with the rotation of the black hole is called \textit{frame dragging}. It turns out that frame dragging has quite dramatic consequences for the motion of particles around a rotating black hole. Consider a photon in the equatorial plane moving only in the $\phi$ direction. Then

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
  \label{eq: frame_dragging_angular_velocity}
  \frac{d\phi}{dt} = 0, \quad \frac{d\phi}{dt} = - \frac{2g_{t\phi}}{g_{\phi \phi}} = \frac{a}{2 + a^2}
\end{equation}
$$

The second solution is just a photon which moves along with the rotation of the black hole. The other solution is a photon which is instantaneously at rest in the BL coordinates. This is to say that at this boundary where $g_{tt} = 0$ no particles can move against the rotation of the black hole. This is also true for the region below this boundary, which is called the *ergosphere*. For $g_{tt}$ to be zero, we need $\Sigma = 2r$, meaning 

$$
\begin{align*}
    &r^2 + a^2 \cos^2\theta - 2r = 0 \\ \\
    \implies &r = 1 \pm \sqrt{1 - a^2 \cos^2\theta}.
\end{align*}
$$

Clearly there are two solutions, so there is both an inner and outer ergosphere boundary just like with the inner and outer event horizons. The following figure illustrates how the event horizons and ergosphere boundaries look with a slice of these surfaces in the $xz$-plane. This illustration is for $a = 0.99$.

![EH and ergosphere](../images/miscellaneous/EH_and_ergosphere.png)

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
    %\label{eq: camera_four_mom_ZAMO}
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

### 6.1 Why does it happen

### 6.2 Gravitational & Doppler redshift

### 6.3 Velocity distribution for the accretion disk
___

## 7 Colors

### 7.1 What colors are hot things

### 7.2 Why a blackbody

### 7.3 Redshift and apparent temperature

### 7.4 Color theory

### 7.5 From temperature to color

### 7.6 Example 1: Redshift

### 7.7 Example 2: Relativistic beaming

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

### 10.1 What

### 10.2 Simple non-relativistic example

### 10.3 Relativistic example

### 10.4 Disk delay angle

### 10.5 Example

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

## 16 References 

[1] Sean M. Carroll. Spacetime and Geometry: An Introduction to General Relativity. Cambridge University Press, 2019

[2] Masud Mansuripur. An exact derivation of the thomas precession rate using the lorentz transformation. In Henri-Jean M. Drouhin, Jean-Eric Wegrowe, and Manijeh Razeghi, editors, Spintronics XIII. SPIE, August 2020