---
title: Rendering a rotating black hole
date: 2023-12-28
categories: [physics, general relativity, black holes]
tags: [black holes]
img_path: '../images'
math: true
---

## Introduction

The aim of this post is to describe in detail how one can create a visualization of a rotating black hole. There are of course many resources online describing how one can do this, but we generally found these to lack good explanations. Here we want to hopefully remedy this. Naturally, since black holes are consequences of general relativity, you need a fairly solid understanding of general relativity in order to understand the derivations of the equations we use. But the key *equations of motion* we use are just standard differential equations. So an understanding of how to solve differential equations numerically might suffice.

Our simulation takes into account the following relativistic effects:

* Redshift/blueshift
* Relativistic beaming
* Light travel time delay
* Relativistic aberration

In addition to this the camera is allowed to follow any arbitrary geodesic, and the effects of redshift and aberration become very apparent when the camera is moving. In order to get realistic colors, also taking into account how the redshift affects the colors, we have had to convert a blackbody spectrum to RGB colors. This is not trivial, and is therefore explained in detail here. It should also be noted that this is exclusively a simulation of the *relativistic effects* related to a black hole, while the accretion disk and jet are not physically simulated. They are just procedural volumes (*correct, Erik, if it's wrong*) where we have chosen velocity distributions. The velocity distribution for the disk assumes the particles in the disk move on circular orbits, which is admittedly not true according to the animations, but it is close. For the jet, however, we chose a completely arbitrary velocity distribution in order to get roughly the desired effect.  

### Notation and terminology

We use

* The metric signature (-, +, +, +)
* Natural units, with $c = G = 1$

### The Kerr metric

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

A relativistic effect which is unique to rotating black holes is *frame dragging*. This effect describes how the rotation of the black hole "drags" spacetime along with it, causing initially radially infalling frames to gain a non-zero angular velocity. In Schwarzschild we typically talk of *shell observers*, which are observers stationary in the Schwarzschild coordinates. But due to frame dragging these shell observers are no longer as natural in Kerr. Here it is more natural to talk of frames co-rotating with the black hole, so called *Zero Angular Momentum Observers* (ZAMO). A result of frame dragging is that there exists a region of spacetime called the *ergosphere* where it is impossible to move against the rotation of the black hole.

### Transformation between Cartesian and Boyer-Lindquist coordinates

We will need to convert vectors between Cartesian and BL coordinates multiple times, so a coordinate transformation between these is necessary. This is most easily done by studying an infinitesimal displacement vector $d\vec{r}$ expressed in the two coordinate systems. This can be written as

$$
\begin{equation}
  \label{eq: infinitesimal_displacement_Cartesian} \tag{5}
  d\vec{r} = dx \; \vec{\hat{x}} + dy \; \vec{\hat{y}} + dz \; \vec{\hat{z}}
\end{equation}
$$

in Cartesian, and

$$
\begin{equation}
  \label{eq: infinitesimal_displacement_spheroidal} \tag{6}
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

We can then determine the scaling parameter $R$ by finding the norm of $R\vec{\hat{r}}$ and using the defining property $\lvert\vec{\hat{r}}\rvert = 1$.

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
  \label{eq: r_unit_vector} \tag{7}
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
  \label{eq: unit_vector_relation} \tag{8}
  \left[\begin{matrix}
          \vec{\hat{r}} \\
          \vec{\hat{\theta}} \\
          \vec{\hat{\phi}}
        \end{matrix}\right] = \left[\begin{matrix}
                                \frac{r\sin\theta \cos\phi}{\sqrt{r^2 + a^2\cos^2\theta}} & \frac{r\sin\theta \sin\phi}{\sqrt{r^2 + a^2\cos^2\theta}} & \frac{\sqrt{r^2 + a^2}\cos\theta}{\sqrt{r^2 + a^2\cos^2\theta}} \\
                                \frac{\sqrt{r^2 + a^2}\cos\theta \cos\phi}{\sqrt{r^2 + a^2\cos^2\theta}} & \frac{\sqrt{r^2 + a^2}\cos\theta \sin\phi}{\sqrt{r^2 + a^2\cos^2\theta}} & -\frac{r\sin\theta}{\sqrt{r^2 + a^2\cos^2\theta}} \\
                                -\sin\phi & \cos\phi & 0
                                \end{matrix}\right] \left[\begin{matrix}
                                                \vec{\hat{x}} \\
                                                \vec{\hat{y}} \\
                                                \vec{\hat{z}}
                                                          \end{matrix}\right]
\end{equation}
$$

Let us call this coordinate transformation matrix $M$. It is fairly easy to check that this is an orthogonal matrix, meaning that $M^{-1} = M^T$ (just compute $M^T M$, which is equal to $\mathbb{I}$ for orthogonal matrices $M$). And in that case the inverse transformation, from BL to Cartesian coordinates, is given by $M^T$. Now, of course this is the coordinate transformation between the Cartesian and BL *basis vectors*. And we want the transformation between \textit{vector components}. But it turns out the coordinate transformation for the vector components is the same as for the basis vectors, as is shown by Lutz Lehmann [here](https://math.stackexchange.com/questions/3493647/do-coordinate-components-transform-in-the-same-or-opposite-way-as-their-bases).

## General relativistic raymarching

### Geodesics

In order to visualize the black hole we need to trace *geodesics* in the Kerr spacetime, which we do by solving the *geodesic equation*

$$
\begin{equation}
  \label{eq: geodesic_equation} \tag{9}
  \frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\rho \sigma}\frac{dx^\rho}{d\lambda}\frac{dx^\sigma}{d\lambda} = 0,
\end{equation}
$$

where $x^\mu$ is the four-position of a photon, and $\lambda$ is an affine parameter for the geodesic, chosen such that $p^\mu \equiv \frac{dx^\mu}{d\lambda}$ is the four-momentum of the photon. $\Gamma^\mu_{\rho \sigma}$ are the *Christoffel symbols* of the metric. Equation $\eqref{eq: geodesic_equation}$ is completely general, and we of course need explicit expressions for the geodesic equation for each value of $\mu$. But computing the Christoffel symbols tends to be very tedious work, and certainly so for the Kerr metric. So instead of deriving the Christoffel symbols by hand we use the Sympy package in Python to derive them, and thus also explicit expressions for each component of the geodesic equation. Code for doing this can be found \href{https://github.com/chrvill/Geodesic_EOM_deriver}{here}.

We of course have the normalization requirement

$$
\begin{equation}
    \label{eq: normalization} \tag{10}
    p_\nu p^\nu = \mu
\end{equation}
$$

with $\mu = 0$ for massless particles like photons and $\mu = -1$ for massive particles (This is really the normalization of the four-*velocity* of a massive particle, not the four-momentum. But the difference is only a factor of $m^2$, with $m$ being the mass of the particle in question). This can be written explicitly as

$$
\begin{equation}
    \label{eq: normalization_explicit} \tag{11}
    g_{tt} \left(p^t\right)^2 + g_{rr}\left(p^r\right)^2 + g_{\theta\theta} \left(p^\theta\right)^2 + g_{\phi \phi}\left(p^\phi\right)^2 + 2g_{t\phi} p^t p^\phi = \mu
\end{equation}
$$

When we choose initial directions for our rays we will essentially provide the spatial components of the four-momentum, so we need to compute $p^t$ from $p^r, p^\theta$ and $p^\phi$. Solving \eqref{eq: normalization_explicit} for $p^t$ gives

$$
\begin{equation}
    \label{eq: p^t expression} \tag{12}
    p^t = -\frac{g_{t\phi}}{g_{tt}} p^\phi \pm \sqrt{\left(\frac{g_{t\phi}}{g_{tt}} p^\phi\right)^2 - \frac{1}{g_{tt}}\left(g_{rr} \left(p^r\right)^2 + g_{\theta\theta} \left(p^\theta\right)^2 + g_{\phi \phi} \left(p^\phi\right)^2 - \mu\right)}
\end{equation}
$$

But we have to find out which sign is appropriate. Consider a case where $p^\phi = 0$, which is only possible outside the ergosphere. Then, since $p^t$ needs to be positive we need to use the $+$ sign outside the ergosphere. Arguing for which sign is appropriate inside the ergosphere proved trickier, but we found that we had to use the $-$ sign inside the ergosphere in order for $p^t$ to increase continuously as you pass the boundary of the ergosphere.

### Numerically integrating the equations of motion

Numerically \eqref{eq: geodesic_equation} is just a completely standard second order differential equation. So in that sense solving it is just like solving any other second order differential equation numerically. However it should be noted that some of the Christoffel symbols will diverge as the photon approaches the event horizon, which needs to be taken into consideration when solving the geodesic equation. Assume for the purposes of this discussion that we have a first order differential equation of the form

$$
\begin{equation}
    \label{eq: general_diff_eq} \tag{13}
  \frac{dy}{dt} = f(t, y).
\end{equation}
$$

Our second order equations of motion can of course be expressed in terms of first order DEs by just by writing the corresponding equations for $\frac{dx^\mu}{d\lambda}$ and $\frac{d^2 x^\mu}{d\lambda^2}$. With an appropriate choice of variable timestep the geodesic equation could most likely be solved using for example the 4th order Runge-Kutta integration scheme. But here, due to an unrelated error that will be discussed later, we chose to use the *Runge-Kutta-Fehlberg* scheme, abbreviated RKF45. This is, as the name implies, also in the family of Runge-Kutta methods. The main benefit of using this scheme is that it allows us to calculate an adaptive timestep $h$ very easily - the scheme itself can choose a small timestep when needed and revert to a bigger timestep when things evolve slowly. The scheme itself is similar to 4th order Runge-Kutta, and is described in detail [here](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method), where we use the table listed under "FORMULA 2". It is a 4th order accurate method where we perform a few more function evaluations in order to obtain a measure of the error associated with the scheme. When the right hand side of the differential equation is not explicitly a function of time, like here, the algorithm can be written as

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

### Relativistic aberration

#### Tetrads

Before we explain aberration we will do a short detour explaining what *tetrads* are.

From the equivalence principle we know that we can always transform to a local inertial frame. And in contrast to Schwarzschild spacetime where shell observers are natural frames we will instead work with ZAMO frames, since these exist even inside the ergosphere. But the question is how we perform this transformation. How can we transform from the global BL coordinates to the local ZAMO coordinates? This is where the concept of a \textit{vierbein}, also called a *frame field*, is extremely useful (An introduction to vierbeins can be found [here](https://jila.colorado.edu/~ajsh/courses/astr5770_21/text.html})). A vierbein is a set of orthonormal axes which form a local inertial frame. So we want the vierbeins to take us from the global coordinate system to the local one, which is to say from the metric $g_{\mu \nu}$ to the Minkowski metric $\eta_{\mu \nu}$. This means that

$$
\begin{equation}
    \label{eq: vierbein_metric_diagonalization} \tag{14}
    g_{\mu \nu} e_m^{\mu} e_n^{\nu} = \eta_{m n},
\end{equation}
$$

where $e_m^{\mu}$ is such a vierbein field or a frame field. Keep in mind that Latin indices are here used to refer to the indices in the local inertial frame. We also have the dual vierbein field $e^{*m}_{\mu}$, also called the co-frame field, defined through

$$
\begin{equation}
  e_m^{\mu} e^{*m}_{\nu} = \delta_{\nu}^{\mu}, \quad e_m^{\mu} e^{*n}_{\mu} = \delta_m^n.
\end{equation}
$$

Using this we can also write

$$
\begin{equation}
    \label{eq: inverse_vierbein_def} \tag{15}
    \eta_{mn} e^{*m}_{\mu} e^{*n}_{\nu} = g_{\mu \nu}.
\end{equation}
$$

From [here](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=c4f5cec0107f3f04c41f87d8130ffe23cc3b3b71) and [here](https://arxiv.org/abs/1506.01473v2) we have that the frame fields for the local ZAMO frame are given by

$$
\begin{align*}
  e_t^{\mu} &= \delta_t^{\mu} \sqrt{\frac{\Lambda}{\Delta \Sigma}} + \delta_\phi^\mu \frac{2ar}{\sqrt{\Lambda \Delta \Sigma}} \label{eq: vierbein_tmu} \tag{16} \hspace{1.1cm}\\ \\
  e_r^{\mu} &= \delta_r^{\mu} \sqrt{\frac{\Delta}{\Sigma}} \label{eq: vierbein_rmu} \tag{17} \\ \\
  e_{\theta}^{\mu} &= \delta_\theta^\mu \frac{1}{\sqrt{\Sigma}} \label{eq: vierbein_thetamu} \tag{18} \\ \\
  e_{\phi}^{\mu} &= \delta_\phi^\mu \sqrt{\frac{\Sigma}{\Lambda}}\frac{1}{\sin \theta}. \label{eq: vierbein_phimu} \tag{19}
\end{align*}
$$

while the co-frame fields are given by

$$
\begin{align*}
  e^{*t}_{\mu} &= \delta_\mu^t \sqrt{\frac{\Delta \Sigma}{\Lambda}} \label{eq: vierbein^tmu} \tag{20} \\ \\
  e^{*r}_{\mu} &= \delta_\mu^r \sqrt{\frac{\Sigma}{\Delta}} \label{eq: vierbein^rmu} \tag{21} \\ \\
  e^{*\theta}_{\mu} &= \delta_\mu^\theta \sqrt{\Sigma} \label{eq: vierbein^thetamu} \tag{22} \\ \\
  e^{*\phi}_{\mu} &= -\delta_\mu^t \frac{2ar\sin\theta}{\sqrt{\Lambda \Sigma}} + \delta_{\mu}^{\phi} \sin\theta \sqrt{\frac{\Lambda}{\Sigma}}. \label{eq: vierbein^phimu} \tag{23}
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

#### Implementing aberration

Our approach to initializing the light-rays emitted by the camera is to emit rays as described in section \ref{sec: initializing_rays}. These initial values for the rays describe the $x$-, $y$- and $z$-components of the initial four-momenta of the rays. We have to convert these to BL coordinates, and then calculate $p^t$ from those components using the requirement that $p^\mu p_\mu = 0$ for photons. One might then assume that we can then just plug in this initial four-momentum $p^\mu$ of each ray into the geodesic equation to calculate the evolution of the position and four-momentum of each ray. But this is not correct in general. If the camera is moving then these initial $p^\mu$ values that we choose for the rays emitted by the camera are expressed in the instantaneous rest frame of the camera. And the geodesic equation is working in the global frame. So we have to perform the transformation from the instantaneous camera rest frame to the global frame before we can proceed to solving the geodesic equation. By implementing this transformation we will also get the effect of *relativistic aberration* for free. Relativistic aberration denotes the changes in angles between different reference frames when the frames are moving at relativistic speeds relative to each other. This aberration causes the distribution of the rays we emit from our camera to change in the global frame depending on the four-velocity of the camera. And this manifests in the simulations through changes in the field of view and an apparent "bending" of the accretion disk.

In order to hopefully alleviate some confusion regarding the approach here we list the steps here before we get into the mathematical details.

1. Compute the local velocity of the camera in the local ZAMO frame
2. Lorentz transform the initial photon four-momenta from the instantaneous camera rest frame to the local ZAMO frame
3. Transform the four-momenta of the photons from the local ZAMO frame to the global frame using the ZAMO tetrads

##### Step 1
___

From above we know how to transform from a local ZAMO frame to the global frame using the ZAMO tetrads. So we will first have to transform from the instantaneous camera rest frame to the local ZAMO frame, in order to then use the tetrads from Appendix \ref{app: tetrads} to finally transform to the global frame. Let the camera have four-velocity $u^\mu$ in the global coordinate frame. Also note for later that this four-velocity will be expressed in BL coordinates. Imagine then a ZAMO frame at the same position as the camera. The four-velocity of the camera in the ZAMO tetrad frame is given by

$$
\begin{equation}
    \label{eq: local_camera_four_vel} \tag{24}
    \tilde{u}^m = e^m_{\mu} u^{\mu},
\end{equation}
$$

Since this four-velocity is measured in a local inertial frame it can also be expressed as

$$
\begin{equation}
  \label{eq: ZAMO_four_velocity_camera} \tag{25}
  \tilde{u}^m = \gamma\left(1, \vec{v}\right),
\end{equation}
$$

where $\gamma = \frac{1}{\sqrt{1 - \vec{v}^2}}$ and $\vec{v}$ is the local velocity of the camera in the ZAMO frame. Then $\gamma = \tilde{u}^t$, and so the components $v^i$ of the local velocity $\vec{v}$ of the camera in the ZAMO frame are given by

$$
\begin{equation}
  \label{eq: ZAMO_local_camera_velocity} \tag{26}
  v^i = \frac{\tilde{u}^i}{\gamma} = \frac{\tilde{u}^i}{\tilde{u}^t}.
\end{equation}
$$

##### Step 2
___

Consider now the rays that we emit from our camera. As mentioned, we first need to transform to the local ZAMO frame. This is just a Lorentz transformation, whose general form is the following (see for example [2])

$$
\begin{equation}
  \Lambda^m_{\;\;n} = \left[\begin{matrix}\gamma && -\gamma v^x && -\gamma v^y && -\gamma v^z\ \\
                          -\gamma v^x && 1 + \left(\gamma - 1\right)\frac{\left(v^x\right)^2}{\vec{v}^2} && \left(\gamma - 1\right) \frac{v^x v^y}{\vec{v}^2} && \left(\gamma - 1\right)\frac{v^x v^z}{\vec{v}^2} \\
                          -\gamma v^y && \left(\gamma - 1\right)\frac{v^x v^y}{\vec{v}^2} && 1 + \left(\gamma - 1\right)\frac{\left(v^y\right)^2}{\vec{v}^2} && \left(\gamma - 1\right)\frac{v^y v^z}{\vec{v}^2} \\
                          -\gamma v^z && \left(\gamma - 1\right)\frac{v^x v^z}{\vec{v}^2} && \left(\gamma - 1\right)\frac{v^y v^z}{\vec{v}^2} && 1 + \left(\gamma - 1\right)\frac{\left(v^z\right)^2}{\vec{v}^2}
                         \end{matrix}\right]
\end{equation}
$$

in Cartesian coordinates. So we first have to convert the velocity vector $\vec{v}$ from BL coordinates to Cartesian coordinates. From Appendix \ref{sec: coordinate_transformation} we know that the components of the velocity vector in Cartesian components are  

$$
\begin{align*}
  v^x &= \frac{1}{\sqrt{r^2 + a^2\cos^2\theta}}\left[r\sin\theta \cos\phi \,v^r + \sqrt{r^2 + a^2}\cos\theta \cos\phi \,v^\theta\right] - \sin\phi \,v^\phi \\ \\
  v^y &= \frac{1}{\sqrt{r^2 + a^2 \cos^2\theta}}\left[r\sin\theta \sin\phi \,v^r + \sqrt{r^2 + a^2}\cos\theta \sin\phi \,v^\theta\right] + \cos\phi \,v^\phi \\ \\
  v^z &= \frac{1}{\sqrt{r^2 + a^2 \cos^2\theta}}\left[\sqrt{r^2 + a^2}\cos\theta \,v^r - r\sin\theta \,v^\theta\right].
\end{align*}
$$

when the components are $v^r, v^\theta$ and $v^\phi$ in BL coordinates. Let us call the four-momentum of the initial photon rays $\tilde{p}{^m}$ in the camera rest frame. In the ZAMO frame the four-momentum is then given by the primed

$$
\begin{equation}
    \label{eq: camera_four_mom_ZAMO} \tag{27}
    \tilde{p}'^m = \Lambda^m_n\left(\vec{v}\right) \tilde{p}^n
\end{equation}
$$

Then we should of course move back into BL coordinates, since that is the coordinate system we work prefer here. And that coordinate transformation is simply given by the matrix in \eqref{eq: unit_vector_relation}.

##### Step 3
___

The final step is then rather simple. The four-momentum of the ray in the global frame is simply

$$
\begin{equation}
    \label{eq: camera_four_mom_global} \tag{28}
    p^\mu = e_m^{\mu} \tilde{p}'^m,
\end{equation}
$$

with the understanding that $\tilde{p}'^m$ is expressed in BL coordinates, as described at the end of Step 2.

Having implemented this we have a fully general approach which takes into account relativistic aberration. And this is implemented purely through a coordinate transformation from the camera's rest frame to the global frame, which is in any case necessary.

## The accretion disk

### The four-velocity of the accretion disk

The accretion disk is in reality moving in a spiral towards the black hole. But when dealing with the physics here we will simplify and assume all points in the disk move in circular orbits. From [3] we know that the four-velocity of a massive particle in a prograde circular orbit in the equatorial plane around a Kerr black hole is given by $u^\mu = \left(u^t, 0, 0, u^\phi\right)$, where

$$
\begin{equation}
  u^\phi = \frac{\sqrt{Mr}}{r\sqrt{r^2 - 3Mr + 2\lvert a\rvert \sqrt{Mr}}}.
\end{equation}
$$

And $u^t$ is of course determined by the condition $u^\mu u_\mu = -1$.

### Redshift

Consider an observer at some point $P$ in spacetime with four-velocity $u^\mu$ in some coordinate system. Assume we also have a photon at the same point $P$ in spacetime with four-momentum $p^\mu = \frac{dx^\mu}{d\lambda}$ in that same coordinate system. The frequency that the observer measures for the photon at $P$ is

$$
\begin{equation}
  \label{eq: redshift_definition} \tag{29}
  \nu = -p^\mu u_\mu = -g_{\mu \rho} p^\mu u^\rho,
\end{equation}
$$

which is of course coordinate invariant. Here redshift is of interest because we have an "observer" at the accretion disk around the black hole that emits light with certain frequencies, and that light is then observed by the camera at another point in spacetime\footnote{Of course in our simulation it is the other way around, but in reality the photon would be emitted by the accretion disk and be observed by the camera.}. And since we have solved the geodesic equation for the rays we already know the four-momentum of the photons at all points along their geodesics. In particular, we know their four-momenta at the disk and at the camera. Consider a point on the disk at coordinates $\left(r, \theta, \phi\right)$ and let the camera be at $\left(R, \Theta, \Phi\right)$. Let the parameter $\lambda$ equal 0 at the camera and $\lambda_1 > 0$ at the point on the disk. We will call the four-velocity of the relevant point on the disk $u^\mu_\mathrm{disk}$ and the four-velocity of the camera $u^\mu_\mathrm{camera}$. Then the redshift $1 + z$ for the light emitted at the disk and observed by our camera is given by

$$
\begin{equation}
  \label{eq: redshift_factor} \tag{30}
  1 + z = \frac{\lambda_\mathrm{camera}}{\lambda_\mathrm{disk}} = \frac{\nu_\mathrm{disk}}{\nu_\mathrm{camera}} = \frac{g_{\mu \nu}\left(r, \theta, \phi\right) p^\mu\left(\lambda = \lambda_1\right) u^\nu_\mathrm{disk}}{g_{\mu \nu}\left(R, \Theta, \Phi\right) p^\mu \left(\lambda = 0\right) u^\nu_\mathrm{camera}}
\end{equation}
$$

where $\nu_\mathrm{disk}$ and $\nu_\mathrm{camera}$ are the frequencies measured at the disk and the camera respectively, and likewise $\lambda_\mathrm{disk}$ and $\lambda_\mathrm{camera}$ are the measured wavelengths at the disk and camera respectively. This general expression has the advantage that it does not assume anything about the motion of neither the disk nor the camera. So this places no restrictions on how we can move our camera around or what kind of velocity distribution the disk can have. Also note that by "redshift" we really mean both redshift and blueshift. We are just referring to both decrease and increase in wavelength under the same umbrella term for simplicity, which is really just a habit from cosmology. It is also important to note that $\eqref{eq: redshift_factor}$ accounts for both the Doppler-like redshift caused by motion and the gravitational redshift. Both of these two kinds of redshift are baked into the metric itself.

### Relativistic beaming

This redshift also causes a change in the brightness that is observed by the camera. Discussion of this can be found in [4]. According to Liouville's theorem the quantity $I_\nu/\nu^3$ is invariant. And the intensity $I_\nu$ is defined in terms of differentials, such that $I_\nu d\nu = - I_\lambda d\lambda$, since $\lambda = \frac{1}{\nu}$. We therefore have

$$
\begin{equation}
  I_\nu = -I_\lambda \frac{d\lambda}{d\nu} = I_\lambda \frac{1}{\nu^2} = I_\lambda \lambda^2
\end{equation}
$$

Therefore the quantity $I_\lambda \lambda^5$ is invariant. This means that the brightness $I_\lambda$ that is observed by the camera is

$$
\begin{equation}
  I_\lambda = \left(\frac{\lambda'}{\lambda}\right)^5 I_\lambda'
\end{equation}
$$

where $\lambda'$ is the emitted wavelength and $\lambda$ the observed one. And $I_\lambda'$ is the intensity at emission. Therefore

$$
\begin{equation}
  \label{eq: beaming} \tag{31}
  I_\lambda = \left(1 + z\right)^{-5} I_\lambda'
\end{equation}
$$

is the observed intensity. This causes an increase(decrease) when the light is blueshifted(redshifted), as one might intuitively expect.

### Light travel time delay

One aspect of the visualizations which might not be obvious is that we have to take into account how long it takes light to travel from each part of the disk to the camera. Since it takes different amounts of time for light from different parts of the disk to reach us the photons that we receive at some time $t$ must have been emitted from the different parts of the disk at different times. So if we for example take our image at global time $t$ and it takes a time $t_\mathrm{A}$ for light from a point A on the disk to reach us, then the light that reaches us at time $t$ must have left A at time $t - t_\mathrm{A}$. So when the travel time $t_\mathrm{A}$ differs between the points on the disk the parts of the disk we see will correspond to different emission times. Let us assume we have a function $f(\vec{x}, t)$ which calculates the density at the position $\vec{x}$ on the disk evaluated at global time $t$. Then the density that the camera actually observes at point $\vec{x}_i$ when it takes an image at time $t$ is $f(\vec{x}_i, t - t_i)$, where $t$ is the time when the image is taken and $t_i$ is the time it takes light to travel from point $\vec{x}_i$ on the disk to the camera.

A convenient aspect of the geodesic tracing here is that we already calculate this travel time in the global coordinate system when we solve the geodesic equation. And since we assume that each point on the disk moves along a circular orbit this correction due to the travel time only requires a correction for the angle that a given point on the disk moves during the travel time. Furthermore, since the four-velocity of the points on the disk is constant this angular correction is particularly simple:

$$
\begin{equation}
    \label{eq: light_travel_delay_angle} \tag{32}
    \Delta \phi = \int_0^{t_\mathrm{travel}} \frac{d\phi}{dt'} dt' = \int_0^{t_\mathrm{travel}} \frac{d\phi}{d\tau} \frac{d\tau}{dt'} dt' = \frac{u^\phi}{u^t}t_\mathrm{travel},
\end{equation}
$$

where $u^\phi = \frac{d\phi}{d\tau}$ is the $\phi$-component of the four-velocity of the relevant point on the disk, and $\tau$ is the proper time of that same point. And $u^t = \frac{dt}{d\tau}$ is the time-component of the four-velocity of the point on the disk. $t$ is again the travel time from the camera to the relevant point on the disk.

## Color theory

Using $\eqref{eq: redshift_factor}$ we can compute the redshift for all points on the disk. But the computer screen of course renders RGB colors, not wavelengths directly. The disk emits photons at all wavelengths, and the distribution of photons emitted as a function of wavelength is not uniform. Here we will assume that the disk is a perfect *blackbody* - meaning it is a perfect absorber of light - the so-called *intensity distribution* is given by

$$
\begin{equation}
  \label{eq: blackbody distribution} \tag{33}
  I_{\lambda}\left(\lambda; T\right) = \frac{2h c^2}{\lambda^5} \frac{1}{e^{hc/\lambda k_B T} - 1},
\end{equation}  
$$

where $h$ and $k_B$ are Planck's constant and Boltzmann's constant respectively. $T$ is the temperature of the disk and $\lambda$ wavelength. This is the *blackbody distribution*, and describes the amount of radiation emitted as a function of wavelength. We now need a way to convert an intensity distribution, like the blackbody distribution, to an RGB color.

The relative amounts of the visible wavelengths that are emitted determine the color of the object. Unfortunately our eyes do not perceive and respond to all visible wavelengths of light equally. So converting an intensity spectrum $I_\lambda$, which represents the intensity distribution that is actually observed, to the color our eyes actually see is far from trivial. The procedure for calculating an RGB color comes in two main steps: converting a spectrum to so-called \textit{tristimulus values}, and then converting the tristimulus values to RGB values. We found [5], [this](https://color.org/chardata/rgb/sRGB.pdf) and [this](https://en.wikipedia.org/wiki/SRGB) to be particularly helpful here. Assume now that we have a blackbody at some temperature $T$. We can find the tristimulus values $X, Y$ and $Z$ using the so-called *color matching functions* $\bar{x}, \bar{y}$ and $\bar{z}$, which describe how our eyes respond to different wavelengths. The $X, Y$ and $Z$ values are calculated in the following way

$$
\begin{align*}
  X &= \int_{\mathrm{380 nm}}^{\mathrm{780 nm}} I_\lambda \left(\lambda; T\right) \bar{x}\left(\lambda\right)d\lambda \label{eq: color_X} \tag{34} \\ \\
  Y &= \int_{\mathrm{380 nm}}^{\mathrm{780 nm}} I_\lambda \left(\lambda; T\right) \bar{y}\left(\lambda\right)d\lambda \label{eq: color_Y} \tag{35} \\ \\
  Z &= \int_{\mathrm{380 nm}}^{\mathrm{780 nm}} I_\lambda \left(\lambda; T\right) \bar{z}\left(\lambda\right)d\lambda. \label{eq: color_Z} \tag{36}
\end{align*}
$$

Approximations to these color matching functions can be found in [6]. But it turns out that at low temperatures these approximations are not accurate enough anymore. So we will instead be using the lookup table given in [7]. In the following figure we have plotted the color matching functions as functions of wavelength. The $x$-axis represents wavelength $\lambda$ in nm, while the $y$-axis represents the output of the color matching functions, which can be taken to have units such that $X, Y$ and $Z$ are dimensionless.

![The color matching functions](miscellaneous/color_matching_functions.png)

## References

[1] Sean M. Carroll. Spacetime and Geometry: An Introduction to General Relativity. Cambridge University Press, 2019

[2] Masud Mansuripur. An exact derivation of the thomas precession rate using the lorentz transformation. In Henri-Jean M. Drouhin, Jean-Eric Wegrowe, and Manijeh Razeghi, editors, Spintronics XIII. SPIE, August 2020

[3] Z. Kh. Kurmakaev. Circular orbits in the Kerr metric. , 18:110, August 1974

[4] Eric Bruneton. Real-time high-quality rendering of non-rotating black holes, 2020

[5] Philippe Colantoni. Color space transformations. 2006

[6] Chris Wyman et al. Simple analytic approximations to the cie xyz color matching functions. 2013

[7] Kasajima Ichiro. Plotting colors on color circle: Interconversion between xyz values and rgb color system. Current Trends in Analytical and Bioanalytical Chemistry, 1, 05 2017
