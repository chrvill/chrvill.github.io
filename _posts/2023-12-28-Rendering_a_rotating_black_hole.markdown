---
title: Rendering a rotating black hole
date: 2023-12-28
categories: [physics, general relativity, black holes]
tags: [black holes]
img_path: '../../images/black_hole_renders'
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

## General relativistic raymarching

### Geodesics

In order to visualize the black hole we need to trace *geodesics* in the Kerr spacetime, which we do by solving the *geodesic equation*

$$
\begin{equation}
  \label{eq: geodesic_equation} \tag{5}
  \frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\rho \sigma}\frac{dx^\rho}{d\lambda}\frac{dx^\sigma}{d\lambda} = 0,
\end{equation}
$$

where $x^\mu$ is the four-position of a photon, and $\lambda$ is an affine parameter for the geodesic, chosen such that $p^\mu \equiv \frac{dx^\mu}{d\lambda}$ is the four-momentum of the photon. $\Gamma^\mu_{\rho \sigma}$ are the *Christoffel symbols* of the metric. Equation $\eqref{eq: geodesic_equation}$ is completely general, and we of course need explicit expressions for the geodesic equation for each value of $\mu$. But computing the Christoffel symbols tends to be very tedious work, and certainly so for the Kerr metric. So instead of deriving the Christoffel symbols by hand we use the Sympy package in Python to derive them, and thus also explicit expressions for each component of the geodesic equation. Code for doing this can be found \href{https://github.com/chrvill/Geodesic_EOM_deriver}{here}.

We of course have the normalization requirement

$$
\begin{equation}
    \label{eq: normalization} \tag{6}
    p_\nu p^\nu = \mu
\end{equation}
$$

with $\mu = 0$ for massless particles like photons and $\mu = -1$ for massive particles (This is really the normalization of the four-*velocity* of a massive particle, not the four-momentum. But the difference is only a factor of $m^2$, with $m$ being the mass of the particle in question). This can be written explicitly as

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

### Numerically integrating the equations of motion

Numerically \eqref{eq: geodesic_equation} is just a completely standard second order differential equation. So in that sense solving it is just like solving any other second order differential equation numerically. However it should be noted that some of the Christoffel symbols will diverge as the photon approaches the event horizon, which needs to be taken into consideration when solving the geodesic equation. Assume for the purposes of this discussion that we have a first order differential equation of the form

$$
\begin{equation}
    \label{eq: general_diff_eq} \tag{9}
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

## Relativistic aberration

### Tetrads

Before we explain aberration we will do a short detour explaining what *tetrads* are.

From the equivalence principle we know that we can always transform to a local inertial frame. And in contrast to Schwarzschild spacetime where shell observers are natural frames we will instead work with ZAMO frames, since these exist even inside the ergosphere. But the question is how we perform this transformation. How can we transform from the global BL coordinates to the local ZAMO coordinates? This is where the concept of a \textit{vierbein}, also called a *frame field*, is extremely useful (An introduction to vierbeins can be found [here](https://jila.colorado.edu/~ajsh/courses/astr5770_21/text.html})). A vierbein is a set of orthonormal axes which form a local inertial frame. So we want the vierbeins to take us from the global coordinate system to the local one, which is to say from the metric $g_{\mu \nu}$ to the Minkowski metric $\eta_{\mu \nu}$. This means that

$$
\begin{equation}
    \label{eq: vierbein_metric_diagonalization} \tag{10}
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
    \label{eq: inverse_vierbein_def} \tag{11}
    \eta_{mn} e^{*m}_{\mu} e^{*n}_{\nu} = g_{\mu \nu}.
\end{equation}
$$

From [here](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=c4f5cec0107f3f04c41f87d8130ffe23cc3b3b71) and [here](https://arxiv.org/abs/1506.01473v2) we have that the frame fields for the local ZAMO frame are given by

$$
\begin{align*}
  e_t^{\mu} &= \delta_t^{\mu} \sqrt{\frac{\Lambda}{\Delta \Sigma}} + \delta_\phi^\mu \frac{2ar}{\sqrt{\Lambda \Delta \Sigma}} \label{eq: vierbein_tmu} \tag{12} \hspace{1.1cm}\\ \\
  e_r^{\mu} &= \delta_r^{\mu} \sqrt{\frac{\Delta}{\Sigma}} \label{eq: vierbein_rmu} \tag{13} \\ \\
  e_{\theta}^{\mu} &= \delta_\theta^\mu \frac{1}{\sqrt{\Sigma}} \label{eq: vierbein_thetamu} \tag{14} \\ \\
  e_{\phi}^{\mu} &= \delta_\phi^\mu \sqrt{\frac{\Sigma}{\Lambda}}\frac{1}{\sin \theta}. \label{eq: vierbein_phimu} \tag{15}
\end{align*}
$$

while the co-frame fields are given by

$$
\begin{align*}
  e^{*t}_{\mu} &= \delta_\mu^t \sqrt{\frac{\Delta \Sigma}{\Lambda}} \label{eq: vierbein^tmu} \tag{16} \\ \\
  e^{*r}_{\mu} &= \delta_\mu^r \sqrt{\frac{\Sigma}{\Delta}} \label{eq: vierbein^rmu} \tag{17} \\ \\
  e^{*\theta}_{\mu} &= \delta_\mu^\theta \sqrt{\Sigma} \label{eq: vierbein^thetamu} \tag{18} \\ \\
  e^{*\phi}_{\mu} &= -\delta_\mu^t \frac{2ar\sin\theta}{\sqrt{\Lambda \Sigma}} + \delta_{\mu}^{\phi} \sin\theta \sqrt{\frac{\Lambda}{\Sigma}}. \label{eq: vierbein^phimu} \tag{19}
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

### Implementing aberration

Our approach to initializing the light-rays emitted by the camera is to emit rays as described in section \ref{sec: initializing_rays}. These initial values for the rays describe the $x$-, $y$- and $z$-components of the initial four-momenta of the rays. We have to convert these to BL coordinates, and then calculate $p^t$ from those components using the requirement that $p^\mu p_\mu = 0$ for photons. One might then assume that we can then just plug in this initial four-momentum $p^\mu$ of each ray into the geodesic equation to calculate the evolution of the position and four-momentum of each ray. But this is not correct in general. If the camera is moving then these initial $p^\mu$ values that we choose for the rays emitted by the camera are expressed in the instantaneous rest frame of the camera. And the geodesic equation is working in the global frame. So we have to perform the transformation from the instantaneous camera rest frame to the global frame before we can proceed to solving the geodesic equation. By implementing this transformation we will also get the effect of *relativistic aberration* for free. Relativistic aberration denotes the changes in angles between different reference frames when the frames are moving at relativistic speeds relative to each other. This aberration causes the distribution of the rays we emit from our camera to change in the global frame depending on the four-velocity of the camera. And this manifests in the simulations through changes in the field of view and an apparent "bending" of the accretion disk.

In order to hopefully alleviate some confusion regarding the approach here we list the steps here before we get into the mathematical details.

1. Compute the local velocity of the camera in the local ZAMO frame
2. Lorentz transform the initial photon four-momenta from the instantaneous camera rest frame to the local ZAMO frame
3. Transform the four-momenta of the photons from the local ZAMO frame to the global frame using the ZAMO tetrads

#### Step 1
___

From above we know how to transform from a local ZAMO frame to the global frame using the ZAMO tetrads.

## References

[1] Sean M. Carroll. Spacetime and Geometry: An Introduction to General Relativity. Cambridge University Press, 2019
