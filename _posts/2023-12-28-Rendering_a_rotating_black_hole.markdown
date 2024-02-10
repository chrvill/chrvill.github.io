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

## References

[1] Sean M. Carroll. Spacetime and Geometry: An Introduction to General Relativity. Cambridge University Press, 2019
