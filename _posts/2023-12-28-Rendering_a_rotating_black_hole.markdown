---
title: Rendering a rotating black hole
date: 2023-12-28
categories: [physics, general relativity, black holes]
tags: [black holes]
math: true
---

![Animation](https://github.com/chrvill/chrvill.github.io/blob/main/images/black_hole_renders/black_hole_animation.gif)

## 1 Introduction

The aim of this post is to describe in detail how one can create a visualization of a rotating black hole. There are of course many resources online describing how one can do this, but we generally found these to lack good explanations. Here we want to hopefully remedy this. Naturally, since black holes are consequences of general relativity, you need a fairly solid understanding of general relativity in order to understand the derivations of the equations we use. But the key *equations of motion* we use are just standard differential equations. So an understanding of how to solve differential equations numerically might suffice.

Our simulation takes into account the following relativistic effects:

* Redshift/blueshift
* Relativistic beaming
* Light travel time delay
* Relativistic aberration

In addition to this the camera is allowed to follow any arbitrary geodesic, and the effects of redshift and aberration become very apparent when the camera is moving. In order to get realistic colors, also taking into account how the redshift affects the colors, we have had to convert a blackbody spectrum to RGB colors. This is not trivial, and is therefore explained in detail here. It should also be noted that this is exclusively a simulation of the *relativistic effects* related to a black hole, while the accretion disk and jet are not physically simulated. They are just procedural volumes where we have chosen velocity distributions. The velocity distribution for the disk assumes the particles in the disk move on circular orbits, which is admittedly not true according to the animations, but it is close. For the jet, however, we chose a completely arbitrary velocity distribution in order to get roughly the desired effect.  

## 2 Fundamental black hole physics

Black holes are a product of the theory of *general relativity* (GR), so we need to discuss the GR in order to describe how one can create a visualization of a black hole. Relativity is commonly divided into special relativity (SR) and GR in pop-science, despite the fact that GR entirely encompasses SR. The distinction is made because SR assumes a spacetime with no matter- or energy-content, while GR allows an arbitrary distribution of matter and energy. Central to both GR and SR is the concept of *spacetime*. Spacetime can be thought of as a combination of the familiar 3 spatial dimensions and the 1 dimension of time, commonly described as a (3 + 1)-dimensional spacetime. 

### Special relativity

Before we talk about GR it is instructional to discuss SR and introduce the relevant mathematical framework in the simpler context of SR. 

When a particle moves through normal 3-dimensional space we can trace out its path to describe its motion through space. Completely analogously we can trace out a particle's path through *spacetime* to describe its motion in spacetime. This trajectory through spacetime is called the *worldline* of the particle. The following figure shows what is called a *spacetime diagram* limited to one spatial dimension, here chosen as the $x$-direction, and the time dimension $t$. The solid curve represents the worldline of some massive particle, and the dashed line represents the worldline of a massless particle like a photon.

![Spacetime diagram](../images/miscellaneous/spacetime_diagram.png)

Photons always move at the speed of light $c$ in all locally inertial reference frames. Here *local* means that we are talking about a small region in spacetime, meaning a small interval of time and/or a small region of space. And *inertial* means non-accelerating. For a photon moving in the $x$-direction its $x$-coordinate is given by $x(t) = ct$, where $t$ is time. Here, however, we use so-called *natural units* in which $c \equiv 1$. Then $x(t) = t$, so the worldline of a photon is linear with a slope of 1, as shown in the spacetime diagram above. Notice that all physical particles (which have to move at or below the speed of light) must have a slope which is greater than or equal to 1 at all points along their worldlines. It should also be noted that this discussion is implicitly assuming that we are talking about the spacetime dealt with in SR - so-called *Minkowski spacetime*. There are a few nuances that make the discussion in a general spacetime more complicated, but the discussion above paints a sufficiently reasonable picture for us here. 

The so-called *line element* (also called the *spacetime interval*) is central to relativity. It describes the geometry of spacetime, specifically how to compute distances. We denote the line element by $ds^2$, and in our normal 3-dimensional space, called *Euclidean* space, it is given by 

$$
\begin{equation}
  \label{eq: Euclidean_metric} \tag{1}
  ds^2 = dx^2 + dy^2 + dz^2
\end{equation}
$$

where $dx, dy$ and $dz$ are infinitesimal displacements along the $x-, y-$ and $z-$directions. Equation $\eqref{eq: Euclidean_metric}$ is just the Pythagorean theorem, thus justifying our interpretation of the line element as a measure of distance. In index notation we can write the line element as

$$
\begin{equation}
  \label{eq: Euclidean_metric_index_notation} \tag{2}
  ds^2 = \sum_{i = 1}^n \sum_{j = 1}^n \delta_{ij} dx^i dx^j,
\end{equation}
$$

where we have collected $dx, dy$ and $dz$ into the vector $dx^i = (dx, dy, dz)$ and where $n$ is the dimension of the Euclidean space. And $\delta_{ij}$ is the Kronecker-delta, defined as $\delta_{ij} = 1$ if $i = j$ and $\delta_{ij} = 0$ if $i \neq j$. We now introduce a new notation that is standard when working with relativity, namely the \textit{Einstein notation convention}. With this convention we omit the summation symbols in expressions like $\eqref{eq: Euclidean_metric_index_notation}$ with the understanding that any time we have a repeated index that index is summed over. We can then write the line element of Euclidean space as

$$
\begin{equation}
  \label{eq: Euclidean_metric_einstein_convention} \tag{3}
  ds^2 = \delta_{ij} dx^i dx^j
\end{equation}
$$

As mentioned earlier, SR deals with Minkowski spacetime. But this spacetime is not just a simple extension of Euclidean space where time is treated the same as the spatial dimensions. Instead the Minkowski line element is given by

$$
\begin{equation}
  \label{eq: Minkowski_metric} \tag{4}
  ds^2 = -dt^2 + dx^2 + dy^2 + dz^2
\end{equation}
$$

where $dt$ is an infinitesimal interval of time. In analogy with $\eqref{eq: Euclidean_metric_einstein_convention}$ we can write the Minkowski line element as

$$
\begin{equation}
  \label{eq: Minkowski_metric_tensor_notation} \tag{5}
  ds^2 = \eta_{\mu \nu} dx^\mu dx^\nu
\end{equation}
$$

where $dx^\mu = (dt, dx, dy, dz)$ and $\eta_{\mu \nu}$ is the Minkowski \textit{metric tensor} and is given by

$$
\begin{equation}
  \label{eq: Minkowski_metric_tensor} \tag[6]
  \eta_{\mu \nu} = \left(\begin{matrix}
                      -1 & 0 & 0 & 0 \\
                      0 & 1 & 0 & 0 \\
                      0 & 0 & 1 & 0 \\
                      0 & 0 & 0 & 1
                   \end{matrix}\right) = \mathrm{diag}\left(-1, 1, 1, 1\right).
\end{equation}
$$

Note that we have transitioned from using Latin letters like $i$ and $j$ to represent indices to using Greek letters like $\mu$ and $\nu$, in accordance with the standard notation in relativity. Similarly to how $\eqref{eq: Euclidean_metric_einstein_convention}$ tells us how to calculate distances in Euclidean space $\eqref{eq: Minkowski_metric_tensor_notation}$ tells us how to calculate distances in Minkowski spacetime. Consider now two events separated by a time $d\tau$ which are located at the same position, meaning $dx = dy = dz = 0$. Then $ds^2 = -d\tau^2$, where $\tau$ is called the \textit{proper time} between the two events. This proper time is always the time measured in the rest frame of the events, where the two events occur at the same location. Events which are such that the line element can be written in terms of a proper time like this are said to have \textit{timelike separation}. A very important property of the line element is that it has the same value in all reference frames. So this means that some other observer measuring the line element between the two events, where $dx, dy$ and $dz$ are not necessarily 0, will measure $dt, dx, dy$ and $dz$ to be in accordance with

$$
\begin{equation}
    ds^2 = -dt^2 + dx^2 + dy^2 + dz^2 \equiv -d\tau^2.
\end{equation}
$$

### Metrics and exact solutions

### Geodesics

### The Schwarzschild solution

### The Kerr solution


## Tracing geodesics in curved spacetime

### Coordinate systems

### Metric and equations of motion

### Initial conditions

### Flat/Schwarzschild/Kerr example


## A simple pin-hole camera model

### Raymarching

### Rendering a sphere

### Reference frames

### Coordinate transformation from Cartesian to BL

### Examples

### 2D thin disk

### Examples

### Problems and limitations


## The non-stationary camera

### What does this mean

### Basic frame of reference transformation

### The ZAMO and ergosphere

### Rest frame to ZAMO to global transformation

### Examples

### Aberration of light 

### Free-falling camera

### Examples


## Redshift 

### Why does it happen

### What do we want

### Gravitational redshift

### Velocity distributions

### Velocity redshift


## Colors

### What colors are hot things

### Why a blackbody

### Redshift and apparent temperature

### Color matching functions

### From temperature to color

### Example 1: Doppler beaming

### Example 2: Dopple shift


## The volumetric accretion disk

### Disk models

### Basic volume rendering 

### Challenges

### Procedural voronoi based disk noise

### Examples

### Temperature distributions

### Examples


## Light travel delay

### What

### Simple non-relativistic example

### Relativistic example

### Disk delay angle

### Example


## The astrophysical jet

### Jet models

### Procedural jet model with special coordinates

### Velocity distribution

### Temperature distribution

### Examples

### Why does it look so weird

### Delaying the jet

### Examples

### Why does it look even weirder


## Improving the renders

### Super sampling

### Anti aliasing

### Physical motion blur


## Test cases

### Example 1

### Example 2 

### Example 3

### Contemporary example comparisons


## Conclusion

### Limitations

### Improvements

### Special thanks

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

### Dimensional analysis

Let us introduce the dimensionless quantities

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

Let us call this coordinate transformation matrix $M$. It is fairly easy to check that this is an orthogonal matrix, meaning that $M^{-1} = M^T$ (just compute $M^T M$, which is equal to $\mathbb{I}$ for orthogonal matrices $M$). And in that case the inverse transformation, from BL to Cartesian coordinates, is given by $M^T$. Now, of course this is the coordinate transformation between the Cartesian and BL *basis vectors*. And we want the transformation between *vector components*. But it turns out the coordinate transformation for the vector components is the same as for the basis vectors, as is shown by Lutz Lehmann [here](https://math.stackexchange.com/questions/3493647/do-coordinate-components-transform-in-the-same-or-opposite-way-as-their-bases).

## General relativistic raymarching

### Geodesics

In normal raymarching one typically just casts straight lines to model how light moves in so-called *flat spacetime* - spacetime with negligible gravity. But in curved spacetime light no longer follows these straight paths through space. Instead it follows the curved-space generalization of "straight paths", namely *geodesics*. In "normal" spaces (meaning not spacetime, simply space. The proper terminology is *Euclidean* space) a geodesic can be thought of as the path which minimizes the length between two points $A$ and $B$. In Euclidean space (think of a plane, for example) this just reduces to a straight line as we expect. But we can get more interesting geodesics if we consider a sphere, for example. This is shown in the following figure, where the red curve shows the geodesic between two points $A$ and $B$ on the sphere. 

![Sphere geodesic](../images/miscellaneous/A-geodesic-on-the-surface-of-a-sphere.png)

This red line, the geodesic, is the shortest possible path between $A$ and $B$ if we are constrained to move on the spherical surface. You might also notice that the geodesic lies along a line of constant longitude, meaning that lines of constant longitude represent "straight" lines on the sphere. What this concept of a path being "straight" here really means is simply that you can follow such a path without having to steer - you can just move in your forwards direction. It might not be obvious, however, that lines of constant latitude are *not* geodesics on the sphere - they are not "straight". This is to say that if you were constrained to the surface of a sphere you could move along lines of constant longitude without steering, but not along lines of constant latitude. The concept of geodesics is the reason why the flight path of a plane can look unreasonable when projected into a two dimensional map, even though it is actually the shortest path.

The mathematics behind geodesics carries over to curved spacetime, but the intuition is harder to grasp. Rather than being the path of shortest spatial distance between two points it is the path of *longest proper time*, which is the time measured by an observer following the path. Let us assume we have two observers, which we creatively call Alice and Bob, moving along two different paths from point $A$ to point $B$ in spacetime. Alice moves along a geodesic between these points, while Bob does not. This then means that Alice will measure a longer time from $A$ to $B$ on her clock than Bob will measure on his clock. 

In relativity the geodesic also has the interpretation of being the path an unaccelerated particle follows. In our previous analogy then, Alice is in unaccelerated free-fall, while Bob accelerates in some way. Now, a possible cause for confusion is the fact that in classical physics we tend to say that an object accelerates due to gravity. But in GR this is not the case anymore - an object moving only under the influence of gravity is said to be unaccelerated. This means that an object moving only under the influence of gravity, so called free-fall will follow a geodesic. That is why we say that in GR gravity is the manifestation of spacetime curvature - a collection of matter curves spacetime, therefore causing geodesics to also be curved. And the effect we perceive as gravity is simply caused by the fact that geodesics deviate from straight lines through space. Our photons are unaccelerated and therefore follow geodesics. 

In order to visualize the black hole we therefore need to trace geodesics in the Kerr spacetime, which we do by solving the *geodesic equation*

$$
\begin{equation}
  \label{eq: geodesic_equation} \tag{9}
  \frac{d^2 x^\mu}{d\lambda^2} + \Gamma^\mu_{\rho \sigma}\frac{dx^\rho}{d\lambda}\frac{dx^\sigma}{d\lambda} = 0,
\end{equation}
$$

where $x^\mu$ is the four-position of a photon, and $\lambda$ is an affine parameter for the geodesic, chosen such that $p^\mu \equiv \frac{dx^\mu}{d\lambda}$ is the four-momentum of the photon. $\Gamma^\mu_{\rho \sigma}$ are the *Christoffel symbols* of the metric. Equation $\eqref{eq: geodesic_equation}$ is completely general, and we of course need explicit expressions for the geodesic equation for each value of $\mu$. But computing the Christoffel symbols tends to be very tedious work, and certainly so for the Kerr metric. So instead of deriving the Christoffel symbols by hand we use the Sympy package in Python to derive them, and thus also explicit expressions for each component of the geodesic equation. Code for doing this can be found [here](https://github.com/chrvill/Geodesic_EOM_deriver).

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
  h_\mathrm{new} = 0.9 h \left(\frac{\epsilon}{TE}\right)^{1/5}
\end{equation}
$$

where $\epsilon$ is a tolerance value we can choose in order to specify the level of accuracy we want to achieve. Then, if $TE > \epsilon$ the error is too big and so we replace $h$ with $h_\mathrm{new}$ and repeat the step. We perform this iteration until $TE < \epsilon$. And then in the next step we use this value of $h_\mathrm{new}$ as the new $h$.

### Relativistic aberration

#### Tetrads

Before we explain aberration we will do a short detour explaining what *tetrads* are.

From the equivalence principle we know that we can always transform to a local inertial frame. And in contrast to Schwarzschild spacetime where shell observers are natural frames we will instead work with ZAMO frames, since these exist even inside the ergosphere. But the question is how we perform this transformation. How can we transform from the global BL coordinates to the local ZAMO coordinates? This is where the concept of a *vierbein*, also called a *frame field*, is extremely useful (An introduction to vierbeins can be found [here](https://jila.colorado.edu/~ajsh/courses/astr5770_21/text.html})). A vierbein is a set of orthonormal axes which form a local inertial frame. So we want the vierbeins to take us from the global coordinate system to the local one, which is to say from the metric $g_{\mu \nu}$ to the Minkowski metric $\eta_{\mu \nu}$. This means that

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

From above we know how to transform from a local ZAMO frame to the global frame using the ZAMO tetrads. So we will first have to transform from the instantaneous camera rest frame to the local ZAMO frame, in order to then use the tetrads from above to finally transform to the global frame. Let the camera have four-velocity $u^\mu$ in the global coordinate frame. Also note for later that this four-velocity will be expressed in BL coordinates. Imagine then a ZAMO frame at the same position as the camera. The four-velocity of the camera in the ZAMO tetrad frame is given by

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

in Cartesian coordinates. So we first have to convert the velocity vector $\vec{v}$ from BL coordinates to Cartesian coordinates. From the section on the coordinate transformation between Cartesian and BL coordinates we know that the components of the velocity vector in Cartesian components are given by

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

Having implemented this we have a fully general approach which takes into account relativistic aberration. And this is implemented purely through a coordinate transformation from the camera's rest frame to the global frame, which is in any case necessary if we want the camera to move freely around.

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

The relative amounts of the visible wavelengths that are emitted determine the color of the object. Unfortunately our eyes do not perceive and respond to all visible wavelengths of light equally. So converting an intensity spectrum $I_\lambda$, which represents the intensity distribution that is actually observed, to the color our eyes actually see is far from trivial. The procedure for calculating an RGB color comes in two main steps: converting a spectrum to so-called *tristimulus values*, and then converting the tristimulus values to RGB values. We found [5], [this](https://color.org/chardata/rgb/sRGB.pdf) and [this](https://en.wikipedia.org/wiki/SRGB) to be particularly helpful here. Assume now that we have a blackbody at some temperature $T$. We can find the tristimulus values $X, Y$ and $Z$ using the so-called *color matching functions* $\bar{x}, \bar{y}$ and $\bar{z}$, which describe how our eyes respond to different wavelengths. The $X, Y$ and $Z$ values are calculated in the following way

$$
\begin{align*}
  X &= \int_{\mathrm{380 nm}}^{\mathrm{780 nm}} I_\lambda \left(\lambda; T\right) \bar{x}\left(\lambda\right)d\lambda \label{eq: color_X} \tag{34} \\ \\
  Y &= \int_{\mathrm{380 nm}}^{\mathrm{780 nm}} I_\lambda \left(\lambda; T\right) \bar{y}\left(\lambda\right)d\lambda \label{eq: color_Y} \tag{35} \\ \\
  Z &= \int_{\mathrm{380 nm}}^{\mathrm{780 nm}} I_\lambda \left(\lambda; T\right) \bar{z}\left(\lambda\right)d\lambda. \label{eq: color_Z} \tag{36}
\end{align*}
$$

Approximations to these color matching functions can be found in [6]. But it turns out that at low temperatures these approximations are not accurate enough anymore. So we will instead be using the lookup table given in [7]. In the following figure we have plotted the color matching functions as functions of wavelength. The $x$-axis represents wavelength $\lambda$ in nm, while the $y$-axis represents the output of the color matching functions, which can be taken to have units such that $X, Y$ and $Z$ are dimensionless.

![The color matching functions](../images/miscellaneous/color_matching_functions.png)

While we have color-coded the different color matching functions, the three values $X, Y$ and $Z$ do not directly correspond to red, green and blue. Instead $X, Y$ and $Z$ live in an abstract color space.

### The chromaticity diagram

When we talk about colors we tend to group them into "classes" of colors that are similar, like greens, reds, blues etc. And we have an intuition that some of those colors are actually the same basic colors, just with different brightnesses, or more precisely *luminances*. And the quantity that is the same between different luminances is called the *chromaticity*.

To explain this in more detail we can study the *chromaticity diagram*. This is a 2d plot which shows the chromaticities corresponding to different points in $XYZ$ space. And the color matching functions give us a map between spectra and the $XYZ$ space. Consider now a spectrum consisting entirely of a single wavelength $\lambda_0$ - which is to say monochromatic light. We can represent this mathematically by $I_\lambda\left(\lambda\right) = \delta\left(\lambda - \lambda_0\right)$ (nevermind the units or the scaling), where $\delta\left(x - y\right)$ is the Dirac-delta "function". Plugging this into \eqref{eq: color_X} - \eqref{eq: color_Z} gives

$$
\begin{align*}
  X &= \bar{x}\left(\lambda_0\right) \\ \\
  Y &= \bar{y}\left(\lambda_0\right) \\ \\
  Z &= \bar{z}\left(\lambda_0\right).
\end{align*}
$$

So if we let $\lambda_0$ vary over the visible range we can trace out the curve in $XYZ$ space made by the wavelengths in the visible spectrum. But in order to visualize this space in 2d we define the new coordinates $x$ and $y$ through

$$
\begin{align*}
  x &= \frac{X}{X + Y + Z} \label{eq: x_cie1931} \tag{37}\\ \\
  y &= \frac{Y}{X + Y + Z} \label{eq: y_cie1931} \tag{38}.
\end{align*}
$$

It also naturally follows that we can define $z = \frac{Z}{X + Y + Z} = 1 - x - y$. We can now plot the $(x, y)$ coordinates of the monochromatic spectrum as we vary $\lambda_0$, and the resulting curve is called the *spectral locus*, and is shown in the following figure. The marked points show where different wavelengths fall on the spectral locus.

![The spectral locus](../images/miscellaneous/chromaticity_diagram_without_planckian_locus.png)

This spectral locus is the boundary of the chromaticity diagram mentioned earlier, and which is shown in [this](https://en.wikipedia.org/wiki/CIE_1931_color_space) Wikipedia article. Here we have colored in the different points along the spectral locus according to which RGB color they correspond to. But we have not yet explained how we can compute these RGB colors from the $x$ and $y$ values, we get to that later.

We can also calculate the curve in the chromaticity diagram which represents the color of blackbodies at different temperatures. This just involves calculating $X, Y$ and $Z$ from the spectrum $I_\lambda\left(\lambda; T\right)$ for varying $T$. The resulting curve in $(x, y)$ coordinates is called the *Planckian locus* (or equivalently the *blackbody locus*). Plotted together with the spectral locus the result we get is shown in the following figure where the thick curve represents the Planckian locus. We have again calculated the RGB color for each point along the Planckian locus.

![The Planckian locus](../images/miscellaneous/chromaticity_diagram_with_planckian.png)

### Transformation from $XYZ$ to RGB

There are many different RGB color spaces designed for outputting colors to a screen, and these different RGB color spaces are defined by their so-called *primaries*. When choosing an RGB color space we choose which primaries red, green and blue to use. These primaries correspond to the vertices in the colored triangle in the following figure. The set of points inside the triangle between the primaries is called the *gamut* of the RGB space, and makes up the colors that the given RGB space can represent. And since different RGB spaces have different primaries the sets of colors that they can represent are therefore also different.

![Chromaticity diagram with sRGB gamut](../images/miscellaneous/SRGB_chromaticity_CIE1931.svg){: width="700" height="400" background-color=white}

The points inside the gamut are linear combinations of the primaries, so the primaries form a basis for the RGB space. Red $R$, green $G$ and blue $B$, being the primaries, are of course written as

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
             \end{matrix}\right].
\end{equation}
$$

in this basis. We will also have an XYZ basis, in which the RGB primaries can be written as

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
             \end{matrix}\right].
\end{equation}
$$

The problem of finding the RGB values for some given set of $XYZ$ values then translates into finding the coordinate transformation between the RGB coordinates and $XYZ$ coordinates. We can write the transformation as

$$
\begin{equation}
  \label{eq: general_XYZ_to_RGB_relation} \tag{39}
  \left[\begin{matrix}
          R \\
          G \\
          B
        \end{matrix}\right] = M\left[\begin{matrix}
                                      X \\
                                      Y \\
                                      Z
                                     \end{matrix}\right]
\end{equation}
$$

where $M$ is the transformation matrix we want to find. The standard way of specifying the primaries of an RGB color space is by giving the $x, y$ and $Y$ values of the primaries. So in order to work with the $X, Y$ and $Z$ values we first have to calculate the unknown $X$ and $Z$ values. From \eqref{eq: x_cie1931} we see that

$$
\begin{equation}
    X = \left(X + Y + Z\right)x = \frac{x}{y}Y.
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
                                     \end{matrix}\right] \equiv MA
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
                                    \end{matrix}\right].
\end{equation}
$$

With these values we get

$$
\begin{equation}
  M = \left[\begin{matrix}
              3.24156456 & -1.53766524 & -0.49870224 \\
              -0.96920119 & 1.87588535 & 0.04155324 \\
              0.05562416 & -0.20395525 & 1.05685902
            \end{matrix}\right]
\end{equation}
$$

It is also typically advised to perform a final non-linear correction after this linear transformation. But our results looked better without this non-linear correction, so we chose to omit it.

Remember that the reason we wanted to find the transformation from $XYZ$ to RGB was that we know how to calculate the $XYZ$ values corresponding to a given spectrum. So with this coordinate transformation in hand we can convert a spectrum into an RGB color. This is of interest to us here because we know the temperature $T$ and redshift $(1 + z)$ of points on the accretion disk. So we can calculate the blackbody spectrum of each point and convert that to an RGB color to display on screen. Incorporating redshift into this is actually very easy. Consider a blackbody spectrum redshifted such that $\lambda_\mathrm{shifted} = \left(1 + z\right)\lambda$. Then

$$
\begin{align*}
    I_{\lambda}\left(\lambda_\mathrm{shifted}; T\right) &= \frac{2h c^2 \left(1 + z\right)^5}{\lambda_\mathrm{shifted}^5} \frac{1}{e^{hc \left(1 + z\right)/\lambda_\mathrm{shifted} k_B T} - 1} \\ \\
    &= \frac{2hc^2 \left(1 + z\right)^5}{\lambda_\mathrm{shifted}^5} \frac{1}{e^{hc/\lambda_\mathrm{shifted} k_B T_\mathrm{shifted}} - 1}, \label{eq: redshifted_spectrum} \tag{40}
\end{align*}
$$

with $T_\mathrm{shifted} \equiv \frac{T}{1 + z}$. Notice now that this is also a blackbody spectrum, evaluated at a "redshifted temperature". And we have an additional vertical scaling by $(1 + z)^5$. From setion earlier we know that relativistic beaming modifies the observed intensity by a factor $(1 + z)^{-5}$, and when applying this to \eqref{eq: redshifted_spectrum} these factors of $(1 + z)$ cancel out. This means that we can actually just use the unmodified blackbody spectrum $I_\lambda$, just with a temperature $T_\mathrm{shifted} = \frac{T}{1 + z}$. So this accounts for relativistic beaming.

We can plot the map showing the blackbody color for different combinations of $T$ and $1 + z$. This is shown in the following figure for temperatures $T \in \left[200, 10000\right]$ K and for $(1 + z) \in \left[0.1, 2\right]$.

![Temperature-redshift map](../images/miscellaneous/temp_redshift_map.png)

## References

[1] Sean M. Carroll. Spacetime and Geometry: An Introduction to General Relativity. Cambridge University Press, 2019

[2] Masud Mansuripur. An exact derivation of the thomas precession rate using the lorentz transformation. In Henri-Jean M. Drouhin, Jean-Eric Wegrowe, and Manijeh Razeghi, editors, Spintronics XIII. SPIE, August 2020

[3] Z. Kh. Kurmakaev. Circular orbits in the Kerr metric. , 18:110, August 1974

[4] Eric Bruneton. Real-time high-quality rendering of non-rotating black holes, 2020

[5] Philippe Colantoni. Color space transformations. 2006

[6] Chris Wyman et al. Simple analytic approximations to the cie xyz color matching functions. 2013

[7] Kasajima Ichiro. Plotting colors on color circle: Interconversion between xyz values and rgb color system. Current Trends in Analytical and Bioanalytical Chemistry, 1, 05 2017
