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

Black holes are product of the theory of *general relativity* (GR), so we need to discuss GR in order to describe how one can create a visualization of a black hole. GR is infamously hard, so we will not do the subject justice in this post, but will have to include quite an extensive amount of math nonetheless. 

### 2.1 Special relativity

Special relativity is based on two postulates:

1. The laws of physics take the same form in all inertial reference frames.
2. The speed of light measured in all inertial frames is always the same

A consequence of these postulates is that there exists no inertial frames moving at the speed of light relative to any other inertial frames. Because then light would be at rest in such a frame, and so the second postulate would be broken. This is to say that questions such as "what would you see moving at the speed of light?" or "what does a photon measure?" are fundamentally flawed, because there does not exist any frames moving at the speed of light. So relativity cannot answer these questions. 

Central to both GR and SR is the concept of *spacetime*. Spacetime can be thought of as a combination of the familiar 3 spatial dimensions and the 1 dimension of time, commonly described as a (3 + 1)-dimensional spacetime. 

Have to add more here....



### 2.2 General relativity and the metric 

General relativity is like the big brother of special relativity. One can describe special relativity with the exact same mathematical framework as is used in general relativity, but in special relativity the math all reduces to much simpler forms. This is because special relativity assumes that spacetime has no matter- or energy-contents, while general relativity allows an arbitrary matter- and energy-distribution. 

General relativity is grounded in very large part on the *equivalence principle*. One formulation of this principle is that *locally* the laws of physics reduce to those of special relativity, and that it is impossible to detect the existence of a gravitational field locally. Locality in this context means that we are considering infinitesimal regions of spacetime, which is to say infinitesimal intervals of time and/or infinitesimal separations in space. So said in a simpler way, if you only consider a small region around yourself and a short interval of time you are fine using only special relativity. And in this small region of space and short interval of time you will not feel gravity. In practice this means that whenever some observer wants to measure something local to them they can use special relativity instead of having to use the full machinery of general relativity. 

Spacetime in special relativity is called *Minkowski spacetime*, which is also called *flat* spacetime. This is spacetime in the absence of gravity. While general relativity deals with the general case (naturally) where spacetime is allowed to have *curvature*. This curvature is set up by the matter- and energy-content which is present, and the curvature and matter/energy content are related through the *Einstein field equations* (roughly). 

In everyday life we are used to the fact that objects move along straight lines through space when they are not accelerating, that we can have lines which are parallel everywhere, and that distances between points can be calculated using the Pythagorean theorem. However, these phenomena and this intuition come from special relativity and do not in general hold in curved spacetime. In particular, the Pythagorean theorem no longer holds in curved spacetime - we instead need a generalization that describes curved spacetime. This is where the *metric* becomes important. Essentially the metric is a matrix that tells you how to compute distances between points in spacetime, and it encodes the geometry of spacetime. 

As an example we can take the metric for Minkowski spacetime, which is to say the metric describing special relativity. The following is the metric (also called the line element)

$$
\begin{equation}
    \label{eq: metric_minkowski}
    ds^2 = -dt^2 + dx^2 + dy^2 + dz^2
\end{equation}
$$

where $dt$ is an infinitesimal interval of time and $dx, dy$ and $dz$ are infinitesimal distanes along the $x$, $y$ and $z$ directions. Notice that if we take $dt = 0$ this reduces to the Pythagorean theorem. But with $dt \neq 0$ this tells us that if we want to compute the spacetime distance between two points in spacetime the contribution to the distance from the infinitesimal interval of time is different from the contributions from the infinitesimal spatial increments, due to the minus sign in front of the $dt^2$-term. This metric is written in *Cartesian coordinates* $x, y, z$. But a fundamental principle of general relativity is that the coordinate system one chooses to express physical quantities in is completely irrelevant - the physics should be the same no matter which coordinate system you use. So let us move to spherical coordinates $r, \theta, \phi$. These are related to the Cartesian coordinates $x, y, z$ through

$$
\begin{align*}
    x &= r\sin\left(\theta\right)\cos\left(\phi\right) \\
    y &= r\sin\left(\theta\right)\sin\left(\phi\right) \\ 
    z &= r\cos\left(\theta\right)
\end{align*}
$$

Now we want to rewrite the metric $ds^2$ in terms of $dr, d\theta$ and $d\phi$. All these quantities $dx, dy, dz$ and $dr, d\theta, d\phi$ are so-called *differentials*, so in order to do the switch $dx, dy, dz \to dr, d\theta, d\phi$ we have to use calculus. But omitting all the unnecessary details here the metric we get in spherical coordinates is 

$$
\begin{equation}
    ds^2 = -dt^2 + dr^2 + r^2\left(d\theta^2 + \sin^2\left(\theta\right) d\phi^2\right)
\end{equation}
$$

The point of rewriting the metric to spherical coordinates was just to introduce spherical coordinates, and to explain that the coordinates that we use are irrelevant to the physics. We are completely free to choose whatever coordinates we want, and furthermore the coordinates do not even need to have a clear physical interpretation. 

#### 2.2.1 The Schwarzschild solution

We mentioned before that the curvature of the spacetime and its matter- and energy-contents are related through the Einstein field equations. That is a slight oversimplification, and it's more accurate to say that the *metric* is related to the matter- and energy-contents through the Einstein field equations. And the curvature is then related to the metric. So we have this equation which tells us that, given an matter/energy distribution, spits out the form that the metric for that spacetime must take. However, the process of actually finding the solution for the metric given the matter/energy distribution is very complicated, because the Einstein field equations are extremely nasty. 

The only hope we have of finding analytical solutions to the Einstein field equations is to assume that the metric, which we're solving for, has some sort of symmetry. Symmetries in physics are extremely important, and roughly describe the fact that under certain circumstances some types of transformations can leave our equations unchanged. For example, let's say you are moving through completely empty space. Then it doesn't matter whether you're 1 meter to the right or 1 meter to the left, there's no reason to expect your motion to be any different in empty space. This is a symmetry of empty space - the fact that you can just add some arbitrary vector to your current position, and it doesn't change the physics. As a sidenote, this is called *translational symmetry* and leads to conservation of momentum. If this is interesting look up *Noether's theorem*. 

As another example assume you're moving in an orbit around the Earth. Now, the gravitational potential energy is only dependent on your distance from the center of the Earth. So in this case it shouldn't matter if you rotated your orbit by say $90^\circ$ - the energy is still the same. In this case we have *rotational symmetry*, also called *axial symmetry*, and this symmetry leads to conservation of angular momentum.

Now comes the actual relevant point to us here. Let us say that we want to model the spacetime around some massive object. What assumptions can we make here that simplify the math? Or put another way, which symmetries does this system possess? The key thing to notice here is that it doesn't matter what angle we view this system from - it looks the same from all angles. The only thing that matters is the distance from the central massive object. This is to say that we have *spherical symmetry*. Another assumption we can make is that the massive object is the same at all times, which is to say it's *static*. These two assumptions, spherical symmetry and a static spacetime, allow us to drastically simplify the possible form that the metric can take. The only thing we still have to resolve here is that even with these two assumptions there is a free parameter in the metric. But the key to fixing that parameter in place is that we know how the gravitational field should look around a massive object in Newtonian physics, which corresponds to the low-velocity, weak-gravity regime. Evaluating the form of the metric in the low-velocity, weak-gravity limit fixes the free parameter, and the free parameter can be identified as the *mass* of the central object. The final form of the metric is 

$$
\begin{equation}
    \label{eq: schwarzschild_metric} \tag{2}
    ds^2 = -\left(1 - \frac{2M}{r}\right)dt^2 + \left(1 - \frac{2M}{r}\right)^{-1} dr^2 + r^2 \left(d\theta^2 + \sin^2\left(\theta\right) d\phi^2\right),
\end{equation}
$$

where $M$ is the mass of the central object. This is called the *Schwarzschild metric*. You might notice that the last parenthesis is the same as in the Minkowski metric. And that's because both the Schwarzschild metric and the Minkowski metric exhibit spherical symmetry. The Schwarzschild metric describes the spacetime around a *non-rotating*, *charge-neutral* black hole. 

#### 2.2.2 The Kerr solution

The Schwarzschild metric is very useful, particularly in analytical work, because it's simple enough that we can do many calculations by hand. But it's still complicated enough to not be trivial and to exhibit some interesting phenomena. So we could've used this metric to visualize a non-rotating black hole. But we wanted to go a bit further and look at a *rotating* black hole instead. Consider now again what assumptions we can make about the rotating black hole which simplifies the math. We can no longer assume spherical symmetry, because the black hole is rotating around some axis, which breaks this symmetry. But the spacetime will be rotationally symmetric around the axis that the black hole rotates around. The spacetime is also not static anymore, because the black hole is rotating. But it's rotating in the same exact same way all the time, which means the spacetime is *stationary*. Exactly what this means isn't that important here, but the point is that even with these weaker assumptions there still exists an analytical solution. This is the *Kerr metric* given by

$$
\begin{equation}
    \label{eq: Kerr_metric} \tag{3}
    ds^2 = g_{tt}dt^2 + g_{rr}dr^2 + g_{\theta\theta}d\theta^2 + g_{\phi\phi}d\phi^2 + g_{t\phi}\left(dtd\phi + d\phi dt\right),
\end{equation}
$$

with 

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
    x &= \sqrt{r^2 + a^2}\sin\theta\cos\phi, \label{eq: x_BL} \tag{4} \\
    y &= \sqrt{r^2 + a^2}\sin\theta\sin\phi, \label{eq: y_BL} \tag{5} \\
    z &= r\cos\theta. \label{eq: z_BL} \tag{6}
\end{align*}
$$

These BL coordinates are similar to the spherical coordinates we encountered earlier, but they aren't exactly the same. This metric is the one we'll use throughout. 

### 2.3 Geodesics

In normal raymarching one typically just casts straight lines to model how light moves in flat spacetime. But in curved spacetime light no longer follows these straight paths through space. Instead it follows the curved-space generalization of "straight paths", namely *geodesics*. In "normal" spaces (meaning not spacetime, simply space. The proper terminology is *Euclidean* space) a geodesic can be thought of as the path which minimizes the length between two points $A$ and $B$. In Euclidean space (think of a plane, for example) this just reduces to a straight line as we expect. But we can get more interesting geodesics if we consider a sphere, for example. This is shown in the following figure, where the red curve shows the geodesic between two points $A$ and $B$ on the sphere. 

![Sphere geodesic](../images/miscellaneous/A-geodesic-on-the-surface-of-a-sphere.png)

This red line, the geodesic, is the shortest possible path between $A$ and $B$ if we are constrained to move on the spherical surface. You might also notice that the geodesic lies along a line of constant longitude, meaning that lines of constant longitude represent "straight" lines on the sphere. What this concept of a path being "straight" here really means is simply that you can follow such a path without having to steer - you can just move in your forwards direction. It might not be obvious, however, that lines of constant latitude are *not* geodesics on the sphere - they are not "straight". This is to say that if you were constrained to the surface of a sphere you could move along lines of constant longitude without steering, but not along lines of constant latitude. The concept of geodesics is the reason why the flight path of a plane can look unreasonable when projected into a two dimensional map, even though it is actually the shortest path.

The mathematics behind geodesics carries over to curved spacetime, but the intuition is harder to grasp. Rather than being the path of shortest spatial distance between two points it is the path of *longest proper time*, which is the time measured by an observer following the path. Let us assume we have two observers, which we creatively call Alice and Bob, moving along two different paths from point $A$ to point $B$ in spacetime. Alice moves along a geodesic between these points, while Bob does not. This then means that Alice will measure a longer time from $A$ to $B$ on her clock than Bob will measure on his clock. 

In relativity the geodesic also has the interpretation of being the path an unaccelerated particle follows. In our previous analogy then, Alice is in unaccelerated free-fall, while Bob accelerates in some way. Now, a possible cause for confusion is the fact that in classical physics we tend to say that an object accelerates due to gravity. But in GR this is not the case anymore - an object moving only under the influence of gravity is said to be unaccelerated. This means that an object moving only under the influence of gravity, so called free-fall will follow a geodesic. That is why we say that in GR gravity is the manifestation of spacetime curvature - a collection of matter curves spacetime, therefore causing geodesics to also be curved. And the effect we perceive as gravity is simply caused by the fact that geodesics deviate from straight lines through space. Our photons are unaccelerated and therefore follow geodesics. 
___

## 3 Tracing geodesics in curved spacetime

### 3.1 Equations of motion

#### 3.1.1 Example of Newtonian gravity

It's instructional to consider *classical mechanics*, and in particular *Newtonian gravity*, and to briefly discuss the equations of motion in that framework. The Newtonian picture of classical mechanics is based on Newton's laws, the 2nd of which being

$$
\begin{equation}
    \sum \vec{F} = m \vec{a},
\end{equation}
$$

with $\sum \vec{F}$ representing the sum of all the forces acting on a body with mass $m$, and $\vec{a}$ being the acceleration of that body. Newton also gave us his famous law of gravity 

$$
\begin{equation}
    \vec{F}_g = -\frac{Gm_1 m_2}{r^2} \vec{\hat{r}}.
\end{equation}
$$

Here $\vec{F}_g$ represents the gravitational force between two bodies with masses $m_1$ and $m_2$ a distance $r$ from each other. $\vec{\hat{r}}$ is the unit vector pointing from one of the masses towards the other and $G$ is the gravitational constant. Consider now the two-body problem, in which case the only force acting on the two bodies is the gravitational force between them. Then 

$$
\begin{align*}
    \sum \vec{F}_1 &= -\frac{G m_1 m_2}{r^2} \vec{\hat{r}} = m_1 \vec{a}_1 \\ \\
    \sum \vec{F}_2 &= \frac{G m_1 m_2}{r^2} \vec{\hat{r}} = m_2 \vec{a}_2
\end{align*}
$$

$\vec{a}_1$ and $\vec{a}_2$ are the accelerations of body 1 and 2 respectively. The acceleration is the second time derivative of the position, so these two equations can be written as 

$$
\begin{align*}
    \frac{d^2 \vec{r}_1}{dt^2} &= -\frac{G m_2}{r^2} \vec{\hat{r}} \label{eq: Newtonian_EOM_1} \tag{7} \\ \\
    \frac{d^2 \vec{r}_2}{dt^2} &= \frac{G m_1}{r^2} \vec{\hat{r}} \label{eq: Newtonian_EOM_2} \tag{8}
\end{align*}
$$

where $\vec{r}_1$ and $\vec{r}_2$ are the position vectors of body 1 and 2 respectively. These two equations are coupled differential equations (DEs), and are called the *equations of motion* for the two-body system. Solving a system's equations of motion tells us how the system evolves over time - in this specific example a solution to the equations of motion tells us how the two bodies move in space. The two-body problem is a particularly simple problem where there actually exists a closed form analytical solution. However, this is not the case for most physical systems. In general we have to use numerical methods to solve the equations of motion. In order to see how this can be done we rewrite the two *second order* DEs above as a set of four *first order* DEs:

$$
\begin{align*}
    \frac{d\vec{r}_1}{dt} &= \vec{v}_1 \label{eq: EOM_1} \tag{9} \\ \\
    \frac{d\vec{v}_1}{dt} &= -\frac{G m_2}{r^2} \vec{\hat{r}} \label{eq: EOM_2} \tag{10} \\ \\
    \frac{d\vec{r}_2}{dt} &= \vec{v}_2 \label{eq: EOM_3} \tag{11} \\ \\
    \frac{d\vec{v}_2}{dt} &= \frac{G m_1}{r^2} \vec{\hat{r}} \label{eq: EOM_4} \tag{12}
\end{align*}
$$

In $\eqref{eq: EOM_1}$ and $\eqref{eq: EOM_3}$ we just define the velocities $\vec{v}_1$ and $\vec{v}_2$ of bodies 1 and 2 as the time derivatives of their positions. And in equations $\eqref{eq: EOM_2}$ and $\eqref{eq: EOM_4}$ we have just plugged these definitions into $\eqref{eq: Newtonian_EOM_1}$ and $\eqref{eq: Newtonian_EOM_2}$. This system of four coupled, first order DEs can be solved fairly easily numerically, and we will go into details on how this can be done later. This procedure of transforming $n$ second order DEs into $2n$ first order DEs is completely general, one just has to use the appropriate right hand side in $\eqref{eq: EOM_1} - \eqref{eq: EOM_}$. 

#### 3.1.2 The geodesic equation 

In order to form the equations of motion for particles moving only under the influence of gravity in classical mechanics we used Newton's 2nd law together with Newton's law of gravity. That gave us the DEs we had to solve - the equations of motion. Here in this project we want to find and solve the equations of motion for particles moving along geodesics in the Kerr spacetime. For that we need the *geodesic equation*

$$
\begin{equation}
    \label{eq: geodesic_eq_general} \tag{13}
    \frac{d^2 x^\mu}{d\lambda^2} + \Gamma^{\mu}_{\rho \sigma} \frac{dx^\rho}{d\lambda} \frac{dx^\sigma}{d\lambda} = 0.
\end{equation}
$$

Here $x^\mu$ is the so-called *four-position* describing the position of some particle in spacetime. It is essentially a vector with four components, and can be written as $x^\mu = \left(t, x, y, z\right)$ in Cartesian coordinates. So it specifies the spatial position of the particle $(x, y, z)$ along with the point in time $t$. $\lambda$ is a so-called *affine parameter* for the geodesic, and can be thought of as a parameter that parametrizes the geodesic - so by varying $\lambda$ you can trace out the geodesic. $\Gamma^\mu_{\rho \sigma}$ are the so-called *Christoffel symbols* of the metric. These encode the geometry at each point in spacetime, and is what allows the curvature of spacetime to influence the path that a freely falling particle falls along. You might notice now that the first term in $\eqref{eq: geodesic_eq_general}$ looks very similar to the terms on the left-hand side of $\eqref{eq: Newtonian_EOM_1}$ and $\eqref{eq: Newtonian_EOM_2}$. However, an aspect that makes the geodesic equation substantially more complicated is that there are implicit sums over $\rho$ and $\sigma$. This means that the geodesic equation is really 

$$
\begin{equation}
    \label{eq: geodesic_eq_explicit_sums} \tag{14}
    \frac{d^2 x^\mu}{d\lambda^2} + \sum_{\rho} \sum_{\sigma} \Gamma^\mu_{\rho \sigma} \frac{dx^\rho}{d\lambda} \frac{dx^\sigma}{d\lambda} = 0,
\end{equation}
$$

where $\rho$ and $\sigma$ take on the values $0, 1, 2, 3$, which is to say that it is an index. So we will for example have a term like $\Gamma^\mu_{00} \frac{dx^0}{d\lambda} \frac{dx^0}{d\lambda}$ for $\rho = 0$, $\sigma = 0$, for example. And the superscript $\mu$ can also take on the values $0, 1, 2, 3$, however there is not a sum over it. Instead it tells you which component of $\frac{d^2 x^\mu}{d\lambda^2}$ we are talking about. For example for $\mu = 0$ we get an equation telling us how $x^0$ evolves. There are 4 values for all the indices, which means that we have 4 different equations of motion, one for each of the components of $x^\mu = (t, x, y, z)$. So the equation for $x^0$, for $\mu = 0$, tells us how $t$ evolves over time. And likewise for the three other values. The double sum in $\eqref{eq: geodesic_eq_explicit_sums}$ therefore have $4^2 = 16$ terms

The Christoffel symbols are given by 

$$
\begin{equation}
    \label{eq: Christoffel_symbol_definition} \tag{14}
    \Gamma^\mu_{\rho \sigma} = \frac{1}{2} g^{\mu \nu} \left(\partial_\rho g_{\sigma \nu} + \partial_\sigma g_{\rho \nu} - \partial_\nu g_{\rho \sigma}\right)
\end{equation}
$$

where $\partial_\mu = \frac{\partial}{\partial x^\mu}$ is the partial derivative with respect to the coordinate $x^\mu$. 

### 3.2 Initial conditions

### 3.3 Flat/Schwarzschild/Kerr example (+Integration schemes)

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

### 5.4 Rest frame to ZAMO to global transformation

#### 5.4.1 Relativistic aberration

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