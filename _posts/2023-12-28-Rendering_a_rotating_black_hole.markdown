---
title: Rendering a rotating black hole
date: 2023-12-28
categories: [physics, general relativity, black holes]
tags: [black holes]     # TAG names should always be lowercase
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
