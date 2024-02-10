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

Our simulation takes the following relativistic effects into account:

* Redshift/blueshift
* Relativistic beaming
* Light travel time delay
* Relativistic aberration

In addition to this the camera is allowed to follow any arbitrary geodesic, and the effects of redshift and aberration become very apparent when the camera is moving. In roder to get realistic colors, also taking into account how redshift and blueshift affects the colors, we had to convert a blackbody spectrum to RGB values. This is a highly non-trivial task, and is therefore explained in detail. But before getting into the details we need to introduce notational choices and also the Kerr metric itself.

### Notation and terminology

* We use the *Einstein notation convention*, so that $\sum_\mu a^\mu b_\mu \rightarrow a^\mu b_\mu$

The aim of this post is to describe in detail how one can create a visualization of a rotating black hole. There are of course many resources online describing how one can do this, but we generally found these to lack good explanations. Here we want to remedy this. Naturally, since black holes are consequences of general relativity, you need a fairly solid understanding of general relativity in order to understand the derivations of the equations we use. But the key *equations of motion* we use are just standard differential equations. So an understanding of how to solve differential equations numerically might suffice.
