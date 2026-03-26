---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(assignment_precip)=
# Assignment 8 -- precipitation.

Downlad link for [assignment8_precipitation.ipynb](https://drive.google.com/file/d/1grmUECBF4ePSNUO_264ibG546GiYKk9y/view?usp=sharing)

## Problem 1 -- Cloud Condensation Nuclei Counter

Hand in the {ref}`worksheet10_problem`

## Precipitation problems

For each of the 3 problems below:  First solve them analytically by hand
then check that solution with python using odeint.  In your notebook compare your analytic equation with the numerical solution.


### Problem 2 - Thompkins


A drop with an initial radius of 100 µm falls through a cloud containing 100 droplets per cubic centimeter that it collects in a continuous manner with a collection efficiency of 0.800. If all the cloud droplets have a radius of 10 µm, how long will it take for the drop to reach a radius of 1 mm? You may assume that for the drops of the size considered in this problem the terminal fall speed v (in $m s^{-1}$) of a drop of radius r (in meters) is given by $v= 8 x 10^3\;r$. Assume that the cloud droplets are stationary and that the updraft velocity in the cloud is negligible.  Hint:  Integrate Thompkins equation 4.28 analytically -- you can also check numerically with python

+++

### Problem 3  -- Wallace and Hobbs problem 6.21  -- updraft

Derive an expression for the height h above cloud base of a droplet at time t that is growing by condensation only in a cloud with a steady updraft velocity w and supersaturation S. 

That is -- integrate

$$
\frac{dh}{dt} = w - \nu
$$
where h is the height, w is the updraft speed and \nu is the terminal droplet fall speed
Use WH equation 6.24 for the terminal fall speed of a droplet together with WH eq 6.21 for the approximate droplet growth ratge

### Problem 4 -- Wallace and Hobbs problem 6.24 -- falling precip

If a raindrop has a radius of 1 mm at cloud base, which is located 5 km above the ground, what will be its radius at the ground and how long will it take to reach the ground if the relative humidity between cloud base and ground is constant at 60%? [Hint: Use (6.21) and the relationship between v and r given in Exercise 6.23. If r is in micrometers, the value of Gl in (6.21) is 100 for cloud droplets, but for the large drop sizes considered in this problem the value of Gl should be taken as 700 to allow for ventilation effects.]
