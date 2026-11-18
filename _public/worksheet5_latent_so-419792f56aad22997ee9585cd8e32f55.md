---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.18.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(worksheet5_latent_solution)=
# Worksheet 5 solution: Latent heat

This is a pencil and paper problem asking you to show that $l_v = h_v - h_l$

+++

## Question

Starting with the first law of thermodynamics, use the definition of the enthalphy $h = u + p \alpha$ to show that for liquid water and water vapor in equilibrium (i.e. constant pressure and temperature) the latent heat is just the difference between the vapor and liquid water enthalpies.

Hint:  This is mostly done by Thompkins on p. 26.


Hand in a scanned pdf with your solution.

## Answer

First law:

$$
du = q\,dt - w\,dt = q\,dt - p\,dv
$$(eq:first)

+++

Definition of enthalpy:

$$
h = u + pv \\
dh = du + p\,dv + v\,dp
$$(eq:enthalpy)

+++

At constant pressure $dp = 0$ so insert {eq}`eq:enthalpy` into {eq}`eq:first` and integrate between the transition from liquid to vapor phases over the time period $t1$ to $t2$ using total integrated heating Q

+++

$$
\int_l^v dh = h_v - h_l = \int_{t1}^{t2} q\,dt = Q = l_v
$$

```{code-cell} ipython3

```
