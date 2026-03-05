---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.1
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
  language: python
---

(assign6_koehler)=
# Assignment 6

## Notebook part

- Hand in {ref}`worksheet8_interp1d`
- Add a plot of vertical velocity vs. CAPE, borrowing the code from {ref}`cape_part1`.  Does the plot agree with Thompkins equation 3.4?
  $$
  w_{max} = \sqrt{2 \,\text{CAPE}}
  $$

## On paper


+++


### Problem 1: Taylor series

Show using two Taylor series expansions (one for $a_w$ and one for
$\exp \left [ \frac{2\sigma}{\rho_l R_v T r} \right ]$) show how to get  Thompkins 4.15:

$S= a_w \exp \left [ \frac{2\sigma}{\rho_l R_v T r} \right ] \approx \left ( 1 + \frac{a}{r} - \frac{b}{r^3} \right )$

+++

### Problem 2: Term comparison  

For the aerosol defined in the kohler.ipynb notebook ($10^{-19}$ kg, ammonium sulphate), inside a droplet of radius $r=1\ \mu m$
      find the size of the smallest term you've kept (either $\frac{a}{r}$ or $\frac{b}{r^3}$ and compare it to
      the size of the largest term you've dropped.   Repeat this for a droplet of radius   $r=0.1\ \mu m$

+++

### Problem 3: surface energy:

Suppose you have $r_l$ =1 g/kg of liquid water spread among N spherical drops of radius 10 $\mu m$.  Compare the surface energy of this
      population (which we've been neglecting) with the enthalpy required to vaporize it:  $l_v r_l$.  Is it negligible in comparison?

```{code-cell} ipython3

```
