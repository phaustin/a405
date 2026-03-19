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

(assign7_aerosols)=
# Assignment 7: aerosols

-  Due Friday March 20 midnight

+++

## Problem 1: cloud chamber

+++

1. Given the critical supersaturation from the kohler notes:

    $$
    SS=S^* - 1= \left ( \frac{4 a^3}{27b} \right )^{1/2}
    $$

show that this implies, for $(NH_4)_2 SO_4$, density $\rho_{aer}$ = 1775
${kg}\,{m^{-3}}$ , van hoft i=3, that:

 $$
S^* -1 \approx 1.54 \times 10^{-12}~ m_{aer}^{-0.5}
 $$

where $m_{aer}$ is the ammonium sulphate aerosol mass in kg.

Note that this is why a cloud chamber can get the aerosol mass distribution from a series of
saturation and light scattering measurements as smaller and smaller aerosols are pushed over
their critical supersaturation and activated.

+++

## Problem 2: Koehler stability

+++

2. Show that the expression for second derivative of the thermodynamic potential derived in the  kohler stability notes:

     $$
     \frac{\delta ^2G}{\delta r^2} = - 4 \pi R_v T \rho_l \left [ 2 a - r^2 \left ( 1 +
           \frac{b}{r^3} \right ) \frac{3b}{r^4}  \right ] + 8 \pi \sigma
     $$

     Changes sign from stable (positive) to unstable (negative) at $r_{crit}$.

     Hint -- first show that the second derivative is zero at the critical radius.  Then show that
     the third derivative is negative above and below the critical radius, which means that
     there has to be a sign change from + to -.

+++

## Problem 3: CCN lifetime

+++

3. Assuming that cloud condensation nuclei (CCN) are removed from the atmosphere by first serving as the centers on which cloud droplets form, and the droplets subsequently grow to form precipitation particles, estimate the residence time of a CCN in a column extending from the surface of the Earth to an altitude of 5 km. Assume that the annual rainfall is 100 cm/year and the cloud liquid water content is 0.30 $g/kg$ .  *Hint:  Assume that all drops in the cloud droplets have  a radii of 10 microns and that every droplet contains exactly 1 CCN.   How many CCN are in 1 kg of air?  About how many kg of air are there in a 5 km column?  About how many CCN are taken out by a rain rate of 1 m/year?  Find the time constant for removal of the form  1/N dN/dt = 1/tau*
