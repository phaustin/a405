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

(assignment_precip_solution)=
# Assignment 8 -- precipitation solution

Downlad link for [assignment8_precipitation.ipynb](https://drive.google.com/file/d/1grmUECBF4ePSNUO_264ibG546GiYKk9y/view?usp=sharing)

## Problem 1 -- Cloud Condensation Nuclei Counter

Hand in the {ref}`worksheet10_problem`

## Precipitation problems

For each of the 3 problems below:  First solve them analytically by hand
then check that solution with python using odeint.  In your notebook compare your analytic equation with the numerical solution.


### Problem 2 - Thompkins

A drop with an initial radius of 100 µm falls through a cloud containing 100 droplets per cubic centimeter that it collects in a continuous manner with a collection efficiency of 0.800. If all the cloud droplets have a radius of 10 µm, how long will it take for the drop to reach a radius of 1 mm? You may assume that for the drops of the size considered in this problem the terminal fall speed v (in $m s^{-1}$) of a drop of radius r (in meters) is given by $v= 8 x 10^3\;r$. Assume that the cloud droplets are stationary and that the updraft velocity in the cloud is negligible.  Hint:  Integrate Thompkins equation 4.28 analytically -- you can also check numerically with python

+++

#### Problem 2 Analytic solution

Equation 4.30 (be sure you can derive this):

$$
  \frac{ dR}{dt}  = \frac{L V \epsilon E}{4\rho_l}
$$(eq:coal)
where $R$ is the radius of the collector drop in meters, $L$ is the liquid water content in $kg\,m^{-3}$, $\rho_l$ is the density of liquid water in $kg\,m^{-3}$, $V=8000R$ is the differential fall speed in m/s, $\epsilon=1$ is the coalescence efficiency and $E=0.8$ is the collection efficiency.

+++

For the $L$, we have 100 droplets $cm^{-3}$ with $\overline{r }= 10^{-3}$ cm, which gives a lwc of $4 \times 10^{-4} \ kg\,m^{-3}$, calculated in the cell below:

+++

#### Calculate the lwc

```{code-cell} ipython3
from a405.thermo.constants import constants as c
import numpy as np
N=100e6  # number/m^3
r = 1.e-5  #meters
L = c.rhol*4./3.*np.pi*N*r**3.
print(f"liquid water content lwc: {L:6.2g} kg/m^3")
```

#### Get the coefficient

+++

Taking
E=0.8 and $V=8000\;R$ cm/s, with R in cm, {eq}`eq:coal` becomes:

$$
\begin{align}
  \frac{ dR}{dt} &= \frac{ 4.2 \times  10^{-4} \times 8000 \times R \times  0.8}{4 \rho_l} \\
\frac{ dR}{R}& = 6.7 \times 10^{-4} dt
\end{align}
$$
where the numbers are calculated in the cell below:

```{code-cell} ipython3
lhs = 4.2e-4*8000*0.8/(4*1000)
print(f"{lhs=:.1e}")
```

#### Integrate both sides

$$
\int_{100 \times 10^{-6}\ m}^{1000 \times 10^{-6}\ m} \frac{dR}{R} dR & = 6.7 \times 10^{-4} \int_0^t dt = 5.04 \times 10^{-4} t \\
\frac{\ln 1000- \ln 100 }{6.7 \times 10^{-4}} &= t \\
t = 3437\ seconds &= 57.3\ minutes
$$

```{code-cell} ipython3
the_time = (np.log(1000) - np.log(100))/6.7e-4
print(f"growth time is {the_time:.1f} seconds")
the_minutes = the_time/60.
print(f"growth time is {the_minutes:.1f} minutes")
```

#### Problem 2 odeint solution

```{code-cell} ipython3
from scipy.integrate import odeint
import pandas as pd

def find_deriv(the_var,the_time):
    dr_dt = 6.7e-4*the_var
    return dr_dt

R0 = 100
the_time = np.linspace(0,3600)
    
sol = odeint(find_deriv,R0, the_time)
df_output = pd.DataFrame.from_records(sol,columns = ['R'])
df_output['time'] = the_time

fig, ax = plt.subplots(1,1)
ax.plot('time','R',data = df_output);
```

### Problem 3  -- Wallace and Hobbs problem 6.21  -- updraft

Derive an expression for the height h above cloud base of a droplet at time t that is growing by condensatyion only in a cloud with a steady updraft velocity w and supersaturation S. 

That is -- integrate

$$
\frac{dh}{dt} = w - \nu
$$
where h is the height, w is the updraft speed and $\nu$ is the terminal droplet fall speed
Use WH equation 6.24 for the terminal fall speed of a droplet together with WH eq 6.21 for the approximate droplet growth ratge

+++ {"jp-MarkdownHeadingCollapsed": true}

$$
\frac{ dh}{dt}  = w - v \label{one}
$$(one)

$$
  v = \frac{ 2 g \rho_l r^2}{9 \eta}
$$(two)

$$
r \frac{ dr}{dt}  = G_l S \\
\Rightarrow \frac{ 1}{2}\frac{dr^2 }{dt}   = G_l S 
$$(three)

Integrate {eq}`three` from $(t^\prime=0,r^\prime=0)$ to $(t^\prime=t,r^\prime = r)$

$$
r^2 = 2 G_l S t
$$(four)

Insert {eq}`four` into {eq}`one`

$$
\frac{ dh}{dt}  = w - \frac{ 2 g \rho_l 2 G_l S t}{9 \eta}
$$

+++ {"jp-MarkdownHeadingCollapsed": true}


$$
    \frac{ dh}{dt}  = w - \frac{ 2 g \rho_l 2 G_l S t}{9 \eta} \label{three} \\
\int_0^h  dh^\prime = \int_0^t \left [ w - \frac{ 2 g \rho_l 2 G_l S t^\prime}{9 \eta} \right ] dt^\prime \\
h= wt - \frac{ 2 g \rho_l  G_l S t^2}{9 \eta}
$$


+++ {"jp-MarkdownHeadingCollapsed": true}

### Problem 4 -- Wallace and Hobbs problem 6.24 -- falling precip

If a raindrop has a radius of 1 mm at cloud base, which is located 5 km above the ground, what will be its radius at the ground and how long will it take to reach the ground if the relative humidity between cloud base and ground is constant at 60%? [Hint: Use (6.21) and the relationship between v and r given in Exercise 6.23. If r is in micrometers, the value of Gl in (6.21) is 100 for cloud droplets, but for the large drop sizes considered in this problem the value of Gl should be taken as 700 to allow for ventilation effects.]
