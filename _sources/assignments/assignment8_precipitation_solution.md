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

Downlad link for [assignment8_precipitation_solution.ipynb](https://drive.google.com/file/d/1FStdQn59V8l3eYsyKWmPwqsPFsqgNg3R/view?usp=sharing)

## Problem 1 -- Cloud Condensation Nuclei Counter

Solution: {ref}`worksheet10_ccnc_solution`



## Precipitation problems

For problem 2:  First solve analytically then
then check that solution with python using odeint.  In your notebook compare your analytic equation with the numerical solution.  Do the same for problem 4, using the analytic equation from problem 3.



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
from matplotlib import pyplot as plt
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
Use WH equation 6.24 for the terminal fall speed of a droplet together with WH eq 6.21 for the approximate droplet growth rate

#### Problem 3 answer



$$
\frac{ dh}{dt}  = w - v \label{one}
$$(one)

$$
  v = \frac{ 2 g \rho_l r^2}{9 \eta}
$$(two)

write the supersaturation as $SS = (S - 1)$


$$
r \frac{ dr}{dt}  = G_l SS \\
\Rightarrow \frac{ 1}{2}\frac{dr^2 }{dt}   = G_l SS 
$$(three)

Integrate {eq}`three` from $(t^\prime=0,r^\prime=r_0)$ to $(t^\prime=t,r^\prime = r)$

$$
r^2 = r_0^2 + 2 G_l SS t
$$(four)

+++ {"jp-MarkdownHeadingCollapsed": true}

Write

$$
\beta =  \frac{ 2 g \rho_l}{9 \eta} \\
\nu = \beta r^2
$$(simple)

Insert {eq}`four` into {eq}`simple` and then use {eq}`one`

$$
\frac{ dh}{dt}  = w - \beta r^2 = w - \beta \left ( r_0^2 + 2 G_l SS t \right )\\
\int_{h0}^h dh^\prime = \int_0^t w dt^\prime - \int_0^t \beta \left ( r_0^2 + 2 G_l SS t \right ) dt^\prime \\
h = h_0 + w t - \beta r_0^2 t - \beta G_l SS t^2
$$(fall1)

+++ {"jp-MarkdownHeadingCollapsed": true}

So if the fall speed was larger than the updraft speed, we
could set $h=0$ and solve the quadratic equation for time $t$ when drop hits the ground.

We would also need a value for the dynamic viscosity $\eta$
for example from [this website](https://www.engineeringtoolbox.com/dynamic-absolute-kinematic-viscosity-d_412.html), which says that at 0 deg C, $\eta=19.8 \times 10^{-6}\ Pa\,s$

+++

#### Check values of $G_l$

Do a numerical check on $G_l$

- Thompkins 4.25

$$
r\frac{d r}{d t} \simeq \frac{D e_s^{\infty}}{\rho_L  R_v T}(S-1)
$$(eq:dropgrow)
where
- $D$ = $2.2 \times 10^{-5}$ $m^2\,s^{-1}$ at 0 deg C
- $R_v$ = 461.5 $J\,kg^{-1}\,K^{-1}$
- $\rho_l$ = 1000 $kg\,m^{-3}$
- $e_s^\infty$ = 611 $Pa$ at 0 deg C
- $T$ = 273.15 K

+++

(assign8_G_l)=
#### Numerical value for $G_l$ at 0 deg C

Confirm $G_l = 100\ m^2/s$ for small drops

```{code-cell} ipython3
D=2.2e-5  #m^2/s
es0 = 611 #Pa
Temp = 273.15 #K
G_l = D*es0/(c.rhol*c.Rv*Temp)
print(f"{G_l:.3g} m^s/s or {G_l*1.e12:.3g} μm^2/s at 0 deg C")
```

Which agrees with Wallace and Hobbs

+++

### Problem 4 Answer -- Wallace and Hobbs problem 6.24 -- falling precip

If a raindrop has a radius of 1 mm at cloud base, which is located 5 km above the ground, what will be its radius at the ground and how long will it take to reach the ground if the relative humidity between cloud base and ground is constant at 60%? [Hint: Use (6.21) and the relationship between v and r given in Exercise 6.23. If r is in micrometers, the value of Gl in (6.21) is 100 for cloud droplets, but for the large drop sizes considered in this problem the value of Gl should be taken as 700 to allow for ventilation effects.]

+++

#### Numerical solution

For extra credit:  

1) Make this work with v = -6000r
2) change fall speed to that used in analytic solution and compare
3) Repeat this using solve_ivp, stopping at h=0

```{code-cell} ipython3
from scipy.integrate import solve_ivp
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


def hit_ground(t, the_vars):
    """Event: when height (the_vars[0]) is zero."""
    return the_vars[0]

hit_ground.terminal = True
hit_ground.direction = -1

G_l = 700e-12
S = 0.6

def find_derivs(the_time,the_vars):
    """
    the_vars[0] = h  (height, m)
    the_vars[1] = r (radius, m)
    """
    h, r = the_vars
    dh_dt = -6000*r #WH prob 6.23, m
    dr_dt =  G_l*(S - 1)/r #WH eq. 6.21
    #print(f"{r=}, {h=}, {the_time=}")
    return [dh_dt,dr_dt]

init_0 = [5000,0.001]  #h, r
tspan=[0,5000]

sol = solve_ivp(
    find_derivs,
    tspan,
    init_0,
    events=[hit_ground],
    dense_output=True # Enable interpolation for plotting smooth curve
)
```

```{code-cell} ipython3
the_times = np.arange(0,sol.t[-1],1)
yplot = sol.sol(the_times)
df_output = pd.DataFrame.from_records(yplot.T,columns = ['height','radius'])
df_output['time']=the_times
df_output
```

```{code-cell} ipython3
fig, ax = plt.subplots(1,1)
ax.grid(True)
ax.plot(df_output['radius']*1.e6,df_output['height'])
ax.set_xlabel("drop radius (microns)")
ax.set_ylabel("height (m)")
```

### Check this answer analytically

$$ 
r^2 &=r_0^2+2 G_l SS t \\
\frac{d h}{d t}&=-6000 r \\ 
d h &=-6000 \sqrt {\left(r_0^2+2 G_l SS t\right) }\, dt
$$(uint)

+++

#### Integrate {eq}`uint` with change of variables

$$
u & =r_0^2+2 G_l SS t \\
d u & =2 G_l SS d t \\
h-h_0 & =-\int_0^t \frac{6000}{2 G_l SS} u^{1 / 2} d u \\
h-h_0 & = -\left.\int_0^t \frac{2}{3} \frac{6000}{2 G_l SS} u^{3 / 2}\right|_0 ^t \\
h-h_0 & = - \left. \frac{2000}{ G_l SS} u^{3 / 2}\right|_0 ^t \\
h-h_0 & = - \frac{2000}{ G_l SS}\left[\left(r_0^2+2 G_l SS t\right)^{3 / 2} - r_0^3\right]
$$

+++

#### Put in some numbers

```{code-cell} ipython3
Gl = 7.e-10
r0 = 0.001
SS = -0.4
time = 985.
r0=0.001
h0=5000.
```

```{code-cell} ipython3
u=r0**2. + 2*Gl*SS*time
h = h0 -2000/(Gl*SS)*(u**1.5 - r0**3.)
print(f"final analytic height at {time=} seconds in {h=:.1f} m")
```

```{raw-cell}
So the analytic and numerical solutions match pretty closely
```

### Repeat numerical solution using fall speed from {eq}`two`

```{code-cell} ipython3
eta = 19.8e-6
beta = 2*c.g0*c.rhol/(9.*eta)

def find_derivs(the_time,the_vars):
    """
    the_vars[0] = h  (height, m)
    the_vars[1] = r (radius, m)
    """
    h, r = the_vars
    dh_dt = -beta*r**2.  #WH eq. 6.24
    dr_dt =  G_l*(S - 1)/r #WH eq. 6.21
    #print(f"{r=}, {h=}, {the_time=}")
    return [dh_dt,dr_dt]
```

```{code-cell} ipython3
init_0 = [5000,0.001]  #h, r
tspan=[0,5000]

sol = solve_ivp(
    find_derivs,
    tspan,
    init_0,
    events=[hit_ground],
    dense_output=True # Enable interpolation for plotting smooth curve
)
```

```{code-cell} ipython3
the_times = np.arange(0,sol.t[-1],1)
yplot = sol.sol(the_times)
df_output = pd.DataFrame.from_records(yplot.T,columns = ['height','radius'])
df_output['time']=the_times
df_output
```

Compare with the analytic result from {eq}`fall1`

```{code-cell} ipython3
time = 46
r0 = 0.001
h0=5000
Gl = 7.e-10
SS = -0.4
h = h0 - beta*r0**2.*time - beta*Gl*SS*time**2
print(f"final analytic height at {time=} seconds in {h=:.1f} m")
```

Again, a close match.  Note how much difference the slower fall that includes drag makes on the answer.

```{code-cell} ipython3

```
