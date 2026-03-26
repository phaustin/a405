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

(assign7_solution)=
# Assignment 7 solution

+++

## Question 1

+++

1. Assuming that cloud condensation nuclei (CCN) are removed from the atmosphere by first serving as the centers on which cloud droplets form, and the droplets subsequently grow to form precipitation particles, estimate the residence time of a CCN in a column extending from the surface of the Earth to an altitude of 5 km. Assume that the annual rainfall is 100 cm/year and the cloud liquid water content is 0.30 $g/kg$ .  *Hint:  Assume that all drops in the cloud droplets have  a radii of 10 microns and that every droplet contains exactly 1 CCN.   How many CCN are in 1 kg of air?  About how many kg of air are there in a 5 km column?  About how many CCN are taken out by a rain rate of 1 m/year?  Find the time constant for removal of the form  1/N dN/dt = 1/tau*

+++

### Question 1 answer

+++

First find $n$, the number of aersols/volume 
by assuming that the mean droplet size is 
$\overline{r}=10\ \mu m$  and
that there is one aerosol in each droplet.

$$
  w_l = 0.3\ \ g\,m^{-3} = \frac{4}{3}\pi \rho_l N\overline{ r}^3  
$$
which gives $n=72\ \ cm^{-3} = 72\ \times 10^6\ \ {kg^{-1}}$ if $\rho_{air} \approx 1\ \ kg\,m^{-3} $

+++

If the column is well mixed, then $n=constant$ between the surface
and 5 km $\approx$ 500 hPa.   Integrate
the hydrostatic equation between those levels to find the mass of dry air:

$$
M=  \int_0^{5000} \rho dz = \frac{1 }{g} \int_{50000}^{100000} dp \nonumber \\
= \frac{ 1}{10} 5 \times 10^4 = 5000\ kg\,m^{-2}
$$
So to get the total CCN in a column multiply $M \times n$ =
$N=3.6 \times 10^{11}$ CCN in a $1 \ m^{2}$ column.

+++

Now what about $dN/dt$?  If the rainfall is 1 m/year
= 1000 $kg\,m^{-2}\,yr^{-1}$ and $\overline{ r}=10\ \mu m$
then the mass of an average drop is 

$$
\frac{ 4}{3}\pi \rho_l N\overline{ r}^3
$$
and the drops removed in a year is:

$$
  \frac{1000 \ kg\,m^{-2}\,yr^{-1} }{\frac{ 4}{3}\pi \rho_l N (10^{-5})^3} = 2.39 \times 10^{14}\ drops/year
$$

So put these together:

$$
  \tau = \frac{N}{dN/dt} = \frac{ 3.6 \times 10 11}{ 2.39 \times 10^{14}} 
= 0.0015\ years = 0.55\ days
$$

+++

## Question 2

+++

A drop with an initial radius of 100 µm falls through a cloud containing 100 droplets per cubic centimeter that it collects in a continuous manner with a collection efficiency of 0.800. If all the cloud droplets have a radius of 10 µm, how long will it take for the drop to reach a radius of 1 mm? You may assume that for the drops of the size considered in this problem the terminal fall speed v (in $m s^{-1}$) of a drop of radius r (in meters) is given by $v= 8 x 10^3\;r$. Assume that the cloud droplets are stationary and that the updraft velocity in the cloud is negligible.  Hint:  Integrate Thompkins equation 4.30 (not 4.28) analytically -- you  can also check numerically with python

+++

### Question 2 answer

+++

Equation 4.30 (be sure you can derive this):

$$
  \frac{ dR}{dt}  = \frac{L V \epsilon E}{4\rho_l}
$$
where $R$ is the radius of the collector drop in meters, $L$ is the liquid water content in $kg\,m^{-3}$, $\rho_l$ is the density of liquid water in $kg\,m^{-3}$, $V=6000R$ is the differential fall speed in m/s, $\epsilon=1$ is the coalescence efficiency and $E0.8$ is the collection efficiency.

+++

 
For the $L$, we have 100 droplets $cm^{-3}$ with $\overline{r }= 10^{-3}$ cm, which gives a lwc of $4 \times 10^{-4} \ kg\,m^{-3}$.

```{code-cell} ipython3
from a405.thermo.constants import constants as c
import numpy as np
N=100e6  # number/m^3
r = 1.e-5  #meters
L = c.rhol*4./3.*np.pi*N*r**3.
print(f"{L:6.2g} kg/m^3")
```

Taking
E=0.8 and $V=6000 r$ cm/s, with r in cm:

$$
\begin{align}
  \frac{ dR}{dt} &= \frac{ 6000 \times R \times 4.2 \times  10^{-4} \times 0.8}{4 \rho_l} \\
\frac{ dR}{R}& = 1.6 \pi 10^{-4} dt
\end{align}
$$

+++

Integrating both sides:

$$
\begin{align}
\int_{100 \times 10^{-6}\ m}^{1000 \times 10^{-6}\ m} \frac{dR}{R} dR & = 5.04 \times 10^{-4} \int_0^t dt = 5.04 \times 10^{-4} t \\
\frac{\ln 1000- \ln 100 }{5.04 \times 10^{-4}} &= t \\
t = 4569\ seconds &= 76.15\ minutes
\end{align}
$$

```{code-cell} ipython3
6000*4.2e-4*0.8/(4.*c.rhol)
```

```{code-cell} ipython3
(np.log(1000) - np.log(100))/5.04e-4
```

```{code-cell} ipython3
4569/60.
```

## Question 3

+++

3. Compare the droplet growth equation in Thompkins equation 4.24 with  with Lohmann 7.28 for 3 micron drop nucleated on a $1 \times 10^{-18}$ kg ammonium sulphate aerosol.  Show numerical values for all the terms, and the total percentage difference in the $dr/dt$ for the two equations.

+++

### Question 3 answer

Reuse some code from worksheets [koehler1](https://phaustin.github.io/a405_2024/notebooks/worksheets/kohler_equilibrium_students.html) and  [koehler2](https://phaustin.github.io/a405_2024/notebooks/koehler2_worksheet.html) to evaluate:

Thompkins 4.24

$$
\frac{d r}{d t}=\frac{D e_s^{\infty}}{\rho_L r R_v T}\left(S-1-\frac{a}{r T}+\frac{b}{r^3}\right)
$$

Lohmann 7.28:

$$
r \frac{d r}{d t}=\frac{(S-1)-a / r+b / r^3}{F_k+F_d} .
$$

$$
F_k^l=\left(\frac{L_v}{R_v T}-1\right) \frac{L_v}{K T} \approx\left(\frac{L_v^2}{K R_v T^2}\right) ;
$$

$$
F_d^l=\frac{1}{D_v \rho_{v s}}=\frac{R_v T}{D_v e_{s, w}(T)}
$$


Note that Thompkins takes the temperature out of his $a$ coefficient, but it's identical to Lohmann if you move temperature into $a$.

```{code-cell} ipython3
from a405.thermo.thermlib import find_esat, find_lv
from a405.utils.helper_funs import make_tuple
from a405.thermo.constants import constants as c
import numpy as np

temp = 280 #K
press = 80000 #Pa
D = 2.8e-5 #m^2/s
K = 0.025 #W/m/K
r = 3.e-6 #m
#
# from the koehler worksheets
#
aero_specs = {
    "Ms": 132,
    "Mw": 18.0,
    "Sigma": 0.075,
    "vanHoff": 3.0,
    "rho": 1775,
    "mass": 1e-18,
    "comments": "ammonum sulfate (NH4)2SO4"
}

aero = make_tuple(aero_specs)
lv = find_lv(temp)

a=(2.*aero.Sigma)/(c.Rv*temp*c.rhol)  #curvature term
b=(aero.vanHoff*aero.Mw)/((4./3.)*np.pi*c.rhol*aero.Ms)*aero.mass
esat = find_esat(temp)

a_over_r = a/r
b_over_r3 = b/r**3.
thompkins_coefficient = D*esat/(c.rhol*c.Rv*temp)
print(f"{(a_over_r, b_over_r3,thompkins_coefficient)=}")
lohmann_FK = lv**2/(K*c.Rv*temp**2.)
lohmann_FD = 1/thompkins_coefficient
lohmann_coefficient = 1/(lohmann_FK + lohmann_FD)
print(f"{lohmann_coefficient/thompkins_coefficient=:.4f}")
```

## Question 4


Derive Lohmann equations 7.31 and 7.32:

$$
\frac{d S}{d t}=Q_1 \frac{d z}{d t}-Q_2 \frac{d r_l}{d t}
$$

$$
\begin{aligned} 
& Q_1=\frac{1}{T}\left(\frac{\epsilon L_v g}{R_d c_p T}-\frac{g}{R_d}\right), \\
& Q_2=\rho\left(\frac{R_d T}{\epsilon e_{s}}+\frac{\epsilon L_v^2}{p T c_p}\right)
\end{aligned}
$$

+++

   Problem 4 hints:
   
   - Start with this approximate relation between the vapor mixing ratio and the saturations

      $$
        r_v = \frac{\rho_v}{\rho_d} = \frac{R_d}{R_v} \frac{e}{p - e} \approx \frac{R_d}{R_v} \frac{e}{p}=
       \epsilon \frac{e}{p} = \frac{S \epsilon e_s}{p} 
      $$
             
   - first show that with the chain rule wind up with this:
   
        $$
        \frac{d r_v}{dt} = S \left [ \frac{-\epsilon e_s}{p^2} 
        \left ( \frac{-g p V}{R_d T} \right ) + \frac{\epsilon}{p} \left ( 
        \frac{\epsilon e_s l_v}{R_d T^2} \right ) \frac{dT}{dt} \right ]
        + \frac{\epsilon e_s}{p} \frac{dS}{dt}
        $$(eq:rv)

     where $V$ is the vertical velocity $dz/dt$ and I've used the Clausius-Clapeyron equation and
     assumed hydrostatic balance:
     
     $$
     \begin{align}
    \frac{de_s}{dT} &= \frac{\epsilon l_v e_s}{R_d T^2}\\
    \frac{dp}{dt} & =  - \rho g \frac{dz}{dt} = -\frac{g p}{R_d T} V
     \end{align}
     $$

   - To get Lohmann's coefficients, recognize that $dr_v/dt = - dr_l/dt$ if total water is conserved

+++

### Question 4 answer

+++

We need to eliminate dT/dt from {eq}`eq:rv`

For $dT/dt$ use the fact that $h_m = c_p T + l_v r_v + gz=constant$ is conserved for adiatbatic ascent take the derivative to get

$$
\frac{d T}{d t}=\frac{l_v}{c_p} \frac{d r_v}{d t}-\frac{g V}{c_p}
$$(eq:hm)

+++

Plug {eq}`eq:hm` into {eq}`eq:rv` and rearrange

$$
\frac{d S}{d t}=-\frac{g V}{R_d T}-\frac{l_v}{R_v T^2}\left[-\frac{l_v}{c_p} \frac{d r_v}{d t}-\frac{g V}{c_p}\right]+\frac{p}{\epsilon e_s} \frac{d r_v}{d t}
$$
where I've made the approximation that $S \approx 1$

+++

$$
\frac{d S}{d t}=\frac{g}{T}\left[\frac{l_v}{c_p R_v T}-\frac{1}{R_d}\right] V-\left[\frac{p}{e_s \epsilon}+\frac{l_v^2}{R_v c_p T^2}\right] \frac{dr_l}{dt}
$$

+++

Finally, eliminate $R_v$ using the defination $\epsilon = R_d/R_v$ to get Lohmann's $Q_1$ and $Q_2$

```{code-cell} ipython3

```
