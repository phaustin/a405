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

(assign7_aerosols_solution)=
# Assignment 7: aerosols: solution

- Download link: [assign7_aerosols_solution](https://drive.google.com/file/d/1a2NHWaATsAaY-U5m78yWG9kCJySwBI2w/view?usp=sharing)

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
 $$(eq:scrit)

where $m_{aer}$ is the ammonium sulphate aerosol mass in kg.

Note that this is why a cloud chamber can get the aerosol mass distribution from a series of
saturation and light scattering measurements as smaller and smaller aerosols are pushed over
their critical supersaturation and activated.

+++

### Problem 1 answer

+++

Given a and b:

$$
a=\frac{2 \sigma}{\rho_l R_v T}
$$
$$
b=\frac{i m M_w}{4 / 3 \pi \rho_s M_s}
$$

+++

Insert the definition of b into {eq}`eq:scrit`:

$$
(S^* - 1 ) = \left (\frac{4a^3 \,4/3 \pi \rho_s M_s }{27 i M_w} \right )^{1/2}  m_{aer}^{-1/2}
$$

```{code-cell} ipython3
from a405.thermo.constants import constants as c
import numpy as np
sigma = 0.075 #N/m
temp = 280 # K
rhos = 1775 #kg/m^3
i = 3
Ms = 132  #kg/kmol
Mw = 18 #kg/kmol

a = 2*sigma/(c.rhol*c.Rv*temp)
numerator = 4*a**3.*4./3.*np.pi*rhos*Ms
denom = 27*i*Mw
coeff = (numerator/denom)**0.5
coeff
```

```{code-cell} ipython3
import numpy as np

Ms=132.  #kg/kmole
Mw=18.  #kg/kmole
Rhol=1.e3  #1000
Sigma=0.075  #N/m
rhoaero=1775.  #kg/m^3
Rv=461.   #J/kg/K
vanhoff=3
Tinit=280.  #K
a=(2.*Sigma)/(Rv*Tinit*Rhol) # units m
bnomass=(vanhoff*Mw)/((4./3.)*np.pi*Rhol*Ms)  #units m^3/kg
Scoeff=((4*a**3)/(27*bnomass))**0.5  #units  kg^{-0.5}
print(f'Scrit coeff without aerosol mass: {Scoeff:6.3e} kg**0.5')
```

#### some trial values

```{code-cell} ipython3
m_aer = 1.e-19
SS = Scoeff*m_aer**(-0.5)
print(f"{m_aer=}, {SS*100=:.2f} %")
m_aer = 1.e-17
SS = Scoeff*m_aer**(-0.5)
print(f"{m_aer=}, {SS*100=:.2f} %")
```

## Problem 2: Koehler stability

+++

2. Show that the expression for second derivative of the thermodynamic potential derived in the  Koehler stability notes:

     $$
     \frac{\delta ^2G}{\delta r^2} = - 4 \pi R_v T \rho_l \left [ 2 a - r^2 \left ( 1 +
           \frac{b}{r^3} \right ) \frac{3b}{r^4}  \right ] + 8 \pi \sigma
     $$(eq:gfull)

     Changes sign from stable (positive) to unstable (negative) at $r_{crit}$.

     Hint -- first show that the second derivative is zero at the critical radius.  Then show that
     the third derivative is negative above and below the critical radius, which means that
     there has to be a sign change from + to -.

+++

### Problem 2 answer

+++

First confirm that 
$\frac{\delta^2 G}{\delta r^2}=0$ 
at the critical radius:

Inserting 
$r_{crit}=\left(\frac{3 b}{a}\right)^{1 / 2}$
into {eq}`eq:gfull`:

$$
\begin{aligned}
\frac{\delta^2 G}{\partial r^2}= & -4 \pi R_v T \rho_l\left(2 a-\frac{3 b a}{3 b}\right)+8 \pi \sigma \\
& =-4 \pi R_v T \rho_l a+8 \pi \sigma
\end{aligned}
$$

but 

$$a=\frac{2\sigma}{R_v T \rho_l}$$

so

$$
\frac{\delta^2 G}{\delta r^2}=-4 \pi R_v T \rho_l \frac{2 \sigma}{R_v T \rho_l}+8 \pi \sigma=-8 \pi \sigma+8 \pi \sigma=0
$$

+++

We walt to confirm that 

$\frac{\delta^2 G}{\delta r^2}<0$ (unstable) for 
$r>r_{crit}$ and 
$\frac{\delta^2 G}{\delta r^2}>0$ (stable)
for $r<r_{crit}$.

That is equivalent to showing that 

$$
\frac{\delta}{\delta r}\left(\frac{\delta^2 G}{\delta r^2}\right)<0 \text{ at } r=r_{crit}
$$

+++

Sure enough:

$$
\begin{aligned}
& \frac{\delta}{\delta r}\left(\frac{\delta^2 G}{\delta r^2}\right)=\frac{\delta}{\delta r}\left(-4 \pi R_V T \rho_l\left(-\frac{3 b}{r^2}\right)\right) \\
& =-4 \pi R_V T \rho_l \frac{3 b}{r^3}<0 \text { at } r=r_{crit}.
\end{aligned}
$$


## Problem 3: CCN lifetime

+++

3. Assuming that cloud condensation nuclei (CCN) are removed from the atmosphere by first serving as the centers on which cloud droplets form, and the droplets subsequently grow to form precipitation particles, estimate the residence time of a CCN in a column extending from the surface of the Earth to an altitude of 5 km. Assume that the annual rainfall is 100 cm/year and the cloud liquid water content is 0.30 $g/kg$ .  *Hint:  Assume that all drops in the cloud droplets have  a radii of 10 microns and that every droplet contains exactly 1 CCN.   How many CCN are in 1 kg of air?  About how many kg of air are there in a 5 km column?  About how many CCN are taken out by a rain rate of 1 m/year?  Find the time constant for removal of the form  1/N dN/dt = 1/tau*

### Problem 3 answer

+++

First find $n$, the number of aersols/volume 
by assuming that the mean droplet size is 
$\overline{r}=10\ \mu m$  and
that there is one aerosol in each droplet.

$$
  w_l = 0.3\ \ g\,m^{-3} = \frac{4}{3}\pi \rho_l N\overline{ r}^3  
$$

which gives 
$n=72\ \ cm^{-3} = 72\ \times 10^6\ \ {kg^{-1}}$ if $\rho_{air} \approx 1\ \ kg\,m^{-3}$

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
  \tau = \frac{N}{dN/dt} = \frac{ 3.6 \times 10^{11}}{ 2.39 \times 10^{14}} 
= 0.0015\ years = 0.55\ days
$$
