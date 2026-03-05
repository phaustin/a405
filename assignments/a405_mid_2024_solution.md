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

(mid_2024_solutions)=
# 2024 midterm solutions

+++

## Question 1: cooling

For air at 700 hPa with 7 g/kg of vapor (saturated) and 1 g/kg of liquid.

-  (4 points)  Find: 
  - The LCL of this air: **755 hPa**
  - The approximate temperature if it was brought adiabatically to a pressure of 1000 hPa.
    - **$\theta$ at LCL = 303 K**
 
- (8 points)  Suppose this air was cooled by 6 degrees C at a constant pressure of 700 hPa.  Find:

  - The amount of liquid water condensed by the cooling (g/kg)
    - **$r_{vafter}$ = 4.5 g/kg so now 3.5 g/kg of liquid, an extra 2.5 g/kg has been condensed**
  - The new LCL, assuming no precipitation: **914 hPa**
  - The amount of energy $\Delta q_{out}$ (J/kg) shed to the environment during the cooling.

    $$
    \begin{align}
    \Delta h &= c_p \Delta T + l_v \Delta r_v\\
             &= -(1004 \times 6) + -(2.5 \times 10^6 \times 2.5 \times 10^{-3}) \\
             &= -6024 + -6250  = -12274\ J/kg
    \end{align}
    $$

+++

```{figure} images/question1_answer.png
---
width: 70%
name: figure1
---
Question 1 cooling problem
```

+++

## Question 2: Taylor series problem (9 points)

- (3 points)  Use a Taylor series to expand the following expression with respect to a temperature perturbation (this should yield 3 terms)
  
$$
d  \left ( \frac{l_v(T) r_s(T)}{c_p T} \right ) 
$$

- (6 points) Use the attached tephigram to estimate the size of each of the three terms for a 2 deg C  temperature perturbation at a pressure of 900 hPa and a temperature of 10 deg C.  To first order, what is error in percent caused by dropping two of those terms in the following approximation?  (Hint: get $dr_s/dT$ from the tephigram, not the C-C equation)

$$
d  \left ( \frac{l_v(T) r_s(T)}{c_p T} \right ) \approx  \left ( \frac{l_v(T) }{c_p T} \right ) dr_s
$$

+++

### Question 2 answer

$$
\begin{aligned}
& f=\frac{l_v r_s}{c_p T} \\
& \frac{d f}{d T}=\frac{l_v r_s}{c_p} \\
& \frac{d f}{d T}=\frac{l_v}{c_p T} \frac{d r_s}{d T}+\frac{r_s}{c_p T} \frac{d l_v}{d T}-\frac{l_v r_s}{c_p T^2}
\end{aligned}
$$

+++

- Put in some numbers:  at 10 deg C and 900 hPa the tephigram says

  - T = 283.15 K
  - $r_s$ = 8.6 g/kg
  - $l_v$ = 2.5e6 J/kg
  - $c_p$ = 1004 J/kg/K
  - $dr_s/dT$ \approx (9.2 - 8) $g\,kg^{-1}$/2 K = 1.2e-3/2 = 0.6e-3 kg/kg/K
  - $dl_v/dT$ = $c_{pv} - c_l$ = 1870 - 4190 = -2320  J/kg/K

+++

Putting in numbers shows the second and third terms contribute about 6% to the total

```{code-cell} ipython3
from a405.thermo.thermlib import find_rsat, find_theta
lv=2.5e6
cp = 1004 
T = 283.15
rs = 8.6e-3
drs_dT = 0.6e-3
dlv_dT = -2320
term1 = lv*drs_dT/cp/T
term2 = rs*dlv_dT/cp/T
term3 = -lv*rs/cp/T**2.
allterms = (term1 + term2 + term3)
term1/allterms*100.,term2/allterms*100,term3/allterms*100
```

```{code-cell} ipython3
find_rsat(283.15,900.e2), find_rsat(284.15,900.e2), find_rsat(282.15,900.e2)
```

## Question 3: Mixing (12 points):

  Surface air at 1000 hPa with a temperature of 20 deg C and a dewpoint of 16 deg C is lifted adiabatically to 800 hPa, where it mixes 50/50 with air that has a $\theta_e$ of 307 K and a mixing ratio of 4 g/kg.  Use the tephigram to find:

  - The $\theta_e$ and LCL of the of the surface air
    - 324.4 K
    - $r_t$ = 7.8 g/kg
    - 942 hPa
  - The $\theta_e$, $r_v$ and $r_l$ of the mixture
    - $\theta_e$ = 315.7 K
    - $r_v$ = 7 g/kg
    - $r_l$ = 0.8 g/kg
  - The LCL of the mixture
    - 850 hPa
  - The temperature of the mixture at 800 hPa
    - 5.3 deg C

+++

```{figure} images/question3_answer.png
---
width: 70%
name: figure1
---
Question 3 mixing problem
```

+++

## Question 1 cooling code

```{code-cell} ipython3
from a405.soundings.wyominglib import read_soundings
from a405.skewT.skewlib import makeSkewDry
from a405.thermo.thermlib import convertTempToSkew, find_lcl, find_thetaet, find_Td, tinvert_thetae
from a405.thermo.thermlib import find_rsat, find_theta
from a405.thermo.constants import constants as c
from a405.skewT.fullskew import makeSkewWet,find_corners,make_default_labels

import datetime
import pytz
import numpy as np
from matplotlib import pyplot as plt
```

```{code-cell} ipython3
def label_fun():
    """
    override the default rs labels with a tighter mesh
    """
    from numpy import arange
    #
    # get the default labels
    #
    tempLabels,rsLabels, thetaLabels, thetaeLabels = make_default_labels()
    #
    # change the temperature and rs grids
    #
    tempLabels = range(-40, 50, 2)
    rsLabels = [0.1, 0.25, 0.5, 1, 2, 3] + list(np.arange(4, 28, 1)) 
    thetaeLabels = np.arange(295,335,3)
    return tempLabels,rsLabels, thetaLabels, thetaeLabels

pa2hPa = 1.e-2
```

```{code-cell} ipython3
fig,ax =plt.subplots(1,1,figsize=(11,11))
skew = 35
corners=[5,20]
ax, skew = makeSkewWet(ax, corners=corners, skew=skew,label_fun=label_fun)
ax.set_title('Cooling problem')

xcorners=find_corners(corners,skew=skew)
ax.set(xlim=xcorners,ylim=[1000,600])
fig.savefig('cooling.pdf')
```

### Find temperature at 700hPa (blue diamond)

Saturated parcel, so plot either temperature or dewpoint

```{code-cell} ipython3
press_700 = 700e2
rv=7e-3
rl = 1.e-3
rt_cloud = rv + rl
#find_thetaet(Td, rt, T, p)
Td_700 = find_Td(rv,press_700)
temp_700 = Td_700
thetae_before = find_thetaet(Td_700, rt_cloud, temp_700,press_700)
print(f"{(temp_700 - c.Tc)=:0.1f} deg C")
#
# plot the dewpoint temperature
#
xplot=convertTempToSkew(Td_700 - c.Tc,press_700*pa2hPa,skew)
bot=ax.plot(xplot, press_700*pa2hPa, 'bd', markersize=14, markerfacecolor='r',label="before cooling")
```

### Find the lcl prior cooling

Plot as a green diamond

```{code-cell} ipython3
press_900 = 900.e2
Temp_900,rv_900,rl_900=tinvert_thetae(thetae_before,rt_cloud,press_900)
# find lcl
Td_900 = find_Td(rv_900, press_900)
T_lcl, press_lcl = find_lcl(Td_900,Temp_900,press_900)
xplot = xplot=convertTempToSkew(T_lcl - c.Tc,press_lcl*pa2hPa,skew)
lcl=ax.plot(xplot, press_lcl*pa2hPa, 'bd', markersize=14, markerfacecolor='g',
           label = "original lcl")
theta_lcl = find_theta(T_lcl,press_lcl)
print(f"LCL potential temperature {theta_lcl:0.1f} K")
print(f"LCL pressure {press_lcl*pa2hPa:0.1f} hPa")
ax.legend()
display(fig)
```

### Check that $\theta_e$ hasn't changed

```{code-cell} ipython3
thetae_check= find_thetaet(Td_900, rt_cloud, Temp_900,press_900)
thetae_check, thetae_before
```

### now cool by 6 K

```{code-cell} ipython3
new_temp = temp_700 - 6
thetae_after = find_thetaet(Td_700, rt_cloud, new_temp,press_700)
Temp_after,rv_after,rl_after=tinvert_thetae(thetae_after,rt_cloud,press_700)
print(f"{(thetae_before - thetae_after)=} K")
print(f"New rv {rv_after*1.e3:0.1f} g/kg")
xplot = xplot=convertTempToSkew(Temp_after - c.Tc,press_700*pa2hPa,skew)
lcl=ax.plot(xplot, press_700*pa2hPa, 'bd', markersize=14, markerfacecolor='b',
           label = "after 6K cooling")
ax.legend()
display(fig)
```

### Find the new lcl

Go down to the surface in case lcl is below 900

```{code-cell} ipython3
rv_after
```

```{code-cell} ipython3
press_1000=1.e5 # Pa
Temp_1000,rv_1000,rl_1000=tinvert_thetae(thetae_after,rt_cloud,press_1000)
Td_1000 = find_Td(rt_cloud, press_1000)
T_lcl, press_lcl = find_lcl(Td_1000,Temp_1000,press_1000)
print(f"{(T_lcl, press_lcl)=}")
xplot = xplot=convertTempToSkew(T_lcl - c.Tc,press_lcl*pa2hPa,skew)
lcl=ax.plot(xplot, press_lcl*pa2hPa, 'cd', markersize=14, markerfacecolor='c',
           label = "new lcl")
ax.legend()
fig.savefig("images/question1_answer.png")
display(fig)
```

## Question 3 mixing code (12 points):

  Surface air at 1000 hPa with a temperature of 20 deg C and a dewpoint of 16 deg C is lifted adiabatically to 800 hPa, where it mixes 50/50 with air that has a $\theta_e$ of 307 K and a mixing ratio of 4 g/kg.  Use the tephigram to find:

  - The $\theta_e$ and LCL of the of the surface air
  - The $\theta_e$, $r_v$ and $r_l$ of the mixture
  - The LCL of the mixture
  - The temperature of the mixture at 800 hPa

+++

### get $\theta_e$ of the surface air

```{code-cell} ipython3
temp_1000 = 20 + c.Tc
press_1000 = 1.e5 
Td_1000 = 16 + c.Tc
rt_cloud = find_rsat(Td_1000,press_1000)
thetae_cloud = find_thetaet(Td_1000, rt_cloud, temp_1000,press_1000)
T_lcl, press_lcl = find_lcl(Td_1000,temp_1000,press_1000)
print(f"thetae: {thetae_cloud:0.1f} K, lcl_press {press_lcl:0.1f} Pa")
```

### Plot T, Td and the lcl

```{code-cell} ipython3
fig,ax =plt.subplots(1,1,figsize=(11,11))
skew = 35
corners=[5,21]
ax, skew = makeSkewWet(ax, corners=corners, skew=skew,label_fun=label_fun)
ax.set_title('mixing problem')
xcorners=find_corners(corners,skew=skew)
ax.set(xlim=xcorners,ylim=[1005,700])
xplot=convertTempToSkew(temp_1000 - c.Tc,press_1000*pa2hPa,skew)
ax.plot(xplot, press_1000*pa2hPa, 'gd', markersize=14, markerfacecolor='g',
           label = "temp_sfc")
xplot=convertTempToSkew(Td_1000 - c.Tc,press_1000*pa2hPa,skew)
ax.plot(xplot, press_1000*pa2hPa,'bd', markersize=14, markerfacecolor='b',
           label = "dewpoint_sfc")
xplot=convertTempToSkew(T_lcl - c.Tc,press_lcl*pa2hPa,skew)
ax.plot(xplot, press_lcl*pa2hPa,'kd', markersize=14, markerfacecolor='k',
           label = "lcl_cloud");
```

### Lift to 800

```{code-cell} ipython3
press_800 = 8.e4 #Pa
Temp_800,rv_800,rl_800=tinvert_thetae(thetae_cloud,rt_cloud,press_800)
xplot=convertTempToSkew(Temp_800 - c.Tc,press_800*pa2hPa,skew)
ax.plot(xplot, press_800*pa2hPa,'rd', markersize=14, markerfacecolor='r',
           label = "cloud_800");
```

```{code-cell} ipython3
display(fig)
```

### add the environment

```{code-cell} ipython3
rt_env = 4.e-3
thetae_env = 307.
Td_env =  find_Td(rt_env, press_800)
Temp_env,rv_env,rl_env=tinvert_thetae(thetae_env,rt_env,press_800)
T_lclenv, press_lclenv = find_lcl(Td_env,Temp_env,press_800)
```

```{code-cell} ipython3
xplot=convertTempToSkew(Temp_env - c.Tc,press_800*pa2hPa,skew)
ax.plot(xplot, press_800*pa2hPa,'gs', markersize=14, markerfacecolor='g',
           label = "temp_env")
xplot=convertTempToSkew(Td_env - c.Tc,press_800*pa2hPa,skew)
ax.plot(xplot, press_800*pa2hPa,'bs', markersize=14, markerfacecolor='b',
           label = "Td_env")
xplot=convertTempToSkew(T_lclenv - c.Tc,press_lclenv*pa2hPa,skew)
ax.plot(xplot, press_lclenv*pa2hPa,'ks', markersize=14, markerfacecolor='k',
           label = "lcl_env")
ax.legend()
display(fig)
```

### add the mixture

```{code-cell} ipython3
thetae_env
```

```{code-cell} ipython3
thetae_mix = 0.5*thetae_env + 0.5*thetae_cloud
rt_mix = 0.5*rt_env + 0.5*rt_cloud
Temp_mix,rv_mix,rl_mix=tinvert_thetae(thetae_mix,rt_mix,press_800)
Temp_1000,rv_1000,rl_1000=tinvert_thetae(thetae_mix,rt_mix,press_1000)
Td_1000 = find_Td(rv_1000,press_1000)
T_lclmix, press_lclmix = find_lcl(Td_1000,Temp_1000,press_1000)
print(f"{thetae_mix=:0.1f} K, {rt_mix=:0.4f} kg/kg")
print(f"{Temp_mix - c.Tc=:0.1f} K, {rv_mix=:0.4f} kg/kg, {rl_mix=:0.4f} kg/kg")
print(f"{press_lclmix=:0.1f} Pa")
```

```{code-cell} ipython3
xplot=convertTempToSkew(Temp_mix - c.Tc,press_800*pa2hPa,skew)
ax.plot(xplot, press_800*pa2hPa,'co', markersize=14, markerfacecolor='c',
           label = "mixture")
xplot=convertTempToSkew(T_lclmix - c.Tc,press_lclmix*pa2hPa,skew)
ax.plot(xplot, press_lclmix*pa2hPa,'cs', markersize=14, markerfacecolor='c',
           label = "mixture lcl")
ax.legend()
fig.savefig("images/question3_answer.png")
display(fig)
```

```{code-cell} ipython3

```
