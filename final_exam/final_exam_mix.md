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

(final_exam_mix)=
# final exam mixing question

+++

## [6)] (12) Mixing

  Surface air at 1000 hPa with a temperature of 20 deg C and a dewpoint of 16 deg C is lifted adiabatically to 800 hPa, where it entrains 70\%  of environmental air that has a $\theta_e$ of 307 K and a mixing ratio of 4 g/kg.  Use the tephigram labeled ``mixing problem'' to find:

  -  The $\theta_e$ and LCL of the of the surface air
  -  The $\theta_e$, $r_v$ and $r_l$ of the mixture
  -  The LCL of the mixture
  -  The temperature of the mixture at 800 hPa

   Is the mixture negatively, positively, or neutrally buoyant with its surrounding environment?  Explain.

```{code-cell} ipython3
from a405.thermo.constants import constants as c
from a405.thermo.thermlib import find_rsat, find_theta, find_Tmoist
from a405.thermo.thermlib import convertTempToSkew, find_lcl, find_thetaet, find_Td, tinvert_thetae
from matplotlib import pyplot as plt
from a405.skewT.fullskew import makeSkewWet,find_corners,make_default_labels
import numpy as np
from pathlib import Path
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
image_dir = Path() / "images"
if not image_dir.exists():
   image_dir.mkdir(parents=True, exist_ok=True)
```

```{code-cell} ipython3
thetae_cloud = 324
thetae_env = 307
rv_cloud = 11.5e-3
rv_env = 4.e-3
thetae_mix = 0.7*thetae_cloud + 0.3*thetae_env
rv_mix = 0.7*rv_cloud + 0.3*rv_env
print(f"!!!!!!!{thetae_mix=:.1f} K, {rv_mix=:.1e} kg/kg")
```

### Question 3 code

+++

#### get $\theta_e$ of the surface air

```{code-cell} ipython3
temp_1000 = 20 + c.Tc
press_1000 = 1.e5 
Td_1000 = 16 + c.Tc
rt_cloud = find_rsat(Td_1000,press_1000)
thetae_cloud = find_thetaet(Td_1000, rt_cloud, temp_1000,press_1000)
T_lcl, press_lcl = find_lcl(Td_1000,temp_1000,press_1000)
print(f"!!!!!!!thetae: {thetae_cloud:0.1f} K, r_T {rt_cloud*1.e3:0.1f} g/kg, lcl_press {press_lcl:0.1f} Pa")
```

#### Plot T, Td and the lcl

```{code-cell} ipython3
fig,ax2 =plt.subplots(1,1,figsize=(11,11))
skew = 35
corners=[5,21]
ax2, skew = makeSkewWet(ax2, corners=corners, skew=skew,label_fun=label_fun)
ax2.set_title('mixing problem')
xcorners=find_corners(corners,skew=skew)
ax2.set(xlim=xcorners,ylim=[1005,700])
xplot=convertTempToSkew(temp_1000 - c.Tc,press_1000*pa2hPa,skew)
ax2.plot(xplot, press_1000*pa2hPa, 'gd', markersize=14, markerfacecolor='g',
           label = "temp_sfc")
xplot=convertTempToSkew(Td_1000 - c.Tc,press_1000*pa2hPa,skew)
ax2.plot(xplot, press_1000*pa2hPa,'bd', markersize=14, markerfacecolor='b',
           label = "dewpoint_sfc")
xplot=convertTempToSkew(T_lcl - c.Tc,press_lcl*pa2hPa,skew)
ax2.plot(xplot, press_lcl*pa2hPa,'kd', markersize=14, markerfacecolor='k',
           label = "lcl_cloud");
```

#### Lift to 800

```{code-cell} ipython3
press_800 = 8.e4 #Pa
Temp_800,rv_800,rl_800=tinvert_thetae(thetae_cloud,rt_cloud,press_800)
xplot=convertTempToSkew(Temp_800 - c.Tc,press_800*pa2hPa,skew)
ax2.plot(xplot, press_800*pa2hPa,'rd', markersize=14, markerfacecolor='r',
           label = "cloud_800");
```

```{code-cell} ipython3
display(fig)
```

#### add the environment

```{code-cell} ipython3
rt_env = 4.e-3
thetae_env = 307.
Td_env =  find_Td(rt_env, press_800)
Temp_env,rv_env,rl_env=tinvert_thetae(thetae_env,rt_env,press_800)
T_lclenv, press_lclenv = find_lcl(Td_env,Temp_env,press_800)
print(f"!!!!{press_lclenv=:.1f} Pa")
```

```{code-cell} ipython3
xplot=convertTempToSkew(Temp_env - c.Tc,press_800*pa2hPa,skew)
ax2.plot(xplot, press_800*pa2hPa,'gs', markersize=14, markerfacecolor='g',
           label = "temp_env")
xplot=convertTempToSkew(Td_env - c.Tc,press_800*pa2hPa,skew)
ax2.plot(xplot, press_800*pa2hPa,'bs', markersize=14, markerfacecolor='b',
           label = "Td_env")
xplot=convertTempToSkew(T_lclenv - c.Tc,press_lclenv*pa2hPa,skew)
ax2.plot(xplot, press_lclenv*pa2hPa,'ks', markersize=14, markerfacecolor='k',
           label = "lcl_env")
ax2.legend()
display(fig)
```

(mixture_answer)=
#### add the mixture

```{code-cell} ipython3
frac_env = 0.7

thetae_mix = 0.3*thetae_env + 0.7*thetae_cloud
rt_mix = 0.3*rt_env + 0.7*rt_cloud
Temp_mix,rv_mix,rl_mix=tinvert_thetae(thetae_mix,rt_mix,press_800)
Temp_1000,rv_1000,rl_1000=tinvert_thetae(thetae_mix,rt_mix,press_1000)
print(rv_1000)
Td_1000 = find_Td(rv_1000,press_1000)
print(Td_1000 - c.Tc)
T_lclmix, press_lclmix = find_lcl(Td_1000,Temp_1000,press_1000)
print(f"{thetae_mix=:0.1f} K, {rt_mix=:0.4f} kg/kg")
print(f"{Temp_mix - c.Tc=:0.1f} K, {rv_mix=:0.4f} kg/kg, {rl_mix=:0.4f} kg/kg")
print(f"{press_lclmix=:0.1f} Pa")
```

```{code-cell} ipython3
xplot=convertTempToSkew(Temp_mix - c.Tc,press_800*pa2hPa,skew)
ax2.plot(xplot, press_800*pa2hPa,'co', markersize=14, markerfacecolor='c',
           label = "mixture")
xplot=convertTempToSkew(T_lclmix - c.Tc,press_lclmix*pa2hPa,skew)
ax2.plot(xplot, press_lclmix*pa2hPa,'cs', markersize=14, markerfacecolor='c',
           label = "mixture lcl")
ax2.legend()
fig.savefig("images/question3_answer.png")
display(fig)
```
