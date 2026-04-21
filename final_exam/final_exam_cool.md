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

(final_exam_cool)=
# final exam cooling question

+++

## [5)] (12) Cooling

  Use the tephigram labeled ``cooling problem'' to calculate the following:

  
  For air at 700 hPa with 6 g/kg of vapor (saturated) and 1 g/kg of liquid.


  -  (4) Find
  
     -  The LCL of this air
     -  The approximate temperature if it was brought adiabatically to a pressure of 1000 hPa.

  -  (8) Suppose this air was cooled by 6 degrees C at a constant pressure of 700 hPa.  Find:

     -  The amount of liquid water condensed by the cooling (g/kg)
     -  The new LCL, assuming no precipitation
     -  The amount of energy $\Delta q_{out}$ (J/kg) shed to the environment during the cooling.

### Question 1 cooling code

```{code-cell} ipython3
from a405.thermo.constants import constants as c

from a405.soundings.wyominglib import read_soundings
from a405.skewT.skewlib import makeSkewDry
from a405.thermo.thermlib import convertTempToSkew, find_lcl, find_thetaet, find_Td, tinvert_thetae
from a405.thermo.thermlib import find_rsat, find_theta, find_Tmoist

from a405.skewT.fullskew import makeSkewWet,find_corners,make_default_labels

import datetime
import pytz
import numpy as np
from matplotlib import pyplot as plt

hPa2pa = 100  # convert hPascals to Pascals
```

```{code-cell} ipython3
delta_T = -6 #K
delta_rv = -2.5*1.e-3  #kg/kg
delta_qout = c.cpd*delta_T + c.lv0*delta_rv
print(f"!!!!!!cooling total q = {delta_qout:.1e} J/kg")
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

```

#### Find temperature at 700 hPa

```{code-cell} ipython3
press_700 = 700e2
rv=7e-3
rl = 1.e-3
rt_cloud = rv + rl

Td_700 = find_Td(rv,press_700)
temp_700 = Td_700
thetae_before = find_thetaet(Td_700, rt_cloud, temp_700,press_700)
print(f"{(temp_700 - c.Tc)=:0.1f} deg C")
```

#### Q1a Find the lcl with $\theta$ prior cooling

Plot as a green diamond

```{code-cell} ipython3
press_900 = 900.e2
Temp_900,rv_900,rl_900=tinvert_thetae(thetae_before,rt_cloud,press_900)
# find lcl
Td_900 = find_Td(rt_cloud, press_900)
T_lcl, press_lcl = find_lcl(Td_900,Temp_900,press_900)
print(f"!!!!!!{press_lcl=}")

fig,ax1 =plt.subplots(1,1,figsize=(11,11))
skew = 35
corners=[5,20]
ax1, skew = makeSkewWet(ax1, corners=corners, skew=skew,label_fun=label_fun)
ax1.set_title('Cooling problem')

xcorners=find_corners(corners,skew=skew)
ax1.set(xlim=xcorners,ylim=[1000,600])
fig.savefig('cooling.pdf')
xplot=convertTempToSkew(Td_700 - c.Tc,press_700*pa2hPa,skew)
bot=ax1.plot(xplot, press_700*pa2hPa, 'bd', markersize=14, markerfacecolor='r',label="before cooling")
xplot = xplot=convertTempToSkew(T_lcl - c.Tc,press_lcl*pa2hPa,skew)
# lcl=ax1.plot(xplot, press_lcl*pa2hPa, 'bd', markersize=14, markerfacecolor='g',
#            label = "original lcl")
theta_lcl = find_theta(T_lcl,press_lcl)
print(f"!!!!!LCL potential temperature {theta_lcl:0.1f} K")
print(f"!!!!!LCL pressure {press_lcl*pa2hPa:0.1f} hPa")
ax1.legend();
```

#### Check that $\theta_e$ hasn't changed

```{code-cell} ipython3
thetae_check= find_thetaet(Td_900, rt_cloud, Temp_900,press_900)
thetae_check, thetae_before
```

#### now cool by 6 K

```{code-cell} ipython3
new_temp = temp_700 - 6
thetae_after = find_thetaet(Td_700, rt_cloud, new_temp,press_700)
Temp_after,rv_after,rl_after=tinvert_thetae(thetae_after,rt_cloud,press_700)
print(f"{(thetae_before - thetae_after)=} K")
print(f"!!!!!New rv {rv_after*1.e3:0.1f} g/kg")
print(f"!!!!!rv_after, rv_change: {rv_after*1.e3:0.1f} g/kg, {(rv_after - rv)*1.e3:0.1f} g/kg")
xplot = xplot=convertTempToSkew(Temp_after - c.Tc,press_700*pa2hPa,skew)
lcl=ax1.plot(xplot, press_700*pa2hPa, 'bd', markersize=14, markerfacecolor='b',
           label = "after 6K cooling")
lcl=ax1.plot(xplot, press_lcl*pa2hPa, 'bd', markersize=14, markerfacecolor='g',
         label = "original lcl")
ax1.legend()
display(fig)
```

```{code-cell} ipython3
rv_after
```

#### Find the new lcl

Go down to the surface in case lcl is below 900

```{code-cell} ipython3
press_1000=1.e5 # Pa
Temp_1000,rv_1000,rl_1000=tinvert_thetae(thetae_after,rt_cloud,press_1000)
Td_1000 = find_Td(rt_cloud, press_1000)
T_lcl, press_lcl = find_lcl(Td_1000,Temp_1000,press_1000)
print(f"!!!!!!!!new lcl temperature, pressure {T_lcl:.1f}, {press_lcl:.1f}")
xplot = xplot=convertTempToSkew(T_lcl - c.Tc,press_lcl*pa2hPa,skew)
lcl=ax1.plot(xplot, press_lcl*pa2hPa, 'cd', markersize=14, markerfacecolor='c',
           label = "new lcl")
ax1.legend()
fig.savefig("images/question1_answer.png")
display(fig)
```
