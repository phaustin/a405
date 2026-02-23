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

(assign4_paper_solution)=
# Assignment 4 part 2  -- solution

These questions will be due the weekend of February 7/8 

Upload a scanned image with your answer to the following questions

You can use the same tephigram (in separate copies) for questions 1/2 and then question 3.  

Here is the link to [assign4b_tephi](https://docs.google.com/presentation/d/1UMxOLz9VCRBHDpUvJOuWqgN6eFeSkVGErlthLPhbDVw/edit?usp=sharing)

+++

## Question 1

-  The linked tephigram shows the temperature (circle) and dewpoint (diamond) for air at 900 hPa and 800 hPa Find (showing your work on the tephigram):

      - The LCL (hPa) of the 900 hPa air
      - The equivalent potential temperature $\theta_e$ (K) 
      - The entropy $\phi$ (J/kg/K) or the 900 hPa air
      - The wet bulb potential temperature $\theta_w$ (K) of the 900 hPa air 

- Is the air at 900 hPa absolutely stable, conditionally unstable, or unstable?  Explain

- Is the layer between 900 hPa and 800 hPa convectively unstable? Explain.

+++

### Question 1 answer

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
    tempLabels = range(-40, 50, 5)
    rsLabels = [0.1, 0.25, 0.5, 1, 2, 3] + list(np.arange(4, 28, 1)) 
    return tempLabels,rsLabels, thetaLabels, thetaeLabels
```

```{code-cell} ipython3
from matplotlib import pyplot as plt
import numpy as np
import a405.thermo.thermlib as tl
from a405.skewT.skewlib import makeSkewDry
from a405.thermo.thermlib import convertTempToSkew
from a405.skewT.fullskew import makeSkewWet,find_corners,make_default_labels
from a405.thermo.constants import constants as c
pa2hPa = 1.e-2
fig,ax =plt.subplots(1,1,figsize=(13,8))
corners = [10, 30]
ax, skew = makeSkewWet(ax, corners=corners, skew=28,label_fun=label_fun)
ax.set_title('Assign4b tephigram')
xcorners=find_corners(corners,skew=skew)
ax.set(xlim=xcorners,ylim=[1000,700])
#
# At 900 hPa set Tdew = 16 K and Temp = 19 K
#
Tdew_900 = 15.8 + c.Tc
Temp_900 = 18.7 + c.Tc
press_900 = 900e2
rv_900=tl.find_rsat(Tdew_900,press_900)
#
# What is the lcl of this air?
#
T900_lcl,p900_lcl = tl.find_lcl(Tdew_900,Temp_900,press_900)
print(f"{T900_lcl=:.1f} K, Question 1a answer: {p900_lcl=:.1f} Pa")
# #
# #  draw these on the sounding at 900 hPa as a red circle and blue diamond
# #
xplot=convertTempToSkew(Temp_900 - c.Tc,press_900*pa2hPa,skew)
bot=ax.plot(xplot, press_900*pa2hPa, 'ro', markersize=14, markerfacecolor='r')
xplot=convertTempToSkew(Tdew_900 - c.Tc,press_900*pa2hPa,skew)
bot=ax.plot(xplot, press_900*pa2hPa, 'bd', markersize=14, markerfacecolor='b')
#
# now do the 800 hPa air with Temp = 16 deg C and Tdew = 4 deg c
#
Tdew_800 =  6.8 + c.Tc
Temp_800 =  16. + c.Tc
press_800 = 800e2
rv_800=tl.find_rsat(Tdew_800,press_900)
T800_lcl,p800_lcl = tl.find_lcl(Tdew_800,Temp_800,press_800)
print(f"{T800_lcl=:.1f} K,{p800_lcl=:.1f} Pa")
#
# What is the lcl of this air?
#
press_700=700.e2
thetae_800=332.  #K
Temp_700=tl.find_Tmoist(thetae_800,press_700)
Tdew_700 = Temp_700
rv_700=tl.find_rsat(Temp_700,press_700)
press_800=800.e2
rv_800=rv_700   #total water is conserved
#Temp_800,rv_800,rl=tl.tinvert_thetae(thetae_800,rv_800,press_800)
#Tdew_800=tl.find_Td(rv_800,press_800)
print(f"{Temp_800-c.Tc=:.1f} K,{Tdew_800-c.Tc=:.1f} Pa")
#
# put these points on the sounding at 800 hPa
#
xplot=convertTempToSkew(Temp_800 - c.Tc,press_800*pa2hPa,skew)
bot=ax.plot(xplot, press_800*pa2hPa, 'ro', markersize=14, markerfacecolor='r')
xplot=convertTempToSkew(Tdew_800 - c.Tc,press_800*pa2hPa,skew)
bot=ax.plot(xplot, press_800*pa2hPa, 'bd', markersize=14, markerfacecolor='b')
xplot = convertTempToSkew(T900_lcl - c.Tc,p900_lcl*pa2hPa,skew)
lcl_900 = ax.plot(xplot, p900_lcl*pa2hPa, 'kv', markersize=14, markerfacecolor='k')
fig.savefig('assign4b_tephi.jpg')
```

## Question 2

Suppose air at 900 hPa is lifted adiabatically to 800 hPa, where it mixes 50/50 with
air at that level.

Using the tephigram what is the LCL of the resulting mixture?  -- Clearly label your work, including the tephigram.


## Question 3

- Use the linked tephigram labeled  to sketch the following thermodynamic
cycle.  You have a sealed cannister containing **1 kg of dry air and 14 grams of water**.  Initially it is at a pressure of **900 hPa and a temperature of 25 deg C** (Point A).  It then is taken through the following stages:

  - adiabatic expansion until the internal pressure is 700 hPa to Point B
  - isothermal compression back to 900 hPa (Point C)
  - isobaric (constant pressure) heating back to 900 hPa  (Point A)

Find the following (you can use $l_v = l_{v0}$ and $c_p = c_{pd}$) :

- The equivalent potential temperatures of point A and point C (Kelvins)
- The heating of the cannister $Q_{in}$ (during C $\rightarrow$ A) and 
cooling $Q_{out}$) during B $\rightarrow$ C  (in Joules/kg)
- The efficiency of this heat engine (in percent)

```{code-cell} ipython3

```
