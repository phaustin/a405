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
---
jupyter:
  outputs_hidden: false
---
from a405.utils.helper_funs import formatvar as f
#
# now make one of these functions with the default 5.2f format
#
n = f()
#
# and another function with exponential format
#
e = f(fmtstring = '{:6.3e}')
#
#  use n() below like this
#
print('here0: {} K'.format(n(273.14156)))
print('here1: {} K'.format(e(273.14156)))
#
#write a regular function to convert from K to deg C
#
from a405.thermo.constants import constants as c
def k2c(temp): return (temp - c.Tc)
print('here2: {} deg C'.format(n(k2c(300))))
```

## back to midterm question 2:  make the tephigram -- with LCLs as black dots

```{code-cell} ipython3
---
jupyter:
  outputs_hidden: false
---
pa2hPa=1.e-2
from a405.skewT.fullskew import makeSkewWet, find_corners
from matplotlib import pyplot as plt
import numpy as np
fig, ax = plt.subplots(1, 1, figsize=(12, 8))
corners=[10,30]
ax, skew = makeSkewWet(ax,corners=corners)
ax.set(ylim=[1000,700])
from a405.thermo.thermlib import find_Tmoist,find_rsat,find_Td,tinvert_thetae,convertTempToSkew,\
              find_lcl
import a405.thermo.thermlib as tl
#
# set the lcl for 900 hPa to 860 hPa and thetae to 338 K
#
press=860.e2
thetae_900=338.  #K
Temp_860=find_Tmoist(thetae_900,press)
rv_860=find_rsat(Temp_860,press)
rv_900 = rv_860  #vapor is conserved
Tdew_860=Temp_860
print("temp,Tdew,rv at LCL press:  {} hPa {} C {} C {} kg/kg"\
      .format(n(press*1.e-2),e(k2c(Temp_860)),e(k2c(Tdew_860)),e(rv_900)))
#
# now descend adiabatically to 900 hPa
#
press=900.e2
Temp_900,rv_900,rl_900=tinvert_thetae(thetae_900,rv_900,press)
Tdew_900=find_Td(rv_900,press)
print("temperature and dewpoint at {} hPa hPa = {} C {} C".format(n(press*1.e-2),e(k2c(Temp_900)), e(k2c(Tdew_900))))
#
#  draw these on the sounding at 900 hPa as a red circle and blue diamond
#
xplot=convertTempToSkew(Temp_900 - c.Tc,press*pa2hPa,skew)
bot=ax.plot(xplot, press*pa2hPa, 'ro', markersize=14, markerfacecolor='r')
xplot=convertTempToSkew(Tdew_900 - c.Tc,press*pa2hPa,skew)
bot=ax.plot(xplot, press*pa2hPa, 'bd', markersize=14, markerfacecolor='b')
#
#now look at an LCL of 700 hPa  -- with thetae = 332 K
#
press=700.e2
thetae_800=332.  #K
Temp_700=find_Tmoist(thetae_800,press)
Tdew_700 = Temp_700
rv_700=find_rsat(Temp_700,press)

print("temp,Tdew, rv at LCL press: {} hPa: {} C  {} C {} kg/kg"\
      .format(n(press*1.e-2),e(k2c(Temp_700)),e(k2c(Tdew_700)),e(rv_700)))
# get the temperature and dewpoint at 800 hPa
#
press=800.e2
rv_800=rv_700   #total water is conserved
Temp_800,rv_800,rl=tinvert_thetae(thetae_800,rv_800,press)
Tdew_800=find_Td(rv_800,press)
print("temperature and dewpoint at {} hPa = {} C {} C".format(n(press*1.e-2),n(k2c(Temp_800)),n(k2c(Tdew_800))))
#
# put these points on the sounding at 800 hPa
#
xplot=convertTempToSkew(Temp_800 - c.Tc,press*pa2hPa,skew)
bot=ax.plot(xplot, press*pa2hPa, 'ro', markersize=14, markerfacecolor='r')
xplot=convertTempToSkew(Tdew_800 - c.Tc,press*pa2hPa,skew)
bot=ax.plot(xplot, press*pa2hPa, 'bd', markersize=14, markerfacecolor='b')
#
# draw the two LCLs as black circles
#
press=860.e2
xplot=convertTempToSkew(Temp_860 - c.Tc,press*pa2hPa,skew)
bot=ax.plot(xplot, press*pa2hPa, 'ko', markersize=14, markerfacecolor='k')
press=702.e2  #add 2 hPa so we can see it
xplot=convertTempToSkew(Temp_700 - c.Tc,press*pa2hPa,skew)
bot=ax.plot(xplot, press*pa2hPa, 'ko', markersize=14, markerfacecolor='k')
fig.savefig('mid-tephi.pdf')


def find_thetaes(temp,press,rv):
    """ input temp in K, press in Pa, rv in kg/kg
        output thetae in K
    """
    thetaes = temp*(c.p0/press)**(c.Rd/c.cpd)*np.exp(c.lv0*rv/(c.cpd*temp))
    return thetaes

press=700.e2
thetaes_700 = find_thetaes(Temp_700,press,rv_700)
accu_thetaes_700 = tl.find_thetaes(Temp_700,press)
print('approx thetaes at {} = {} K'.format(press*1.e-2,n(thetaes_700)))
print('accurate thetaes at {} = {} K'.format(press*1.e-2,n(accu_thetaes_700)))

press=860.e2
thetaes_860 = find_thetaes(Temp_860,press,rv_860)
accu_thetaes_860 = tl.find_thetaes(Temp_860,press)
print('approx thetaes at {} = {} K'.format(press*1.e-2,n(thetaes_860)))
print('accu thetaes at {} = {} K'.format(press*1.e-2,n(accu_thetaes_860)))
print('entropy for 860 hPa = {}'.format(e(c.cpd*np.log(thetaes_860))))

#
#wet bulb temp potential temperature for 900 hPa  -- bring air to 1000 hPa
#along a moist adiabat
#
press = 1.e5
Temp_1000=find_Tmoist(thetae_900,press)
print('wet bulb potential temperature for 900 hPa air = {} C'.format(n(k2c(Temp_1000))))
#
# add as a green circle
#
press=998.e2  #subtrack 2 hPa so we can see it
xplot=convertTempToSkew(Temp_1000 - c.Tc,press*pa2hPa,skew)
bot=ax.plot(xplot, press*pa2hPa, 'go', markersize=14, markerfacecolor='g')
xcorners=find_corners(corners,skew=skew)
ax.set(xlim=xcorners,ylim=[1000,700])
```

```{code-cell} ipython3
help(find_corners)
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
