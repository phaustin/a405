---
jupytext:
  formats: md:myst,ipynb
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

(assign4_rootfind_solution)=
# Assignment 4a rootfinding

Download link to [rootfind.ipynb](https://drive.google.com/file/d/1HBjqe9Y8WbUPuOts2RQtia2d7w8VE96A/view?usp=sharing)

Using a rootfinder to solve implicit equations

+++

We'll need to use rootfinding to solve equations like finding the lifting condensation level, i.e. the level
where the dewpoint=temperature.  This notebook works through a simple example and sets up question 1, which asks
you to use a rootfinder to find the dewpoint temperature for a given saturation vapor pressure

+++

Documentation is here [scipy.optimize.brentq](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html)

+++

The cell below shows how to use a rootfinder to solve the equation

$$
cos(x) - 0.75 = 0
$$

The algorithm evaluates the function for increasing values of x until in sees a sign change.  It then backtracks to the other side of the sign change and starts increaing x again until it sees the sign change.  It will stop when the change in x is below a selected tolerance

```{code-cell} ipython3
---
jupyter:
  outputs_hidden: false
---
from scipy import optimize
from matplotlib import pyplot as plt
import numpy as np
plt.style.use('ggplot')


xvals=np.linspace(0,10.)
fig,ax = plt.subplots(1,1,figsize=(10,8))
ax.plot(xvals,np.cos(xvals))
straight_line=np.ones_like(xvals)
ax.plot(xvals,straight_line*0.75,'b-')

def root_function(x):
    """Function we want to find the root of
       input: x value
       output: y value -- should be zero when x is a root
    """
    return np.cos(x) - 0.75

root1 = optimize.brentq(root_function,0,2)
root2 = optimize.brentq(root_function,4,6)
root3 = optimize.brentq(root_function,6,8)
xvals=np.array([root1,root2,root3])
yvals=np.cos(xvals)
ax.plot(xvals,yvals,'ro',markersize=20);
```

## Problem 1 (part of assignment 4)

Giventhe Clausius Clapeyron equation and a saturation vapor pressure  $esat$ in Pa, find the dewpoint temperature in Kelvin, i.e. the temperature where:

    esat = find_esat(Tdewpoint)

```{code-cell} ipython3
def find_esat(temp):
    """
    Thompkins 2.15 version of Clausius-Clapeyron equation
    
    Parameters
    ----------

    temp: temperature (K)

    Returns
    -------

    esatOut: saturation vapor pressure (Pa)
    """
    Tc = temp - 273.15
    esatOut = 611.2 * np.exp(17.67 * Tc / (Tc + 243.5))
    return esatOut
```

In the cell below, use the rootfinder to write this function:

    def find_dewpoint(esat):
       ...
       return Tdew

and use it to make a plot of dewpoint vs temperature from 0 deg C to 30 deg C

```{code-cell} ipython3
# you code here
def root_function(Tdew,esat_target):
    return esat_target - find_esat(Tdew)
    
def find_dewpoint(esat):
    Tdew = optimize.brentq(root_function,200,350,args=(esat,))
    return Tdew
    
find_dewpoint(1200)   
```

```{code-cell} ipython3
from a405.thermo.constants import constants as c
import numpy as np
esat_vals = np.linspace(500, 4000, 50)
dewpoint = [find_dewpoint(esat) for esat in esat_vals]
dewpoint_tc = np.array(dewpoint) - c.Tc
esat_vals_mb = np.array(esat_vals)*0.01

plt.style.use('ggplot')
fig, ax = plt.subplots(1,1,figsize = [8,6])

# Plot values and convert the saturation vapour pressure to hPa
ax.plot(esat_vals_mb,dewpoint_tc)
ax.set_xlabel("Saturation vapor pressure (hPa)")
ax.set_ylabel("Dewpoint temperature (°C)")
plt.grid(True)
```

```{code-cell} ipython3

```
