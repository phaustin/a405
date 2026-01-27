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

(assign2_hydro_solution)=
# Assign 2: Solution - Scale heights for typical atmospheric soundings 

Here is the [hydrostatic_scaling_solution.ipynb](https://drive.google.com/file/d/1r2VmOC0yWLP_OyRPLmau9o-8tLrrUt8r/view?usp=sharing)  download link.

+++

## Tasks

1. Get a unique sounding for your assigned region and season
2. Plot the dewpoint and temperature soundings using metpy
3. Write a function to calculate the pressure scale height
4. Write a funcion to calculate the  density scale height
5. Plot the  vertical pressure profile vs. the hydrostatic pressure profile for your pressure scale height

+++

### 1. Getting a sounding

In the cell below, change the region, year, month and station to get a unique sounding for analysis (i.e. I'd like to spread out the soundings among regions and seasons).  The url for the Wyoming page is:  [https://weather.uwyo.edu/upperair/sounding.html](https://weather.uwyo.edu/upperair/sounding.html)

```{code-cell} ipython3
from a405.soundings.wyominglib import write_soundings
import numpy as np
from pathlib import Path
current_dir = Path()
sounding_dir = (current_dir / 'sounding_dir').absolute()
print(f"{str(sounding_dir)=}")
write = False
if write:
    region = 'samer'
    year = '2013'
    month= '7'
    start = '0100'
    stop = '0110'
    station = '71109'
    values=dict(region=region,year=year,month=month,start=start,stop=stop,station=station)
    write_soundings(values, sounding_dir)
```

```{code-cell} ipython3
from a405.soundings.wyominglib import read_soundings
wyoming_dict = read_soundings(sounding_dir)
sounding_dict = wyoming_dict['sounding_dict']
the_sounding = sounding_dict[(2013, 7, 1, 0)]
press = the_sounding['pres'].to_numpy() #hPa
dewpoint = the_sounding['dwpt'].to_numpy() #degC
temp = the_sounding['temp'].to_numpy() # degC
height = the_sounding['hght'].to_numpy() #m
tempK = temp + 273.15
pressPa = press*100.
height_km = height/1000.
```

```{code-cell} ipython3
the_sounding
```

###  2. Plot the sounding

```{code-cell} ipython3
from metpy.plots import SkewT
from metpy.units import units
from matplotlib import pyplot as plt
import numpy as np
fig,ax =plt.subplots(1,1,figsize=(8,8))
fig.clf()
skew_plot = SkewT(fig)
skew_plot.ax.set_title("metpy example")
skew_plot.ax.set(xlim=(0,25),ylim=(1000,600))
theta = np.array([0,10,20,30,40,50,60]) + 273.15
theta = theta*units("K")
skew_plot.plot_dry_adiabats(t0=theta)
skew_plot.plot_moist_adiabats()
skew_plot.plot_mixing_lines()
skew_plot.plot(press,temp,'r')
skew_plot.plot(press,dewpoint,'g');
```

### 3. Calculate the pressure scale height

+++

Here is equation 14 of the [hydrostatic balance notes](https://clouds.eos.ubc.ca/~phil/courses/atsc405/docs/hydro.pdf)

+++

$$\frac{ 1}{\overline{H_p}} =  \overline{ \left ( \frac{1 }{H} \right )} = \frac{\int_{0 }^{z}\!\frac{1}{H} dz^\prime  }{z-0} $$

where

$$H=R_d T/g$$

+++

a. Turn this into a python function called calc_scale_height that takes sounding vectors of temperature, pressure and height and returns
   the pressure scale height in meters

b.  Use this to find the pressure scale height in meters  as a scalar variable named Hbar

```{code-cell} ipython3
#
#  in this cell define a function called
#
#  calc_scale_height(T,p,z)
#  which takes vertical profiles of temperature, pressure and height and calculates the
#  pressure scale height in meters
#

g=9.8  #don't worry about g(z) for this exercise
Rd=287.  #kg/m^3
def calc_scale_height(T,p,z):
    """
    Calculate the pressure scale height H_p
    
    Parameters
    ----------
    
    T: vector (float)
      temperature (K)
      
    p: vector (float) of len(T)
      pressure (pa)
      
    z: vector (float) of len(T
      height (m)
      
    Returns
    -------
    
    Hbar: vector (float) of len(T)
      pressure scale height (m)
    
    """
    dz=np.diff(z)
    TLayer=(T[1:] + T[0:-1])/2.
    oneOverH=g/(Rd*TLayer)
    Zthick=z[-1] - z[0]
    oneOverHbar=np.sum(oneOverH*dz)/Zthick
    Hbar = 1/oneOverHbar
    return Hbar
```

### 3. Answer

```{code-cell} ipython3
Hbar = calc_scale_height(tempK,pressPa,height)
print(f"{Hbar=:.0f} meters")
```

### 4. Calculate the density scale height

Similarly, equation (23) of the [hydrostatic balance notes](https://clouds.eos.ubc.ca/~phil/courses/atsc500/docs/hydro.pdf)
is:

$$\frac{d\rho }{\rho} = - \left ( \frac{1 }{H} + 
                   \frac{1 }{T} \frac{dT }{dz} \right ) dz \equiv - \frac{dz }{H_\rho} $$
                   
Which leads to 

$$\frac{ 1}{\overline{H_\rho}}  = \frac{\int_{0 }^{z}\!\left [ \frac{1}{H} + \frac{1 }{T} \frac{dT }{dz} \right ] dz^\prime  }{z-0} $$

and the following python function:

```{code-cell} ipython3
#
#  in this cell define a function called
#
#  calc_dense_height(T,p,z)
#  which takes vertical profiles of temperature, pressure and height and calculates the
#  density scale height in meters
#

def calc_dense_height(T,p,z):
    """
    Calculate the density scale height H_rho
    
    Parameters
    ----------
    
    T: vector (float)
      temperature (K)
      
    p: vector (float) of len(T)
      pressure (pa)
      
    z: vector (float) of len(T
      height (m)
      
    Returns
    -------
    
    Hbar: vector (float) of len(T)
      density scale height (m)
    """
    dz=np.diff(z)
    TLayer=(T[1:] + T[0:-1])/2.
    dTdz=np.diff(T)/np.diff(z)
    oneOverH=g/(Rd*TLayer) + (1/TLayer*dTdz)
    Zthick=z[-1] - z[0]
    oneOverHbar=np.sum(oneOverH*dz)/Zthick
    Hbar = 1/oneOverHbar
    return Hbar
```

### 4. Answer

```{code-cell} ipython3
Hrho = calc_dense_height(tempK,pressPa,height)
print(f"{Hrho=:.0f} meters")
```

### 5. How does the hydrostatic profile compare to the observed pressure sounding?

+++

Now check the hydrostatic approximation by plotting the pressure column against

$$p(z) = p_0 \exp \left (-z/\overline{H_p} \right )$$

vs. the actual sounding p(T):

+++

### 5. Answer

```{code-cell} ipython3
#
#  In this cell make a plot with two lines, 
#
#  line 1: press (hPa) vs. height_km (vertical axis)
#  line 2: press p(z) from equation 5 (hPa) vs. height_km
#
#  include a legend that clear identifies the two lines
#  set the y axis limits from 0 to 10 km
#
fig,theAx=plt.subplots(1,1)
hydroPress=press[0]*np.exp(-height/Hbar)
theAx.plot(press,height_km,label='sounding')
theAx.plot(hydroPress,height_km,label='hydrostat approx')
theAx.set_title('height vs. pressure')
theAx.set_xlabel('pressure (hPa)')
theAx.set_ylabel('height (km)')
theAx.set_xlim([200,1000])
theAx.set_ylim([0,10])
theAx.grid(True)
theAx.legend(loc='best');
```

```{code-cell} ipython3

```
