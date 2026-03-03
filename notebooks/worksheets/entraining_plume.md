---
jupytext:
  formats: md:myst,ipynb
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

(entrain_cloud)=
# Modeling an entraining cloud updraft

This notebook calculates the time evolution of four variables:

\[velocity, height, $\theta_{ecld}$, $r_{Tcloud}$ \]

in a rising, entraining cloud with constant entrainment rate.  The environment is specified from a Wyoming sounding,
and interpolated each timestep using [scipy.interpolate.interp1d](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html) 

The variables are intergrated with respect to time using [scipy.integrate.RK45](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html)

```{code-cell} ipython3
import numpy as np
import pandas as pd
from functools import partial
from pprint import pformat
from a405.thermo.constants import constants as c
from a405.thermo.thermlib import find_Tmoist,find_rsat,find_thetaep
from scipy.interpolate import interp1d
from scipy.integrate import RK45
from a405.soundings.wyominglib import write_soundings, read_soundings
import json
import xarray as xr

import matplotlib.pyplot as plt
from a405.skewT.nudge import nudge
```

## Find the derivatives wrt time of each of the 4 variables

See [entrain.pdf](https://drive.google.com/file/d/1VgY5LWMSmavpoHFQYFE-wkPxGnD5x44O/view?usp=sharing)  notes.  We want to integrate to find $[\theta_e, r_T, w_{vel}, z]$ for the ascending cloud parcel.

The equations:

$$
\frac{d\,\theta_e}{dt} &= \lambda \left (\theta_{e\,e} - \theta_{e,c} \right )\\
\frac{d\,r_T}{dt} &= \lambda \left (r_{T\,e} - r_{T\,c} \right ) \\
\frac{ d\,w_{vel}}{dt} & = \text{buoyant acceleration} \\
\frac{dz}{dt} & = w_{vel}
$$
Note that we really should be mixing $\log \theta_e$, not $\theta_e$ (why?).  We'll see whether that makes a difference later.

For the entrainment calculation, we need environmental temperature, dewpoint and pressure at arbitrary heights.  Use [scipy.interpolate.interp1d](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html) for this

```{code-cell} ipython3
def derivs(t, y, entrain_rate=None, tinterp=None, tdinterp = None, pinterp=None):
    """Function that computes the derivative vector for the ode integrator

    Parameters
    ----------
    
    t: float
       time (s)
    y: vector
       4-vector containing wvel (m/s), height (m), thetae (K), rt (kg/kg)
    entrain_rate: float
                  1/m dm/dt (s-1)
    tinterp: func
                interp1d function for environmental temperature T(z) 
    tdinterp: func
                interp1d function for environmental dewpoint temperature Td(z)
    pinterp: func
                interp1d function for presusure  p(z)

    Returns
    -------

    yp: vector
       4-vector containing time derivatives of wvel (m/s^2), height (m/s), thetae (K/s), rt (kg/kg/s)
    """
    #print(f"inside derivs")
    #print(f"{y=}")
    yp = np.zeros((4,),dtype=float)
    velocity = y[0]
    height = y[1]
    thetae_cloud = y[2]
    rt_cloud = y[3]
    #
    # fill the yp (yprime) vector with the derivatives
    #
    # yp[0] is the acceleration, in this case the buoyancy 
    #
    yp[0] = calcBuoy(height, thetae_cloud, tinterp, tdinterp, pinterp)
    press = pinterp(height)*100. #Pa
    Tdenv = tdinterp(height) + c.Tc #K
    Tenv = tinterp(height) + c.Tc #K
    rtenv = find_rsat(Tdenv, press) #kg/kg
    thetaeEnv = find_thetaep(Tdenv,Tenv, press)
    #
    # yp[1] is the rate of change of height
    #
    yp[1] = velocity
    #
    # yp[2] is the rate of change of thetae_cloud
    #
    yp[2] = entrain_rate*(thetaeEnv - thetae_cloud)
    #
    # yp[3] is the rate of change of rt_cloud
    #
    yp[3] = entrain_rate*(rtenv - rt_cloud)
    #print(f" derivs returning")
    #print(f"{yp=}")
    return yp
```

## Find the buoyancy from the cloud and environment $\theta_e$ and $r_T$

```{code-cell} ipython3
def calcBuoy(height, thetae0, tinterp, tdinterp, pinterp):
    """function to calculate buoyant acceleration for an ascending saturated parcel
       this version neglects liquid water loading
    
    Parameters
    ----------
    
    height: float
            parcel height (m)
    thetae0: float
            parcel thetae (K)

    tinterp: func
                interp1d function for environmental temperature T(z) 
    tdinterp: func
                interp1d function for environmental dewpoint temperature Td(z)
    pinterp: func
                interp1d function for presusure  p(z)

    Returns
    -------

    buoy: float
          buoyancy (m/s/s)
    """
    #input: height (m), thetae0 (K), plus function handles for
    #T,Td, press soundings
    #output: Bout = buoyant acceleration in m/s^2
    #neglect liquid water loading in the virtual temperature
    
    press=pinterp(height)*100.#%Pa
    Tcloud=find_Tmoist(thetae0,press) #K
    rvcloud=find_rsat(Tcloud,press); #kg/kg
    Tvcloud=Tcloud*(1. + c.eps*rvcloud)
    Tenv=tinterp(height) + c.Tc
    Tdenv=tdinterp(height) + c.Tc
    #print(f"inside buoy {(height,Tenv,Tdenv,press)=}")
    rvenv=find_rsat(Tdenv,press); #kg/kg
    Tvenv=Tenv*(1. + c.eps*rvenv)
    TvDiff=Tvcloud - Tvenv
    buoy = c.g0*(TvDiff/Tvenv)
    return buoy
```

## Integrator 


Use [scipy.integrate.RK45](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.RK45.html) to integrate our system of 4 odes

```{code-cell} ipython3
def integ_entrain(df_sounding,entrain_rate,wvel_init,press_init,integ_time=2500,time_interval=10):
    """integrate an ascending parcel given a constant entrainment rate
       this version hardwired to start parcel at 800 hPa with cloud base
       values of environment at 900 hPa

    Parameters
    ----------

    df_sounding: pandas dataframe 
               : cloumns are temperature, dewpoint, height, press

    entrain_rate: float
                  1/m dm/dt (s-1)

    wvel_init: float
             initial vertical velocity (m/s)

    press_init: float
             initial parcel pressure level (hPa)

    integ_time: float (optional, default=2500 s)
             total integration time (s)
             
    time_interval: float (optional, default = 10 s)
          recording interval for output

    Returns
    -------

    df_out: dataframe
          dataframe containing wvel (m/s) ,cloud_height (m) , thetae (K), rt (kg/kg) for assending parcel

   interpPress: func
              interp1d function for presusure  p(z) (used for plotting)
    """
    press = df_sounding['pres'].values
    height = df_sounding['hght'].values
    temp = df_sounding['temp'].values
    dewpoint = df_sounding['dwpt'].values
    #
    # 
    #
    envHeight= nudge(height)

    interpTenv = interp1d(envHeight,temp)
    interpTdEnv = interp1d(envHeight,dewpoint)
    interpPress = interp1d(envHeight,press)

    args=dict(entrain_rate = entrain_rate,
              tinterp=interpTenv,
              tdinterp=interpTdEnv,
              pinterp=interpPress)
    #
    # use functools.partial to supply the extra arguments to the derivs function
    # this changes the derivs function signature from
    #
    # derivs(t, y, entrain_rate=None, tinterp=None, tdinterp = None, pinterp=None)
    #
    # to 
    #
    # the_derivs(t,y) 
    #
    # as required by the RK45 integrator
    # 
    #
    the_derivs = partial(derivs,**args)
    #
    #  find the initial toral water mixing ratio and thetae at press_init
    #
    pinit_level = len(press) - np.searchsorted(press[::-1],press_init)
    rtcloud = find_rsat(dewpoint[pinit_level] + c.Tc, press[pinit_level]*100.)
    thetaeVal=find_thetaep(dewpoint[pinit_level] + c.Tc,
                           temp[pinit_level] + c.Tc,press[pinit_level]*100.)
    #
    # start the parcel with initial vertical velocity wvel_init
    #
    #
    height_init=height[pinit_level]
    yinit = [wvel_init, height_init, thetaeVal, rtcloud]  
    tinit = 0  #seconds
    dt = time_interval #seconds
    output_times = np.arange(tinit, integ_time, dt)
    #
    # want to integrate derivs using dopr15 runge kutta described at
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.ode.html
    #
    init_vals = (yinit, tinit)
    rk45 = RK45(the_derivs, tinit, yinit, t_bound=integ_time, max_step=dt)
    yout = []
    while rk45.status == "running":
        #
        # stop if wvel < 0
        #
        if rk45.y[0] < 0:
            break
        press = interpPress(rk45.y[1])*100.
        yvals= [rk45.t,press]
        yvals.extend(rk45.y)
        yout.append(yvals)
        rk45.step()
    #
    # convert the output into a pandas dataframe
    #
    colnames=['time','press','wvel','cloud_height','thetae_cloud','rt_cloud']
    df_out=pd.DataFrame.from_records(yout,columns=colnames)
    df_out = df_out.set_index('time',drop=False)
    return df_out
```

## Read in a sounding to set the environment

```{code-cell} ipython3
from a405.soundings.wyominglib import write_soundings
from pathlib import Path
curr_dir = Path()
stem = 'little_rock_july'
sounding_dir = (curr_dir / stem ).absolute()
station = '72340'
year = '2012'
month = '7'
values=dict(region='naconf',year=year,month=month,start='0100',stop='3000',station=station)
write = False
if write:
    write_soundings(values, sounding_dir)
    soundings= read_soundings(sounding_dir)
else:
    soundings= read_soundings(sounding_dir)
```

```{code-cell} ipython3
day = 9
hour = 0
the_time=(2012,7,day,hour)
output_filename = curr_dir / f"little_rock_july_{day:d}.nc"
sounding=soundings['sounding_dict'][the_time]
```

## Do the integration

```{code-cell} ipython3
entrain_rate = 8.e-4  #s^{-1}
wvel_init = 10
press_init = 900
df = integ_entrain(sounding,entrain_rate,wvel_init,press_init,integ_time=3500,time_interval=3)
```

```{code-cell} ipython3
df
```

## Convert the dataframe to xarray

```{code-cell} ipython3
the_ds = df.to_xarray()
the_ds
```

## Add units to the variables plus dataset attributes

```{code-cell} ipython3
the_ds = the_ds.set_coords(['cloud_height','press'])
the_ds['press'] = the_ds['press'].assign_attrs(units = 'Pa')
the_ds['cloud_height'] = the_ds['cloud_height'].assign_attrs(units = 'm')
the_ds['wvel'] = the_ds['wvel'].assign_attrs(units = 'm/s')
the_ds['thetae_cloud'] = the_ds['thetae_cloud'].assign_attrs(units = 'K')
the_ds['rt_cloud'] = the_ds['rt_cloud'].assign_attrs(units = 'kg/kg')
tephigram_skew = 35
the_ds.attrs = {'history': ' written by entraining_plume.ipynb',
                'entrainment_rate':entrain_rate,
                'entrainment_units':'s^{-1}',
                'tephigram_skew' :tephigram_skew,
                'sounding_dir' : str(sounding_dir),
                'sounding_time':the_time,
                'station':station}
```

```{code-cell} ipython3
the_ds
```

## Add the environment variables

```{code-cell} ipython3
env_height = xr.DataArray(data = sounding['hght'], dims={'envlevs':len(sounding)})
env_height = env_height.assign_attrs(units = 'm')
the_ds["env_height"] = env_height

env_thetae = xr.DataArray(data = sounding['thte'], dims={'envlevs':len(sounding)})
env_thetae = env_thetae.assign_attrs(units = 'K')
the_ds["env_thetae"] = env_thetae

env_temp = xr.DataArray(data = sounding['temp'], dims={'envlevs':len(sounding)})
env_temp = env_temp.assign_attrs(units = 'degC')
the_ds["env_temp"] = env_temp

env_dewpoint = xr.DataArray(data = sounding['dwpt'], dims={'envlevs':len(sounding)})
env_dewpoint = env_dewpoint.assign_attrs(units = 'degC')
the_ds["env_dewpoint"] = env_dewpoint
```

```{code-cell} ipython3
the_ds
```

## write the dataset to netcdf

```{code-cell} ipython3
the_ds.to_netcdf(output_filename)
```

## Plot the result on a tephigram

```{code-cell} ipython3
from a405.thermo.constants import constants as c
from a405.thermo.thermlib import convertSkewToTemp, convertTempToSkew
from a405.skewT.fullskew import makeSkewWet,find_corners,make_default_labels
```

```{code-cell} ipython3

```
