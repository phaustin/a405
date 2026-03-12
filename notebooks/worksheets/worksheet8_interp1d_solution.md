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

(worksheet8_interp1d_solution)=
# Worksheet 8: Interp1d Solution

This worksheet is part 1 of Assignment 6 and is due with the rest of the assignment on March 12.

Download [worksheet8_interp1d.ipynb](https://drive.google.com/file/d/1KzMEkSRzKD61yhLW01-hBK8DUjkCEAIR/view?usp=sharing)

Modify a copy of the [entraining_plume.ipynb](https://drive.google.com/file/d/1rLSqA6c6DqP_DVNW2ie2HfERtL4_t002/view?usp=sharing) notebook to add a new function called  `calc_interp`

that contains the following lines from [integ_entrain](https://phaustin.github.io/a405/notebooks/worksheets/entraining_plume.html#define-the-integrator-function)

+++

    press = df_sounding['pres'].values
    height = df_sounding['hght'].values
    temp = df_sounding['temp'].values
    dewpoint = df_sounding['dwpt'].values
    #
    # the nudge function moves any identical heights slightly up or down
    # to avoid breaking the interpolation
    
    envHeight= nudge(height)
    
    interpTenv = interp1d(envHeight,temp)
    interpTdEnv = interp1d(envHeight,dewpoint)
    interpPress = interp1d(envHeight,press)

+++

The signature should look like this:

     def calc_interp(df_sounding):

and it should return the three interperlators interpTenv, interpTdEnv, interpPress.

+++

Check the accuracy of the interpolaters by adding the interpolated soundings onto the environmental sounding plot as a series of large dots.  Are the interpolaters accurately mapping the jagged changes in the temperature and depoint with height?

+++

## Solution

Since I have the sounding data in the netcdf file written by the entraining_plume notebook, I'll use that below.

```{code-cell} ipython3

```

```{code-cell} ipython3
from pathlib import Path
import xarray as xr
from a405.skewT.nudge import nudge
from scipy.interpolate import interp1d
import numpy as np

def calc_interp(entrain_ds):
    press = entrain_ds['env_press'].values
    envHeight = entrain_ds['env_height'].values
    temp = entrain_ds['env_temp'].values
    dewpoint = entrain_ds['env_dewpoint'].values
    emvHeight = nudge(envHeight)
    interpTenv = interp1d(envHeight,temp)
    interpTdEnv = interp1d(envHeight,dewpoint)
    interpPress = interp1d(envHeight,press)
    return interpTenv,interpTdEnv,interpPress
    
```

```{code-cell} ipython3
in_file = Path() / "little_rock_july_9.nc"
if not in_file.exists():
    raise Exception(f"missing input file {in_file.absolute()}")
entrain_ds = xr.open_dataset(in_file)
entrain_ds
```

```{code-cell} ipython3
interpTenv,interpTdEnv,interpPress = calc_interp(entrain_ds)
```

```{code-cell} ipython3
even_height = np.arange(200,9500,10)
even_press=interpPress(even_height)
even_temp = interpTenv(even_height)
even_dew = interpTdEnv(even_height)
```

```{code-cell} ipython3
press = entrain_ds['env_press'].values
temp = entrain_ds['env_temp'].values
dewpoint = entrain_ds['env_dewpoint'].values
tephigram_skew = entrain_ds.attrs['tephigram_skew']
```

```{code-cell} ipython3
from a405.thermo.constants import constants as c
from a405.thermo.thermlib import convertSkewToTemp, convertTempToSkew
from a405.skewT.fullskew import makeSkewWet,find_corners,make_default_labels
from matplotlib import pyplot as plt

fig,ax =plt.subplots(1,1,figsize=(10,10))
corners = [0, 35]
ax, skew = makeSkewWet(ax, corners=corners, skew=tephigram_skew)
xcorners=find_corners(corners,skew=skew)
ax.set(xlim=xcorners,ylim=[1000,200])
ax.set_title("entrainment example")

xcoord_envtemp = convertTempToSkew(temp,press,tephigram_skew)
xcoord_envdew = convertTempToSkew(dewpoint,press,tephigram_skew)
xcoord_even_temp = convertTempToSkew(even_temp,even_press,tephigram_skew)
xcoord_even_dew = convertTempToSkew(even_dew,even_press,tephigram_skew)
ax.plot(xcoord_envtemp,press,'g-',linewidth=3)
ax.plot(xcoord_envdew,press,'b-',linewidth=3)
ax.plot(xcoord_even_temp,even_press,'ro')
ax.plot(xcoord_even_dew,even_press,'co')
```
