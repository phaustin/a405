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

(assign6_cape_solution)=
# Assignment 6 CAPE solution

## Notebook question

### Original problem statement

- Hand in {ref}`worksheet8_interp1d`
- Add a plot of vertical velocity vs. CAPE, borrowing the code from {ref}`cape_part1`.  Does the plot agree with Thompkins equation 3.4?
  $$
  w_{max} = \sqrt{2 \,\text{CAPE}}
  $$

+++

### Answer

+++

The solution below reads from the netcdf file produced by {ref}`entrain_cloud` and calculates and plots the CAPE as a function of height, using Thompkins equation 3.2.  As Thompkins says, the CAPE is just the integrated work between an initial height (we can start at $z^\prime = 0$) and the observation height $z$, if all that buoyant accleration is converted to kinetic energy with $E_k = \frac{1}{2} m w^2$ and for m=1 kg  then the velocity $w$ at height $z$ is:

$$
B &= g \left(\frac{T_{v\,\text{cloud}}-T_{v\,\text{env}}}{T_{v\,\text{env}}}\right) 
$$
$$
\text{CAPE}&=\int_0^{z} B dz^ \prime
$$

$$
\frac{1}{2} w^2 = \text{CAPE}
$$

$$
w(z) = \sqrt{2 \, \text{CAPE}}
$$

+++



```{code-cell} ipython3
import numpy as np

from a405.thermo.constants import constants as c
from a405.skewT.nudge import nudge
```

### From Assignb 6 part 1: use `calc_interp` to construct interpolators

```{code-cell} ipython3
from pathlib import Path
import xarray as xr
from matplotlib import pyplot as plt
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

### Read in the entrainment run

Read in the file produced by {ref}`entrain_cloud`

```{code-cell} ipython3
in_file = Path() / "little_rock_july_9.nc"
if not in_file.exists():
    raise Exception(f"missing input file {in_file.absolute()}")
entrain_ds = xr.open_dataset(in_file)
entrain_ds
```

### copy get_buoy from the cape notebook

```{code-cell} ipython3
from numpy.typing import ArrayLike as array
import a405.thermo.thermlib as tl
from a405.thermo.constants import constants as c


def get_buoy(T_env: array, Td_env: array, T_cloud: array, press: array) -> array:
    """
    Given an atmospheric sounding (Temperature and dew point profiles),
    calculate the buoyant force per unit mass (m/s2) of cloud air previously
    calculated from the entraining parcel notebook

    Parameters
    ----------
    Tenv, Tdenv, Tcloud and press all need to be the same length

    Tenv:  environmental emperature profile (K)
    Tdenv:  environmental dew point temperature profile (K)
    pres: Pressures at which to calculate B (Pa)

    Returns
    -------

    B: Buoyancy per unit mass (m/s2)

    """
    rv_cloud = tl.find_rsat(T_cloud,press)
    Tv_cloud = tl.find_Tv(T_cloud, rv_cloud)
   
    rv_env = tl.find_rsat(Td_env, press)
    Tv_env = tl.find_Tv(T_env, rv_env)
    
    # Thompkins 1.61
    buoy = c.g0 * (Tv_cloud - Tv_env) / Tv_env
    return buoy
```

### interpolate to the cloud pressure levels

```{code-cell} ipython3
interpTenv,interpTdEnv,interpPress = calc_interp(entrain_ds)
even_height=entrain_ds['cloud_height'].values
even_press=interpPress(even_height)*c.hPa2pa
even_temp = interpTenv(even_height) + c.Tc
even_dew = interpTdEnv(even_height) + c.Tc
wvel = entrain_ds['wvel'].values
cloud_temp = entrain_ds['cloud_temp'].values

buoy = get_buoy(even_temp, even_dew, cloud_temp, even_press)
```

### plot temperatures, buoyancy and vertical velocity

```{code-cell} ipython3
fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(10,6))
ax1.plot(cloud_temp,even_press,label="cloud")
ax1.plot(even_temp,even_press,label="env")
ax1.set_xlabel("temperature (K)")
ax1.set_ylabel("pressure (Pa)")
ax1.invert_yaxis()
ax1.grid(True)
ax1.set_title("environmental and cloud temperture")
ax1.legend()
ax2.invert_yaxis()
ax2.plot(buoy,even_press)
ax2.axes.get_yaxis().set_visible(False)
ax2.grid(True)
ax2.set_xlabel("B (m/s^2)")
ax2.set_title("Buoyancy")
ax3.plot(wvel,even_press)
ax3.axes.get_yaxis().set_visible(False)
ax3.set_xlabel("$w_{vel}$ (m/s)")
ax3.set_title("vertical velocity")
```

```{code-cell} ipython3
from matplotlib import pyplot as plt
```

### Calculate the cape and $\sqrt{2*\text{CAPE}}$

```{code-cell} ipython3
print(f"initial wvel at cloud base is {wvel[0]:.1f} m/s")
```

Note the negative buoyancy at cloud base that quickly slows the parcel down to close to 0 m/s, making the agreement between wvel and CAPE much better than it would be otherwise.

```{code-cell} ipython3
mid_buoy = (buoy[:-1] + buoy[1:])/2.
dz = np.diff(even_height)
cumulative_cape = np.cumsum(mid_buoy*dz)
wvel_estimate=(2*cumulative_cape)**0.5 
```

### Agreement is pretty good

```{code-cell} ipython3
fig, ax = plt.subplots(1,1)
ax.plot(wvel[:-1],even_height[:-1],label='wvel')
ax.plot(wvel_estimate,even_height[:-1],label='cape')
ax.set_xlabel('velocity (m/s)')
ax.set_ylabel('height (m)')
ax.legend();
```

```{code-cell} ipython3

```
