---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

## Parcel model with 30 aerosol masses, lognormal distribution

```{code-cell} ipython3
import json
import a405.utils
from pathlib import Path
import numpy as np
from a405.dropgrow.aerolib import lognormal,create_koehler
from a405.utils.helper_funs import make_tuple, find_centers
from a405.thermo.thermlib import find_esat, find_lv
from a405.thermo.rootfinder import find_interval, fzero
from a405.dropgrow.drop_grow import find_diff, rlcalc, find_derivs, Scalc
from a405.thermo.constants import constants as c
from scipy.integrate import odeint
import pandas as pd
from matplotlib import pyplot as plt
import datetime
import importlib_resources as ir 
import pprint
pp = pprint.PrettyPrinter(indent=4)
```



+++

## Read in the json file and set the koehler function for this aerosol

```{code-cell} ipython3
aerosol_specs = {
    "Ms": 114,  #molecular weight of aerosol
    "Mw": 18.0,  #molecular weight of water
    "Sigma": 0.075,  # surface tension N/m^2
    "vanHoff": 2.0,  
    "rhoaero": 1778, #aerosol density, kg/m^3
    "themean": 2e-17,  #mean mass kg
    "sd": 1.7,  #standard deviation (kg)
    "totmass": 1.5e-09  #kg/m^3 
}
```

```{code-cell} ipython3
initial_conditions = {
        "Tinit": 280.0, #K
        "Zinit": 1000.0, #m
        "Pinit": 90000.0, #Pa
        "Sinit": 0.995,
        "wvel": 0.5  #m/s
}
```

```{code-cell} ipython3
aero=make_tuple(aerosol_specs)
parcel=make_tuple(initial_conditions)

koehler_fun = create_koehler(aero,parcel)
    
```

## initialize the lognormal mass and number distributions for 30 bins

```{code-cell} ipython3
#
#set the edges of the mass bins
#31 edges means we have 30 droplet bins
#
numrads = 30
mass_vals = np.linspace(-20,-16,numrads+1) 
mass_vals = 10**mass_vals  #aerosol mass in kg
mu=aero.themean
sigma = aero.sd
totmass = aero.totmass
mdist = totmass*lognormal(mass_vals,np.log(mu),np.log(sigma))
mdist = find_centers(mdist)*np.diff(mass_vals)  #kg/m^3 of aerosol in each bin
center_mass = find_centers(mass_vals)
ndist = mdist/center_mass  #number/m^3 of aerosol in each bin
#save these in an ordered dictionary to pass to functions
cloud_vars = dict()
cloud_vars['mdist'] = mdist
cloud_vars['ndist'] = ndist
cloud_vars['center_mass'] = center_mass
cloud_vars['koehler_fun'] = koehler_fun
```

## find the equilibrium radius for each bin at saturation Sinit

```{code-cell} ipython3
S_target = parcel.Sinit
logr_start = np.log(0.1e-6)

initial_radius = []
dry_radius = []
for mass in center_mass:
    brackets = np.array(find_interval(find_diff,logr_start,S_target,mass,koehler_fun))
    left_bracket, right_bracket = np.exp(brackets)*1.e6  #get brackets in microns for printing
    equil_rad = np.exp(fzero(find_diff,brackets,S_target,mass,koehler_fun))

    initial_radius.append(equil_rad)
    dry_rad = (mass/(4./3.*np.pi*aero.rhoaero))**(1./3.)
    dry_radius.append(dry_rad)

    print('mass = {mass:6.3g} kg'.format_map(locals()))
    print('equlibrium radius at S={} is {:5.3f} microns\n'.format(S_target,equil_rad*1.e6))
```

## now add the intial conditions to the cloud_vars dictionary and make it a namedtuple

the vector var_vec holds 30 droplet radii plus three extra variables at the
end of the vector: the temperature, pressure and height.

```{code-cell} ipython3
cloud_vars['initial_radiius'] = initial_radius
cloud_vars['dry_radius'] = dry_radius
cloud_vars['masses'] = center_mass
numrads = len(initial_radius)
var_vec = np.empty(numrads + 3)
for i in range(numrads):
    var_vec[i] = initial_radius[i]

#
# temp, press and height go at the end of the vector
#
var_vec[-3] = parcel.Tinit
var_vec[-2] = parcel.Pinit
var_vec[-1] = parcel.Zinit

cloud_tup = make_tuple(cloud_vars)
#calculate the total water (kg/kg)
rl=rlcalc(var_vec,cloud_tup);
e=parcel.Sinit*find_esat(parcel.Tinit);
rv=c.eps*e/(parcel.Pinit - e)
#save total water
cloud_vars['rt'] = rv + rl
cloud_vars['wvel'] = parcel.wvel
cloud_vars['wvel'] = 1.5
#
# pass this to the find_derivs function
#
cloud_tup= make_tuple(cloud_vars)
```

## use odeint to integrate the variable in var_vec from tinit to tfin with outputs every dt seconds

```{code-cell} ipython3
var_out = []
time_out =[]

tinit=input_dict['integration']['dt']
dt = input_dict['integration']['dt']
tfin = input_dict['integration']['tend']

t = np.arange(0,tfin,dt)
sol = odeint(find_derivs,var_vec, t, args=(cloud_tup,))
```

## create a dataframe with 33 columns to hold the data

```{code-cell} ipython3
colnames = ["r{}".format(item) for item in range(30)]
colnames.extend(['temp','press','z'])
df_output = pd.DataFrame.from_records(sol,columns = colnames)
```

## store the dataframe in an csv file, including a copy of the input dictionary for future reference

```{code-cell} ipython3
if input_dict['dump_output']:
    outfile_name = f'{input_dict["output_file"]}.csv'
    with open(outfile_name,'w') as store:
       df_output.to_csv(store)

    metadata_name = f'{input_dict["output_file"]}.json'
    date=datetime.datetime.now().strftime('%Y-%M-%d')
    with open(metadata_name,'w') as meta:
        history ="file produced by drop_grow.py on {}".format(date)
        print('history: ',history)
        input_dict['history']=history
        json.dump(input_dict,meta,indent=4)
```

```{code-cell} ipython3
fig, ax = plt.subplots(1,1,figsize=[10,8])
for i in colnames[:-3]:
    ax.plot(df_output[i]*1.e6,df_output['z'],label=i)
out=ax.set(ylim=[1000,1040],xlim=[0,6],
       xlabel='radii (microns)',ylabel='height (m)',
              title='radii vs. height in a {} m/s updraft'.format(cloud_tup.wvel))
```

```{code-cell} ipython3
Svals = []
for index,row in df_output.iterrows():
    var_vec = row.values
    Svals.append(Scalc(var_vec,cloud_tup))
fig,ax = plt.subplots(1,1,figsize=[10,8])
ax.plot(Svals,df_output['z'])
out=ax.set(ylim=[1000,1050],title='Saturation in a {} m/s updraft'.format(cloud_tup.wvel))
```

## Predict $dSS/dt$ at 1010m

Define a function that implements equation 6 from the equilibrium supersaturation notes:

$$
\frac{dr_v}{dt} = (1 + SS)\left[\frac{-\epsilon e_s}{p^2}\left(\frac{-gpV}{R_dT}\right) + \frac{\epsilon}{p}\left(\frac{\epsilon Le_s}{R_dT^2}\right)\frac{dT}{dt}\right] + \frac{\epsilon e_s}{p}\frac{dSS}{dt}\tag{PA 6}
$$

Re-arrange to solve for $dSS/dt$:

$$
\frac{dSS}{dt} = \frac{p}{\epsilon e_s}\left(\frac{dr_v}{dt} - (1 + SS)\left[\frac{-\epsilon e_s}{p^2}\left(\frac{-gpV}{R_dT}\right) + \frac{\epsilon}{p}\left(\frac{\epsilon Le_s}{R_dT^2}\right)\frac{dT}{dt}\right]\right)
$$

To get $dT/dt$, assume a static temperature at each level and constant updraft speed of 1.5m/s, yielding:

$$
\frac{dT}{dt} = \frac{dT}{dz}\frac{dz}{dt} = V\frac{dT}{dz}
$$

Similarly find $dr_v/dt$:

$$
r_v = \frac{\epsilon e_s}{p - e_s}
$$

$$
\frac{dr_v}{dt} = \frac{dr_v}{dz}\frac{dz}{dt} = V\frac{dr_v}{dz}
$$

```{code-cell} ipython3
# pull out variables for PA6
SS = np.array(Svals)
press = df_output.press
z = np.array(df_output.z)
T = df_output.temp
es = find_esat(T)
V = 1.5  # m/s

# pack inputs into a tuple
inputs = (SS, press, z, T, es, V)
```

```{code-cell} ipython3
def find_dSSdt(SS, press, z, T, es, V):
    """
    finds time derivative of SS given an updraft speed and other params
    """
    # get the height deriv of T and interpolate
    dTdz_lyr = np.diff(T) #/ np.diff(z)
    dTdz_plus = np.append(dTdz_lyr, dTdz_lyr[-1])
    dTdz_minus = np.append(dTdz_lyr[0], dTdz_lyr)
    Tgrad = (dTdz_plus + dTdz_minus) / 2
    
    # same for rv
    rv = c.eps * es / (press - es)
    drvdz_lyr = np.diff(rv) #/ np.diff(z)
    drvdz_plus = np.append(drvdz_lyr, drvdz_lyr[-1])
    drvdz_minus = np.append(drvdz_lyr[0], drvdz_lyr)
    rv_grad = (drvdz_plus + drvdz_minus) / 2
    
    # solve for dSS/dt
    A = press / (c.eps * es)
    B = 1 + SS
    C = (c.eps * es / press ** 2) * ((c.g0 * press * V) / (c.Rd * T))
    D = (c.eps / press) * (c.eps * find_lv(T) * es * V * Tgrad) / (c.Rd * T ** 2)
    dSSdt = A * (rv_grad - B * (C + D))
    return dSSdt
```

```{code-cell} ipython3
dSSdt = find_dSSdt(*inputs)

ht = 1010 # m
idx = (np.abs(z - ht)).argmin()
print(f"dSS/dt at {ht} m: {dSSdt[idx]}")
```
