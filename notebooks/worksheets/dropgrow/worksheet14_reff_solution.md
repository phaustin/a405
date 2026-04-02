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

(worksheet14_reff_solution)=
# Worksheet 14 solution: liquid water and effective radius

- notebook download link:  [worksheet14_reff_solution.ipynb](https://drive.google.com/file/d/1yZ5JPgDTnd1tNj18Dt_9sSe8-FknoOkK/view?usp=sharing)

This notebook demonstates how to use [pandas.DataFrame.apply](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.apply.html) to add new
columns to a dataframe by iterating over dataframe rows.  It asks you to add columns for the liquid water content and the effective radius, segregated
by updraft vs. downdraft

```{code-cell} ipython3
from pathlib import Path
import numpy as np
from a405.dropgrow.aerolib import lognormal,create_koehler
from a405.utils.helper_funs import make_tuple, find_centers
from a405.thermo.thermlib import find_esat, find_rsat, find_thetaet, tinvert_thetae, find_Td
from a405.thermo.rootfinder import find_interval, fzero
from a405.dropgrow.drop_grow import find_diff
from a405.thermo.constants import constants as c
from scipy.integrate import odeint
import pandas as pd
from matplotlib import pyplot as plt
```

## Library functions from the week 10 koehler worksheet

New functions from [the a405 libraries](https://phaustin.github.io/a405_lib/full_listing.html)

- a405.utils.helper_funs.make_tuple
- a405.utils.helper_funs.find_centers
- a405.dropgrow.aerolib.lognormal
- a405.dropgrow.aerolib.create_koehler
- a405.dropgrow.drop_grow.find_diff

+++

## new functions for droplet growth

Calculate the liquid water content, supersaturation and the derivatives for the drop sizes

+++

#### Define the wvel function

Make it go through 0 $\rightarrow$ $2\pi$ over nsteps

```{code-cell} ipython3
def find_wvel(the_time, period):
    wvel = 0.5*np.sin(the_time*period)
    return wvel

nsteps = 200
period = 2*np.pi/nsteps
times = np.linspace(0,200,nsteps)
wvel_out = [find_wvel(the_time, period) for the_time in times]
plt.plot(times,wvel_out);
```

(rlcalc)=
#### liquid water mixing ratio $r_l$

```{code-cell} ipython3
def rlcalc(var_vec,cloud_tup,numrads):
    """
    calculate the liquid water for the distribution

    the last  3 variables in var_vec are temperature, pressure and height (mks)

    Parameters
    ----------

    var_vec: vector(float)
           vector of values to be integrated
    cloud_top: namedtuple
           tuple of necessary coefficients
    numrads: int
           number of radii in var_vec
    """
    rl=cloud_tup.ndist*(var_vec[:numrads]**3.)
    rl=np.sum(rl)
    rl=rl*4./3.*np.pi*c.rhol
    return rl
```

#### Parcel supersaturation

```{code-cell} ipython3
def Scalc(var_vec,cloud_tup, numrads):
    """
    calculate the environmental saturation using conservation
    of total water mixing ratio cloud_top.rt and the current
    value of the liquid water mixing ratio rl

    That is diagnosing the vapor mixing ratio as the difference
    between the total water in the parcel and the liquid water
    in the droplets

    the last  3 variables in var_vec are temperature, pressure and height (mks)

    Parameters
    ----------

    var_vec: vector(float)
        vector of values to be integrated

    cloud_top: namedtuple
           tuple of necessary coefficients
    numrads: int
           number of radii in var_vec

    Returns
    -------

    Sout: float
       environmental saturation

    """
    temp,press,height = var_vec[numrads:]
    rl=rlcalc(var_vec,cloud_tup,numrads)
    rv=cloud_tup.rt - rl
    e=rv*press/(c.eps + rv)
    Sout=e/find_esat(temp)
    return Sout
```

#### derivative of $r_l$

```{code-cell} ipython3
def rlderiv(var_vec,deriv_vec,cloud_tup,numrads):
    """
    calculate the time derivative of the liquid water content
    using drop_grow.pdf eqn 21b
    the last  3 variables in var_vec are temperature, pressure and height (mks)
    
    Parameters
    ----------

    var_vec: vector(float)
        vector of values to be integrated

    deriv_vec: vector(float)
         derivatives of each of var_vec members

    cloud_tup: namedtuple
           tuple of input coefficients

    nrads:  int
        number of droplet radii variables
        
    Returns
    -------

    drldt: float
         rate of change of rl
    """
    #
    # the slice [:numrads] gives only the droplet radii
    #
    rlderiv=(var_vec[:numrads])**2.
    rlderiv=cloud_tup.ndist*rlderiv
    rlderiv=rlderiv*deriv_vec[:numrads]
    drldt = np.sum(rlderiv)*4.*np.pi*c.rhol
    return drldt
```

### Calculate the other derivatives

```{code-cell} ipython3
def find_derivs(var_vec,the_time,cloud_tup,numrads):
    """
    calcuate derivatives of var_vec 

    a

    Parameters
    ----------

    var_vec: vector(float)
        vector of values to be integrated

    the_time: float
       timestep 

    cloud_tup: namedtuple
           tuple of necessary coefficients
           
    nrads:  int
        number of droplet radii variables
    

    Returns
    -------

    deriv_vec: vector(float)
         derivatives of each of var_vec
    
    """
    # def find_wvel(the_time,period):
    #     wvel = 0.5*np.sin(the_time*period)
    #     return wvel

    period = cloud_tup.period
    the_wvel = find_wvel(the_time,period)
    #print('inside: ',var_vec)
    #print(f"{the_time=}")
    temp,press,height = var_vec[numrads:]
    dry_radius = cloud_tup.dry_radius
    rho=press/(c.Rd*temp)
    #
    # find the evironmental S by water balance
    #
    S=Scalc(var_vec,cloud_tup,numrads)
    deriv_vec=np.zeros_like(var_vec)
    #dropgrow notes equation 18 
    for i in range(numrads):
        m=cloud_tup.masses[i]
        if var_vec[i] < dry_radius[i]:
            var_vec[i] = dry_radius[i]
        Seq=cloud_tup.koehler_fun(var_vec[i],m)  
        rhovr=(Seq*find_esat(temp))/(c.Rv*temp)
        rhovinf=S*find_esat(temp)/(c.Rv*temp)
        #day 25 drop_grow.pdf eqn. 18
        deriv_vec[i]=(c.D/(var_vec[i]*c.rhol))*(rhovinf - rhovr)
    #
    # moist adiabat
    #
    deriv_vec[-3]=find_lv(temp)/c.cpd*rlderiv(var_vec,deriv_vec,cloud_tup,numrads) - c.g0/c.cpd*the_wvel
    #
    # hydrostatic balance  dp/dt = -rho g dz/dt
    #
    deriv_vec[-2]= -1.*rho*c.g0*the_wvel
    #
    # how far up have we traveled?
    #
    deriv_vec[-1] = the_wvel
    return deriv_vec
```

## Get the aerosol and initial conditions

```{code-cell} ipython3
---
jupyter:
  outputs_hidden: false
---
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

initial_conditions = {
        "Tinit": 280.0, #K
        "Zinit": 1000.0, #m
        "Pinit": 90000.0, #Pa
        "Sinit": 0.995,
        "wvel": 0.5  #m/s
}

integration = {
    "tstart": 0,
    "tend": 200,
    "nsteps": 200,
    "numrads":30
}

aero=make_tuple(aerosol_specs)
parcel=make_tuple(initial_conditions)
integration = make_tuple(integration)
koehler_fun = create_koehler(aero,parcel)
```

### initialize the lognormal mass and number distributions for 30 bins

```{code-cell} ipython3
#
#set the edges of the mass bins
#31 edges means we have 30 droplet bins
#
mass_vals = np.linspace(-20,-16,integration.numrads+1) 
mass_vals = 10**mass_vals  #aerosol mass in kg
mu=aero.themean
sigma = aero.sd
totmass = aero.totmass
mdist = totmass*lognormal(mass_vals,np.log(mu),np.log(sigma))
mdist = find_centers(mdist)*np.diff(mass_vals)  #kg/m^3 of aerosol in each binlou
center_mass = find_centers(mass_vals)
ndist = mdist/center_mass  #number/m^3 of aerosol in each bin
#save these in an ordered dictionary to pass to functions
cloud_vars = dict()
cloud_vars['mdist'] = mdist
cloud_vars['ndist'] = ndist
cloud_vars['center_mass'] = center_mass
cloud_vars['koehler_fun'] = koehler_fun
```

```{code-cell} ipython3
def find_lv(temp):
    """
    Calculates the temperature dependent
    enthalpy of evaporation

    Parameters
    ----------

    temp : float or array_like
           Temperature of parcel (K).

    Returns
    -------

    lv : float or list
         entalpy of evaporation


    Examples
    --------

    >>> find_lv(293.)
    2547052.0

    References
    ----------

    Day 13 moist static energy notes
    """
    lv = c.lv0 - (c.cpv - c.cl) * (temp - c.Tc)
    return lv
```

### find the equilibrium radius for each bin at saturation Sinit

```{code-cell} ipython3
---
jupyter:
  outputs_hidden: false
---
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

### now add the intial conditions to the cloud_vars dictionary and make it a namedtuple

the vector var_vec holds 30 droplet radii plus three extra variables at the
end of the vector: the temperature, pressure and height.

```{code-cell} ipython3

the_times = np.linspace(integration.tstart,
                        integration.tend,
                        integration.nsteps)


cloud_vars['initial_radiius'] = initial_radius
cloud_vars['dry_radius'] = dry_radius
cloud_vars['masses'] = center_mass
#
# period for up/down oscillation
#
cloud_vars['period'] = 2*np.pi/integration.nsteps

var_vec = np.empty(integration.numrads + 3)
for i in range(integration.numrads):
    var_vec[i] = initial_radius[i]

#
# temp, press and height go at the end of the vector
#
var_vec[-3] = parcel.Tinit
var_vec[-2] = parcel.Pinit
var_vec[-1] = parcel.Zinit

cloud_tup = make_tuple(cloud_vars)
#calculate the total water (kg/kg)
wl=rlcalc(var_vec,cloud_tup,integration.numrads);
e=parcel.Sinit*find_esat(parcel.Tinit);
wv=c.eps*e/(parcel.Pinit - e)
#save total water
cloud_vars['rt'] = wv + wl
cloud_vars['wvel'] = parcel.wvel
#
# pass this to the find_derivs function
#
cloud_tup= make_tuple(cloud_vars)
cloud_tup.period
```

### use odeint to integrate the variable in var_vec from tinit to tfin with outputs every dt seconds

```{code-cell} ipython3
var_out = []
time_out =[]


sol = odeint(find_derivs,var_vec, the_times, args=(cloud_tup,integration.numrads))
```

### create a dataframe with 33 columns to hold the data

```{code-cell} ipython3
colnames = ["r{}".format(item) for item in range(30)]
colnames.extend(['temp','press','z'])
df_output = pd.DataFrame.from_records(sol,columns = colnames)
```

```{code-cell} ipython3
---
jupyter:
  outputs_hidden: false
---
plt.style.use('ggplot')
plt.close('all')
fig, ax = plt.subplots(1,1,figsize=[10,8])
for i in colnames[:-4]:
    ax.plot(df_output[i]*1.e6,df_output['z'],label=i)
out=ax.set(ylim=[1000,1040],xlim=[0,6],
       xlabel='radii (microns)',ylabel='height (m)',
              title='radii vs. height in a {} m/s updraft'.format(cloud_tup.wvel))
```

## Plot the parcel supersaturation

```{code-cell} ipython3
---
jupyter:
  outputs_hidden: false
---
Svals = []
for index,row in df_output.iterrows():
    out_vec = row.values
    Svals.append(Scalc(out_vec,cloud_tup,integration.numrads))
fig,ax = plt.subplots(1,1,figsize=[10,8])
ax.plot(Svals,df_output['z'])
out=ax.set(ylim=[1000,1050],title='Saturation in a {} m/s updraft'.format(cloud_tup.wvel))
```

## Worksheet problem 1 -- calculating $r_l$

+++

Use [pandas.DataFrame.apply](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.apply.html) to add two new columns to the dataframe

+++

### simple example

1. we want to apply a function along columns (axis=1) rather than rows (axis=0)
2. we want to return the result as a new column in the dataframe, rather than a pandas Series

Downdraft drops: -- write a function that takes a row of the dataframe and any other information needed.
Inside that function, add a new member to the row and return the modified row.

Here's an example that returns the number of droplets with radii $>$ 1 $\mu m$ in descending air

+++

#### Add the time and wvel vectors as new columns

```{code-cell} ipython3
df_output['time'] = the_times
df_output['wvel'] = find_wvel(the_times,cloud_tup.period)
```

```{code-cell} ipython3
plt.plot(df_output['wvel']);
```

```{code-cell} ipython3
def test_fun(row,ndist,numrads):
    """
    sum those members of ndist corresponding to droplet sizes
    larger than 1 micron

    the data frame row contains 30 drop radii followed by temperature, pressure and vertical velocity

    Parameters
    ----------

    row: pd.Series
         Dataframe passed by pandas.DataFrame.apply
    ndist: ndarray -- number of aerosols in each row position
    numrads: int
           number of radii in row
    """
    #print(f"{row.name=}, {row=}")
    #print(f"{row=}")
    full_vec = row.to_numpy()
    drop_radii = full_vec[:numrads]
    hit = drop_radii > 1.e-6
    big_drops = ndist[hit]
    num_sum = np.sum(big_drops)
    row["num_bigdrops"] = num_sum
    return row
```

### apply test_fun to df_output

```{code-cell} ipython3
new_df =df_output.apply(test_fun,axis=1,args = (ndist,integration.numrads) )
new_df.columns
```

#### Plot updraft and downdraft separately

```{code-cell} ipython3
df_updraft = new_df.loc[new_df['wvel'] >= 0.01]
df_downdraft = new_df.loc[new_df['wvel'] < -0.01]
```

```{code-cell} ipython3
df_updraft['time']
```

```{code-cell} ipython3
fig, ax = plt.subplots(1,1)
ax.plot('num_bigdrops','z','b-',data=df_updraft,label='updraft')
ax.plot('num_bigdrops','z','c-',data=df_downdraft,label='downdraft')
fig.legend();
```

## Worksheet problem 2: solution

now add the $r_l$ and $r_{eff}$ column 

1) use pandas apply with the {ref}`rlcalc` function above to add a new dataframe column for total liquid water.  Plot this for the updraft and downdraft
2) write a new function that caculates $r_{eff}$, add the column to the dataframe using apply, and plot for updraft and downdraft

+++

### two new functions

```{code-cell} ipython3
def rl_fun(row,ndist,numrads):
    rl=ndist*row[:numrads]**3.
    rl=np.sum(rl)
    rl=rl*4./3.*np.pi*c.rhol
    row['rl'] = rl
    return row

def reff_fun(row,ndist,numrads):
    numerator = np.sum(ndist*row[:numrads]**3.)
    denominator = np.sum(ndist*row[:numrads]**2.)
    reff = numerator/denominator
    row['reff'] = reff
    return row
```

```{code-cell} ipython3
ws_df = new_df.apply(rl_fun,axis=1,args = (ndist,integration.numrads) )
ws_df = ws_df.apply(reff_fun, axis=1,args = (ndist,integration.numrads) )
ws_df.columns
```

```{code-cell} ipython3
ws_updraft = ws_df.loc[new_df['wvel'] >= 0.01]
ws_downdraft = ws_df.loc[new_df['wvel'] < -0.01]
```

### plots

```{code-cell} ipython3
fig, ax = plt.subplots(1,1)
ax.plot('rl','z','b-',data=ws_updraft,label='updraft')
ax.plot('rl','z','c-',data=ws_downdraft,label='downdraft')
ax.set_ylabel('height (m)')
ax.set_xlabel('$r_l$ (kg/kg)')
ax.set_title('liquid water mixing ratio vs. height (kg/kg)')
ax.legend(loc='best');
```

```{code-cell} ipython3
fig, ax = plt.subplots(1,1)
ax.plot('reff','z','b-',data=ws_updraft,label='updraft')
ax.plot('reff','z','c-',data=ws_downdraft,label='downdraft')
ax.set_ylabel('height (m)')
ax.set_xlabel('$r_{eff}$ (m)')
ax.set_title('effective ratius (m)')
ax.legend(loc='best');
```

```{code-cell} ipython3

```
