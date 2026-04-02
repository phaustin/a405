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

(worksheet13_odeint)=
# Worksheet 13 odeint practice

Simple integration practice for assignment 8, using odeint as in 
{ref}`worksheet11_dropgrow_solution`

Download [worksheet13_odeint.md](https://drive.google.com/file/d/14Aw1Pr-rP-W4N5S0QsK4xx5tgnZI1jVL/view?usp=sharing)

+++

Problem: integrate the following system of odes:

$$
\frac{dx}{dt} &= y \\
\frac{dy}{dt} &= -x
$$
with intial conditions [x, y] = [0, 1] using first [odeint](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html) and then
[solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)

+++

## Using odeint

Step 1 supply the derivatives, initial condition, and time vector

```{code-cell} ipython3
from scipy.integrate import odeint
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# your code here
#  define find_derivs, int_0 and the_time
#

#
# run the integrator
#
sol = odeint(find_derivs,init_0, the_time)
#
# turn the solution object into a dataframe
#
df_output = pd.DataFrame.from_records(sol,columns = ['xvar','yvar'])
df_output['time'] = the_time
```

### Plot the result

```{code-cell} ipython3
fig, ax = plt.subplots(1,1)
ax.plot('time','xvar',data=df_output)
ax.plot('time','yvar',data=df_output)
ax.legend();
```

## Repeat using [solve_ivp](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)

set up an event to stop when $x< 0$

```{code-cell} ipython3
from scipy.integrate import solve_ivp

#
# define an event called "hit_ground" that will return 0 when variable x =0
#


#
# your code here
# 

# Set hit_ground to be a terminal event (stops integration)
# and only trigger when the value is going from positive to negative (direction=-1)
# to avoid triggering at the start (t=0).
hit_ground.terminal = True
hit_ground.direction = -1

sol = solve_ivp(
    find_derivs,
    tspan,
    init_0,
    events=[hit_ground],
    dense_output=True # Enable interpolation for plotting smooth curve
)
```

```{code-cell} ipython3
dir(sol)
```

```{code-cell} ipython3
sol.y_events
```

## Interpolate to uniform grid

Use the sol object to interpolate onto a uniform set of times

Note that the shape is now [2,numrows], where above sol was [numrows, 2], you you need to 
transpose yplot to create the dataframe

```{code-cell} ipython3
the_times = np.arange(0,sol.t[-1],0.01)
yplot = sol.sol(the_times)
df_output = pd.DataFrame.from_records(yplot.T,columns = ['xvar','yvar'])
df_output['time']=the_times
```

```{code-cell} ipython3
ig, ax = plt.subplots(1,1)
ax.plot('time','xvar',data=df_output)
ax.plot('time','yvar',data=df_output)
ax.legend();
```

```{code-cell} ipython3

```
