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

# assignment 8  solution

+++

### Problem 4 -- Wallace and Hobbs problem 6.24 -- falling precip

If a raindrop has a radius of 1 mm at cloud base, which is located 5 km above the ground, what will be its radius at the ground and how long will it take to reach the ground if the relative humidity between cloud base and ground is constant at 60%? [Hint: Use (6.21) and the relationship between v and r given in Exercise 6.23. If r is in micrometers, the value of Gl in (6.21) is 100 for cloud droplets, but for the large drop sizes considered in this problem the value of Gl should be taken as 700 to allow for ventilation effects.]

+++

From WH 6.23: 

$$
v=6 \times 10^3 r
$$

From wH Eq. 6.21:

$$
r \frac{d r}{d t}=G_{\ell} S
$$

```{code-cell} ipython3

```

```{code-cell} ipython3
from scipy.integrate import solve_ivp

def hit_ground(t, the_vars):
    """Event: when height (the_vars[0]) is zero."""
    return the_vars[0]

def find_derivs(the_time,the_vars):
    dx_dt = the_vars[1]
    dy_dt = -the_vars[0]
    return [dx_dt,dy_dt]

# Set hit_ground to be a terminal event (stops integration)
# and only trigger when the value is going from positive to negative (direction=-1)
# to avoid triggering at the start (t=0).
hit_ground.terminal = True
hit_ground.direction = -1
tspan=[0,6]

sol = solve_ivp(
    find_derivs,
    tspan,
    init_0,
    events=[hit_ground],
    dense_output=True # Enable interpolation for plotting smooth curve
)
```
