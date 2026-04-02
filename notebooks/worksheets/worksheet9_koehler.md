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

(worksheet9_koehler)=
# Worksheet9: Koehler Equilibrium

```{code-cell} ipython3
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from a405.thermo.constants import constants as c
```

This worksheet is designed to get you familiar with the Koehler equation for a specific aerosol.

The download link for [worksheet9_koehler.ipynb](https://drive.google.com/file/d/1fIWW1Zo1-2Qm0ROP6ZqgG6rKcrC_ntuy/view?usp=sharing)

+++

## Question 1

+++

Plot the Kohler curve (Thompkins equation 4.15, and Fig. 4.13) for an aerosol:

$$
S = \frac{e_\chi}{e_{s\,\infty}} \approx 1 + \frac{a}{r} - \frac{b}{r^3}\tag{AT 4.15}
$$
where $e_\chi$ is the equilibrium vapor pressure over a curved surface with water and ions, and $e_{s\, \infty}$ is the equilibrium vapor pressure over pure water with a flat surface (i.e. the Clausius-Clapeyron equation)


$$
a = \frac{2\sigma}{\rho_l R_v T}
$$

$$
b = \frac{imM_w}{4/3\pi \rho_s M_s}
$$

+++

Download this [aerosol_file](https://drive.google.com/file/d/1BFG11eITRO_6qLblE_U1g457RC5JYvde/view?usp=sharing) into the same folder as  your notebook

```{code-cell} ipython3
# load the properties of ammonium sulphate from JSON file
filename =  Path("ammonium_sulphate.json")
with open(filename) as infile:
    ammonium_sulphate = json.load(infile)
print(ammonium_sulphate)
```

```{code-cell} ipython3
def find_S(r, aerosol, T):
    """
    calculates supersaturation S given an aerosol dictionary,
    temperature T, and droplet radius r 
    uses Thompkins 4.15

    Parameters
    ----------
    r: radius (me)
    aerosol: dictionary with keys: Sigma (N/m^2), vanHoff (int), Mw (kg/kmole), Ms (kg/kmole), mass (kg)
    T: temperature (K)

    Returns
    -------
    S: supersaturation (unitless)
    
    """
    # your code here
  
    return S
```

```{code-cell} ipython3
# In this cell create a vector of radii r and use find_S to get the corresponding S values
```

```{code-cell} ipython3
# In this cell make an (x,y) plot of (r,S)
```

## Question 2 -- rootfinding the equilibrium radius

Use our rootfinder to find the equilibrium radius for a haze particle at a relative humidity of 90% and a temperature of 15 deg C

```{code-cell} ipython3
#In this cell write a function find_r that takes the arguments to find_S and a relative humidity RH
# and changes sign when S crosses over that relative humidity 
# then use optimize.brentq to find the radius where S(r) = RH
```

```{code-cell} ipython3

```
