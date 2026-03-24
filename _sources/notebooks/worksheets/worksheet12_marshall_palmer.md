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

(worksheet12_marshall_palmer)=
# Worksheet 12: Marshall-Palmer Distribution

Download link: [worksheet12_marshall_palmer](https://drive.google.com/file/d/1f95GnUpYVnHY7vtQmttrx8t1_Pmh8zNI/view?usp=sharing)

+++

## Code and questions about the Marshall-Palmer distribution

See Thompkins p. 81 and Wallace and Hobbs page 232

```{code-cell} ipython3
import numpy as np
from matplotlib import pyplot as plt
plt.style.use('ggplot')

def marshallpalmer(R):
    """
    marshall palmer size distribution
    given rainrate R in mm/hr, return
    n(D), the number concentration of drops with
    diameter D

    Parameters
    ----------
    R: float
        rainrate (mm/hr)

    Returns
    -------

    d: vector (float)
      drop diameters (cm)

    n: vector (float)
     the number distribution n(d) #m^{-3} mm^{-1}

    """
    D=np.arange(0,8,0.01)
    Dmm=D
    Dcm=D*0.1
    N0=0.08*1.e6*1.e-1 #m**{-3} mm^{-1}
    theLambda=41*R**(-0.21)
    n=N0*np.exp(-theLambda*Dcm)
    return Dcm,n

curve_dict={}
Rvals = [1,5,25]
for R in Rvals:
    diam,ndist = marshallpalmer(R)
    curve_dict[R] = ndist
fig, ax = plt.subplots(1,1,figsize=(10,8))
for R in Rvals:
    ax.semilogy(diam,curve_dict[R],label='{} mm/hr'.format(R))
ax.set_xlabel('Drop diameter (mm)')
ax.set_ylabel('$n(D) m^{-3} mm^{-1}$')
ax.set_title('Marshall Palmer distribution for three rain rates')
out=ax.legend(loc='best')
```

##  Question 1: Mean diameter

- Confirm that the mean diameter of the Marshall Palmer distribution is $1/\Lambda$

+++

## Question 2: Rain rate (part of assignment 9)


- Calculate the precipitation flux (mm/hr) by integrating the total volumen in the Marshall Palmer size distribution and with the fall speed of Villermaux and Bossa (2009): $w = - \sqrt{\rho_l/\rho_{air} * g *D}$ where
D is the drop diameter and $\rho_l,\rho_{air}$ are the liquid and air densities.  Show
that you get about R=15 mm/hour back from this calculation if you use the $\Lambda$ that's appropriate for  R=15 mm/hour

```{code-cell} ipython3
# your code here
```
