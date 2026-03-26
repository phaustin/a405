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
    D=np.arange(0,8,0.001)
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

### Answer 1: calculus

Given a probability density $p(x)$ the definition of the mean is just the weighted average ("first moment") of the distribution:

$$
\left < x \right > = \frac{\int_0^\infty x\,p(x) dx}{\int_0^\infty p(x) dx}
$$

So for the Marshall Palmer distribution that becomes:

$$
\left < D \right > = \frac{\int_0^\infty D\,N_0 \exp(-\Lambda D) dD}{\int_0^\infty N_0 \exp(-\Lambda D) dD}
$$

The denominator is easy:

$$
\int_0^\infty N_0 \exp(-\Lambda D) dD = -\frac{N_0}{\Lambda} \exp(-\Lambda D) \bigg \rvert_0^\infty = \frac{N_0}{\Lambda}
$$(denominator)

For the numerator, the secret is the first line in [this list of integrals](https://en.wikipedia.org/wiki/List_of_integrals_of_exponential_functions):

$$
\int x e^{c x} d x=e^{c x}\left(\frac{c x-1}{c^2}\right)
$$
In our case $c = -\Lambda$ and $x=D$ so:

$$
N_0 \int_0^\infty D e^{-\Lambda D} d D= N_0 e^{-\Lambda D}\left(\frac{-\Lambda D -1}{\Lambda^2}\right) \bigg \rvert_0^\infty = \frac{N_0}{\Lambda^2}
$$(numerator)

Divide {eq}`numerator` by {eq}`denominator` and get:

$$
\left < D \right > = \frac{1}{\Lambda}
$$

+++

### Answer 2: Python

```{code-cell} ipython3
R = 15 #mm/hr
theLambda=41*R**(-0.21)
print(f'mean diameter = {1./theLambda:6.3g} cm')
diam,ndist = marshallpalmer(R)  #cm, m^{-3} mm^{-1}
binwidth = np.diff(diam)[0]*10  #bin width in mm
#
# the python version of equation 2
#
approx_diam = np.sum(diam*ndist*binwidth)/np.sum(ndist*binwidth)
print(f'approx diameter = {approx_diam:6.3g} cm')
```

(mp_ques2)=
## Question 2: Rain rate (part of assignment 9)


- Calculate the precipitation flux (mm/hr) by integrating the total volumen in the Marshall Palmer size distribution and with the fall speed of Villermaux and Bossa (2009): $w = - \sqrt{\rho_l/\rho_{air} * g *D}$ where
D is the drop diameter and $\rho_l,\rho_{air}$ are the liquid and air densities.  Show
that you get about R=15 mm/hour back from this calculation if you use the $\Lambda$ that's appropriate for  R=15 mm/hour

```{code-cell} ipython3
# your code here
```
