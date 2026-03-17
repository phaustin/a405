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

# Assignment 6 paper/pencil problems

+++

# Kohler problem solutions

## 1. Taylor series

Show using two Taylor series expansions how to get  Lohmann 6.24:

$S= a_w \exp \left [ \frac{2\sigma}{\rho_l R_v T r} \right ] \approx \left ( 1 + \frac{a}{r} - \frac{b}{r^3} \right )$

### 1. Answer

#### First taylor series

$a_w = \left ( \frac{n_w}{n_w + n_s} \right ) = \left ( \frac{1}{1 + n_s/n_w} \right ) = (1 + n_s/n_w)^{-1}$

If $n_s/n_w  = x \ll 1$  then expand $(1 + x)^{-1}$ in a Taylor series about x=0:

$f^\prime (0) = -(1 + 0 )^{-2} = -1$, $f^{\prime \prime} (0) = 2 ( 1 + 0 )^{-3} = 2$

so $(1 + x)^{-1} \approx 1 - x + 2x^2/2 + \ldots$

#### Second taylor series

$\exp \left [ \frac{2\sigma}{\rho_l R_v T r} \right ] = \exp \left ( \frac{a}{r} \right ) = \exp(y)$

if $y \ll 1$ then expand exp in a taylor series about y=0:

$f^\prime (0) =\exp(0) = 1$

$f^{\prime \prime}(0) = 1$

so $\exp(y) \approx 1 + y  + y^2/2 + \ldots$

#### Combining:

$a_w \exp \left [ \frac{2\sigma}{\rho_l R_v T r} \right ] \approx (1 - x + x^2 )(1 + y + y^2/2) = 1 + y - x + y^2/2 + x^2 + \ldots$

and keeping only first order terms:

$S = 1 + y - x = 1 + \frac{a}{r} - n_s/n_w$

and since $n_w \propto r^3$: $S = 1 + y - x = 1 + \frac{a}{r} - \frac{b}{r^3}$

where we've dropped all terms of order $xy$ or higher, i.e. we've dropped

$$
\text{dropped\_term} =\frac{ab}{r^4}
$$

+++

## 2. Term Comparison

For the aerosol defined in the kohler.ipynb notebook ($10^{-19}$ kg, ammonium sulphate), inside a droplet of radius $r=1\ \mu m$
      find the size of the smallest term you've kept (either $\frac{a}{r}$ or $\frac{b}{r^3}$ and compare it to
      the size of the largest term you've dropped.   Repeat this for a droplet of radius   $r=0.1\ \mu m$
### 2. Answer

```{code-cell} ipython3
from a405.thermo.constants import constants as c
import numpy as np

def find_terms(r):
    Tinit=c.Tc + 15 #Temperature K
    Ms=132 #ammonium sulphae kg/Kmole
    Mw=18. #water kg/Kmole
    Sigma=0.075  #N/m
    vanHoff=3. #van hoff for ammonium bisulphate
    #calculate kohler coefficients:
    a=(2.*Sigma)/(c.Rv*Tinit*c.rhol)  #curvature term
    mass_aero = 1.e-19  #kg
    b=(vanHoff*Mw)/((4./3.)*np.pi*c.rhol*Ms)*mass_aero   #solution term
    kv=a/r
    rt=b/(r**3.)
    dropped= kv*rt
    return kv,rt ,dropped

def template(kelvin_term,raoult_term,radius,dropped):
    out=(f'kelvin term = {kelvin_term:8.1g} vs '
       f'raoult term {raoult_term:8.1g} for ' 
       f'radius= {radius*1.e6:5.1f} microns\n'
       f'largest term dropped is: {dropped:8.1e}')
    return out


for radius in [1.e-6,1.e-7,1.e-8]:
    print(f"\nradius = {radius*1.e6:.3f} microns\n")
    kelvin_term,raoult_term, dropped=find_terms(radius)
    print(template(kelvin_term,raoult_term,radius, dropped))
```

Conclusion: don't use the Taylor series expansion for very small droplets

+++

## Surface energy vs. latent heat

Suppose you have $r_l$ =1 g/kg of liquid water spread among N spherical drops of radius 10 $\mu m$.  Compare the surface energy of this
      population (which we've been neglecting) with the enthalpy required to vaporize it:  $l_v r_l$.  Is it negligible in comparison?

```{code-cell} ipython3
sigma = 0.075 # N/m^2
rhol=1000  #kg/m^3
rl=1.e-3 #g/kg
r=1.e-5  #m
N = rl/(4./3.*np.pi*r**3.*rhol)
surface_energy = N*sigma*4*np.pi*r**2.
bulk_enthalpy = c.lv0*rl
print(f'there are {N=:6.1g} drops/kg')
print(f'{surface_energy=:6.3f} J/kg, bulk enthalpy of evaporation {bulk_enthalpy:.1f} J/kg')
error = surface_energy/bulk_enthalpy*100
print(f"we are neglecting {error:.1e}% of the energy")
```
