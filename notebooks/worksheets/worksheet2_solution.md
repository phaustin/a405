---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.18.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

(worksheet2_solution)=
# Worksheet2 -- Thermodynamics worksheet solution

+++

Download link [worksheet2_solution.ipynb](https://drive.google.com/file/d/1xtJfUjoR4FDEch_6zeDEANwa2ZkVvVXl/view?usp=sharing)

## Working with the A405 thermodynamics libraries 

### Installation 

+++

The source code for the libraries is in the [A405 github repository](https://github.com/phaustin/a405_2024/tree/main/src/a405).  To install on your laptop, use this [requirements.txt](https://github.com/phaustin/a405/blob/main/requirements.txt
) with pip  (note only one dash in -r, but two dashes in --update) in the a405 environment.

+++

    pip install -r requirements.txt --upgrade

+++

or just use the line inside the requirements.txt file directly from the prompt

+++

     pip install git+https://github.com/phaustin/a405.git --upgrade

+++

### Documentation

+++

The A405 library modules are documented at [https://phaustin.github.io/a405_lib/full_listing.html](https://phaustin.github.io/a405_lib/full_listing.html)

+++

## Practice problems

+++

1)  Find the potential temperature of air at a pressure of 450 hPa and a temperature of 300 K

+++

### 1. Answer

```{code-cell} ipython3
import numpy as np
import a405
from a405.thermo import thermlib
temp=300
press=45000
answer = a405.thermo.thermlib.find_theta(temp, press, rv=0)
print(f"{answer=:.2f} K")
```

2. Two jars of liquid water are placed in an insulated vacumn chamber with completely reflecting surfaces in the longwave.  They exchange longwave radiation until they come into equilibrium.  Neglect the water vapor in the chamber and the container material, and assume that Jar A contains 1 kg of liquid and Jar B contains 2 kg of liquid.

    a. I the 1 kg jar A starts at 350 K and the 2 kg Jar B starts at 280 K, what is the final equilibrium temperature of each jar?  
    b. What is the total change in entropy for the system composed of the two jars?  (liquid water is incompressible)

+++

### 2a. Answer

```{code-cell} ipython3
from a405.thermo.constants import constants as c
```

```{code-cell} ipython3
massA=1
massB=2
tempA=350
tempB=280
H_A=massA*c.cl*tempA
H_B=massB*c.cl*tempB
H_tot = H_A + H_B
h_new = (H_A + H_B)/(massA + massB)
temp_new = h_new/c.cl
print(f"{temp_new=:0.1f} K")
```

### 2b. Answer

Combine the first and second laws to get the liquid water entropy.  Use the fact that liquid water is incompressible so that the work term $w_t$ in the first law is zero.  We also need to replace the atmospheric heat capacity $c_v$ with the liquid water heat capacity $c_l$

- First law

$$
q_t   = \frac{du}{dt} + w_t = \frac{du}{dt}  + 0  = c_l \frac{dT}{dt} \\
q_t dt = c_l dT
$$

- Second law

$$
d\phi \ge \frac{q_t dt}{T}
$$

- Combined for reversible (quasi-equilibrium) process


$$
d\phi = c_l \frac{dT}{T} = c_l\, d\ln T
$$

- Integrate from $T^\prime=0$ to $T^\prime = T$

$$
\phi = c_l \ln T
$$

```{code-cell} ipython3
massA=1
massB=2
tempA=350
tempB=280
Phi_A=massA*c.cl*np.log(tempA)
Phi_B=massB*c.cl*np.log(tempB)
Phi_start = Phi_A + Phi_B
Phi_new = (massA + massB)*c.cl*np.log(temp_new)
```

```{code-cell} ipython3
print(f"Entropy before = {Phi_start:.0f} J/K")
print(f"Entropy after = {Phi_new:.0f} J/K")
print(f"Total change in Entropy = {Phi_new - Phi_start:.0f} J/K")
```

## Study question

Use the second law of thermodynamics to show that energy must always flow from warmer to cooler objects.

```{code-cell} ipython3

```
