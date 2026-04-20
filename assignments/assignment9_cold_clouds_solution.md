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

(assignment9_cold_solution)=
# Assignment 9 solution -- cold clouds

+++

Questions 1: Hand in you answer to {ref}`mp_ques2` from the Marshall-Palmer worksheet.


Question 2: Lohmann problem 4 page 250:

  - Mixed-phase clouds contain ice crystals as well as liquid droplets. Consider such a cloud at a temperature $T=-4^{\circ} \mathrm{C}$, pressure $p=800 \mathrm{hPa}$ and a humidity which corresponds to supersaturations with respect to ice of $5 \%$ and with respect to water of $1 \%$. In the cloud, an ice crystal and a droplet both grow by diffusion, each starting from mass $m_0=10^{-8} \mathrm{~g}$. The ice crystal is a thin hexagonal plate, so that its capacitance can be approximated by $C=2 r_i / \pi$, where $r_i$ is the radius of the ice crystal.

  (a) Determine the times it takes for the droplet and the ice crystal to grow to a total mass $m_1=1.1 \times 10^{-8} \mathrm{~g}$. You can neglect solution and curvature effects when calculating the droplet growth. For the ice crystal, you can assume that its mass $m_i$ and diameter $d_i$ are related by $m_i=\alpha d_i^3$, with $\alpha=1.9 \times 10^{-2} \mathrm{~g} \mathrm{~cm}^{-3}$.

  (b) Which of the two cloud particles grows faster? Explain the main reason for the difference in growth speed.

  (c) How would the situation in the mixed-phase cloud change for supersaturation with respect to ice but subsaturation with respect to water? Explain qualitatively in a few sentences.

+++

## Two equations

+++

Droplet growth:

Use WH 6.21:

$$
r\frac{d r}{d t} \simeq G_l S
$$
where $G_l \simeq 1 \times 10^{-10}\ m^2/s$ as we saw in {ref}`assign8_G_l`

For a droplet of mass $m$ this means that 

$$
\frac{d m}{d t}=4 \pi \rho_l r^2 \frac{d r}{d t} = 4\pi \rho_l r G_l S 
$$

Crystal growth:

$$

$$

```{code-cell} ipython3
import numpy as np
from a405.thermo.constants import constants as c
m0 = 1.1e-11 #kg
r0 = (m0/(4./3.*np.pi*c.rhol))**0.3333
r0
S = 0.01
G_l = 1.e-10
dmdt = 4*np.pi*c.rhol*r0*G_l*S
dmdt
```

Wallace and Hobbs 6.37

$$\frac{d M}{d t}=\frac{C}{\varepsilon_0} G_i S_i
$$

```{code-cell} ipython3

```
