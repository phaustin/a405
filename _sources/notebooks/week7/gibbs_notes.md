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

(gibbs_notes)=
# Gibbs Free Energy notes

## 1. Show that for phase change at constant temperature and pressure, the Gibbs free energy can only decrease or stay the same.

Start with the definition of the specific Gibbs free energy:

$$
g=h-T \phi
$$(eq:g)
where $h$ is the enthalpy, $T$ the temperature and $\phi$ the entropy.

We also have the second law:

$$
d \phi \geq  \frac{qdt}{T}
$$(eq:second)
where $q$ is the heating rate applied for time $dt$.

Take the differential of {eq}`eq:g`

$$
d g= -T d \phi \\
$$(eq:dg)
where we have used the fact that $dh = 0$ at constant T.

Finally, note that the approach to equilibrium is occuring without any heating, so {eq}`eq:second` becomes:

$$
T d\phi \geq 0
$$
which means from  {eq}`eq:dg` that 

$$
dg \leq 0
$$

So the Gibbs free energy can never increase during phase change.

+++

## 2. Show that in phase equilibrium the vapor and liquid water Gibbs free energies are equal

Given a closed system with vapor mass $m_v$ and Gibbs free energy $g_v$ and liquid water mass $m_l$ and vapor mixing ratio $g_l$, the total Gibbs free energy is:

$$
 G=m_v g_v + m_l g_l
$$

Since the system is closed the total mass of water is constant:

$$
m_v &+ m_l = C \\
dm_v &+ dm_l = 0\\
dm_v &= -dm_l
$$(eq:constant)

At constant temperture and pressure, both $g_v$ and $g_l$ are constant so:

$$
dG =  g_v dm_v + g_l dm_l
$$

And using {eq}`eq:constant`


$$
dG =  g_v dm_v + g_l (-dm_v) = (g_v - g_l) dm_v
$$

But we know from question 1. that in equilibrium $dG = 0$ so

$$
dG = 0 = (g_v - g_l) dm_v 
$$

Which reduces to

$$
g_v = g_l
$$

```{code-cell} ipython3

```
