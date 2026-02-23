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

(worksheet4_taylor)=
# Worksheet 4: Taylor series

This is a pencil and paper problem about expanding the equation of state in a Taylor series

+++

## Tayor series examples

### One dimension

From the equation sheet:

$$
  f(x)  &=& f(x_0) + f^\prime(x_0)(x - x_0) \nonumber\\ 
        &+&  \frac{f^{\prime\prime}(x_0)}{2}(x-x_0)^2 +  \ldots
$$

#### Example:

Show that to first order, if x is close to 0:

$$
\frac{1}{1-x} \approx 1 + x
$$

#### Answer

$$
f(x_0) = \frac{1}{(1 - 0)} = 1 \\
f^\prime(x_0) = -(1 - x_0)^{-2} \times (-1) = 1 \\
f^{\prime \prime}(x_0) = (-2) \times (1 - x_0)^{-3} \times (-1) = -2 \\
f(x) = 1 + x - \frac{2}{2} x^2 + \ldots = 1 + x - x^2 + \ldots
$$

+++

### Two dimensions: Taylor's series expansion of the equation of state

+++

Consider this equation from the
[Wikipedia entry on Taylor series](http://en.wikipedia.org/wiki/Taylor_series)

$$
\begin{align}
f(x,y) & \approx f(a,b) +(x-a)\, f_x(a,b) +(y-b)\, f_y(a,b)\nonumber \\
&  + \frac{1}{2!}\left[ (x-a)^2\,f_{xx}(a,b) + 2(x-a)(y-b)\,f_{xy}(a,b) +(y-b)^2\, f_{yy}(a,b) \right]
\end{align}
$$ (eq:taylor)
where $f_{xy} = \frac{ \partial^2 f}{\partial x \partial y }$, etc. If you expand $f=p=\rho R_d T$ about the point 
$p_0(z) = \rho_0(z) R_d T_0(z)$ where $p_0,\ \rho_0,\ T_0$ are the pressure,
density and temperature at height $z$ for a hydrostatic atmosphere.  Using {eq}`eq:taylor` with $a=\rho_0$ and $b=T_0$ and $f(a,b) = p_0$ you should be able to show that
to second order: 


$$
\frac{\Delta p}{p_0} = \frac{\Delta T}{T_0} + \frac{\Delta \rho}{\rho_0} + \frac{\Delta T \Delta \rho}{T_0 \rho_0}
$$ (eq:full)

Where $\Delta p = p - p_0$, etc. Note that $\Delta p$, $\Delta T$, and $\Delta \rho$ are all functions of
(t,x,y,z).

If the atmosphere is close to hydrostatic balance, then we can expect the $\Delta$ differences to be small if $p_0$ is the hydrostatic pressure, and we can drop the
$\frac{ \Delta T \Delta \rho}{T_0 \rho_0}$ term and write

$$
\frac{\Delta p}{p_0} = \frac{\Delta T}{T_0} + \frac{\Delta \rho}{\rho_0}
$$(eq:hydro)

We will show later that away from active convection we also can expect $\frac{\Delta p }{p_0}$ to be small.

+++

## Problem -- on a piece of paper, show that {eq}`eq:hydro` is true

Upload a scanned image to canvas

```{code-cell} ipython3

```
