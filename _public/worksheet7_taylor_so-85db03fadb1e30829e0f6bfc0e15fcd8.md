---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.18.1
kernelspec:
  name: python3
  display_name: Python 3 (ipykernel)
  language: python
---

(worksheet7_taylor_solution)=
# Worksheet 7 -- Taylor series problem solution

The link to [worksheet7_taylor_solution.ipynb](https://drive.google.com/file/d/1lD93cY3wWHKVC6ie7WJ_ghZXjlieoOwM/view?usp=sharing)

+++

1) Use a Taylor series to expand the following expression with respect to a temperature perturbation (this should yield 3 terms)
  
$$
d  \left ( \frac{l_v(T) r_s(T)}{c_p T} \right ) 
$$

Use the [assignment 4b tephigram](https://docs.google.com/presentation/d/18DG-wXVmWE2lj38h2tS_axL7-uBQvpiINoUlARlmiaE/edit?usp=sharing) to estimate the size of each of the three terms for a 2 deg C  temperature perturbation at a pressure of 900 hPa and a temperature of 10 deg C.  To first order, what is error in percent caused by dropping two of those terms in the following approximation?  (Hint: get $dr_s/dT$ from the tephigram, not the C-C equation)

$$
d  \left ( \frac{l_v(T) r_s(T)}{c_p T} \right ) \approx  \left ( \frac{l_v(T) }{c_p T} \right ) dr_s
$$

+++

## 1 Answer

Use the chain rule on the three temperature dependent variables:

$$
d\left(\frac{l_v r_s}{c_p T}\right)=\frac{1}{c_p}\left [ \underbrace{\left(\frac{r_s}{T}\right) \frac{dl_v}{dT} d T }_a +
\underbrace{\left(\frac{l_v}{T}\right) \frac{d r_s}{\partial T} d T}_b \\
- \underbrace{\left(\frac{l_v r_s}{T^2}\right) d T}_c \right ] \approx \underbrace{\left(\frac{l_v}{T}\right) \frac{d r_s}{\partial T} d T}_b  = \underbrace{\left(\frac{l_v}{T} dr_s\right) }_b 
$$(eq:approx)

+++

We know the following:

$$
\frac{dl_v}{dT} = c_{pv} - c_l 
$$
$$
\frac{dr_s}{dT} dT = dr_x
$$

+++

put in some numbers from the tephigram

- $r_s$(14 deg) \approx 11.5 g/kg
- $r_s$(16 deg) \approx 12.5 g/gk
- $l_v$ \approx 2.5 $\times 10^6$ J/kg
- $c_{pv}$ = 1870 J/kg/K
- $c_l$ = 4190 J/kg/K

```{code-cell} ipython3
drs = 1.e-3 #kg/kg
dT = 2  # K
rs = 12.e-3 # kg/kg
cpv = 1870 #J/kg/K
cl = 4190 # J/kg/K
lv = 2.5e6 # J/kg
temp = 273.15 + 15 # K
term_a = (rs/temp)*(cpv - cl)*dT
term_b = (lv/temp)*drs
term_c = (lv*rs/temp**2.)*dT
print(f"{term_a=:.1f}, {term_b=:.1f}, {term_c=:.1f}")
```

Conclusion, terms a and c can be dropped if we can live with an error of about 10%

+++

2) Using this approximation, integrate the moist entropy equation to show that



$$
\phi = c_p \log \theta_e 
$$

## 2 Answer

i) Start with the first law of thermodynamics for a saturated process:

$$
q\,dt = dh - \alpha dp =  c_p\, dT\ + l_v\,dr_s - \alpha dp =   c_p\, dT\ + l_v\,dr_s - R_d T \frac{dp}{p}
$$(eq:first_taylor)

where the specific volume $\alpha = 1/\rho = R_d T/p$

+++

ii) Recall the second law:

$$
d\phi = \frac{q\,dt}{T}
$$

So divide {eq}`eq:first_taylor` by T:

$$
d\phi = c_p\, \frac{dT}{T}\ + \frac{l_v}{T}\,dr_s - R_d  \frac{dp}{p}
$$

Now use the approximation {eq}`eq:approx`.

$$
d\phi = c_p\, d \log T\ + d \left ( \frac{l_v}{T}\,r_s \right ) - R_d  d \log p
$$

and using (11) in the midterm equation sheet we can write that as 

$$
d\phi = c_p\, d \log \theta + d \left ( \frac{l_v}{T}\,r_s \right )
$$(eq:dentropy)

So how do we integrate this  along an adiabat?  Set $d\phi = 0$ and integrate from the bottom of the atmosphere where

$$
\theta^\prime = \theta \\
\left ( \frac{l_v}{T}\,r_s \right )^\prime = \left ( \frac{l_v}{T}\,r_s \right )
$$

to the top of the atmosphere where 

$$
\theta^\prime = \theta_e \\
\left ( \frac{l_v}{T}\,r_s \right )^\prime = 0
$$

This gives us:

$$
 - c_p \log \theta + c_p \log \theta_e =  0 - \left ( \frac{l_v}{T}\,r_s \right )
$$

Finally, rearranging and taking the $\exp$ of both sides to get:


+++

$$
\theta_e = \theta \exp \left ( \frac{l_v r_s}{c_p T} \right )
$$(eq:thetae)

+++

Taking the differential of {eq}`eq:thetae` and comparing it to {eq}`eq:dentropy` gives:

$$
d\phi = c_p d\log \theta_e =  c_p\, d \log \theta + d \left ( \frac{l_v}{T}\,r_s \right )
$$

So changes in $c_p \log \theta_e$ are the same as changes in $\phi$

```{code-cell} ipython3

```
