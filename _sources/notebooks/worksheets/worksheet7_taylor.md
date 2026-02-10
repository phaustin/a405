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

(worksheet7_taylor)=
# Worksheet 7 -- Taylor series problem

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

2) Using this approximation, integrate the moist entropy equation to show that

$$
\phi = c_p \log \theta_e 
$$

```{code-cell} ipython3

```
