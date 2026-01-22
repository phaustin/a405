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

(worksheet2)=
# Week 2 worksheet

+++

This notebook has separate problems for Tuesday and Thursday.  Hand in both problems after class on Thursday.

The download link is [worksheet2.ipynb](https://drive.google.com/file/d/1NzrJy3mG9O9CGORxv3Vtujvd9CAQ-jwZ/view?usp=sharing)

+++ {"user_expressions": []}

## Working with the A405 thermodynamics libraries

### Installation

+++

The source code for the libraries is in the [A405 github repository](https://github.com/phaustin/a405/tree/main/src/a405).  To install on your laptop, use this
[requirements.txt](https://github.com/phaustin/a405/blob/main/requirements.txt) with pip  (note only one dash in -r, but two dashes in --update) in the a405 environment.

+++ {"user_expressions": []}

   pip install -r requirements.txt --uprade

+++

or just copy/paste the line inside the requirements.txt file at the prompt

+++

     pip install  git+https://github.com/phaustin/a405.git --upgrade

+++ {"user_expressions": []}

### Documentation

+++ {"user_expressions": []}

The A405 library modules are documented at [https://phaustin.github.io/a405_lib/full_listing.html](https://phaustin.github.io/a405_lib/full_listing.html)

+++ {"user_expressions": []}

## Practice problems

Use the library functions to solve the following problems

+++ {"user_expressions": []}

1)  Find the potential temperature of air at a pressure of 450 hPa and a temperature of 300 K

```{code-cell} ipython3
# your code here (week 2 Tuesday)
```

2. For Thursday: Two jars of liquid water are placed in an insulated vacumn chamber which has completely reflecting surfaces in the longwave, so that the chamber doesn't affect their temperatures. They exchange longwave radiation until they come into equilibrium.  Neglect the water vapor in the chamber and the container material, and assume that Jar A contains 1 kg of liquid and Jar B contains 2 kg of liquid.

    a. I the 1 kg jar A starts at 350 K and the 2 kg Jar B starts at 280 K, what is the final equilibrium temperature of each jar?  
    b. What is the total change in entropy for the system composed of the two jars?  (liquid water is incompressible)

```{code-cell} ipython3
# your code here (week 2 Thursday)
```

```{code-cell} ipython3

```
