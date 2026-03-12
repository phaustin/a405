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

(worksheet8_interp1d)=
# Worksheet 8: Interp1d

This worksheet is part 1 of Assignment 6 and is due with the rest of the assignment on March 12.

Download [worksheet8_interp1d.ipynb](https://drive.google.com/file/d/1KzMEkSRzKD61yhLW01-hBK8DUjkCEAIR/view?usp=sharing)

Modify a copy of the [entraining_plume.ipynb](https://drive.google.com/file/d/1rLSqA6c6DqP_DVNW2ie2HfERtL4_t002/view?usp=sharing) notebook to add a new function called  `calc_interp`

that contains the following lines from [integ_entrain](https://phaustin.github.io/a405/notebooks/worksheets/entraining_plume.html#define-the-integrator-function)

+++

    press = df_sounding['pres'].values
    height = df_sounding['hght'].values
    temp = df_sounding['temp'].values
    dewpoint = df_sounding['dwpt'].values
    #
    # the nudge function moves any identical heights slightly up or down
    # to avoid breaking the interpolation
    
    envHeight= nudge(height)
    
    interpTenv = interp1d(envHeight,temp)
    interpTdEnv = interp1d(envHeight,dewpoint)
    interpPress = interp1d(envHeight,press)

+++

The signature should look like this:

     def calc_interp(df_sounding):

and it should return the three interperlators interpTenv, interpTdEnv, interpPress.

+++

Check the accuracy of the interpolaters by adding the interpolated soundings onto the environmental sounding plot as a series of large dots.  Are the interpolaters accurately mapping the jagged changes in the temperature and depoint with height?

```{code-cell} ipython3

```
