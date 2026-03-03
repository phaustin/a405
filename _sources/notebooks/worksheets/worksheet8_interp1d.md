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

# Worksheet 8: Interp1d

Modify a copy of the [entraining_plume.ipynb](https://drive.google.com/file/d/1rLSqA6c6DqP_DVNW2ie2HfERtL4_t002/view?usp=sharing) to add a new function called  `calc_interp`

that contains the following lines from

    press = df_sounding['pres'].values
    height = df_sounding['hght'].values
    temp = df_sounding['temp'].values
    dewpoint = df_sounding['dwpt'].values
    #
    # 
    #
    envHeight= nudge(height)

    interpTenv = interp1d(envHeight,temp)
    interpTdEnv = interp1d(envHeight,dewpoint)
    interpPress = interp1d(envHeight,press)

```{code-cell} ipython3

```
