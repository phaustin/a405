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

# ATSC 405 2024 final exam solution

+++

## Instructions 

- Do all questions (note weights) showing all you work.  Don't forget to put your name on your tephigram.

+++

## Q1 (10) Cooling 

A cloud at a pressure of 600 hPa and temperature of 5 deg C has 12 g/kg of total water.  Suppose it cools by 10 K and rains so that at -5 deg C it has 1 g/kg of liquid water remaining.  Find

1. The total precipitated water, in g/kg
2. The total change in moist static energy $h_m$ in J/kg
3. The lifting condensation level before and after the cooling, in hPa

+++

## Q2 (8) Taylor series

  Recall that buoyancy for a parcel is defined

$$
  buoyancy = g \frac{\rho_{env} - \rho_{parcel}}{\rho_{env}} \approx g \frac{T_{vparcel} - T_{venv}}{T_{venv}}
$$


1. Use a Taylor series expansion to derive \eqref{eq:buoy} showing the second order terms you drop.  (Hint, expand in the small parameter $\Delta T_v/T_{venv}$, where $\Delta T_v = T_{vparcel} - T_{venv}$)

2. Suppose a parcel at $T_v$=281 K maintains a 1 degree temperature difference over an environment at $T_{venv}$=280 K for an extended period of time.  How many minutes would it take for the parcel to reach a velocity of 10 m/s?

+++

## Q3  (8) Drizzle

Suppose you know that cloud drops have a size distribution given
    by 

$$
    N(r)=N_0 \exp (- \chi r)
$$

where $N(r)$ ($\mathrm m^{-3}$) is the number of drops per unit volume with radius $>
r$ (for radii $50\ \mathrm{\mu m} < r < 200\ \mathrm{\mu m}$ ) and $N_0$
and $\chi$ are constants and $r$ is in $\mu m$.  Suppose also that the drop fall speed
$v(r)\ (\mathrm{m\,s^{-1}})$ 
depends on linearly on radius:
$$
  v(r)=Jr
$$
where $J= 6000\ s^{-1}$ and now $r$ is in meters.

Derive equations (including unit conversions if necessary) for:


- The number density $n(D)$ ($m^{-3}\,\mu m^{-1}$)
- The liquid water content ($\mathrm{kg\,m^{-3}}$)
- The precipitation rate ($\mathrm{mm\,hr^{-1}}$)

+++

## Q4 (10) Collision/coalescence

   Derive with the help of a sketch the simple model of the growth of rain drops represented by
   equation (\ref{eq:collec}), defining all your symbols.  If the fall speed is given by $J= 6000r \ m\,s^{-1}$ as in (\ref{eq:collec}), roughly how long would it take a 50 \mum radius drop to double to 100 \mum radius, given $w_l = 0.3\ g\,m^{-3}$?  What processes need to be added to
    (\ref{eq:collec}) to make the description of rain formation more
    accurate? Explain.

```{code-cell} ipython3

```

## Q5 (11) Koehler curve


```{figure} ./m18.png
---
width: 70%
name: directive-fig
alt: pha
---
Koehler curve
```


The figure shows two droplets with radii of 0.2 $\mu m$ and
1 $\mu m$ at a environmental relative humidity of about 1.001.

- (3) Both droplets lie on a curved line called the Koehler curve.  What
do points on this line represent?  i.e.  what is it that all points lying on
this line have in common?

- (4) Suppose the relative humidity increases from 1.001 to 1.002.  Describe
qualitatively what happens to the two drops.


- (4) Quantitatively, what is the initial evaporation/growth rate of droplet
2 in microns/second, assuming a RH of 1.002, temperature T=280, esat(T)=10 hPa.

+++

## Q6 (10) Adiabatic ascent

  Suppose someone asks you to provide the temperature, pressure, water vapor mixing ratio and liquid water mixing ratio as a function of height for surface air at press= 1000 hPa, temperature T and relative humidity RH as it ascends adiabatically through the atmosphere.  Using pseudo code sketch out a python script that would do this assuming a hydrostatic balance.  Assume you can use odeint and the rootfinder as needed.   Note:  1) the equations you would need to solve and 2) how you would solve them.

+++

```{figure} ./tephi_final_2024.jpg
---
width: 90%
name: fig_tephi
alt: pha
---
tephigram
```

```{code-cell} ipython3

```
