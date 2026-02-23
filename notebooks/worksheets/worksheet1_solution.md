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

(worksheet1_solution)=
# Week 1 worksheet solutions

Three problems from Lohmann Chapter 1

Download link for this file [worksheet1_solution.ipynb](https://drive.google.com/file/d/1dZ-ZcxA2w4rVilm_Gcm9ojpxvRI69Rh4/view?usp=sharing)

+++

4) Assume a low level cloud with the following properties: vertical extent 1.5 km, horizontal
extent 1 km, cloud droplet number concentration 90 $cm^{−3}$ air and mean diameter 10
μm. Estimate the mass of the cloud by assuming that all the cloud droplets are spherical
and have the same size. Take 1000 $kg\,m^{−3}$ to be the density of water.

+++

4. Answer

```{code-cell} ipython3
cloud_diameter=1000 #m
cloud_height = 1500 #m
drop_density = 90e6 # /m^3
drop_diameter = 10.e-6  #m
rho = 1000  #kg/m^3
pi = 3.1416
cloud_volume = cloud_height*pi*(cloud_diameter/2.)**2.
print(f"{cloud_volume*1.e-6=}")
drop_volume = 4./3.*pi*(drop_diameter/2.)**3.
total_number = cloud_volume*drop_density
print(f"{total_number=}")
total_mass = total_number*drop_volume*rho #kg
total_mass_tonnes = total_mass/1000.
print(f"Answer: {total_mass_tonnes:0.1f} tonnes")
```

5) Imagine a pure water cloud which consists of cloud droplets of uniform size with $R_d$
=5 μm. Assume a uniform spacing of the cloud droplets, with a number concentration of
Nc = 170 $cm^{-3}$

   (a) Calculate the average distance between the cloud droplets
   
   (b) What would happen to the distance if the number concentration was trippled?

+++

5. Answer

Assume that each drop gets its own cube inside the 1  $cm^{3}$ of air.  How many cubes does it take to pack a cm?

```{code-cell} ipython3
N = 170.
drop_cube_volume = 1/N
cube_edge_size = drop_cube_volume**0.3333
print(f"{cube_edge_size=:0.2f} cm")
```

So if the drop diameter is 10 microns = $10^{-3}$ cm, then there are about 0.18/0.001 = 180 drop diameters between droplets.

b) If the concentration is trippled, then the drop volume is smaller by a factor of 3, and the edge is shorter by $1/3^{1/3}$ = 0.7

+++

3. The average pressure at Earth’s surface is 985 hPa. Knowing the radius of the Earth
$r_{Earth}$ = 6370 km and using the definition of pressure, estimate the total mass of the
atmosphere.

+++

3. Answer

- pressure = force of gravity per unit area
- force of gravity = mass $\times$ acceleration

So find the mass M of a 1 $m^2$ column 

$$M= \int_0^\infty \rho dz$$

by integrating the hydrostatic equation

$$
dp = -\rho g dz \\
\int_{psfc}^0 dp = 0 - p_{sfc} = \int_0^\infty -\rho g dz \approx -g M\\
M = p_{sfc}/g
$$

```{code-cell} ipython3
pressure = 98500  #pascals
g = 9.8 # m/s/s
mass_m2 = pressure/g
r_earth = 6370e3  #m
earth_area = 4*pi*r_earth**2.
mass_atm = mass_m2*earth_area #kg
mass_atm_gtonne = mass_atm/1.e12
print(f"{mass_atm_gtonne:0.1f} gtonnes")
```

```{code-cell} ipython3

```
