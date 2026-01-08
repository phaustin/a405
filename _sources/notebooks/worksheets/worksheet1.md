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

(worksheet1)=
# Week 1 worksheet

Download the ipynb file for notebook:  [worksheet1 ipynb link](https://drive.google.com/file/d/1rlzSExCyxtXlVCjo2HlxQZ6hWzrZpahL/view?usp=drive_link)

Credit/no credit (i.e. try your best)

Three problems from Lohmann Chapter 1: 4, 5, 6

Upload your finished ipynb file to canvas.  Due Friday midnight.

+++

4. Assume a low level cloud with the following properties: vertical extent 1.5 km, horizontal
extent 1 km, cloud droplet number concentration 90 $cm^{−3}$ air and mean diameter 10
μm. Estimate the mass of the cloud by assuming that all the cloud droplets are spherical
and have the same size. Take 1000 $kg\,m^{−3}$ to be the density of water.

+++

4. Answer

```{code-cell} ipython3
# your code here
```

5) Imagine a pure water cloud which consists of cloud droplets of uniform size with radius $R_d$
=5 μm. Assume a uniform spacing of the cloud droplets, with a number concentration of
Nc = 170 $cm^{-3}$

   (a) Calculate the average distance between the cloud droplets, both in microns and in drop diameters.
   
   (b) What would happen to the distance if the number concentration was trippled?

+++

5. Answer

Hint: Assume that each drop gets its own cube inside the 1  $cm^{3}$ of air.  How many cubes does it take to pack a cm?

```{code-cell} ipython3
# your code here
```

So if the drop diameter is 10 microns, then there are about 180 drop diameters between droplets

+++

6. The average pressure at Earth’s surface is 985 hPa. Knowing the radius of the Earth
$r_{Earth}$ = 6370 km and using the definition of pressure, estimate the total mass of the
atmosphere.

+++

6. Answer

Hint:  use the hydrostatic equation

- pressure = force of gravity per unit area
- force of gravity = mass $\times$ acceleration per unit area

You'll need to find the total mass of a 1 $m^2$ column of the atmosphere by integrating the hydrostatic equation (Stull 1.25c) or

```{code-cell} ipython3
# your code here
```

```{code-cell} ipython3

```
