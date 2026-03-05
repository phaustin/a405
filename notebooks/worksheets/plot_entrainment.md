---
jupytext:
  formats: md:myst,ipynb
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

(entrain_plot)=
# Plotting the entraining plume

- Run the [entraining plume](https://phaustin.github.io/a405_2024/notebooks/entraining_plume.html) notebook three times to produce three netcdf files with entrainment rates of $\lambda = 0$, $\lambda = 2 \times 10^{-4}$ and $\lambda = 4 \times 10^{-4}$ $S^{-1}$

- Do you see a difference in the cloud top height and maximum updraft speed for the three runs?

- Make a plot, with a legend, that shows $\theta_e$ vs. height for the three cases

- Repeat, but with vertical velocity

```{code-cell} ipython3
import xarray as xr
from matplotlib import pyplot as plt
from pathlib import Path
```

```{code-cell} ipython3
curr_dir = Path()
nc_files = list(curr_dir.glob("*.nc"))
print(f"{nc_files=}")
the_file = nc_files[0]
```

```{code-cell} ipython3
the_ds = xr.open_dataset(the_file)
the_ds
```

```{code-cell} ipython3
fig, ax = plt.subplots(1,1)
ax.plot('thetae_cloud','cloud_height',data = the_ds)
ax.plot('env_thetae','env_height',data = the_ds)
ax.grid(True)
```

```{code-cell} ipython3
fig, ax = plt.subplots(1,1)
ax.plot('thetae_cloud','cloud_height',data = the_ds)
```

```{code-cell} ipython3

```
