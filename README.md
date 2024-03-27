# EMCC_models
Modelling of Electromicrobiological Concentration Cells (EMCCs)

These are simple models used to generate plots for the paper entitled "Electromicrobiological concentration cells are an overlooked potential energy conservation mechanism for subsurface microorganisms" by Ian P.G. Marshall. The code is written in the [Julia programming language](https://www.julialang.org) using the [Gadfly](http://gadflyjl.org) and [DataFrames](https://dataframes.juliadata.org/) packages.

The concentration profile is defined based on the function `twopart_concentration_profile`, for example:

```
twopart_concentration_profile.(1E-3, 1E-9, 1, 5, 0.5)
```

specifies a concentration profile with a starting concentration of 1E-3 that decreases linearly to 1E-9 over a distance of 1 mm, with the lower concentration from 1-5 mm, calculated for the depth of 0.5 mm (this is calculated for every depth at an interval of 0.01 mm in the examples provided here).

Other important parameters include `z` for the number of electrons transferred per substrate molecule and `voltage_drop_mm` for the energy loss per mm of conductive structure.
