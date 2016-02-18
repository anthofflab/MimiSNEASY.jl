# Sneasy.jl

This is an implementation of SNEASY in Julia. Everything is experimental at this point and this is not the official SNEASY code.

## Requirements

The minimum requirement to run SNEASY.jl is [Julia](http://julialang.org/) and the [Mimi](https://bitbucket.org/davidanthoff/mimi.jl) package. To plot you should also install the [Gadfly](https://github.com/dcjones/Gadfly.jl) package. To run the example IJulia notebook file you need to install [IJulia](https://github.com/JuliaLang/IJulia.jl).

## Getting started

``examples/sneasy_demo.ipynb`` contains a simple example of how to run SNEASY.jl. It requires a working IJulia installation to run.

``src/main.jl`` runs SNEASY.jl once and is a good starting point to understand the code.

## Folder structure

The ``src`` folder has all the model code.

The ``calibration`` folder has the MCMC assimilation code. The code there can assimilate both the original Fortran implementation of SNEASY and SNEASY.jl.

The ``data`` folder has data files that are needed to run SNEASY. It does *not* contain the data files that are used as observational constraints, those are located in the ``calibration/data`` folder.

The ``examples`` folder has some example code that uses SNEASY.jl.

The ``test`` folder has code that cross checks the julia implementation of SNEASY with the Fortran implementation.
