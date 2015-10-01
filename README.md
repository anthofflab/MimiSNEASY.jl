# Sneasy.jl

This is an implementation of SNEASY in Julia. Everything is experimental at this point and this is not the official SNEASY code.

## Requirements

The minimum requirement to run SNEASY is [Julia](http://julialang.org/) and the [Mimi](https://bitbucket.org/davidanthoff/mimi.jl) package. To plot you should also install the [Gadfly](https://github.com/dcjones/Gadfly.jl) package. To run the example IJulia notebook file you need to install [IJulia](https://github.com/JuliaLang/IJulia.jl).

### Installing Julia

You can download Julia from [http://julialang.org/downloads/](http://julialang.org/downloads/). You should use the latest v0.3.x version and install it.

### Installing the Mimi package

Start Julia and enter the following command on the Julia prompt:

````julia
Pkg.add("Mimi")
````

### Installing Gadfly (optional)

[Gadfly](https://github.com/dcjones/Gadfly.jl) is one of many Julia plotting packages. To install it, start Julia and enter the following command on the Julia prompt:

````julia
Pkg.add("Gadfly")
````

### Installing IJulia (optional)

[IJulia](https://github.com/JuliaLang/IJulia.jl) provides the excellent [IPython](http://ipython.org/) interface for the Julia language. Installation is a two step procedure:

First, you need a working IPython installation. The Anaconda distribution [(http://continuum.io/downloads)](http://continuum.io/downloads) is currently the easiest way to install IPython.

Second, start Julia, then enter the following command on the Julia prompt:

````julia
Pkg.add("IJulia")
````

To start the IJulia notebook, open a normal command prompt (not Julia), change to the directory with the SNEASY source code and enter the following command:

````
ipython notebook --profile julia
````

### Keeping requirements up-to-date (optional)

Many of these requirements are regularily updated. To make sure you have the latest versions, periodically execute the following command on the Julia prompt:

````jl
Pkg.update()
````
