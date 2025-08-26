# PostCFD

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://yuiwa-TK.github.io/PostCFD.jl/dev/)
<!-- [![Build Status](https://travis-ci.com/yuiwa-TK/PostCFD.jl.svg?branch=main)](https://travis-ci.com/yuiwa-TK/PostCFD.jl)
[![Coverage](https://codecov.io/gh/yuiwa-TK/PostCFD.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yuiwa-TK/PostCFD.jl) -->

## Install

```julia
using Pkg;Pkg.add(https://github.com/yuiwa-TK/PostCFD.jl.git)
```

## Quick reference
To use the functions of PostCFD.jl, run
```julia
using PostCFD
```

#### To read a grid/flow/function/restart file (FileReader)

```julia
# Grid data
file_grid="path/to/gridfile"
G = FileReader.read_grid(file_grid); # Array of [jmax,kmax,lmax,3]
# or 
G = read_grid(file_grid)

# Flow data
file_flow="path/to/flowfile"
F = FileReader.read_flow(file_flow) # Array of [jmax,kmax,lmax,5]

# If you want to read a function/restart file,
F = FileReader.read_function(filepath)
F = FileReader.read_restart(filepath)

# If you set verbose=0, no message is shown.
G = FileReader.read_grid(file_grid; verbose=0)
```

#### To write a grid/flow/function file (FileWriter)
```julia
# For a grid file (w/o record marker)
FileWriter.write_grid(filename, G; precision="single") # single precision
# or
FileWriter.write_grid(filename, G; precision="double") # double precision

# Similarly, 
write_function(filename, F; precision="single") # double precision
# - you can read this file by `fv`. 

# For a flow file, `params` is required as the third argument.
params = [fvmach,alpha,re,time]; # `params =zeros(4)` works, if these params are not important
write_flow(filename, F, params; precision="single", verbose=0)  #single precision;
# - you can read this file by `fv`. 
# - verbose also works.

# To write a restart file, 'params' of [fvmach,alpha,time] and 'nc' is also required.
params = [fvmach,alpha,time] # Float
nc = nstep_restart # Int 
write_restart(filename, F, params, nc; verbose=0) 
```

#### To take derivatives (Mathlib)

- 1st order sided
- 2nd order central
- 6th-order compact differencing


```julia

fx(xx) = 0.5.*xx.*xx; # 0.5x^2
x = collect(0.0:0.01:1.0);
ff= fx(x);

df = derivative_2ndcentral(ff,x);
# julia> df≈x
#  true

# 6th-order compact diff is also implimented.
df = derivative_compact_6th(ff,x)
```

#### To compute Jacobian and metrics for curvilinear grids (Geometry)
```julia

func_for_deriv = Mathlib.derivative_2ndcentral
J, Met = metrics(G, func_for_deriv)

# size(Met) = jmax,kmax,lmax,3,3
#       --             --
#      |xi_x eta_x zeta_x|   [1,1] [1,2] [1,3]
# Met= |xi_y eta_y zeta_y|   [2,1]  ...    :
#      |xi_z eta_z zeta_z|   [3,1]  ...  [3,3]
#       --             --

# # Fast-computing implementation by reducing memory allocation ... (in development)
J = Jacobian_fast(G, func_for_deriv; )
J, Met = metrics_fast(G, func_for_deriv; )

# For fast-computing implementation with 6th-order compact differencing
J = Jacobian_fast_compact6th(G) 
```