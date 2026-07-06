""" 
  write_grid(filename::AbstractString,q::AbstractArray{T,4};precision::AbstractString, verbose=2) where {T<:Real,S<:Real}

keyword 'precision' ∈ [ "single", "double" ].
If precision == "single", an output file can be visualized by 'FieldView'.

"""
function write_grid2d(filename::AbstractString,xyz::AbstractArray{T}
                      ;precision::AbstractString, verbose=2) where T<:Real
  if verbose>=1; @show filename; end
  if verbose>=2
    jmax,kmax = size(xyz)
    @info jmax,kmax
  end

  if precision == "single"
    write_grid2d_single(filename,xyz;verbose=0)
  elseif precision == "double"
    write_grid2d_double(filename,xyz;verbose=0)
  end
  return filename
end

"""
  function  write_grid2d_single(filename::AbstractString,xyz::AbstractArray{T,3})
writes 2D grid file with single precision.
  -  Arg1: filename
  -  Arg2: grid data 
"""
function  write_grid2d_single(filename::AbstractString,xyz::AbstractArray{T,3};verbose=2) where T<:Real
  if verbose>=1; @show filename; end
  jmax,kmax,n = size(xyz)
  if verbose>=2; @info jmax,kmax; end
  open(filename,"w") do f
    write(f,Int32(jmax))
    write(f,Int32(kmax))
    write(f,Float32.(xyz))
  end
end

"""
  function  write_grid2d_double(filename::AbstractString,xyz::AbstractArray{T,3})
writes 2D grid file with single precision.
  -  Arg1: filename
  -  Arg2: grid data 
"""
function  write_grid2d_double(filename::AbstractString,xyz::AbstractArray{T,3};verbose=2) where T<:Real
  if verbose>=1; @show filename; end
  jmax,kmax,n = size(xyz)
  if verbose>=2; @info jmax,kmax; end
  open(filename,"w") do f
    write(f,Int32(jmax))
    write(f,Int32(kmax))
    write(f,Float64.(xyz))
  end
end