""" 
  write_grid(filename::AbstractString,q::AbstractArray{T,4};precision::AbstractString, verbose=2) where {T<:Real,S<:Real}

keyword 'precision' ∈ [ "single", "double" ].
If precision == "single", an output file can be visualized by 'FieldView'.

"""
function write_grid(filename::AbstractString,xyz::AbstractArray{T}
                      ;precision::AbstractString, verbose=2) where T<:Real
  if verbose>=1; @show filename; end
  if verbose>=2
    jmax,kmax,lmax = size(xyz)
    @info jmax,kmax,lmax
  end

  if precision == "single"
    write_grid_single(filename,xyz;verbose=0)
  elseif precision == "double"
    write_grid_double(filename,xyz;verbose=0)
  end
  return filename
end

"""
  function  write_grid_single(filename::AbstractString,xyz::AbstractArray{T})
writes grid file with single precision.
  -  Arg1: filename
  -  Arg2: grid data 
"""
function  write_grid_single(filename::AbstractString,xyz::AbstractArray{T};verbose=2) where T<:Real
  if verbose>=1; @show filename; end
  jmax,kmax,lmax,n = size(xyz)
  if verbose>=2; @info jmax,kmax,lmax; end
  open(filename,"w") do f
    write(f,Int32(jmax))
    write(f,Int32(kmax))
    write(f,Int32(lmax))
    write(f,Float32.(xyz))
  end
end

"""
  write grid file for fieldview

  function  write_grid_fv(filename::AbstractString,xyz::AbstractArray{T})
  -  Arg1: filename
  -  Arg2: grid data 
"""
write_grid_fv = write_grid_single

"""
  function  write_grid_double(filename::AbstractString,xyz::AbstractArray{T})
writes grid file with double precision.
  -  Arg1: filename
  -  Arg2: grid data 
"""
function  write_grid_double(filename::AbstractString,xyz::AbstractArray{T};verbose=2) where T<:Real
  if verbose>=1; @show filename; end
  jmax,kmax,lmax,n = size(xyz)
  if verbose>=2; @info jmax,kmax,lmax; end
  open(filename,"w") do f
    write(f,Int32(jmax))
    write(f,Int32(kmax))
    write(f,Int32(lmax))
    write(f,Float64.(xyz))
  end
end