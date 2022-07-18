"""
  write grid file for fieldview

  function  write_grid_single(filename::String,xyz::Array{T})
    Arg1: filename
    Aeg2: grid data 
"""
function  write_grid_single(filename::String,xyz::Array{T}) where T<:Real
  @show filename
  jmax,kmax,lmax,n = size(xyz)
  @show jmax,kmax,lmax
  open(filename,"w") do f
    write(f,Int32(jmax))
    write(f,Int32(kmax))
    write(f,Int32(lmax))
    write(f,Float32.(xyz))
  end
end
write_grid_fv = write_grid_single

function  write_grid_double(filename::String,xyz::Array{T}) where T<:Real
  @show filename
  jmax,kmax,lmax,n = size(xyz)
  @show jmax,kmax,lmax
  open(filename,"w") do f
    write(f,Int32(jmax))
    write(f,Int32(kmax))
    write(f,Int32(lmax))
    write(f,Float64.(xyz))
  end
end