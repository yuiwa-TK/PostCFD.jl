""" 
  write pl3d file to visualize by field view.

"""
function write_flow_single(filename::String,q::Array{T},params::Vector{S}) where {T<:Real,S<:Real}
  @show filename
  println("fvmach,alpha,re,time=",params)
  jmax,kmax,lmax,nmax = size(q)
  open(filename,"w") do f
    write(f,Int32(jmax))
    write(f,Int32(kmax))
    write(f,Int32(lmax))
    write(f,Float32.(params))
    write(f,Float32.(q))
  end
  return filename
end
write_flow_fv = write_flow_single

"""
  write restart flow 

  write_restart_double(filename::String,q::Array{T},params::Array{T},nc::Int32) where T
  +Arg1: filename 
  +Arg2: q
  +Arg3: params (mach , AoA, time)
  +Arg4: nc
"""
function write_restart_double(filename::String,q::Array{T,4},params::Vector{S},nc::Int) where {T<:Real,S<:Real}
  @show filename
  println("fvmach,alpha,time,nc=",params,nc)
  jmax,kmax,lmax,nmax = size(q)

  open(filename,"w") do f
    write(f,Int32(jmax))
    write(f,Int32(kmax))
    write(f,Int32(lmax))
    write(f,Float64.(params))
    write(f,Int32(nc))
    write(f,Float64.(q))
  end
end