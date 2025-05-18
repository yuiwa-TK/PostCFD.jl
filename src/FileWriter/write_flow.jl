""" 
  write_flow(filename::AbstractString,q::Array{T,4},params::Vector{S};precision::AbstractString) where {T<:Real,S<:Real}

keyword 'precision' ∈ [ "single", "double" ].
If precision == "single", an output file can be visualized by 'FieldView'.

3rd argument 'params' includes (mach , AoA, time, Re)
"""
function write_flow(filename::AbstractString,q::AbstractArray{T,4},params::Vector{S}
                      ;precision::AbstractString) where {T<:Real,S<:Real}
  @show filename
  println("fvmach,alpha,re,time=",params)
  jmax,kmax,lmax,nmax = size(q)

  if precision == "single"
    open(filename,"w") do f
      write(f,Int32(jmax))
      write(f,Int32(kmax))
      write(f,Int32(lmax))
      write(f,Float32.(params))
      write(f,Float32.(q))
    end
  elseif precision == "double"
    open(filename,"w") do f
      write(f,Int32(jmax))
      write(f,Int32(kmax))
      write(f,Int32(lmax))
      write(f,Float64.(params))
      write(f,Float64.(q))
    end
  else
    @error "invalid" precision
  end
  return filename
end

"""
  write_restart(filename::String,q::Array{T,4},params::Vector{S},nc::Int;verbose=2) where {T<:Real,S<:Real}

+Arg1: filename 
+Arg2: q
+Arg3: params (mach , AoA, time)
+Arg4: nc
"""
function write_restart(filename::String,q::AbstractArray{T,4},params::Vector{S},nc::Int;verbose=2) where {T<:Real,S<:Real}
  if verbose>=1;  @show filename;end
   jmax,kmax,lmax,nmax = size(q)

  if verbose==2
    println("jmax,kmax,lmax",jmax,kmax,lmax)
    println("fvmach,alpha,time,nc=",params,nc)
  end

  open(filename,"w") do f
    write(f,Int32(jmax))
    write(f,Int32(kmax))
    write(f,Int32(lmax))
    write(f,Float64.(params))
    write(f,Int32(nc))
    write(f,Float64.(q))
  end
end