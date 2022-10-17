""" 
  write_function(filename::AbstractString,q::Array{T,4},;precision::AbstractString) where {T<:Real,S<:Real}

keyword 'precision' ∈ [ "single", "double" ].
If precision == "single", an output file can be visualized by 'FieldView'.

"""
function write_function(filename::AbstractString,q::Array{T,4}
                      ;precision::AbstractString) where {T<:Real}
  @show filename
  jmax,kmax,lmax,nvar = size(q)

  if precision == "single"
    open(filename,"w") do f
      write(f,Int32(jmax))
      write(f,Int32(kmax))
      write(f,Int32(lmax))
      write(f,Int32(nvar))
      write(f,Float32.(q))
    end
  elseif precision == "double"
    open(filename,"w") do f
      write(f,Int32(jmax))
      write(f,Int32(kmax))
      write(f,Int32(lmax))
      write(f,Int32(nvar))
      write(f,Float64.(q))
    end
  end
  return filename
end