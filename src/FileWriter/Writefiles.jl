module FileWriter

include("./write_flow.jl")
include("./write_grid.jl")
include("./write_function.jl")

export writefiles

"""
  writefiles(filename, data...; mode) 

is the wrapper fucntion for writting files related to HPC_TEMPLATE
mode ∈ [grid_fv(=grid_single),
          flow_fv(=pl3d, =flow_single),
          grid_double,
          restart]
"""
function writefiles(filename::String,xyz::AbstractArray{T}; mode::String) where T<:Real
  @show mode

  if mode == "grid_fv" || mode == "grid_single"
    return write_grid_single(filename,xyz)
  elseif mode == "grid_double"
    return write_grid_double(filename,xyz)
  end
end

function writefiles(filename::String,q::AbstractArray{T},params::AbstractArray{T}; mode::String) where T<:Real
  @show mode

  if mode == "flow_fv" || mode == "pl3d" || mode =="flow_single" 
    return write_flow_fv(filename,q,params)
  end
end

function writefiles(filename::String,q::AbstractArray{T},params::AbstractArray{T},nc::Int; mode::String) where T
  @show mode
  if mode == "restart"
    return write_restart(filename,q,params,nc)
  end
end

end