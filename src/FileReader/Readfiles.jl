module FileReader

include("./read_grid.jl")
include("./read_flow.jl")
include("./read_divided.jl")
include("./read_header.jl")

export readfiles

"""
    readfiles(filename::String; mode::String [, ids]) 

is a wrapper fucntion for reading files
+readfiles(filename, mode)     :: read files
+readfiles(compdir, mode, ids) :: read diveded grid/restart files 
                                    in directory "compdir"

mode ∈ [grid (=grid_double), 
        grid_single,
        flow_fv (=pl3d. flow_single),
        restart,
        dims (=size),
        header,
        ---------------
[!NOTE!]divided_xxx is assumed to be double precision[!NOTE!]
        divided_grid(=grid),
        divided_restart,
        divided_restart_old,
        divided_flow,
        divided_flow_old, ]

[!Note!]
        divided_flow_slice
"""
function readfiles(filename::String; mode::String)
    @show mode
    
    # Grid data
    if     mode == "grid" || mode == "grid_double"
        return read_grid(filename)
    elseif mode == "grid_single" 
        return read_grid_single(filename)

    # Flow data
    elseif mode == "flow" || mode == "flow_double"
        return read_flow_double(filename)
    elseif mode == "pl3d"|| mode == "flow_single" || mode =="flow_fv"
        return read_flow_single(filename)
    
    # Restart flow data
    elseif mode == "restart"
        return read_restart(filename)

    elseif mode == "dims" || mode == "size"
        return read_grid_dims(filename)

    # normal .dat file with header
    elseif mode == "header"
        return read_header(filename)
    end
end
function readfiles(compdir::String, ids::Union{OrdinalRange,Vector{Int}};  mode::String)
    if mode == "grid" ||  mode == "divided_grid"
        return read_dividedgrid(compdir,ids)
    elseif mode == "restart" || mode=="divided_restart" 
        return  read_dividedrestart(compdir,ids,"restart.")
    elseif mode == "restart_old" || mode=="divided_restart_old" 
        return  read_dividedrestart(compdir,ids,"restart_old.")
    elseif mode == "flow" || mode == "divided_flow"
        return  read_dividedrestart(compdir,ids,"flow_z") 
    elseif mode == "flow_old" || mode == "divided_flow_old"
        return  read_dividedrestart(compdir,ids,"flow_old_z") 
    end
end
function readfiles(compdir::String, ids::Union{OrdinalRange,Vector{Int}},filebase::String;  mode::String)
    if mode == "divided_flow_slice"
        return read_dividedrestart(compdir,ids,filebase)
    end
end

end
