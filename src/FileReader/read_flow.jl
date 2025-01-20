"""
    qall = read_flow_double(filename::AbstractString)

reads a file with pl3d format written in *little-endian* & *stream*.  
and returns *qall=Array{Float64}(jmax,kmax,lmax,5)*.
"""
function read_flow_double(filename::AbstractString; verbose=2, endian="little")
    if verbose>=1
        @info filename
    end
    # settings ==============================================
    dims = Array{Int32}(undef,(3))
    nov  = 5

    dims    = Array{Int32}(undef,(3))
    params  = Array{Float64}(undef,(4))
    qall    = 0.0

    open(filename,"r") do io 
        read!(io,dims)
        read!(io,params)
        if endian in ["b", "B", "big", "big-endian"]
            dims=ntoh.(dims)
            params=ntoh.(params)
        end
        if verbose>=2
            @show dims
            @show params
        end
        jmax = dims[1]
        kmax = dims[2]
        lmax = dims[3]
        qall = Array{Float64}(undef,(jmax,kmax,lmax,nov))
        read!(io,qall)
    end
    if endian in ["b", "B", "big", "big-endian"]
        return ntoh.(qall)
    else
        return qall
    end
end

"""
    qall = read_flow_single(filename::AbstractString) or
    qall = read_flow_fv(filename::AbstractString)

reads a file with pl3d format written in *little-endian* & *stream*.  
and returns *qall=Array{Float32}(jmax,kmax,lmax,5)*.
"""
function read_flow_single(filename::AbstractString;verbose=2,endian="little")
    if verbose>=1
        @info filename
    end
    # settings ==============================================
    dims = Array{Int32}(undef,(3))
    nov  = 5

    dims    = Array{Int32}(undef,(3))
    params  = Array{Float32}(undef,(4))
    qall    = 0.0
   
    open(filename,"r") do io 
        read!(io,dims)
        read!(io,params)
        if endian in ["b", "B", "big", "big-endian"]
            dims=ntoh.(dims)
            params=ntoh.(params)
        end
        if verbose>=2
            @show dims
            @show params
        end
        jmax = dims[1]
        kmax = dims[2]
        lmax = dims[3]
        qall = Array{Float32}(undef,(jmax,kmax,lmax,nov))
        read!(io,qall)
    end
    if endian in ["b", "B", "big", "big-endian"]
        return ntoh.(qall)
    else
        return qall
    end
end
read_flow_fv = read_flow_single

"""
    read_restart(filename::AbstractString)
read a file with *reatart* format written in *little-endian" & "stream*.  
this function returns *qrestart=Array{Float64}(jmax,kmax,lmax,5)*.
"""
function read_restart(filename::AbstractString;verbose=2,endian="little")
    if verbose>=1
        @info filename
    end
    # settings ==============================================
    dims= Array{Int32}(undef,(3))
    nov  = 5

    dims    = Array{Int32}(undef,(3))
    params  = Array{Float64}(undef,(3))
    nc      = Array{Int32}(undef,1)
    qall    = 0.0
   
    open(filename,"r") do io 
        read!(io,dims)
        read!(io,params)
        read!(io,nc) 
        if endian in ["b", "B", "big", "big-endian"]
            dims=ntoh.(dims)
            params=ntoh.(params)
            nc = ntoh(nc)
        end
        if verbose>=2
            @show dims
            @show params
            @show nc
        end
        jmax = dims[1]
        kmax = dims[2]
        lmax = dims[3]
        qall = Array{Float64}(undef,(jmax,kmax,lmax,nov))
        read!(io,qall)
    end
    if endian in ["b", "B", "big", "big-endian"]
        return ntoh.(qall)
    else
        return qall
    end
end

"""
    dim,param = read_flow_params(filename::AbstractString; type)

returns the dimension and parameter of flowfile written in pl3d format.

"type" : choose from "single", "double", and "restart"
"""
function read_flow_params(filename::AbstractString; type)
    println("Reading flow params")
    @show filename

    if type == "single"
        dims    = Array{Int32}(undef,(3))
        params  = Array{Float32}(undef,(4))
    
        open(filename,"r") do io
            read!(io,dims)
            read!(io,params)
        end
        println("Mach:", params[1], "| AoA:",params[2],"| Time:",params[3],"| Re",params[4])
        return dims,params
        
    elseif type=="double"
        dims    = Array{Int32}(undef,(3))
        params  = Array{Float32}(undef,(3))
        nc      = Array{Int32}(undef,(1))
        open(filename,"r") do io
            read!(io,dims)
            read!(io,params)
            read!(io,nc)
            append!(params, nc);
        end
        println("Mach:", params[1], "| AoA:",params[2],"| Time:",params[3],"| nc",params[4])
        return dims, params
    elseif type=="restart"
        dims    = Array{Int32}(undef,(3))
        params  = Array{Float64}(undef,(3))
        nc      = Array{Int32}(undef,(1))
        open(filename,"r") do io
            read!(io,dims)
            read!(io,params)
            read!(io,nc)
            append!(params, nc);
        end
        println("Mach:", params[1], "| AoA:",params[2],"| Time:",params[3],"| nc",params[4])
        return dims, params
    end
end

"""
Read grid size from a flow file.
This program return Tuple:(jmax,kmax,lmax)
"""
function read_flow_dims(filename::AbstractString, endian="little")
    dims = Array{Int32}(undef,(3))
    open(filename,"r") do io
        read!(io,dims)
    end
    if endian in ["b", "B", "big", "big-endian"]
        return ntoh.(dims)
    else
        return dims
    end
end

"""
    read_flow_auto(filename::AbstractString)
automatically determines the file type written in pl3d format.
"""
function read_flow_auto(filename::AbstractString; verbose=2, endian="little")
    Nb_INT32 = 4
    NBF_FLOAT64 = 8
    NBF_FLOAT32 = 4
    Npoints = prod(read_flow_dims(filename))
    Nb_file = filesize(filename)

    if Nb_file == 3*Nb_INT32 + 4*NBF_FLOAT32 +5*Npoints*NBF_FLOAT32
        return read_flow_single(filename; verbose= verbose, endian=endian)
    elseif Nb_file == 3*Nb_INT32 + 4*NBF_FLOAT64 + 5*Npoints*NBF_FLOAT64
        return read_flow_double(filename; verbose= verbose, endian=endian)
    elseif Nb_file == 3*Nb_INT32 + 3*NBF_FLOAT64 + Nb_INT32 + 5*Npoints*NBF_FLOAT64
        return read_restart(filename;  verbose= verbose, endian=endian)
    else
        @error println("$filename is not written in pl3d format or written with record marker .")
        return NaN
    end
end
read_flow=read_flow_auto