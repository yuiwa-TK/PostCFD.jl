using Printf
"""
"DOUBE PRECISION GRID"
read a file with xyz format written in "little-endian" & "stream".  
This function returns xyz=Array{Float64}(jmax,kmax,lmax,3).
"""
function read_grid_double(filename::AbstractString; verbose=2)
    if verbose>=1
        @info filename
    end
    xyz = 0
    dims= Array{Int32}(undef,(3))
    open(filename,"r") do io 
        read!(io,dims)
        if verbose>=2
            @show dims
        end
        jmax = dims[1]
        kmax = dims[2]
        lmax = dims[3]
        xyz = Array{Float64}(undef,(jmax,kmax,lmax,3))
        read!(io,xyz)
    end
    return xyz
end

"""
"SIGLE PRECISION GRID"
read a file with xyz format written in "little-endian" & "stream".  
This function returns xyz=Array{Float32}(jmax,kmax,lmax,3).
"""
function read_grid_single(filename::AbstractString;verbose=2)
    if verbose>=1
        @info filename
    end
    xyz = 0
    dims= Array{Int32}(undef,(3))
    open(filename,"r") do io 
        read!(io,dims)
        if verbose>=2
            @show dims
        end
        jmax = dims[1]
        kmax = dims[2]
        lmax = dims[3]
        xyz = Array{Float32}(undef,(jmax,kmax,lmax,3))
        read!(io,xyz)
    end
    return xyz
end
read_grid_fv = read_grid_single


"""
Read grid size.
This program return Tuple:(jmax,kmax,lmax)
"""
function read_grid_dims(filename::AbstractString;verbose=2)
    if verbose>=1
        @info filename
    end
    dims = Array{Int32}(undef,(3))
    open(filename,"r") do io
        read!(io,dims)
    end
    return dims
end

"""
    read_grid_auto(filename::AbstractString)
automatically determines the file type written in pl3d format.
"""
function read_grid_auto(filename::AbstractString;verbose=2)
    Nb_INT32 = 4
    NBF_FLOAT64 = 8
    NBF_FLOAT32 = 4

    Npoints = prod(read_grid_dims(filename;verbose=verbose))
    Nb_file = filesize(filename)

    if Nb_file == 3*Nb_INT32 + 3*Npoints*NBF_FLOAT32
        return read_grid_single(filename;verbose=verbose)
    elseif Nb_file == 3*Nb_INT32 + 3*Npoints*NBF_FLOAT64
        return read_grid_double(filename;verbose=verbose)
    else
        @error println("$filename is not written in pl3d format or written with record marker .")
        return NaN
    end
end
read_grid=read_grid_auto


function typeof_gridfile(filename::AbstractString; verbose=2)
    Nb_INT32 = 4
    NBF_FLOAT64 = 8
    NBF_FLOAT32 = 4
    Nvars = prod(read_grid_dims(filename; verbose=0))
    Nb_file = filesize(filename)

    if Nb_file == 3*Nb_INT32 + 3*Nvars*NBF_FLOAT32
        if verbose>=1
            println("$filename is single format")
        end
        return "single"
    elseif Nb_file ==3*Nb_INT32 + 3*Nvars*NBF_FLOAT64
        if verbose>=1
            println("$filename is double format")
        end
        return "double"
    else
        @error println("$filename is not written in pl3d format or written with record marker .")
    end
    return 0
end
function read_grid_specifying_xyz(filename::String,iddir::Int; verbose=2)
    NBF_FLOAT64 = 8
    NBF_FLOAT32 = 4
    if verbose>=1
        @show filename
    end
    tp = typeof_gridfile(filename)
    dims = Array{Int32}(undef,(3))

    if tp=="single"
        io   = open(filename,"r") 
        read!(io,dims)
        qvar = Array{Float32}(undef,(dims[1],dims[2],dims[3]))
        Nb_skip = prod(@view dims[1:2])*(iddir-1)*NBF_FLOAT32
        skip(io,Nb_skip)
        read!(io,qvar)
        close(io)
        return qvar
    elseif tp=="double"
        io   = open(filename,"r") 
        read!(io,dims)
        qvar = Array{Float64}(undef,(dims[1],dims[2],dims[3]))
        Nb_skip = prod(@view dims[1:2])*(iddir-1)*NBF_FLOAT64
        skip(io,Nb_skip)
        read!(io,qvar)
        close(io)
        return qvar
    else
        return nothing
    end
end

"""
read_grid_specifying_xyz_rect(filename::String,iddir::Int; verbose=2) returns the vector containing the data of specified direction assuming the rectangular cell.
"""
function read_grid_specifying_xyz_rect(filename::String,iddir::Int; verbose=2)
    NBF_FLOAT64 = 8
    NBF_FLOAT32 = 4
    if verbose>=1
        @show filename
    end
    tp = typeof_gridfile(filename)
    dims = Array{Int32}(undef,(3))

    if tp=="single"
        io   = open(filename,"r")
        read!(io,dims)
        idgs=similar(dims);trimed_dims=similar(dims)
        if iddir == 1
            trimed_dims= [dims[1], 1, 1]
            idgs=[1:dims[1],1,1]
        elseif iddir == 2
            trimed_dims= [dims[1], dims[2], 1]
            idgs=[1,1:dims[2],1]
        elseif iddir == 3
            trimed_dims= dims
            idgs=[1,1,1:dims[3]]
        end
        Nb_skip = prod(@view dims[1:2])*(iddir-1)*NBF_FLOAT32
        skip(io,Nb_skip)
        qvar = Array{Float32}(undef,trimed_dims...)
        read!(io,qvar)
        close(io)
        return Base.getindex(qvar, idgs...)
    elseif tp=="double"
        io   = open(filename,"r")
        read!(io,dims)
        idgs=similar(dims);trimed_dims=similar(dims)
        if iddir == 1
            trimed_dims= [dims[1], 1, 1]
            idgs=[1:dims[1],1,1]
        elseif iddir == 2
            trimed_dims= [dims[1], dims[2], 1]
            idgs=[1,1:dims[2],1]
        elseif iddir == 3
            trimed_dims= dims
            idgs=[1,1,1:dims[3]]
        end
        Nb_skip = prod(@view dims[1:2])*(iddir-1)*NBF_FLOAT64
        skip(io,Nb_skip)
        qvar = Array{Float64}(undef,trimed_dims...)
        read!(io,qvar)
        close(io)
        return Base.getindex(qvar, idgs...)
    else
        return nothing
    end
end