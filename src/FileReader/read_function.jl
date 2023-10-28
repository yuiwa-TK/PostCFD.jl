function read_function_single(filename::String)
    @show filename
    # settings ==============================================
    dims = Array{Int32}(undef,(4))
    qall = 0

    open(filename,"r") do io 
        @show read!(io,dims)
        jmax = dims[1]
        kmax = dims[2]
        lmax = dims[3]
        nvar = dims[4]
        qall = Array{Float32}(undef,(jmax,kmax,lmax,nvar))
        read!(io,qall)
    end
    return qall
end

function read_function_double(filename::String)
    @show filename
    # settings ==============================================
    dims = Array{Int32}(undef,(4))
    qall = 0

    open(filename,"r") do io 
        @show read!(io,dims)
        jmax = dims[1]
        kmax = dims[2]
        lmax = dims[3]
        nvar = dims[4]
        qall = Array{Float64}(undef,(jmax,kmax,lmax,nvar))
        read!(io,qall)
    end
    return qall
end

function read_function_dims(filename::String)
    dims = Array{Int32}(undef,(4))
    open(filename,"r") do io 
        @show read!(io,dims)
    end
    return dims
end

"""
    read_flow_auto(filename::AbstractString)
automatically determines the file type written in pl3d format.
"""
function read_function_auto(filename::AbstractString)
    Nb_INT32 = 4
    NBF_FLOAT64 = 8
    NBF_FLOAT32 = 4
    Nvars = prod(read_function_dims(filename))
    Nb_file = filesize(filename)

    if Nb_file == 4*Nb_INT32 + Nvars*NBF_FLOAT32
        return read_function_single(filename)
    elseif Nb_file == 4*Nb_INT32 + Nvars*NBF_FLOAT64
        return read_function_double(filename)
    else
        @error println("$filename is not written in pl3d format or written with record marker .")
        return NaN
    end
end
read_function=read_function_auto