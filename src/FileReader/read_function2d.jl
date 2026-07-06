function read_function2d_single(filename::String; verbose=2,endian="little")
    if verbose>=1
        @show filename
    end
    # settings ==============================================
    dims = Array{Int32}(undef,(3))
    qall = 0

    open(filename,"r") do io 
        read!(io, dims)
        jmax = dims[1]
        kmax = dims[2]
        nvar = dims[3]
        if endian!="little"
            jmax,kmax,nvar=ntoh.([jmax,kmax,nvar])
        end
        if verbose>=1 ; @show jmax,kmax,nvar; end
        qall = Array{Float32}(undef,(jmax,kmax,nvar))
        read!(io,qall)
    end
    if endian!="little"
        return ntoh.(qall)
    else
        return qall
    end
end

function read_function2d_double(filename::String; verbose=2,endian="little")
    if verbose>=1
        @show filename
    end
    # settings ==============================================
    dims = Array{Int32}(undef,(3))
    qall = 0

    open(filename,"r") do io 
        read!(io, dims)
        jmax = dims[1]
        kmax = dims[2]
        nvar = dims[3]
        if endian!="little"
            jmax,kmax,nvar=ntoh.([jmax,kmax,nvar])
        end
        if verbose>=1 ; @show jmax,kmax,nvar; end
        qall = Array{Float64}(undef,(jmax,kmax,nvar))
        read!(io,qall)
    end
    if endian!="little"
        return ntoh.(qall)
    else
        return qall
    end
end

function read_function2d_dims(filename::String; verbose=2, endian="little")
    dims = Array{Int32}(undef,(3))
    open(filename,"r") do io 
        read!(io,dims)
    end
    if verbose>=2
        @info dims
    end
    if endian!="little"
        ntoh.(dims)
    else
        return dims
    end
end


function typeof_functionfile2d(filename::AbstractString; verbose=2,endian="little")
    Nb_INT32 = 4
    NBF_FLOAT64 = 8
    NBF_FLOAT32 = 4
    Nvars = prod(read_function2d_dims(filename; verbose=verbose,endian=endian))
    Nb_file = filesize(filename)

    if Nb_file == 3*Nb_INT32 + Nvars*NBF_FLOAT32
        if verbose>=1
            println("$filename is single format")
        end
        return "single"
    elseif Nb_file == 3*Nb_INT32 + Nvars*NBF_FLOAT64
        if verbose>=1
            println("$filename is double format")
        end
        return "double"
    else
        @error println("$filename is not written in pl3d format or written with record marker .")
    end
    return 0
end

"""
    read_flow_auto(filename::AbstractString; verbose=2, endian="little")
automatically determines the file type written in pl3d format.

# input
- endian = "little"/"big"
- verbose = 2(show filename), 0(no output)
"""
function read_function2d_auto(filename::AbstractString; verbose=2,endian="little")
    Nb_INT32 = 4
    NBF_FLOAT64 = 8
    NBF_FLOAT32 = 4
    Nvars = prod(read_function2d_dims(filename;verbose=0))
    if endian!="little"
        Nvars = prod(ntoh.(read_function2d_dims(filename;verbose=0)))
    end
    Nb_file = filesize(filename)

    if Nb_file == 3*Nb_INT32 + Nvars*NBF_FLOAT32
        return read_function2d_single(filename; verbose=verbose,endian=endian)
    elseif Nb_file == 3*Nb_INT32 + Nvars*NBF_FLOAT64
        return read_function2d_double(filename; verbose=verbose,endian=endian)
    else
        @error println("$filename is not written in pl3d format or written with record marker .")
        return NaN
    end
end
read_function2d=read_function2d_auto

"""
    read_function2d_specifyingvariable(filename::String,idvar::Int; verbose=2,endian="little")
# input
- endian = "little"/"big"
"""
function read_function2d_specifyingvariable(filename::String,idvar::Int; verbose=2,endian="little")
    NBF_FLOAT64 = 8
    NBF_FLOAT32 = 4
    if verbose>=1
        @show filename
    end
    tp = typeof_functionfile2d(filename; verbose=0,endian=endian)
    dims = Array{Int32}(undef,(3))

    if tp=="single"
        io   = open(filename,"r") 
        read!(io,dims)
        if endian!="little"
            dims=ntoh.(dims)
        end
        qvar = Array{Float32}(undef,(dims[1],dims[2]))
        Nb_skip = prod(@view dims[1:2])*(idvar-1)*NBF_FLOAT32
        skip(io,Nb_skip)
        read!(io,qvar)
        close(io)
        return qvar
    elseif tp=="double"
        io   = open(filename,"r") 
        read!(io,dims)
        if endian!="little"
            dims=ntoh.(dims)
        end
        qvar = Array{Float64}(undef,(dims[1],dims[2]))
        Nb_skip = prod(@view dims[1:2])*(idvar-1)*NBF_FLOAT64
        skip(io,Nb_skip)
        read!(io,qvar)
        close(io)
        return qvar
    else
        return nothing
    end
end