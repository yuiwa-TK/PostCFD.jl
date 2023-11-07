function read_function_single(filename::String; verbose=2)
    if verbose<=1
        @show filename
    end
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

function read_function_double(filename::String; verbose=2)
    if verbose<=1
        @show filename
    end
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

function read_function_dims(filename::String; verbose=2)
    dims = Array{Int32}(undef,(4))
    open(filename,"r") do io 
        read!(io,dims)
    end
    if verbose>=2
        @info dims
    end
    return dims
end


function typeof_functionfile(filename::AbstractString; verbose=2)
    Nb_INT32 = 4
    NBF_FLOAT64 = 8
    NBF_FLOAT32 = 4
    Nvars = prod(read_function_dims(filename; verbose=verbose))
    Nb_file = filesize(filename)

    if Nb_file == 4*Nb_INT32 + Nvars*NBF_FLOAT32
        if verbose>=1
            println("$filename is single format")
        end
        return "single"
    elseif Nb_file == 4*Nb_INT32 + Nvars*NBF_FLOAT64
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
    read_flow_auto(filename::AbstractString)
automatically determines the file type written in pl3d format.
"""
function read_function_auto(filename::AbstractString; verbose=2)
    Nb_INT32 = 4
    NBF_FLOAT64 = 8
    NBF_FLOAT32 = 4
    Nvars = prod(read_function_dims(filename;verbose=0))
    Nb_file = filesize(filename)

    if Nb_file == 4*Nb_INT32 + Nvars*NBF_FLOAT32
        return read_function_single(filename; verbose=verbose)
    elseif Nb_file == 4*Nb_INT32 + Nvars*NBF_FLOAT64
        return read_function_double(filename; verbose=verbose)
    else
        @error println("$filename is not written in pl3d format or written with record marker .")
        return NaN
    end
end
read_function=read_function_auto


function read_function_specifyingvaribale(filename::String,idvar::Int; verbose=2)
    NBF_FLOAT64 = 8
    NBF_FLOAT32 = 4
    if verbose>=1
        @show filename
    end
    tp = typeof_functionfile(filename)
    dims = Array{Int32}(undef,(4))

    if tp=="single"
        io   = open(filename,"r") 
        read!(io,dims)
        qvar = Array{Float32}(undef,(dims[1],dims[2],dims[3]))
        Nb_skip = prod(@view dims[1:3])*(idvar-1)*NBF_FLOAT32
        skip(io,Nb_skip)
        read!(io,qvar)
        close(io)
        return qvar
    elseif tp=="double"
        io   = open(filename,"r") 
        read!(io,dims)
        qvar = Array{Float64}(undef,(dims[1],dims[2],dims[3]))
        Nb_skip = prod(@view dims[1:3])*(idvar-1)*NBF_FLOAT64
        skip(io,Nb_skip)
        read!(io,qvar)
        close(io)
        return qvar
    else
        return nothing
    end
end

"""
    fuction read_function_specifying_l_and_variable(filename::String,lid::Int, idvar::Int; verbose=2)
reads the value of specified veriable and at l index.

## i.e.)
    read_function_specifying_l_and_variable(file, 10, 1) == read_function_auto(file)[:,:,10,1]
"""
function read_function_specifying_l_and_variable(filename::String,lid::Int, idvar::Int; verbose=2)
    NBF_FLOAT64 = 8
    NBF_FLOAT32 = 4
    if verbose>=1
        @show filename
    end
    tp = typeof_functionfile(filename)
    dims = Array{Int32}(undef,(4))

    if tp=="single"
        io   = open(filename,"r") 
        read!(io,dims)
        qvar = Array{Float32}(undef,(dims[1],dims[2],1))
        
        # skip bytes for variable
        Nb_skip = prod(@view dims[1:3])*(idvar-1)*NBF_FLOAT32
        skip(io,Nb_skip)

        # skip bytes for ldim
        Nb_skip2 = prod(@view dims[1:2])*(lid-1)*NBF_FLOAT32
        skip(io,Nb_skip2)
        read!(io,qvar)
        close(io)
        return qvar
    elseif tp=="double"
        io   = open(filename,"r") 
        read!(io,dims)
        qvar = Array{Float64}(undef,(dims[1],dims[2],dims[3]))

        # skip bytes for variable
        Nb_skip = prod(@view dims[1:3])*(idvar-1)*NBF_FLOAT64
        skip(io,Nb_skip)
      
        # skip bytes for ldim
        Nb_skip2 = prod(@view dims[1:2])*(lid-1)*NBF_FLOAT64
        skip(io,Nb_skip2)
        read!(io,qvar)
        close(io)
        return qvar
    else
        return nothing
    end
end
