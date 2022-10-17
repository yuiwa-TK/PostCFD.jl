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
    @show filename,nvar
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