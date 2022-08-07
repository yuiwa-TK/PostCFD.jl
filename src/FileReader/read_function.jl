function read_function_single(filename::String,nvar::Int)
    @show filename
    # settings ==============================================
    dims = Array{Int32}(undef,(4))
    params  = Array{Float32}(undef,(4))
    qall    = 0.0

    open(filename,"r") do io 
        @show read!(io,ltoh(dims))
        jmax = dims[1]
        kmax = dims[2]
        lmax = dims[3]
        qall = Array{Float32}(undef,(jmax,kmax,lmax,nvar))
        read!(io,ltoh(qall))
    end
    return qall
end