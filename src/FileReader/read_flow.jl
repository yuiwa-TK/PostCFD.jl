"""
read a file with pl3d format written in *little-endian* & *stream*.  
This function returns *qall=Array{Float64}(jmax,kmax,lmax,5)*.
"""
function read_flow_double(filename::String)
    @show filename
    # settings ==============================================
    dims = Array{Int32}(undef,(3))
    nov  = 5

    dims    = Array{Int32}(undef,(3))
    params  = Array{Float64}(undef,(4))
    qall    = 0.0

    open(filename,"r") do io 
        @show read!(io,ltoh(dims))
        @show read!(io,ltoh(params))
        jmax = dims[1]
        kmax = dims[2]
        lmax = dims[3]
        qall = Array{Float64}(undef,(jmax,kmax,lmax,nov))
        read!(io,ltoh(qall))
    end
    return qall
end

"""
read a file with pl3d format written in *little-endian* & *stream*.  
This function returns *qall=Array{Float32}(jmax,kmax,lmax,5)*.
"""
function read_flow_single(filename::String)
    @show filename
    # settings ==============================================
    dims = Array{Int32}(undef,(3))
    nov  = 5

    dims    = Array{Int32}(undef,(3))
    params  = Array{Float32}(undef,(4))
    qall    = 0.0
   
    open(filename,"r") do io 
        @show read!(io,dims)
        @show read!(io,params)
        jmax = dims[1]
        kmax = dims[2]
        lmax = dims[3]
        qall = Array{Float32}(undef,(jmax,kmax,lmax,nov))
        read!(io,qall)
    end
    return qall
end
read_flow_fv = read_flow_single

"""
read a file with "reatart" format written in *little-endian" & "stream*.  
this function returns *qall=Array{Float64}(jmax,kmax,lmax,5)*.
"""
function read_restart(filename)
    @show filename
    # settings ==============================================
    dims= Array{Int32}(undef,(3))
    nov  = 5

    dims    = Array{Int32}(undef,(3))
    params  = Array{Float64}(undef,(3))
    nc      = Array{Int32}(undef,1)
    qall    = 0.0
   
    open(filename,"r") do io 
        @show read!(io,dims)
        @show read!(io,params)
        @show read!(io,nc)
        
        jmax = dims[1]
        kmax = dims[2]
        lmax = dims[3]
        qall = Array{Float64}(undef,(jmax,kmax,lmax,nov))
        read!(io,qall)
    end
    return qall
end


"""
    dim,param = read_flow_params(filename::AbstractString)

returns the dimension and parameter of flowfile written in pl3d format.
"""
function read_flow_params(filename::AbstractString; precision)
    println("Reading flow params")
    @show filename

    if precision == "single"
        dims    = Array{Int32}(undef,(3))
        params  = Array{Float32}(undef,(4))
    
        open(filename,"r") do io
            read!(io,dims)
            read!(io,params)
        end
        return dims,Float32.(params)

    elseif precision=="double"
        dims    = Array{Int32}(undef,(3))
        params  = Array{Float64}(undef,(3))
        nc      = Array{Int32}(undef,(1))
        open(filename,"r") do io
            read!(io,dims)
            read!(io,params)
            read!(io,nc)
            params=append!(params, nc)
        end
        return dims, Float32.(params)
    end
end