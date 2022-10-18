"""
    read_dividedrestart(Compdir::AbstractString,flowids::Union{OrdinalRange,Vector{Int}})

read divided restart files.
This function returns xyz_vector=Vector{Array{Float64}(jmax,kmax,lmax,3)}.
Arg1 : name of directory that contains the divided restart files.
Arg2 : flow rank(number) indicated by range (eg. 1:2:19) or array (e.g. [1,3,5,⋯,19]) where you will read.
"""
function read_dividedrestart(Compdir::AbstractString,flowids::Union{OrdinalRange,Vector{Int}})
    # OrdinalRange includes typeof(ns:ne) and typeof(ns:nskip:ne) 
    # flowids=ns:nskip:ne
    
    vec_flow = Vector{Array{AbstractFloat}}(undef,length(flowids))
    
    count = 0
    for n in flowids
        divide_num  = @sprintf "%5.5i"  n
        flowfile    = Compdir*"flow_z"*divide_num
        println("read file ",flowfile)
        count+=1
        vec_flow[count]=read_restart(flowfile)
    end
    return vec_flow
end
function read_dividedrestart(Compdir::AbstractString,flowids::Union{OrdinalRange,Vector{Int}},filebase::AbstractString)
    vec_flow = Vector{Array{AbstractFloat}}(undef,length(flowids))
    count = 0
    for n in flowids
        divide_num  = @sprintf "%5.5i"  n
        flowfile    = Compdir*filebase*divide_num
        println("read file ",flowfile)
        count+=1
        vec_flow[count]=read_restart(flowfile)
    end
    return vec_flow
end

function read_dividedsnapshot(Compdir::AbstractString,flowids::Union{OrdinalRange,Vector{Int}},steps::Union{OrdinalRange,Vector{Int}},
                            filebase::AbstractString,suffix::AbstractString;precision::AbstractString)

    vec_flow =  Vector{Array{AbstractFloat}}(undef,length(flowids)*length(steps))
    if precision=="double"
        my_read = read_flow_double
    elseif precision=="single"
        my_read = read_flow_single
    end
    
    count = 0
    for s     in steps
    for irank in flowids
        @info  s, irank
        count      += 1
        srank       = @sprintf "%5.5i"  irank
        sstep       = @sprintf "%8.8i"  s 
        flowfile    = Compdir*filebase*srank*"_"*sstep*suffix
        vec_flow[count]=my_read(flowfile)
    end
    end
    return vec_flow
end

"""
read divided grid files. (eg. grid files in dir of "/Grid/Work.divide/comp/" in HPC Template. )
This function returns xyz_vector=Vector{Array{Float64}(jmax,kmax,lmax,3)}.
Arg1 : name of directory that contains the divided grid files.
Arg2 : grid rank(number) indicated by range (eg. 1:2:19) or array (e.g. [1,3,5,⋯,19]) where you will read.
"""
function read_dividedgrid(Compdir::AbstractString,gridids::Union{OrdinalRange,Vector{Int}})
    # OrdinalRange includes typeof(ns:ne) and typeof(ns:nskip:ne) 
    # gridids=ns:nskip:ne
    vec_xyz = Vector{Array{AbstractFloat}}(undef,length(gridids))
    count = 0
    for n in gridids
        divide_num  = @sprintf "%5.5i"  n
        gfile = Compdir*"grid."*divide_num
        count+=1
        vec_xyz[count]=read_grid(gfile)
    end
    return vec_xyz
end