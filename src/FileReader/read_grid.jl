using Printf
"""
"DOUBE PRECISION GRID"
read a file with xyz format written in "little-endian" & "stream".  
This function returns xyz=Array{Float64}(jmax,kmax,lmax,3).
"""
function read_grid_double(filename::AbstractString)
    @show filename
    xyz = 0
    dims= Array{Int32}(undef,(3))
    open(filename,"r") do io 
        read!(io,dims)
        @show dims
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
function read_grid_single(filename::AbstractString)
    @show filename
    xyz = 0
    dims= Array{Int32}(undef,(3))
    open(filename,"r") do io 
        read!(io,dims)
        @show dims
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
function read_grid_dims(filename::AbstractString)
    @show filename
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
function read_grid_auto(filename::AbstractString)
    Nb_INT32 = 32
    NBF_FLOAT64 = 64
    NBF_FLOAT32 = 32

    Npoints = prod(read_grid_dims(filename))
    Nb_file = filesize(filename)

    if Nb_file == 3*Nb_INT32 + Npoints*NBF_FLOAT32
        return read_grid_single(filename)
    elseif Nb_file == 3*Nb_INT32 + Npoints*NBF_FLOAT64
        return read_grid_double(filename)
    else
        @error println("$filename is not written in pl3d format.")
        return NaN
    end
end
read_grid=read_grid_auto