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
read_grid = read_grid_double #alias

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
    @show dims
    return dims
end
