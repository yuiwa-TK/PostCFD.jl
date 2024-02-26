module FortranFileReader
using FortranFiles
# URL for the document (https://traktofon.github.io/FortranFiles.jl/stable/index.html)
# using DelimitedFiles
# https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/


"""
    fread.Grid(filename::String,myEndian::String,myMarker::String,myPrecision::String)

    -Read grid files (plot3d) written by Fortran code shown below
    --integer*4 :: i,j,k,imax,jmax,kmax
    --open(10,file='bin_double_big.dat',form="unformatted")
    --write(10)imax,jmax,kmax
    --write(10)(((q(i,j,k),i=1,imax),j=1,jmax),k=1,kmax)
    --close(10)
    -input variables
    --filename (T:String)   : "the name of file"
    --myEndian (T:String)   : "little-endian" or "big-endian"
    --myMarker (T: ???)     : RECMRK4B or  RECMRK8B or "stream"
    --myPrecision(T:String) : "single" or "double"
"""
function grid(filename::String; myEndian::String,myMarker::String,myPrecision::String)
    @show filename
    if myEndian == "big" || myEndian =="b" || myEndian =="big-endian"
        myEndian = "big-endian"
    elseif myEndian == "little" || myEndian =="l" ||  myEndian =="little-endian"
        myEndian = "little-endian"
    else
        println("input error : myEndian")
    end
    if myPrecision == "single" || myPrecision =="s"
        myPrecision = Float32
    elseif myPrecision == "double" || myPrecision =="d"
        myPrecision = Float64
    else
        println("input error : myPrecision")
    end

    if myMarker == "RECMRK4B" || myMarker == "4"
        myMarker = RECMRK4B
    elseif myMarker == "RECMRK8B" || myMarker == "8"
        myMarker = RECMRK8B
    elseif myMarker == "stream" || myMarker == "s"
        println("not supported")
    else
        println("not supported")
    end
    # @show myEndian,myMarker,myPrecision
    io    = FortranFile(filename,convert=myEndian,marker=myMarker)
    sizes = read(io,(Int32,3))
    n     = Tuple(sizes)
    x,y,z = read(io,(myPrecision,n),(myPrecision,n),(myPrecision,n))

    return sizes, x,y,z
end


"""
    flow(filename::String,myEndian::String,myMarker::String,myPrecision::String,flag::Bool)

    returns sizes, info,q

-Read grid files (plot3d) written by Fortran code shown below
--integer*4 :: i,j,k,imax,jmax,kmax
--open(10,file=file_field,form="unformatted",access="myMarker")
--read(10)imax,jmax,kmax
--read(10)fvmach,alpha,re,time
--read(10)((((q_read(i,j,k,kv),i=1,imax),j=1,jmax),k=1,kmax),kv=1,5)
--close(10)

-input variables
--filename (T:String)   : "the name of file"
--myEndian (T:String)   : "little-endian" or "big-endian"
--myMarker (T: ???)     : RECMRK4B or  RECMRK8B or "stream"
--myPrecision(T:String) : "single" or "double"
"""
function flow(filename::String; myEndian::String,myMarker::String,myPrecision::String)
    @show filename
    if myEndian == "big" || myEndian =="b" || myEndian =="big-endian"
        myEndian = "big-endian"
    elseif myEndian == "little" || myEndian =="l" ||  myEndian =="little-endian"
        myEndian = "little-endian"
    else
        println("input error : myEndian")
    end
    if myPrecision == "single" || myPrecision =="s"
        myPrecision = Float32
    elseif myPrecision == "double" || myPrecision =="d"
        myPrecision = Float64
    else
        println("input error : myPrecision")
    end

    if myMarker == "RECMRK4B" || myMarker == "4"
        myMarker = RECMRK4B
    elseif myMarker == "RECMRK8B" || myMarker == "8"
        myMarker = RECMRK8B
    elseif myMarker == "stream" || myMarker == "s"
        println("not supported")
    else
        println("not supported")
    end

    # @show myEndian,myMarker,myPrecision,flag
    io    = FortranFile(filename,convert=myEndian,marker=myMarker)
    sizes = read(io,(Int32,3))                       #imax,jmax,kmax
    n     = Tuple(push!(sizes,5))
    q     = read(io,(myPrecision,n))                 #q
    return sizes, info,q
end



"""
    func(filename::String,myEndian::String,myMarker::String,myPrecision::String)

    returns sizes,q

-Read grid files (plot3d) written by Fortran code shown below
--integer*4 :: i,j,k,imax,jmax,kmax
--open(10,file=file_field,form="unformatted")
--read(10)imax,jmax,kmax,nvar
--read(10)((((q_read(i,j,k,kv),i=1,imax),j=1,jmax),k=1,kmax),kv=1,nov)
--close(10)
-input variables
--filename (T:String)   : "the name of file"
--myEndian (T:String)   : "little-endian" or "big-endian"
--myMarker (T: ???)     : RECMRK4B or  RECMRK8B or "stream"
--myPrecision(T:String) : "single" or "double"
"""
function func(filename::String,myEndian::String,myMarker::String,myPrecision::String)
    @show filename
    if myEndian == "big" || myEndian =="b"
        myEndian = "big-endian"
    elseif myEndian == "little" || myEndian =="l"
        myEndian = "little-endian"
    else
        println("input error : myEndian")
    end
    if myPrecision == "single" || myPrecision =="s"
        myPrecision = Float32
    elseif myPrecision == "double" || myPrecision =="d"
        myPrecision = Float64
    else
        println("input error : myPrecision")
    end

    if myMarker == "RECMRK4B" || myMarker == "4"
        myMarker = RECMRK4B
    elseif myMarker == "RECMRK8B" || myMarker == "8"
        myMarker = RECMRK8B
    elseif myMarker == "stream" || myMarker == "s"
        println("not supported")
    else
        println("not supported")
    end

    # @show myEndian,myMarker,myPrecision,flag
    io    = FortranFile(filename,convert=myEndian,marker=myMarker)
    sizes = read(io,(Int32,4))                       #imax,jmax,kmax,nov
    @show sizes
    q     = read(io,(myPrecision,Tuple(sizes)))                 #q

    return sizes, q

end

end #module

# HOW TO USE --------------------------------------------------------------------
# filename = "./FortranFiles/bin_grid_single_big.dat"
# sizes,x,y,z = readGrid(filename, "big","RECMRK4B","single")
# filename = "./FortranFiles/bin_q_single_big.dat"
# sizes,info,q = readQ(filename,"b","4","s",5)
#------------------------------------------------------------------------------

module FortranFileWriter
using FortranFiles

"""
    write_grid(filename::String,x::Real,y::Real,z::Real,myEndian::String,myMarker::String,myPrecision::String)
"""
function write_grid(filename::String,x::Real,y::Real,z::Real,myEndian::String,myMarker::String,myPrecision::String)

    filename = filename*"_"*myEndian*"_"*myMarker*"_"*myPrecision*".xyz"
    io = FortranFile(filename,convert=myEndian,marker=myMarker)
    size = size(x)
    write(io,size[1],size[2],size[3])
    write(io,x,y,z)
    
end

"""
    write_flow(filename::String,q::Real,myEndian::String,myMarker::String,myPrecision::String)
"""
function write_flow(filename::String,q::Real,myEndian::String,myMarker::String,myPrecision::String)
    #=Write grid files (plot3d) written by Fortran code shown below
     ++++ integer*4 :: i,j,k,imax,jmax,kmax
     ++++ open(10,file=file_field,form="unformatted",access="stream")
     ++++ write(10)imax,jmax,kmax
     ++++ Write(10)fvmach,alpha,re,time
     ++++ Write(10)((((q_Write(i,j,k,kv),i=1,imax),j=1,jmax),k=1,kmax),kv=1,nov)
     ++++ close(10)
    =#
    size = size(q) # size = (imax,jmax,kmax,nov)
    filename = filename*"_"*myEndian*"_"*maMarker*"_"*myPrecision*"pl3d"
    io = FortranFile(filename,convert=myEndian,marker=myMarker)
    write(size[1],size[2],size[3])
    write(io,q)

end

end # module