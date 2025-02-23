module FortranFileReader 
using FortranFiles
# URL for the document (https://traktofon.github.io/FortranFiles.jl/stable/index.html)
# using DelimitedFiles
# https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/


function assign_fileformat(myEndian,myPrecision,myMarker;verbose=verbose)

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
    
    if verbose>1
        @info myEndian,myMarker,myPrecision
   end

   return myEndian,myPrecision,myMarker

end

"""
    fread.Grid(filename::AbstractString,myEndian::AbstractString,myMarker::AbstractString,myPrecision::AbstractString,verbose::Int=0)

    -Read grid files (plot3d) written by Fortran code shown below
    --integer*4 :: i,j,k,imax,jmax,kmax
    --open(10,file='bin_double_big.dat',form="unformatted")
    --write(10)imax,jmax,kmax
    --write(10)(((q(i,j,k),i=1,imax),j=1,jmax),k=1,kmax)
    --close(10)
    -input variables
    --filename (T:AbstractString)   : "the name of file"
    --myEndian (T:AbstractString)   : "little-endian" or "big-endian"
    --myMarker (T:AbstractString)     : "4" (="RECMRK4B") or "8"(="RECMRK8B") or "stream"
    --myPrecision(T:AbstractString) : "single" or "double"
"""
function grid(filename::AbstractString; myEndian="big"::AbstractString,myMarker="4"::AbstractString,myPrecision="single"::AbstractString,verbose=0)
    if verbose>0
        @show filename
    end
    myEndian,myPrecision,myMarker = assign_fileformat(myEndian,myPrecision,myMarker;verbose=verbose)
    io    = FortranFile(filename,convert=myEndian,marker=myMarker)
    sizes = read(io,(Int32,3))
    if verbose>1
        @info sizes
    end
    n  = Tuple(push!(sizes,3))
    G  = read(io,(myPrecision,n))

    return sizes, G
end


"""
    flow(filename::AbstractString,myEndian::AbstractString,myMarker::AbstractString,myPrecision::AbstractString,verbose::Int=0)

    returns sizes, info,q

-Read grid files (plot3d) written by Fortran code shown below
--integer*4 :: i,j,k,imax,jmax,kmax
--open(10,file=file_field,form="unformatted",access="myMarker")
--read(10)imax,jmax,kmax
--read(10)fvmach,alpha,re,time
--read(10)((((q_read(i,j,k,kv),i=1,imax),j=1,jmax),k=1,kmax),kv=1,5)
--close(10)

-input variables
--filename (T:AbstractString)   : "the name of file"
--myEndian (T:AbstractString)   : "little-endian" or "big-endian"
--myMarker (T:AbstractString)     : "4" (="RECMRK4B") or "8"(="RECMRK8B") or "stream"
--myPrecision(T:AbstractString) : "single" or "double"
"""
function flow(filename::AbstractString; myEndian="big"::AbstractString,myMarker="4"::AbstractString,myPrecision="single"::AbstractString, verbose::Int=0)
    if verbose>0
        @show filename
    end
    myEndian,myPrecision,myMarker = assign_fileformat(myEndian,myPrecision,myMarker;verbose=verbose)

    io    = FortranFile(filename,convert=myEndian,marker=myMarker)
    sizes = read(io,(Int32,3))                       #imax,jmax,kmax
    if verbose>1
        @info sizes
    end
    n     = Tuple(push!(sizes,5))
    info  = read(io,(myPrecision,4))
    if verbose>1
        @info "mach,alpha,re,time", info
    end
    q     = read(io,(myPrecision,n))                 #q
    return sizes, info,q
end

"""
flow_params(filename::AbstractString,myEndian::AbstractString,myMarker::AbstractString,myPrecision::AbstractString,verbose::Int=0)

returns sizes, info

-input variables
--filename (T:AbstractString)   : "the name of file"
--myEndian (T:AbstractString)   : "little-endian" or "big-endian"
--myMarker (T:AbstractString)     : "4" (="RECMRK4B") or "8"(="RECMRK8B") or "stream"
--myPrecision(T:AbstractString) : "single" or "double"
"""
function flow_params(filename::AbstractString; myEndian="big"::AbstractString,myMarker="4"::AbstractString,myPrecision="single"::AbstractString, verbose::Int=0)
    if verbose>0
        @show filename
    end
    myEndian,myPrecision,myMarker = assign_fileformat(myEndian,myPrecision,myMarker;verbose=verbose)

    io    = FortranFile(filename,convert=myEndian,marker=myMarker)
    sizes = read(io,(Int32,3))                       #imax,jmax,kmax
    if verbose>1
        @info sizes
    end
    info  = read(io,(myPrecision,4))
    if verbose>1
        @info "mach,alpha,re,time", info
    end
    return info
end

"""
restart_params(filename::AbstractString,myEndian::AbstractString,myMarker::AbstractString,myPrecision::AbstractString,verbose::Int=0)

returns sizes..., info (i.e., mach, aoa, re, nc)

-input variables
--filename (T:AbstractString)   : "the name of file"
--myEndian (T:AbstractString)   : "little-endian" or "big-endian"
--myMarker (T:AbstractString)     : "4" (="RECMRK4B") or "8"(="RECMRK8B") or "stream"
--myPrecision(T:AbstractString) : "single" or "double"
"""
function restart_params(filename::AbstractString; myEndian="big"::AbstractString,myMarker="4"::AbstractString,myPrecision="single"::AbstractString, verbose::Int=0)
    if verbose>0
        @show filename
    end
    myEndian,myPrecision,myMarker = assign_fileformat(myEndian,myPrecision,myMarker;verbose=verbose)

    io    = FortranFile(filename,convert=myEndian,marker=myMarker)
    sizes = read(io,(Int32,3))                       #imax,jmax,kmax
    if verbose>1
        @info sizes
    end
    info,nstep  = read(io,(myPrecision,3), (Int32,1))
    if verbose>1
        @info "mach,alpha,re,nstep", info,nstep
    end
    return info...,nstep[1]
end

"""
restart(filename::AbstractString,myEndian::AbstractString,myMarker::AbstractString,myPrecision::AbstractString,verbose::Int=0)

returns q

-input variables
--filename (T:AbstractString)   : "the name of file"
--myEndian (T:AbstractString)   : "little-endian" or "big-endian"
--myMarker (T:AbstractString)     : "4" (="RECMRK4B") or "8"(="RECMRK8B") or "stream"
--myPrecision(T:AbstractString) : "single" or "double"
"""
function restart(filename::AbstractString; myEndian="big"::AbstractString,myMarker="4"::AbstractString,myPrecision="single"::AbstractString, verbose::Int=0)
    if verbose>0
        @show filename
    end
    myEndian,myPrecision,myMarker = assign_fileformat(myEndian,myPrecision,myMarker;verbose=verbose)

    io    = FortranFile(filename,convert=myEndian,marker=myMarker)
    sizes = read(io,(Int32,3))                       #imax,jmax,kmax
    if verbose>1
        @info sizes
    end
    
    info,nstep  = read(io,(myPrecision,3), (Int32,1))
    if verbose>1
        @info "mach,alpha,re,nstep", info,nstep
    end
    n = Tuple(push!(sizes,5))
    q = read(io,(myPrecision,n))                 #q
    return q
end


"""
    func(filename::AbstractString,myEndian::AbstractString,myMarker::AbstractString,myPrecision::AbstractString)

    returns sizes,q

-Read grid files (plot3d) written by Fortran code shown below
--integer*4 :: i,j,k,imax,jmax,kmax
--open(10,file=file_field,form="unformatted")
--read(10)imax,jmax,kmax,nvar
--read(10)((((q_read(i,j,k,kv),i=1,imax),j=1,jmax),k=1,kmax),kv=1,nov)
--close(10)
-input variables
--filename (T:AbstractString)   : "the name of file"
--myEndian (T:AbstractString)   : "little-endian" or "big-endian"
--myMarker (T:AbstractString)     : "4" (="RECMRK4B") or "8"(="RECMRK8B") or "stream"
--myPrecision(T:AbstractString) : "single" or "double"
"""
function func(filename::AbstractString,myEndian::AbstractString,myMarker::AbstractString,myPrecision::AbstractString)
    if verbose>0
         @show filename
    end
    myEndian,myPrecision,myMarker = assign_fileformat(myEndian,myPrecision,myMarker;verbose=verbose)

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
    write_grid(filebase::AbstractString,x,y,z,myEndian::AbstractString,myMarker::AbstractString,myPrecision::AbstractString)
"""
function write_grid(filebase::AbstractString,x,y,z,myEndian::AbstractString,myMarker::AbstractString,myPrecision::AbstractString)
    if myEndian == "big" || myEndian =="b"
        myEndian = "big-endian"
    elseif myEndian == "little" || myEndian =="l"
        myEndian = "little-endian"
    else
        println("input error : myEndian")
    end
    filename = filebase*"_"*myEndian*"_"*myMarker*"_"*myPrecision*".xyz"
    io = FortranFile(filename,convert=myEndian,marker=myMarker,"w")
    size = size(x)
    write(io,size[1],size[2],size[3])
    write(io,x,y,z)
end

"""
    write_grid(filebase::AbstractString,G;myEndian="little"::AbstractString,myMarker="4"::AbstractString,myPrecision="single"::AbstractString)
"""
function write_grid(filebase::AbstractString,G;myEndian="little-endian"::AbstractString,myMarker="4"::AbstractString,myPrecision="single"::AbstractString)
    if myEndian == "big" || myEndian =="b"
        myEndian = "big-endian"
    elseif myEndian == "little" || myEndian =="l"
        myEndian = "little-endian"
    else
        println("input error : myEndian")
    end
    @show typeof(G)

    filename = filebase*"_"*myEndian*"_"*myMarker*"_"*myPrecision*".xyz"
    # io = FortranFile(filename,"w",convert=myEndian,marker=myMarker)
    io = FortranFile(filename,"w")
    size = size(G)
    write(io,size[1],size[2],size[3])
    write(io,G)
    
end

"""
    write_flow(filename::AbstractString,q;myEndian="little"::AbstractString,myMarker="4"::AbstractString,myPrecision="single"::AbstractString)
"""
function write_flow(filename::AbstractString,q;myEndian="little"::AbstractString,myMarker="4"::AbstractString,myPrecision="single"::AbstractString)
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