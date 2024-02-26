function read_uns(filename::AbstractString; endian="little")
    ngridpoints=Vector{Int32}(undef,1)
    uidum   = Vector{UInt8}(undef,80)
    idum1   = Vector{Int32}(undef,1)
    idum2   = Vector{Int32}(undef,2)
    idum4   = Vector{Int32}(undef,4)
    rdum4   = Vector{Float32}(undef,4)
    nvar=Vector{Int32}(undef,1)
    qall = 0


    open(filename,"r") do io
        read!(io,idum1)
        a=read!(io,uidum)
        String(a)
        read!(io,idum4)
        read!(io,rdum4)
        read!(io,idum1)
        @show read!(io,nvar)
        dumstr = Vector{UInt8}(undef,80*nvar[1]) #namelist of function
        aa=read!(io,dumstr)
        read!(io,idum2)
        @show read!(io,ngridpoints)
        if nvar[1]==0
            # nvar==0ならdumstr(namelist of functions)をskipする
            println("nvar=0 is detected. nvar is reset as 5 (to read q file)")
            nvar[1]=5
        end
        qall = Array{Float32}(undef,(ngridpoints[1],nvar[1]))
        read!(io,qall)
        read!(io,idum1)
    end
    if endian=="little"
        return qall
    else
        return ntoh.(qall)
    end
end


# #IN Fortran 
# open(32, file= trim(File_field), form = "unformatted",access="stream")
# read(32) FV_MAGIC
# read(32) FV_txt
# read(32) dum1, dum2, FV_RESULTS_FILE, dum3
# print*, dum1,dum2,FV_RESULTS_FILE,dum3
# read(32) rdum1,rdum2,rdum3,rdum4
# print*, rdum1,rdum2,rdum3,rdum4
# read(32) ngrid
# read(32) nvar
# allocate(FV_var(nvar))
# do i=1,nvar
#    read(32) FV_var(i)
# enddo
# read(32)FV_NODES, npt
# read(32) FV_VARIABLES
# allocate(qp(nvar,FV_VARIABLES))
# read(32) ((qp(i,n), i=1,nvar),n=1,FV_VARIABLES)
# read(32) FV_BNDRY_VARS
# close(32)