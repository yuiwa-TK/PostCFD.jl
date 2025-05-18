# 3D flow
function prim2conv(qp,dim::Tuple{Int,Int,Int,Int};verbose=0)
    @assert dim == size(qp)
    if verbose>=2
      println("3D flow @ function prim2conv")
    end
    jmax,kmax,lmax,nvar = dim
    gamma1 = 0.4
    qc= Array{AbstractFloat}(undef,(jmax,kmax,lmax,5))
    qc[:,:,:,1] = qp[:,:,:,1]
    qc[:,:,:,2] = qp[:,:,:,1].*qp[:,:,:,2]
    qc[:,:,:,3] = qp[:,:,:,1].*qp[:,:,:,3]
    qc[:,:,:,4] = qp[:,:,:,1].*qp[:,:,:,4] 
    qc[:,:,:,5] = qp[:,:,:,6]/gamma1 + 0.5*qp[:,:,:,1].*(qp[:,:,:,2].*qp[:,:,:,2] + qp[:,:,:,3].*qp[:,:,:,3] + qp[:,:,:,4].*qp[:,:,:,4])
    return qc
end

# 2D flow
function prim2conv(qp,dim::Tuple{Int,Int,Int};verbose=2)
    @assert dim == size(qp)
    if verbose>=2
      println("2D flow @ function prim2conv")
    end
    jmax,lmax,nvar = dim
    gamma1 = 0.4
    qc= Array{AbstractFloat}(undef,(jmax,lmax,5))
    qc[:,:,1] = qp[:,:,1]
    qc[:,:,2] = qp[:,:,1].*qp[:,:,2]
    qc[:,:,3] = qp[:,:,1].*qp[:,:,3]
    qc[:,:,4] = qp[:,:,1].*qp[:,:,4] 
    qc[:,:,5] = qp[:,:,6]/gamma1 + 0.5*qp[:,:,1].*(qp[:,:,2].*qp[:,:,2] + qp[:,:,3].*qp[:,:,3] + qp[:,:,4].*qp[:,:,4])
    return qc
end

# 1D flow
function prim2conv(qp,dim::Tuple{Int,Int};verbose=2)
    @assert dim == size(qp)
    if verbose>=2
      println("1D flow @ function prim2conv")
    end
    lmax,nvar = dim
    gamma1 = 0.4
    qc= Array{AbstractFloat}(undef,(lmax,5))
    qc[:,1] = qp[:,1]
    qc[:,2] = qp[:,1].*qp[:,2]
    qc[:,3] = qp[:,1].*qp[:,3]
    qc[:,4] = qp[:,1].*qp[:,4] 
    qc[:,5] = qp[:,6]/gamma1 + 0.5*qp[:,1].*(qp[:,2].*qp[:,2] + qp[:,3].*qp[:,3] + qp[:,4].*qp[:,4])
    return qc
end
  