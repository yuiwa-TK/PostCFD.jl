# 3D flow
function conv2prim(qc,dim::Tuple{Int,Int,Int,Int})
  @assert dim == size(qc)
  println("3D flow")
  jmax,kmax,lmax,nvar = dim

  gamma  = 1.4
  gamma1 = 0.4
  gasc   = 1/gamma

  r = qc[:,:,:,1]
  ru= qc[:,:,:,2]
  rv= qc[:,:,:,3]
  rw= qc[:,:,:,4]
  re= qc[:,:,:,5]

  qp= Array{AbstractFloat}(undef,(jmax,kmax,lmax,6))
  u = ru./r
  v = rv./r
  w = rw./r
  ke= 0.5.*r.*(u.*u + v.*v + w.*w)
  p = gamma1.*(re-ke)
  t = p./(r*gasc)

  qp[:,:,:,1] = r 
  qp[:,:,:,2] = u 
  qp[:,:,:,3] = v 
  qp[:,:,:,4] = w 
  qp[:,:,:,5] = t
  qp[:,:,:,6] = p
  return qp
end

# 2D flow
function conv2prim(qc,dim::Tuple{Int,Int,Int})
  @assert dim == size(qc)
  println("2D flow")
  jmax,lmax,nvar = dim

  gamma  = 1.4
  gamma1 =0.4
  gasc   = 1/gamma

  r = qc[:,:,1]
  ru= qc[:,:,2]
  rv= qc[:,:,3]
  rw= qc[:,:,4]
  re= qc[:,:,5]

  qp= Array{AbstractFloat}(undef,(jmax,lmax,6))
  u = ru./r
  v = rv./r
  w = rw./r
  ke= 0.5.*r.*(u.*u + v.*v + w.*w)
  p = gamma1.*(re-ke)
  t = p./(r*gasc)

  qp[:,:,1] = r 
  qp[:,:,2] = u 
  qp[:,:,3] = v 
  qp[:,:,4] = w 
  qp[:,:,5] = t
  qp[:,:,6] = p
  return qp
end

# 1D flow
function conv2prim(qc,dim::Tuple{Int,Int})
  @assert dim == size(qc)
  println("1D flow")
  jmax,nvar = dim

  gamma  = 1.4
  gamma1 = 0.4
  gasc   = 1.0/gamma

  r = qc[:,1]
  ru= qc[:,2]
  rv= qc[:,3]
  rw= qc[:,4]
  re= qc[:,5]

  qp= Array{AbstractFloat}(undef,(jmax,6))
  u = ru./r
  v = rv./r
  w = rw./r
  ke= 0.5*r.*(u.*u + v.*v + w.*w)
  p = gamma1.*(re-ke)
  t = p./(r*gasc)

  qp[:,1] = r 
  qp[:,2] = u 
  qp[:,3] = v 
  qp[:,4] = w 
  qp[:,5] = t
  qp[:,6] = p
  return qp
end
