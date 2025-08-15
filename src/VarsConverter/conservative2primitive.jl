# 1D flow
function conv2prim(qc,dim::Tuple{Int,Int};verbose=2)
  @assert dim == size(qc)
  if verbose>=2
    println("1D flow")
  end
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

  qp[:,1] .= r 
  qp[:,2] .= u 
  qp[:,3] .= v 
  qp[:,4] .= w 
  qp[:,5] .= t
  qp[:,6] .= p
  return qp
end

# 3D flow
function conv2prim(qc,dim::Tuple{Int,Int,Int,Int};verbose=2)
  @assert dim == size(qc)
  if verbose>=2
    println("3D flow")
  end
  jmax,kmax,lmax,nvar = dim

  gamma  = 1.4
  gamma1 = 0.4
  gasc   = 1/gamma

  @views begin
    r = qc[:,:,:,1]
    ru= qc[:,:,:,2]
    rv= qc[:,:,:,3]
    rw= qc[:,:,:,4]
    re= qc[:,:,:,5]
  end

  # u = ru./r
  # v = rv./r
  # w = rw./r
  # ke= 0.5.*r.*(u.*u + v.*v + w.*w)
  # p = gamma1.*(re-ke)
  # t = p./(r*gasc)

  qp= Array{AbstractFloat}(undef,(jmax,kmax,lmax,6))
  @views qp[:,:,:,1] .= r 
  @views qp[:,:,:,2] .= ru./r 
  @views qp[:,:,:,3] .= rv./r 
  @views qp[:,:,:,4] .= rw./r 
  @views qp[:,:,:,6] .= gamma1.* (re - 0.5.*(ru.*ru + rv.*rv + rw.*rw)./r)
  @views qp[:,:,:,5] .= qp[:,:,:,6]./(r*gasc)
  return qp
end

# 2D flow
function conv2prim(qc,dim::Tuple{Int,Int,Int};verbose=2)
  @assert dim == size(qc)
  if verbose>=2
    println("2D flow")
  end
  jmax,lmax,nvar = dim

  gamma  = 1.4
  gamma1 =0.4
  gasc   = 1/gamma

  @views begin
    r = qc[:,:,1]
    ru= qc[:,:,2]
    rv= qc[:,:,3]
    rw= qc[:,:,4]
    re= qc[:,:,5]
  end

  qp= Array{AbstractFloat}(undef,(jmax,lmax,6))
  u = ru./r
  v = rv./r
  w = rw./r
  ke= 0.5.*r.*(u.*u + v.*v + w.*w)
  p = gamma1.*(re-ke)
  t = p./(r*gasc)

  @views qp[:,:,1] .= r 
  @views qp[:,:,2] .= u 
  @views qp[:,:,3] .= v 
  @views qp[:,:,4] .= w 
  @views qp[:,:,5] .= t
  @views qp[:,:,6] .= p
  return qp
end