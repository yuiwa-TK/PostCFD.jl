"""
    vorticity(v::Matrix{S},w::Matrix{S},y::Vector{S},z::Vector{S}) where S
compute the vortices in xi direction, assuming in the rectangular grid.
    In this function, 4th-order central derivative is used.
    dwdy - dvdy is returned
"""
function vorticity(v::Matrix{S},w::Matrix{S},y::Vector{S},z::Vector{S}) where S
    kmax = length(y)
    lmax = length(z)
    @assert kmax==size(v,1)
    @assert lmax==size(v,2)
    
    dydet=zeros(S,kmax)
    for k=2:kmax-1
        dydet[k] = 0.5*(y[k+1]-y[k-1])
    end
    dydet[1] = (-3y[1]+4y[2]-y[3])*0.5 
    dydet[end] = (3y[end]-4y[end-1]+y[end-2])*0.5

    dzdzt=zeros(S,lmax)
    for l=2:lmax-1
        dzdzt[l] = 0.5*(z[l+1]-z[l-1])
    end
    dzdzt[1] = (-3z[1]+4z[2]-z[3])*0.5
    dzdzt[end]=(3z[end]-4z[end-1]+z[end-2])*0.5

    dwdy = ones(S,kmax,lmax)
    dvdz = ones(S,kmax,lmax)

    # diff for eta dir
    for l=1:lmax
        for k=3:kmax-2
            dwdy[k,l] = (8(w[k+1,l]-w[k-1,l])-(w[k+2,l]-w[k-2,l]))/12/dydet[k]
        end
        dwdy[1,l] = (-3w[1,l]+4w[2,l]-w[3,l])/ (2.0*dydet[1])
        dwdy[2,l] = 0.5*(w[3,l]-w[1,l])/(dydet[2])
        dwdy[end-1,l] = 0.5*(w[end,l]-w[end-2,l])/(dydet[end-1])
        dwdy[end,l] = (3w[end,l]-4w[end-1,l]+w[end-2,l])/(2.0*dydet[end])
    end
    
    # diff for zeta dir
    for k=1:kmax
        for l=3:lmax-2
            dvdz[k,l] = (8(v[k,l+1]-v[k,l+1]) - (v[k,l+2]-v[k,l+2]))/12/dzdzt[l]
        end
        dvdz[k,1]  = (-3v[k,1]+4v[k,2]-v[k,3])*0.5/(dzdzt[1])
        dvdz[k,2]  = 0.5*(v[k,2+1]-v[k,2-1])/(dzdzt[2])
        dvdz[k,end-1]  = 0.5*(v[k,end]-v[k,end-2])/(dzdzt[end-1])
        dvdz[k,end]= (3v[k,end]-4v[k,end-1]+v[k,end-2])*0.5/(dzdzt[end])
    end
    # return dwdy-dvdz,dwdy,dvdz
    return dwdy-dvdz
end