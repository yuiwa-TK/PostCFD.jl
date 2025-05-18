
"""
     derivative_2ndcentral(f::AbstractVector{T},x::AbstractVector{S}) where {S,T}
returns df/dx
"""
function derivative_2ndcentral(f::AbstractVector{T}, x::AbstractVector{S}) where {S,T}
    dxdxi = similar(x)
    nmax = length(dxdxi)
    for n in 2:nmax-1
        dxdxi[n] = (x[n+1] - x[n-1]) / 2.0
    end
    dxdxi[1] = (-3x[1] + 4x[2] - x[3]) / (2.0)
    dxdxi[end] = (3x[end] - 4x[end-1] + x[end-2]) / (2.0)

    dfdx = similar(f)
    @inbounds for n in 2:nmax-1
        dfdx[n] = (f[n+1] - f[n-1]) / (2.0 * dxdxi[n])
    end
    dfdx[1] = (-3f[1] + 4f[2] - f[3]) / (2.0 * dxdxi[1])
    dfdx[end] = (3f[end] - 4f[end-1] + f[end-2]) / (2.0 * dxdxi[end])
    return dfdx
end

"""
   derivative_1stsided(f::AbstractVector{T},x::AbstractVector{S}) where {S,T}
returns df/dx
"""
function derivative_1stsided(f::AbstractVector{T}, x::AbstractVector{S}) where {S,T}
    dxdxi = similar(x)
    nmax = length(dxdxi)
    for n in 1:nmax-1
        dxdxi[n] = (x[n+1] - x[n])
    end
    dxdxi[end] = (3x[end] - 4x[end-1] + x[end-2]) / (2.0)

    dfdx = similar(f)
    @inbounds for n in 1:nmax-1
        dfdx[n] = (f[n+1] - f[n]) / (dxdxi[n])
    end
    dfdx[end] = (3f[end] - 4f[end-1] + f[end-2]) / (2.0 * dxdxi[end])
    return dfdx
end

"""
    derivative_2ndcentral(f::AbstractVector{S}) where S
"""
function derivative_2ndcentral(f::AbstractVector{S}) where S
    nmax = length(f)
    df = similar(f)

    @inbounds for j = 2:nmax-1
        df[j] = 0.5 * ((f[j+1] - f[j]) + (f[j] - f[j-1]))
    end
    df[1] = (-3f[1] + 4f[2] - f[3]) * 0.5
    df[end] = (3f[end] - 4f[end-1] + f[end-2]) * 0.5
    return df
end


"""
    derivative_compact_6th(g::AbstractVector,x::AbstractVector{S}) computes dg/dx
"""
function derivative_compact_6th(g::AbstractVector{T}, x::AbstractVector{S}) where {S,T}
    dxdxi = derivative_compact_6th(x)
    dgdxi = derivative_compact_6th(g)
    dfdx = dgdxi ./ dxdxi
    return dfdx
end


"""
    derivative_compact_6th(f::AbstractVector{S}) where S
"""
function derivative_compact_6th(f::AbstractVector{S}) where S
    aa, bb, cc = 1 / 3, 1.0 , 1 / 3
    div1, coe1, coe2 = 1 / 36, 1.0, 28
    div2, coe21, coe22 = 1.0, 0.5, 0.5

    nmax = length(f)
    NS, NE = 3, nmax - 2
    df = similar(f)

    # rhs
    df[1] = (-3f[1] + 4f[2] - f[3]) * 0.5
    df[end] = (3f[end] - 4f[end-1] + f[end-2]) * 0.5
    df[NS-1] = (coe21 * (f[2] - f[1]) + coe22 * (f[3] - f[2])) * div2
    df[NE+1] = -(coe21 * (f[nmax-1] - f[nmax]) + coe22 * (f[nmax-2] - f[nmax-1])) * div2
    @inbounds for j = NS:NE
        df[j] = (coe1 * (f[j+2] - f[j-2]) + coe2 * (f[j+1] - f[j-1])) * div1
    end
    df[NS] = df[NS] - aa * df[NS-1]
    df[NE] = df[NE] - cc * df[NE+1]
    df[NS:NE] = kernel_tridiagonal(aa, bb, cc, NS, NE, df)
    return df
end

"""
    function kernel_tridiagonal(a,b,c,ns,ne,rhs) returns f
    
## example
    --             --
    | b c         | | f[ns]   |       | rhs[ns] | 
    | a b c       | | f[ns+1] |       |         |
    |     :       | |    :    |  =    |         |
    |     a b c   | | f[ne-1] |       |         |
    |       a b c | | f[ne]   |       | rhs[ne] | 
"""
function kernel_tridiagonal(a::S, b::S, c::S, ns::T, ne::T, rhs::AbstractVector{U}) where {S,T,U}
    gamma = similar(rhs) #upper
    rsvec = rhs
    # removing lower part (forward sweep)
    rsvec[ns] = rsvec[ns] / b
    gamma[ns] = c / b
    @inbounds for n in ns+1:ne
        beta_inv = 1 / (b - a * gamma[n-1])
        rsvec[n] = (rsvec[n] - a * rsvec[n-1]) * beta_inv #update rsvec
        gamma[n] = c * beta_inv
    end
    # removing upper part (backward sweep)
    @inbounds for n in ne-1:-1:ns
        rsvec[n] = rsvec[n] - gamma[n] * rsvec[n+1]
    end
    return rsvec[ns:ne]
end