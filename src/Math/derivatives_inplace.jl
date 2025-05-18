
function derivative_compact_6th!(f::AbstractVector{S},df::AbstractVector{T},kernel_tridiagonal!::Function) where {S,T}
    aa, bb, cc = 1 / 3, 1.0 , 1 / 3
    div1, coe1, coe2 = 1 / 36, 1.0, 28
    div2, coe21, coe22 = 1.0, 0.5, 0.5

    nmax = length(f)
    NS, NE = 3, nmax - 2

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
    df[NS:NE] = kernel_tridiagonal!(aa, bb, cc, NS, NE, df)
    return df
end

function kernel_tridiagonal!(a::S, b::S, c::S, ns::T, ne::T, rsvec::AbstractVector{U}, workm::AbstractVector{V}) where {S,T,U,V}
    fill!(workm,0.0)
    # removing lower part (forward sweep)
    rsvec[ns] = rsvec[ns] / b
    workm[ns] = c / b
    @inbounds @simd for n in ns+1:ne
        beta_inv = 1 / (b - a * workm[n-1])
        rsvec[n] = (rsvec[n] - a * rsvec[n-1]) * beta_inv #update rsvec
        workm[n] = c * beta_inv
    end
    # removing upper part (backward sweep)
    @inbounds @simd for n in ne-1:-1:ns
        rsvec[n] = rsvec[n] - workm[n] * rsvec[n+1]
    end
    return rsvec[ns:ne]
end


"""
    derivative_2ndcentral(f::AbstractVector{S},df::AbstractVector{T}) where {S,T}
"""
function derivative_2ndcentral!(f::AbstractVector{S},df::AbstractVector{T}) where {S,T}
    nmax = length(f)

    @inbounds for j = 2:nmax-1
        df[j] = 0.5 * ((f[j+1] - f[j]) + (f[j] - f[j-1]))
    end
    df[1] = (-3f[1] + 4f[2] - f[3]) * 0.5
    df[end] = (3f[end] - 4f[end-1] + f[end-2]) * 0.5
    return df
end