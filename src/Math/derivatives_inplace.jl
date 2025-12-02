
function derivative_compact_6th!(f::AbstractVector{<:AbstractFloat},df::AbstractVector{<:AbstractFloat},workm::AbstractVector{<:AbstractFloat})
    aa, bb, cc = 1 / 3, 1.0 , 1 / 3
    div1, coe1, coe2 = 1 / 36, 1.0, 28
    div2, coe21, coe22 = 1.0, 0.5, 0.5

    nmax = length(f)
    NS, NE = 3, nmax - 2

    # rhs
    df[1] = (-3f[1] + 4f[2] - f[3]) * 0.5
    df[end] = (3f[end] - 4f[end-1] + f[end-2]) * 0.5
    df[NS-1] =  (coe21 * (f[2] - f[1]) + coe22 * (f[3] - f[2])) * div2
    df[NE+1] = -(coe21 * (f[nmax-1] - f[nmax]) + coe22 * (f[nmax-2] - f[nmax-1])) * div2
    @inbounds for j = NS:NE
        df[j] = (coe1 * (f[j+2] - f[j-2]) + coe2 * (f[j+1] - f[j-1])) * div1
    end
    df[NS] = df[NS] - aa * df[NS-1]
    df[NE] = df[NE] - cc * df[NE+1]
    @views begin
    df[NS:NE] .= kernel_tridiagonal!(aa, bb, cc, NS, NE, df,workm)
    end
    return df
end

function kernel_tridiagonal!(a::AbstractFloat, b::AbstractFloat, c::AbstractFloat, ns::Integer, ne::Integer, rsvec::AbstractVector{<:AbstractFloat}, workm::AbstractVector{<:AbstractFloat})
    # fill!(workm,0.0)
    # removing lower part (forward sweep)
    rsvec[ns] = rsvec[ns] / b
    workm[ns] = c / b
    @inbounds @simd for n in ns+1:ne
        # beta_inv = 1 / (b - a * workm[n-1])
        rsvec[n] = (rsvec[n] - a * rsvec[n-1]) / (b - a * workm[n-1]) #update rsvec
        workm[n] = c / (b - a * workm[n-1])
    end
    # removing upper part (backward sweep)
    @inbounds @simd for n in ne-1:-1:ns
        rsvec[n] = rsvec[n] - workm[n] * rsvec[n+1]
    end
    return rsvec[ns:ne]
end


"""
    derivative_2ndcentral!(f::AbstractVector{<:AbstractFloat},df::AbstractVector{<:AbstractFloat})
"""
function derivative_2ndcentral!(f::AbstractVector{<:AbstractFloat},df::AbstractVector{<:AbstractFloat}) 
    @assert length(f)==length(df)
    nmax = length(f)
    c11, c12 = 1.5, -0.5

    # @inbounds for j = 2:nmax-1
    for j = 2:nmax-1
        df[j] = 0.5 * ((f[j+1] - f[j]) + (f[j] - f[j-1]))
    end
    # df[1] = (-3f[1] + 4f[2] - f[3]) * 0.5
    # df[end] = (3f[end] - 4f[end-1] + f[end-2]) * 0.5

    df[1] =   c11 * (f[2] - f[1]) + c12 * (f[3] - f[2])
    df[end] = c11 * (- f[end-1] + f[end]) + c12 * (-f[end-2] + f[end-1])
    return df
end


"""
    derivative_4thcentral!(f::AbstractVector{<:AbstractFloat},df::AbstractVector{<:AbstractFloat})
"""
function derivative_4thcentral!(f::AbstractVector{<:AbstractFloat}, df::AbstractVector{<:AbstractFloat})
    @assert length(f)==length(df)
    nmax = length(f)

    # interior: 4th-order central
    @inbounds for j in 3:nmax-2
        df[j] = (-f[j+2] + 8f[j+1] - 8f[j-1] + f[j-2]) / 12
    end

    # boundary: 2nd-order (same style as your code)
    c11, c12 = 1.5, -0.5

    # j = 1, 2 (forward)
    df[1] = c11*(f[2]-f[1]) + c12*(f[3]-f[2])
    df[2] = c11*(f[3]-f[2]) + c12*(f[4]-f[3])

    # j = end, end-1 (backward)
    df[end]     = c11*(-f[end-1]+f[end])     + c12*(-f[end-2]+f[end-1])
    df[end-1]   = c11*(-f[end-2]+f[end-1])   + c12*(-f[end-3]+f[end-2])

    return df
end
