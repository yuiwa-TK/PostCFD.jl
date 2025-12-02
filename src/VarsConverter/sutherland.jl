"""
    sutherland(T::S, T0::S=one(S), μ0::S=one(S), C::S=117*one(S)) where S

Compute the dynamic viscosity using Sutherland's law.

# Arguments
- `T::S`: Temperature (same type S for all arguments)
- `T0::S=one(S)`: Reference temperature (default: 1 of type S)
- `μ0::S=one(S)`: Reference viscosity (default: 1 of type S)
- `C::S=117*one(S)`: Sutherland's constant (default: 117 of type S)

# Returns
- Dynamic viscosity at temperature `T`, with the same type as `T`.

# Example
```julia
μ = sutherland(300.0, 273.15, 1.716e-5, 110.4)
```
"""
function sutherland(T::AbstractFloat, T0=1.0, μ0=1.0, C=117)
    return μ0 * (T / T0)^(1.5) * (T0 + C) / (T + C)
end

function sutherland(T::AbstractArray{<:AbstractFloat},T0=1.0,μ0=1.0, C=117)
    mu = similar(T)
    mu .= μ0.*(T./T0).^(1.5).*(T0+C)./(T.+C)
    return mu
end