function sutherland(T::S,T0::S=1.0,μ0::S=1.0, C::S=117) where S
    return μ0*(T/T0)^(1.5)*(T0+C)/(T+C)
end
function sutherland(T::Vector{S},T0::S=1.0,μ0::S=1.0, C::S=117) where S<:AbstractFloat
    return μ0.*(T./T0).^(1.5)*(T0+C)./(T.+C)
end