function ∫fdy(f,y; intrange, type="trapezoid_extp")
    if type=="trapezoid"
        return ∫fdy_tp(f,y; intrange)
    elseif type=="trapezoid_extp"
        return ∫fdy_tpex(f,y; intrange)
    end
end

"""
∫fdy_tp(f,y; intrange)
computes the integral function value `f` within `intrange` with trapezoid approximation.
"""
function ∫fdy_tp(f,y; intrange)
    zs   = filter(y-> intrange[1]<=y<=intrange[2], y)
    if length(zs)<2
        @warn "length(zs<2). 0 is returned."
        return 0
    else
        dz   = diff(zs)
        data = f[findall(y-> intrange[1]<=y<=intrange[2], y)]
        sumv = 0.0
        for i in eachindex(dz)
            fm = data[i]
            fp = data[i+1]
            sumv += dz[i]*(fp+fm)/2
        end
        return sumv
    end
end


"""
∫fdy_tpex(f,y; intrange)
computes the integral function value `f` within `intrange` with trapezoid approximation, padding the loss the gap between intrange[2] and y[le] where le is the greatest integer satisfying y[le]<=intrange[2] using a linear intporation.
"""
function ∫fdy_tpex(f,y; intrange)
    zs   = filter(y-> intrange[1]<=y<=intrange[2], y)
    if length(zs)<2
        @warn "length(zs<2). 0 is returned."
        return 0
    else
        dz   = diff(zs)
        data = f[findall(y-> intrange[1]<=y<=intrange[2], y)]
        sumv = 0.0
        for i in eachindex(dz)
            fm = data[i]
            fp = data[i+1]
            sumv += dz[i]*(fp+fm)/2
        end
        le = findfirst(j -> intrange[2] < y[j+1], 1:length(y)) # y[le]<=intrange[2]
        if le<length(y)
            extp = (intrange[2]-y[le]) / (y[le+1]-y[le])* (y[le+1]-y[le]) * 0.5*(f[le]+f[le+1])
            return sumv + extp
        else
            return sumv
        end
    end
end