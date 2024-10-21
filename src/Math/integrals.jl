function ∫fdy(f,y; intrange, type="trapezoid")
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
