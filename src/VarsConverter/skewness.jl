"""
skewness (f, ff, fff; is_normalize=true)
computes the cross skewness ⟨f'f'f'⟩/⟨f'f'⟩^{1.5} (if is_normalize is true)
"""
function skewness(f, ff, fff; is_normalize=true)
    if is_normalize
        variance = ff .- f.*f
        d = variance.*sqrt.(abs.(variance))
        d.+=eps()
    else
        d=1.0
    end
    return (fff .- 3.0.*ff.*f .+ 2.0.*f.*f.*f)./d
end

"""
skewness (f, g, fg, ff, ffg)
computes the cross skewness ⟨f'f'g'⟩
"""
function skewness(f, g, fg, ff, fgh)
    return fgh .- 2.0.*f.*fg .- g.*ff .+ 2.0.*f.*g.*f
end

"""
skewness (f, g, h, fg, gh, hf, fgh)
computes the cross skewness ⟨f'g'h'⟩
"""
function skewness(f, g, h, fg, gh, hf, fgh)
    return fgh .- f.*gh .-g.*hf .- h.*fg .+ 2.0.*f.*g.*h
end

"""
skewness (r, rf, rff, rfff; is_normalize=true)
computes the Favre-averaged skewness ⟨rf'f'f'⟩/⟨rf'f'⟩^{1.5} (if is_normalize is true).
"""
function skewness(r, rf, rff, rfff; is_normalize=true)

    f_fave = rf./r

    if is_normalize
        variance = (rff .- rf.*rf)./r # ⟨rf''f''⟩/⟨r⟩ = {f''f''}
        d = variance.*sqrt.(abs.(variance))
        d .+=eps()
        return (rfff .- 3.0.*rff.*f_fave .+ 3.0.*rf.*f_fave .- r.*f_fave.*f_fave.*f_fave)./r./d #{f''f''f''}/{f''f''}^{1.5}
    else
        return (rfff .- 3.0.*rff.*f_fave .+ 3.0.*rf.*f_fave .- r.*f_fave.*f_fave.*f_fave)./r
    end
    return (fff .- 3.0.*ff.*f + 2.0.*f.*f.*f)./d
end

"""
skewness (r, rf, rff, rfff; is_normalize=true)
computes the Favre-averaged skewness {f''f''g''} = ⟨rf'f'g'⟩/⟨r⟩
"""
function skewness(r, rf, rg, rff, rfg, rffg)

    rinv = 1.0./r
    f_fave = rf.*rinv
    g_fave = rg.*rinv
    v1 = rinv.*rffg
    v2 = rinv.*rfg.*f_fave
    v3 = rinv.*rff.*g_fave
    v4 = rinv.*rf.*f_fave.*g_fave

    return v1 .- 2.0.*v2 .- v3 .+2.0.*v4
end
