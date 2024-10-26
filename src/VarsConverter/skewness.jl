"""
skewness (f, ff, fff; is_normalize=true)
computes the cross skewness ⟨f'f'f'⟩/⟨f'f'⟩^{1.5} (if is_normalize is true)
"""
function skewness(f, ff, fff; is_normalize=true)
    if is_normalize
        variance = ff .- f.*f
        d = variance.*sqrt.(abs.(variance))
    else
        d=1.0
    end
    return (fff .- 3.0.*ff.*f + 2.0.*f.*f.*f)./d
end

"""
skewness (f, g, h, fg, gh, hf, fgh)
computes the cross skewness ⟨f'g'h'⟩
"""
function skewness(f, g, h, fg, gh, hf, fgh)
    return fgh .- f.*gh .-g.*hf .- h.*fg + 2.0.*f.*g.*h
end