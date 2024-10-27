

function flatness(f,ff,fff,ffff; is_normalize=true, eps=1f-8)
    
    if is_normalize
        d = ff .- f.*f
        d .+= eps
    else
        d = 1.0
    end
    return (ffff .- 4.0.*fff.*f .+ 6.0.*ff.*f.*f .- 3.0.*f.*f.*f.*f)./d

end


function flatness(r,rf,rff,rfff,rffff; is_normalize=true, eps=1f-8)
    
    if is_normalize
        variance = (rff .- rf.*rf./r)./r # ⟨rf''f''⟩/⟨r⟩ = {f''f''}
        d = variance.*variance 
        d.+= eps
    else
        d = 1.0
    end
    f_fave = rf./r
    return (rffff .- 4.0.*rfff.*f_fave .+ 6.0.*rff.*f_fave.*f_fave .- 3.0.*r.*f_fave.*f_fave.*f_fave.*f_fave)./d

end

kurtosis=flatness