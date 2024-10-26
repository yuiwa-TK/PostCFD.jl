

function flatness(f,ff,fff,ffff; is_normalize=true)
    
    if is_normalize
        d = ff .- f.*f
        d .+= eps()
    else
        d = 1.0
    end
    return (ffff .- 4.0.*fff.*f .+ 6.0.*ff.*f.*f .- 3.0.*f.*f.*f.*f)./d

end


function flatness(r,rf,rff,rfff,rffff; is_normalize=true)
    
    f_fave = rf./r
    if is_normalize
        d = rff .- r.*f_fave.*f_fave
        d .+= eps()
    else
        d = 1.0
    end
    return (rffff .- 4.0.*rfff.*f_fave .+ 6.0.*rff.*f_fave.*f_fave .- 3.0.*r.*f_fave.*f_fave.*f_fave.*f_fave)./d

end