@doc raw"""
reynolds_stress(rho, rf, rff)
computes the Reynods stress 

`` \overline{\rho} \widetilde{f''f''} ``
"""
function reynolds_stress(rho, rf, rff)
    ffave= rf./rho
    return rff.-rho.*ffave.*ffave # ρu''iu''j
end

@doc raw"""
reynolds_stress(rho, rf, rg, rff)
computes the Reynods stress 

`` \overline{\rho} \widetilde{f''f''} ``
"""
function reynolds_stress(rho, rf,rg, rfg)
    ffave= rf./rho
    gfave= rg./rho
    return rfg.-ffave.*gfave.*rg # ρu''iu''j
end

function reynolds_stress_fave(rho, ffave,gfave, rfg)
    return rfg.-rho.*ffave.*gfave # ρu''iu''j
end