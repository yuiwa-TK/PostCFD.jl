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
reynolds_stress(rho, rf, rg, rfg)
computes the Reynods stress 

`` \overline{\rho} \widetilde{f''f''} ``
"""
function reynolds_stress(rho, rf,rg, rfg)
    return rfg.-rf.*rg./rho # ρu''iu''j
end

@doc raw"""
reynolds_stress(rho, ffave, gfave, rfg)
computes the Reynods stress 

`` \overline{\rho} \widetilde{f''f''} ``
"""
function reynolds_stress_fave(rho, ffave,gfave, rfg)
    return rfg.-rho.*ffave.*gfave # ρu''iu''j
end