"""
    recovery_temperature(Mach, r = sqrt(0.72), gamma=1.4)

computes the recovery temperature, computed as 

 T_r/T_inf = 1 + 0.5* r * (γ-1) * M^2

"""
function recovery_temperature(Mach; r=sqrt(0.72), gamma=1.4)
    return 1 + 0.5 * r * (gamma - 1) * Mach * Mach
end