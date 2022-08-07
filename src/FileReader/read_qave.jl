"""
 read_qave(filename::AbstractString; mode)

returns statics values Qstat.

---
>  -fave1-   -fave2-   -fave3-   -fave4-     -fave5-      -fave7-
>  1: r      30: ru    34: rr    148: ruu    158: uuu     171: uuuu
>  2: u      31: rv    35: uu    149: rvv    159: vvv     172: vvvv
>  3: v      32: rw    36: vv    150: rww    160: www     173: wwww
>  4: w      33: rT    37: ww    151: ruv
>  5: p                38: uv    152: rvw    -fave6-      -fave8-
>  6: t                39: vw    153: rwu    161: ruuu    174: ruuuu
>  7: mu               40: wu    154: rtt    162: rvvv    175: rvvvv
>  8: LAD              41: pp    155: rut    163: rwww    176: rwwww
>  9: tmusgs           42: tt    156: rvt    164: ruuv
> 10: tprsgs           43: ut    157: rwt    165: ruvv    -fave9-
> 11: turmu            44: vt                

166: rvvw    177    : Ar
> 12-17: tau(lam)      45: wt                167: rvww    178-180: u*Ar
> 18-23: tau(tub)   46-48: up                168: rwwu    181-186: uu*Ar
> 24-26: qx(lam)    49-57: uxp               169: rwuu    
> 27-29: qx(trb)    58-75: u*tau(lam)        170: ruvw    -fave10-
>                   76-93: u*tau(trb)                     186-189: Aru
>                  94-120: ux*tau(lam)                    190-198: u*Aru
>                 121-147: ux*tau(trb)                    199    : Ae
>              
---         
"""
function read_statistics(filename::AbstractString; mode::Number)
    @show filename

    dims    = Array{Int32}(undef,(3)) # jmax,kmax,lmax
    ndtotal = Array{Int32}(undef,(1))

    open(filename) do f
        read!(f,dims)
        read!(f,ndtotal)
        @show dims,ndtotal
    end
    jmax,kmax,lmax=dims
    Nstats=set_num_of_statistics_variables(mode)

    if mode ==1
        println("\t not supported yet")
    elseif mode==2
        ReAve       =  [Array{Float64}(undef,(jmax,kmax,lmax)) for _=1:Nstats[1]]
        FaAve       =  [Array{Float64}(undef,(jmax,kmax,lmax)) for _=1:Nstats[2]]
        ReFluc_2nd  =  [Array{Float64}(undef,(jmax,kmax,lmax)) for _=1:Nstats[3]]
        FaFluc_2nd  =  [Array{Float64}(undef,(jmax,kmax,lmax)) for _=1:Nstats[4]]
        ReFluc_3rd  =  [Array{Float64}(undef,(jmax,kmax,lmax)) for _=1:Nstats[5]]
        FaFluc_3rd  =  [Array{Float64}(undef,(jmax,kmax,lmax)) for _=1:Nstats[6]]
        ReFluc_4th  =  [Array{Float64}(undef,(jmax,kmax,lmax)) for _=1:Nstats[7]]
        FaFluc_4th  =  [Array{Float64}(undef,(jmax,kmax,lmax)) for _=1:Nstats[8]]
        Diffusion_Cmpct_ContinuousEq    =  [Array{Float64}(undef,(jmax,kmax,lmax)) for _=1:Nstats[9]]
        Diffusion_Cmpct_MomentEnergyEq  =  [Array{Float64}(undef,(jmax,kmax,lmax)) for _=1:Nstats[10]]

        open(filename) do f
            read!(f,dims)
            read!(f,ndtotal)
            for i=1:Nstats.ReAve; read!(f,ReAve[i]); end
            for i=1:Nstats.FaAve; read!(f,FaAve[i]); end
            for i=1:Nstats.ReFluc_2nd; read!(f,ReFluc_2nd[i]); end
            for i=1:Nstats.FaFluc_2nd; read!(f,FaFluc_2nd[i]); end
            for i=1:Nstats.ReFluc_3rd; read!(f,ReFluc_3rd[i]); end
            for i=1:Nstats.FaFluc_3rd; read!(f,FaFluc_3rd[i]); end
            for i=1:Nstats.ReFluc_4th; read!(f,ReFluc_4th[i]); end
            for i=1:Nstats.FaFluc_4th; read!(f,FaFluc_4th[i]); end
            for i=1:Nstats.Diffusion_Cmpct_ContinuousEq
                read!(f,Diffusion_Cmpct_ContinuousEq[i])
            end
            for i=1:Nstats.Diffusion_Cmpct_MomentEnergyEq
                read!(f,Diffusion_Cmpct_MomentEnergyEq[i])
            end
        end
        return (;ReAve,FaAve, ReFluc_2nd,FaFluc_2nd, ReFluc_3rd, FaFluc_3rd, ReFluc_4th, FaFluc_4th, Diffusion_Cmpct_ContinuousEq, Diffusion_Cmpct_MomentEnergyEq)
    end
end

function set_num_of_statistics_variables(mode)
    ReAve         = 29   #fave1
    FaAve         = 4    #fave2
    ReFluc_2nd    = 114  #...
    FaFluc_2nd    = 10
    ReFluc_3rd    = 3
    FaFluc_3rd    = 10
    ReFluc_4th    = 3
    FaFluc_4th    = 3
    Diffusion_Cmpct_ContinuousEq   = 10   #...
    Diffusion_Cmpct_MomentEnergyEq = 13 #fave10

    if mode==1
        return [11]
    else
        return (;ReAve, FaAve, ReFluc_2nd, FaFluc_2nd,
                ReFluc_3rd, FaFluc_3rd, ReFluc_4th, FaFluc_4th,
                Diffusion_Cmpct_ContinuousEq, Diffusion_Cmpct_MomentEnergyEq)
    end
end