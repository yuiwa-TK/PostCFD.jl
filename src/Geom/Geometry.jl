module Geometry
using DelimitedFiles
include("../Math/derivatives_inplace.jl")
include("fast_jacobian_metrics.jl")

export Jacobian,Jacobian_symmetric,metrics,metrics_symmetric

"""
    J,metrics =metrics_symmetric(Grid::AbstractArray{T,4},  Func_Deriv::Function) where T

computes the metrics of grid &
returns met::Array{S}(3,3,jmax,kmax,lmax)

## Input
- Grid : 4dimensional array (i.e., size(Grid)=jmax,kmax,lmax,3)

"""
function metrics_symmetric(Grid::AbstractArray{T,4},  Func_Deriv::Function) where T
    jmax,kmax,lmax,dum1= size(Grid)

    xc = view(Grid,:,:,:,1)
    yc = view(Grid,:,:,:,2)
    zc = view(Grid,:,:,:,3)

    @assert xc ==Grid[:,:,:,1]

    xξ = Array{T}(undef,(jmax,kmax,lmax))
    yξ = Array{T}(undef,(jmax,kmax,lmax))
    zξ = Array{T}(undef,(jmax,kmax,lmax))
    xη = Array{T}(undef,(jmax,kmax,lmax))
    yη = Array{T}(undef,(jmax,kmax,lmax))
    zη = Array{T}(undef,(jmax,kmax,lmax))
    xζ = Array{T}(undef,(jmax,kmax,lmax))
    yζ = Array{T}(undef,(jmax,kmax,lmax))
    zζ = Array{T}(undef,(jmax,kmax,lmax))

    met=Array{T}(undef,(jmax,kmax,lmax,3,3))

    # compute on k=const.,l=const. plane
    # x_xi, y_xi, z_xi
    for l=1:lmax
    for k=1:kmax
        xξ[:,k,l] = Func_Deriv( view(xc, :, k, l) )
        yξ[:,k,l] = Func_Deriv( view(yc, :, k, l) ) 
        zξ[:,k,l] = Func_Deriv( view(zc, :, k, l) )
    end
    end

    # compute on j=const.,l=const. plane
    # x_eta, y_eta, z_eta
    for l=1:lmax
    for j=1:jmax
        xη[j,:,l] = Func_Deriv( view(xc, j, :, l) )
        yη[j,:,l] = Func_Deriv( view(yc, j, :, l) )
        zη[j,:,l] = Func_Deriv( view(zc, j, :, l) )
    end
    end

    # compute on k=const. plane
    # zeta-diff on j=const
    # x_zeta, y_zeta, z_zeta 
    for k=1:kmax
    for j=1:jmax
        xζ[j,k,:] = Func_Deriv( view(xc, j, k, :) )
        yζ[j,k,:] = Func_Deriv( view(yc, j, k, :) )
        zζ[j,k,:] = Func_Deriv( view(zc, j, k, :) )
    end
    end

    # for l=1:lmax
    # for k=1:kmax
    # for j=1:jmax
    #   @show  xξ[j,k,l]
    #   @show  yη[j,k,l]
    #   @show  zζ[j,k,l]     
    # end
    # end
    # end 

    # ================================================================= 
    # calculating symmetric-conservative spatial metrics
    #                                       using explicit scheme
    # ================================================================= 
    # ------------------- ξn (met[j,k,l,1,n]) -------------------------------------
    for j=1:jmax
        ff=zeros(kmax,lmax,3); 
        gg=zeros(kmax,lmax,3); 
        df1vec=zeros(kmax,lmax,3)
        df2vec=zeros(kmax,lmax,3)
        for l=1:lmax,k=1:kmax
            ff[k,l,1] = yη[j,k,l]*zc[j,k,l] - zη[j,k,l]*yc[j,k,l] #y_eta*z - z_eta*y
            ff[k,l,2] = zη[j,k,l]*xc[j,k,l] - xη[j,k,l]*zc[j,k,l] #z_eta*x - x_eta*z
            ff[k,l,3] = xη[j,k,l]*yc[j,k,l] - yη[j,k,l]*xc[j,k,l] #x_eta*y - y_eta*x

            gg[k,l,1] = yζ[j,k,l]*zc[j,k,l] - zζ[j,k,l]*yc[j,k,l] #y_zeta*z - z_zeta*y
            gg[k,l,2] = zζ[j,k,l]*xc[j,k,l] - xζ[j,k,l]*zc[j,k,l] #z_zeta*x - x_zeta*z
            gg[k,l,3] = xζ[j,k,l]*yc[j,k,l] - yζ[j,k,l]*xc[j,k,l] #x_zeta*y - y_zeta*x
        end
        # zeta diff of ff
        for k=1:kmax
            df1vec[k,:,1] = Func_Deriv( view(ff,k,:,1) ) #(y_eta*z-z_eta*y)_zeta
            df1vec[k,:,2] = Func_Deriv( view(ff,k,:,2) ) #(z_eta*x-x_eta*z)_zeta
            df1vec[k,:,3] = Func_Deriv( view(ff,k,:,3) ) #(x_eta*y-y_eta*x)_zeta
        end

        # eta diff of gg
        for l=1:lmax 
            df2vec[:,l,1] = Func_Deriv( view(gg,:,l,1) )  
            df2vec[:,l,2] = Func_Deriv( view(gg,:,l,2) )  
            df2vec[:,l,3] = Func_Deriv( view(gg,:,l,3) )  
        end

        # eta diff of gg & symmetric summasion
        for l=1:lmax
        for k=1:kmax
            #ξx/J   = {(y_eta*z-z_eta*y)_zeta - (y_zeta*z-z_zeta*y)_eta}/2
            met[j,k,l,1,1] = 0.5*( df1vec[k,l,1] - df2vec[k,l,1] ) 
            #ξy/J   = {(z_eta*x-x_eta*z)_zeta - (z_zeta*x-x_zeta*z)_eta}/2
            met[j,k,l,2,1] = 0.5*( df1vec[k,l,2] - df2vec[k,l,2] ) 
            #ξz/J   = {(x_eta*y-y_eta*x)_zeta - (x_zeta*y-y_zeta*x)_eta}/2
            met[j,k,l,3,1] = 0.5*( df1vec[k,l,3] - df2vec[k,l,3] ) 
        end
        end
    end

    # ------------------- ηn (met[j,k,l,2,n]) -------------------------------------
    for k=1:kmax
        ff=zeros(jmax,lmax,3); 
        gg=zeros(jmax,lmax,3); 
        df1vec=zeros(jmax,lmax,3)
        df2vec=zeros(jmax,lmax,3)
        for l=1:lmax,j=1:jmax
            ff[j,l,1] = yζ[j,k,l]*zc[j,k,l] - zζ[j,k,l]*yc[j,k,l] #y_zeta*z-z_zeta*y
            ff[j,l,2] = zζ[j,k,l]*xc[j,k,l] - xζ[j,k,l]*zc[j,k,l] #z_zeta*x-x_zeta*z
            ff[j,l,3] = xζ[j,k,l]*yc[j,k,l] - yζ[j,k,l]*xc[j,k,l] #x_zeta*y-y_zeta*x

            gg[j,l,1] = yξ[j,k,l]*zc[j,k,l] - zξ[j,k,l]*yc[j,k,l] #(y_xi*z-z_xi*y)
            gg[j,l,2] = zξ[j,k,l]*xc[j,k,l] - xξ[j,k,l]*zc[j,k,l] #(z_xi*x-x_xi*z)
            gg[j,l,3] = xξ[j,k,l]*yc[j,k,l] - yξ[j,k,l]*xc[j,k,l] #(x_xi*y-y_xi*x)

        end
        # xi diff of ff
        for l=1:lmax
            df1vec[:,l,1] = Func_Deriv( view(ff,:,l,1) )
            df1vec[:,l,2] = Func_Deriv( view(ff,:,l,2) )
            df1vec[:,l,3] = Func_Deriv( view(ff,:,l,3) )
        end

        for j=1:jmax
            df2vec[j,:,1] = Func_Deriv( view(gg,j,:,1))
            df2vec[j,:,2] = Func_Deriv( view(gg,j,:,2))
            df2vec[j,:,3] = Func_Deriv( view(gg,j,:,3))
        end

        # zeta diff of ff and symmetric summasion
        for l=1:lmax,j=1:jmax
            # {(y_zeta*z-z_zeta*y)_xi-(y_xi*z-z_xi*y)_zeta}/2=eta_x/J
            met[j,k,l, 1,2] = 0.5*(df1vec[j,l,1] - df2vec[j,l,1] )
            # !{(z_zeta*x-x_zeta*z)_xi-(z_xi*x-x_xi*z)_zeta}/2=eta_y/J
            met[j,k,l, 2,2] = 0.5*(df1vec[j,l,2] - df2vec[j,l,2] )
            # !{(x_zeta*y-y_zeta*x)_xi-(x_xi*y-y_xi*x)_zeta}/2=eta_z/J
            met[j,k,l, 3,2] = 0.5*(df1vec[j,l,3] - df2vec[j,l,3] )
        end
    end

    # ------------------- ζn/J (met[j,k,l,3,n]) -------------------------------------
    for l=1:lmax
        ff=zeros(jmax,kmax,3); df1vec=zeros(jmax,kmax,3)
        gg=zeros(jmax,kmax,3); df2vec=zeros(jmax,kmax,3)
        for k=1:kmax,j=1:jmax
            ff[j,k,1] = yξ[j,k,l]*zc[j,k,l] - zξ[j,k,l]*yc[j,k,l] #y_xi*z-z_xi*y
            ff[j,k,2] = zξ[j,k,l]*xc[j,k,l] - xξ[j,k,l]*zc[j,k,l] #z_xi*x-x_xi*z
            ff[j,k,3] = xξ[j,k,l]*yc[j,k,l] - yξ[j,k,l]*xc[j,k,l] #x_xi*y-y_xi*x

            gg[j,k,1] = yη[j,k,l]*zc[j,k,l] - zη[j,k,l]*yc[j,k,l] #(y_eta*z-z_eta*y)
            gg[j,k,2] = zη[j,k,l]*xc[j,k,l] - xη[j,k,l]*zc[j,k,l] #(z_eta*x-x_eta*z)
            gg[j,k,3] = xη[j,k,l]*yc[j,k,l] - yη[j,k,l]*xc[j,k,l] #(x_eta*y-y_eta*x)
        end
        # eta diff of ff on j=const. line
        for j=1:jmax
            df1vec[j,:,1] = Func_Deriv( view(ff,j,:,1) ) # (y_xi*z-z_xi*y)_eta
            df1vec[j,:,2] = Func_Deriv( view(ff,j,:,2) ) # (z_xi*x-x_xi*z)_eta
            df1vec[j,:,3] = Func_Deriv( view(ff,j,:,3) ) # (x_xi*y-y_xi*x)_eta
        end
        
        for k=1:kmax
            df2vec[:,k,1] = Func_Deriv( view(gg,:,k,1) )
            df2vec[:,k,2] = Func_Deriv( view(gg,:,k,2) )
            df2vec[:,k,3] = Func_Deriv( view(gg,:,k,3) )
        end

        # xi diff of ff on k=const. and symmetric sum
        for k=1:kmax,j=1:jmax
            met[j,k,l,1,3] = 0.5*( df1vec[j,k,1] - df2vec[j,k,1] )
            met[j,k,l,2,3] = 0.5*( df1vec[j,k,2] - df2vec[j,k,2] )
            met[j,k,l,3,3] = 0.5*( df1vec[j,k,3] - df2vec[j,k,3] )
        end
    end

    # -------------------------------------------------------------------------------
    # Jacobian 
    # Ref. 
    # 1. Abe, 2011, 25th CFD sympo (Eq.39)
    # -------------------------------------------------------------------------------
    Jacobi_inv = zeros(jmax,kmax,lmax)
    # Jacobi_inv/3 = [ zeta_z/J * z + zeta_x/J * x + zeta_y/J * y ]_zeta (: zeta-Func_deriv part)
    #               +[ eta_z/J  * z + eta_x/J  * x + eta_y/J  * y ]_eta  (: eta-Func_deriv part)
    #               +[ xi_z/J   * z + xi_x/J   * x + xi_y/J   * y ]_xi   (: xi-Func_deriv part)

    # zeta-Func_deriv part
    for k=1:kmax
        for j=1:jmax
            ff=zeros(lmax)
            for l=1:lmax
                ff[l] = met[j,k,l,3,3]*zc[j,k,l] # zeta_z/J * z
                      + met[j,k,l,1,3]*xc[j,k,l] # zeta_x/J * x
                      + met[j,k,l,2,3]*yc[j,k,l] # zeta_y/J * y
            end
            Jacobi_inv[j,k,:]    += Func_Deriv(ff) # defiv for zeta-direction
        end
    end
    
    # eta -Func_deriv part
    for l=1:lmax
        for j=1:jmax
            ff=zeros(kmax)
            for k=1:kmax
                ff[k] = met[j,k,l,3,2] * zc[j,k,l]
                      + met[j,k,l,1,2] * xc[j,k,l]
                      + met[j,k,l,2,2] * yc[j,k,l]
            end
            Jacobi_inv[j,:,l] +=Func_Deriv(ff)
        end
    end

    # xi-Func_deriv part
    for l=1:lmax
        for k=1:kmax
            ff=zeros(jmax)
            for j=1:jmax
                ff[j] = met[j,k,l,3,1] * zc[j,k,l]
                      + met[j,k,l,2,1] * xc[j,k,l]
                      + met[j,k,l,1,1] * yc[j,k,l]
            end
            Jacobi_inv[:,k,l] +=Func_Deriv(ff)
        end
    end

    # divided by 3
    Jacobi_inv = Jacobi_inv/3

    # ---------------------------------------------------------
    # compute symmetric Jaccobian & Metrics
    # ---------------------------------------------------------
    negatives=[];
    for l=1:lmax
        for k=1:kmax
            for j=1:jmax
                Jacobi_inv[j,k,l] = 1/Jacobi_inv[j,k,l]

                met[j,k,l,1,1] = met[j,k,l,1,1] * Jacobi_inv[j,k,l]
                met[j,k,l,2,1] = met[j,k,l,2,1] * Jacobi_inv[j,k,l]
                met[j,k,l,3,1] = met[j,k,l,3,1] * Jacobi_inv[j,k,l]
                met[j,k,l,1,2] = met[j,k,l,1,2] * Jacobi_inv[j,k,l]
                met[j,k,l,2,2] = met[j,k,l,2,2] * Jacobi_inv[j,k,l]
                met[j,k,l,3,2] = met[j,k,l,3,2] * Jacobi_inv[j,k,l]
                met[j,k,l,1,3] = met[j,k,l,1,3] * Jacobi_inv[j,k,l]
                met[j,k,l,2,3] = met[j,k,l,2,3] * Jacobi_inv[j,k,l]
                met[j,k,l,3,3] = met[j,k,l,3,3] * Jacobi_inv[j,k,l]
            end
        end
    end
    # if length(negatives) !=0
    #     open("negative_Jacobian.dat","w+") do io
    #         writedlm(io,reshape(negatives,4,:)')
    #     end
    #     @error "negative_Jacobian"
    #     return 1
    # end
        
    if any(x->x<=0, J)
        @error "negativeJacobian"
        return 1
    end

    # [Note] variable "Jacobi_inv" is "not" the inverse value of Jacobian.
    return Jacobi_inv, met
end

"""
    J = Jacobian_symmetric(Grid::AbstractArray{T,4},  Func_Deriv::Function) where T

computes the Jacobian

## Input
- Grid : 4dimensional array (i.e., size(Grid)=jmax,kmax,lmax,3)
- Func_Deriv
"""
function Jacobian_symmetric(Grid::AbstractArray{T,4},  Func_Deriv::Function) where T
    jmax,kmax,lmax,dum1= size(Grid)

    xc = view(Grid,:,:,:,1)
    yc = view(Grid,:,:,:,2)
    zc = view(Grid,:,:,:,3)

    @assert xc ==Grid[:,:,:,1]

    xξ = Array{T}(undef,(jmax,kmax,lmax))
    yξ = Array{T}(undef,(jmax,kmax,lmax))
    zξ = Array{T}(undef,(jmax,kmax,lmax))
    xη = Array{T}(undef,(jmax,kmax,lmax))
    yη = Array{T}(undef,(jmax,kmax,lmax))
    zη = Array{T}(undef,(jmax,kmax,lmax))
    xζ = Array{T}(undef,(jmax,kmax,lmax))
    yζ = Array{T}(undef,(jmax,kmax,lmax))
    zζ = Array{T}(undef,(jmax,kmax,lmax))

    met=Array{T}(undef,(jmax,kmax,lmax,3,3))

    # compute on k=const.,l=const. plane
    # x_xi, y_xi, z_xi
    for l=1:lmax
    for k=1:kmax
        xξ[:,k,l] = Func_Deriv( view(xc, :, k, l) )
        yξ[:,k,l] = Func_Deriv( view(yc, :, k, l) ) 
        zξ[:,k,l] = Func_Deriv( view(zc, :, k, l) )
    end
    end

    # compute on j=const.,l=const. plane
    # x_eta, y_eta, z_eta
    for l=1:lmax
    for j=1:jmax
        xη[j,:,l] = Func_Deriv( view(xc, j, :, l) )
        yη[j,:,l] = Func_Deriv( view(yc, j, :, l) )
        zη[j,:,l] = Func_Deriv( view(zc, j, :, l) )
    end
    end

    # compute on k=const. plane
    # zeta-diff on j=const
    # x_zeta, y_zeta, z_zeta 
    for k=1:kmax
    for j=1:jmax
        xζ[j,k,:] = Func_Deriv( view(xc, j, k, :) )
        yζ[j,k,:] = Func_Deriv( view(yc, j, k, :) )
        zζ[j,k,:] = Func_Deriv( view(zc, j, k, :) )
    end
    end

    # ================================================================= 
    # calculating symmetric-conservative spatial metrics
    #                                       using explicit scheme
    # ================================================================= 
    # ------------------- ξn (met[j,k,l,1,n]) -------------------------------------
    for j=1:jmax
        ff=zeros(kmax,lmax,3); 
        gg=zeros(kmax,lmax,3); 
        df1vec=zeros(kmax,lmax,3)
        df2vec=zeros(kmax,lmax,3)
        for l=1:lmax,k=1:kmax
            ff[k,l,1] = yη[j,k,l]*zc[j,k,l] - zη[j,k,l]*yc[j,k,l] #y_eta*z - z_eta*y
            ff[k,l,2] = zη[j,k,l]*xc[j,k,l] - xη[j,k,l]*zc[j,k,l] #z_eta*x - x_eta*z
            ff[k,l,3] = xη[j,k,l]*yc[j,k,l] - yη[j,k,l]*xc[j,k,l] #x_eta*y - y_eta*x

            gg[k,l,1] = yζ[j,k,l]*zc[j,k,l] - zζ[j,k,l]*yc[j,k,l] #y_zeta*z - z_zeta*y
            gg[k,l,2] = zζ[j,k,l]*xc[j,k,l] - xζ[j,k,l]*zc[j,k,l] #z_zeta*x - x_zeta*z
            gg[k,l,3] = xζ[j,k,l]*yc[j,k,l] - yζ[j,k,l]*xc[j,k,l] #x_zeta*y - y_zeta*x
        end
        # zeta diff of ff
        for k=1:kmax
            df1vec[k,:,1] = Func_Deriv( view(ff,k,:,1) ) #(y_eta*z-z_eta*y)_zeta
            df1vec[k,:,2] = Func_Deriv( view(ff,k,:,2) ) #(z_eta*x-x_eta*z)_zeta
            df1vec[k,:,3] = Func_Deriv( view(ff,k,:,3) ) #(x_eta*y-y_eta*x)_zeta
        end

        # eta diff of gg
        for l=1:lmax 
            df2vec[:,l,1] = Func_Deriv( view(gg,:,l,1) )  
            df2vec[:,l,2] = Func_Deriv( view(gg,:,l,2) )  
            df2vec[:,l,3] = Func_Deriv( view(gg,:,l,3) )  
        end

        # eta diff of gg & symmetric summasion
        for l=1:lmax
        for k=1:kmax
            #ξx/J   = {(y_eta*z-z_eta*y)_zeta - (y_zeta*z-z_zeta*y)_eta}/2
            met[j,k,l,1,1] = 0.5*( df1vec[k,l,1] - df2vec[k,l,1] ) 
            #ξy/J   = {(z_eta*x-x_eta*z)_zeta - (z_zeta*x-x_zeta*z)_eta}/2
            met[j,k,l,2,1] = 0.5*( df1vec[k,l,2] - df2vec[k,l,2] ) 
            #ξz/J   = {(x_eta*y-y_eta*x)_zeta - (x_zeta*y-y_zeta*x)_eta}/2
            met[j,k,l,3,1] = 0.5*( df1vec[k,l,3] - df2vec[k,l,3] ) 
        end
        end
    end

    # ------------------- ηn (met[j,k,l,2,n]) -------------------------------------
    for k=1:kmax
        ff=zeros(jmax,lmax,3); 
        gg=zeros(jmax,lmax,3); 
        df1vec=zeros(jmax,lmax,3)
        df2vec=zeros(jmax,lmax,3)
        for l=1:lmax,j=1:jmax
            ff[j,l,1] = yζ[j,k,l]*zc[j,k,l] - zζ[j,k,l]*yc[j,k,l] #y_zeta*z-z_zeta*y
            ff[j,l,2] = zζ[j,k,l]*xc[j,k,l] - xζ[j,k,l]*zc[j,k,l] #z_zeta*x-x_zeta*z
            ff[j,l,3] = xζ[j,k,l]*yc[j,k,l] - yζ[j,k,l]*xc[j,k,l] #x_zeta*y-y_zeta*x

            gg[j,l,1] = yξ[j,k,l]*zc[j,k,l] - zξ[j,k,l]*yc[j,k,l] #(y_xi*z-z_xi*y)
            gg[j,l,2] = zξ[j,k,l]*xc[j,k,l] - xξ[j,k,l]*zc[j,k,l] #(z_xi*x-x_xi*z)
            gg[j,l,3] = xξ[j,k,l]*yc[j,k,l] - yξ[j,k,l]*xc[j,k,l] #(x_xi*y-y_xi*x)

        end
        # xi diff of ff
        for l=1:lmax
            df1vec[:,l,1] = Func_Deriv( view(ff,:,l,1) )
            df1vec[:,l,2] = Func_Deriv( view(ff,:,l,2) )
            df1vec[:,l,3] = Func_Deriv( view(ff,:,l,3) )
        end

        for j=1:jmax
            df2vec[j,:,1] = Func_Deriv( view(gg,j,:,1))
            df2vec[j,:,2] = Func_Deriv( view(gg,j,:,2))
            df2vec[j,:,3] = Func_Deriv( view(gg,j,:,3))
        end

        # zeta diff of ff and symmetric summasion
        for l=1:lmax,j=1:jmax
            # {(y_zeta*z-z_zeta*y)_xi-(y_xi*z-z_xi*y)_zeta}/2=eta_x/J
            met[j,k,l, 1,2] = 0.5*(df1vec[j,l,1] - df2vec[j,l,1] )
            # !{(z_zeta*x-x_zeta*z)_xi-(z_xi*x-x_xi*z)_zeta}/2=eta_y/J
            met[j,k,l, 2,2] = 0.5*(df1vec[j,l,2] - df2vec[j,l,2] )
            # !{(x_zeta*y-y_zeta*x)_xi-(x_xi*y-y_xi*x)_zeta}/2=eta_z/J
            met[j,k,l, 3,2] = 0.5*(df1vec[j,l,3] - df2vec[j,l,3] )
        end
    end

    # ------------------- ζn/J (met[j,k,l,3,n]) -------------------------------------
    for l=1:lmax
        ff=zeros(jmax,kmax,3); df1vec=zeros(jmax,kmax,3)
        gg=zeros(jmax,kmax,3); df2vec=zeros(jmax,kmax,3)
        for k=1:kmax,j=1:jmax
            ff[j,k,1] = yξ[j,k,l]*zc[j,k,l] - zξ[j,k,l]*yc[j,k,l] #y_xi*z-z_xi*y
            ff[j,k,2] = zξ[j,k,l]*xc[j,k,l] - xξ[j,k,l]*zc[j,k,l] #z_xi*x-x_xi*z
            ff[j,k,3] = xξ[j,k,l]*yc[j,k,l] - yξ[j,k,l]*xc[j,k,l] #x_xi*y-y_xi*x

            gg[j,k,1] = yη[j,k,l]*zc[j,k,l] - zη[j,k,l]*yc[j,k,l] #(y_eta*z-z_eta*y)
            gg[j,k,2] = zη[j,k,l]*xc[j,k,l] - xη[j,k,l]*zc[j,k,l] #(z_eta*x-x_eta*z)
            gg[j,k,3] = xη[j,k,l]*yc[j,k,l] - yη[j,k,l]*xc[j,k,l] #(x_eta*y-y_eta*x)
        end
        # eta diff of ff on j=const. line
        for j=1:jmax
            df1vec[j,:,1] = Func_Deriv( view(ff,j,:,1) ) # (y_xi*z-z_xi*y)_eta
            df1vec[j,:,2] = Func_Deriv( view(ff,j,:,2) ) # (z_xi*x-x_xi*z)_eta
            df1vec[j,:,3] = Func_Deriv( view(ff,j,:,3) ) # (x_xi*y-y_xi*x)_eta
        end
        
        for k=1:kmax
            df2vec[:,k,1] = Func_Deriv( view(gg,:,k,1) )
            df2vec[:,k,2] = Func_Deriv( view(gg,:,k,2) )
            df2vec[:,k,3] = Func_Deriv( view(gg,:,k,3) )
        end

        # xi diff of ff on k=const. and symmetric sum
        for k=1:kmax,j=1:jmax
            met[j,k,l,1,3] = 0.5*( df1vec[j,k,1] - df2vec[j,k,1] )
            met[j,k,l,2,3] = 0.5*( df1vec[j,k,2] - df2vec[j,k,2] )
            met[j,k,l,3,3] = 0.5*( df1vec[j,k,3] - df2vec[j,k,3] )
        end
    end

    # -------------------------------------------------------------------------------
    # Jacobian 
    # Ref. 
    # 1. Abe, 2011, 25th CFD sympo (Eq.39)
    # -------------------------------------------------------------------------------
    Jacobi_inv = zeros(jmax,kmax,lmax)
    # Jacobi_inv/3 = [ zeta_z/J * z + zeta_x/J * x + zeta_y/J * y ]_zeta (: zeta-Func_deriv part)
    #               +[ eta_z/J  * z + eta_x/J  * x + eta_y/J  * y ]_eta  (: eta-Func_deriv part)
    #               +[ xi_z/J   * z + xi_x/J   * x + xi_y/J   * y ]_xi   (: xi-Func_deriv part)

    # zeta-Func_deriv part
    for k=1:kmax
        for j=1:jmax
            ff=zeros(lmax)
            for l=1:lmax
                ff[l] = met[j,k,l,3,3]*zc[j,k,l] # zeta_z/J * z
                      + met[j,k,l,1,3]*xc[j,k,l] # zeta_x/J * x
                      + met[j,k,l,2,3]*yc[j,k,l] # zeta_y/J * y
            end
            Jacobi_inv[j,k,:]    += Func_Deriv(ff) # defiv for zeta-direction
        end
    end
    
    # eta -Func_deriv part
    for l=1:lmax
        for j=1:jmax
            ff=zeros(kmax)
            for k=1:kmax
                ff[k] = met[j,k,l,3,2] * zc[j,k,l]
                      + met[j,k,l,1,2] * xc[j,k,l]
                      + met[j,k,l,2,2] * yc[j,k,l]
            end
            Jacobi_inv[j,:,l] +=Func_Deriv(ff)
        end
    end

    # xi-Func_deriv part
    for l=1:lmax
        for k=1:kmax
            ff=zeros(jmax)
            for j=1:jmax
                ff[j] = met[j,k,l,3,1] * zc[j,k,l]
                      + met[j,k,l,2,1] * xc[j,k,l]
                      + met[j,k,l,1,1] * yc[j,k,l]
            end
            Jacobi_inv[:,k,l] +=Func_Deriv(ff)
        end
    end

    # divided by 3
    Jacobi_inv = Jacobi_inv/3

    # ---------------------------------------------------------
    # compute symmetric Jaccobian & Metrics
    # ---------------------------------------------------------
    negatives=[];
    for l=1:lmax
        for k=1:kmax
            for j=1:jmax
                if Jacobi_inv[j,k,l]<=0
                    # @error "Negative Jacobian", Jacobi_inv[j,k,l],"@(j,k,l)=",j,k,l
                    append!(negatives,j,k,l,Jacobi_inv[j,k,l])
                end
                Jacobi_inv[j,k,l] = 1/Jacobi_inv[j,k,l]
            end
        end
    end
    if length(negatives) !=0
        open("negative_Jacobian.dat","w+") do io
            writedlm(io,reshape(negatives,4,:)')
        end
        @error "negative_Jacobian"
        return 1
    end

    # [Note] variable "Jacobi_inv" is "not" the inverse value of Jacobian.
    return Jacobi_inv
end

"""
    metrics(Grid::AbstractArray{T,4}, Func_Deriv1dvec::Function) where T
returns Jaccobian and metrics
"""
function metrics(Grid::AbstractArray{T,4}, Func_Deriv1dvec::Function) where T
    jmax,kmax,lmax,dum1= size(Grid)

    xξ = Array{T}(undef,(jmax,kmax,lmax))
    yξ = Array{T}(undef,(jmax,kmax,lmax))
    zξ = Array{T}(undef,(jmax,kmax,lmax))
    xη = Array{T}(undef,(jmax,kmax,lmax))
    yη = Array{T}(undef,(jmax,kmax,lmax))
    zη = Array{T}(undef,(jmax,kmax,lmax))
    xζ = Array{T}(undef,(jmax,kmax,lmax))
    yζ = Array{T}(undef,(jmax,kmax,lmax))
    zζ = Array{T}(undef,(jmax,kmax,lmax))
    met= Array{T}(undef,(jmax,kmax,lmax,3,3))
    
    # compute on k=const.,l=const. plane
    # x_xi, y_xi, z_xi
    if jmax > 1
        for l in axes(Grid,3)
            for k in axes(Grid,2)
                xξ[:,k,l] = Func_Deriv1dvec( view(Grid, :, k, l,1) )
                yξ[:,k,l] = Func_Deriv1dvec( view(Grid, :, k, l,2) ) 
                zξ[:,k,l] = Func_Deriv1dvec( view(Grid, :, k, l,3) )
            end
        end
    else
        fill!(xξ, 1.0)
        fill!(yξ, 0.0)
        fill!(zξ, 0.0)
    end

    # compute on j=const.,l=const. plane
    # x_eta, y_eta, z_eta
    if kmax >1
        for l in axes(Grid,3)
            for j in axes(Grid,1)
                xη[j,:,l] = Func_Deriv1dvec( view(Grid, j, :, l,1) )
                yη[j,:,l] = Func_Deriv1dvec( view(Grid, j, :, l,2) )
                zη[j,:,l] = Func_Deriv1dvec( view(Grid, j, :, l,3) )
            end
        end
    else
        fill!(xη, 0.0)
        fill!(yη, 1.0)
        fill!(zη, 0.0)
    end
    
    # compute on k=const. plane
    # zeta-diff on j=const
    # x_zeta, y_zeta, z_zeta 
    if lmax >1
        for k in axes(Grid,2)
            for j in axes(Grid,1)
                xζ[j,k,:] = Func_Deriv1dvec( view(Grid, j, k, :,1) )
                yζ[j,k,:] = Func_Deriv1dvec( view(Grid, j, k, :,2) )
                zζ[j,k,:] = Func_Deriv1dvec( view(Grid, j, k, :,3) )
            end
        end
    else
        fill!(xζ, 0.0)
        fill!(yζ, 0.0)
        fill!(zζ, 1.0)
    end

    # ================================================================
    #               Computing metrics
    # ================================================================
    xix=(yη.*zζ .- yζ.*zη)
    etx=(yζ.*zξ .- yξ.*zζ)
    ztx=(yξ.*zη .- yη.*zξ)
    J = xξ.*xix + xη.*etx + xζ.*ztx; J=1.0./J;
    
    # met (銀本の転置)
    met[:,:,:,1,1] = J.*xix                 #xi_x
    met[:,:,:,2,1] = J.*(zη.*xζ .- zζ.*xη)  #xi_y
    met[:,:,:,3,1] = J.*(xη.*yζ .- xζ.*yη)  #xi_z
    met[:,:,:,1,2] = J.*etx                 #eta_x 
    met[:,:,:,2,2] = J.*(zζ.*xξ .- zξ.*xζ)  #eta_y
    met[:,:,:,3,2] = J.*(xζ.*yξ .- xξ.*yζ)  #eta_z
    met[:,:,:,1,3] = J.*ztx                 #zeta_x
    met[:,:,:,2,3] = J.*(zξ.*xη .- zη.*xξ)  #zeta_y
    met[:,:,:,3,3] = J.*(xξ.*yη .- xη.*yξ)  #zeta_z

    if any(x->x<=0, J)
        @error "negativeJacobian"
        return 1
    end
    return J,met
end


function Jacobian(Grid::AbstractArray{T,4}, Func_Deriv1dvec::Function) where T
    jmax,kmax,lmax,dum1= size(Grid)
    xξ = Array{T}(undef,(jmax,kmax,lmax))
    yξ = Array{T}(undef,(jmax,kmax,lmax))
    zξ = Array{T}(undef,(jmax,kmax,lmax))
    xη = Array{T}(undef,(jmax,kmax,lmax))
    yη = Array{T}(undef,(jmax,kmax,lmax))
    zη = Array{T}(undef,(jmax,kmax,lmax))
    xζ = Array{T}(undef,(jmax,kmax,lmax))
    yζ = Array{T}(undef,(jmax,kmax,lmax))
    zζ = Array{T}(undef,(jmax,kmax,lmax))

    # compute on k=const.,l=const. plane
    # x_xi, y_xi, z_xi
    if jmax > 1
        @inbounds begin 
        for l in axes(Grid,3)
            for k in axes(Grid,2)
                xξ[:,k,l] .= Func_Deriv1dvec( view(Grid, :, k, l,1) )
                yξ[:,k,l] .= Func_Deriv1dvec( view(Grid, :, k, l,2) ) 
                zξ[:,k,l] .= Func_Deriv1dvec( view(Grid, :, k, l,3) )
            end
        end
        end
    else
        fill!(xξ, 1.0)
        fill!(yξ, 0.0)
        fill!(zξ, 0.0)
    end

    # compute on j=const.,l=const. plane
    # x_eta, y_eta, z_eta
    if kmax >1
        @inbounds begin 
        for l in axes(Grid,3)
            for j in axes(Grid,1)
                xη[j,:,l] .=  Func_Deriv1dvec( view(Grid, j, :, l,1) )
                yη[j,:,l] .=  Func_Deriv1dvec( view(Grid, j, :, l,2) )
                zη[j,:,l] .=  Func_Deriv1dvec( view(Grid, j, :, l,3) )
            end
        end
        end
    else
        fill!(xη, 0.0)
        fill!(yη, 1.0)
        fill!(zη, 0.0)
    end
    
    # compute on k=const. plane
    # zeta-diff on j=const
    # x_zeta, y_zeta, z_zeta 
    if lmax >1
        @inbounds begin 
        for k in axes(Grid,2)
            for j in axes(Grid,1)
                xζ[j,k,:] .=  Func_Deriv1dvec( view(Grid, j, k, :,1) )
                yζ[j,k,:] .=  Func_Deriv1dvec( view(Grid, j, k, :,2) )
                zζ[j,k,:] .=  Func_Deriv1dvec( view(Grid, j, k, :,3) )
            end
        end
        end
    else
        fill!(xζ, 0.0)
        fill!(yζ, 0.0)
        fill!(zζ, 1.0)
    end

    # ================================================================
    #               Computing metrics
    # ================================================================
    xix=(yη.*zζ .- yζ.*zη)
    etx=(yζ.*zξ .- yξ.*zζ)
    ztx=(yξ.*zη .- yη.*zξ)
    J = xξ.*xix + xη.*etx + xζ.*ztx; 
    J = 1.0./J;
    
    if any(x->x<=0, J)
        @error "negativeJacobian"
        return 1
    end
    return J
end

end # module