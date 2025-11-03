"""
    derivative_curvilinear_inplace!(df, f::AbstractArray{<:AbstractFloat},metrics::AbstractArray{T},func_deriv::Function,idirection::Int) where {<:AbstractFloat, T}

Compute derivative in curvilinear coordinate.

# Arguments
- `f::AbstractArray{<:AbstractFloat,3}`: array of variable 
- `metrics::AbstractArray{T,5}` : array of metrics whose size is [3,3, jmax, kmax, lmax]. 
- `func_deriv::Function` : function for derivative, whose argument is 1D vector and which returns the vector's derivative.
- `idirection::Int` : choose from [1,2,3] the direction 

# Details of arguments
 The metrics should be defined in the way shown below:
```text

                   --             --
                  |xi_x eta_x zeta_x|   [1,1] [1,2] [1,3]
 met[j,k,l,:,:] = |xi_y eta_y zeta_y|   [2,1]  ...    :
                  |xi_z eta_z zeta_z|   [3,1]  ...  [3,3]
                   --             --
```
"""
function derivative_curvilinear_inplace!(df, f::AbstractArray{<:AbstractFloat,3},metrics::AbstractArray{<:AbstractFloat,5},func_deriv!::Function,idirection::Int)
    JD, KD, LD = size(f);
    dfdξ, dfdη, dfdζ = similar(f), similar(f), similar(f)
    fill!(dfdξ,300)
    fill!(dfdη,300)
    fill!(dfdζ,300)

    @assert sum(isnan.(f))==0

    if JD >3 && KD > 3 && LD >3
        # df/dξ, df/dη, df/dζ
        for l in axes(f,3), k in axes(f,2)
            func_deriv!(@view(f[:,k,l]), @view(dfdξ[:,k,l]))
        end
        for l in axes(f,3), j in axes(f,1)
            func_deriv!(@view(f[j,:,l]), @view(dfdη[j,:,l]))
        end
        for k in axes(f,2), j in axes(f,1)
            func_deriv!(@view(f[j,k,:]), @view(dfdζ[j,k,:]))
        end
    elseif JD >3 && KD == 3 && LD >3
        # df/dξ, df/dη, df/dζ
        k=2;
        for l in axes(f,3)
            func_deriv!(@view(f[:,k,l]), @view(dfdξ[:,k,l]) ) # copyではなくview(配列そのもの）を渡す
        end
        for j in axes(f,1)
            func_deriv!(@view(f[j,k,:]), @view(dfdζ[j,k,:]) )
        end

        dfdξ[:,1,:] .= dfdξ[:,k,:]
        dfdξ[:,3,:] .= dfdξ[:,k,:]
        fill!(dfdη,0.0)
        dfdζ[:,1,:] .= dfdζ[:,k,:]
        dfdζ[:,3,:] .= dfdζ[:,k,:]
    else
        @error "not supported sized array", JD,KD,LD
        return 1
    end

    if idirection==1 # df/dx
        # df = similar(f);
        df .= metrics[:,:,:,1,1].*dfdξ + metrics[:,:,:,1,2].*dfdη + metrics[:,:,:,1,3].*dfdζ
        return df
    elseif idirection==2 # df/dy
        # df = similar(f);
        df .=metrics[:,:,:,2,1].*dfdξ + metrics[:,:,:,2,2].*dfdη + metrics[:,:,:,2,3].*dfdζ
        return df

    elseif idirection==3 # df/dz
        # df = similar(f)
        df .=metrics[:,:,:,3,1].*dfdξ + metrics[:,:,:,3,2].*dfdη + metrics[:,:,:,3,3].*dfdζ
        return df

    else #[df/dx; df/dy; df/dz]
        # JD,KD,LD = size(f)
        # df = similar(f,JD,KD,LD,3)
        @views begin
            df[:,:,:,1] .= metrics[:,:,:,1,1].*dfdξ + metrics[:,:,:,1,2].*dfdη + metrics[:,:,:,1,3].*dfdζ
            df[:,:,:,2] .= metrics[:,:,:,2,1].*dfdξ + metrics[:,:,:,2,2].*dfdη + metrics[:,:,:,2,3].*dfdζ
            df[:,:,:,3] .= metrics[:,:,:,3,1].*dfdξ + metrics[:,:,:,3,2].*dfdη + metrics[:,:,:,3,3].*dfdζ
        end
        return df
    end

end
