"""
    derivative_curvilinear(f::AbstractArray{<:AbstractFloat},metrics::AbstractArray{T},func_deriv::Function,idirection::Int) where {<:AbstractFloat, T}

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
function derivative_curvilinear(f::AbstractArray{<:AbstractFloat,3},metrics::AbstractArray{<:AbstractFloat,5},func_deriv::Function,idirection::Int)

    # df/dξ, df/dη, df/dζ
    dfdξ, dfdη, dfdζ = similar(f), similar(f), similar(f)
    for l in axes(f,3), k in axes(f,2)
        @views dfdξ[:,k,l] .= func_deriv(view(f,:,k,l))
    end
    for l in axes(f,3), j in axes(f,1)
        @views dfdη[j,:,l] .= func_deriv(view(f,j,:,l))
    end
    for k in axes(f,2), j in axes(f,1)
        @views dfdζ[j,k,:] .= func_deriv(view(f,j,k,:))
    end

    if idirection==1 # df/dx
        df = similar(f);
        df .= metrics[:,:,:,1,1].*dfdξ + metrics[:,:,:,1,2].*dfdη + metrics[:,:,:,1,3].*dfdζ
        return df
    elseif idirection==2 # df/dy
        df = similar(f);
        df .=metrics[:,:,:,2,1].*dfdξ + metrics[:,:,:,2,2].*dfdη + metrics[:,:,:,2,3].*dfdζ
        return df

    elseif idirection==3 # df/dz
        df = similar(f)
        df .=metrics[:,:,:,3,1].*dfdξ + metrics[:,:,:,3,2].*dfdη + metrics[:,:,:,3,3].*dfdζ
        return df

    else #[df/dx; df/dy; df/dz]
        JD,KD,LD = size(f)
        df = similar(f,JD,KD,LD,3)
        @views begin
            df[:,:,:,1] .= metrics[:,:,:,1,1].*dfdξ + metrics[:,:,:,1,2].*dfdη + metrics[:,:,:,1,3].*dfdζ
            df[:,:,:,2] .= metrics[:,:,:,2,1].*dfdξ + metrics[:,:,:,2,2].*dfdη + metrics[:,:,:,2,3].*dfdζ
            df[:,:,:,3] .= metrics[:,:,:,3,1].*dfdξ + metrics[:,:,:,3,2].*dfdη + metrics[:,:,:,3,3].*dfdζ
        end
        return df
    end
end

function derivative_curvilinear(f::AbstractArray{<:AbstractFloat,2},metrics::AbstractArray{<:AbstractFloat,4},func_deriv::Function,idirection::Int)

    # df/dξ, df/dη
    dfdξ, dfdη = similar(f), similar(f)
    for k in axes(f,2)
        @views dfdξ[:,k] .= func_deriv(view(f,:,k))
    end
    for j in axes(f,1)
        @views dfdη[j,:] .= func_deriv(view(f,j,:))
    end
    
    if idirection==1 # df/dx
        df = similar(f);
        df .= metrics[:,:,1,1].*dfdξ + metrics[:,:,1,2].*dfdη
        return df

    elseif idirection==2 # df/dy
        df = similar(f);
        df .=metrics[:,:,2,1].*dfdξ + metrics[:,:,2,2].*dfdη
        return df

    else #[df/dx; df/dy]
        JD, KD = size(f)
        df = similar(f,JD,KD,2)
        @views df[:,:,1].= metrics[:,:,1,1].*dfdξ + metrics[:,:,1,2].*dfdη
        @views df[:,:,2].= metrics[:,:,2,1].*dfdξ + metrics[:,:,2,2].*dfdη
        return  df
    end
end

derivative_curvilinear2D(f::AbstractArray{<:AbstractFloat,2},metrics::AbstractArray{<:AbstractFloat,4},func_deriv::Function,idirection::Int) =derivative_curvilinear(f::AbstractArray{<:AbstractFloat,2},metrics::AbstractArray{<:AbstractFloat,4},func_deriv::Function,idirection::Int)