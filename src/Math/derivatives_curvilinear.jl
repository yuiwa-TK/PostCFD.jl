"""
    f::AbstractArray{S},metrics::AbstractArray{T},func_deriv::Function,idirection::Int
"""
function derivative_curvilinear(f::AbstractArray{S},metrics::AbstractArray{T},func_deriv::Function,idirection::Int) where {S, T}

    # df/dξ, df/dη, df/dζ
    dfdξ, dfdη, dfdζ = similar(f), similar(f), similar(f)
    for l in axes(f,3), k in axes(f,2)
        dfdξ[:,k,l] = func_deriv(view(f,:,k,l))
    end
    for l in axes(f,3), j in axes(f,1)
        dfdη[j,:,l] = func_deriv(view(f,j,:,l))
    end
    for k in axes(f,2), j in axes(f,1)
        dfdζ[j,k,:] = func_deriv(view(f,j,k,:))
    end

    if idirection==1 # df/dx
        return metrics[1,1,:,:,:].*dfdξ + metrics[1,2,:,:,:].*dfdη + metrics[1,3,:,:,:].*dfdζ
    elseif idirection==2 # df/dy
        return metrics[2,1,:,:,:].*dfdξ + metrics[2,2,:,:,:].*dfdη + metrics[2,3,:,:,:].*dfdζ
    elseif idirection==3 # df/dz
        return metrics[3,1,:,:,:].*dfdξ + metrics[3,2,:,:,:].*dfdη + metrics[3,3,:,:,:].*dfdζ
    else #[df/dx; df/dy; df/dz]
        return  metrics[1,1,:,:,:].*dfdξ + metrics[1,2,:,:,:].*dfdη + metrics[1,3,:,:,:].*dfdζ, 
                metrics[2,1,:,:,:].*dfdξ + metrics[2,2,:,:,:].*dfdη + metrics[2,3,:,:,:].*dfdζ,
                metrics[3,1,:,:,:].*dfdξ + metrics[3,2,:,:,:].*dfdη + metrics[3,3,:,:,:].*dfdζ
    end
end
