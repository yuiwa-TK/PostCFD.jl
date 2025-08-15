"""
    vorticity(velfield::AbstractArray{<:AbstractFloat,4},metrics::AbstractArray{<:AbstractFloat,5},func_deriv::Function) 
    vorticity(velfield::AbstractArray{<:AbstractFloat,4},metrics::AbstractArray{<:AbstractFloat,5},func_deriv::Function,idir) 
    ---
    vorticity(velfield::AbstractArray{<:AbstractFloat,3},metrics::AbstractArray{<:AbstractFloat,4},func_deriv::Function) 
    vorticity(velfield::AbstractArray{<:AbstractFloat,3},metrics::AbstractArray{<:AbstractFloat,4},func_deriv::Function,idir::Int) 
    vorticity(u2Dfield::AbstractArray{<:AbstractFloat,2},v2Dfield::AbstractArray{<:AbstractFloat,2},metrics::AbstractArray{<:AbstractFloat,4},func_deriv::Function) 

Compute vorticity of 3D/ 2D velocity field.


# Details of arguments
 The metrics should be defined in the way shown below:

 - For 3D case (metrics::AbstractArray{<:AbstractFloat,5})
```text

                   --             --
                  |xi_x eta_x zeta_x|   [1,1] [1,2] [1,3]
 met[j,k,l,:,:] = |xi_y eta_y zeta_y|   [2,1]  ...    :
                  |xi_z eta_z zeta_z|   [3,1]  ...  [3,3]
                   --             --
```

- For 2D case (metrics::AbstractArray{<:AbstractFloat,4})
```text

                   --        --
                  |xi_x eta_x|   [1,1] [1,2] 
  met[j,k, :,:] = |xi_y eta_y|   [2,1] [2,2] 
                   --        --
```

# example 

- When you have 3D velocity vector and want vorticity vector.

```julia
    using PostCFD
    Q = FileReader.read_flow("flowfile")
    Qp = VarsConverter.conv2prim(Q,size(Q))
    @views Uarray = Qp[:,:,:,2:4];

    your_deriv_func = PostCFD.MathLib.derivative_2ndcentral;

    G = FileReader.read_grid("gridfile")
    J, metrics = Geometry.metrics(G, your_deriv_func);

    vor = vorticity(Uarray, metrics, your_deriv_func)
```

- When you want omega_y = dw/dx - du/dz
```julia
    @views uwmat = Qp[:,1,:, [2,4]];
    @views submet = metrics[:, 1, :, [1,3],[1,3]]

    omega_y = vorticity(uwmat, submet, your_deriv_func)
    # equivalently,
    omega_y = vorticity(Qp[:,1,:,1], Qp[:,1,:,3], submet, your_deriv_func) # 1

    omega_y = vorticity(Uarray, metrics, your_deriv_func, (1,3))

```

"""
function vorticity(velfield::AbstractArray{<:AbstractFloat,4}, metrics::AbstractArray{<:AbstractFloat,5}, func_deriv::Function)
    om = similar(velfield)
    @views begin
        ufield = velfield[:, :, :, 1]
        vfield = velfield[:, :, :, 2]
        wfield = velfield[:, :, :, 3]
    end

    ~, dudy, dudz = MathLib.derivative_curvilinear(ufield, metrics, func_deriv, -1)
    dvdx, ~, dvdz = MathLib.derivative_curvilinear(vfield, metrics, func_deriv, -1)
    dwdx, dwdy, ~ = MathLib.derivative_curvilinear(wfield, metrics, func_deriv, -1)

    @views begin
        om[:, :, :, 1] .= dwdy - dvdz
        om[:, :, :, 2] .= dwdx - dudz
        om[:, :, :, 3] .= dvdx - dudy
    end
    return om
end

function vorticity(velfield::AbstractArray{<:AbstractFloat,4}, metrics::AbstractArray{<:AbstractFloat,5}, func_deriv::Function, idir)
    jmax, kmax, lmax, ~ = size(velfield)
    om = Array{3}(undef, jmax, kmax, lmax)
    d2 = similar(om)

    @views begin
        ufield = velfield[:, :, :, 1]
        vfield = velfield[:, :, :, 2]
        wfield = velfield[:, :, :, 3]
    end

    if idir == 1 # dwdy - dvdz
        om .= MathLib.derivative_curvilinear(wfield, metrics, func_deriv, 2)
        d2 .= MathLib.derivative_curvilinear(vfield, metrics, func_deriv, 3)
        om .-= d2

    elseif idir == 2 # dwdx - dudz
        om .= MathLib.derivative_curvilinear(wfield, metrics, func_deriv, 1)
        d2 .= MathLib.derivative_curvilinear(ufield, metrics, func_deriv, 3)
        om .-= d2

    elseif idir == 3 # dvdx - dudy
        om .= MathLib.derivative_curvilinear(vfield, metrics, func_deriv, 1)
        d2 .= MathLib.derivative_curvilinear(ufield, metrics, func_deriv, 2)
        om .-= d2
    else
        @error "Choose idir from [1, 2, 3]"
    end
    return om
end


function vorticity(velfield::AbstractArray{<:AbstractFloat,3}, metrics::AbstractArray{<:AbstractFloat,4}, func_deriv::Function)
    @assert (size(metrics, 3), size(metrics, 4)) == (2, 2)
    JD, KD, nv = size(velfield)
    @assert nv == 2

    om = Array{2}(undef, JD, KD)
    @views begin
        ufield = velfield[:, :, 1]
        vfield = velfield[:, :, 2]
    end

    om .= MathLib.derivative_curvilinear(vfield, metrics, func_deriv, 1) # dvdx
    om .-= MathLib.derivative_curvilinear(ufield, metrics, func_deriv, 2) # dudy

    return om # dvdx-dudy
end


function vorticity(u2Dfield::AbstractArray{<:AbstractFloat,2}, v2Dfield::AbstractArray{<:AbstractFloat,2}, metrics::AbstractArray{<:AbstractFloat,4}, func_deriv::Function)
    @assert (size(metrics, 3), size(metrics, 4)) == (2, 2)
    @assert size(u2Dfield) == size(v2Dfield)

    om = similar(u2Dfield)
    om .= MathLib.derivative_curvilinear(v2Dfield, metrics, func_deriv, 1) # dvdx
    om .-= MathLib.derivative_curvilinear(u2Dfield, metrics, func_deriv, 2) # dudy
    return om # dvdx-dudy
end

function vorticity(velfield::AbstractArray{<:AbstractFloat,4}, metrics::AbstractArray{<:AbstractFloat,5}, func_deriv::Function, idplane::Tuple{Int,Int})
    JD, KD, LD, nv = size(velfield)
    @assert nv == 3

    if idplane ∈ [(1, 2), (2, 1)]  # x-y plane
        @views begin
            idp = div(1 + LD, 2)
            ufield = velfield[:, :, idp, 1]
            vfield = velfield[:, :, idp, 2]
            submet = metrics[:, :, idp, [1, 2], [1, 2]]
        end
        return vorticity2D_2Dfield(ufield, vfield, submet, func_deriv)

    elseif idplane ∈ [(2, 3), (3, 2)] # y-z plane
        @views begin
            idp = div(1 + JD, 2)
            vfield = velfield[idp, :, :, 2]
            wfield = velfield[idp, :, :, 3]
            submet = metrics[idp, :, :, [2, 3], [2, 3]]
        end
        return vorticity2D_2Dfield(vfield, wfield, submet, func_deriv)

    elseif idplane ∈ [(3, 1), (1, 3)] # z-x plane
        @views begin
            idp = div(1 + KD, 2)
            ufield = velfield[:, idp, :, 1]
            wfield = velfield[:, idp, :, 3]
            submet = metrics[:, idp, :, [1, 3], [1, 3]]
        end
        return vorticity2D_2Dfield(ufield, wfield, submet, func_deriv)
    else
        @error "Arguments error, idplane", idplane
    end

end


"""
    alias for vorticity(velfield::AbstractArray{<:AbstractFloat,4},metrics::AbstractArray{<:AbstractFloat,5},func_deriv::Function) 
"""
vorticity3D(velfield::AbstractArray{<:AbstractFloat,4}, metrics::AbstractArray{<:AbstractFloat,5}, func_deriv::Function) = vorticity(velfield::AbstractArray{<:AbstractFloat,4}, metrics::AbstractArray{<:AbstractFloat,5}, func_deriv::Function)



"""
    alias for vorticity(velfield::AbstractArray{<:AbstractFloat,3},metrics::AbstractArray{<:AbstractFloat,4},func_deriv::Function) 
"""
vorticity2D(velfield::AbstractArray{<:AbstractFloat,3}, metrics::AbstractArray{<:AbstractFloat,4}, func_deriv::Function) = vorticity(velfield::AbstractArray{<:AbstractFloat,3}, metrics::AbstractArray{<:AbstractFloat,4}, func_deriv::Function)


"""
    alias for vorticity(velfield::AbstractArray{<:AbstractFloat,4},metrics::AbstractArray{<:AbstractFloat,5},func_deriv::Function,idir) 
"""
vorticity3D(velfield::AbstractArray{<:AbstractFloat,4}, metrics::AbstractArray{<:AbstractFloat,5}, func_deriv::Function, idir) = vorticity(velfield::AbstractArray{<:AbstractFloat,4}, metrics::AbstractArray{<:AbstractFloat,5}, func_deriv::Function, idir)

"""
    alias for vorticity(u2Dfield::AbstractArray{<:AbstractFloat,2},v2Dfield::AbstractArray{<:AbstractFloat,2},metrics::AbstractArray{<:AbstractFloat,4},func_deriv::Function) 
"""
vorticity2D_2Dfield(u2Dfield::AbstractArray{<:AbstractFloat,2}, v2Dfield::AbstractArray{<:AbstractFloat,2}, metrics::AbstractArray{<:AbstractFloat,4}, func_deriv::Function) = vorticity(u2Dfield::AbstractArray{<:AbstractFloat,2}, v2Dfield::AbstractArray{<:AbstractFloat,2}, metrics::AbstractArray{<:AbstractFloat,4}, func_deriv::Function)



"""
    vorticity_4thcentral_linearmesh(v::Matrix{<:AbstractFloat},w::Matrix{<:AbstractFloat},y::Vector{<:AbstractFloat},z::Vector{<:AbstractFloat})
compute the vortices in xi direction, assuming in the rectangular grid.
    In this function, 4th-order central derivative is used.
    dwdy - dvdy is returned
"""
function vorticity_4thcentral_linearmesh(v::Matrix{<:AbstractFloat}, w::Matrix{<:AbstractFloat}, y::Vector{<:AbstractFloat}, z::Vector{<:AbstractFloat}) 
    kmax = length(y)
    lmax = length(z)
    @assert kmax == size(v, 1)
    @assert lmax == size(v, 2)

    dydet = similar(y, kmax)
    for k = 2:kmax-1
        dydet[k] = 0.5 * (y[k+1] - y[k-1])
    end
    dydet[1] = (-3.0*y[1] + 4.0*y[2] - y[3]) * 0.5
    dydet[end] = (3.0*y[end] - 4.0*y[end-1] + y[end-2]) * 0.5

    dzdzt = similar(z, lmax)
    for l = 2:lmax-1
        dzdzt[l] = 0.5 * (z[l+1] - z[l-1])
    end
    dzdzt[1] = (-3*z[1] + 4*z[2] - z[3]) * 0.5
    dzdzt[end] = (3*z[end] - 4*z[end-1] + z[end-2]) * 0.5

    dwdy = similar(w, kmax, lmax)
    dvdz = similar(v, kmax, lmax)

    # diff for eta dir
    for l = 1:lmax
        for k = 3:kmax-2
            dwdy[k, l] = (8.0*(w[k+1, l] - w[k-1, l]) - (w[k+2, l] - w[k-2, l])) / 12.0 / dydet[k]
        end
        dwdy[1, l] = (-3w[1, l] + 4w[2, l] - w[3, l]) / (2.0 * dydet[1])
        dwdy[2, l] = 0.5 * (w[3, l] - w[1, l]) / (dydet[2])
        dwdy[end-1, l] = 0.5 * (w[end, l] - w[end-2, l]) / (dydet[end-1])
        dwdy[end, l] = (3*w[end, l] - 4*w[end-1, l] + w[end-2, l]) / (2.0 * dydet[end])
    end

    # diff for zeta dir
    for k = 1:kmax
        for l = 3:lmax-2
            dvdz[k, l] = (8.0*(v[k, l+1] - v[k, l-1]) - (v[k, l+2] - v[k, l-2])) / 12.0 / dzdzt[l]
        end
        dvdz[k, 1] = (-3*v[k, 1] + 4*v[k, 2] - v[k, 3]) * 0.5 / (dzdzt[1])
        dvdz[k, 2] = 0.5 * (v[k, 3] - v[k, 1]) / (dzdzt[2])
        dvdz[k, end-1] = 0.5 * (v[k, end] - v[k, end-2]) / (dzdzt[end-1])
        dvdz[k, end] = (3*v[k, end] - 4*v[k, end-1] + v[k, end-2]) * 0.5 / (dzdzt[end])
    end
    # return dwdy-dvdz,dwdy,dvdz
    return dwdy - dvdz
end