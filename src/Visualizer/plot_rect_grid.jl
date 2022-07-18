"""
##plot_grid_rectangular(xyz, plane::AbstractString)

"""
function plot_grid_rectangular(xyz::Array{T},plane::AbstractString) where T<:Real
    xx = 1; yy = 2; zz= 3;

    if plane == "xy"
        v1 = xyz[:,2,1,xx]
        v2 = xyz[1,:,1,yy]
    elseif plane == "xz"
        v1 = xyz[:,2,1,xx]
        v2 = xyz[1,2,:,zz]
    elseif plane == "yz"
        v1 = xyz[1,:,1,yy]
        v2 = xyz[1,2,:,zz]
    end

    xv = repeat(v1,inner=(length(v2)))
    yv = repeat(v2,outer=(length(v1)))

    p=scatter(xv,yv,label=false)
end

function plot_grid_rectangular(xyz::Array{T},plane::S,outfile::S) where {T<:Real,S<:AbstractString}
    xx = 1; yy = 2; zz= 3;

    if plane == "xy"
        v1 = xyz[:,2,1,xx]
        v2 = xyz[1,:,1,yy]
    elseif plane == "xz"
        v1 = xyz[:,2,1,xx]
        v2 = xyz[1,2,:,zz]
    elseif plane == "yz"
        v1 = xyz[1,:,1,yy]
        v2 = xyz[1,2,:,zz]
    end

    xv = repeat(v1,inner=(length(v2)))
    yv = repeat(v2,outer=(length(v1)))

    p=scatter(xv,yv,label=false)
    savefig(p,outfile)
end
