"""
plot 2D profile of primitive variables.
input: grid_data, flow_data
keywords : filename,gridrange
"""
function plot_2dprofile_prim(grid,flow; filename::String,
                             x_range::Union{Tuple,Nothing},
                             y_range::Union{Tuple,Nothing},
                             c_range::Union{Tuple,Nothing})
    @assert size(grid[:,:,1])==size(flow[:,:,1])
    @assert size(grid,3) ==2

    nv = size(flow,3)
    if     nv==5
        fp  = VarsConverter.conv2prim(flow,size(flow))
    elseif nv==6
        fp  = flow
    end

    p1 = heatmap(grid[:,1,1],grid[1,:,2],fp[:,:,1]',xlabel="ρ")
    p2 = heatmap(grid[:,1,1],grid[1,:,2],fp[:,:,2]',xlabel="u")
    p3 = heatmap(grid[:,1,1],grid[1,:,2],fp[:,:,3]',xlabel="v")
    p4 = heatmap(grid[:,1,1],grid[1,:,2],fp[:,:,4]',xlabel="w")
    p5 = heatmap(grid[:,1,1],grid[1,:,2],fp[:,:,5]',xlabel="T")
    p6 = heatmap(grid[:,1,1],grid[1,:,2],fp[:,:,6]',xlabel="p")
    p  = plot(p1,p2,p3,p4,p5,p6,layout=(6,1),size=(1200,1200),legend=false,
        xrange=x_range,yrange=y_range,
        color=:lightrainbow)
    # display(p)
    savefig(p,filename)
end
