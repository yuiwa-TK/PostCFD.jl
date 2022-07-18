"""
    plot 1D profile of primitive variables.
    input: grid_data, flow_data
    keywords : filename,gridrange
"""
function plot_1dprofile_prim(grid,flow; filename::String)
  @assert length(grid)==length(flow[:,1])

  if size(flow,2)==5
    fp  = conv2prim(flow,size(flow))
  elseif size(flow,2)==6
    fp  = flow
  end
  p1  = plot(fp[:,1],grid,xlabel="ρ")
  p2  = plot(fp[:,2],grid,xlabel="u")
  p3  = plot(fp[:,3],grid,xlabel="v")
  p4  = plot(fp[:,4],grid,xlabel="w")
  p5  = plot(fp[:,5],grid,xlabel="T")
  p6  = plot(fp[:,6],grid,xlabel="P")
  p   = plot(p1,p2,p3,p4,p5,p6,layout=(2,3),size=(1000,600),legend=false)
  savefig(p,filename)
end

function plot_1dprofile_prim(grid,flow; filename::String,gridrange::Tuple)
  @assert length(grid)==length(flow[:,1])

  if size(flow,2)==5
    fp  = conv2prim(flow,size(flow))
  elseif size(flow,2)==6
    fp  = flow
  end
  
  p1  = plot(fp[:,1],grid,xlabel="ρ")
  p2  = plot(fp[:,2],grid,xlabel="u")
  p3  = plot(fp[:,3],grid,xlabel="v")
  p4  = plot(fp[:,4],grid,xlabel="w")
  p5  = plot(fp[:,5],grid,xlabel="T")
  p6  = plot(fp[:,6],grid,xlabel="P",xrange=(0.71,0.72))
  p   = plot(p1,p2,p3,p4,p5,p6,layout=(2,3),yrange=gridrange, size=(1000,600),legend=false)
  savefig(p,filename)
end

"""
    plot 1D profile of conservative variables.
    input: grid_data, flow_data
    keywords : filename,gridrange
"""
function plot_1dprofile_conv(grid,flow; filename::String)
  @assert length(grid)==length(flow[:,1])

  p1  = plot(flow[:,1],grid,xlabel="ρ")
  p2  = plot(flow[:,2],grid,xlabel="ρu")
  p3  = plot(flow[:,3],grid,xlabel="ρv")
  p4  = plot(flow[:,4],grid,xlabel="ρw")
  p5  = plot(flow[:,5],grid,xlabel="ρE")
  pn  = plot()
  p   = plot(p1,p2,p3,p4,p5,pn,layout=(2,3),size=(1000,600),legend=false)
  savefig(p,filename)
end
function plot_1dprofile_conv(grid,flow; filename::String,gridrange::Tuple)
  @assert length(grid)==length(flow[:,1])

  p1  = plot(flow[:,1],grid,xlabel="ρ")
  p2  = plot(flow[:,2],grid,xlabel="ρu")
  p3  = plot(flow[:,3],grid,xlabel="ρv")
  p4  = plot(flow[:,4],grid,xlabel="ρw")
  p5  = plot(flow[:,5],grid,xlabel="ρE")
  pn  = plot()
  p   = plot(p1,p2,p3,p4,p5,pn,layout=(2,3),yrange=gridrange, size=(1000,600),legend=false)
  savefig(p,filename)
end