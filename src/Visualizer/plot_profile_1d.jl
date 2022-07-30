"""
    plot_1dprofile_prim(xyz,flow,[gridrange]; filename)

    plot 1D profile of primitive(conservative) variables.
    input: xyz_data, flow_data
    keywords : filename,
"""
function plot_1dprofile_prim(xyz,flow; filename::String)
  @assert length(xyz)==length(flow[:,1])

  if size(flow,2)==5
    fp  = conv2prim(flow,size(flow))
  elseif size(flow,2)==6
    fp  = flow
  end
  p1  = plot(fp[:,1],xyz,xlabel="ρ")
  p2  = plot(fp[:,2],xyz,xlabel="u")
  p3  = plot(fp[:,3],xyz,xlabel="v")
  p4  = plot(fp[:,4],xyz,xlabel="w")
  p5  = plot(fp[:,5],xyz,xlabel="T")
  p6  = plot(fp[:,6],xyz,xlabel="P")
  p   = plot(p1,p2,p3,p4,p5,p6,layout=(2,3),size=(1000,600),legend=false)
  savefig(p,filename)
end

function plot_1dprofile_prim(xyz,flow,gridrange::Tuple; filename::String)
  @assert length(xyz)==length(flow[:,1])

  if size(flow,2)==5
    fp  = conv2prim(flow,size(flow))
  elseif size(flow,2)==6
    fp  = flow
  end
  
  p1  = plot(fp[:,1],xyz,xlabel="ρ")
  p2  = plot(fp[:,2],xyz,xlabel="u")
  p3  = plot(fp[:,3],xyz,xlabel="v")
  p4  = plot(fp[:,4],xyz,xlabel="w")
  p5  = plot(fp[:,5],xyz,xlabel="T")
  p6  = plot(fp[:,6],xyz,xlabel="P")
  p   = plot(p1,p2,p3,p4,p5,p6,layout=(2,3),yrange=gridrange, size=(1000,600),legend=false)
  savefig(p,filename)
end

"""
  plot_1dprofile_conv(xyz,flow,[gridrange]; filename)

  plot 1D profile of conservative variables.
  input: xyz_data, flow_data
    keywords : filename,gridrange
"""
function plot_1dprofile_conv(xyz,flow; filename::String)
  @assert length(xyz)==length(flow[:,1])

  p1  = plot(flow[:,1],xyz,xlabel="ρ")
  p2  = plot(flow[:,2],xyz,xlabel="ρu")
  p3  = plot(flow[:,3],xyz,xlabel="ρv")
  p4  = plot(flow[:,4],xyz,xlabel="ρw")
  p5  = plot(flow[:,5],xyz,xlabel="ρE")
  pn  = plot()
  p   = plot(p1,p2,p3,p4,p5,pn,layout=(2,3),size=(1000,600),legend=false)
  savefig(p,filename)
end
function plot_1dprofile_conv(xyz,flow,gridrange::Tuple; filename::String)
  @assert length(xyz)==length(flow[:,1])

  p1  = plot(flow[:,1],xyz,xlabel="ρ")
  p2  = plot(flow[:,2],xyz,xlabel="ρu")
  p3  = plot(flow[:,3],xyz,xlabel="ρv")
  p4  = plot(flow[:,4],xyz,xlabel="ρw")
  p5  = plot(flow[:,5],xyz,xlabel="ρE")
  pn  = plot()
  p   = plot(p1,p2,p3,p4,p5,pn,layout=(2,3),yrange=gridrange, size=(1000,600),legend=false)
  savefig(p,filename)
end