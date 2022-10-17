"""
    slice_flow(flow, samle1         ; mode) 
    slice_flow(flow, samle1, sample2; mode) !not implemented yet

returns sliced 2d flow field and 1d flowfield, respectively.
mode ∈ [x, y, z,(1d) xy, yz, zx(2d)]
"""
function slice_flow(flow, gsample; mode)
    idsample = 1;
    jmax,kmax,lmax,n = size(flow)
    if mode=="xy" || mode=="yx"
        idsample = argmin(abs.(flow[1,1,:,3].-gsample))
        @info "size of sliced flow::",jmax,kmax,n
        return flow[:,:,idsample,:]
    elseif mode =="yz" || mode=="zy"
        idsample = argmin(abs.(flow[:,1,1,1].-gsample))
        @info "size of sliced flow::",kmax,lmax,n
        return flow[idsample,:,:,:]
    elseif mode =="zx" || mode=="xz"
        idsample = argmin(abs.(flow[1,:,1,2].-gsample))
        @info "size of sliced flow::",jmax,lmax,n
        return flow[:,idsample,:,:]
    end
end
    