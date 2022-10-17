"""
    slice_grid(xyz, samle1         ; mode) 
    slice_grid(xyz, samle1, sample2; mode) !not implemented yet

returns sliced 2d grid and 1d grid, respectively.
mode ∈ [x, y, z,(1d) xy, yz, zx(2d)]
"""
function slice_grid(xyz, gsample; mode)
    id =[];
    idsample = 1;
    jmax,kmax,lmax,n = size(xyz)
    if mode=="xy" || mode=="yx"
        id = [1,2]
        idsample = argmin(abs.(xyz[1,1,:,3].-gsample))
        @info "size of sliced grid::",jmax,kmax,2
        return xyz[:,:,idsample,id]
    elseif mode =="yz" || mode=="zy"
        id = [2,3]
        idsample = argmin(abs.(xyz[:,1,1,1].-gsample))
        @info "size of sliced grid::",kmax,lmax,2
        return xyz[idsample,:,:,id]
    elseif mode =="zx" || mode=="xz"
        id = [1,3]
        idsample = argmin(abs.(xyz[1,:,1,2].-gsample))
        @info "size of sliced grid::",jmax,lmax,2
        return xyz[:,idsample,:,id]
    end
end
    