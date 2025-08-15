module VarsConverter

using ..MathLib

include("./conservative2primitive.jl")
include("./primitive2conservative.jl")
include("./vorticity.jl")
include("./reynolds_stress.jl")
include("./skewness.jl")
include("./flatness.jl")
include("./sutherland.jl")

include("./slice_grid.jl")
include("./slice_flow.jl")

export  vorticity, vorticity2D,vorticity3D,vorticity2D_2Dfield,
        reynolds_stress, reynolds_stress_fave, 
        skewness, flatness,
        prim2conv,conv2prim,
        sutherland

end #module