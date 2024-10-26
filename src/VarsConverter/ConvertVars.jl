module VarsConverter

include("./conservative2primitive.jl")
include("./primitive2conservative.jl")
include("./velocity2vorticity.jl")
include("./reynolds_stress.jl")
include("./skewness.jl")

include("./slice_grid.jl")
include("./slice_flow.jl")

export vorticity, 
       skewness, reynolds_stress, reynolds_stress_fave, 
       prim2conv,conv2prim

"""
    convertvars(Inputs; options, mode::AbstractString)

converts certain flow variables (e.g., a set of convervative variables ) into some new variables (e.g. a set of primitive variables )
and returns "Outputs".
"""
function convertvars(Inputs; options, mode::AbstractString)

    if mode=="conv2prim"
       Outputs= conv2prim(Inputs,options)
    end
    return Outputs
end

end #module