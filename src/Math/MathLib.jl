module MathLib # MathLib

export derivative_1stsided, derivative_2ndcentral,derivative_compact_6th,derivative_curvilinear, ∫fdy

include("derivatives.jl")
include("derivatives_inplace.jl")
include("derivatives_curvilinear.jl")
include("integrals.jl")

end # module MathLib
