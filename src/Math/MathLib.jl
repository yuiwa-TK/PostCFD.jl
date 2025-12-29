module MathLib # MathLib

export derivative_1stsided, derivative_2ndcentral, derivative_4thcentral, 
        derivative_compact_6th,derivative_curvilinear,
        derivative_curvilinear_inplace!, derivative_curvilinear_inplace_compact!,
        derivative_2ndcentral!,derivative_4thcentral!,derivative_compact_6th!,generate_compact6th!,
        ∫fdy

include("derivatives.jl")
include("derivatives_inplace.jl")
include("derivatives_curvilinear.jl")
include("integrals.jl")
include("derivatives_curvilinear_inplace.jl")

end # module MathLib
