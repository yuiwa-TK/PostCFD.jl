module PostCFD

# Write your package code here.
"""Hello"""
function helloworld()
    println("helloworld")
end


include("./sub_readfiles/Readfiles.jl")
include("./sub_writefiles/Writefiles.jl")

end # module
