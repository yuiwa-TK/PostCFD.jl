"""
    Visualizer

Flow & Grid Visualization
"""
module Visualizer

import Plots

export 

# 
include("../BasicFuncsContainer/conservative2primitive.jl") #
include("./plot_profile_1d.jl")
include("./plot_profile_2d.jl")

#
include("./plot_rect_grid.jl")

end