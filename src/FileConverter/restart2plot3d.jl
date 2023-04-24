
include("../FileReader/Readfiles.jl"); using ..FileReader
include("../FileWriter/Writefiles.jl");using ..FileWriter

function restart2plot3d(input_filename::AbstractString,output_filename::AbstractString)
    data = FileReader.read_restart(input_filename);
    # dims,params = FileReader.read_flow_params(input_filename);
    dims,params = FileReader.read_flow_params(input_filename; precision="double")

    print("Dimension:"); println(dims)
    print("Parameter:"); println(params)

    FileWriter.write_flow(output_filename,data, params ; precision="single")
end

