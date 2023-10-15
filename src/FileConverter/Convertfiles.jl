module FileConverter

include("../FileReader/Readfiles.jl"); import .FileReader
include("../FileWriter/Writefiles.jl");import .FileWriter

function restart2plot3d(input_filename::AbstractString,output_filename::AbstractString)
    data = FileReader.read_restart(input_filename);
    # dims,params = FileReader.read_flow_params(input_filename);
    dims,params = FileReader.read_flow_params(input_filename; precision="restart")

    print("Dimension:"); println(dims)
    print("Parameter:"); println(params)

    FileWriter.write_flow(output_filename,data, params ; precision="single")
end

function data2restart(input_filename::AbstractString,output_filename::AbstractString,params, nc)
    data = FileReader.read_flow_auto(input_filename);
    FileWriter.write_restart(output_filename, data, params, nc)
end
function data2restart(input_filename::AbstractString,output_filename::AbstractString; filetype)
    dim, params= FileReader.read_flow_params(input_filename; type=filetype)
    nc = 1
    data2restart(input_filename,output_filename,params,Int32(nc))
end

end