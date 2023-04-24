using PostCFD
using Test

@testset "PostCFD.jl" begin
    # Write your tests here.

    # using Revise
    using Random
    @testset "FileConverter" begin
        jmax=2;kmax=3;lmax=4
        fsmach  = 0.2;
        aoa     = 0.0;
        totime  = 10.0;
        nc      = 100;
        q       = rand(0.0:0.1:1.0,jmax,kmax,lmax,5)
        q       = Float32.(q)

        filename = "dataset/test.restart"
        fileout  = "dataset/test.pl3d"

        println("====writing restart file")
        FileWriter.write_restart(filename,q,[fsmach,aoa,totime],nc)
        
        println("==== ConvertFiles")
        FileConverter.restart2plot3d(filename,fileout)
        println("==== reading pl3d file")
        qa       = FileReader.read_flow_single(fileout)
        
        @test prod(q.==qa)
    end

end
