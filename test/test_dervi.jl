using Revise
using PostCFD
using BenchmarkTools
using Test
x = 0.0:1.0:150;
y = x.*x;

ref = collect(2x);

d1=similar(x);
# @testset "derivatives" begin
@btime derivative_2ndcentral!(y,d1);
@test ref==d1
@btime derivative_4thcentral!(y,d1);
@test ref==d1
# end

# k1(a,b,c,ns,ne,rsvec)=kernel_tridiagonal!(a,b,c,ns,ne,rsvec, similar(x));

@btime derivative_compact_6th!(y,d1,similar(x));
@test ref≈ d1



