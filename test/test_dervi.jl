using Revise
using PostCFD
using BenchmarkTools
x = 0.0:1.0:150;
y = 0.5.*x.*x

d1=similar(x);
@btime derivative_2ndcentral!(y,d1);
@btime derivative_4thcentral!(y,d1);

# k1(a,b,c,ns,ne,rsvec)=kernel_tridiagonal!(a,b,c,ns,ne,rsvec, similar(x));
@btime derivative_compact_6th!(y,d1,similar(x));
