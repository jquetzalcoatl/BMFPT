using Roots, SpecialFunctions
using Random, Statistics, DelimitedFiles, Arpack, Distributions, Images, CSV, LinearAlgebra, SparseArrays, PlotlyJS
using ORCA

af(z,Q,dd,n,D)=z^2 * Q/(4*dd*n)-log(erf((sqrt(Q/D)*z)/2))
f(z,Q,dd,n,D)= (Q*z)/(2*dd*n)-exp(-(Q*z^2)/(4*D))*sqrt(Q/D)/(sqrt(pi)*erf((sqrt(Q/D)*z)/2))
f2(z,Q,dd,n,D)= (Q)/(2*dd*n)+exp(-(Q*z^2)/(4*D))*sqrt(Q/D)*Q*z/(2*D*sqrt(pi)*erf((sqrt(Q/D)*z)/2))+exp(-(Q*z^2)/(2*D))*Q/(pi*D*erf((sqrt(Q/D)*z)/2)^2)

ff(z)=f(z,0.5,1,1,2)

find_zero(ff,1)

function taylor(y,Q,dd,n,D)
    ff1(z)=f(z,Q,dd,n,D)
    #ff2(z)=f2(z,Q,dd,n,D)
    roo=find_zero(ff1,1.0)
    er=af(roo,Q,dd,n,D)+0.5*(y-roo)^2 * f2(roo,Q,dd,n,D)
    return er
end

function meanErfc(Q,dd,n,D)
    ff1(z)=f(z,Q,dd,n,D)
    #ff2(z)=f2(z,Q,dd,n,D)
    roo=find_zero(ff1,1.0)
    er=1-sqrt(Q/2*dd)*exp(-af(roo,Q,dd,n,D))/sqrt(f2(roo,Q,dd,n,D))
    return er
end

function prob(Q,dd,n,D)
    a=0
    if n==0
        a=erfc(10*sqrt(Q)/(2*sqrt(2)))
    elseif n==1
        a=meanErfc(Q,dd,n,D)*erf(10*sqrt(Q)/(2*sqrt(2)))
    else
        a=meanErfc(Q,dd,n,D)*erf(10*sqrt(Q)/(2*sqrt(2)))*prod(1-meanErfc(Q,dd,i,D) for i=1:n-1)
    end
end

plot(scatter(x=collect(1:10),y=[meanErfc(0.5,1,n,2) for n=1:10]))
plot([scatter(x=collect(1:10),y=[prob(0.01,1,n,2) for n=1:8]), scatter(x=collect(1:10),y=[prob(0.1,1,n,2) for n=1:8])])

plot([scatter(x=collect(0:15),y=[sum(prob(0.01,1,m,2) for m=0:n) for n=0:8]), scatter(x=collect(0:15),y=[sum(prob(0.1,1,m,2) for m=0:n) for n=0:8]), scatter(x=collect(0:15),y=[sum(prob(1,1,m,2) for m=0:n) for n=0:8])])

plot(scatter(x=collect(1:15),y=[sum(prob(0.01,1,m,2)*m for m=0:n) for n=1:15]))

plot(scatter(x=[Q for Q in Qlist],y=[sum(prob(Q,1,m,2)*m for m=0:6) for Q in Qlist]))

Qlist=[0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 10, 100, 500]

yy=[sum(prob(Q,1,m,2) for m=0:6) for Q in Qlist]
xx=[Q for Q in Qlist]

writedlm("saddlepoint.dat", [xx yy])

yy1=[sum(prob(0.01,1,m,2) for m=0:n) for n=0:8]
yy2=[sum(prob(0.1,1,m,2) for m=0:n) for n=0:8]
yy3=[sum(prob(1,1,m,2) for m=0:n) for n=0:8]
yy4=[sum(prob(10,1,m,2) for m=0:n) for n=0:8]
yy5=[sum(prob(100,1,m,2) for m=0:n) for n=0:8]
xx=[Q for Q=0:8]
writedlm("saddlepoint.dat", [xx yy1 yy2 yy3 yy4 yy5])

##################
plot([scatter(x=collect(0:0.1:10),y=[af(z,0.05,1,1,2) for z=0:0.1:10],mode="lines"), scatter(x=collect(0:0.1:10),y=[taylor(z,0.05,1,1,2) for z=0:0.1:10],mode="lines")])
