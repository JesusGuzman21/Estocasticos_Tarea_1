
#Flat top laser beam (Ejercicio 4)

using StatsPlots

n=1000; 
l=200;

function f(x) #Values between a and b uniformily divided in n parts
    y=sqrt.(1 .- exp.(-x.^2)); #Define function 
    #y=10 .+ 5*sin.(5*x);
    return y
end

function IntegrationMC(n,x) #MonteCarlo 

    A=minimum(x)
    B=maximum(x)
    C=minimum(f(x))
    D=maximum(f(x))

    rand_x= rand(Uniform(A,B),n,1)
    rand_y= rand(Uniform(C,D),n,1) #Random numbers uniformily distributed between A and B
    #hits_outside= rand_y .> f(x)
    #hits_inside_x= rand_x[rand_y .<= f(rand_x)]
    hits_inside_y= rand_y[rand_y .<= f(rand_x)]
    total_hits=size(hits_inside_y,1)
    MC=(B-A)*maximum(f(x))*total_hits/n
    #MCC=round(MC,digits=4)

    #plot(x, f(x), label="ϕ(r)=(1-exp(-r²))",linewith=3,linecolor=:black,title="Monte Carlo [N=$n,I=$MCC]")
    #scatter!(rand_x,rand_y,markersize=4,color=:red,label="Hits outside")
    #scatter!(hits_inside_x,hits_inside_y,markersize=4,color=:green,label="Hits inside")
    return MC
end

function Φ(n,l) #Flat top laser beam

 a=0

 ϕ=zeros(1,l);

 r=collect(range(.1,2,l));
 r=r';

 size(r,2)

 for i=1:size(r,2)
    b=sqrt(2)*r[1,i]
    x= LinRange(a,b,n)
    ϕ[1,i]=IntegrationMC(n,x)
 end

 I=round(sum(ϕ),digits=3)

 scatter(r,ϕ,label=false,xlabel="r",ylabel="ϕ(r)",title="Área total=$I",titlefont=font(12,"Computer Modern"))

end

Φ(n,l)