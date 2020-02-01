using Plots, Random, Statistics, DelimitedFiles, Arpack, Distributions
plotlyjs()
#using Plotly
Random.seed!(123)

t=Levy(0,100/4)
#x=Normal(0,sqrt(2*10))
#r=Uniform()
#tau=-log(rand(r))/0.01
list=rand(t,10000)
writedlm("l2.dat", list)

#newposition=0.0
function firstPassage(Q)
    T=0;                #Total time
    steps=0;
    newposition=0;
    limit=1000000;           #limit to compare times
    for i=1:limit
        #println("increase limit")
        #println("s "*string(i))
        if i==1
            newposition=10
        end
        #println(newposition)
        tt=newposition==0 ? 0 : rand(Levy(0,newposition^2/4));  #rand(Levy(0,1/4));      #First passage time Random Variable
        tau2=-log(rand(Uniform()))/Q            #Transition time random variable
        if tt<tau2
            T=T+tt
            steps=i-1;
            break
        else
            #println(T)
            #T=T+tau2
            newposition=rand(Normal(0,sqrt(2*T)))       # T==0 ? 0 : rand(Levy(0,sqrt(2*T)))            #rand(Normal(0,sqrt(2*T)))       #New position random variable
            T=T+tau2
        end
        if i==limit
            println("increase limit")
            break
        end
    end

    return  T, steps;
end

firstPassage(0.001)
tim=[0.0]
steps=[0];
QQ=[0.00001 0.00005 0.0001 0.0005 0.001 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 200 500]
meantime=[0.0]
meanres=[0.0]
stdres=[0.0]
for s in QQ
    global tim=[0.0]
    global steps=[0]
    lim=15000000         #sample size
    for i=1:lim
        T, stps=firstPassage(s);
        append!(tim, T)
        append!(steps, stps)
    end
    popfirst!(tim)
    popfirst!(steps)
    append!(meantime,sum(tim)/length(tim))
    append!(meanres,mean(steps))
    append!(stdres,std(steps))
end
popfirst!(meantime)
popfirst!(meanres)
popfirst!(stdres)

writedlm("FPT2Reset.dat",[QQ' meantime meanres stdres]')

plot(QQ, [i for i in meantime]', seriestype=:scatter, xaxis=(:log),yaxis=(:log))

println([QQ' meantime]')

writedlm("FPT.dat",[QQ' meantime]')
