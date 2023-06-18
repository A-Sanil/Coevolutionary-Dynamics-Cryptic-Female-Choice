#Code to run/solve adaptive dynamics model
#for Kustra and Alonzo "The coevolutionary dynamics of cryptic female choice"
#please send any questions to mkustra@ucsc.edu
using Plots, SymPy, Latexify, LambertW

#declaring varibles/paramters of model
#a = shape parameter
#b= center of sigmoid
#q= probability female mates twice
#x_e = equilibrium sperm
#x_m = sperm mutant
@syms a, b, q, x_e, x_m,x

#expected number of mates for mutant relative to equilibrium male
n_xm_xe=(1-(1/(1+exp(-a*(x_m-b)))))/(1-(1/(1+exp(-a*(x_e-b)))))
#println(latexify(n_xm_xe))

#expected fertilization sucess
v_xm_xe=(1-q)+2*q*(x_m/(x_e + x_m))
#latexify(v_xm_xe))

#overall fitness is product of expected matings and expected success of each mating
w_xm_xe=n_xm_xe*v_xm_xe

#calculate derivative
D_w_xm_xe=diff(w_xm_xe,x_m)

#subsittue
D_w_x=D_w_xm_xe.subs(x_m,x).subs(x_e,x)

#solve
solution=solve(D_w_x,x)
print(solution)

#println(latexify(D_w_xm_xe))
#println(latexify(simplify(D_w_x)))
#println(latexify(solution))

#code up the function of solution for plotting
function sol1(q,a,b)
    (q/2 + lambertw(q*exp(a*b - q/2)/2))/a
end

#Get specific numbers to compare to simulation
risks = [0.25,0.5,0.75,1]
#solution for no tradeoff
sol1.(risks,1/20,50)

#solution for tradeoff total investment
sol1.(risks,1/1000,2500)


#solution for tradeoff sperm number assuming equal investment
sol1.(risks,1/1000,2500).^0.5
#solution for tradeoff sperm number assuming 50 for ejaculate trait
sol1.(risks,1/1000,2500)./50


###Figures
using ColorSchemes
#function plotting

function mate(x,a,b)
    1-(1/(1+exp(-a*(x-b))))
  end
#Make the plots for a
pa1=plot(mate.(collect(1:150),1/5,50),label="α=1/5",linewidth=2.5,legend=:topright,color=:thistle2)
plot!(pa1,mate.(collect(1:150),1/10,50),label="α=1/10",xlabel="Ejaculate allocation",ylabel="Prob. mating success (n)",linewidth=2.5,color=:violet)
plot!(pa1,mate.(collect(1:150),1/20,50),label="α=1/20",xlabel="Ejaculate allocation",ylabel="Prob. mating success (n)",linewidth=2.5,color=:mediumorchid3)
plot!(pa1,mate.(collect(1:150),1/50,50),label="α=1/50",xlabel="Ejaculate allocation",ylabel="Prob. mating success (n)",linewidth=2.5,color=:darkorchid4)
      
#make the plots for b
pb1=plot(mate.(collect(1:150),1/20,25),label="β=25",linewidth=2.5,color=:khaki)
plot!(pb1,mate.(collect(1:150),1/20,50),label="β=50",xlabel="Ejaculate allocation (x)",ylabel="Prob. mating success (n)",linewidth=2.5,color=:burlywood2)
plot!(pb1,mate.(collect(1:150),1/20,75),label="β=75",xlabel="Ejaculate allocation (x)",ylabel="Prob. mating success (n)",linewidth=2.5,color=:orange2)
plot!(pb1,mate.(collect(1:150),1/20,100),label="β=100",xlabel="Ejaculate allocation (x)",ylabel="Prob. mating success (n)",linewidth=2.5,color=:darkorange3)
  

#making plots of various things
pa2=plot(collect(0:0.001:1),sol1.(collect(0:0.001:1),1/5,50),label="α=1/5", linewidth=2.5,legend=:topleft,color=:thistle2)
plot!(pa2,collect(0:0.001:1),sol1.(collect(0:0.001:1),1/10,50),label="α=1/10",linewidth=2.5,color=:violet)
plot!(pa2,collect(0:0.001:1),sol1.(collect(0:0.001:1),1/20,50),label="α=1/20",linewidth=2.5,color=:mediumorchid3)
plot!(pa2,collect(0:0.001:1),sol1.(collect(0:0.001:1),1/50,50),label="α=1/50",linewidth=2.5,color=:darkorchid4)
plot!(pa2,xlabel="Risk of sperm competiiton(q)",ylabel="Ejaculate allocation ESS (x_m)",linewidth=2.5)


#make heatmap
x=collect(0:0.001:1)
y=collect(5:0.05:55)
z=zeros((length(y),length(x)))
for i in 1:length(y)
    for j in 1:length(x)
        z[i,j]=sol1(x[j],1/y[i],50)
    end
end

#heatmap for a
pa3=heatmap(z,xticks=([0,250, 500,750, 1000],[0,0.25,0.5,0.75,1]),yticks=([0,250,500,750,1000],[5,17.5,30,42.5,55]),xlabel="Risk of Sperm Competition",ylabel="1/a",c= cgrad(:thermal,rev=true), aspect_ratio=:equal)


#ess plot for b
pb2=plot(collect(0:0.001:1),sol1.(collect(0:0.001:1),1/20,25),label="β=25",linewidth=2.5,legend=:topleft,color=:khaki)
plot!(pb2,collect(0:0.001:1),sol1.(collect(0:0.001:1),1/20,50),label="β=50",linewidth=2.5,color=:burlywood2)
plot!(pb2,collect(0:0.001:1),sol1.(collect(0:0.001:1),1/20,75),label="β=75",linewidth=2.5,color=:orange2)
plot!(pb2,collect(0:0.001:1),sol1.(collect(0:0.001:1),1/20,100),label="β=100",linewidth=2.5,color=:darkorange3)
plot!(pb2,xlabel="Risk of Sperm competiiton (q)",ylabel="Ejaculate allocation ESS(x_m)")


#making heatmap for b
xb=collect(0:0.001:1)
yb=collect(25:0.1:125)
zb=zeros((length(yb),length(xb)))
for i in 1:length(yb)
    for j in 1:length(xb)
        zb[i,j]=sol1(xb[j],1/20,yb[i])
    end
end

pb3=heatmap(zb,xticks=([0,250, 500,750, 1000],[0,0.25,0.5,0.75,1]),yticks=([0,250,500,750,1000],[25,50,75,100,125]),xlabel="Risk of Sperm Competition",ylabel="b",c= cgrad(:thermal,rev=true), aspect_ratio=:equal)
#heatmaps not shown in paper but good to see what they look like

#pall=plot(pa1,pa2,pa3,pb1,pb2,pb3,size=(1000,600))

pall=plot(pa1,pa2,pb1,pb2,size=(1000,800),margin=5Plots.mm)

#Save figure
savefig(pall,"Analytical.pdf")
