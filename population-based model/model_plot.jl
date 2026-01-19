using DifferentialEquations, GLMakie, DataFrames,CSV, RecursiveArrayTools, Distributions,Random,Colors, Chain, DataFrameMacros
using LinearAlgebra

bg = :white
fig_font = "Arial"
gg_theme = Theme(
    backgroundcolor = bg,
    Figure = (backgroundcolor = bg,),
    Scene = (backgroundcolor = bg,),
    Axis = (
        titlealign = :left,
        xgridvisible   = false,
        ygridvisible   = false,
        ygridcolor   = (:grey,0.3),
        xgridcolor   = (:grey,0.3),
        xgridwidth   = 2,
        ygridwidth   = 2,
        spinewidth = 4,
        rightspinecolor = :grey,
        leftspinecolor = :grey,
        topspinecolor = :grey,
        bottomspinecolor = :grey,
        xlabelsize = 20,
        titlefont = fig_font,
        titlesize = 26,
        xticklabelfont = fig_font,
        xticklabelsize = 16,
        yticklabelfont = fig_font,
        ylabelfont = fig_font,
        xlabelfont = fig_font,
        ylabelsize = 20,
        yticklabelsize = 16,
        backgroundcolor= bg,
        topspinevisible = false,
    bottomspinevisible = false
    ),
    Axis3 = (
        titlealign = :left,
        xgridvisible   = false,
        ygridvisible   = false,
        xlabelsize = 26,
        titlefont = fig_font,
        titlesize = 26,
        xticklabelfont = fig_font,
        xticklabelsize = 20,
        yticklabelfont = fig_font,
        zticklabelfont = fig_font,
        ylabelfont = fig_font,
        xlabelfont = fig_font,
        zlabelfont = fig_font,
        zlabelsize = 26,
        backgroundcolor= bg
    ),
    Legend = (
    backgroundcolor = bg,
      titlefont = fig_font,
      labelfont = fig_font,
      labelsize = 20,
      titlesize = 20,
    ),
    Colorbar = (
      titlefont = fig_font,
      labelfont = fig_font,
      labelsize = 16,
      titlesize = 16,
    )
)

set_theme!(gg_theme);


Random.seed!(2)

function affect!(integrator)
    terminate!(integrator)
end
function condition(u,t,integrator)
  t> 2500.0 && all(abs.(integrator.uprev .- u) .< 0.00001)
  end
cb =DiscreteCallback(condition,affect!)


U(R,Vₘₐₓ,Kₘ) = (Vₘₐₓ*R)/(Kₘ+R) 

function multi_levin(du,u,p,t)

    R,B₀,Bₚ = u.x
    ρ,c,d,G₀,λ,γ,Vₘₐₓ,ω,Kₘ = p
    H,M,N = length(R), length(B₀), size(Bₚ,1)
    for i in 1:M # for each strain i
        #assign plasmid free deaths 
        du.x[2][i] = -d[i]*B₀[i]*(sum(Bₚ[:,i])+B₀[i])
        for j in 1:N # for each plasmid j
            #assign plasmid carrying deaths 
            du.x[3][j,i] = -d[i]*Bₚ[j,i]*(sum(Bₚ[:,i])+B₀[i])
            for k in 1:M #for each strain k
                # subtract conjugation from plasmid-free strain
                du.x[2][i] += -γ[k,j,i]*B₀[i]*Bₚ[j,k]
                # add conjugation to plasmid-carrying strain
                du.x[3][j,i] += γ[k,j,i]*B₀[i]*Bₚ[j,k] 
            end
        end
    end

    for h in 1:H
        du.x[1][h] =  ρ[h]*R[h]*(c[h]-R[h]) 
        for i in 1:M # for each strain i
            # subtract plasmid free consumption from R
            du.x[1][h] -= U(R[h], Vₘₐₓ[i,h],Kₘ[i,h])*B₀[i]
            # add plasmid-free births 
            du.x[2][i] += G₀[i]*U(R[h], Vₘₐₓ[i,h],Kₘ[i,h])*B₀[i]
            for j in 1:N # for each plasmid j
                # subtract plasmid-carrying consumption from R
                du.x[1][h] -= U(R[h], Vₘₐₓ[i,h]*ω[j,i],Kₘ[i,h])*Bₚ[j,i] 
                # add plasmid-carrying births and losses from segregation
                du.x[3][j,i] += (1-λ[j,i])*G₀[i]*ω[j,i]*U(R[h], Vₘₐₓ[i,h]*ω[j,i],Kₘ[i,h])*Bₚ[j,i] 
                # add segregational loss to plasmid-free strain
                du.x[2][j] += λ[j,i]*G₀[i]*ω[j,i]*U(R[h], Vₘₐₓ[i,h]*ω[j,i],Kₘ[i,h])*Bₚ[j,i]
            end
        end
    end
end

H,M,N = 3,3,3
B = Array{Float64}(undef,0,8+M+2N)
kₘ = 1.0
vₘₐₓ = kₘ * 6e-10
vₘₐₓlims =  (4e-10,8e-10) .* kₘ
G₀ = 8e+8
d = 1e-4
λ = 1e-9
γ = 1e-6 
Γ =fill(γ,M,N,M)
for i in 1:M
    Γ[i,:,i] .*=1
end
w = fill(0.985,N,M)

w[1,1]+= √0.007
w[1,2:3] .-= 0.007

p = [fill(1.5,H),fill(1.0,H),rand(truncated(Normal(d,d),0.0,Inf),M),rand(truncated(Normal(G₀,100),4.8e+8,1.2e+9),M),
rand(truncated(Normal(λ,λ),0,Inf),N,M),Γ,
rand(truncated(Normal(vₘₐₓ,1e-10*kₘ),vₘₐₓlims...),M,H),w,fill(kₘ,M,H)]

prob = ODEProblem(multi_levin,ArrayPartition(fill(1.0,H),fill(500.0,M),fill(50.0,N,M)),(0,10000),p)

s = solve(prob, Tsit5(), callback = cb);
so = hcat(vcat([u.x[2]' for u in s.u]...),vcat([vec(u.x[3])' for u in s.u]...),[u.x[1][1] for u in s.u])
col = [:yellow,:red,:blue]
c1 ="#b2d5eb"
c2 ="#41b199"
c3 ="#f5c12f"
c4 ="#e13228"
mut_col = colorant"#ffcc00"
bac_col = colorant"#e3212f"
plas_cola = colorant"#0ea29f"
plas_colb = colorant"#1b708c"
col = [mut_col,plas_cola,plas_colb]
hosts = ["A","B","C"]
fig = Figure()

format(values) =  [i>1 ? "$(float(value))×10⁴" : "$(float(value))"  for (i,value) in enumerate(values)]
g3 = GridLayout(fig[1,1])

axs = [Axis(g3[i,1],yticklabelsize = 12, xlabel = "Generations",ylabelrotation = 0, ylabel = "Host " * hosts[i] * "\ndensity", ytickformat =format) for i in 1:3]
axs[1].topspinevisible = true
axs[3].bottomspinevisible = true
x = s.t[1:60]
sol = so[1:60,:] ./1000
for i in 1:M
    hlines!(axs[i],[0], color = (:grey,0.4), linestyle = :dash)
    if i <3
        hidexdecorations!(axs[i])
    end
    y = sol[:,i] 
    lower = [Point2f(x[i], 0) for i in 1:length(x)]
    upper = [Point2f(x[i], y[i]) for i in 1:length(x)]
    band!(axs[i],lower, upper; color = (bac_col, 0.99), label =("Plasmid free"))
    for j in 1:N
        k = j + M + N*(i-1)
        lower = [Point2f(x[i], y[i]) for i in 1:length(x)]
        y .+= sol[:,k]
        upper = [Point2f(x[i], y[i]) for i in 1:length(x)]
       band!(axs[i],lower, upper; color = (col[j], 0.99), label = "Plasmid $j") 
    end
    lines!(axs[i],x,y, color = :black)
end
rowgap!(g3, 0)

hideydecorations!(ax3, ticklabels = false)
hidexdecorations!(ax3,label = false, ticklabels = false, ticks = false)
ax3.yticks = (100:100:300, ["Host A","Host B","Host C"])
Legend(fig[2,1], ax3, tellheight = false, tellwidth = false,merge = true)
Label(g3[1, 1, TopLeft()], "a)",
        fontsize = 20,
        padding = (5, 5, 10, 5),
        halign = :right)



w = fill(0.985,N,M)

w .+= diagm(fill(√0.007,3))

p = [fill(1.5,H),fill(1.0,H),rand(truncated(Normal(d,d),0.0,Inf),M),rand(truncated(Normal(G₀,100),4.8e+8,1.2e+9),M),
rand(truncated(Normal(λ,λ),0,Inf),N,M),Γ,
rand(truncated(Normal(vₘₐₓ,1e-10*kₘ),vₘₐₓlims...),M,H),w,fill(kₘ,M,H)]

prob = ODEProblem(multi_levin,ArrayPartition(fill(1.0,H),fill(500.0,M),fill(50.0,N,M)),(0,10000),p)

s = solve(prob, Tsit5(), callback = cb);
so = hcat(vcat([u.x[2]' for u in s.u]...),vcat([vec(u.x[3])' for u in s.u]...),[u.x[1][1] for u in s.u])

g3 = GridLayout(fig[1,3])
axs = [Axis(g3[i,1], yticklabelsize = 12,xlabel = "Generations", ytickformat =format) for i in 1:3]
axs[1].topspinevisible = true
axs[3].bottomspinevisible = true
x = s.t[1:60]
sol = so[1:60,:] ./1000
for i in 1:M
    hlines!(axs[i],[0], color = (:grey,0.4), linestyle = :dash)
    if i <3
        hidexdecorations!(axs[i])
    end
    y = sol[:,i] 
    lower = [Point2f(x[i], 0) for i in 1:length(x)]
    upper = [Point2f(x[i], y[i]) for i in 1:length(x)]
    band!(axs[i],lower, upper; color = (bac_col, 0.99), label =("Plasmid free"))
    for j in 1:N
        k = j + M + N*(i-1)
        lower = [Point2f(x[i], y[i]) for i in 1:length(x)]
        y .+= sol[:,k]
        upper = [Point2f(x[i], y[i]) for i in 1:length(x)]
       band!(axs[i],lower, upper; color = (col[j], 0.99), label = "Plasmid $j") 
    end
    lines!(axs[i],x,y, color = :black)
end
rowgap!(g3, 0)

ax3.yticks = (100:100:300, ["Host A","Host B","Host C"])
Label(g3[1, 1, TopLeft()], "c)",
        fontsize = 20,
        padding = (5, 5, 10, 5),
        halign = :right)

w = fill(0.985,N,M)

w .+= diagm(fill(0.01,3))

p = [fill(1.5,H),fill(1.0,H),rand(truncated(Normal(d,d),0.0,Inf),M),rand(truncated(Normal(G₀,100),4.8e+8,1.2e+9),M),
rand(truncated(Normal(λ,λ),0,Inf),N,M),Γ,
rand(truncated(Normal(vₘₐₓ,1e-10*kₘ),vₘₐₓlims...),M,H),w,fill(kₘ,M,H)]

prob = ODEProblem(multi_levin,ArrayPartition(fill(1.0,H),fill(500.0,M),fill(50.0,N,M)),(0,10000),p)

s = solve(prob, Tsit5(), callback = cb);
so = hcat(vcat([u.x[2]' for u in s.u]...),vcat([vec(u.x[3])' for u in s.u]...),[u.x[1][1] for u in s.u])

g3 = GridLayout(fig[1,2])
axs = [Axis(g3[i,1],yticklabelsize = 12, xlabel = "Generations", ytickformat =format) for i in 1:3]
axs[1].topspinevisible = true
axs[3].bottomspinevisible = true
x = s.t[1:60]
sol = so[1:60,:] ./1000
for i in 1:M
    hlines!(axs[i],[0], color = (:grey,0.4), linestyle = :dash)
    if i <3
        hidexdecorations!(axs[i])
    end
    y = sol[:,i] 
    lower = [Point2f(x[i], 0) for i in 1:length(x)]
    upper = [Point2f(x[i], y[i]) for i in 1:length(x)]
    band!(axs[i],lower, upper; color = (bac_col, 0.99), label =("Plasmid free"))
    for j in 1:N
        k = j + M + N*(i-1)
        lower = [Point2f(x[i], y[i]) for i in 1:length(x)]
        y .+= sol[:,k]
        upper = [Point2f(x[i], y[i]) for i in 1:length(x)]
       band!(axs[i],lower, upper; color = (col[j], 0.99), label = "Plasmid $j") 
    end
    lines!(axs[i],x,y, color = :black)
end
rowgap!(g3, 0)
Label(g3[1, 1, TopLeft()], "b)",
        fontsize = 20,
        padding = (5, 5, 10, 5),
        halign = :right)

        all_data = CSV.read("all_data.csv",DataFrame)
        degree = all_data.x4
        all_flip = CSV.read("all_flip.csv",DataFrame)
        degree_flip = all_flip.x4
        
        z = Array{Float64}(undef,4,4)
        for i in 0:3
            for j in 0:3
                z[i+1,j+1] = sum((degree .==i) .& (degree_flip .== j))
            end
            z[i+1,:] ./= sum(z[i+1,:])
        end
        g3 = GridLayout(fig[2,2])
        
        ax = Axis(g3[2,1],leftspinevisible=false,rightspinevisible=false, xlabel = "Parasite connectance", ylabel = "Mutualist connectance", xticks = (1:4,["0","1/3","2/3","1"]), yticks = (1:4,["0","1/3","2/3","1"]))
        hm = heatmap!(ax,z, colormap=:thermal)
        Colorbar(g3[1,1], hm,vertical = false, label = "Frequency of change in connectance")
        
        for i in 1:4
            for j in 1:4
                col = z[i,j] > 0.5 ? :black : :white
                text!(ax,Point2f(i,j),text =string(round(z[i,j],digits=3)), color = col, align = (:center, :center))
            end
        end
        Label(g3[1, 1, TopLeft()], "d)",
            fontsize = 20,
                padding = (5, 5, 50, 5),
                halign = :right)
        



all_data = CSV.read("one_data.csv",DataFrame)
degree = all_data.x4[1:3:end]
all_flip = CSV.read("one_flip.csv",DataFrame)
degree_flip = all_flip.x4[1:3:end]
z = Array{Float64}(undef,4,4)
for i in 0:3
    for j in 0:3
        z[i+1,j+1] = sum((degree .==i) .& (degree_flip .== j))
    end
    z[i+1,:] ./= sum(z[i+1,:])
end

g3 = GridLayout(fig[2,3])
axs[1].topspinevisible = true
axs[3].bottomspinevisible = true
ax = Axis(g3[2,1], leftspinevisible=false,rightspinevisible=false,xlabel = "Parasite connectance", ylabel = "Mutualist connectance", xticks = (1:4,["0","1/3","2/3","1"]), yticks = (1:4,["0","1/3","2/3","1"]))
hm = heatmap!(ax,z, colorrange = (0,1),colormap=:thermal)
Colorbar(g3[1,1], hm,vertical = false, label = "Frequency of change in connectance")

for i in 1:4
    for j in 1:4
        col = z[i,j] > 0.5 ? :black : :white
        text!(ax,Point2f(i,j),text =string(round(z[i,j],digits=3)), color = col, align = (:center, :center))
    end
end
Label(g3[1, 1, TopLeft()], "e)",
        fontsize = 20,
        padding = (5, 5, 50, 5),
        halign = :right)

