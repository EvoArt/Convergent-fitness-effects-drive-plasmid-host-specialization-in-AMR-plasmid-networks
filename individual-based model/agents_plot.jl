using GLMakie, Distributions,CSV,DataFrames
using Chain,DataFrameMacros, StatsBase, LinearAlgebra
function flip(x::Array,cost_c::Int64,ben_c::Int64)
    y = x[[x[:,1] .== cost_c]...,2]
    if length(y) >= 1

        n = sum(.!isnan.(y))
        return sum(y .== ben_c)//n
    else
        println([cost_c ben_c])
        return 0//1
    end
end

function flip(X,cost_cs,ben_cs)
    Y = Array{Rational}(undef,length(ben_cs),length(cost_cs))
    for (i,b) in enumerate(ben_cs)
        for (j,c) in enumerate(cost_cs)
            Y[i,j] = flip(X,c,b)
        end
    end
    Y
end

csvs = filter(contains("3_host_3_plas_mid_ab_agents_results"),readdir())
results = CSV.read(csvs,DataFrame)
n_plas = 3
n_host = 3
cost_function_order = [1,3,2,4]
n_outcome =4
fig3 = Figure()
Label(fig3[1, 5], "Single resistant", rotation = pi/2, font = :bold,tellheight = false)
Label(fig3[2, 5], "All resistant", rotation = pi/2, font = :bold,tellheight = false)
Colorbar(fig3[3,1:4], limits = (0,1), colormap = :thermal,vertical = false,
tellheight = true, tellwidth = false,label = "Frequency of change in connectance")
for j in 1:4
    Label(fig3[0, j],
    ["Single cost","Sublinear cost","Linear cost","Incompatible plasmids"][j] ,
    font =:bold, tellwidth = false)
    res =results[results.cost_func .==cost_function_order[j],:]
    for i in 1:2
        k = 4*(i-1)+j
        println([i,j])
        res_array = Array(res[res.resistance .==i,:])[:,[1,3]]
        z = rotr90(flip(res_array,0:n_host,0:n_host)[end:-1:begin,:])
        ax = Axis(fig3[i,j],leftspinevisible=false,rightspinevisible=false, 
        xticks = (1:n_outcome,["0","⅓","⅔","1"]), 
        yticks = (1:n_outcome,["0","⅓","⅔","1"]))
        hm = heatmap!(ax,float.(z), colormap=:thermal,colorrange = (0,1))  
        ax.title = "abcdefgh"[k] * ")"
        ax.titlefont = :regular
        ax.titlealign = :left
        if i == 2
            ax.xlabel = "Costly connectance"
        end
        if j ==1     
        ax.ylabel = "Beneficial connectance"
        end
        for k in 1:n_outcome
            for l in 1:n_outcome
                col = z[k,l] > 0.5 ? :black : :white
                text!(ax,Point2f(k,l),text =string(round(z[k,l],digits=3)), color = col, align = (:center, :center))
            end
        end
    end
end
