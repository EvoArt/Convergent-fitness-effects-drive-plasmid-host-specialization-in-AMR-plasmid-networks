using DataFrames,CSV,GLMakie, Turing, Random, QuadGK, LinearAlgebra,Optim, KernelDensity
fig_font = "Arial"
bg = :white #background colour for plotting
include("plotting.jl")
Random.seed!(1234)
set_theme!(gg_theme);

# Functions for inferring growth rate params
growth(r ,k ,u ,t ) = (k * u) /(u + (k -u) * exp(-r * t)) # logistic growth func
growth(r ,k ,u ,lag,t) = t>lag ? growth(r ,k ,u ,t-lag ) : u # logistic growth with lag

@model function get_growth_params(y,t, ::Type{T}=Float64) where {T}
    x = Vector{T}(undef,length(t))
    r ~ Exponential()
    k ~ Exponential()
    u ~ Exponential()
    σ ~ Exponential()
    α =0.04
    lag ~ Exponential()
    x .= growth.(r ,k ,u,lag,t) .+α
    y ~ MvNormal(x, σ^2 *I)
end

# function to estimate the area under the curve, minus th optical density of the media
function get_AUC(t,x,n_rep =8)
    if maximum(x) < 0.07 # if no growth, do not attempt to fit function to data
        return 0.0
    else
        res = Vector{Any}(undef,n_rep)
        Threads.@threads for i in 1:8
            res[i] =optimize(get_growth_params(x,t), MAP(), Optim.NelderMead(),Optim.Options(iterations=100_000,allow_f_increases=true))
        end
        best_res = res[argmax([r.lp] for r in res)]
        est = best_res.values
        func(x) = growth(est[:r] ,est[:k] ,est[:u],est[:lag] ,x )
        return (quadgk(func,0.0,maximum(t))[1])
    end
end

# Model for inferring impact of each plasmid
@model function plasmid_effects(y)
    σ ~ Exponential()
    α ~ Normal(0,10)
    β ~ filldist(Normal(0,10),3)

    y[1] .~ Normal( α, σ)
    for i in 1:3
        y[i+1] .~ Normal(β[i] + α, σ)
    end
end

# Estimate mean OD for each well
data = CSV.read("growth.csv", DataFrame)[1:end,:]
t_dict = Dict("Bord " => 36, "Och " => 24, "Pse " => 96)
hosts = ["Bord ","Bord ","Bord ","Bord ","Och ","Och ","Och ","Och ","Pse ","Pse ","Pse ","Pse "]
plas = ["no plas","pb5","pb12","pkjk5","no plas","pb5","pb12","pkjk5","no plas","pb5","pb12","pkjk5"]
trt = ["none","weak","strong","none","weak","strong","none","weak"]
dat = data[data.time .<= 1200,2:end]
t = data.time[data.time .<= 1200]
T = String[]
H = String[]
P = String[]
mean_OD = Float64[]
for i in 1:12
    t_mask = 8 .<= data.time .<= t_dict[hosts[i]]
    mask = [parse(Int,match(r"[0-9]+",name).match) ==i for name in names(dat)]
    est = get_AUC.([t[t_mask] .- minimum(t[t_mask])],eachcol(dat[t_mask,mask])) ./ t_dict[hosts[i]]
    push!(mean_OD, est...)
    push!(T, trt...)
    push!(H, fill(hosts[i],8)...)
    push!(P, fill(plas[i],8)...)
end
df = DataFrame(host = H, plas = P, trt = T, mean_OD = mean_OD)

# Estimate change in mean OD affected by each plasmid in each treatment
chns = Chains[]
for h in unique(H)
    for t in unique(T)
        mask = (h .== H) .& (T .==t) 
        sub_df = df[mask,:]
        y = [sub_df.mean_OD[sub_df.plas .== p] for p in unique(P)]
        md = plasmid_effects(y)
        chn = sample(md, NUTS(),MCMCThreads(),4000,8)
        push!(chns,chn)
    end
end

# Plot results
fig= Figure()
no_tet_plot!(chns[1],fig,[1,1], "a) Bordetella")
no_tet_plot!(chns[4],fig,[1,2],"b) Ochrobactrum")
no_tet_plot!(chns[7],fig,[1,3],"c) Pseudomonas")
tet_plot!(chns[[2,3]],fig.content[1])
tet_plot!(chns[[5,6]],fig.content[2])
tet_plot!(chns[[8,9]],fig.content[3])
[ax.xticks= -0.025:0.025:0.05 for ax in fig.content[1:3]]
linkaxes!(fig.content...)
axislegend(position = :rb)
