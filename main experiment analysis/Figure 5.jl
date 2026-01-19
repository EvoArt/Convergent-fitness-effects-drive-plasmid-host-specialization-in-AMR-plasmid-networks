using DataFrames, CSV,Turing, CategoricalArrays, Random, GLMakie
Random.seed!(123)
box_clrs = [colorant"#ABD9E9",colorant"#FDAE61",colorant"#D73027"]
bg = :white
fig_font = "Arial"
include("plotting.jl")
set_theme!(gg_theme)
df = CSV.read("velvet.csv",DataFrame, missingstring ="NA")
df = df[df.gloop .==0,:]
df.trt_week .= df.treatment .* string.(df.week)
df.trt_rep .= df.treatment .* string.(df.replicate)
n_reps = length(unique(df.trt_rep))
reps = levelcode.(categorical(df.trt_rep))
trts = levelcode.(categorical(df.treatment))
levs = levelcode.(categorical(df.trt_week))
trt_weeks = ["n1","w1","s1","n6","w6","s6"]
n,m = size(df)

hosts = split("bop","")
plass = split("0EGP","")
p =[]
o = []
b = []
ps = [Float64[],Float64[],Float64[]]
grps= [Int64[],Int64[],Int64[]]
xs= [Int64[],Int64[],Int64[]]

for i in 1:n
    for (h,host) in enumerate(hosts)
        N = df[i,host]
        for (p,plas) in enumerate(plass)
            comb = host * plas
            data = df[i,comb]
            if !ismissing(data)
                push!(xs[h],i)
                push!(ps[h],data/N)
                push!(grps[h],p)
            end
        end
    end
end

hs = zeros(Int,n,3)
Ps = zeros(Int,n,3)
hosts = split("bop","")
plass = split("EGP","")
for i in 1:n
    for (h,host) in enumerate(hosts)
        for (p,plas) in enumerate(plass)
            comb = host * plas
            data = df[i,comb]
            if !ismissing(data)
                hs[i,h] += data
                Ps[i,p] += data
            end
        end
    end
end

@model function multi(y)
    k = sum(y,dims = 2)
    (n,m) = size(y)
    α ~ Exponential()
    p ~ Dirichlet(m,α)
    conc ~ Exponential()

    for i in 1:n
        y[i,:] ~ DirichletMultinomial(k[i],p .*conc)
    end
end

# no plasmid chains
o_df = CSV.read("plating.csv",DataFrame)
o_df = o_df[(o_df.gloop .==0) .& (o_df.plasmids .==0),:]
tw = o_df.treatment .* string.(o_df.week) 
tws = ["n01","w01","s01","n06","w06","s06"]
x_chns = [sample(multi(Array(o_df[tw .== w,7:9]) .÷ 10),NUTS(),MCMCThreads(),1000,8) for w in tws]
# with plasmid chains
o_df = CSV.read("plating.csv",DataFrame)
o_df = o_df[(o_df.gloop .==0) .& (o_df.plasmids .==1),:]
tw = o_df.treatment .* string.(o_df.week) 
tws = ["n11","w11","s11","n16","w16","s16"]
xp_chns = [sample(multi(Array(o_df[tw .== w,7:9]) .÷ 10),NUTS(),MCMCThreads(),1000,8) for w in tws]
# plasmid proportion chains
p_chns = [sample(multi(Ps[df.trt_week .== "$tw",:]),NUTS(),MCMCThreads(),1000,8) for tw in trt_weeks]


markers = [:rect,:circle,:diamond]
markers3 = [:rect,:circle,:diamond,:rect,:circle,:diamond,:rect,:circle,:diamond]
markers2 = [:utriangle,:utriangle,:utriangle]
cmap = ["#002f4a", "#d62a29", "#fe6625"]
ts = ["No Tet","Low Tet","High Tet"]
ws = ["1","6"]
fig = Figure()
for i in 1:3
    t = ts[i]
    for j in 1:3
        chns = [x_chns,xp_chns,p_chns][j]
        chn1 = chns[i]
        chn6 = chns[i+3]
        P = [3,3,3,9][j]
        ax = Axis(fig[i,j],title = "abcdefghi"[3*(j-1)+i] * ")")
        ax.xticks = ([1,6] .+ P/48,string.([1,6]))
        if j==2 
            hideydecorations!(ax,label = false,grid = false)
        end
        if i < 3
            hidexdecorations!(ax)
        else
            hidexdecorations!(ax, ticks = false, ticklabels = false,label = false)
        end
        if j == 3
            ax.ylabel = ts[i]
            ax.ylabelcolor = box_clrs[i]
            ax.yaxisposition = :right
        end
        if j == 2
        ax.xlabel = "Week"
        end
        if (i ==2) & (j ==1)
        ax.ylabel = "Relative abundance"
        end
        clrs = [[:black,:black,:black],[:black,:black,:black],cmap][j]
         marker= [markers,markers,markers2,markers3][j]
        for p in 1:P
            off = p/8
            clr = clrs[p]
            lines!(ax,[1,6] .+off,[median(chn1["p[$p]"]),median(chn6["p[$p]"])], linewidth = 3, color = (clr,0.5))  
            scatter!(ax,[1] .+off,[median(chn1["p[$p]"])],strokewidth = 3,color = (:white,1.0), strokecolor = clr,marker = marker[p])
            scatter!(ax,[6] .+off,[median(chn6["p[$p]"])], strokewidth = 3,color = (:white,1.0), strokecolor = clr,marker = marker[p])     
            lines!(ax,[1,1] .+off,quantile(vec(chn1["p[$p]"]),[0.025,0.975]), color = :black)
            lines!(ax,[6,6] .+off,quantile(vec(chn6["p[$p]"]),[0.025,0.975]), color = :black)       
        end
    end
end
linkyaxes!(fig.content...)
fig.content[1].title = "Host\na)"
fig.content[2].title = "Host\nd)"
fig.content[3].title = "Plasmid\ng)"

mk_els =[MarkerElement(color = :black, marker = mk, markersize = 15) for mk in markers]
clr_els = [MarkerElement(color = clr, marker = :rect, markersize = 15) for clr in cmap]

fig[1:3,4] =Legend(fig,
[mk_els, clr_els],
[["Bordetella","Ochorobactrum","Pseudomonas"], ["pB12","pB5","pKJK5"]],
["Host", "Plasmid"])
