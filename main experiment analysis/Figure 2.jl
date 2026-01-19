using GLMakie, CSV, DataFrames,CategoricalArrays

bg = :white
fig_font = "Arial"
box_clrs = [colorant"#ABD9E9",colorant"#FDAE61",colorant"#D73027"]
include("plotting.jl")
set_theme!(gg_theme);

df = CSV.read("velvet.csv",DataFrame, missingstring ="NA")
df = df[(df.gloop .==0),:]
df.trt_week .= df.treatment .* string.(df.week)
df.trt_rep .= df.treatment .* string.(df.replicate)
n_reps = length(unique(df.trt_rep))
trts = levelcode.(categorical(df.treatment))
levs = levelcode.(categorical(df.trt_week))
trt_weeks = ["none 1","w1","s1","n6","w6","s6"]
n,m = size(df)



for host in split("bop","")
    sub_df = Array(df[!, host .* split("0EGP","")])
    sub_df[ismissing.(sub_df)] .=0
    df[!,host] .= sum(sub_df,dims = 2)
end

df.plas_prop .= 1 .-  sum(Array(df[!,split("bop","") .*"0"]),dims = 2) ./ sum(Array(df[!,split("bop","")]),dims = 2) 

hosts = split("bop","")
plass = split("0EGP","")
p =[]
o = []
b = []
ps = [Float64[],Float64[],Float64[]]
grps= [Int64[],Int64[],Int64[]]
xs= [Int64[],Int64[],Int64[]]

bar_trts = ["n1","n6","w1","w6","s1","s6"]

j = 0
for trt in bar_trts
    sub_df = df[df.trt_week .== trt,:]
    sub_n = size(sub_df,1)
    j+=1
    for i in 1:sub_n
        for (h,host) in enumerate(hosts)
            N = sub_df[i,host]
            for (p,plas) in enumerate(plass)
                comb = host * plas
                data = sub_df[i,comb]
                if !ismissing(data)
                    push!(xs[h],j)
                    push!(ps[h],data/N)
                    push!(grps[h],p)
                end
            end
        end
        j +=1
    end
end



X = Array(df[!,split("bop","")])
props = (X .+1) ./ sum(X .+1,dims = 2)


fig = Figure()
ax = Axis(fig[1:3,1], xticks = (1:6,split(repeat("16",3),"")), xlabel = "Week",
ylabel = "Total plasmid prevalence", title = "a)")

labs = ["No Tet","Low Tet","High Tet"]
for i in [1,3,2]
    if i == 2
        j = 1
    elseif i ==3 
        j = -1
    else
        j = 0
    end
    boxplot!(2j.+ levs[trts .==i],df.plas_prop[trts .==i], color = (box_clrs[i+j],0.7), strokewidth = 2,
show_outliers = false,medianlinewidth = 3,label = labs[i+j])
scatter!(ax, 2j.+ levs[trts .==i] .+ randn(sum(trts .==i)) ./10,df.plas_prop[trts .==i],color = :white,
strokewidth = 3,markersize = 10,strokecolor = box_clrs[i+j])
end
axislegend(; merge = true, position = :rb, framecolor = :black)

cmap = ["#002f4a", "#d62a29", "#fe6625", "#fdc04c", "#e9e3b7"]
titles = ["b) Bordatella","c) Ochrabactrum","d) Pseudomonas"]
    for i in 1:3
    ax = Axis(fig[i,2:4], title = titles[i],xgridvisible   = false,ygridvisible   = false)
    ax.xticklabelrotation = π/2
    cmap = ["#002f4a", "#d62a29", "#fe6625", "#fdc04c", "#e9e3b7"]
        tbl = (x = xs[i],
        height = ps[i],
        grp = grps[i],
        )
        barplot!(ax,tbl.x, tbl.height,
                stack = tbl.grp,
                color = tbl.grp,
                offset = 15,
                colormap = cmap
                )
                if i ==79
                    ax.xticks = (1:n,df.trt_week)
                else
                    hidexdecorations!(ax)
                end
end


Legend(fig[2,5],
[MarkerElement(color = cmap[c], marker = :rect, markersize = 15) for c in [1,2,4,5]],
["None","pB12","pB5","PKJK5"],"Plasmid")

