gg_theme = Theme(
	    backgroundcolor = bg,
	    Figure = (backgroundcolor = bg,),
	    Scene = (backgroundcolor = bg,),
	    Axis = (
	        titlealign = :left,
	        xgridvisible   = true,
	        ygridvisible   = true,
	        ygridcolor   = (:grey,0.3),
	        xgridcolor   = (:grey,0.3),
	        xgridwidth   = 2,
	        ygridwidth   = 2,
	        xminorticks= IntervalsBetween(10),
	        spinewidth = 4,
	        rightspinecolor = :grey,
	        leftspinecolor = :grey,
	        topspinecolor = :grey,
	        bottomspinecolor = :grey,
	        xlabelsize = 26,
	        titlefont = fig_font,
	        titlesize = 26,
	        xticklabelfont = fig_font,
	        xticklabelsize = 20,
	        yticklabelfont = fig_font,
	        ylabelfont = fig_font,
	        xlabelfont = fig_font,
	        ylabelsize = 26,
	        yticklabelsize = 26,
	        backgroundcolor= bg
	    ),
	    Axis3 = (
	        titlealign = :left,
	        xgridvisible   = false,
	        ygridvisible   = false,
	        xlabelsize = 26,
	        titlefont = fig_font,
	        titlesize = 26,
	        xticklabelfont = fig_font,
	        xticklabelsize = 26,
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
	            labelsize = 20,
	      titlesize = 20,
	    )
	)


clrs =[colorant"#ABD9E9",colorant"#FDAE61",colorant"#D73027"]

function ci(x,y=0.0;orientation = 1,p = [0.025,0.975])
    k = kde(x)
    q = quantile(x,p)
    mask = q[1] .< k.x .< q[2]
    lower = Point2f.(k.x[mask], y)
    upper = Point2f.(k.x[mask], orientation*k.density[mask]/2maximum(k.density[mask]) .+y)
    orientation==-1 ? (lower,upper) = (upper,lower) : nothing
    lower,upper
end
function no_tet_plot!(chn,fig,loc,title)
    ax = Axis(fig[loc[1],loc[2]])
    xlims!(ax, (-0.025,0.075))
    pb5 = vec(chn["β[1]"])
    pb12 = vec(chn["β[2]"])
    pkjk5 = vec(chn["β[3]"])
    
    band!(ax,ci(pb5,1,orientation = -1,p=[0,1])...,color = (clrs[1],0.3))
    band!(ax,ci(pb5,1,orientation = -1)...,color = (clrs[1],0.6))
    band!(ax,ci(pb12,2,orientation = -1,p=[0,1])...,label = "No Tet",color = (clrs[1],0.3))
    band!(ax,ci(pb12,2,orientation = -1)...,color = (clrs[1],0.6))
    band!(ax,ci(pkjk5,3,orientation = -1,p=[0,1])...,color = (clrs[1],0.3))
    band!(ax,ci(pkjk5,3,orientation = -1)...,color = (clrs[1],0.6))
    vlines!(ax,[0], linestyle = :dash, color = :black, linewidth = 4)

    if loc[2] == 1
        ax.yticks = (1:3, ["pB5","pB12","pKJK5"])
         hidedecorations!(ax,ticks = false, ticklabels = false, grid = false)
    else
        hidexdecorations!(ax,ticks = false, ticklabels = false,grid = false,label = false)
        hideydecorations!(ax, grid = false)
    end
    if loc[2] ==2
        ax.xlabel = "Change in mean optical density"
    end
    ax.title = title
end
function tet_plot!(chns,ax)
    labs = ["Low Tet","High Tet"]
    for i in 1:2
      chn = chns[i]
      clr = clrs[2:3][i]
      pb5 = vec(chn["β[1]"])
      println(mean(pb5))
      pb12 = vec(chn["β[2]"])
      pkjk5 = vec(chn["β[3]"])
      band!(ax,ci(pb5,1,orientation = 1,p=[0,1])...,color = (clr,0.3))
      band!(ax,ci(pb5,1,orientation = 1)...,color = (clr,0.6))
      band!(ax,ci(pb12,2,orientation = 1,p=[0,1])...,color = (clr,0.3))
      band!(ax,ci(pb12,2,orientation = 1)...,color = (clr,0.6),label = labs[i])
      band!(ax,ci(pkjk5,3,orientation = 1,p=[0,1])...,color = (clr,0.3))
      band!(ax,ci(pkjk5,3,orientation = 1)...,color = (clr,0.6))
      for (j,x) in enumerate([pb5,pb12,pkjk5])
            if quantile(x,0.025) >0
              t = string(round(median(x)*100,sigdigits = 2)) *"×10⁻²"
            if i ==1
              y = j+0.5
              align = (:center, :bottom)
            else
              y = j
              align = (:center, :top)
            end
            text!(ax,(median(x),y),text = t, align = align,font = "Arial")
          end
        end
    end
end
