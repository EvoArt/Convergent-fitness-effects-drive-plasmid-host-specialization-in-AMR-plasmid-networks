using Pkg
Pkg.activate(".")
using Agents, Random,  Distributions, DataFrames,CSV
using Chain,DataFrameMacros, StatsBase, LinearAlgebra

conjugation_rate_multiplier = parse(Int64,ARGS[1])*0.01

@agent struct Host(NoSpaceAgent{2})
    Species::Int64
    Plasmids::Vector{Bool}
    Plasmid_effect:: Float64
    Resistant:: Bool
end

function initialize_model(;
    cf = maximum,
    inc = false,
    specialised = true,
    AB = 0.0,
    env_size = 1000,
    n_host = 3,
    n_plas = 3,
    cost =  0.5 .* rand(n_host,n_plas),
    stability = rand(n_plas,n_host),
    resistance = falses(n_plas),
    dens = 0.2,
    conjugation_rate = 0.01 .* rand(n_host,n_host,n_plas) .+  0.1I(n_host),
    combinations = zeros(1,n_host,n_plas),
    alpha = [0.5,0.1],
    seed = 23182,
    )

    rng = MersenneTwister(seed)
    properties = (
    specialised = specialised,
    inc = cf == false ? true : false,
    cost_func = cf == false ?  x-> isempty(x) ? 0.0 : maximum(x) : x-> isempty(x) ? 0.0 : cf(x),
    AB = AB,
    env_size = env_size,
    dens = dens,
    alpha = alpha,
    n_host = n_host,
    n_plas = n_plas,
    has_plas = falses(n_plas),
    cost = cost,
    stability = stability,
    resistance = resistance,
    combinations = combinations,
    host_dens = zeros(n_host),
    host_count = zeros(n_host),
    conjugation_rate= conjugation_rate,
    interactions = zeros(n_host)
    )
    model = StandardABM(Host;
        agent_step! = Host_step!, model_step! = env_step!,
        properties, rng, scheduler = Schedulers.Randomly(), warn = false,
        agents_first = false
    )
    # Add agents
    for i in 1:n_host
        for j in 1:2dens*env_size÷(3n_host)
            add_agent!(Host, model,i,falses(n_plas),0.0, false)
            model.host_count[i]+=1
        end

        if model.specialised
            for j in 1:dens*env_size÷(9n_host)
                for k in 1:n_plas
                    plas = falses(n_plas)
                    plas[k] = true
                    add_agent!(Host, model,i,plas,model.cost[i,k],model.resistance[k])
                    model.host_count[i]+=1
                    model.combinations[1,i,k] +=1
                end
            end
        else
            plasmid_costs = model.cost[i,:]
            plasmid_cost = model.cost_func(plasmid_costs)
            for j in 1:dens*env_size÷(3n_host)
            add_agent!(Host, model,i,trues(n_plas),plasmid_cost,any(model.resistance))
            model.host_count[i]+=1
            model.combinations[1,i,:] .+=1
        end
        end

    end
    model.host_dens .= model.host_count ./ model.env_size
    model.interactions .=  model.alpha[1] .* model.host_dens .+ model.alpha[2] .* sum((1 .-I(model.n_host)) .*model.host_dens',dims=2)
    return model
end

function Host_step!(bac::Host, model)
    i = bac.Species
    base_fitness = 1 - model.AB * !bac.Resistant
    p = Distributions.logistic( base_fitness-model.interactions[i]- bac.Plasmid_effect)

    if rand(abmrng(model)) ≤ p
        model.has_plas .= bac.Plasmids
         for plas in eachindex(model.has_plas)
            if model.has_plas[plas]
                model.has_plas[plas] = rand(Bernoulli(model.stability[plas,i]))
            end
        end
        replicate!(bac, model,
        Plasmids = model.has_plas, 
        Plasmid_effect =model.cost_func(model.cost[i,model.has_plas]),
        Resistant = any(model.resistance[model.has_plas]))
        model.host_count[i]+=1
        model.combinations[1,i,:] .+= model.has_plas
    elseif rand(abmrng(model)) > p
        remove_agent!(bac, model)
        model.host_count[i]-=1
        model.combinations[1,i,:] .-= bac.Plasmids
    end
end


function env_step!(model)
    model.host_dens .= model.host_count ./ model.env_size
    model.interactions .=  model.alpha[1] .* model.host_dens .+ model.alpha[2] .* sum((1 .-I(model.n_host)) .*model.host_dens',dims=2)
    effective_rate = effective_conjugation_rate(model)
    if model.inc
        for bac in allagents(model)
            if !any(bac.Plasmids)
                for plas in shuffle(eachindex(bac.Plasmids))
                    conj = rand(Bernoulli(effective_rate[bac.Species,plas]))
                    if conj
                        bac.Plasmids[plas] = true
                        bac.Plasmid_effect = model.cost[bac.Species,plas]
                        model.combinations[1,bac.Species,plas] +=1
                        bac.Resistant = any(model.resistance[bac.Plasmids])
                        break
                    end
                end
            end
        end

    else
        for bac in allagents(model)
            for plas in eachindex(bac.Plasmids)
                if !bac.Plasmids[plas]
                    conj = rand(Bernoulli(effective_rate[bac.Species,plas]))
                    if conj
                        bac.Plasmids[plas] = true
                        bac.Plasmid_effect = model.cost_func(model.cost[bac.Species,bac.Plasmids])
                        model.combinations[1,bac.Species,plas] +=1
                        bac.Resistant = any(model.resistance[bac.Plasmids])
                    end
                end
            end
        end
    end
end


function effective_conjugation_rate(combinations,conjugation_rate,env_size)

    p_conj_event = conjugation_rate ./ env_size
    p_no_conj = (1 .- p_conj_event) .^ combinations
    p_conj= 1 .- prod(p_no_conj, dims = 1)
    p_conj[1,:,:]

end
effective_conjugation_rate(model) = effective_conjugation_rate(model.combinations,model.conjugation_rate,model.env_size)

n_rep = 200
let n_host = parse(Int64,ARGS[2])
    results = Array{Any}(undef,0,7)
n_plas = parse(Int64,ARGS[3])
n_outcome = n_host+1
resistances = [falses(n_plas),trues(n_plas)]
resistances[1][1] = true
for rep in 1:n_rep
    @label loop_failed
    sub_res = Array{Any}(undef,8,7)
    row_num = 0
    env_size = 1000
    cost = rand(n_host,n_plas)
    stability = rand(n_plas,n_host)
    dens = 2.0
    conjugation_rate = conjugation_rate_multiplier .* rand(n_host,n_host,n_plas)
    conjugation_rate .+=  conjugation_rate .* I(n_host)
    alpha = [1.0,0.1]
    seed = rand(Poisson(100))

    for (cf_type,cf) in enumerate([maximum, sum, x -> mean([maximum(x), sum(x)]),false])
        for (res_type,resistance) in enumerate(resistances)
            row_num +=1
            for (AB_pres,AB) in enumerate([0.0,1.0,2.0])
                adata = [(a ->(a.Species == i), count) for i in 1:n_host]
                mdata = [m -> m.combinations[1,:,1]]
                sidmod = initialize_model(
                    cf = cf,
                    AB = AB,
                    env_size = env_size,
                    dens = dens,
                    alpha = alpha,
                    n_host = n_host,
                    n_plas = n_plas,
                    cost = cost,
                    stability = stability,
                    resistance = resistance,
                    conjugation_rate= conjugation_rate,
                    seed = seed        
                )
                df,mdf = run!(sidmod, 5000; adata,mdata,when = 5000,showprogress=true)
                if sum(Array(df[end,2:end]) .> 0) < n_host
                    @goto loop_failed
                end
                connectance = sum((mdf[end,2] ./ Array(df[end,2:end])) .> 0.01)
                sub_res[row_num,AB_pres] = connectance
                sub_res[row_num,4] = cf_type
                sub_res[row_num,5] = res_type
                sub_res[row_num,6] = round(conjugation_rate_multiplier,sigdigits=3)
                sub_res[row_num,7] = sum(Array(df[end,2:end]) .> 0)
            end
        end
    end
    results = vcat(results,sub_res)
    res = DataFrame(results,["no_tet","weak_tet","strong_tet","cost_func","resistance","rate_multiplier","host_rich"])
    CSV.write("$(ARGS[2])_host_$(ARGS[3])_plas_mid_ab_agents_results_$(ARGS[1]).csv",res)
end
end
