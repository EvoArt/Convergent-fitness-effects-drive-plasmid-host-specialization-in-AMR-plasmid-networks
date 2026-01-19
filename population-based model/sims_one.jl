using DifferentialEquations,DataFrames,CSV, RecursiveArrayTools, Distributions, DataFrameMacros,Chain


function affect!(integrator)
    terminate!(integrator)
end
function condition(u,t,integrator)
  t> 2500.0 && all(abs.(integrator.uprev .- u) .< 0.00001)
  end
cb =DiscreteCallback(condition,affect!)


U(R,Vₘₐₓ,Kₘ) = (Vₘₐₓ*R)/(Kₘ+R) 


function simple_multi_levin(du,u,p,t)

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
            du.x[2][i] += G₀[i]*U(R[h], Vₘₐₓ[i,h],Kₘ[i,h])*B₀[i] #- d[i]*B₀[i]*B
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
B = Array{Float64}[]
Bflip = Array{Float64}[]
B2 = Array{Float64}[]

#kₘ = 1.0
#vₘₐₓ = kₘ * 6e-10
#vₘₐₓlims =  (4e-10,8e-10) .* kₘ
#G₀ = 8e+8
d = 1e-4
for j in 1:4
    for kₘ in [1.0]
        vₘₐₓ = kₘ * 6e-10
        vₘₐₓlims =  (4e-10,8e-10) .* kₘ
        for gmod in -2:2:2
            G₀ = 8e+8 *10.0^gmod
            for seg in 2:2:6 
                λ = 1e-8 * 10^seg
                for conj in 2:2:6
                    γ = 1e-10 * 10^conj
                    Γ =rand(Uniform(0.0,2γ),M,N,M)
                    for i in 1:M
                        Γ[i,:,i] .*=10
                    end
                    for w in 0.007:0.007:0.007
                        A = Array{Float64}(undef,24,8+3M)
                        A[:,end] .= w
                        A[:,end-1] .= γ
                        A[:,end-2] .= λ
                        Aflip = copy(A)
                        Threads.@threads for i in 1:8
                            W = rand(truncated(Normal(0.985,√w),0,1),N,M)
                            p = [fill(1.5,H),fill(1.0,H),rand(truncated(Normal(d,d),0.0,Inf),M),rand(truncated(Normal(G₀,100),4.8e+8,1.2e+9),M),
                            rand(truncated(Normal(λ,λ),0,Inf),N,M),Γ,
                            rand(truncated(Normal(vₘₐₓ,1e-10*kₘ),vₘₐₓlims...),M,H),W,fill(kₘ,M,H)]
                            diff = (1 - minimum(W)) + (√w)/2
                            for flip in [false, true]
                                if flip
                                    p[end-1][1,:] = W[1,:] .+ diff
                                end
                                prob = ODEProblem(simple_multi_levin,ArrayPartition(fill(1.0,H),fill(500.0,M),fill(50.0,N,M)),(0,10000),p)
                                s = solve(prob, Tsit5(), callback = cb);
                                bac = s.u[end].x[2]'
                                plas = s.u[end].x[3]
                                tot = sum(plas,dims = 1) .+ bac
                                network = plas .> tot ./ 100
                                degree = sum(network,dims = 2)
                                w =prob.p[8][:,:,1]
                                muts = sum(w .>1, dims = 2)
                                pos =  sum(w , dims = 2)
                                pro = prod(w,dims = 2)
                                if flip
                                Aflip[(i-1)*N+1:i*N,1] .= muts
                                Aflip[(i-1)*N+1:i*N,2] .= pos
                                Aflip[(i-1)*N+1:i*N,3] .= pro
                                Aflip[(i-1)*N+1:i*N,4] .= degree
                                Aflip[(i-1)*N+1:i*N,5] .= i
                                Aflip[(i-1)*N+1:i*N,6:5+M] .= tot .* ones(N)
                                Aflip[(i-1)*N+1:i*N,6+M:5+2M] .= bac .* ones(N)
                                Aflip[(i-1)*N+1:i*N,6+2M:5+3M] .= plas
                                else                  
                                A[(i-1)*N+1:i*N,1] .= muts
                                A[(i-1)*N+1:i*N,2] .= pos
                                A[(i-1)*N+1:i*N,3] .= pro
                                A[(i-1)*N+1:i*N,4] .= degree
                                A[(i-1)*N+1:i*N,5] .= i
                                A[(i-1)*N+1:i*N,6:5+M] .= tot .* ones(N)
                                A[(i-1)*N+1:i*N,6+M:5+2M] .= bac .* ones(N)
                                A[(i-1)*N+1:i*N,6+2M:5+3M] .= plas
                                end
                            end
                        end
                        push!(B,A)
                        push!(Bflip,Aflip)
                        CSV.write("one_data.csv",DataFrame(vcat(B...),:auto))
                        CSV.write("one_flip.csv",DataFrame(vcat(Bflip...),:auto))
                        println([j kₘ gmod seg conj])
                    end
                end
            end
        end
    end
end