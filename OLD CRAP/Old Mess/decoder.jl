using Printf
using PastaQ
using ITensors
using Random
import PastaQ: gate

#single qubit gate failure location
gate(::GateName"pz"; γ::Number) = 
  [1-γ  γ
   γ    1-γ]

# CNOT gate failure -- first qubit is control, second qubit is target
gate(::GateName"cnot"; pz₁::Number, pz₂::Number, pz₁₂::Number) = 
  [1-pz₁-pz₂-pz₁₂   pz₂   pz₁   pz₁₂
   pz₁₂   pz₁   pz2   1-pz₁-pz₂-pz₁₂
   pz₁   pz₁₂   1-pz₁-pz₂-pz₁₂   pz₂
   pz₂   1-pz₁-pz₂-pz₁₂   pz₁₂   pz₁]

gate(::GateName"Π"; σ::Int) =
  [ 1-σ 0
    0   σ]


function parameters(p::Number; kwargs...)
  pXm = get(kwargs, :pXm, 0.004) 
  κ₂  = get(kwargs, :κ₂, 1e7)   
  α²  = get(kwargs, :α², 8)   
 
  params = Dict()
  
  params["pXm"] = pXm
  params["κ₂"]  = κ₂
  params["α²"]  = α²
  
  params["p₁"] = 10 * p 
  params["p₂"] = 0.299 * p^2
  params["p₃"] = p * κ₂ * α² * (350*1e-9 + 10 / (κ₂*α²))
  
  params["pz₁"]  = 0.845 * √p
  params["pz₂"]  = 0.133 * √p
  params["pz₁₂"] = 0.133 * √p
  
  params["pₐᵢ"]  = 15 * p / 2
  
  return params
end


function generate_error(nqubits::Int, nreps::Int, params::Dict)
  pXm  = params["pXm"] 
  p₁   = params["p₁"] 
  p₂   = params["p₂"] 
  p₃   = params["p₃"] 
  pz₁  = params["pz₁"] 
  pz₂  = params["pz₂"]
  pz₁₂ = params["pz₁₂"]
  pₐᵢ  = params["pₐᵢ"]
  
  Qdata = zeros(Int64, nqubits)
  Adata = [zeros(Int64, nqubits-1) for _ in 1:nreps]

  #Adata = zeros(Int64, nreps, nqubits-1)
  # loop over repetitions
  for r ∈ 1:nreps
    
    # t = 0
    for q ∈ 1:nqubits-1
      r == 1     && rand()<p₁ && (Qdata[q] ⊻= 1)
      rand()<pₐᵢ && (Adata[r][q] ⊻= 1)
    end
    r == 1 && rand()<p₁  && (Qdata[nqubits] ⊻= 1)

    # t = 1,2
    for ξ ∈ (0,1)
      for q ∈ 1:nqubits-1
        rn = rand() 
        Adata[r][q] ⊻= (rn < (pz₁ + pz₂ + pz₁₂) ? Qdata[q+ξ] ⊻ 1 : Qdata[q+ξ])
        if (rn > pz₁) && (rn < (pz₁ + pz₂))
          Qdata[q] = (Qdata[q]+1)%2
        elseif (rn > (pz₁ + pz₂)) && (rn < (pz₁ + pz₂ + pz₁₂))
          Qdata[q] ⊻= 1 
        end
      end
      rand() < p₂ && (Qdata[nqubits] ⊻= 1)
    end
    # t = 3
    for q ∈ 1:nqubits-1
      rand() < p₃  && (Qdata[q] ⊻= 1)
      rand() < pXm && (Adata[r][q] ⊻= 1)
    end
    rand() < p₃  && (Qdata[nqubits] ⊻= 1) 
  end
  return Qdata, Adata

end

generate_error(nqubits::Int, nreps::Int, p::Number,nsamples::Int) =
  [generate_error(nqubits, nreps, p) for _ in 1:nsamples]

logical_sector(qubit_data::Vector{Int}) = 
  qubit_data[(length(qubit_data)-1) ÷ 2 + 1]

logical_sector(qubit_data::Vector{Vector{Int}}) = 
  [logical_sector(data) for data in qubit_data]




function generate_circuit(nqubits::Int, nreps::Int, mbits::Vector{Vector{Int}}, params::Dict)
  pXm  = params["pXm"] 
  p₁   = params["p₁"] 
  p₂   = params["p₂"] 
  p₃   = params["p₃"] 
  pz₁  = params["pz₁"] 
  pz₂  = params["pz₂"]
  pz₁₂ = params["pz₁₂"]
  pₐᵢ  = params["pₐᵢ"]
  gates = []
  for r ∈ 1:nreps
    # first time step
    for q ∈ 1:nqubits
      r == 1 && push!(gates, ("pz",2*q-1,(p = p₁,)))
    end
    for q ∈ 1:nqubits
      if r == 1
        push!(gates, ("pz",2*q,(p = pₐᵢ,)))
      else
        @show q,r-1,mbits[r-1]
        push!(gates, ("Π",2*q,(σ = mbits[r-1][q],))) 
      end
    end

  #  # second time step
  #  for q ∈ 1:nqubits-1
  #    push!(gates,("cnot",(2*q,2*q-1),(pz₁=pz₁,pz₂=pz₂,pz₁₂=pz₁₂))) 
  #  end
  #  push!(gates,("pz",2*nqubits-1,(p=p₂,)))
  #  
  #  # third time step
  #  for q ∈ 1:nqubits-1
  #    push!(gates,("cnot",(2*(q-1),2*q-1),(pz₁=pz₁,pz₂=pz₂,pz₁₂=pz₁₂))) 
  #  end
  #  push!(gates,("pz",nqubits,(p=p₂,)))

  #  # fourth tiem step
  #  for q ∈ 1:nqubits-1
  #    push!(gates,("pz",2*q,(p=pXm,)))
  #  end
  #  for q ∈ 1:nqubits
  #    push!(gates,("pz",2*q-1,(p=p₃,)))
  #  end
  end
   
  return gates
end



Random.seed!(1234)

nqubits  = 3
nreps    = 2
maxdim   = 100
cutoff   = 1e-10
nsamples = 1000 

p       = 0.01
δ       = 0.05

params = parameters(p)

Random.seed!(1234)
for i in 1:nsamples
  Qdata, Adata = generate_error(nqubits, nreps, params)
  #gates = generate_circuit(nqubits, nreps, Adata, params)
#  #decode(Qdata, Adata; cutoff = cutoff, maxdim = maxdim)
end

