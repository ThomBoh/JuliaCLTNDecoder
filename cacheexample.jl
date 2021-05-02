using Random
using Printf
using ITensors
using PastaQ


Random.seed!(1234)

N = 17
depth = 5
repetitions = 2

# this is a single block
circuit = randomcircuit(N, depth; twoqubitgates = "CX", onequbitgates = ["Rx","Rz"])

# generate MPO for the block
U = runcircuit(circuit; process = true)

# remember to create the state passing the MPO, so the site indices match
ψ = productstate(U)

for r in 1:repetitions
  # run MPO-MPS contraction (this is a different algorithm
  # and is more precise than contracting each single gate separately :)
  ψ = noprime(U * ψ)

  # Here you should add the measurements.
  # I think you can just create the "measurement circuit" with
  # the correct projectors, run that and that's it.

end
