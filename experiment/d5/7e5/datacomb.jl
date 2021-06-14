using DelimitedFiles
#using StatsBase

function comb()

dx = 3
dz = 5
p = "0.00007"
num_sim = 4

nfails = []
ntrials = []
failure_rate = []


failtot=0
trialtot=0


for sim_id in 1:num_sim
  fname = "experiment_dz$(dz)_dx$(dx)_p$(p)_id$(sim_id).txt"
  fin = open(fname,"r")
  data = readdlm(fin)
  n = Int(data[1])
  f = Int(data[2])
  failtot = Int(failtot+f)
  trialtot = Int(trialtot+n)
  push!(failure_rate, f/n)
  close(fin)
end

fname="combdata_dz$(dz)_dx$(dx)_p$(p)_id.txt"
fout=open(fname,"w")
writedlm(fout,[failtot,trialtot,failtot/trialtot])
close(fout)

end

comb()

#@show mean(failure_rate)
#@show sem(failure rate

