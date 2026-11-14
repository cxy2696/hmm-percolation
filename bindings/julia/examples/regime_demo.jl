# Run from the repo root: julia bindings/julia/examples/regime_demo.jl
include(joinpath(@__DIR__, "..", "HP.jl"))
using .HP
using Random

Random.seed!(7)
# Two latent regimes: low-vol positive drift and high-vol negative drift.
T = 1200
regime = zeros(Int, T)
r = zeros(T)
state = 1
for t in 1:T
    global state
    if t > 1
        stay = state == 1 ? 0.98 : 0.96
        state = rand() < stay ? state : 3 - state
    end
    regime[t] = state
    mu = state == 1 ? 0.10 / 252 : -0.20 / 252
    sig = state == 1 ? 0.14 / sqrt(252) : 0.34 / sqrt(252)
    r[t] = mu + sig * randn()
end
x = hcat(r .* 100.0)  # one feature column here; the app uses return + vol

m = GaussianHMM(2, 1)
rep = train!(m, x; max_iter=120, restarts=3, seed=42)
println("train: ", rep)
println("state means: ", [round(state_mean(m, j), digits=3) for j in 0:1])

path, lp = viterbi(m, x)
acc = count(path .== (regime .- 1)) / T
println("viterbi accuracy (labels may swap): ", round(max(acc, 1 - acc), digits=3))

# Percolation sanity: 2D bond square lattice, p_c = 0.5.
pc, _ = perc_crossing(64; kind=:bond, reps=80, seed=11, iters=12)
println("bond L=64 crossing: ", round(pc, digits=4))
