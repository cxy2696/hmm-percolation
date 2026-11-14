using Test
using Random
include(joinpath(@__DIR__, "..", "HP.jl"))
using .HP

@testset "HP bindings" begin
    @testset "version" begin
        @test occursin(r"^\d+\.\d+\.\d+$", lib_version())
    end

    @testset "uniform HMM" begin
        h = DiscreteHMM(2, 2)
        _, ll = forward(h, [0, 1, 0])
        @test isapprox(ll, log(1 / 8); atol=1e-12)
        path, lp = viterbi(h, [0, 1, 0])
        @test length(path) == 3
        @test isfinite(lp)
    end

    @testset "RNG determinism" begin
        a, b = RNG(123), RNG(123)
        @test [uniform(a) for _ in 1:5] == [uniform(b) for _ in 1:5]
        @test [normal(a) for _ in 1:5] != zeros(5)
    end

    @testset "GHMM recovery" begin
        Random.seed!(3)
        x = vcat(randn(400, 1) .* 0.6 .- 2.0, randn(400, 1) .* 0.6 .+ 2.0)
        m = GaussianHMM(2, 1)
        rep = train!(m, x; max_iter=100, restarts=3, seed=3)
        @test isfinite(rep.loglik)
        @test state_mean(m, 0) < -1.5
        @test state_mean(m, 1) > 1.5
        post, ll = posterior(m, x)
        @test size(post) == (2, 800)
        @test all(isapprox.(sum(post; dims=1), 1.0; atol=1e-10))
        @test isapprox(ll, rep.loglik; atol=1e-6)
        save_model(m, joinpath(tempdir(), "hp_jl_model.txt"))
        m2 = load_ghmm(joinpath(tempdir(), "hp_jl_model.txt"))
        @test isapprox(state_mean(m2, 0), state_mean(m, 0); atol=1e-10)
    end

    @testset "lattice" begin
        lat = Lattice(48, :bond)
        bernoulli!(lat, 0.6, RNG(5))
        st = analyze(lat)
        @test st.n_sites == 48 * 48
        @test st.spanning
        @test sum(cluster_sizes(lat)) == 48 * 48
        pc, _ = perc_crossing(32; reps=40, seed=9, iters=10)
        @test abs(pc - 0.5) < 0.08
    end

    @testset "corr and backtest" begin
        n = 20
        corr = [i == j ? 1.0 : ((i <= 10) == (j <= 10) ? 0.8 : 0.0) for i in 1:n, j in 1:n]
        c = corr_critical(corr)
        @test abs(c.theta_half - 0.8) < 0.03
        g = graph_giant(corr, 0.5)
        @test g.n_components == 2
        @test isapprox(g.giant_fraction, 0.5; atol=1e-12)

        ret = [NaN; zeros(4); 0.01; 0.0; 0.0]
        sig = [0, 0, 0, 0, 1, 0, 0, 0]
        bt = backtest(ret, sig; cost_bps=0.0)
        @test isapprox(bt.equity[6], 1.01; atol=1e-12)
        @test isapprox(bt.total_return, 0.01; atol=1e-12)
    end
end
