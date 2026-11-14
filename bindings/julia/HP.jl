"""
    HP

Julia `ccall` bindings for the C core `libhp` (HMM + percolation + backtest).
Build the shared library first: `make shared` in the project root, or set
`ENV["HP_LIB"]` to the dylib/so path.

Arrays follow C layout at the boundary; Julia-facing results are converted to
ordinary Julia matrices/vectors (states x time for posteriors).
"""
module HP

export HPError, RNG, uniform, normal,
       DiscreteHMM, GaussianHMM, train!, forward, posterior, filtered, viterbi,
       save_model, load_hmm, load_ghmm, params, state_mean,
       Lattice, bernoulli!, analyze, cluster_sizes,
       perc_crossing, perc_spanning_prob, er_giant, bond_threshold, config_giant,
       corr_matrix, corr_critical, graph_giant,
       backtest

const _ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const LIBHP = get(ENV, "HP_LIB",
    Sys.isapple() ? joinpath(_ROOT, "build", "libhp.dylib") :
                    joinpath(_ROOT, "build", "libhp.so"))
isfile(LIBHP) || error("libhp not found at $LIBHP; run `make shared` or set HP_LIB")

struct HPError <: Exception
    status::Cint
    msg::String
end
Base.showerror(io::IO, e::HPError) = print(io, "HPError(", e.status, "): ", e.msg)

function check(st::Integer)
    st == 0 && return nothing
    msg = unsafe_string(ccall((:hp_status_message, LIBHP), Cstring, (Cint,), st))
    throw(HPError(Cint(st), msg))
end

# ------------------------------------------------------------------ RNG
mutable struct RNG
    s::Vector{UInt64}
    function RNG(seed::Integer=42, stream::Integer=0)
        s = zeros(UInt64, 4)
        ccall((:hp_rng_seed_stream, LIBHP), Cvoid,
              (Ptr{UInt64}, UInt64, UInt64), pointer(s), seed, stream)
        new(s)
    end
end
uniform(r::RNG) = ccall((:hp_rng_uniform, LIBHP), Cdouble, (Ptr{UInt64},), pointer(r.s))
normal(r::RNG) = ccall((:hp_rng_normal, LIBHP), Cdouble, (Ptr{UInt64},), pointer(r.s))

# ------------------------------------------------------------ HMM structs
struct HMMData
    n_states::Csize_t
    n_symbols::Csize_t
    pi::Ptr{Cdouble}
    a::Ptr{Cdouble}
    b::Ptr{Cdouble}
end

struct GHMMData
    n_states::Csize_t
    dim::Csize_t
    n_mix::Csize_t
    pi::Ptr{Cdouble}
    a::Ptr{Cdouble}
    mix::Ptr{Cdouble}
    mean::Ptr{Cdouble}
    var::Ptr{Cdouble}
end

struct HMMOpts
    max_iter::Cuint
    tol::Cdouble
    n_restarts::Cuint
    seed::UInt64
    alpha::Cdouble
    init::Cint
    verbose::Cint
end

struct GHMMOpts
    max_iter::Cuint
    tol::Cdouble
    n_restarts::Cuint
    seed::UInt64
    alpha::Cdouble
    var_floor::Cdouble
    mix_floor::Cdouble
    init::Cint
    verbose::Cint
end

struct HMMReport
    loglik::Cdouble
    bic::Cdouble
    aic::Cdouble
    n_params::Csize_t
    iterations::Cuint
    converged::Cint
    restarts::Cuint
end

const INIT_QUANTILE = Cint(0)
const INIT_RANDOM = Cint(1)

# ------------------------------------------------------------------ HMM
mutable struct DiscreteHMM
    ptr::Ptr{Cvoid}
    function DiscreteHMM(n_states::Integer, n_symbols::Integer)
        p = ccall((:hp_hmm_new, LIBHP), Ptr{Cvoid}, (Csize_t, Csize_t),
                  n_states, n_symbols)
        p == C_NULL && error("hp_hmm_new failed")
        obj = new(p)
        finalizer(o -> ccall((:hp_hmm_free, LIBHP), Cvoid, (Ptr{Cvoid},), o.ptr), obj)
        obj
    end
end

mutable struct GaussianHMM
    ptr::Ptr{Cvoid}
    function GaussianHMM(n_states::Integer, dim::Integer, n_mix::Integer=1)
        p = ccall((:hp_ghmm_new, LIBHP), Ptr{Cvoid}, (Csize_t, Csize_t, Csize_t),
                  n_states, dim, n_mix)
        p == C_NULL && error("hp_ghmm_new failed")
        obj = new(p)
        finalizer(o -> ccall((:hp_ghmm_free, LIBHP), Cvoid, (Ptr{Cvoid},), o.ptr), obj)
        obj
    end
end

function Base.close(h::DiscreteHMM)
    h.ptr == C_NULL && return
    ccall((:hp_hmm_free, LIBHP), Cvoid, (Ptr{Cvoid},), h.ptr)
    h.ptr = C_NULL
    nothing
end

function Base.close(m::GaussianHMM)
    m.ptr == C_NULL && return
    ccall((:hp_ghmm_free, LIBHP), Cvoid, (Ptr{Cvoid},), m.ptr)
    m.ptr = C_NULL
    nothing
end

hmm_data(h::DiscreteHMM) = unsafe_load(Ptr{HMMData}(h.ptr))
ghmm_data(m::GaussianHMM) = unsafe_load(Ptr{GHMMData}(m.ptr))

function params(h::DiscreteHMM)
    c = hmm_data(h)
    n, m = Int(c.n_states), Int(c.n_symbols)
    pi = copy(unsafe_wrap(Array, c.pi, n))
    a = permutedims(reshape(copy(unsafe_wrap(Array, c.a, n * n)), n, n))
    b = permutedims(reshape(copy(unsafe_wrap(Array, c.b, n * m)), m, n))
    (pi=pi, a=a, b=b)
end

function params(m::GaussianHMM)
    c = ghmm_data(m)
    n, d, k = Int(c.n_states), Int(c.dim), Int(c.n_mix)
    pi = copy(unsafe_wrap(Array, c.pi, n))
    a = permutedims(reshape(copy(unsafe_wrap(Array, c.a, n * n)), n, n))
    mix = permutedims(reshape(copy(unsafe_wrap(Array, c.mix, n * k)), k, n))
    mean = permutedims(reshape(copy(unsafe_wrap(Array, c.mean, n * k * d)), d, k, n), (3, 2, 1))
    var = permutedims(reshape(copy(unsafe_wrap(Array, c.var, n * k * d)), d, k, n), (3, 2, 1))
    (pi=pi, a=a, mix=mix, mean=mean, var=var)
end

state_mean(m::GaussianHMM, state::Integer, dim::Integer=1) = begin
    p = params(m)
    sum(p.mix[state + 1, c] * p.mean[state + 1, c, dim] for c in 1:Int(size(p.mix, 2)))
end

function train!(h::DiscreteHMM, obs::AbstractVector{<:Integer};
                max_iter::Integer=200, tol::Real=1e-6, restarts::Integer=4,
                seed::Integer=42, alpha::Real=1.0, init::Integer=INIT_QUANTILE,
                verbose::Integer=0)
    o = Vector{UInt32}(obs)
    opts = HMMOpts(max_iter, tol, restarts, seed, alpha, Cint(init), Cint(verbose))
    rep = Ref{HMMReport}()
    c = hmm_data(h)
    check(ccall((:hp_hmm_train, LIBHP), Cint,
        (Ptr{UInt32}, Csize_t, Csize_t, Csize_t, Ref{HMMOpts}, Ptr{Cvoid}, Ref{HMMReport}),
        o, length(o), Int(c.n_states), Int(c.n_symbols), opts, h.ptr, rep))
    (loglik=rep[].loglik, bic=rep[].bic, aic=rep[].aic,
     n_params=Int(rep[].n_params), iterations=Int(rep[].iterations),
     converged=rep[].converged != 0, restarts=Int(rep[].restarts))
end

function train!(m::GaussianHMM, x::AbstractMatrix{<:Real};
                max_iter::Integer=200, tol::Real=1e-6, restarts::Integer=4,
                seed::Integer=42, alpha::Real=1.0, var_floor::Real=1e-4,
                mix_floor::Real=1e-8, init::Integer=INIT_QUANTILE,
                verbose::Integer=0)
    T, d = size(x)
    c = ghmm_data(m)
    Int(c.dim) == d || error("expected dim $(Int(c.dim)), got $d")
    xf = Vector{Float64}(vec(permutedims(Matrix{Float64}(x))))
    opts = GHMMOpts(max_iter, tol, restarts, seed, alpha, var_floor, mix_floor,
                    Cint(init), Cint(verbose))
    rep = Ref{HMMReport}()
    check(ccall((:hp_ghmm_train, LIBHP), Cint,
        (Ptr{Cdouble}, Csize_t, Csize_t, Csize_t, Csize_t, Ref{GHMMOpts},
         Ptr{Cvoid}, Ref{HMMReport}),
        xf, T, d, Int(c.n_states), Int(c.n_mix), opts, m.ptr, rep))
    (loglik=rep[].loglik, bic=rep[].bic, aic=rep[].aic,
     n_params=Int(rep[].n_params), iterations=Int(rep[].iterations),
     converged=rep[].converged != 0, restarts=Int(rep[].restarts))
end

function forward(h::DiscreteHMM, obs::AbstractVector{<:Integer})
    o = Vector{UInt32}(obs)
    n = Int(hmm_data(h).n_states)
    T = length(o)
    la = Vector{Float64}(undef, T * n)
    ll = Ref{Cdouble}(0.0)
    check(ccall((:hp_hmm_forward_log, LIBHP), Cint,
        (Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{Cdouble}, Ref{Cdouble}),
        h.ptr, o, T, la, ll))
    (reshape(la, n, T), ll[])
end

function viterbi(h::DiscreteHMM, obs::AbstractVector{<:Integer})
    o = Vector{UInt32}(obs)
    T = length(o)
    path = Vector{UInt32}(undef, T)
    lp = Ref{Cdouble}(0.0)
    check(ccall((:hp_hmm_viterbi, LIBHP), Cint,
        (Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{UInt32}, Ref{Cdouble}),
        h.ptr, o, T, path, lp))
    (Int.(path), lp[])
end

function posterior(h::DiscreteHMM, obs::AbstractVector{<:Integer})
    o = Vector{UInt32}(obs)
    n = Int(hmm_data(h).n_states)
    T = length(o)
    g = Vector{Float64}(undef, T * n)
    ll = Ref{Cdouble}(0.0)
    check(ccall((:hp_hmm_posterior, LIBHP), Cint,
        (Ptr{Cvoid}, Ptr{UInt32}, Csize_t, Ptr{Cdouble}, Ref{Cdouble}),
        h.ptr, o, T, g, ll))
    (reshape(g, n, T), ll[])
end

function forward(m::GaussianHMM, x::AbstractMatrix{<:Real})
    T, d = size(x)
    c = ghmm_data(m)
    Int(c.dim) == d || error("expected dim $(Int(c.dim)), got $d")
    xf = Vector{Float64}(vec(permutedims(Matrix{Float64}(x))))
    n = Int(c.n_states)
    la = Vector{Float64}(undef, T * n)
    ll = Ref{Cdouble}(0.0)
    check(ccall((:hp_ghmm_forward_log, LIBHP), Cint,
        (Ptr{Cvoid}, Ptr{Cdouble}, Csize_t, Ptr{Cdouble}, Ref{Cdouble}),
        m.ptr, xf, T, la, ll))
    (reshape(la, n, T), ll[])
end

function posterior(m::GaussianHMM, x::AbstractMatrix{<:Real})
    T, d = size(x)
    c = ghmm_data(m)
    Int(c.dim) == d || error("expected dim $(Int(c.dim)), got $d")
    xf = Vector{Float64}(vec(permutedims(Matrix{Float64}(x))))
    n = Int(c.n_states)
    g = Vector{Float64}(undef, T * n)
    ll = Ref{Cdouble}(0.0)
    check(ccall((:hp_ghmm_posterior, LIBHP), Cint,
        (Ptr{Cvoid}, Ptr{Cdouble}, Csize_t, Ptr{Cdouble}, Ref{Cdouble}),
        m.ptr, xf, T, g, ll))
    (reshape(g, n, T), ll[])
end

function filtered(m::GaussianHMM, x::AbstractMatrix{<:Real})
    T, d = size(x)
    c = ghmm_data(m)
    Int(c.dim) == d || error("expected dim $(Int(c.dim)), got $d")
    xf = Vector{Float64}(vec(permutedims(Matrix{Float64}(x))))
    n = Int(c.n_states)
    f = Vector{Float64}(undef, T * n)
    ll = Ref{Cdouble}(0.0)
    check(ccall((:hp_ghmm_filter, LIBHP), Cint,
        (Ptr{Cvoid}, Ptr{Cdouble}, Csize_t, Ptr{Cdouble}, Ref{Cdouble}),
        m.ptr, xf, T, f, ll))
    (reshape(f, n, T), ll[])
end

function viterbi(m::GaussianHMM, x::AbstractMatrix{<:Real})
    T, d = size(x)
    c = ghmm_data(m)
    Int(c.dim) == d || error("expected dim $(Int(c.dim)), got $d")
    xf = Vector{Float64}(vec(permutedims(Matrix{Float64}(x))))
    path = Vector{UInt32}(undef, T)
    lp = Ref{Cdouble}(0.0)
    check(ccall((:hp_ghmm_viterbi, LIBHP), Cint,
        (Ptr{Cvoid}, Ptr{Cdouble}, Csize_t, Ptr{UInt32}, Ref{Cdouble}),
        m.ptr, xf, T, path, lp))
    (Int.(path), lp[])
end

save_model(h::DiscreteHMM, path::AbstractString) =
    check(ccall((:hp_hmm_save, LIBHP), Cint, (Ptr{Cvoid}, Cstring), h.ptr, path))
save_model(m::GaussianHMM, path::AbstractString) =
    check(ccall((:hp_ghmm_save, LIBHP), Cint, (Ptr{Cvoid}, Cstring), m.ptr, path))

function load_hmm(path::AbstractString)
    p = ccall((:hp_hmm_load, LIBHP), Ptr{Cvoid}, (Cstring,), path)
    p == C_NULL && error("cannot load HMM from $path")
    h = DiscreteHMM(2, 2)  # registers the single finalizer
    close(h)               # release the placeholder
    h.ptr = p
    h
end

function load_ghmm(path::AbstractString)
    p = ccall((:hp_ghmm_load, LIBHP), Ptr{Cvoid}, (Cstring,), path)
    p == C_NULL && error("cannot load GHMM from $path")
    m = GaussianHMM(2, 1, 1)  # registers the single finalizer
    close(m)                  # release the placeholder
    m.ptr = p
    m
end

# -------------------------------------------------------------- lattice
struct LatticeStatsC
    n_sites::Csize_t
    n_occupied::Csize_t
    n_clusters::Csize_t
    max_cluster::Csize_t
    s_max::Cdouble
    susceptibility::Cdouble
    occupancy::Cdouble
    spanning_x::Cint
    spanning_y::Cint
    spanning::Cint
end

mutable struct Lattice
    ptr::Ptr{Cvoid}
    L::Int
    kind::Symbol
    periodic::Bool
    function Lattice(L::Integer, kind::Symbol=:bond, periodic::Bool=false)
        kind in (:bond, :site) || error("kind must be :bond or :site")
        k = kind == :bond ? 1 : 0
        p = ccall((:hp_lattice_new, LIBHP), Ptr{Cvoid},
                  (Cint, Cint, Cint), L, k, periodic ? 1 : 0)
        p == C_NULL && error("bad lattice")
        obj = new(p, Int(L), kind, periodic)
        finalizer(o -> ccall((:hp_lattice_free, LIBHP), Cvoid, (Ptr{Cvoid},), o.ptr), obj)
        obj
    end
end

function bernoulli!(lat::Lattice, p::Real, rng::RNG)
    check(ccall((:hp_lattice_bernoulli, LIBHP), Cint,
        (Ptr{Cvoid}, Cdouble, Ptr{UInt64}), lat.ptr, p, pointer(rng.s)))
    lat
end

function analyze(lat::Lattice)
    st = Ref{LatticeStatsC}()
    check(ccall((:hp_lattice_analyze, LIBHP), Cint,
        (Ptr{Cvoid}, Ref{LatticeStatsC}), lat.ptr, st))
    s = st[]
    (n_sites=Int(s.n_sites), n_occupied=Int(s.n_occupied),
     n_clusters=Int(s.n_clusters), max_cluster=Int(s.max_cluster),
     s_max=s.s_max, susceptibility=s.susceptibility, occupancy=s.occupancy,
     spanning=s.spanning != 0, spanning_x=s.spanning_x != 0,
     spanning_y=s.spanning_y != 0)
end

function cluster_sizes(lat::Lattice)
    n = Ref{Csize_t}(0)
    ptr = ccall((:hp_lattice_cluster_sizes, LIBHP), Ptr{Csize_t},
                (Ptr{Cvoid}, Ref{Csize_t}), lat.ptr, n)
    ptr == C_NULL && return Int[]
    out = copy(unsafe_wrap(Array, ptr, Int(n[])))
    ccall((:hp_free, LIBHP), Cvoid, (Ptr{Cvoid},), ptr)
    out
end

struct PercExpC
    L::Cint
    kind::Cint
    periodic::Cint
    reps::Csize_t
    seed::UInt64
end

function perc_crossing(L::Integer; kind::Symbol=:bond, periodic::Bool=false,
                       reps::Integer=100, seed::Integer=1,
                       p_lo::Real=0.35, p_hi::Real=0.65, iters::Integer=14)
    exp = PercExpC(L, kind == :bond ? 1 : 0, periodic ? 1 : 0, reps, seed)
    pc = Ref{Cdouble}(0.0)
    sp = Ref{Cdouble}(0.0)
    check(ccall((:hp_perc_crossing, LIBHP), Cint,
        (Ref{PercExpC}, Cdouble, Cdouble, Cuint, Ref{Cdouble}, Ref{Cdouble}),
        exp, p_lo, p_hi, iters, pc, sp))
    (pc[], sp[])
end

function perc_spanning_prob(L::Integer, p::Real; kind::Symbol=:bond,
                            periodic::Bool=false, reps::Integer=100,
                            seed::Integer=1)
    exp = PercExpC(L, kind == :bond ? 1 : 0, periodic ? 1 : 0, reps, seed)
    ccall((:hp_perc_spanning_prob, LIBHP), Cdouble, (Ref{PercExpC}, Cdouble), exp, p)
end

er_giant(c::Real) = ccall((:hp_perc_er_giant_theory, LIBHP), Cdouble, (Cdouble,), c)
bond_threshold(mu::Real) = ccall((:hp_perc_bond_threshold, LIBHP), Cdouble, (Cdouble,), mu)

function config_giant(pk::AbstractVector{<:Real}; transmissibility::Real=1.0)
    g = Ref{Cdouble}(0.0)
    r0 = Ref{Cdouble}(0.0)
    mean = Ref{Cdouble}(0.0)
    check(ccall((:hp_perc_config_giant, LIBHP), Cint,
        (Ptr{Cdouble}, Csize_t, Cdouble, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
        Vector{Float64}(pk), length(pk) - 1, transmissibility, g, r0, mean))
    (giant=g[], r0=r0[], mean_degree=mean[])
end

# ------------------------------------------------------- correlation
function corr_matrix(x::AbstractMatrix{<:Real})
    T, N = size(x)
    xf = Vector{Float64}(vec(permutedims(Matrix{Float64}(x))))
    raw = Vector{Float64}(undef, N * N)
    ccall((:hp_corr_matrix, LIBHP), Cvoid,
          (Ptr{Cdouble}, Csize_t, Csize_t, Ptr{Cdouble}), xf, T, N, raw)
    permutedims(reshape(raw, N, N))
end

struct CorrCriticalC
    theta_c::Cdouble
    theta_half::Cdouble
    chi_max::Cdouble
    giant_at_c::Cdouble
    censored_high::Cint
    censored_low::Cint
end

function corr_critical(corr::AbstractMatrix{<:Real}; theta_min::Real=0.1,
                       theta_max::Real=0.99, steps::Integer=45)
    n = size(corr, 1)
    cf = Vector{Float64}(vec(permutedims(Matrix{Float64}(corr))))
    c = Ref{CorrCriticalC}()
    check(ccall((:hp_corr_critical, LIBHP), Cint,
        (Ptr{Cdouble}, Csize_t, Cdouble, Cdouble, Csize_t, Ref{CorrCriticalC}),
        cf, n, theta_min, theta_max, steps, c))
    (theta_peak=c[].theta_c, theta_half=c[].theta_half, chi_max=c[].chi_max,
     giant_at_peak=c[].giant_at_c, censored_high=c[].censored_high != 0,
     censored_low=c[].censored_low != 0)
end

struct ComponentsC
    n::Csize_t
    n_components::Csize_t
    giant::Csize_t
    giant_fraction::Cdouble
    susceptibility::Cdouble
    labels::Ptr{Csize_t}
    sizes::Ptr{Csize_t}
end

function graph_giant(corr::AbstractMatrix{<:Real}, theta::Real)
    n = size(corr, 1)
    cf = Vector{Float64}(vec(permutedims(Matrix{Float64}(corr))))
    g = ccall((:hp_graph_from_corr, LIBHP), Ptr{Cvoid},
              (Ptr{Cdouble}, Csize_t, Cdouble), cf, n, theta)
    g == C_NULL && error("graph_from_corr failed")
    c = Ref{ComponentsC}()
    try
        check(ccall((:hp_graph_components, LIBHP), Cint,
            (Ptr{Cvoid}, Cdouble, Ref{ComponentsC}), g, -Inf, c))
        sizes = copy(unsafe_wrap(Array, c[].sizes, Int(c[].n_components)))
        (n_edges=Int(ccall((:hp_graph_edge_count, LIBHP), Csize_t, (Ptr{Cvoid},), g)),
         n_components=Int(c[].n_components), giant=Int(c[].giant),
         giant_fraction=c[].giant_fraction, susceptibility=c[].susceptibility,
         sizes=sizes)
    finally
        ccall((:hp_components_free, LIBHP), Cvoid, (Ref{ComponentsC},), c)
        ccall((:hp_graph_free, LIBHP), Cvoid, (Ptr{Cvoid},), g)
    end
end

# ----------------------------------------------------------- backtest
struct BTOptsC
    cost_bps::Cdouble
    apply_costs::Cint
end

struct BTMetricsC
    total_return::Cdouble
    cagr::Cdouble
    ann_vol::Cdouble
    sharpe::Cdouble
    max_drawdown::Cdouble
    calmar::Cdouble
    avg_turnover::Cdouble
    exposure::Cdouble
    hit_rate::Cdouble
    final_equity::Cdouble
    n_days::Csize_t
end

struct BTResultC
    equity::Ptr{Cdouble}
    strategy_ret::Ptr{Cdouble}
    position::Ptr{Cdouble}
end

function backtest(asset_ret::AbstractVector{<:Real}, signal::AbstractVector{<:Real};
                  cost_bps::Real=5.0, apply_costs::Bool=true)
    T = length(asset_ret)
    length(signal) == T || error("signal length mismatch")
    r = Vector{Float64}(asset_ret)
    s = Vector{Float64}(signal)
    opts = BTOptsC(cost_bps, apply_costs ? 1 : 0)
    res = Ref{BTResultC}()
    m = Ref{BTMetricsC}()
    check(ccall((:hp_backtest, LIBHP), Cint,
        (Ptr{Cdouble}, Csize_t, Ptr{Cdouble}, Ref{BTOptsC},
         Ref{BTResultC}, Ref{BTMetricsC}), r, T, s, opts, res, m))
    out = (total_return=m[].total_return, cagr=m[].cagr, ann_vol=m[].ann_vol,
           sharpe=m[].sharpe, max_drawdown=m[].max_drawdown, calmar=m[].calmar,
           avg_turnover=m[].avg_turnover, exposure=m[].exposure,
           hit_rate=m[].hit_rate, final_equity=m[].final_equity,
           n_days=Int(m[].n_days),
           equity=copy(unsafe_wrap(Array, res[].equity, T)),
           strategy_ret=copy(unsafe_wrap(Array, res[].strategy_ret, T)),
           position=copy(unsafe_wrap(Array, res[].position, T)))
    ccall((:hp_bt_result_free, LIBHP), Cvoid, (Ref{BTResultC},), res)
    out
end

end # module
