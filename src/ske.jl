struct SKEstimator
    M::Int
    N::Int
    d::Float64

    u1::Float64
    u2::Float64
    u3::Float64
    u4::Float64
end

function SKEstimator(M::Integer, N::Integer=1, d::AbstractFloat=1)
    M >= 2 || error("M ($M) must be > 2")
    N >= 1 || error("N ($N) must be > 1")
    d > 0 || error("d ($d) must be > 0")
try
    Nd = N*d
    MNd = M*Nd

    log_Nd = log(Nd)
    log_Nd1 = log(Nd+1)
    log_M = log(M)
    log_M1 = log(M-1)
    log_MNd = log(M)+log(Nd)

    loggammadiv_MNd24 = MNd < 6 ? loggamma(MNd+2)-loggamma(MNd+4) : loggammadiv(2, MNd+2)
    loggammadiv_MNd26 = MNd < 6 ? loggamma(MNd+2)-loggamma(MNd+6) : loggammadiv(4, MNd+2)
    loggammadiv_MNd28 = MNd < 6 ? loggamma(MNd+2)-loggamma(MNd+8) : loggammadiv(6, MNd+2)

    u1 = 1

    u2 = exp(
        log(2) + log_Nd + log_Nd1 + 2log_M + loggammadiv_MNd24
        -
        log_M1
    )

    u3 = exp(
        log(8) + log_Nd + log_Nd1 + 3log_M + loggammadiv_MNd26
        -
        2log_M1
        +
        log((Nd+4)*MNd - 5Nd - 2)
    )

    #=
    u4 = exp(
        log(12) + log_Nd + log_Nd1 + 4log_M + loggammadiv_MNd28
        -
        3log_M1
        +
        log(
            M^3*Nd^4 + 3M^2*Nd^4 + MNd^3 + 68M^2*Nd^3
            - 93M*Nd^3 + 125MNd^2 - 245M*Nd^2
            + 84Nd^2 - 32MNd + 48Nd + 24
        )
    )
    =#

    logu4q = (
        log(12) + log_Nd + log_Nd1 + 4log_M + loggammadiv_MNd28
        -
        3log_M1
    )

    u4 = mapreduce(((s,x),)->s*exp(logu4q+x), +,
        (
            (+1, 3log_M+4log_Nd), (+1, log(3)+2log_M+4log_Nd), (+1, 3log_MNd), (+1, log(68)+2log_M+3log_Nd),
            (-1, log(93)+log_M+3log_Nd), (+1, log(125)+2log_MNd), (-1, log(245)+log_M+2log_Nd),
            (+1, log(84)+2log_Nd), (-1, log(32)+log_MNd), (+1, log(48)+log_Nd), (+1, log(24))
        )
    )

    SKEstimator(intM, intN, d, u1, u2, u3, u4)
catch
    @show M N d
    rethrow()
end    
end    

function B1(u2, u3)
    u3^2 / u2^3
end

function B1(ske::SKEstimator)
    B1(ske.u2, ske.u3)
end

function B2(u2, u4)
    u4 / u2^2
end

function B2(ske::SKEstimator)
    B2(ske.u2, ske.u4)
end

function g1(ske::SKEstimator)
    sqrt(B1(ske))
end

function g2(ske::SKEstimator)
    B2(ske) - 3
end
