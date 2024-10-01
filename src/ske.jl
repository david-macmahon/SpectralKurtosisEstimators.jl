struct SKEstimator
    M::Int
    N::Int
    d::Float64

    u1::Float64
    u2::Float64
    u3::Float64
    u4::Float64
end

function SKEstimator(M::Float64, N::Float64=1.0, d::Real=1)
    isinteger(M) || error("value of M ($M) must be an integer")
    isinteger(N) || error("value of N ($N) must be an integer")
    isinteger(M*N*d) || error("value of M*N*d ($(M*N*d)) must be an integer")

    M >= 2 || error("M ($M) must be >= 2")
    N >= 1 || error("N ($N) must be >= 1")
    d > 0 || error("d ($d) must be > 0")
try
    Nd = N*d
    MNd = M*Nd

    u1=1

    u2=(
        (2M^2 * Nd * (1+Nd))
        /
        ((M-1) * (6 + 5M*Nd + M^2*Nd^2))
    )

    u3=(
        (8M^3 * Nd * (Nd+1))
        * ((-2 + Nd * (-5 + M *(4 + Nd))))
        / ((M-1)^2 * (MNd+5) * (MNd+4) * (MNd+3) * (MNd+2))
    )

    u4=(
        (12M^4 * Nd * (Nd+1))
        * (24 + Nd*(48 + 84Nd + M*(-32 + Nd*(-245 - 93Nd + M*(125 + Nd*(68 + M + (3 + M)*Nd))))))
        /
        ((M-1)^3 * (MNd+7) * (MNd+6) * (MNd+5) * (MNd+4) * (MNd+3) * (MNd+2))
    )

    SKEstimator(M, N, d, u1, u2, u3, u4)
catch
    @show M N d
    rethrow()
end
end

function SKEstimator(M::Integer, N::Integer=1, d::Real=1)
    SKEstimator(float(M), float(N), d)
end

function g1(ske::SKEstimator)
    ske.u3 / sqrt(skw.u2)^3
end

function g2(ske::SKEstimator)
    ske.u4 / ske.u2^2 - 3
end
