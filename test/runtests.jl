using Test
using Random
using Statistics
using Distributions
using SpectralKurtosisEstimators
using SpectralKurtosisEstimators: relative_error, thresholds, skewness

# Unit tests reconstruct some of the entries of Table 1 of
# [Nita [2010b]](https://doi.org/10.1111/j.1745-3933.2010.00882.x)

M = 300
N = 10
d = 1
nsigma = 3

lowervi_expected = 0.766_48
uppervi_expected = 1.283_13
# Table 1 has 0.18% of u4errvi_expected, but our errors are not percent
u4errvi_expected = -0.001_8

loweriii_expected = 0.767_54
upperiii_expected = 1.282_12
u4erriii_expected = -0.007_1

loweriv_expected = 0.766_13
upperiv_expected = 1.283_45

# Create SKEstimator for the given paramters
ske = SKEstimator(M, N, d)

# Get PearsonTypeVI distribution for our SKEsitmator along with a measure of the
# error in the fourth central moment.
pdvi = PearsonTypeVI(ske)
u4errvi = relative_error(pdvi, ske)

# Get PearsonTypeIII distribution for our SKEsitmator along with a measure of
# the error in the fourth central moment.
pdiii = PearsonTypeIII(ske)
u4erriii = relative_error(pdiii, ske)

# Get PearsonTypeIV distribution for our SKEsitmator
pdiv = PearsonTypeIV(ske)

# Compute lower and upper thresholds for our "probability curves" corresponding
# to ±3 sigma of a standard normal distribution.
lowervi, uppervi = thresholds(pdvi, nsigma)
loweriii, upperiii = thresholds(pdiii, nsigma)
loweriv, upperiv = thresholds(pdiv, nsigma)

# Test that everything agrees with the Table 1 entries to the given precision
@testset "PearsonTypeVI " begin
    @test lowervi ≈ lowervi_expected atol=0.000_005
    @test uppervi ≈ uppervi_expected atol=0.000_005
    @test u4errvi ≈ u4errvi_expected atol=0.000_05
end

@testset "PearsonTypeIII" begin
    @test loweriii ≈ loweriii_expected atol=0.000_005
    @test upperiii ≈ upperiii_expected atol=0.000_005
    @test u4erriii ≈ u4erriii_expected atol=0.000_05
end

 getm(::PearsonTypeIV{m,ν,a,λ}) where {m,ν,a,λ} = m
 getν(::PearsonTypeIV{m,ν,a,λ}) where {m,ν,a,λ} = ν
 geta(::PearsonTypeIV{m,ν,a,λ}) where {m,ν,a,λ} = a
 getλ(::PearsonTypeIV{m,ν,a,λ}) where {m,ν,a,λ} = λ

@testset "PearsonTypeIV " begin
    @test loweriv ≈ loweriv_expected atol=0.000_005
    @test upperiv ≈ upperiv_expected atol=0.000_005

    # Test calculated PearsonTypeIV fields for various values of M (and N=1,
    # d=1) using examples from figure 5 of:
    # [Nita [2010a]](https://doi.org/10.1086/652409)
    @testset "M=$M" for (M,  expected) in (
        (32,   (; m=  5.760, ν= -21.847, a=0.389, λ=0.108)),
        (512,  (; m= 20.125, ν= -32.442, a=0.410, λ=0.652)),
        (4096, (; m=132.094, ν=-211.441, a=0.393, λ=0.683)),
        (8192, (; m=260.092, ν=-416.381, a=0.392, λ=0.685)),
    )
        piv = PearsonTypeIV(SKEstimator(M))
        @test getm(piv) ≈ expected.m atol=0.000_5
        @test getν(piv) ≈ expected.ν atol=0.000_5
        @test geta(piv) ≈ expected.a atol=0.000_5
        @test getλ(piv) ≈ expected.λ atol=0.000_5
    end
end

# -----------------------------------------------------------------------------
# Tests for the extension to non-integer M*N*d
#
# The moment formulas of SKEstimator are a closed-form evaluation of the
# gamma-function ratios of equation 9 of Nita & Gary (2010), which are valid
# for any real MNd > 0, not just integer MNd.  The tests below verify this
# against exact rational arithmetic and against Monte Carlo simulation.
# -----------------------------------------------------------------------------

# Exact (Rational{BigInt}) evaluation of the moment formulas, used as an
# infinite-precision reference: for rational inputs every operation below is
# exact, so any transcription error in the formulas is caught without any
# floating point oracle
function reference_moments(M::Real, Nd::Real, MNd::Real)
    M = Rational{BigInt}(M)
    Nd = Rational{BigInt}(Nd)
    MNd = Rational{BigInt}(MNd)
    u2 = (2M^2 * Nd * (1+Nd)) /
         ((M-1) * (6 + 5M*Nd + M^2*Nd^2))
    u3 = (8M^3 * Nd * (Nd+1)) * (-2 + Nd * (-5 + M * (4 + Nd))) /
         ((M-1)^2 * (MNd+5) * (MNd+4) * (MNd+3) * (MNd+2))
    u4 = (12M^4 * Nd * (Nd+1)) *
         (24 + Nd*(48 + 84Nd + M*(-32 + Nd*(-245 - 93Nd + M*(125 + Nd*(68 + M + (3 + M)*Nd)))))) /
         ((M-1)^3 * (MNd+7) * (MNd+6) * (MNd+5) * (MNd+4) * (MNd+3) * (MNd+2))
    (u2, u3, u4)
end

@testset "non-integer MNd" begin
    @testset "golden values" begin
        # With rational inputs, the moment formulas involve only rational
        # operations, so the expected moments below are exact rationals
        @testset "M=$M, N=$N, d=$d (MNd=$(M*N*d))" for (M, N, d, expected) in (
            (3, 1, 0.5,  (u2 = 3//7,       u3 = 162//1001,     u4 = 7209//17017)),
            (3, 1, 1.5,  (u2 = 9//13,      u3 = 3294//4199,    u4 = 232065//96577)),
            (7, 2, 0.75, (u2 = 49//135,    u3 = 132398//364095, u4 = 4126633//4005045)),
        )
            ske = SKEstimator(M, N, d)
            @test ske.u1 == 1
            @test ske.u2 ≈ expected.u2 rtol=1e-12
            @test ske.u3 ≈ expected.u3 rtol=1e-12
            @test ske.u4 ≈ expected.u4 rtol=1e-12
        end

        # Generic (non half-integer) non-integer MNd, expected values computed
        # with the product forms at high (BigFloat) precision
        ske = SKEstimator(5, 2, 1.23)  # MNd = 12.3
        @test ske.u1 == 1
        @test ske.u2 ≈ 0.4862882215823392 rtol=1e-12
        @test ske.u3 ≈ 0.5618207727554533 rtol=1e-12
        @test ske.u4 ≈ 1.8001004095621955 rtol=1e-12
    end

    @testset "sweep vs exact reference" begin
        # Dense sweep over d (hitting integer, half-integer, and quarter-integer
        # MNd values, plus random points) must agree with the exact rational
        # evaluation of the moment formulas everywhere.  The reference is
        # computed at exactly the Float64 values the package sees (converted to
        # exact rationals), so the comparison isolates floating point error.
        ds = vcat(0.05:0.001:2.1, 0.05 .+ 2 .* rand(MersenneTwister(42), 100))
        @testset "M=$M, N=$N" for (M, N) in ((3, 1), (2, 3), (6, 2))
            for d in ds
                ske = SKEstimator(M, N, d)
                u2r, u3r, u4r = reference_moments(M, N*d, M*N*d)
                @test ske.u2 ≈ Float64(u2r) rtol=1e-12
                @test ske.u3 ≈ Float64(u3r) rtol=1e-12
                @test ske.u4 ≈ Float64(u4r) rtol=1e-12
            end
        end
    end

    @testset "guard rails" begin
        @test_throws ErrorException SKEstimator(2.5, 1.0, 1.0)  # non-integer M
        @test_throws ErrorException SKEstimator(2.0, 1.5, 1.0)  # non-integer N
        @test_throws ErrorException SKEstimator(1, 1, 1)        # M < 2
        @test_throws ErrorException SKEstimator(2, 0, 1)        # N < 1
        @test_throws ErrorException SKEstimator(2, 1, 0)        # d <= 0
        @test_throws ErrorException SKEstimator(2, 1, -0.5)     # d <= 0

        # Non-integer M*N*d is now allowed
        ske = SKEstimator(3, 1, 0.5)
        @test ske.M == 3 && ske.N == 1 && ske.d == 0.5

        # Large MNd stays finite
        largeske = SKEstimator(1e6, 1, 1)
        @test isfinite(largeske.u2) && isfinite(largeske.u3) && isfinite(largeske.u4)
        @test largeske.u2 > 0
    end
end

@testset "Rational path" begin
    # Rational d produces SKEstimator{Rational{BigInt}} with exact moments;
    # everything else (including plain Integer d) stays on the Float64 path
    @testset "eltype selection" begin
        @test SKEstimator(3, 2, 1//3) isa SKEstimator{Rational{BigInt}}
        @test SKEstimator(3, 2, 0.5) isa SKEstimator{Float64}
        @test SKEstimator(300, 10, 1) isa SKEstimator{Float64}
        @test SKEstimator(M=3, N=2, d=1//3) isa SKEstimator{Rational{BigInt}}
    end

    @testset "exact moments" begin
        # Same parameters as the Float64 golden case (3, 1, 0.5), now exact
        ske = SKEstimator(3, 1, 1//2)
        @test ske.u1 == 1//1
        @test ske.u2 == 3//7
        @test ske.u3 == 162//1001
        @test ske.u4 == 7209//17017
        @test ske.d == 1//2

        # Hand-derived exact values
        ske = SKEstimator(3, 2, 1//3)
        @test ske.u1 == 1//1
        @test ske.u2 == 1//2
        @test ske.u3 == 2//7
        @test ske.u4 == 5//7
        @test mean(ske) == 1//1
        @test var(ske) == 1//2
        @test kurtosis(ske) == -1//7
    end

    @testset "exact at flagship parameters" begin
        # The case that rules out Rational{Int64} storage: u4's reduced exact
        # form needs 20+23 digits.  Must construct without overflow and equal
        # the exact reference exactly
        ske = SKEstimator(300, 10, 1//1)
        u2r, u3r, u4r = reference_moments(300, 10//1, 3000//1)
        @test ske.u1 == 1//1
        @test ske.u2 == u2r
        @test ske.u3 == u3r
        @test ske.u4 == u4r
    end

    @testset "exact Pearson criterion" begin
        ske = SKEstimator(300, 10, 1//1)
        κ = pearson_criterion(ske)
        @test κ isa Rational{BigInt}
        @test 0 < κ < 1  # exact Type IV region classification
    end

    @testset "skhat with Rational ske" begin
        ske = SKEstimator(3, 2, 1//3)
        s1 = [1.0, 2.0, 3.0]
        s2 = [1.5, 6.0, 12.0]
        sk = skhat(s1, s2, ske)
        # Rational{BigInt} parameters promote Float64 data to BigFloat —
        # Julia promotion semantics preserve maximum precision
        @test sk isa AbstractVector{BigFloat}
        @test all(isfinite, sk)
    end

    @testset "guard rails" begin
        @test_throws ErrorException SKEstimator(7//2, 1, 1)    # non-integer M
        @test_throws ErrorException SKEstimator(3, 3//2, 1//1) # non-integer N
        @test_throws ErrorException SKEstimator(3, 2, -1//3)   # d <= 0
        @test_throws ErrorException SKEstimator(1//1, 1, 1//1) # M < 2
    end
end

@testset "pearson_distribution errors" begin
    # κ < 0 (the Pearson Type I region) is not supported and throws
    ske = SKEstimator(3, 1, 0.5)
    @test pearson_criterion(ske) < 0
    @test_throws ErrorException pearson_distribution(ske)
end

@testset "non-integer MNd Monte Carlo" begin
    # For per-sample powers with shape d, the accumulated sample of N addends
    # is Gamma(N*d) distributed for any real shape, so sampling Gamma(N*d)
    # directly exercises the non-integer M*N*d extension against real data
    rng = MersenneTwister(20260925)
    T = 1_000_000
    @testset "M=$M, N=$N, d=$d (MNd=$(M*N*d))" for (M, N, d, high) in (
        (3, 1, 0.5, false),  # MNd = 1.5 (extreme): mean/var only
        (7, 1, 1.5, true),   # MNd = 10.5
        (6, 2, 1.1, true),   # MNd = 13.2
    )
        ske = SKEstimator(M, N, d)
        xs = rand(rng, Gamma(N*d, 1.0), T, M)
        s1 = vec(sum(xs; dims=2))
        s2 = vec(sum(abs2, xs; dims=2))
        sk = skhat(s1, s2, M, N, d)
        @test mean(sk) ≈ ske.u1 rtol=0.01
        @test var(sk) ≈ var(ske) rtol=(high ? 0.02 : 0.05)
        if high
            @test skewness(sk) ≈ skewness(ske) rtol=0.1
        end
    end
end

# -----------------------------------------------------------------------------
# Tests for the s1s2/s1s2! fused S1/S2 kernels
#
# s1s2 computes the s1 and s2 inputs of skhat with fused reductions that never
# materialize a full-size intermediate array.  The reference below is the
# original three-reduction implementation, against which bit-identical (==)
# agreement is asserted: s1s2 uses the same underlying reduction kernels
# (sum! and Base.mapreducedim!) as sum(...; dims), so results must match exactly.
# -----------------------------------------------------------------------------

function _s1s2_reference(A::AbstractArray, M, N; dims::Integer=ndims(A))
    T = size(A, dims) ÷ (M*N)
    T > 0 || error("$(size(A, dims)) is too few samples for M=$(M) and N=$(N)")
    axesA = axes(A)
    Ausable = M*N*T == size(A, dims) ? A :
        view(A, axesA[1:dims-1]..., 1:M*N*T, axesA[dims+1:end]...)
    Anmt = reshape(Ausable, axesA[1:dims-1]..., N, M, T, axesA[dims+1:end]...)
    s0 = dropdims(sum(Anmt; dims); dims)
    s1 = dropdims(sum(s0; dims); dims)
    s2 = dropdims(sum(abs2, s0; dims); dims)
    (s1, s2)
end

@testset "s1s2" begin
    @testset "matches reference" begin
        @testset "size=$sz dims=$dims (M,N)=($M,$N)" for (sz, dims, M, N) in (
            ((16, 200),   2, 2, 1),
            ((16, 200),   2, 8, 1),
            ((16, 200),   2, 4, 3),
            ((16, 200),   2, 5, 4),
            ((16, 199),   2, 4, 3),   # M*N does not divide size(A, dims): truncation
            ((4, 5, 120), 3, 4, 3),   # 3-D batch, reducing dims=3
            ((4, 5, 120), 3, 5, 4),
            ((16, 50),    1, 4, 3),   # reducing dims=1
        )
            A = rand(MersenneTwister(hash((sz, dims, M, N))), sz...)
            (; s1, s2) = s1s2(A, M, N; dims)
            s1r, s2r = _s1s2_reference(A, M, N; dims)
            @test s1 == s1r
            @test s2 == s2r
            outsize = (sz[1:dims-1]..., sz[dims] ÷ (M*N), sz[dims+1:end]...)
            @test size(s1) == outsize && size(s2) == outsize

            # SKEstimator convenience method agrees
            ske = SKEstimator(M, N, 1)
            s1s = s1s2(A, ske; dims)
            @test s1s.s1 == s1r && s1s.s2 == s2r
        end
    end

    @testset "collapsed call" begin
        # Omitting M and N collapses the entire dims dimension into a single
        # accumulation per slice: M = size(A, dims), N = 1
        A = rand(MersenneTwister(20260928), 16, 200)
        r = s1s2(A)
        s1r, s2r = s1s2(A, 200, 1)
        @test r.s1 == s1r && r.s2 == s2r
        @test size(r.s1) == (16, 1) && size(r.s2) == (16, 1)

        # non-default dims: M follows dims
        A3 = rand(MersenneTwister(20260928), 4, 5, 120)
        r3 = s1s2(A3; dims=3)
        s1r3, s2r3 = s1s2(A3, 120, 1; dims=3)
        @test r3.s1 == s1r3 && r3.s2 == s2r3
        @test size(r3.s1) == (4, 5, 1) && size(r3.s2) == (4, 5, 1)

        # in-place collapsed call: garbage buffers fully overwritten
        b1 = fill!(Array{Float64}(undef, 16, 1),  1e300)
        b2 = fill!(Array{Float64}(undef, 16, 1), -1e300)
        result = s1s2!(b1, b2, A)
        @test result.s1 === b1 && result.s2 === b2
        @test b1 == s1r && b2 == s2r
        s1s2!(b1, b2, A)
        @test b1 == s1r && b2 == s2r

        # one SK value per channel via the explicit estimator
        ske = SKEstimator(200, 1, 1)
        sk = skhat.(s1r, s2r, ske)
        @test size(sk) == (16, 1)
        @test all(isfinite, sk)

        # degenerate: a single sample along dims gives M = 1, T = 1
        A1 = rand(MersenneTwister(20260928), 4, 1)
        r1 = s1s2(A1)
        @test r1.s1 == A1
        @test size(r1.s1) == (4, 1)
    end

    @testset "s1s2! overwrites garbage buffers" begin
        A = rand(MersenneTwister(20260928), 16, 200)
        M, N = 4, 3
        s1r, s2r = _s1s2_reference(A, M, N)
        buf1 = fill!(Array{Float64}(undef, 16, 16),  1e300)
        buf2 = fill!(Array{Float64}(undef, 16, 16), -1e300)
        result = s1s2!(buf1, buf2, A, M, N)
        @test result.s1 === buf1 && result.s2 === buf2
        @test buf1 == s1r && buf2 == s2r

        # idempotent: repeated calls into the same buffers give identical results
        s1s2!(buf1, buf2, A, M, N)
        @test buf1 == s1r && buf2 == s2r

        # SKEstimator convenience method agrees
        buf3 = similar(buf1); buf4 = similar(buf2)
        s1s2!(buf3, buf4, A, SKEstimator(M, N, 1))
        @test buf3 == s1r && buf4 == s2r
    end

    @testset "errors" begin
        A = rand(MersenneTwister(20260928), 16, 200)
        good = Array{Float64}(undef, 16, 16)
        @test_throws DimensionMismatch s1s2!(zeros(15, 16), zeros(16, 16), A, 4, 3)
        @test_throws DimensionMismatch s1s2!(zeros(16, 16), zeros(15, 16), A, 4, 3)
        @test_throws ArgumentError s1s2!(good, good, A,  0, 1)
        @test_throws ArgumentError s1s2!(good, good, A, -2, 1)
        @test_throws ArgumentError s1s2!(good, good, A,  4, 0)
        @test_throws ArgumentError s1s2!(good, good, A,  4, -1)
        @test_throws ErrorException s1s2(A, 201, 1)         # M*N > size(A, dims)
        @test_throws ErrorException s1s2(A, 60, 4)          # M*N > size(A, dims)
        @test_throws ErrorException s1s2!(good, good, A, 60, 4)
    end

    @testset "skhat(A) consistency" begin
        A = rand(MersenneTwister(20260928), 16, 200)
        ske = SKEstimator(4, 3, 1)
        (; s1, sk) = skhat(A, ske)
        s1r, s2r = _s1s2_reference(A, 4, 3)
        @test s1 == s1r
        @test sk == skhat(s1r, s2r, ske)

        # same result via the documented allocation-free pipeline
        buf1 = similar(s1); buf2 = similar(s1); skbuf = similar(sk)
        s1s2!(buf1, buf2, A, ske)
        skbuf .= skhat.(buf1, buf2, ske)
        @test buf1 == s1 && skbuf == sk
    end
end
