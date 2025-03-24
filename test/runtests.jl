using Test
using SpectralKurtosisEstimators
using SpectralKurtosisEstimators: relative_error, thresholds

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
