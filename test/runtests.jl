using SpectralKurtosisEstimators
using Test

# Unit tests reconstruct some of the entries of Table 1 from this paper:
#
# Nita, G. M. & Gary, D. E. [2010b] MNRAS 406, L60,
# doi:10.1111/j.1745-3933.2010.00882.x.

M = 300
N = 10
d = 1
nsigma = 3

lowervi_expected = 0.766_48
uppervi_expected = 1.283_13
# Table 1 has 0.18% of u4errvi_expected, but our errors are not percent
u4errvi_expected = 0.001_8

loweriii_expected = 0.767_54
upperiii_expected = 1.282_12
# Table 1 has -0.71% of u4erriii_expected, but our errors are absolute value
u4erriii_expected = 0.007_1

ske = SKEstimator(M, N, d)

skdvi, u4errvi = @test_logs (:warn,"pearson critereon < 1") pearson_type_vi(ske)
skdiii, u4erriii = @test_logs (:warn,"pearson critereon < 1") pearson_type_iii(ske)

lowervi, uppervi = SpectralKurtosisEstimators.thresholds(skdvi, nsigma)
loweriii, upperiii = SpectralKurtosisEstimators.thresholds(skdiii, nsigma)

@test lowervi ≈ lowervi_expected atol=0.000_005
@test uppervi ≈ uppervi_expected atol=0.000_005
@test u4errvi ≈ u4errvi_expected atol=0.000_05

@test loweriii ≈ loweriii_expected atol=0.000_005
@test upperiii ≈ upperiii_expected atol=0.000_005
@test u4erriii ≈ u4erriii_expected atol=0.000_05
