function pearson_critereon(B1, B2)
    (
        (B1 * (B2+3)^2)
        /
        (4(4B2-3B1) * (2B2-3B1-6))
    )
end

function pearson_critereon(ske::SKEstimator)
    B1 = ske.u3^2 / ske.u2^3
    B2 = ske.u4 / ske.u2^2

    pearson_critereon(B1, B2)
end

function pearson_critereon(M, N, d)
    pearson_critereon(SKEstimator(M, N, d))
end

function pearson_critereon(M)
    pearson_critereon(SKEstimator(M, 1, 1))
end

# Pearson type III

function pearson_type_iii(u2, u3)
    # The parameter naming convention of the 2010c SK paper differs from the
    # parameter naming convention of Distributions.jl for the Gamma
    # distribution.  The paper uses `β` for the *shape* parameter and `α` for
    # the *scale* parameter whereas `Distributions.Gamma` uses `α` for the shape
    # parameter and `θ` for the scale parameter.  To avoid (more) confusion, we
    # will refer to the shape and scale parameters here as `shape` and `scale`.
    # The *location* parameter called `δ` in the paper will be called `location`
    # here.
    shape = 4u2^3 / u3^2
    scale = u3 / 2u2
    location = 1 - 2u2^2 / u3

    Gamma(shape, scale) + location
end

function pearson_type_iii(ske::SKEstimator)
    pearson_critereon(ske) > 1 || @warn("pearson critereon not > 1")

    pearson_type_iii(ske.u2, ske.u3)
end

function pearson_type_iii(A::AbstractArray)
    pearson_type_iii(moment(A, 2), moment(A, 3))
end

function pearson_type_iii(M, N, d)
    pearson_type_iii(SKEstimator(M, N, d))
end

# Pearson type VI

function pearson_type_vi(u2, u3)
    radical = sqrt(16u2^4 + 4u3^2*u2 + u3^2)

    a = (
        32u2^5 - 4u3*u2^3 + 8u3^2*u2^2 + u3^2*u2 - u3^3
        + (8u2^3 - u3*u2 + u3^2) * radical
    ) / u3^3

    B = 2u2 * (4u2^2 + radical) / u3^2 + 3

    d = (B-a-1)/(B-1)

    BetaPrime(a, B) + d
end

function pearson_type_vi(ske::SKEstimator)
    pearson_critereon(ske) > 1 || @warn("pearson critereon not > 1")

    pearson_type_vi(ske.u2, ske.u3)
end

function pearson_type_vi(A::AbstractArray)
    pearson_type_vi(moment(A, 2), moment(A, 3))
end

function pearson_type_vi(M, N, d)
    pearson_type_vi(SKEstimator(M, N, d))
end
