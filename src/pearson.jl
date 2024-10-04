function pearson_criterion(B1, B2)
    (
        (B1 * (B2+3)^2)
        /
        (4(4B2-3B1) * (2B2-3B1-6))
    )
end

function pearson_criterion(ske::SKEstimator)
    B1 = ske.u3^2 / ske.u2^3
    B2 = ske.u4 / ske.u2^2

    pearson_criterion(B1, B2)
end

function pearson_criterion(M, N, d)
    pearson_criterion(SKEstimator(M, N, d))
end

function pearson_criterion(M)
    pearson_criterion(SKEstimator(M, 1, 1))
end

