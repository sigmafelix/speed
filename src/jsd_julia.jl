## Jensen-Shannon divergence
using Base.MathConstants


function scale_minmax(vect)
    vect = (vect .- minimum(vect)) ./ (maximum(vect) - minimum(vect))
    return vect
end

function scale_minmax_matrix(mat, colrange)
    for i in colrange
        mat[:,i] = scale_minmax(mat[:,i])
    end
    return mat
end

function check_matrix(mat)
    res = sum(mat .< 0)
    return res
end

function jensen_shannon(p, q)
    n = length(p)
    m = zeros(Float64, n)
    sum = 0.0
    for i in 1:n
        m[i] = 0.5 * (p[i] + q[i])
        sum += 0.5 * (p[i] * log2(2.0 * p[i] / (p[i] + q[i])) + q[i] * log2(2.0 * q[i] / (p[i] + q[i])))
    end
    return sum
end

function jensen_shannon_divergence(mat_tr, mat_co = nothing, scaling = true)
    if mat_co == nothing
        mat_co = mat_tr
    end

    nrow_tr, ncol_tr = size(mat_tr)
    nrow_co, ncol_co = size(mat_co)
    res = zeros(Float64, nrow_tr, nrow_co)

    for i in 1:nrow_tr
        for j in 1:nrow_co
            p = mat_tr[i,:]
            q = mat_co[j,:]
            res[i, j] = jensen_shannon(p, q)
        end
    end
    return res
end

