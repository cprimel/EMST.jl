
using Random

"""
  verify_emst(x::Array{Float64,2},e::Array{Int64,2},n_checks::Int64)

Checks whether the computed EMST contains the nearest neighbor edges.
"""
function verify_emst(x::Array{Float64,2}, e::Array{Int64,2}, n_checks::Int)
    n = min(n_checks, size(x, 2))
    idx_rand = randperm(size(x, 2))[1:n]
    return verify_emst(x, e, idx_rand)
end

"""
  verify_emst(x::Array{Float64,2},e::Array{Int64,2},idx_check::Vector{Int64})

Checks whether the computed EMST contains the nearest neighbor edges for the
given indeces.
"""
function verify_emst(x::Array{Float64,2}, e::Array{Int64,2}, idx_check::Vector{Int})
    # create edge set
    es::Set{Set{Int64}} = Set{Set{Int64}}()
    for zi = 1:size(e, 1)
        push!(es, Set([e[zi, 1], e[zi, 2]]))
    end

    print("\nVerify: ")
    ok = true
    for zi in 1:length(idx_check)
        ix = idx_check[zi]
        xd = sum((broadcast(-, x[:, ix], x)) .^ 2, dims=1)
        xd[ix] = Inf
        i_min = findmin(xd)[2]
        if (sum(xd .== xd[i_min]) > 1)
            println("WARNING: Equidistant vertices discovered..")
            ok = false
            iequi = find(any(xd .== xd[i_min]))
            for iie in iequi
                ok = isempty(setdiff(Set([Set([ix, iie])]), es))
                if (ok)
                    break
                end
            end
            if (!ok)
                println("emst kaputt! $(ix) : edge to one of $(iequi) should be in emst!")
                break
            end
        else
            i_min = i_min[1]
            ok = isempty(setdiff(Set([Set([ix, i_min])]), es))
            if (!ok)
                println("emst kaputt! $(ix) : edge to $(i_min) should be in emst!")
                break
            end
        end
        if (mod(zi, ceil(length(idx_check) / 10.0)) == 0)
            print(".")
        end
    end
    if (ok)
        print(" ok!\n")
    end
    return ok
end
