using DataStructures
using Distances

export dtb

"""
    Edge(a, b[, w])

Implements undirected edge with optional weight. Paired vertices are ordered by id number.
...
# Fields:
    - `a`: the id of start vertex 
    - `b`: the id of end vertex
    - `w=nothing`: optional weight parameter. 
...
"""
struct Edge
    a::Int64
    b::Int64
    w::Float64
    Edge(pa::Int64, pb::Int64) = new(pa::Int64, pb::Int64, nothing)
    Edge(pa::Int64, pb::Int64, w::Float64) = (pa < pb) ? new(pa, pb, w::Float64) : new(pb, pa, w::Float64)
end
function isequal(a::Edge, b::Edge)
    return a.a == b.a && a.b == b.b
end

"""
    write_edgelist(edges)

Creates edge list and weight list from vector of Edge structs, ordered by weight. 

Note: Seperate lists needed since different types.
"""
function write_edgelist(ee::Vector{Edge})
    ee = sort!(ee, by = e -> e.w)
    m::Array{Int64,2} = Array{Int64}(undef, length(ee), 2)
    w::Array{Float64} = Array{Float64}(undef, length(ee))
    for zi = 1:length(ee)
        m[zi, 1] = ee[zi].a
        m[zi, 2] = ee[zi].b
        w[zi] = ee[zi].w
    end
    return (m, w)
end



"""
Wraps main dual_tree_boruvka function. 
"""
function dtb(q::KDNode, e::IntDisjointSets)
    edges = dual_tree_boruvka(q, e)
    edges, weights = write_edgelist(collect(edges))
    return edges, weights
end

"""
    dual_tree_boruvka(q::KDNode,e::IntDisjointSets)

Implements the dual-tree Boruvka algorithm which computes the Euclidean minimum spanning tree for the given KD-tree.
"""
function dual_tree_boruvka(q::KDNode, e::IntDisjointSets)
    edges = Set{Edge}()

    while (e.ngroups > 1)
        ngroups = e.ngroups
        println("--> ngroups: $ngroups")

        # Initialize dictionary for candidate edges
        C_dcq = Dict{Int64,Float64}()
        C_e = Dict{Int64,Edge}()
        roots = unique(e.parents)

        for ri in roots
            C_dcq[ri] = Inf
        end

        # Compute component neighbors
        find_cn(q, q, e, C_dcq, C_e)
        # and now add the edges..
        for ne::Edge in values(C_e)
            union!(e, ne.a, ne.b)
            push!(edges, ne)
        end
    end
    return edges
end


"""
    find_cn(
        q::KDNode,
        r::KDNode,
        e::IntDisjointSets, 
        C_dcq::Dict{Int64,Float64},
        C_e::Dict{Int64,Edge}
    )

Find Component Neighbors

...
# Arguments:
- `q`: root of first component
- `r`: root of second component
- `e`: union-find datastructure for tracking edges
- `C_dcq`: Component distances to candate edges
- `C_e`: Component candidate edges (the candidate edge of component `i` is C_e[:,i])
"""
function find_cn(q::KDNode, r::KDNode, e::IntDisjointSets, C_dcq::Dict{Int64,Float64}, C_e::Dict{Int64,Edge})

    # Check that all are in same component
    onecomp::Bool = true
    joined = [q.subset; r.subset]
    for ji in joined
        if (~in_same_set(e, joined[1], ji))
            onecomp = false
            break
        end
    end
    if (onecomp)
        return
    end

    # Check that d(Q,R) > d(Q)
    dqr = computeDQR(q.box_lb, q.box_ub, r.box_lb, r.box_ub)
    if (dqr > q.dQ) 
        return
    end

    # Check if R and Q in a leaf node
    if (q.left == q && r.left == r)

        n_dQ::Float64 = q.dQ

        all_d_qr = Distances.pairwise(Euclidean(), q.data, r.data, dims=2)

        for iq = 1:size(q.subset, 1)
            for ir = 1:size(r.subset, 1)
                qq = q.subset[iq]
                rr = r.subset[ir]
                if (in_same_set(e, qq, rr))
                    continue
                end

                cq = find_root!(e, rr) # component of q

                # check distance:
                dist_qr = all_d_qr[iq, ir]
                if (dist_qr < C_dcq[cq])
                    C_dcq[cq] = dist_qr
                    C_e[cq] = Edge(qq, rr, dist_qr) #(qq,rr)
                    # and dQ !
                    n_dQ = max(n_dQ, dist_qr)
                    #println(n_dQ)
                    #n_dQ = dist_qr
                    #println(n_dQ)
                end
            end
        end
        q.dQ = n_dQ # always Inf

        return
    end

    # do recursions..
    find_cn(q.left, r.left, e, C_dcq, C_e)
    find_cn(q.right, r.left, e, C_dcq, C_e)
    find_cn(q.left, r.right, e, C_dcq, C_e)
    find_cn(q.right, r.right, e, C_dcq, C_e)

    q.dQ = max(q.left.dQ, q.right.dQ)
end


"""
  computeDQR( q_lb::Array{Float64,1} , q_ub::Array{Float64,1} , r_lb::Array{Float64,1} , r_ub::Array{Float64,1} )

compute min dist. between bounding boxes, i.e. between rectangular boxes Q/R with
bounds given by q_lb/q_ub and r_lb/r_ub

NOTE: midpoint split leads to deep / infinite recursions related to this function.
"""
function computeDQR(q_lb::Array{Float64,1}, q_ub::Array{Float64,1}, r_lb::Array{Float64,1}, r_ub::Array{Float64,1})
    d::Int64 = length(q_lb)
    rdists::Array{Float64} = zeros(d)
    for zd = 1:d
        rdists[zd] = max(max(q_lb[zd] - r_ub[zd], r_lb[zd] - q_ub[zd]), 0)
    end
    dqr::Float64 = sqrt(sum(rdists .^ 2))
    return dqr
end
