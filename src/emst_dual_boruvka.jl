
using DataStructures
using Distances

import Base.isequal
import Base.hash
import Statistics.median


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

Creates (unweighted) edge list and weight list from vector of Edge structs. 

Note: Seperate lists needed since different types.
"""
function write_edgelist(ee::Vector{Edge})
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
struct DTBProfiling
    n_comparisons::Int64
    n_break_boxdist::Int64
    n_break_onecomp::Int64
end
function DTBProfiling()
    return DTBProfiling(0,0,0)
end
"""

"""
    KDNode

Implements node of a KDtree.

...
# Fields
- `id`:
- `data`: data in 2-dimensional array (row=sample, column=feature)
- `subset`: indices of samples contained in node
- `box_lb`: lower bound of node
- `box_ub`: upper bound of node
- `dQ`: length of candidate edge
- `left`: reference to left child node
- `right`: reference to right child node
...
"""
mutable struct KDNode
    id::Int64 # identifies the node ; root has 00..00 , then
    # levels (right/left) are encoded as 01 or 10
    # i.e. node root->"right"->"right"->"left" ends with:
    # 10 01 01 00
    #  l r  r  R

    data::Array{Float64,2}
    subset::Array{Int64}

    box_lb::Array{Float64,1}
    box_ub::Array{Float64,1}

    dQ::Float64

    left::KDNode
    right::KDNode

    function KDNode(id::Int64, data::Array{Float64,2}, subset::Array{Int64}, box_lb::Array{Float64,1}, box_ub::Array{Float64,1}, dQ::Float64)
        n = new(id, data, subset, box_lb, box_ub, dQ)
        n.left = n
        n.right = n
        return n
    end
end
function isequal(a::KDNode, b::KDNode)
    return a.id == b.id
end
function is_leaf(n::KDNode)
    return n == n
end

"""
    kdtree(data)

Initializes root node of KD tree

Note:: To build a KD tree:
            kdt_root = kdtree( data )
            kdtree_split!( kdt_root , 10 )
"""
function kdtree(xx::Array{Float64,2})
    root = KDNode(Int64(0), xx, collect(Int64(1):Int64(size(xx, 2))), fill!(ones(size(xx, 1)), -Inf), fill!(ones(size(xx, 1)), Inf), Inf)
    return root
end


"""
    kdtree_split!(node, leafSize)

Computes splits of KD tree. Stop splitting when node contains number of samples equal to leafSize.
"""
function kdtree_split!(node::KDNode, nmin::Int64)

    if (size(node.data, 2) <= nmin)
        return node
    end
    #println(size(node.data, 2))
    if (length(node.data) < 1)
        return node
    end
    mind = minimum(node.data, dims=2)
    maxd = maximum(node.data, dims=2)
    s = maxd - mind
    ds = findmax(s)
    ds = ds[2][1]
    vs = median(node.data[ds, :])
    bx = node.data[ds, :] .<= vs
    range_a = node.subset[bx]
    range_b = node.subset[.~bx]

    data_a = node.data[:, bx]
    data_b = node.data[:, .~bx]

    box_lb_a = copy(node.box_lb)
    box_ub_a = copy(node.box_ub)
    box_ub_a[ds] = vs
    box_lb_b = copy(node.box_lb)
    box_lb_b[ds] = vs
    box_ub_b = copy(node.box_ub)

    id::Int64 = node.id
    id_depth = Int(ceil((64 - leading_zeros(id)) / 2)) + 1 # in this cell we are, i.e. we have to shift id_depth times left by two bits..
    id_l = id | (1) << (2 * id_depth)
    #id_l = id_l |  (1)<<(63)
    id_r = id | (2) << (2 * id_depth)
    #id_r = id_r |  (1)<<(63)
    #println("id $(id) -> r $(id_l) , l $(id_r)")

    node_left = kdtree_split!(KDNode(id_l, data_a, range_a, box_lb_a, box_ub_a, Inf) ,nmin)
    node_right = kdtree_split!(KDNode(id_r, data_b, range_b, box_lb_b, box_ub_b, Inf), nmin)

    node.left = node_left
    node.right = node_right
    return node
end

"""
    compute_emst(data[,leafSize])

Computes EMST for the given data (where columns are samples).
leafSize is the max number of elements in kd-tree node.
"""
function compute_emst(data::Array{Float64,2}; leafSize::Int64=64)
    root = kdtree(data)
    kdtree_split!(root, leafSize)
    oldfromnew = Vector{Int64}()
    getleafs(root, oldfromnew)
    edges = dtb(root, IntDisjointSets(size(data, 2)))
    e_out, w_out = EMST.write_edgelist(collect(edges))
    return e_out, w_out, oldfromnew
end


function getleafs(root::KDNode, leafsubsets::Vector{Int64})
    if root == root.left && root == root.right
        append!(leafsubsets,root.subset)
        return
    end
    getleafs(root.left, leafsubsets)
    getleafs(root.right, leafsubsets)
    return
end

"""
    dtb(q::KDNode,e::IntDisjointSets)

Implements the dual-tree Boruvka algorithm which computes the Euclidean minimum spanning tree for the given KD-tree.
"""
function dtb(q::KDNode, e::IntDisjointSets)
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

                cq = find_root!(e, rr) # compoment of q

                # check distance:
                dist_qr = all_d_qr[iq, ir]
                if (dist_qr < C_dcq[cq])
                    C_dcq[cq] = dist_qr
                    C_e[cq] = Edge(qq, rr, dist_qr) #(qq,rr)
                    # and dQ !
                    n_dQ = max(n_dQ, dist_qr)
                end
            end
        end
        q.dQ = n_dQ

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
