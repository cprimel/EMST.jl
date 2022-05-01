import Statistics.median


export kdtree, KDNode

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
function getleaves(root::KDNode, leafsubsets::Vector{Int64})
    if root == root.left && root == root.right
        append!(leafsubsets,root.subset)
        return
    end
    getleaves(root.left, leafsubsets)
    getleaves(root.right, leafsubsets)
    return
end


"""
Initializes and builds k-d tree. Returns root of tree.
"""
function kdtree(xx::Array{Float64,2}, leafSize::Int64; splitType::String="median")
    box_lb = vec(minimum(xx,dims=2))
    box_ub = vec(maximum(xx, dims=2))
    
    root = KDNode(Int64(0), xx, collect(Int64(1):Int64(size(xx, 2))), box_lb, box_ub, Inf)
    kdtree_split!(root, leafSize, splitType)
    return root
end

"""
Computes splits of KD tree. Stop splitting when node contains number of samples equal to leafSize.
"""
function kdtree_split!(node::KDNode, leafSize::Int64, splitType::String)

    if (size(node.data, 2) <= leafSize)
        return node
    end
    
    # Split by median
    if splitType == "median" 
        mind = minimum(node.data, dims=2)
        maxd = maximum(node.data, dims=2)
        s = maxd - mind
        ds = findmax(s)
        ds = ds[2][1]
        vs = median(node.data[ds, :])
    elseif splitType == "midpoint"
    # Split by midpoint
        s = node.box_ub - node.box_lb
        ds = findmax(s)[2][1]
        vs = node.box_ub[ds] - s[ds] / 2.0
    end

    bx = node.data[ds, :] .<= vs
    
    # If all points are less than or equal to split, 
    if all(bx)
        bx = node.data[ds, :] .< vs
        if all(.~bx)
            return node
        end
    end
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
    id_r = id | (2) << (2 * id_depth)


    node_left = kdtree_split!(KDNode(id_l, data_a, range_a, box_lb_a, box_ub_a, Inf) ,leafSize, splitType)
    node_right = kdtree_split!(KDNode(id_r, data_b, range_b, box_lb_b, box_ub_b, Inf), leafSize, splitType)

    node.left = node_left
    node.right = node_right
    return node
end
