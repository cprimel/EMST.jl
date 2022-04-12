
using DataStructures
using Distances

import Base.isequal
import Base.hash
import Statistics.median


"""
  DTBSplit

implements split of KD-Tree
"""
struct DTBSplit
    dim::Int64
    val::Float64
    sgn::Bool # 0: smaller, 1 : greater
end

"""
  MyE

implements undirected edge
"""
struct MyE
    a::Int64
    b::Int64
    w::Float64
    MyE(pa::Int64,pb::Int64,w::Float64) = (pa<pb) ? new(pa,pb,w::Float64) : new(pb,pa,w::Float64)
end
function isequal(a::MyE,b::MyE)
    return a.a==b.a && a.b==b.b
end
function hash(a::MyE)
    return xor( hash(a.a) , hash(a.b) )
end
function to_matrices(ee::Vector{MyE})
    m::Array{Int64,2} = Array{Int64}(undef,length(ee),2)
    w::Array{Float64} = Array{Float64}(undef,length(ee))
    for zi=1:length(ee)
         m[zi,1] = ee[zi].a
         m[zi,2] = ee[zi].b
         w[zi] = ee[zi].w
    end
    return (m, w)
end


struct DTBProfiling
    n_comparisons::Int64
    n_break_boxdist::Int64
    n_break_onecomp::Int64
end
function DTBProfiling()
    return DTBProfiling(0,0,0)
end


"""
  DTBNode

implements node of a KD tree
"""
mutable struct DTBNode
    id::Int64 # identifies the node ; root has 00..00 , then
              # levels (right/left) are encoded as 01 or 10
              # i.e. node root->"right"->"right"->"left" ends with:
              # 10 01 01 00
              #  l r  r  R

    # data / indidces of this kd tree node
    data::Array{Float64,2}
    subset::Array{Int64}
    #splits::Array{DTBSplit,1}

    # kd tree node bounds
    box_lb::Array{Float64,1}
    box_ub::Array{Float64,1}

    # max_{v in this node} ( length of candidate edge of v)
    dQ::Float64

    left::DTBNode
    right::DTBNode

    #DTBNode(data::Array{Float64,2}) = (x = new(data); x.left=x; x.right=x)
    #function DTBNode(data::Array{Float64,2},subset::Array{Int64},splits::Array{DTBSplit,1},dQ::Float64)
    function DTBNode(id::Int64,data::Array{Float64,2},subset::Array{Int64},box_lb::Array{Float64,1},box_ub::Array{Float64,1},dQ::Float64)
        n = new(id,data,subset,box_lb,box_ub,dQ)
        n.left=n
        n.right=n
        return n
    end
end
function isequal(a::DTBNode,b::DTBNode)
    return a.id==b.id
end
function hash(a::DTBNode)
    #a.id
    hash(a.id) # we could also drop this hash.. should not matter..
end
function is_leaf(n::DTBNode)
    return n==n
    #return ( n.id & (1<<63) ) == (1<<63)
end

"""
  kdtree(xx::Array{Float64,2})

initializes root node of KD tree (but does not compute the splits)
To build the tree use:
  kdt_root = kdtree( data )
  kdtree_split!( kdt_root , 10 ) # 10 is the node size where we stop splitting
"""
function kdtree(xx::Array{Float64,2})
    #root = DTBNode(xx,-1,NaN,nothing,nothing)
    #root = DTBNode(Int64(0),xx,collect(1:size(xx,2)),fill!(ones(size(xx,1)),-Inf),fill!(ones(size(xx,1)),Inf),Inf)
    root = DTBNode(Int64(0),xx,collect(Int64(1):Int64(size(xx,2))),fill!(ones(size(xx,1)),-Inf),fill!(ones(size(xx,1)),Inf),Inf)
    return root
end


"""
  kdtree_split!(node::DTBNode, nmin::Int64)

computes splits of KD tree. nmin indicates the node size where we stop
splitting nodes.
"""
function kdtree_split!( node::DTBNode, nmin::Int64)

    # ensure that parent node is in the tree.
    #t[node.id] = node

    if( size(node.data,2) <= nmin )
        return node
    end

    mind = minimum(node.data,dims=2)
    maxd = maximum(node.data, dims=2)
    s = maxd - mind
    # println(s)
    ds   = findmax(s)
    # print("$ds")
    ds   = ds[2][1]
    #println("$ds")
    #println(node.data[ds,:])
    vs   = median(node.data[ds,:])
    # println("$vs")
    bx   = node.data[ds,:].<= vs
    # println(bx)
    range_a = node.subset[bx]
    range_b = node.subset[.~bx]

    data_a = node.data[:,bx]
    data_b = node.data[:,.~bx]

    box_lb_a = copy(node.box_lb)
    box_ub_a = copy(node.box_ub)
    box_ub_a[ds] = vs
    box_lb_b = copy(node.box_lb)
    box_lb_b[ds] = vs
    box_ub_b     = copy(node.box_ub)

    id::Int64 = node.id
    id_depth  = Int( ceil( (64-leading_zeros(id)) / 2) ) + 1 # in this cell we are, i.e. we have to shift id_depth times left by two bits..
    id_l = id   |  (1)<<(2*id_depth)
    #id_l = id_l |  (1)<<(63)
    id_r = id   |  (2)<<(2*id_depth)
    #id_r = id_r |  (1)<<(63)
    #println("id $(id) -> r $(id_l) , l $(id_r)")

    node_left  = kdtree_split!( DTBNode(id_l,data_a,range_a, box_lb_a   , box_ub_a  , Inf) , nmin )
    node_right = kdtree_split!( DTBNode(id_r,data_b,range_b, box_lb_b   , box_ub_b  , Inf) , nmin )

    node.left  = node_left
    node.right = node_right
    #println("SPLIT") #println(data_a) #println(data_b)
    #println(splits_l)
    return node
end

"""
  compute_emst(data::Array{Float64,2};nmin::Int64=64)

Computes EMST for the given data (where columns are samples).
nmin is the max number of elements in kd-tree node.
"""
function compute_emst(data::Array{Float64,2};nmin::Int=64)
    nmin64 = Int64(nmin)
    root = kdtree(data)
    kdtree_split!(root,nmin64)
    edges = dtb(root,IntDisjointSets(size(data,2)))
    return EMST.to_matrices( collect(edges) )
end


"""
  dtb(q::DTBNode,e::IntDisjointSets)

implements the dual-tree Boruvka algorithm which computes the EMST for the
given KD-tree.
"""
function dtb(q::DTBNode,e::IntDisjointSets)
    edges = Set{MyE}()

    #hashmap_QvsR = Dict{Tuple{DTBNode,DTBNode},Array{Float64,2}}()

    while(e.ngroups>1)
        ngroups = e.ngroups;
        println("--> ngroups: $ngroups")

        # prepare dicts for candidate edges..
        C_dcq = Dict{Int64,Float64}()
        C_e   = Dict{Int64,MyE}()
        # init themz
        roots = unique(e.parents)
        for ri in roots
            C_dcq[ri] = Inf
        end

        # compute..
        #fcn(q,q,e,C_dcq,C_e , hashmap_QvsR )
        fcn(q,q,e,C_dcq,C_e)

        # and now add the edges..
        for ne::MyE in values(C_e)
            union!(e,ne.a,ne.b)
            push!(edges,ne)
        end
    end
    return edges
end


"""
  fcn(q::DTBNode,r::DTBNode,e::IntDisjointSets , C_dcq::Dict{Int64,Float64} , C_e::Dict{Int64,MyE} )
Find Component Neighbors

Components are identified by the root of the component
C_dcq : component distances to candidate edges
C_e   : component candidate edges (comp i candidate edge is C_e[:,i])
"""
function fcn(q::DTBNode,r::DTBNode,e::IntDisjointSets , C_dcq::Dict{Int64,Float64} , C_e::Dict{Int64,MyE} )
#function fcn(q::DTBNode,r::DTBNode,e::IntDisjointSets , C_dcq::Dict{Int64,Float64} , C_e::Dict{Int64,Tuple{Int64,Int64}} )
#function fcn(q::DTBNode,r::DTBNode,e::IntDisjointSets , C_dcq::Dict{Int64,Float64} , C_e::Dict{Int64,MyE} , hashmap_QvsR::Dict{Tuple{DTBNode,DTBNode},Array{Float64,2}})

    # check all in same component
    onecomp::Bool  = true
    joined         = [q.subset;r.subset]
    #println(joined)
    for ji in joined
        if(~in_same_set(e,joined[1],ji))
            onecomp=false
            break
        end
    end
    if(onecomp)
        return
    end

    #check d(Q,R) > d(Q)
    dqr = computeDQR( q.box_lb,q.box_ub , r.box_lb,r.box_ub )
    #println("dq : $(dqr) , q.dQ : $(q.dQ)")
    if( dqr > q.dQ)
        #println(".")
        return
    end

    #check if R and Q in a leaf node
    if( q.left==q && r.left==r ) # that's not very elegant but probably what we have to do?!?
    #if( is_leaf(q) && is_leaf(r) )
        #println("LEAF NODES")
        #n_dQ = -Inf
        n_dQ::Float64 = q.dQ

        all_d_qr = Distances.pairwise(Euclidean(),q.data,r.data,dims=2)

        # --> hashing these values does not really speed up things..
        #     and it uses tons of memory, such that fancier handling of
        #     hash map entries would be required
        #if(!haskey(hashmap_QvsR,(q,r)))
        #    if(!haskey(hashmap_QvsR,(r,q)))
        #        all_d_qr = Distances.pairwise(Euclidean(),q.data,r.data)
        #        #hashmap_QvsR[(q,r)] = all_d_qr
        #    else
        #        all_d_qr = transpose(hashmap_QvsR[(r,q)])
        #    end
        #    hashmap_QvsR[(r,q)] = transpose(all_d_qr)
        #else
        #    all_d_qr = hashmap_QvsR[(q,r)]
        #    #print('.')
        #end

        for iq = 1:size(q.subset,1)
            for ir = 1:size(r.subset,1)
                qq = q.subset[iq]
                rr = r.subset[ir]
                if( in_same_set(e , qq , rr ) )
                    continue
                end

                cq = find_root!(e,rr) # compoment of q

                # check distance:
                # dist_qr::Float64 = sqrt( sum( (q.data[:,iq]-r.data[:,ir]).^2 ) )
                dist_qr = all_d_qr[iq,ir]

                if( dist_qr < C_dcq[ cq ] )
                    C_dcq[ cq ] = dist_qr
                    C_e[ cq ]   = MyE(qq,rr, dist_qr) #(qq,rr)
                    # and dQ !
                    n_dQ = max(n_dQ,dist_qr)
                    #n_dQ = dist_qr
                end
            end
        end
        q.dQ = n_dQ

        return
    end

    # do recursions..
    #fcn( q.left,r.left,e ,  C_dcq::Dict{Int64,Float64} , C_e::Dict{Int64,Tuple{Int64,Int64}} )
    #q_left  = is_leaf(q) ? q : t[q.left]
    #q_right = is_leaf(q) ? q : t[q.right]
    #r_left  = is_leaf(r) ? r : t[r.left]
    #r_right = is_leaf(r) ? r : t[r.right]
    fcn( q.left,r.left    , e ,  C_dcq , C_e ) #, hashmap_QvsR )
    fcn( q.right,r.left   , e ,  C_dcq , C_e ) #, hashmap_QvsR )
    fcn( q.left,r.right   , e ,  C_dcq , C_e ) #, hashmap_QvsR )
    fcn( q.right,r.right  , e ,  C_dcq , C_e ) #, hashmap_QvsR )
    q.dQ = max( q.left.dQ , q.right.dQ )
end


"""
  computeDQR( q_lb::Array{Float64,1} , q_ub::Array{Float64,1} , r_lb::Array{Float64,1} , r_ub::Array{Float64,1} )

compute min dist. between bounding boxes, i.e. between rectangular boxes Q/R with
bounds given by q_lb/q_ub and r_lb/r_ub
"""
function computeDQR( q_lb::Array{Float64,1} , q_ub::Array{Float64,1} , r_lb::Array{Float64,1} , r_ub::Array{Float64,1} )
    d::Int64  = length(q_lb)#size(q_lb,1)
    rdists::Array{Float64}    = zeros(d)
    for zd = 1:d
        rdists[zd] = max(  max( q_lb[zd] - r_ub[zd] , r_lb[zd] - q_ub[zd] ) , 0 )
    end
    dqr::Float64 = sqrt( sum( rdists.^2 ))
    return dqr
end
