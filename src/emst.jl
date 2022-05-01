
using DataStructures

import Base.isequal
import Base.hash

export compute_emst


"""
    compute_emst(data[,leafSize])

Computes EMST for the given data (where columns are samples).
leafSize is the max number of elements in kd-tree node.
"""
function compute_emst(data::Array{Float64,2}; leafSize::Int64=1)
    root = kdtree(data, leafSize)
    """
        TODO: compute mapping during construction
    """
    oldfromnew = Vector{Int64}()
    getleaves(root, oldfromnew)
    
    edges, weights = dtb(root, IntDisjointSets(size(data, 2)))
"""
    # NOTE: Might not make sense to do this processing here
    edges::Array{Int64,2} = Array{Int64}(undef,size(e_out)) 
    for i = 1:size(e_out, 1)
        indexA::Int64 = oldfromnew[e_out[i, :][1]]
        indexB::Int64 = oldfromnew[e_out[i, :][2]]

        if indexA < indexB
            edges[i,1] = indexA
            edges[i,2] = indexB
        else
            edges[i,1] = indexB
            edges[i,2] = indexA
        end
    end
"""
    return edges, weights, oldfromnew
end



