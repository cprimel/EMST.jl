using EMST
using Test
using MLDatasets

function test_iris()
    #Iris.download(; i_accept_the_terms_of_use=true);
    iris_features = Iris.features();
    x = unique(iris_features, dims=2)
    x = Array{Float64,2}(x)
    (x_emst, _, _) = EMST.compute_emst(x)

    good_emst = EMST.verify_emst(x, x_emst,size(x,2))
    @test good_emst
end

function test_emst_with_uniform(d,n,l;check_verification=false)
    x = rand(d,n)
    # x = rand(2,2000)

    (x_emst, _, _)    = EMST.compute_emst(x;leafSize=l)
    good_emst = EMST.verify_emst(x,x_emst,size(x,2))
    @test good_emst

    if(check_verification)
        # and show that verify works..
        x_emst_kaputt = copy(x_emst)
        # change n edges..
        for zi in 1:min(n-1,5) ; x_emst_kaputt[ rand( 1:( size(x,2))-1) , rand(1:2) ] = rand(1:size(x,2)); end
        emst_verification_working =  ~verify_emst(x,x_emst_kaputt,size(x,2))
        @test emst_verification_working
    end
end


function test_kdtree()
    iris_features = Iris.features();
    x = unique(iris_features, dims=2)
    x = Array{Float64,2}(x)

    #check_equality(x)

    root = kdtree(x)
    kdtree_split!(root, 1)
    oldfromnew = Vector{Int64}()
    getleaves(root, oldfromnew)

end

@testset "Tests" begin
    #@testset "Test with verification" begin
    #    for zi=1:2
    #        test_emst_with_uniform(8,,1;check_verification=true)
    #    end
    #end
    #@testset "With uniform" begin
    #    for zi=1:2
    #        test_emst_with_uniform(8,1000,1;check_verification=false)
    #    end
    #end
    @testset "Test kdtree" begin
        test_kdtree()
    end
    @testset "With Iris" begin
        test_iris()
    end
end


#using Plots
#using Gadfly
#function plot_emst_2d(x,edges)
#    Plots.gr()
#    Plots.plot()
#    Plots.scatter!(x[1,:],x[2,:],ms=2.0,marker=(stroke(0,:gray)))
#    for zr in 1:size(edges,1)
#        Plots.plot!( x[1,edges[zr,:]] , x[2,edges[zr,:]] )
#    end
#    Plots.scatter!() # for some reason with this the plot shows im juno.. :)
#end
#x_test  = rand(2,2000)
#tic(); e_test  = EMST.compute_emst(x_test;nmin=64); toc();
#EMST.verify_emst(x_test,e_test,200)
#plot_emst_2d(x_test,e_test)


##
