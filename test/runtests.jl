
##
#reload("EMST.jl")
using EMST
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

function test_emst_with_uniform(d,n;check_verification=false)
    x = rand(d,n)
    # x = rand(2,2000)

    (x_emst, _)    = EMST.compute_emst(x)
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





@testset "EMST Tests A" begin
    #for zi=1:100
     #   test_emst_with_uniform(2,500)
    #end
    #for zi=1:100
     #   test_emst_with_uniform(17,500)
    #end
    # and check that our verification actually works.. :)
    for zi=1:20
        test_emst_with_uniform(4,1000;check_verification=true)
    end
end


@testset "EMST Tests B" begin
    test_emst_with_uniform(2,1000)
    test_emst_with_uniform(123,1000)
    test_emst_with_uniform(10,4000)
    test_emst_with_uniform(10,8000)
    test_emst_with_uniform(2,16000)
    test_emst_with_uniform(44,16000)
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
