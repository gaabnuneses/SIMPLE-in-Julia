using Plots

let
    plot(leg=false,aspect_ratio=:equal)
    X = [x[1]*cos(j) for j in y]
    Y = [x[1]*sin(j) for j in y]
    plot!(X,Y,X*0 .+z[1],color=:black,lw=1)
    plot!(X,Y,X*0 .+z[end],color=:black,lw=1)
    # plot!([X';X'],[Y';Y'],[(X'*0 .+z[1]);(X'*0 .+z[end])],color=:black,lw=1)

    X = [x[end]*cos(j) for j in y]
    Y = [x[end]*sin(j) for j in y]
    plot!(X,Y,X*0 .+z[1],color=:black,lw=1)
    plot!(X,Y,X*0 .+z[end],color=:black,lw=1)
    plot!([X';X'],[Y';Y'],[(X'*0 .+z[1]);(X'*0 .+z[end])],color=:black,lw=1)

    # Inner cylinder
    X = [.4*cos(j) for j in y]
    Y = [.4*sin(j) for j in y]
    plot!(X,Y,X*0 .+.2,color=:black,lw=1)
    plot!(X,Y,X*0 .+z[end],color=:black,lw=1)
    plot!([X';X'],[Y';Y'],[(X'*0 .+.2);(X'*0 .+z[end])],color=:black,lw=1)


end

plot!(collect(x),collect(x).*0,P[:,2,:],st=:surface)
