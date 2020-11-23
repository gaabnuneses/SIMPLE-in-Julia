using Plots

function Plot2DVector(x,y,V,arr=.15)
    plot!(x,y,lw=1,color=:cyan)#line_z=V)
    seta_X = [(1-arr);1;(1-arr)]
    seta_Y = [-arr;0;arr]
    S = sqrt((x[1]-x[2])^2+(y[1]-y[2])^2)
    Θ = acos((x[2]-x[1])/S)
    # Θ = atan((y[2]-y[1])/(x[2]-x[1]))
    R = [cos(Θ) -sin(Θ);sin(Θ) cos(Θ)]

    New = R*[seta_X seta_Y]'
    # display(New)
    seta_X=x[1].+S*New[1,:]
    if (y[1].+S*New[2,2])!=y[2]
        S = -S
    end
    seta_Y=y[1].+S*New[2,:]
    plot!(seta_X,seta_Y,aspect_ratio=:equal,lw=1,color=:cyan)#line_z=V)
end

function PlotVelocityField(x,y,u,v)
    maxVel = maximum(sqrt.(u.^2+v.^2))
    p=plot!(leg=false)
    for i in 2:2:nx,j in 2:2:ny
        ux  = .1*u[i,j]/maxVel
        vy  = .1*v[i,j]/maxVel
        V = sqrt(u[i,j]^2+v[i,j]^2)
        # p=plot!([x[i];x[i]+ux],[y[j];y[j]+vy],arrow=true,leg=false)
        p=Plot2DVector([x[i];x[i]+ux],[y[j];y[j]+vy],V)
    end
    return p
end

function PlotSolution(u,v,P)
    p=plot()
    # Pressure Field
    p=plot!(x[2:end-1],y[2:end-1],P[2:end-1,2:end-1]',st=:heatmap,cam=(0,90))
    # Velocity Field
    p=PlotVelocityField(x,y,u,v)
    p=plot!(xlim=(minimum(x),maximum(x)),ylim=(minimum(y),maximum(y)))
    display(p)
end
