using Plots
using LinearAlgebra

function plot3dVec(x,y,z,Um,arr=.15)
    u = [x[2]-x[1]; y[2]-y[1]; z[2]-z[1]]
    S = sqrt(sum(u.^2))
    u = u/S
    p=plot!(x,y,z,lw=1,line_z=Um,leg=false)

    u′ = [x[2]-x[1]; y[2]-y[1]; 0]
    S′ = sqrt(sum(u′.^2))
    u′ = u′/S′
    θ = acos(sum([1;0;0].*u′))
    Θ = [cos(θ) sin(θ) 0;-sin(θ) cos(θ) 0;0 0 1]
    if abs(sum(u.*u′))<1
        ϕ = acos(sum(u.*u′))
    else
        ϕ = acos(1)
    end
    Φ = [cos(ϕ) 0 sin(ϕ);0 1 0;-sin(ϕ) 0 cos(ϕ)]

    seta_X = [(1-arr);1;(1-arr)]
    seta_Y = [-arr;0;arr]
    seta_Z = [0;0;0]

    New = Θ*Φ*[seta_X seta_Y seta_Z]'
    seta_X=x[1].+S*New[1,:]
    if sign(S*New[2,2])==sign(y[2]-y[1])
        seta_Y=y[1].+S*New[2,:]
    else
        seta_Y=y[1].-S*New[2,:]
    end
    if sign(S*New[3,2])==sign(z[2]-z[1])
        seta_Z=z[1].+S*New[3,:]
    else
        seta_Z=z[1].-S*New[3,:]
    end

    p=plot!(seta_X,seta_Y,seta_Z,lw=1,line_z=Um)
    return p
    # Rz = [cos(Θ) -sin(Θ) 0;sin(Θ) cos(Θ) 0;]
end

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


anim=@animate for i=1:size(P,1)
    # plot(x,y,P[:,:,i]',clim=(0,maximum(P[2:nx,2:ny,i])),st=:heatmap)
    p=plot()

    plot()
    # i=length(x)
    for j in 2:ny
        for k in 2:nz
            x1=x[i]
            y1=y[j]
            z1=z[k]
            x2=x[i] + .4*u[i,j,k]
            y2=y[j] + .4*v[i,j,k]
            z2=z[k] + .4*w[i,j,k]
            # Um = P[i,j,k]
            Um = sqrt(u[i,j,k]^2+v[i,j,k]^2+w[i,j,k]^2)
            p=plot3dVec([x1;x2],[y1;y2],[z1;z2],Um)
            # p=scatter!([x[i]],[y[j]],[z[k]],marker_z=P[i,j,k])
        end
    end

    x′ = [0;3;3;0;0;NaN;0;3;3;0;0;NaN;0;0;NaN;0;0;NaN;3;3;NaN;3;3;]
    y′ = [0;0;1;1;0;NaN;0;0;1;1;0;NaN;0;0;NaN;1;1;NaN;0;0;NaN;1;1;]
    z′ = [0;0;0;0;0;NaN;1;1;1;1;1;NaN;0;1;NaN;0;1;NaN;0;1;NaN;0;1;]
    p=plot!(x′,y′,z′,xlim=(0,4),ylim=(-2,3),zlim=(0,1),clim=(0,maximum(u)),lw=2,color=:black)
    p=plot!(cbar=false,cam=(30,60))
end


gif(anim,"teste3d.gif",fps=10)

let
    X = -15 .+ 30*rand(2000)#-3:.1:3
    Y = 1*rand(2000)#0:.05:1
    Z = 1*rand(2000)#0:.05:1
    # X,Y = meshgrid(X,Y)
    anim = @animate for k in 1:100
        p=plot(leg=false)
        x′ = [0;3;3;0;0;NaN;0;3;3;0;0;NaN;0;0;NaN;0;0;NaN;3;3;NaN;3;3;]
        y′ = [0;0;1;1;0;NaN;0;0;1;1;0;NaN;0;0;NaN;1;1;NaN;0;0;NaN;1;1;]
        z′ = [0;0;0;0;0;NaN;1;1;1;1;1;NaN;0;1;NaN;0;1;NaN;0;1;NaN;0;1;]
        p=plot!(x′,y′,z′,xlim=(0,4),ylim=(-2,3),zlim=(0,1),lw=2,color=:black)
        p=plot!(cbar=false,cam=(30,60))
        dt = .1
        Xold = X[:]; Yold= Y[:]; Zold=Z[:]
        for i in 1:length(X)
            ix = findmin(abs.(X[i] .-x))[2]
            iy = findmin(abs.(Y[i] .-y))[2]
            iz = findmin(abs.(Z[i] .-z))[2]
            Ux = u[ix,iy,iz]
            Uy = v[ix,iy,iz]
            Uz = w[ix,iy,iz]
            if X[i]<0
                Ux=1
            end
            X[i] = X[i] + dt*Ux
            Y[i] = Y[i] + dt*Uy
            Z[i] = Z[i] + dt*Uz
            if Ux==0 && Uy ==0 && Uz ==0
                X[i] = NaN
                Y[i] = NaN
                Z[i] = NaN
            end
        end
        vest = sqrt.((Xold-X).^2+(Yold-Y).^2+(Zold-Z).^2 )/dt
        scatter!(X,Y,Z,marker_z=vest,leg=false)
    end
    gif(anim,"Flow.gif",fps=15)

end
