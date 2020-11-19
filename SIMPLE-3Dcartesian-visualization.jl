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
    for i in 2:1:nx,j in 2:1:ny
        ux  = x[2]*u[i,j]/maxVel
        vy  = x[2]*v[i,j]/maxVel
        V = sqrt(u[i,j]^2+v[i,j]^2)
        # p=plot!([x[i];x[i]+ux],[y[j];y[j]+vy],arrow=true,leg=false)
        p=Plot2DVector([x[i];x[i]+ux],[y[j];y[j]+vy],V)
    end
    return p
end

function PlotSolution(u,v,P)
    p=plot()
    # Pressure Field
    p=plot!(x[1:end],y[1:end],P[1:end,1:end]',st=:heatmap,cam=(0,90))
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
    y′ = [0;0;2;2;0;NaN;0;0;2;2;0;NaN;0;0;NaN;2;2;NaN;0;0;NaN;2;2;]
    z′ = [0;0;0;0;0;NaN;2;2;2;2;2;NaN;0;2;NaN;0;2;NaN;0;2;NaN;0;2;]
    p=plot!(x′,y′,z′,xlim=(0,4),ylim=(-2,4),zlim=(0,2),clim=(0,maximum(u)),lw=2,color=:black)

    # x′ = [1;2;2;1;1;NaN;1;2;2;1;1;NaN;1;1;NaN;1;1;NaN;2;2;NaN;2;2;]
    # y′ = [0;0;1;1;0;NaN;0;0;1;1;0;NaN;0;0;NaN;1;1;NaN;0;0;NaN;1;1;]
    # z′ = [0;0;0;0;0;NaN;.5;.5;.5;.5;.5;NaN;0;.5;NaN;0;.5;NaN;0;.5;NaN;0;.5;]
    # p=plot!(x′,y′,z′,xlim=(0,4),ylim=(-1,2),zlim=(0,1),clim=(0,maximum(u)),lw=2,color=:black)

    phi = 0:pi/10:2*pi;
    theta = 0:pi/10:pi;
    xs = [.25*cos(t)*sin(p) for t in theta, p in phi] .+ 1.5 ;
    ys = [.25*sin(t)*sin(p) for t in theta, p in phi] .+ 1 ;
    zs = [.25*cos(p) for t in theta, p in phi] .+ 1 ;
    plot!(xs,ys,zs,color=:black)
    plot!(xs',ys',zs',color=:black)
    p=plot!(cbar=false,cam=(30,60))
end

gif(anim,"teste3d.gif",fps=15)

let
    X = -15 .+ 30*rand(7000)#-3:.1:3
    Y = 2*rand(7000)#0:.05:1
    Z = 2*rand(7000)#0:.05:1
    # X,Y = meshgrid(X,Y)
    anim = @animate for k in 1:400
        p=plot(leg=false)
        x′ = [0;3;3;0;0;NaN;0;3;3;0;0;NaN;0;0;NaN;0;0;NaN;3;3;NaN;3;3;]
        y′ = [0;0;2;2;0;NaN;0;0;2;2;0;NaN;0;0;NaN;2;2;NaN;0;0;NaN;2;2;]
        z′ = [0;0;0;0;0;NaN;2;2;2;2;2;NaN;0;2;NaN;0;2;NaN;0;2;NaN;0;2;]
        p=plot!(x′,y′,z′,xlim=(0,4),ylim=(-2,4),zlim=(0,2),lw=2,color=:black)
        # x′ = [1;2;2;1;1;NaN;1;2;2;1;1;NaN;1;1;NaN;1;1;NaN;2;2;NaN;2;2;]
        # y′ = [0;0;1;1;0;NaN;0;0;1;1;0;NaN;0;0;NaN;1;1;NaN;0;0;NaN;1;1;]
        # z′ = [0;0;0;0;0;NaN;.5;.5;.5;.5;.5;NaN;0;.5;NaN;0;.5;NaN;0;.5;NaN;0;.5;]
        # p=plot!(x′,y′,z′,xlim=(0,4),ylim=(-1,2),zlim=(0,1),clim=(0,maximum(u)),lw=2,color=:black)
        phi = 0:pi/10:2*pi;
        theta = 0:pi/10:pi;
        xs = [.5*cos(t)*sin(p) for t in theta, p in phi] .+ 1.5 ;
        ys = [.5*sin(t)*sin(p) for t in theta, p in phi] .+ 1 ;
        zs = [.5*cos(p) for t in theta, p in phi] .+ 1 ;
        plot!(xs,ys,zs,color=:black)
        plot!(xs',ys',zs',color=:black)

        p=plot!(cbar=false,cam=(30,60))
        dt = .05
        Vel = zeros(length(X))
        for i in 1:length(X)
            if isnan(X[i])==false
                ix = findmin(abs.(X[i] .-x))[2]
                iy = findmin(abs.(Y[i] .-y))[2]
                iz = findmin(abs.(Z[i] .-z))[2]
                Ux = u[ix,iy,iz]
                Uy = v[ix,iy,iz]
                Uz = w[ix,iy,iz]
                if X[i]<0
                    Ux=1
                    Uy=0
                    Uz=0
                end
                X[i] = X[i] + dt*Ux
                Y[i] = Y[i] + dt*Uy
                Z[i] = Z[i] + dt*Uz
                Vel[i] = sqrt(Ux^2+Uy^2+Uz^2)

                if Ux==0 && Uy ==0 && Uz ==0
                    X[i] = 0
                    Y[i] = 2*rand()
                    Z[i] = 2*rand()
                    # Vel[i] = Inf
                end
                if X[i]>3.5 || Y[i]>2 || Y[i]<0 || Z[i]<0 || Z[i]>2
                    X[i] = 0
                    Y[i] = 2*rand()
                    Z[i] = 2*rand()
                    # Vel[i] = Inf
                end
            end
        end
        coef=((X.<0).-1).^(-1).+2
        scatter!(coef.*X,coef.*Y,coef.*Z,m = (3.0, :blue, 0.5, stroke(0, :blue)),leg=false)
    end
    gif(anim,"Flow.gif",fps=15)

end

u[1,:,:]=u[2,:,:]

let
    p=plot(leg=false)
    x′ = [0;30;30;0;0;NaN;0;30;30;0;0;NaN;0;0;NaN;0;0;NaN;30;30;NaN;30;30;]
    y′ = [0;0;20;20;0;NaN;0;0;20;20;0;NaN;0;0;NaN;20;20;NaN;0;0;NaN;20;20;]
    z′ = [0;0;0;0;0;NaN;20;20;20;20;20;NaN;0;20;NaN;0;20;NaN;0;20;NaN;0;20;]
    p=plot!(x′,y′,z′,xlim=(0,40),ylim=(-20,40),zlim=(0,20),lw=2,color=:black)
    # x′ = [1;2;2;1;1;NaN;1;2;2;1;1;NaN;1;1;NaN;1;1;NaN;2;2;NaN;2;2;]
    # y′ = [0;0;1;1;0;NaN;0;0;1;1;0;NaN;0;0;NaN;1;1;NaN;0;0;NaN;1;1;]
    # z′ = [0;0;0;0;0;NaN;.5;.5;.5;.5;.5;NaN;0;.5;NaN;0;.5;NaN;0;.5;NaN;0;.5;]
    # p=plot!(x′,y′,z′,xlim=(0,4),ylim=(-1,2),zlim=(0,1),clim=(0,maximum(u)),lw=2,color=:black)
    phi = 0:pi/10:2*pi;
    theta = 0:pi/10:pi;
    xs = [5*cos(t)*sin(p) for t in theta, p in phi] .+ 15 ;
    ys = [5*sin(t)*sin(p) for t in theta, p in phi] .+ 10 ;
    zs = [5*cos(p) for t in theta, p in phi] .+ 10 ;
    plot!(xs,ys,zs,color=:black)
    plot!(xs',ys',zs',color=:black)

    p=plot!(cbar=false,cam=(30,60))

    dt = .01
    T = 2500
    n = 100
    X = zeros(1,n)#-3:.1:3
    Y = 20*rand(1,n)#0:.05:1
    Z = 20*rand(1,n)#0:.05:1
    Pre = zeros(T,n+1)
    for t = 1:T
        xn = []; yn = []; zn =[]; vn =[]
        for i in 1:n
            ix = 1; iy=1; iz=1
            for j in 1:nx
                if X[end,i]>=x[j] && X[end,i]<x[j+1]
                    ix=j
                    break
                else X[end,i]>=x[j+1]
                    ix=j
                end
            end
            for j in 1:ny
                if Y[end,i]>=y[j] && Y[end,i]<y[j+1]
                    iy=j
                    break
                else Y[end,i]>=y[j+1]
                    iy=j
                end
            end
            for j in 1:nz
                if Z[end,i]>=z[j] && Z[end,i]<z[j+1]
                    iz=j
                    break
                else Z[end,i]>=z[j+1]
                    iz=j
                end
            end
            Ux = u[ix,iy,iz]
            Uy = v[ix,iy,iz]
            Uz = w[ix,iy,iz]

            Pre[t,i] = sqrt(Ux^2+Uy^2+Uz^2)#P[ix,iy,iz]
            if X[i]<0
                Ux=1
                Uy=0
                Uz=0
            end
            xn = [xn; X[end,i] + dt*Ux]
            yn = [yn; Y[end,i] + dt*Uy]
            zn = [zn; Z[end,i] + dt*Uz]
            vn = [vn; sqrt(Ux^2+Uy^2+Uz^2)]
            # if Ux==0 && Uy ==0 && Uz ==0
            #     X[i] = NaN
            #     Y[i] = NaN
            #     Z[i] = NaN
            #     Vel[i] = NaN
            # end
            # if X[i]>3.5 || Y[i]>1 || Y[i]<0 || Z[i]<0 || Z[i]>1
            #     X[i] = NaN
            #     Y[i] = NaN
            #     Z[i] = NaN
            #     Vel[i] = NaN
            # end
        end
        X = [X;xn']
        Y = [Y;yn']
        Z = [Z;zn']
        # Pre = [Pre;Pr']
    end
    display(size(X))
    display(size(Pre))
    # display(Vel)
    display(Pre)
    for i in 1:size(X,2)
        p=plot!(X[:,i],Y[:,i],Z[:,i],color=:inferno,line_z=Pre[:,i],lw=2,clim=(minimum(Pre),maximum(Pre)))
    end
    scatter!([X[1,:];X[end,:]],[Y[1,:];Y[end,:]],[Z[1,:];Z[end,:]],m=2,color=:inferno,marker_z=[Pre[1,:];Pre[end,:]])
    p
end
plot!(xs,ys,zs,color=:black,ls=:dash)
plot!(xs',ys',zs',color=:black,ls=:dash)
savefig("corrente.png")
P[:,1,:]=P[:,2,:]
P[:,:,1]=P[:,:,2]
let
    wall = zeros(0,3)
    for i = 1:nx, j = 1:ny, k = 1:nz
        if u[i,j,k]==0 && v[i,j,k]==0 && w[i,j,k]==0
            wall = [wall;i j k]
        end
    end
    wall = Int.(wall)
    display(wall)
    ix = wall[:,1]
    iy = wall[:,2]
    iz = wall[:,3]
    scatter(x[ix],y[iy],z[iz])
end
