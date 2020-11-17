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

    # p=plot!(seta_X,seta_Y,seta_Z,aspect_ratio=:equal,lw=1,line_z=Um)
    return p
    # Rz = [cos(Θ) -sin(Θ) 0;sin(Θ) cos(Θ) 0;]
end


anim=@animate for i=1:size(P,3)
    # plot(x,y,P[:,:,i]',clim=(0,maximum(P[2:nx,2:ny,i])),st=:heatmap)
    p=plot()

    plot()
    # i=length(x)
    for j in 2:ny
        for k in 2:nz
            x1=x[i]
            y1=y[j]
            z1=z[k]
            x2=x[i] + u[i,j,k]
            y2=y[j] + v[i,j,k]
            z2=z[k] + w[i,j,k]
            Um = P[i,j,k]#sqrt(u[i,j,k]^2+v[i,j,k]^2+w[i,j,k]^2)
            p=plot3dVec([x1;x2],[y1;y2],[z1;z2],Um)
            # p=scatter!([x[i]],[y[j]],[z[k]],marker_z=P[i,j,k])
        end
    end

    x′ = [0;5;5;0;0;NaN;0;5;5;0;0;NaN;0;0;NaN;0;0;NaN;5;5;NaN;5;5;]
    y′ = [0;0;1;1;0;NaN;0;0;1;1;0;NaN;0;0;NaN;1;1;NaN;0;0;NaN;1;1;]
    z′ = [0;0;0;0;0;NaN;1;1;1;1;1;NaN;0;1;NaN;0;1;NaN;0;1;NaN;0;1;]
    p=plot!(x′,y′,z′,xlim=(0,7),ylim=(-2,3),zlim=(0,1),clim=(0,maximum(P)),lw=2,color=:black)
    p=plot!(cbar=false,cam=(30,60))
end


gif(anim,"teste3d.gif",fps=10)

p=plot!()
