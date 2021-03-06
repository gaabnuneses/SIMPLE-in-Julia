include("SIMPLE-3D-cylindrical-solver.jl")
using Plots

function setConditions(u,v,w,P)
    ω = 1/.4
    for i in 2:length(x), k in 2:length(z)
        if x[i]<=.4 && z[k]>= .2
            v[i,:,k] = fill(x[i]*ω,size(v[i,:,k])) # West
            w[i,:,k] = fill(0,size(v[i,:,k])) # West
            u[i,:,k] = fill(0,size(v[i,:,k])) # West
        end
    end

    v[:,ny+1,:]=v[:,ny,:]
    u[:,ny+1,:]=u[:,ny,:]
    w[:,ny+1,:]=w[:,ny,:]
    # P[:,ny+1,:]=P[:,ny,:]

    v[:,1,:]=v[:,ny+1,:]
    u[:,1,:]=u[:,ny+1,:]
    w[:,1,:]=w[:,ny+1,:]
    P[:,1,:]=P[:,ny+1,:]


    v[:,2,:]=v[:,1,:]
    u[:,2,:]=u[:,1,:]
    w[:,2,:]=w[:,1,:]
    # P[:,2,:]=P[:,1,:]

    return u,v,w,P
end

# Domain Discretization
nx = 20; ny = 30;  nz = 20
xmax = 1; ymax = 2π;  zmax = 1
xmin = 0
dx = (xmax-xmin)/nx; dy = ymax/(ny+1);  dz = zmax/nz
x = xmin:dx:xmax; y = 0:dy:(ymax-dy);  z = 0:dz:zmax

# Fluid Properties
ρ = 1
μ = 0.01

# Underelaxation properties
Ωu = .5
Ωv = .5
Ωw = .5
Ωp = .1
Ωpp = 1.7
β = 0.95

# u,v,p initialization
u = zeros(nx+1,ny+1,nz+1)
v = zeros(nx+1,ny+1,nz+1)
w = zeros(nx+1,ny+1,nz+1)
P = zeros(nx+1,ny+1,nz+1)
v[2:nx,2:ny,2:nz] = fill(1.0,size(v[2:nx,2:ny,2:nz]))
using TimerOutputs
const to = TimerOutput();
# u,v,p initialization
u,v,w,P=setConditions(u,v,w,P)
# u[2,:,:] = fill(0.004,size(u[2,:,:])) # West
# u[2,:,:] = fill(0.004,size(u[2,:,:])) # West

@timeit to "Solução" begin
    u,v,w,P,ϵ = Solve(u,v,w,P,10,10,20)
end
# Plot iteration convergence
plot(ϵ,m=4,yaxis=:log)

plot(y,v[6,:,2],m=2)
plot()
for k in 1:length(z)
    plot!(x,v[:,6,k],ylim=(0,1),m=2)
end
plot!()

iz = Int(round(nz/2))
plot(x,z,v[:,iz,:]',st=:heatmap)
for iz = 1:nz
    # PlotSolution(u[:,:,iz],v[:,:,iz],u[:,:,iz])
end
k = Int(round(nz/2))
PlotSolution(u[:,:,k],v[:,:,k],P[:,:,k])

let
    X = [];     Z = [];     Y = [];
    dX = [];    dZ = [];    dY = [];
    for i in 1:nx, j in 1:ny+1
        X = [X; x[i]*cos(y[j])]
        Y = [Y; x[i]*sin(y[j])]
        Z = [Z; P[i,j,iz]]
        dX = [dX;u[i,j,iz]*cos(y[j])-v[i,j,iz]*sin(y[j])]
        dY = [dY;u[i,j,iz]*sin(y[j])+v[i,j,iz]*cos(y[j])]
    end
    plot(X,Y,Z,st=:surface,cam=(0,90),aspect_ratio=:equal,leg=false)
###
    # plot!(Shape(X,Y),color=:black,xlim=(-1,1),ylim=(-1,1))
    # for i in 1:nx, j in 1:ny+1
    #     X = [X; x[i]*cos(y[j])]
    #     Y = [Y; x[i]*sin(y[j])]
    #     Z = [Z; P[i,j,iz]]
    #     dX = [dX;u[i,j,iz]*cos(y[j])-v[i,j,iz]*sin(y[j])]
    #     dY = [dY;u[i,j,iz]*sin(y[j])+v[i,j,iz]*cos(y[j])]
    # end

###
    for i in 1:length(X)
        if sqrt(X[i]^2+Y[i]^2)>.4
            Plot2DVector([X[i];X[i]+.1*dX[i]],[Y[i];Y[i]+.1*dY[i]])
        end
    end
    X = [.4*cos(i) for i in 0:y[2]:y[end]]
    Y = [.4*sin(i) for i in 0:y[2]:y[end]]
    # plot!(Shape(X,Y),color=:white,xlim=(-1,1),ylim=(-1,1))
    plot!(xlim=(-1,1),ylim=(-1,1),cbar=true)
end

let
    X = [];     Z = [];     Y = [];
    dX = [];    dZ = [];    dY = [];
    for i in 1:length(x), k in 1:length(z)
        X = [X; x[i]]
        Y = [Y; z[k]]
        Z = [Z; P[i,5,k]]
        dX = [dX;u[i,5,k]]
        dY = [dY;w[i,5,k]]
    end
    plot(X,Y,Z,st=:surface,cam=(0,90),leg=false)
    for i in 1:length(X)
        Plot2DVector([X[i];X[i]+.0010*dX[i]],[Y[i];Y[i]+.0010*dY[i]])
    end

    X = [];     Z = [];     Y = [];
    dX = [];    dZ = [];    dY = [];
    for i in 1:length(x), k in 1:length(z)
        X = [X; x[i]]
        Y = [Y; z[k]]
        Z = [Z; P[i,20,k]]
        dX = [dX;u[i,20,k]]
        dY = [dY;w[i,20,k]]
    end
    plot!(-X,Y,Z,st=:surface,cam=(0,90),leg=false)
    # plot(X,Y,Z,st=:surface,cam=(0,90),leg=false)
    for i in 1:length(X)
        Plot2DVector(-[X[i];X[i]+.001*dX[i]],[Y[i];Y[i]+.001*dY[i]])
    end

    #
    # X = [x[1]*cos(i) for i in 0:y[2]:y[end]]
    # Y = [x[1]*sin(i) for i in 0:y[2]:y[end]]
    # plot!(Shape(X,Y),color=:white)
    X = [-.4;.4;.4;-.4]
    Y = [.2;.2;z[end];z[end]]
    plot!(Shape(X,Y),color=:white)

end

plot(xlim=(-1,1),ylim=(-1,1))
for i in 1:length(y)
    annotate!(cos(y[i]),sin(y[i]),"$i")
end
plot!(aspect_ratio=:equal)
