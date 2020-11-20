include("SIMPLE-3D-cylindrical-solver.jl")
using Plots

function setConditions(u,v,w,P)
    ω = 1/.4
    for i in 1:length(x), k in 1:length(z)
        if x[i]<.4 && z[k]> .2
            v[i,:,k] = fill(x[i]*ω,size(v[i,:,k])) # West
        end
        # v[i,:,:] = x[i]fill(1,size(u[1,:,:])) # West
    end

    # w[2,:,:] = fill(0.0,size(u[1,:,:])) # West
    # v[2,:,:] = fill(0.0,size(u[1,:,:])) # West
    # P[nx,:,:] = zeros(size(u[nx+1,:,:])) .+0 # LESTE
    # u[10:11,:,1:12] = zeros(size(u[10:11,:,1:12]))
    # v[10:11,:,1:12] = zeros(size(u[10:11,:,1:12]))
    # w[10:11,:,1:12] = zeros(size(u[10:11,:,1:12]))
    #
    u[:,2,:]=u[:,end,:]
    v[:,2,:]=v[:,end,:]
    w[:,2,:]=w[:,end,:]
    # P[:,1,:]=P[:,end,:]
    # u[:,:,1]=u[:,:,end]
    # v[:,:,1]=v[:,:,end]
    # w[:,:,1]=w[:,:,end]
    # for j in 2:ny, k in 2:nz
    #     u[nx+1,j,k] = u[nx,j,k] #.+ dx*dz/dy*(v[nx,j,k] - v[nx,j+1,k]) .+ dx*dy/dz*(w[nx,j,k] - w[nx,j,k+1])
    # end
    return u,v,w,P
end

plot(y,uvo[5,:,5],m=2)
# Domain Discretization
nx = 10; ny = 30;  nz = 10
xmax = 1; ymax = 2π;  zmax = 1
xmin = 0
dx = (xmax-xmin)/nx; dy = ymax/ny;  dz = zmax/nz
x = xmin:dx:xmax; y = 0:dy:ymax;  z = 0:dz:zmax

# Fluid Properties
ρ = 1
μ = 0.01

# Underelaxation properties
Ωu = .5
Ωv = .5
Ωw = .5
Ωp = .3
Ωpp = 1.7
β = 0.95

# u,v,p initialization
u = zeros(nx+1,ny+1,nz+1)
v = zeros(nx+1,ny+1,nz+1)
w = zeros(nx+1,ny+1,nz+1)
P = zeros(nx+1,ny+1,nz+1)

using TimerOutputs
const to = TimerOutput();
# u,v,p initialization
u,v,w,P=setConditions(u,v,w,P)
# u[2,:,:] = fill(0.004,size(u[2,:,:])) # West
# u[2,:,:] = fill(0.004,size(u[2,:,:])) # West

@timeit to "Solução" begin
    u,v,w,P,ϵ = Solve(u,v,w,P,25,100,100)
end
# Plot iteration convergence
plot(ϵ,m=4,yaxis=:log)

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
    for i in 1:length(x), j in 1:length(y)
        X = [X; x[i]*cos(y[j])]
        Y = [Y; x[i]*sin(y[j])]
        Z = [Z; P[i,j,iz]]
        dX = [dX;u[i,j,k]*cos(y[j])-v[i,j,k]*sin(y[j])]
        dY = [dY;u[i,j,k]*sin(y[j])+v[i,j,k]*cos(y[j])]
    end
    plot(X,Y,Z,st=:surface,cam=(0,90),aspect_ratio=:equal,leg=false)
    for i in 1:length(X)
        Plot2DVector([X;X+.0000001*dX],[Y;Y+.0000001*dY])
    end

    X = [x[1]*cos(i) for i in 0:y[2]:y[end]]
    Y = [x[1]*sin(i) for i in 0:y[2]:y[end]]
    plot!(Shape(X,Y),color=:white)

end

let
    X = [];     Z = [];     Y = [];
    dX = [];    dZ = [];    dY = [];
    for i in 1:length(x), k in 1:length(z)
        X = [X; x[i]]
        Y = [Y; z[k]]
        Z = [Z; P[i,13,k]]
        dX = [dX;u[i,13,k]]
        dY = [dY;w[i,13,k]]
    end
    plot(X,Y,Z,st=:surface,cam=(0,90),leg=false)

    X = [];     Z = [];     Y = [];
    dX = [];    dZ = [];    dY = [];
    for i in 1:length(x), k in 1:length(z)
        X = [X; x[i]]
        Y = [Y; z[k]]
        Z = [Z; v[i,2,k]]
        dX = [dX;u[i,1,k]]
        dY = [dY;w[i,1,k]]
    end
    plot!(-X,Y,Z,st=:surface,cam=(0,90),leg=false)
    # for i in 1:length(X)
    #     Plot2DVector([X;X+.0000001*dX],[Y;Y+.0000001*dY])
    # end
    #
    # X = [x[1]*cos(i) for i in 0:y[2]:y[end]]
    # Y = [x[1]*sin(i) for i in 0:y[2]:y[end]]
    # plot!(Shape(X,Y),color=:white)

end
