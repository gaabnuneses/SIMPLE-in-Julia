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
    # u[:,2,:]=u[:,ny,:]
    # u[:,2,:]=u[:,ny+1,:]
    # w[:,2,:]=w[:,ny+1,:]
    v[:,2,:]=v[:,ny+1,:]
    # P[:,2,:]=P[:,ny+1,:]



    # for i in 3:nx, j in 2:ny, k in 2:nz
    #     ae = (x[i+1]+x[i])/2 * (y[j+1]-y[j-1])/2*(z[k+1]-z[k-1])/2
    #     aw = (x[i-1]+x[i])/2 * (y[j+1]-y[j-1])/2*(z[k+1]-z[k-1])/2
    #     u[i,j,k]=u[i,j,k] + 1/Apu[i,j,k]*(aw*Pp[i-1,j,k]-ae*Pp[i,j,k])
    # end
    # # Corrigindo v
    # for i in 2:nx, j in 3:ny, k in 2:nz
    #     an = (x[i+1]-x[i-1])/2*(z[k+1]-z[k-1])/2
    #     as = (x[i+1]-x[i-1])/2*(z[k+1]-z[k-1])/2
    #     v[i,j,k]=v[i,j,k] + 1/Apv[i,j,k]*(as*Pp[i,j-1,k]-an*Pp[i,j,k])
    # end
    # # Corrigindo w
    # for i in 2:nx, j in 2:ny, k in 3:nz
    #     at = (y[j+1]-y[j-1])/2*(((x[i+1]-x[i])/2)^2-((x[i]-x[i-1])/2)^2)/2
    #     ab = (y[j+1]-y[j-1])/2*(((x[i+1]-x[i])/2)^2-((x[i]-x[i-1])/2)^2)/2
    #     w[i,j,k]=w[i,j,k] + 1/Apw[i,j,k]*(ab*Pp[i,j,k-1]-at*Pp[i,j,k])
    # end

    # w[:,2,:]=w[:,ny,:]
    # P[:,1,:]=P[:,end,:]
    # u[:,:,1]=u[:,:,end]
    # v[:,:,1]=v[:,:,end]
    # w[:,:,1]=w[:,:,end]
    # for j in 2:ny, k in 2:nz
    #     u[nx+1,j,k] = u[nx,j,k] #.+ dx*dz/dy*(v[nx,j,k] - v[nx,j+1,k]) .+ dx*dy/dz*(w[nx,j,k] - w[nx,j,k+1])
    # end
    return u,v,w,P
end

# Domain Discretization
nx = 20; ny = 30;  nz = 20
xmax = 1; ymax = 2π;  zmax = 1
xmin = .4
dx = (xmax-xmin)/nx; dy = ymax/ny;  dz = zmax/nz
x = xmin:dx:xmax; y = 0:dy:ymax;  z = 0:dz:zmax

# Fluid Properties
ρ = 1
μ = 0.01

# Underelaxation properties
Ωu = .05
Ωv = .05
Ωw = .05
Ωp = .03
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
    u,v,w,P,ϵ = Solve(u,v,w,P,25,10,120)
end
# Plot iteration convergence
plot(ϵ,m=4,yaxis=:log)

plot(y,v[6,:,5],m=2)

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
    for i in 1:nx, j in 1:ny
        X = [X; x[i]*cos(y[j])]
        Y = [Y; x[i]*sin(y[j])]
        Z = [Z; P[i,j,iz]]
        dX = [dX;u[i,j,iz]*cos(y[j])-v[i,j,iz]*sin(y[j])]
        dY = [dY;u[i,j,iz]*sin(y[j])+v[i,j,iz]*cos(y[j])]
    end
    plot(X,Y,Z,st=:surface,cam=(0,90),aspect_ratio=:equal,leg=false)
    # for i in 1:length(X)
    #     Plot2DVector([X;X+.0000001*dX],[Y;Y+.0000001*dY])
    # end

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
        Z = [Z; P[i,16,k]]
        dX = [dX;u[i,16,k]]
        dY = [dY;w[i,16,k]]
    end
    plot(X,Y,Z,st=:surface,cam=(0,90),leg=false)

    X = [];     Z = [];     Y = [];
    dX = [];    dZ = [];    dY = [];
    for i in 1:length(x), k in 1:length(z)
        X = [X; x[i]]
        Y = [Y; z[k]]
        Z = [Z; P[i,1,k]]
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
