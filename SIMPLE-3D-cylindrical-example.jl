include("SIMPLE-3D-cylindrical-solver.jl")
using Plots

function setConditions(u,v,w,P)
    # w[2,:,:] = fill(0.0,size(u[1,:,:])) # West
    # v[2,:,:] = fill(0.0,size(u[1,:,:])) # West
    # P[nx,:,:] = zeros(size(u[nx+1,:,:])) .+0 # LESTE
    # u[10:11,:,1:12] = zeros(size(u[10:11,:,1:12]))
    # v[10:11,:,1:12] = zeros(size(u[10:11,:,1:12]))
    # w[10:11,:,1:12] = zeros(size(u[10:11,:,1:12]))
    for j in 2:ny, k in 2:nz
        u[nx+1,j,k] = u[nx,j,k] #.+ dx*dz/dy*(v[nx,j,k] - v[nx,j+1,k]) .+ dx*dy/dz*(w[nx,j,k] - w[nx,j,k+1])
    end
    return u,v,w,P
end

# Domain Discretization
nx = 10; ny = 30;  nz = 10
xmax = .02; ymax = 2π;  zmax = .1
xmin = .01
dx = (xmax-xmin)/nx; dy = ymax/ny;  dz = ymax/ny
x = xmin:dx:xmax; y = 0:dy:ymax;  z = 0:dy:ymax

######### VISUALIZAÇÃO
# phi = 0:pi/10:2*pi;
# theta = 0:pi/10:pi;
let
    X = []
    Y = []
    Z = []
    for i in 1:nx+1, j in 1:ny+1, k in 1:nz+1
        X = [X;x[i]*cos(y[j])]
        Y = [Y;x[i]*sin(y[j])]
        Z = [Z;z[k]]
    end
    scatter(X,Y,Z,cam=(30,60),aspect_ratio=:equal)
end
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
v[2,:,:] = fill(1,size(u[1,:,:])) # West
v[1,:,:] = fill(1,size(u[1,:,:])) # West

@timeit to "Solução" begin
    u,v,w,P,ϵ = Solve(u,v,w,P,100,100,100)
end
# Plot iteration convergence
plot(ϵ,m=4,yaxis=:log)
