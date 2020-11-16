include("SIMPLE-3Dcartesian-solver.jl")
include("SIMPLE-3Dcartesian-visualization.jl")

function setConditions(u,v,w,p)

    # Flow in a channel
    u[:,1,:] = fill(0,size(u[:,1,:]))# South
    u[:,end,:] = fill(0,size(u[:,end,:])) # North
    u[2,:,:] = fill(1,size(u[2,:,:])) # West
    u[:,:,1] = fill(0,size(u[:,:,1]))# Bottom
    u[:,:,end] = fill(0,size(u[:,:,end])) # Top

    v[:,1,:] = fill(0,size(v[:,1,:])) # SUL
    v[:,end,:] = fill(0,size(v[:,end,:])) # Norte
    v[1,:,:] = fill(0,size(v[1,:,:])) # West
    v[end,:,:] = fill(0,size(v[end,:,:])) # East
    u[:,:,1] = fill(0,size(u[:,:,1]))# Bottom
    u[:,:,end] = fill(0,size(u[:,:,end])) # Top
    # Outlet boundary condition (pressure)
    P[end-1,:,:] = zeros(size(u[end,:,:])) .+0 # LESTE
    u[end,(2:end-1),(2:end-1)] = u[end-1,(2:end-1),(2:end-1)] .+ dx/dy*(v[end-1,(2:end-1),(2:end-1)] - v[end-1,(3:end),(2:end-1)]) .+ dx/dz*(w[end-1,(2:end-1),(2:end-1)] - w[end-1,(2:end-1),(3:end)])

    ## Optional block
    # u[5:10,1:10] = zeros(length(u[5:10,1:10])) .+0
    # v[5:10,1:10] = zeros(length(u[5:10,1:10])) .+0
    # v[40:60,end] = zeros(length(u[40:60,end])) .-1

    return u,v,w,p
end

# Domain Discretization
nx = 20; ny = 20;  nz = 20
xmax = 5; ymax = 1;  zmax = 1
dx = xmax/nx; dy = ymax/ny;  dz = ymax/ny
x = 0:dx:xmax; y = 0:dy:ymax;  z = 0:dy:ymax

# Fluid Properties
ρ = 1
μ = 0.01

# Underelaxation properties
Ωu = .05
Ωv = .05
Ωw = .05
Ωp = .05
Ωpp = 1.2
β = 0.95

# u,v,p initialization
u = zeros(nx+1,ny+1,nz+1)
v = zeros(nx+1,ny+1,nz+1)
w = zeros(nx+1,ny+1,nz+1)
P = zeros(nx+1,ny+1,nz+1)

# u,v,p initialization
u,v,w,P=setConditions(u,v,w,P)

u,v,w,P,ϵ = Solve(u,v,w,P,50,10,100)

# Plot iteration convergence
plot(ϵ,m=4,yaxis=:log)

PlotSolution(u,v,P)
