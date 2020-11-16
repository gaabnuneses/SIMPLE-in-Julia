include("SIMPLE-2Dcartesian-solver.jl")
include("SIMPLE-2Dcartesian-visualization.jl")
include("Unsteady-SIMPLE-2Dcartesian-solver.jl.jl")
include("Unsteady-SIMPLE-2Dcartesian-GIF.jl")

function setConditions(u,v,p)

    # Flow in a channel
    u[:,1] = fill(0,length(u[:,1]))# South
    u[:,end] = fill(0,length(u[:,end])) # North
    u[2,:] = fill(1,length(u[2,:])) # West

    v[:,1] = fill(0,length(v[:,1])) # SUL
    v[:,end] = fill(0,length(v[:,end])) # Norte
    v[1,:] = fill(0,length(v[1,:])) # West
    v[end,:] = fill(0,length(v[end,:])) # East

    # Outlet boundary condition (pressure)
    P[end-1,:] = zeros(length(u[end,:])) .+0 # LESTE
    u[end,(2:end-1)] = u[end-1,(2:end-1)] .+ dx/dy*(v[end-1,(2:end-1)] - v[end-1,(3:end)])

    ## Optional block
    # u[5:10,1:10] = zeros(length(u[5:10,1:10])) .+0
    # v[5:10,1:10] = zeros(length(u[5:10,1:10])) .+0
    # v[40:60,end] = zeros(length(u[40:60,end])) .-1

    return u,v,p
end

# Domain Discretization
nx = 100; ny = 20
xmax = 5; ymax = 1
dx = xmax/nx; dy = ymax/ny
x = 0:dx:xmax; y = 0:dy:ymax
# Time Discretization
nt = 1000
tmax = 1
dt = tmax/nt
t = 0:dt:tmax

# Fluid Properties
ρ = 1
μ = 0.01

# Underelaxation properties
Ωu = .1
Ωv = .1
Ωp = .1
β = 0.95

# u,v,p initialization
u = zeros(nx+1,ny+1)
v = zeros(nx+1,ny+1)
P = zeros(nx+1,ny+1)

# u,v,p initialization
u,v,P=setConditions(u,v,P)

u,v,P,ϵ = Solve(u,v,P,500,10,100)
Uk,Vk,Pk = Unsteady(u,v,P,5,5,40)
makeGIF(x,y,Uk,Vk,Pk)
