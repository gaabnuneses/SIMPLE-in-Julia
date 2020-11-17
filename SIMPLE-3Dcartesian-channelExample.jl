include("SIMPLE-3Dcartesian-solver.jl")
include("SIMPLE-3Dcartesian-visualization.jl")

function setConditions(u,v,w,p)


    # Flow in a channel
    u[:,1,:] = fill(0.0,size(u[:,1,:]))# South
    u[:,end,:] = fill(0.0,size(u[:,end,:])) # North
    u[2,:,:] = fill(1.0,size(u[2,:,:])) # West
    u[:,:,1] = fill(0.0,size(u[:,:,1]))# Bottom
    u[:,:,end] = fill(0.0,size(u[:,:,end])) # Top

    # v[:,1,:] = fill(0,size(v[:,1,:])) # SUL
    # v[:,end,:] = fill(0,size(v[:,end,:])) # Norte
    # v[1,:,:] = fill(0,size(v[1,:,:])) # West
    # v[end,:,:] = fill(0,size(v[end,:,:])) # East
    # v[:,:,1] = fill(0,size(u[:,:,1]))# Bottom
    # v[:,:,end] = fill(0,size(u[:,:,end])) # Top
    #
    # w[:,1,:] = fill(0,size(v[:,1,:])) # SUL
    # w[:,end,:] = fill(0,size(v[:,end,:])) # Norte
    # w[1,:,:] = fill(0,size(v[1,:,:])) # West
    # w[end,:,:] = fill(0,size(v[end,:,:])) # East
    # w[:,:,1] = fill(0,size(u[:,:,1]))# Bottom
    # w[:,:,end] = fill(0,size(u[:,:,end])) # Top

    # Outlet boundary condition (pressure)
    P[2,:,:] = zeros(size(P[2,:,:]))  # LESTE
    P[nx,:,:] = zeros(size(u[nx+1,:,:])) .+0 # LESTE
    # u[nx,2:ny,2:nz] = u[nx,2:ny,2:nz] .+ dx*dz/dy*(v[nx,2:ny,2:nz] - v[nx,3:(ny+1),3:(nz+1)]) .+ dx*dy/dz*(w[nx,2:ny,2:nz] - w[nx,3:(ny+1),3:(nz+1)])
    u[nx+1,2:ny,2:nz] = u[nx,2:ny,2:nz] .+ dx/dy*(v[nx,2:ny,2:nz] - v[nx,3:(ny+1),2:nz]) .+ dx/dz*(w[nx,2:ny,2:nz] - w[nx,2:ny,3:(nz+1)])
    ## Optional block
    # u[5:10,1:10] = zeros(length(u[5:10,1:10])) .+0
    # v[5:10,1:10] = zeros(length(u[5:10,1:10])) .+0
    # v[40:60,end] = zeros(length(u[40:60,end])) .-1

    return u,v,w,p
end

# Domain Discretization
nx = 90; ny = 30;  nz = 30
xmax = 3; ymax = 1;  zmax = 1
dx = xmax/nx; dy = ymax/ny;  dz = ymax/ny
x = 0:dx:xmax; y = 0:dy:ymax;  z = 0:dy:ymax

# Fluid Properties
ρ = 1
μ = 0.01

# Underelaxation properties
Ωu = .1
Ωv = .1
Ωw = .1
Ωp = .1
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
# u,v,w,P=setConditions(u,v,w,P)
u[2,:,:] = fill(1.0,size(u[2,:,:])) # West

@timeit to "Solução" begin
    u,v,w,P,ϵ = Solve(u,v,w,P,150,10,100)
end
# Plot iteration convergence
plot(ϵ,m=4,yaxis=:log)

anim = @animate for k in 1:size(u,2)
    PlotSolution(u[:,k,:],w[:,k,:],P[:,k,:])
    plot!(clim=(minimum(P),maximum(P)))
end
gif(anim,"uw.gif",fps=2)

anim = @animate for k in 1:size(u,2)
    PlotSolution(u[:,:,k],v[:,:,k],P[:,:,k])
    plot!(clim=(minimum(P),maximum(P)))
end
gif(anim,"uv.gif",fps=2)

PlotSolution(v[16,:,:],w[16,:,:],P[16,:,:])

plot(y,z,u[end,:,:],aspect_ratio=:equal,xlim=(0,1),ylim=(0,1),title="Perfil de velocidades",xlabel="y[m]",ylabel="z[m]",cam=(60,30),st=:heatmap)
savefig("perfilVelocidade.png")

plot(y,z,u[end,:,:],cam=(60,30),st=:wireframe)
