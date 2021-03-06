include("SIMPLE-3Dcartesian-solver.jl")
include("SIMPLE-3Dcartesian-visualization.jl")
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
nx = 2*30; ny = 2*10;  nz = 2*10
xmax = 3; ymax = 1;  zmax = 1
dx = xmax/nx; dy = ymax/ny;  dz = ymax/ny
x = 0:dx:xmax; y = 0:dy:ymax;  z = 0:dy:ymax

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
u[2,:,:] = fill(1,size(u[1,:,:])) # West
u[1,:,:] = fill(1,size(u[1,:,:])) # West

@timeit to "Solução" begin
    u,v,w,P,ϵ = Solve(u,v,w,P,100,100,100)
end
# Plot iteration convergence
plot(ϵ,m=4,yaxis=:log)

iy = Int(round(ny/2))
PlotSolution(u[:,iy,:],w[:,iy,:],P[:,iy,:])
k = Int(round(nz/2))
PlotSolution(u[:,:,k],v[:,:,k],P[:,:,k])
anim = @animate for k in 1:size(u,2)
    PlotSolution(u[:,k,:],w[:,k,:],P[:,k,:])
    plot!(clim=(minimum(P),maximum(P)))
end
gif(anim,"uw.gif",fps=5)

anim = @animate for k in 1:size(u,2)
    PlotSolution(u[:,:,k],v[:,:,k],P[:,:,k])
    plot!(clim=(minimum(P),maximum(P)))
end
gif(anim,"uv.gif",fps=10)

PlotSolution(v[16,:,:],w[16,:,:],P[16,:,:])

plot(y,z,u[3,:,:]',aspect_ratio=:equal,xlim=(0,1),ylim=(0,1),title="Perfil de velocidades",xlabel="y[m]",ylabel="z[m]",cam=(60,30),st=:heatmap)
savefig("perfilVelocidade.png")

plot(y,z,u[end,:,:],cam=(60,30),st=:wireframe)
maximum(u[7,:,:]')
maximum(u[end,:,:]')


to
