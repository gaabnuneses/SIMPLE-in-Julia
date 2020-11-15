# SIMPLE-in-Julia
In computational fluid dynamics (CFD), the SIMPLE algorithm is a widely used numerical procedure to solve the Navier–Stokes equations. SIMPLE is an acronym for Semi-Implicit Method for Pressure Linked Equations.

![alt text](https://github.com/gaabnuneses/SIMPLE-in-Julia/blob/main/Cavidade.png?raw=true)


## Initialization

### Domain discretization
```julia
nx = 100; ny = 20
xmax = 5; ymax = 1
dx = xmax/nx; dy = ymax/ny
x = 0:dx:xmax; y = 0:dy:ymax
```
### Fluid Properties
```julia
ρ = 1
μ = 0.01
```
### Numerical Method Properties
```julia
Ωu = .1
Ωv = .1
Ωp = .1
β = 0.95
```
### u,v,p initialization
```julia
u = zeros(nx+1,ny+1)
v = zeros(nx+1,ny+1)
P = zeros(nx+1,ny+1)
```
### Boundary Conditions
```julia
function setConditions(u,v,p)
    # Flow in a Channel
    # Wall condition
    u[:,1] = zeros(length(u[:,1])) .+0 # SUL
    u[:,end] = zeros(length(u[:,end])) .+0 # Norte
    v[:,1] = zeros(length(u[:,1])) # SUL
    v[:,end] = zeros(length(u[:,end])) # Norte
    v[1,:] = zeros(length(u[1,:])) # OESTE
    v[end,:] = zeros(length(u[end,:])) # LESTE

    # Inlet Speed
    u[2,:] = zeros(length(u[1,:])) .+1 # OESTE

    # Optional Adding a block
    # u[5:10,1:10] = zeros(length(u[5:10,1:10])) .+0
    # v[5:10,1:10] = zeros(length(u[5:10,1:10])) .+0
    # v[40:60,end] = zeros(length(u[40:60,end])) .-1

    # Pressure boundary condition
    P[end-1,:] = zeros(length(u[end,:])) .+0 # LESTE
    u[end,(2:end-1)] = u[end-1,(2:end-1)] .+ dx/dy*(v[end-1,(2:end-1)] - v[end-1,(3:end)])

    return u,v,p
end

u,v,P=setConditions(u,v,P)
```

### Functions that solve the domain for X-Momentum and Y-Momentum
```julia
function MomentoX(u_old,v_old,u,v,P,nIT_vel)
    Ae = zeros(nx+1,ny+1)
    Aw = zeros(nx+1,ny+1)
    An = zeros(nx+1,ny+1)
    As = zeros(nx+1,ny+1)
    Apu = zeros(nx+1,ny+1)
    Dcu = zeros(nx+1,ny+1)
    Dcc = zeros(nx+1,ny+1)
    for i in 3:nx,j in 2:ny
       me = .5*ρ*dy*(u_old[i+1,j]+u_old[i,j])
       mw = .5*ρ*dy*(u_old[i-1,j]+u_old[i,j])
       mn = .5*ρ*dx*(v_old[i,j+1]+v_old[i,j])
       ms = .5*ρ*dx*(v_old[i,j]+v_old[i,j-1])
        Ae[i,j] = maximum([-me,0])
        Aw[i,j] = maximum([mw,0])
        An[i,j] = maximum([-mn,0])
        As[i,j] = maximum([ms,0])
        Apu[i,j] = Ae[i,j]+Aw[i,j]+An[i,j]+As[i,j]
        Dcu[i,j] = -(Ae[i,j]*u[i+1,j]+Aw[i,j]*u[i-1,j]+An[i,j]*u[i,j+1]+As[i,j]*u[i,j-1])+Apu[i,j]*u[i,j]
        Dcc[i,j] = .5*(me*(u[i+1,j]+u[i,j])-mw*(u[i,j]+u[i-1,j])+mn*(u[i,j+1]+u[i,j])-ms*(u[i,j]+u[i,j-1]) )
    end

    Ae = Ae .+ μ*dy/dx
    Aw = Aw .+ μ*dy/dx
    An = An .+ μ*dx/dy
    As = As .+ μ*dx/dy
    for i in 3:nx
       mn = .5*ρ*dx*(v_old[i,ny+1]+v_old[i-1,ny+1])
       ms = .5*ρ*dx*(v_old[i,2]+v_old[i-1,2])
        An[i,ny] = maximum([-mn,0])+ μ*dx/(dy/2)
        As[i,2] = maximum([ms,0])+ μ*dx/(dy/2)
    end


    Apu = Ae+Aw+An+As

    Apu = Apu/Ωu

    for nit = 1:nIT_vel
        for i = 3:nx
            for j = 2:ny
                u[i,j] = (1 - Ωu)*u_old[i,j]+1/Apu[i,j]*(Ae[i,j]*u[i+1,j]+Aw[i,j]*u[i-1,j]+An[i,j]*u[i,j+1]+As[i,j]*u[i,j-1]-β*(Dcc[i,j]-Dcu[i,j])+dy*(P[i-1,j]-P[i,j]))
            end
        end
    end

    return u, Apu
end

function MomentoY(u_old,v_old,u,v,P,nIT_vel)
    Ae = zeros(nx+1,ny+1)
    Aw = zeros(nx+1,ny+1)
    An = zeros(nx+1,ny+1)
    As = zeros(nx+1,ny+1)
    Apv = zeros(nx+1,ny+1)
    Dcu = zeros(nx+1,ny+1)
    Dcc = zeros(nx+1,ny+1)
    # u_old,v_old = u[:,:],v[:,:]
    for i in 2:nx,j in 3:ny
       me = .5*ρ*dy*(u_old[i+1,j]+u_old[i,j])
       mw = .5*ρ*dy*(u_old[i-1,j]+u_old[i,j])
       mn = .5*ρ*dx*(v_old[i,j+1]+v_old[i,j])
       ms = .5*ρ*dx*(v_old[i,j]+v_old[i,j-1])
        Ae[i,j] = maximum([-me,0])
        Aw[i,j] = maximum([mw,0])
        An[i,j] = maximum([-mn,0])
        As[i,j] = maximum([ms,0])
        Apv[i,j] = Ae[i,j]+Aw[i,j]+An[i,j]+As[i,j]
        Dcu[i,j] = -(Ae[i,j]*v[i+1,j]+Aw[i,j]*v[i-1,j]+An[i,j]*v[i,j+1]+As[i,j]*v[i,j-1])+Apv[i,j]*v[i,j]
        Dcc[i,j] = .5*(me*(v[i+1,j]+v[i,j])-mw*(v[i,j]+v[i-1,j])+mn*(v[i,j+1]+v[i,j])-ms*(v[i,j]+v[i,j-1]) )
    end

    Ae = Ae .+ μ*dy/dx
    Aw = Aw .+ μ*dy/dx
    An = An .+ μ*dx/dy
    As = As .+ μ*dx/dy

    for j in 3:ny
       me = .5*ρ*dy*(u_old[nx+1,j]+u_old[nx+1,j-1])
       mw = .5*ρ*dy*(u_old[2,j]+u_old[2,j-1])
        Ae[nx,j] = maximum([-me,0])+ μ*dy/(dx/2)
        Aw[2,j] = maximum([mw,0])+ μ*dy/(dx/2)
    end

    Apv = Ae+Aw+An+As
    Apv = Apv/Ωv

    for nit = 1:nIT_vel
        for i = 2:nx
            for j = 3:ny
                v[i,j] = (1 - Ωv)*v_old[i,j]+1/Apv[i,j]*(Ae[i,j]*v[i+1,j]+Aw[i,j]*v[i-1,j]+An[i,j]*v[i,j+1]+As[i,j]*v[i,j-1]-β*(Dcc[i,j]-Dcu[i,j])+dx*(P[i,j-1]-P[i,j]))
            end
        end
    end
    return v, Apv
end
```

### Pressure Correction and Solution for Mass Conservation
```julia
function Pressao(u,v,P,Apu,Apv,nIT_P)
    Ae = zeros(nx+1,ny+1)
    Aw = zeros(nx+1,ny+1)
    An = zeros(nx+1,ny+1)
    As = zeros(nx+1,ny+1)

    for i in 2:nx, j in 2:ny
        Ae[i,j] = ρ*dy^2/Apu[i+1,j]
        Aw[i,j] = ρ*dy^2/Apu[i,j]
        An[i,j] = ρ*dx^2/Apv[i,j+1]
        As[i,j] = ρ*dx^2/Apv[i,j]
    end

    Ae[nx,:]=zeros(length(Ae[nx,:]))
    Aw[2,:]=zeros(length(Aw[1,:]))
    An[:,ny]=zeros(length(An[:,ny]))
    As[:,2]=zeros(length(As[:,1]))

    App = Ae+Aw+An+As
    # App[2,2] = 1e30
    App[nx,:]= fill(1e30,length(App[nx,:]))
    Pp = zeros(nx+1,ny+1)

    Source = zeros(nx+1,ny+1)
    for i in 2:nx,j in 2:ny
        Source[i,j] = ρ*dy*(u[i+1,j]-u[i,j]) + ρ*dx*(v[i,j+1]-v[i,j])
    end

    # Total = sqrt(sum(Source.^2))
    # println("Source = $Total")
    # erro = [erro;Total]
    for nit in 1:nIT_P
        for j in 2:ny
            for i in 2:nx
                Pp[i,j] = Pp[i,j] + 1.7/App[i,j]*(Ae[i,j]*Pp[i+1,j]+Aw[i,j]*Pp[i-1,j]+An[i,j]*Pp[i,j+1]+As[i,j]*Pp[i,j-1]-Source[i,j]-Pp[i,j]*App[i,j])
            end
        end
    end

    # Cprrigindo P
    for i in 2:nx, j in 2:ny
        P[i,j]=P[i,j] + Ωp*Pp[i,j]
    end
    return P, Pp
end

function ConservacaoMassa(u,v,P,Pp,Apv,Apu)
    # COrrigindo u
    for i in 3:nx, j in 2:ny
        u[i,j]=u[i,j] + dy/Apu[i,j]*(Pp[i-1,j]-Pp[i,j])
    end
    # COrrigindo v
    for i in 2:nx, j in 3:ny
        v[i,j]=v[i,j] + dx/Apv[i,j]*(Pp[i,j-1]-Pp[i,j])
    end
    return u,v
end

```

### Mass Source (called as error) calculation
```julia
function ErroSource(u,v,erro)
    Source = zeros(nx+1,ny+1)
    for i in 2:nx,j in 2:ny
        Source[i,j] = ρ*dy*(u[i+1,j]-u[i,j]) + ρ*dx*(v[i,j+1]-v[i,j])
    end
    Total = sqrt(sum(Source.^2))

    erro = [erro;Total]
    return erro
end
```


### Solver that run previous functions iteractively 
```julia
function Solve(u,v,P,nIt,nIT_vel,nIT_P)
    erro = []
    for OuterLoop in 1:nIt
        println("Número da Iteração: $OuterLoop / $nIt")
        u_old,v_old = u[:,:],v[:,:]

        ## Momento em X
        u,Apu = MomentoX(u_old,v_old,u,v,P,nIT_vel)

        ## Momento em Y
        v,Apv = MomentoY(u_old,v_old,u,v,P,nIT_vel)

        ## Pressure correction
        P,Pp = Pressao(u,v,P,Apu,Apv,nIT_P)

        u,v = ConservacaoMassa(u,v,P,Pp,Apv,Apu)

        u,v,P=setConditions(u,v,P)

        erro = ErroSource(u,v,erro)

        ## Critérios de parada
        if isnan(erro[end])
            break
        end

        if OuterLoop>2
            if (abs(erro[end-1]-erro[end])/abs(erro[end-1]))<1e-12
                break
            end
        end
        # u[end,:] = u[end-1,:]
    end
    return u,v,P,erro
end

u,v,P,ϵ = Solve(u,v,P,500,10,100)
```
