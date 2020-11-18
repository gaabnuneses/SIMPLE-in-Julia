function MomentoX(u_old,v_old,w_old,u,v,w,P,nIT_vel)
    Ae = zeros(nx+1,ny+1,nz+1)
    Aw = zeros(nx+1,ny+1,nz+1)
    An = zeros(nx+1,ny+1,nz+1)
    As = zeros(nx+1,ny+1,nz+1)
    At = zeros(nx+1,ny+1,nz+1)
    Ab = zeros(nx+1,ny+1,nz+1)
    Apu = fill(1.0,nx+1,ny+1,nz+1)
    Dcu = zeros(nx+1,ny+1,nz+1)
    Dcc = zeros(nx+1,ny+1,nz+1)
    for i in 3:nx,j in 2:ny,k in 2:nz
       me = .5*ρ*dy*dz*(u_old[i+1,j,k]+u_old[i,j,k])
       mw = .5*ρ*dy*dz*(u_old[i-1,j,k]+u_old[i,j,k])
       mn = .5*ρ*dx*dz*(v_old[i-1,j+1,k]+v_old[i,j+1,k])
       ms = .5*ρ*dx*dz*(v_old[i-1,j,k]+v_old[i,j,k])
       mt = .5*ρ*dx*dy*(w_old[i,j,k+1]+w_old[i-1,j,k+1])
       mb = .5*ρ*dx*dy*(w_old[i,j,k]+w_old[i-1,j,k])
        Ae[i,j,k] = maximum([-me,0])
        Aw[i,j,k] = maximum([mw,0])
        An[i,j,k] = maximum([-mn,0])
        As[i,j,k] = maximum([ms,0])
        At[i,j,k] = maximum([-mt,0])
        Ab[i,j,k] = maximum([mb,0])
        Apu[i,j,k] = Ae[i,j,k]+Aw[i,j,k]+An[i,j,k]+As[i,j,k]+At[i,j,k]+Ab[i,j,k]
        Dcu[i,j,k] = -(Ae[i,j,k]*u[i+1,j,k]+Aw[i,j,k]*u[i-1,j,k]+An[i,j,k]*u[i,j+1,k]+As[i,j,k]*u[i,j-1,k]+At[i,j,k]*u[i,j,k+1]+Ab[i,j,k]*u[i,j,k-1])+Apu[i,j,k]*u[i,j,k]
        Dcc[i,j,k] = .5*(me*(u[i+1,j,k]+u[i,j,k])-mw*(u[i,j,k]+u[i-1,j,k])+mn*(u[i,j+1,k]+u[i,j,k])-ms*(u[i,j,k]+u[i,j-1,k]) +mt*(u[i,j,k+1]+u[i,j,k])-mb*(u[i,j,k]+u[i,j,k-1]) )
    end

    Ae = Ae .+ μ*dy*dz/dx
    Aw = Aw .+ μ*dy*dz/dx
    An = An .+ μ*dx*dz/dy
    As = As .+ μ*dx*dz/dy
    At = At .+ μ*dx*dy/dz
    Ab = Ab .+ μ*dx*dy/dz
    for i in 3:nx, k in 2:nz
        mn = .5*ρ*dx*dz*(v_old[i-1,ny+1,k]+v_old[i,ny+1,k])
        ms = .5*ρ*dx*dz*(v_old[i-1,2,k]+v_old[i,2,k])
        An[i,ny,k] = maximum([-mn,0])+ μ*dx*dz/(dy/2)
        As[i,2,k] = maximum([ms,0])+ μ*dx*dz/(dy/2)
    end

    for i in 3:nx, j in 2:ny
        mt = .5*ρ*dx*dy*(w_old[i-1,j,nz+1]+w_old[i,j,nz+1])
        mb = .5*ρ*dx*dy*(w_old[i-1,j,2]+w_old[i,j,2])
        At[i,j,nz] = maximum([-mt,0])+ μ*dx*dy/(dz/2)
        Ab[i,j,2] = maximum([mb,0])+ μ*dx*dy/(dz/2)
    end

    Apu = Ae+Aw+An+As+At+Ab

    Apu = Apu/Ωu

    for nit = 1:nIT_vel
        for i = 3:nx
            for j = 2:ny
                for k = 2:nz
                    Valor = Ae[i,j,k]*u[i+1,j,k] + Aw[i,j,k]*u[i-1,j,k] + An[i,j,k]*u[i,j+1,k] + As[i,j,k]*u[i,j-1,k] + At[i,j,k]*u[i,j,k+1] + Ab[i,j,k]*u[i,j,k-1] - β*(Dcc[i,j,k]-Dcu[i,j,k]) + dy*dz*(P[i-1,j,k]-P[i,j,k])
                    u[i,j,k] = (1-Ωu)*u_old[i,j,k]+Valor/Apu[i,j,k]
                end
            end
        end
    end

    return u, Apu
end

function MomentoY(u_old,v_old,w_old,u,v,w,P,nIT_vel)
    Ae = zeros(nx+1,ny+1,nz+1)
    Aw = zeros(nx+1,ny+1,nz+1)
    An = zeros(nx+1,ny+1,nz+1)
    As = zeros(nx+1,ny+1,nz+1)
    At = zeros(nx+1,ny+1,nz+1)
    Ab = zeros(nx+1,ny+1,nz+1)
    Apv = fill(1.0,nx+1,ny+1,nz+1)
    Dcu = zeros(nx+1,ny+1,nz+1)
    Dcc = zeros(nx+1,ny+1,nz+1)
    # u_old,v_old = u[:,:],v[:,:]
    for i in 2:nx,j in 3:ny,k in 2:nz
        me = .5*ρ*dy*dz*(u_old[i+1,j,k]+u_old[i+1,j-1,k])
        mw = .5*ρ*dy*dz*(u_old[i,j-1,k]+u_old[i,j,k])
        mn = .5*ρ*dx*dz*(v_old[i,j+1,k]+v_old[i,j,k])
        ms = .5*ρ*dx*dz*(v_old[i,j,k]+v_old[i,j-1,k])
        mt = .5*ρ*dx*dy*(w_old[i,j,k+1]+w_old[i,j-1,k+1])
        mb = .5*ρ*dx*dy*(w_old[i,j,k]+w_old[i,j-1,k])
         Ae[i,j,k] = maximum([-me,0])
         Aw[i,j,k] = maximum([mw,0])
         An[i,j,k] = maximum([-mn,0])
         As[i,j,k] = maximum([ms,0])
         At[i,j,k] = maximum([-mt,0])
         Ab[i,j,k] = maximum([mb,0])
         Apv[i,j,k] = Ae[i,j,k]+Aw[i,j,k]+An[i,j,k]+As[i,j,k]+At[i,j,k]+Ab[i,j,k]
         Dcu[i,j,k] = -(Ae[i,j,k]*v[i+1,j,k]+Aw[i,j,k]*v[i-1,j,k]+An[i,j,k]*v[i,j+1,k]+As[i,j,k]*v[i,j-1,k]+At[i,j,k]*v[i,j,k+1]+Ab[i,j,k]*v[i,j,k-1])+Apv[i,j,k]*v[i,j,k]
         Dcc[i,j,k] = .5*(me*(v[i+1,j,k]+v[i,j,k])-mw*(v[i,j,k]+v[i-1,j,k])+mn*(v[i,j+1,k]+v[i,j,k])-ms*(v[i,j,k]+v[i,j-1,k]) +mt*(v[i,j,k+1]+v[i,j,k])-mb*(v[i,j,k]+v[i,j,k-1]) )
    end

    Ae = Ae .+ μ*dy*dz/dx
    Aw = Aw .+ μ*dy*dz/dx
    An = An .+ μ*dx*dz/dy
    As = As .+ μ*dx*dz/dy
    At = At .+ μ*dx*dy/dz
    Ab = Ab .+ μ*dx*dy/dz

    for j in 3:ny, k in 2:nz
        me = .5*ρ*dy*dz*(u_old[nx+1,j,k]+u_old[nx+1,j-1,k])
        mw = .5*ρ*dy*dz*(u_old[2,j-1,k]+u_old[2,j,k])
        Ae[nx,j,k] = maximum([-me,0])+ μ*dy*dz/(dx/2)
        Aw[2,j,k] = maximum([mw,0])+ μ*dy*dz/(dx/2)
    end

    for i in 2:nx, j in 3:ny
        mt = .5*ρ*dx*dy*(w_old[i,j-1,nz+1]+w_old[i,j,nz+1])
        mb = .5*ρ*dx*dy*(w_old[i,j-1,2]+w_old[i,j,2])
        At[i,j,nz] = maximum([-mt,0])+ μ*dx*dy/(dz/2)
        Ab[i,j,2] = maximum([mb,0])+ μ*dx*dy/(dz/2)
    end

    Apv = Ae+Aw+An+As+At+Ab

    Apv = Apv/Ωv

    for nit = 1:nIT_vel
        for i = 2:nx
            for j = 3:ny
                for k = 2:nz
                    Valor = Ae[i,j,k]*v[i+1,j,k] + Aw[i,j,k]*v[i-1,j,k] + An[i,j,k]*v[i,j+1,k] + As[i,j,k]*v[i,j-1,k] + At[i,j,k]*v[i,j,k+1] + Ab[i,j,k]*v[i,j,k-1] - β*(Dcc[i,j,k]-Dcu[i,j,k]) + dx*dz*(P[i,j-1,k]-P[i,j,k])
                    v[i,j,k] = (1-Ωv)*v_old[i,j,k]+Valor/Apv[i,j,k]
                end
            end
        end
    end
    return v, Apv
end


function MomentoZ(u_old,v_old,w_old,u,v,w,P,nIT_vel)
    Ae = zeros(nx+1,ny+1,nz+1)
    Aw = zeros(nx+1,ny+1,nz+1)
    An = zeros(nx+1,ny+1,nz+1)
    As = zeros(nx+1,ny+1,nz+1)
    At = zeros(nx+1,ny+1,nz+1)
    Ab = zeros(nx+1,ny+1,nz+1)
    Apw = fill(1.0,nx+1,ny+1,nz+1)
    Dcu = zeros(nx+1,ny+1,nz+1)
    Dcc = zeros(nx+1,ny+1,nz+1)
    for i in 2:nx,j in 2:ny,k in 3:nz
        me = .5*ρ*dy*dz*(u_old[i+1,j,k]+u_old[i+1,j,k-1])
        mw = .5*ρ*dy*dz*(u_old[i,j,k-1]+u_old[i,j,k])
        mn = .5*ρ*dx*dz*(v_old[i,j+1,k]+v_old[i,j+1,k-1])
        ms = .5*ρ*dx*dz*(v_old[i,j,k]+v_old[i,j,k-1])
        mt = .5*ρ*dx*dy*(w_old[i,j,k+1]+w_old[i,j,k])
        mb = .5*ρ*dx*dy*(w_old[i,j,k]+w_old[i,j,k-1])
         Ae[i,j,k] = maximum([-me,0])
         Aw[i,j,k] = maximum([mw,0])
         An[i,j,k] = maximum([-mn,0])
         As[i,j,k] = maximum([ms,0])
         At[i,j,k] = maximum([-mt,0])
         Ab[i,j,k] = maximum([mb,0])
         Apw[i,j,k] = Ae[i,j,k]+Aw[i,j,k]+An[i,j,k]+As[i,j,k]+At[i,j,k]+Ab[i,j,k]
         Dcu[i,j,k] = -(Ae[i,j,k]*w[i+1,j,k]+Aw[i,j,k]*w[i-1,j,k]+An[i,j,k]*w[i,j+1,k]+As[i,j,k]*w[i,j-1,k]+At[i,j,k]*w[i,j,k+1]+Ab[i,j,k]*w[i,j,k-1])+Apw[i,j,k]*w[i,j,k]
         Dcc[i,j,k] = .5*(me*(w[i+1,j,k]+w[i,j,k])-mw*(w[i,j,k]+w[i-1,j,k])+mn*(w[i,j+1,k]+w[i,j,k])-ms*(w[i,j,k]+w[i,j-1,k]) +mt*(w[i,j,k+1]+w[i,j,k])-mb*(w[i,j,k]+w[i,j,k-1]) )
    end

    Ae = Ae .+ μ*dy*dz/dx
    Aw = Aw .+ μ*dy*dz/dx
    An = An .+ μ*dx*dz/dy
    As = As .+ μ*dx*dz/dy
    At = At .+ μ*dx*dy/dz
    Ab = Ab .+ μ*dx*dy/dz

    for i in 2:nx, k in 3:nz
        mn = .5*ρ*dx*dz*(v_old[i,ny+1,k]+v_old[i,ny+1,k-1])
        ms = .5*ρ*dx*dz*(v_old[i,2,k]+v_old[i,2,k-1])
        An[i,ny,k] = maximum([-mn,0])+ μ*dx*dz/(dy/2)
        As[i,2,k] = maximum([ms,0])+ μ*dx*dz/(dy/2)
    end

    for j in 2:ny, k in 3:nz
        me = .5*ρ*dy*dz*(u_old[nx+1,j,k]+u_old[nx+1,j,k-1])
        mw = .5*ρ*dy*dz*(u_old[2,j,k]+u_old[2,j,k-1])
        Ae[nx,j,k] = maximum([-me,0])+ μ*dy*dz/(dx/2)
        Aw[2,j,k] = maximum([mw,0])+ μ*dy*dz/(dx/2)
    end

    Apw = Ae+Aw+An+As+At+Ab

    Apw = Apw/Ωw

    for nit = 1:nIT_vel
        for i = 2:nx
            for j = 2:ny
                for k = 3:nz
                    Valor = Ae[i,j,k]*w[i+1,j,k] + Aw[i,j,k]*w[i-1,j,k] + An[i,j,k]*w[i,j+1,k] + As[i,j,k]*w[i,j-1,k] + At[i,j,k]*w[i,j,k+1] + Ab[i,j,k]*w[i,j,k-1] - β*(Dcc[i,j,k]-Dcu[i,j,k]) + dx*dy*(P[i,j,k-1]-P[i,j,k])
                    w[i,j,k] = (1-Ωw)*w_old[i,j,k]+Valor/Apw[i,j,k]
                end
            end
        end
    end
    return w, Apw
end


function Pressao(u,v,w,P,Apu,Apv,Apw,nIT_P)
    Ae = zeros(nx+1,ny+1,nz+1)
    Aw = zeros(nx+1,ny+1,nz+1)
    An = zeros(nx+1,ny+1,nz+1)
    As = zeros(nx+1,ny+1,nz+1)
    At = zeros(nx+1,ny+1,nz+1)
    Ab = zeros(nx+1,ny+1,nz+1)
    for i in 2:nx, j in 2:ny, k in 2:nz
        Ae[i,j,k] = ρ*(dz*dy)^2/Apu[i+1,j,k]
        Aw[i,j,k] = ρ*(dz*dy)^2/Apu[i,j,k]
        An[i,j,k] = ρ*(dz*dx)^2/Apv[i,j+1,k]
        As[i,j,k] = ρ*(dz*dx)^2/Apv[i,j,k]
        At[i,j,k] = ρ*(dx*dy)^2/Apw[i,j,k+1]
        Ab[i,j,k] = ρ*(dx*dy)^2/Apw[i,j,k]
    end

    Ae[nx,:,:]=zeros(size(Ae[nx,:,:]))
    Aw[2,:,:]=zeros(size(Aw[1,:,:]))
    An[:,ny,:]=zeros(size(An[:,ny,:]))
    As[:,2,:]=zeros(size(As[:,1,:]))
    At[:,:,nz]=zeros(size(At[:,:,nz]))
    Ab[:,:,2]=zeros(size(Ab[:,:,2]))

    App = Ae+Aw+An+As+At+Ab
    # App[2,2,2] = 1e30
    App[nx,:,:]= fill(1.0e30,size(App[nx,:,:]))
    Pp = zeros(nx+1,ny+1,nz+1)
    Source = zeros(nx+1,ny+1,nz+1)
    for i in 2:nx,j in 2:ny,k in 2:nz
        Source[i,j,k] = ρ*dy*dz*(u[i+1,j,k]-u[i,j,k]) + ρ*dx*dz*(v[i,j+1,k]-v[i,j,k]) + ρ*dx*dy*(w[i,j,k+1]-w[i,j,k])
    end

    for nit in 1:nIT_P
        for k in 2:nz
            for j in 2:ny
                for i in 2:nx
                    Pp[i,j,k] = Pp[i,j,k] + Ωpp/App[i,j,k]*(Ae[i,j,k]*Pp[i+1,j,k]+Aw[i,j,k]*Pp[i-1,j,k]+An[i,j,k]*Pp[i,j+1,k]+As[i,j,k]*Pp[i,j-1,k]+At[i,j,k]*Pp[i,j,k+1]+Ab[i,j,k]*Pp[i,j,k-1]-Source[i,j,k]-Pp[i,j,k]*App[i,j,k])
                end
            end
        end
    end

    # Corrigindo P
    for i in 2:nx, j in 2:ny, k in 2:nz
        P[i,j,k]=P[i,j,k] + Ωp*Pp[i,j,k]
    end
    return P, Pp
end

function ConservacaoMassa(u,v,w,P,Pp,Apv,Apu,Apw)
    # Corrigindo u
    for i in 3:nx, j in 2:ny, k in 2:nz
        u[i,j,k]=u[i,j,k] + dy*dz/Apu[i,j,k]*(Pp[i-1,j,k]-Pp[i,j,k])
    end
    # Corrigindo v
    for i in 2:nx, j in 3:ny, k in 2:nz
        v[i,j,k]=v[i,j,k] + dx*dz/Apv[i,j,k]*(Pp[i,j-1,k]-Pp[i,j,k])
    end
    # Corrigindo w
    for i in 2:nx, j in 2:ny, k in 3:nz
        w[i,j,k]=w[i,j,k] + dx*dy/Apw[i,j,k]*(Pp[i,j,k-1]-Pp[i,j,k])
    end
    return u,v,w
end

function ErroSource(u,v,w,erro)
    Source = zeros(nx+1,ny+1,nz+1)
    for i in 2:nx,j in 2:ny,k in 2:ny
        Source[i,j,k] = ρ*dy*dz*(u[i+1,j,k]-u[i,j,k]) + ρ*dx*dz*(v[i,j+1,k]-v[i,j,k])  + ρ*dx*dz*(w[i,j,k+1]-w[i,j,k])
    end
    Total = sqrt(sum(Source.^2))

    erro = [erro;Total]
    return erro
end

function Solve(u,v,w,P,nIt,nIT_vel,nIT_P)
    erro = []
    for OuterLoop in 1:nIt
        println("Número da Iteração: $OuterLoop / $nIt")
        u_old,v_old,w_old = u[:,:,:],v[:,:,:],w[:,:,:]

        ## Momento em X
        u,Apu = MomentoX(u_old,v_old,w_old,u,v,w,P,nIT_vel)

        ## Momento em Y
        v,Apv = MomentoY(u_old,v_old,w_old,u,v,w,P,nIT_vel)

        ## Momento em Z
        w,Apw = MomentoZ(u_old,v_old,w_old,u,v,w,P,nIT_vel)

        ## Pressure correction
        P,Pp = Pressao(u,v,w,P,Apu,Apv,Apw,nIT_P)

        u,v,w = ConservacaoMassa(u,v,w,P,Pp,Apv,Apu,Apw)

        # u,v,w,P=setConditions(u,v,w,P)

        for j in 2:ny, k in 2:nz
            u[nx+1,j,k] = u[nx,j,k] .+ dx*dz/dy*(v[nx,j,k] - v[nx,j+1,k]) .+ dx*dy/dz*(w[nx,j,k] - w[nx,j,k+1])
        end
        # u[nx+1,:,:] = fill(1.0,size(u[nx+1,:,:]))

        erro = ErroSource(u,v,w,erro)
        ## Critérios de parada
        if isnan(erro[end])
            break
        end

        if OuterLoop>2
            if (abs(erro[end-1]-erro[end])/abs(erro[end-1]))<1e-12
                break
            end
        end

        for i in 1:nx+1
            for j in 1:ny+1
                for k in 1:nz+1
                    if x[i]>=1&&x[i]<=2
                        if y[j]>=0&&y[j]<=1
                            if z[k]>=0&&z[k]<=.5
                                u[i,j,k] = 0
                                v[i,j,k] = 0
                                w[i,j,k] = 0
                            end
                        end
                    end
                end
            end
        end

        if erro[end]<1e-5
            break
        end
        # u[end,:] = u[end-1,:]
    end
    return u,v,w,P,erro
end
