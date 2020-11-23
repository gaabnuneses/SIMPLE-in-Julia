function MomentoX_Uns(u0,u_old,v_old,u,v,P,nIT_vel)
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
       mn = .5*ρ*dx*(v_old[i,j+1]+v_old[i-1,j+1])
       ms = .5*ρ*dx*(v_old[i,j]+v_old[i-1,j])
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


    Apu = Ae+Aw+An+As .+ρ*dx*dy/dt

    Apu = Apu/Ωu

    for nit = 1:nIT_vel
        for i = 3:nx
            for j = 2:ny
                u[i,j] = (1 - Ωu)*u_old[i,j]+1/Apu[i,j]*(Ae[i,j]*u[i+1,j]+Aw[i,j]*u[i-1,j]+An[i,j]*u[i,j+1]+As[i,j]*u[i,j-1]-β*(Dcc[i,j]-Dcu[i,j])+dy*(P[i-1,j]-P[i,j])+ ρ*dx*dy*u0[i,j]/dt)
            end
        end
    end
    return u, Apu
end

function MomentoY_Uns(v0,u_old,v_old,u,v,P,nIT_vel)
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
       mn = .5*ρ*dx*(v_old[i,j+1]+v_old[i-1,j+1])
       ms = .5*ρ*dx*(v_old[i,j]+v_old[i-1,j])
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

    Apv = Ae+Aw+An+As .+ρ*dx*dy/dt
    Apv = Apv/Ωv

    for nit = 1:nIT_vel
        for i = 2:nx
            for j = 3:ny
                v[i,j] = (1 - Ωv)*v_old[i,j]+1/Apv[i,j]*(Ae[i,j]*v[i+1,j]+Aw[i,j]*v[i-1,j]+An[i,j]*v[i,j+1]+As[i,j]*v[i,j-1]-β*(Dcc[i,j]-Dcu[i,j])+dx*(P[i,j-1]-P[i,j])+ρ*dx*dy*v0[i,j]/dt)
            end
        end
    end

    return v, Apv
end

function Pressao_Uns(u,v,P,Apu,Apv,nIT_P)
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

function ConservacaoMassa_Uns(u,v,P,Pp,Apv,Apu)
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

function Unsteady(u,v,P,nIt,nIT_vel,nIT_P)
    Uk = zeros(nx+1,ny+1,Int(round(nt/10))+2)
    Vk = zeros(nx+1,ny+1,Int(round(nt/10))+2)
    Pk = zeros(nx+1,ny+1,Int(round(nt/10))+2)
    Uk[:,:,1] = u[:,:]
    Vk[:,:,1] = v[:,:]
    Pk[:,:,1] = P[:,:]
    counter = 0
    k = 2
    erro = []
    for ItTime = 1:nt
        u0 = u[:,:]
        v0 = v[:,:]

        counter += 1
        if counter == 10
            Uk[:,:,k] = u[:,:]
            Vk[:,:,k] = v[:,:]
            Pk[:,:,k] = P[:,:]
            k+=1
            counter = 0
        end

        for OuterLoop in 1:nIt
            u_old,v_old = u[:,:],v[:,:]

            ## Momento em X
            u,Apu = MomentoX_Uns(u0,u_old,v_old,u,v,P,nIT_vel)

            ## Momento em Y
            v,Apv = MomentoY_Uns(v0,u_old,v_old,u,v,P,nIT_vel)

            ## Pressure correction
            P,Pp = Pressao_Uns(u,v,P,Apu,Apv,nIT_P)

            u,v = ConservacaoMassa_Uns(u,v,P,Pp,Apv,Apu)

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
    end
    return Uk,Vk,Pk
end
