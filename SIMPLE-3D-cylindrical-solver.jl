# x,i → r
# y,j → θ
# z,k → z
# Na realidade, momento em R
function MomentoX(u_old,v_old,w_old,u,v,w,P,nIT_vel)
    Ae = zeros(nx+1,ny+1,nz+1)
    Aw = zeros(nx+1,ny+1,nz+1)
    An = zeros(nx+1,ny+1,nz+1)
    As = zeros(nx+1,ny+1,nz+1)
    At = zeros(nx+1,ny+1,nz+1)
    Ab = zeros(nx+1,ny+1,nz+1)
    Apu = fill(0.0,nx+1,ny+1,nz+1)
    Dcu = zeros(nx+1,ny+1,nz+1)
    Dcc = zeros(nx+1,ny+1,nz+1)
    for i in 3:nx,j in 2:ny,k in 2:nz
        # areas
        ae = (x[i]+x[i+1])/2 * (y[j+1]-y[j])*(z[k]+z[k+1])/2
        aw = (x[i]+x[i-1])/2 * (y[j+1]-y[j])*(z[k]+z[k+1])/2
        an = (x[i+1]-x[i+1])/2*(z[k]+z[k+1])/2
        as = (x[i+1]-x[i+1])/2*(z[k]+z[k+1])/2
        at = (y[j+1]-y[j])*(((x[i]+x[i+1])/2)^2-((x[i]+x[i-1])/2)^2)/2
        ab = (y[j+1]-y[j])*(((x[i]+x[i+1])/2)^2-((x[i]+x[i-1])/2)^2)/2
        # vazoes
        me = .5*ρ*ae*(u_old[i+1,j,k]+u_old[i,j,k])
        mw = .5*ρ*aw*(u_old[i-1,j,k]+u_old[i,j,k])
        mn = .5*ρ*an*(v_old[i-1,j+1,k]+v_old[i,j+1,k])
        ms = .5*ρ*as*(v_old[i-1,j,k]+v_old[i,j,k])
        mt = .5*ρ*at*(w_old[i,j,k+1]+w_old[i-1,j,k+1])
        mb = .5*ρ*ab*(w_old[i,j,k]+w_old[i-1,j,k])
        # coeficientes sistema
        Ae[i,j,k] = maximum([-me,0])
        Aw[i,j,k] = maximum([mw,0])
        An[i,j,k] = maximum([-mn,0])
        As[i,j,k] = maximum([ms,0])
        At[i,j,k] = maximum([-mt,0])
        Ab[i,j,k] = maximum([mb,0])
        Apu[i,j,k] = Ae[i,j,k]+Aw[i,j,k]+An[i,j,k]+As[i,j,k]+At[i,j,k]+Ab[i,j,k]
        Dcu[i,j,k] = -(Ae[i,j,k]*u[i+1,j,k]+Aw[i,j,k]*u[i-1,j,k]+An[i,j,k]*u[i,j+1,k]+As[i,j,k]*u[i,j-1,k]+At[i,j,k]*u[i,j,k+1]+Ab[i,j,k]*u[i,j,k-1])+Apu[i,j,k]*u[i,j,k]
        Dcc[i,j,k] = .5*(me*(u[i+1,j,k]+u[i,j,k])-mw*(u[i,j,k]+u[i-1,j,k])+mn*(u[i,j+1,k]+u[i,j,k])-ms*(u[i,j,k]+u[i,j-1,k]) +mt*(u[i,j,k+1]+u[i,j,k])-mb*(u[i,j,k]+u[i,j,k-1]) )
        Ae[i,j,k] += μ*ae/(x[i+1]-x[i])
        Aw[i,j,k] += μ*aw/(x[i]-x[i-1])
        An[i,j,k] += μ*an/((y[j+1]-y[j])*x[i])
        As[i,j,k] += μ*as/((y[j]-y[j-1])*x[i])
        At[i,j,k] += μ*at/(z[k+1]-z[k])
        Ab[i,j,k] += μ*ab/(z[k]-z[k-1])
    end

    # Contorno
    for i in 3:nx, k in 2:nz
        an = (x[i+1]-x[i+1])/2*(z[k]+z[k+1])/2
        as = (x[i+1]-x[i+1])/2*(z[k]+z[k+1])/2
        mn = .5*ρ*an*(v_old[i-1,ny+1,k]+v_old[i,ny+1,k])
        ms = .5*ρ*as*(v_old[i-1,2,k]+v_old[i,2,k])
        An[i,ny,k] = maximum([-mn,0])+ μ*an/((y[ny+1]-y[ny])*x[i]/2)
        As[i,2,k] = maximum([ms,0])+ μ*as/((y[2]-y[2-1])*x[i]/2)
    end

    for i in 3:nx, j in 2:ny
        at = (y[j+1]-y[j])*(((x[i]+x[i+1])/2)^2-((x[i]+x[i-1])/2)^2)/2
        ab = (y[j+1]-y[j])*(((x[i]+x[i+1])/2)^2-((x[i]+x[i-1])/2)^2)/2
        mt = .5*ρ*dx*dy*(w_old[i,j-1,nz+1]+w_old[i,j,nz+1])
        mb = .5*ρ*dx*dy*(w_old[i,j-1,2]+w_old[i,j,2])
        At[i,j,nz] = maximum([-mt,0])+ μ*at/((z[k+1]-z[k])/2)
        Ab[i,j,2] = maximum([mb,0])+ μ*ab/((z[k]-z[k-1])/2)
    end

    # Cálculo
    Apu = Ae+Aw+An+As+At+Ab

    Apu = Apu/Ωu

    for nit = 1:nIT_vel
        for i = 3:nx
            for j = 2:ny
                for k = 2:nz
                    dV = -(y[j-1]-y[j+1])*(z[k-1]-z[k+1])*(2*x[i]+x[i-1]+x[i+1])*(x[i-1]-x[i+1])/32
                    dX = (x[i+1]-x[i-1])/2
                    Valor = Ae[i,j,k]*u[i+1,j,k] + Aw[i,j,k]*u[i-1,j,k] + An[i,j,k]*u[i,j+1,k] + As[i,j,k]*u[i,j-1,k] + At[i,j,k]*u[i,j,k+1] + Ab[i,j,k]*u[i,j,k-1] - β*(Dcc[i,j,k]-Dcu[i,j,k]) + dV/dX*(P[i-1,j,k]-P[i,j,k])
                    u[i,j,k] = (1-Ωu)*u_old[i,j,k]+Valor/Apu[i,j,k]
                end
            end
        end
    end

    return u, Apu
end

# Na realidade, momento em Θ
function MomentoY(u_old,v_old,w_old,u,v,w,P,nIT_vel)
    Ae = zeros(nx+1,ny+1,nz+1)
    Aw = zeros(nx+1,ny+1,nz+1)
    An = zeros(nx+1,ny+1,nz+1)
    As = zeros(nx+1,ny+1,nz+1)
    At = zeros(nx+1,ny+1,nz+1)
    Ab = zeros(nx+1,ny+1,nz+1)
    Apv = fill(0.0,nx+1,ny+1,nz+1)
    Dcu = zeros(nx+1,ny+1,nz+1)
    Dcc = zeros(nx+1,ny+1,nz+1)
    for i in 2:nx,j in 3:ny,k in 2:nz
        # areas
        ae = x[i+1]*(y[j+1]-y[j-1])/2*(z[k]+z[k+1])/2
        aw = x[i]*(y[j+1]-y[j-1])/2*(z[k]+z[k+1])/2
        an = (x[i+1]-x[i])*(z[k]+z[k+1])/2
        as = (x[i+1]-x[i])*(z[k]+z[k+1])/2
        at = (y[j+1]-y[j-1])/2*((x[i+1])^2-(x[i-1])^2)/2
        ab = (y[j+1]-y[j-1])/2*((x[i+1])^2-(x[i-1])^2)/2
        # vazoes
        me = .5*ρ*ae*(u_old[i+1,j,k]+u_old[i+1,j-1,k])
        mw = .5*ρ*aw*(u_old[i,j-1,k]+u_old[i,j,k])
        mn = .5*ρ*an*(v_old[i,j+1,k]+v_old[i,j,k])
        ms = .5*ρ*as*(v_old[i,j,k]+v_old[i,j-1,k])
        mt = .5*ρ*at*(w_old[i,j,k+1]+w_old[i,j-1,k+1])
        mb = .5*ρ*ab*(w_old[i,j,k]+w_old[i,j-1,k])
        # coeficientes
        Ae[i,j,k] = maximum([-me,0])
        Aw[i,j,k] = maximum([mw,0])
        An[i,j,k] = maximum([-mn,0])
        As[i,j,k] = maximum([ms,0])
        At[i,j,k] = maximum([-mt,0])
        Ab[i,j,k] = maximum([mb,0])
        Apv[i,j,k] = Ae[i,j,k]+Aw[i,j,k]+An[i,j,k]+As[i,j,k]+At[i,j,k]+Ab[i,j,k]
        Dcu[i,j,k] = -(Ae[i,j,k]*v[i+1,j,k]+Aw[i,j,k]*v[i-1,j,k]+An[i,j,k]*v[i,j+1,k]+As[i,j,k]*v[i,j-1,k]+At[i,j,k]*v[i,j,k+1]+Ab[i,j,k]*v[i,j,k-1])+Apv[i,j,k]*v[i,j,k]
        Dcc[i,j,k] = .5*(me*(v[i+1,j,k]+v[i,j,k])-mw*(v[i,j,k]+v[i-1,j,k])+mn*(v[i,j+1,k]+v[i,j,k])-ms*(v[i,j,k]+v[i,j-1,k]) +mt*(v[i,j,k+1]+v[i,j,k])-mb*(v[i,j,k]+v[i,j,k-1]) )
        Ae[i,j,k] += μ*ae/(x[i+1]-x[i])
        Aw[i,j,k] += μ*aw/(x[i]-x[i-1])
        An[i,j,k] += μ*an/((y[j+1]-y[j])*x[i])
        As[i,j,k] += μ*as/((y[j]-y[j-1])*x[i])
        At[i,j,k] += μ*at/(z[k+1]-z[k])
        Ab[i,j,k] += μ*ab/(z[k]-z[k-1])
    end

    for j in 3:ny, k in 2:nz
        ae = x[nx+1]*(y[j+1]-y[j-1])/2*(z[k]+z[k+1])/2
        aw = x[2]*(y[j+1]-y[j-1])/2*(z[k]+z[k+1])/2
        me = .5*ρ*ae*(u_old[nx+1,j,k]+u_old[nx+1,j-1,k])
        mw = .5*ρ*aw*(u_old[2,j-1,k]+u_old[2,j,k])
        Ae[nx,j,k] = maximum([-me,0])+ μ*ae/((x[nx+1]-x[nx])/2)
        Aw[2,j,k] = maximum([mw,0])+ μ*aw/((x[2]-x[2-1])/2)
    end

    for i in 2:nx, j in 3:ny
        at = (y[j+1]-y[j-1])/2*((x[i+1])^2-(x[i-1])^2)/2
        ab = (y[j+1]-y[j-1])/2*((x[i+1])^2-(x[i-1])^2)/2
        mt = .5*ρ*at*(w_old[i,j-1,nz+1]+w_old[i,j,nz+1])
        mb = .5*ρ*ab*(w_old[i,j-1,2]+w_old[i,j,2])
        At[i,j,nz] = maximum([-mt,0])+ μ*at/((z[k+1]-z[k])/2)
        Ab[i,j,2] = maximum([mb,0])+ μ*ab/((z[k]-z[k-1])/2)
    end

    Apv = Ae+Aw+An+As+At+Ab

    Apv = Apv/Ωv

    for nit = 1:nIT_vel
        for i = 2:nx
            for j = 3:ny
                for k = 2:nz
                    dV = -(y[j-1]-y[j+1])*(z[k-1]-z[k+1])*(2*x[i]+x[i-1]+x[i+1])*(x[i-1]-x[i+1])/32
                    dY = (y[j+1]-y[j-1])/2
                    Valor = Ae[i,j,k]*v[i+1,j,k] + Aw[i,j,k]*v[i-1,j,k] + An[i,j,k]*v[i,j+1,k] + As[i,j,k]*v[i,j-1,k] + At[i,j,k]*v[i,j,k+1] + Ab[i,j,k]*v[i,j,k-1] - β*(Dcc[i,j,k]-Dcu[i,j,k]) + dV/dY*(P[i,j-1,k]-P[i,j,k])
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
    Apw = fill(0.0,nx+1,ny+1,nz+1)
    Dcu = zeros(nx+1,ny+1,nz+1)
    Dcc = zeros(nx+1,ny+1,nz+1)
    for i in 2:nx,j in 2:ny,k in 3:nz
        # areas
        ae = x[i+1] * (y[j+1]-y[j])*(z[k+1]-z[k-1])/2
        aw = x[i] * (y[j+1]-y[j])*(z[k+1]-z[k-1])/2
        an = (x[i+1]-x[i])*(z[k+1]-z[k-1])/2
        as = (x[i+1]-x[i])*(z[k+1]-z[k-1])/2
        at = (y[j+1]-y[j])*((x[i+1])^2-(x[i])^2)/2
        ab = (y[j+1]-y[j])*((x[i+1])^2-(x[i])^2)/2
        #vazoes
        me = .5*ρ*ae*(u_old[i+1,j,k]+u_old[i+1,j,k-1])
        mw = .5*ρ*aw*(u_old[i,j,k-1]+u_old[i,j,k])
        mn = .5*ρ*an*(v_old[i,j+1,k]+v_old[i,j+1,k-1])
        ms = .5*ρ*as*(v_old[i,j,k]+v_old[i,j,k-1])
        mt = .5*ρ*at*(w_old[i,j,k+1]+w_old[i,j,k])
        mb = .5*ρ*ab*(w_old[i,j,k]+w_old[i,j,k-1])
        Ae[i,j,k] = maximum([-me,0])
        Aw[i,j,k] = maximum([mw,0])
        An[i,j,k] = maximum([-mn,0])
        As[i,j,k] = maximum([ms,0])
        At[i,j,k] = maximum([-mt,0])
        Ab[i,j,k] = maximum([mb,0])
        Apw[i,j,k] = Ae[i,j,k]+Aw[i,j,k]+An[i,j,k]+As[i,j,k]+At[i,j,k]+Ab[i,j,k]
        Dcu[i,j,k] = -(Ae[i,j,k]*w[i+1,j,k]+Aw[i,j,k]*w[i-1,j,k]+An[i,j,k]*w[i,j+1,k]+As[i,j,k]*w[i,j-1,k]+At[i,j,k]*w[i,j,k+1]+Ab[i,j,k]*w[i,j,k-1])+Apw[i,j,k]*w[i,j,k]
        Dcc[i,j,k] = .5*(me*(w[i+1,j,k]+w[i,j,k])-mw*(w[i,j,k]+w[i-1,j,k])+mn*(w[i,j+1,k]+w[i,j,k])-ms*(w[i,j,k]+w[i,j-1,k]) +mt*(w[i,j,k+1]+w[i,j,k])-mb*(w[i,j,k]+w[i,j,k-1]) )
        Ae[i,j,k] += μ*ae/(x[i+1]-x[i])
        Aw[i,j,k] += μ*aw/(x[i]-x[i-1])
        An[i,j,k] += μ*an/((y[j+1]-y[j])*x[i])
        As[i,j,k] += μ*as/((y[j]-y[j-1])*x[i])
        At[i,j,k] += μ*at/(z[k+1]-z[k])
        Ab[i,j,k] += μ*ab/(z[k]-z[k-1])
    end

    for i in 2:nx, k in 3:nz
        an = (x[i+1]-x[i])*(z[k+1]-z[k-1])/2
        as = (x[i+1]-x[i])*(z[k+1]-z[k-1])/2
        mn = .5*ρ*an*(v_old[i,ny+1,k]+v_old[i,ny+1,k-1])
        ms = .5*ρ*as*(v_old[i,2,k]+v_old[i,2,k-1])
        An[i,ny,k] = maximum([-mn,0])+ μ*an/((y[ny+1]-y[ny])*x[i]/2)
        As[i,2,k] = maximum([ms,0])+ μ*as/((y[2]-y[2-1])*x[i]/2)
    end

    for j in 2:ny, k in 3:nz
        ae = x[nx+1]/2 * (y[j+1]-y[j])*(z[k+1]-z[k-1])/2
        aw = x[2]/2 * (y[j+1]-y[j])*(z[k+1]-z[k-1])/2
        me = .5*ρ*ae*(u_old[nx+1,j,k]+u_old[nx+1,j,k-1])
        mw = .5*ρ*aw*(u_old[2,j,k]+u_old[2,j,k-1])
        Ae[nx,j,k] = maximum([-me,0])+μ*ae/(x[nx+1]-x[nx]/2)
        Aw[2,j,k] = maximum([mw,0])+μ*aw/(x[2]-x[2-1]/2)
    end

    Apw = Ae+Aw+An+As+At+Ab

    Apw = Apw/Ωw

    for nit = 1:nIT_vel
        for i = 2:nx
            for j = 2:ny
                for k = 3:nz
                    dV = -(y[j-1]-y[j+1])*(z[k-1]-z[k+1])*(2*x[i]+x[i-1]+x[i+1])*(x[i-1]-x[i+1])/32
                    dZ = (z[k+1]-z[k-1])/2
                    Valor = Ae[i,j,k]*w[i+1,j,k] + Aw[i,j,k]*w[i-1,j,k] + An[i,j,k]*w[i,j+1,k] + As[i,j,k]*w[i,j-1,k] + At[i,j,k]*w[i,j,k+1] + Ab[i,j,k]*w[i,j,k-1] - β*(Dcc[i,j,k]-Dcu[i,j,k]) + dV/dZ*(P[i,j,k-1]-P[i,j,k])
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
        ae = (x[i+1]+x[i])/2 * (y[j+1]-y[j-1])/2*(z[k+1]-z[k-1])/2
        aw = (x[i-1]+x[i])/2 * (y[j+1]-y[j-1])/2*(z[k+1]-z[k-1])/2
        an = (x[i+1]-x[i-1])/2*(z[k+1]-z[k-1])/2
        as = (x[i+1]-x[i-1])/2*(z[k+1]-z[k-1])/2
        at = (y[j+1]-y[j-1])/2*(((x[i+1]-x[i])/2)^2-((x[i]-x[i-1])/2)^2)/2
        ab = (y[j+1]-y[j-1])/2*(((x[i+1]-x[i])/2)^2-((x[i]-x[i-1])/2)^2)/2
        Ae[i,j,k] = ρ*(ae)^2/Apu[i+1,j,k]
        Aw[i,j,k] = ρ*(aw)^2/Apu[i,j,k]
        An[i,j,k] = ρ*(an)^2/Apv[i,j+1,k]
        As[i,j,k] = ρ*(as)^2/Apv[i,j,k]
        At[i,j,k] = ρ*(at)^2/Apw[i,j,k+1]
        Ab[i,j,k] = ρ*(ab)^2/Apw[i,j,k]
    end
    #
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
        ae = (x[i+1]+x[i])/2 * (y[j+1]-y[j-1])/2*(z[k+1]-z[k-1])/2
        aw = (x[i-1]+x[i])/2 * (y[j+1]-y[j-1])/2*(z[k+1]-z[k-1])/2
        an = (x[i+1]-x[i-1])/2*(z[k+1]-z[k-1])/2
        as = (x[i+1]-x[i-1])/2*(z[k+1]-z[k-1])/2
        at = (y[j+1]-y[j-1])/2*(((x[i+1]-x[i])/2)^2-((x[i]-x[i-1])/2)^2)/2
        ab = (y[j+1]-y[j-1])/2*(((x[i+1]-x[i])/2)^2-((x[i]-x[i-1])/2)^2)/2
        Source[i,j,k] = ρ*(ae*u[i+1,j,k]-aw*u[i,j,k]) + ρ*(an*v[i,j+1,k]-as*v[i,j,k]) + ρ*(at*w[i,j,k+1]-ab*w[i,j,k])
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
        ae = (x[i+1]+x[i])/2 * (y[j+1]-y[j-1])/2*(z[k+1]-z[k-1])/2
        aw = (x[i-1]+x[i])/2 * (y[j+1]-y[j-1])/2*(z[k+1]-z[k-1])/2
        u[i,j,k]=u[i,j,k] + 1/Apu[i,j,k]*(aw*Pp[i-1,j,k]-ae*Pp[i,j,k])
    end
    # Corrigindo v
    for i in 2:nx, j in 3:ny, k in 2:nz
        an = (x[i+1]-x[i-1])/2*(z[k+1]-z[k-1])/2
        as = (x[i+1]-x[i-1])/2*(z[k+1]-z[k-1])/2
        v[i,j,k]=v[i,j,k] + 1/Apv[i,j,k]*(as*Pp[i,j-1,k]-an*Pp[i,j,k])
    end
    # Corrigindo w
    for i in 2:nx, j in 2:ny, k in 3:nz
        at = (y[j+1]-y[j-1])/2*(((x[i+1]-x[i])/2)^2-((x[i]-x[i-1])/2)^2)/2
        ab = (y[j+1]-y[j-1])/2*(((x[i+1]-x[i])/2)^2-((x[i]-x[i-1])/2)^2)/2
        w[i,j,k]=w[i,j,k] + 1/Apw[i,j,k]*(ab*Pp[i,j,k-1]-at*Pp[i,j,k])
    end
    return u,v,w
end

function ErroSource(u,v,w,erro)
    Source = zeros(nx+1,ny+1,nz+1)
    for i in 2:nx,j in 2:ny,k in 2:ny
        ae = (x[i+1]+x[i])/2 * (y[j+1]-y[j-1])/2*(z[k+1]-z[k-1])/2
        aw = (x[i-1]+x[i])/2 * (y[j+1]-y[j-1])/2*(z[k+1]-z[k-1])/2
        an = (x[i+1]-x[i-1])/2*(z[k+1]-z[k-1])/2
        as = (x[i+1]-x[i-1])/2*(z[k+1]-z[k-1])/2
        at = (y[j+1]-y[j-1])/2*(((x[i+1]-x[i])/2)^2-((x[i]-x[i-1])/2)^2)/2
        ab = (y[j+1]-y[j-1])/2*(((x[i+1]-x[i])/2)^2-((x[i]-x[i-1])/2)^2)/2
        Source[i,j,k] = ρ*(ae*u[i+1,j,k]-aw*u[i,j,k]) + ρ*(an*v[i,j+1,k]-as*v[i,j,k]) + ρ*(at*w[i,j,k+1]-ab*w[i,j,k])
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

        u,v,w,P=setConditions(u,v,w,P)

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


        if (erro[end]/erro[1])<1e-3
            break
        end
        # u[end,:] = u[end-1,:]
    end
    return u,v,w,P,erro
end
