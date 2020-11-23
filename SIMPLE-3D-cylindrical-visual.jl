using Plots

let
    plot(leg=false,aspect_ratio=:equal)
    X = [x[1]*cos(j) for j in y]
    Y = [x[1]*sin(j) for j in y]
    plot!(X,Y,X*0 .+z[1],color=:black,lw=1)
    plot!(X,Y,X*0 .+z[end],color=:black,lw=1)
    # plot!([X';X'],[Y';Y'],[(X'*0 .+z[1]);(X'*0 .+z[end])],color=:black,lw=1)

    X = [x[end]*cos(j) for j in y]
    Y = [x[end]*sin(j) for j in y]
    plot!(X,Y,X*0 .+z[1],color=:black,lw=1)
    plot!(X,Y,X*0 .+z[end],color=:black,lw=1)
    plot!([X';X'],[Y';Y'],[(X'*0 .+z[1]);(X'*0 .+z[end])],color=:black,lw=1)

    # Inner cylinder
    X = [.4*cos(j) for j in y]
    Y = [.4*sin(j) for j in y]
    plot!(X,Y,X*0 .+.2,color=:black,lw=1)
    plot!(X,Y,X*0 .+z[end],color=:black,lw=1)
    plot!([X';X'],[Y';Y'],[(X'*0 .+.2);(X'*0 .+z[end])],color=:black,lw=1)


end

plot!(collect(x),collect(x).*0,P[:,2,:],st=:surface)

function getXYZ()
    X = [];Y = [];ZETA = [];
    xc = [];yc = [];zc = [];
    for i in [2;nx]
        for j in 2:ny
            for k in 2:nz
                R = x[i]
                Θ = y[j]
                Z = z[k]
                δre = (-x[i]+x[i+1])/2
                δrw = (x[i]-x[i-1])/2
                δθn = (-y[j]+y[j+1])/2
                δθs = (y[j]-y[j-1])/2
                δzt = (-z[k]+z[k+1])/2
                δzb = (z[k]-z[k-1])/2
                xc = [xc;R*cos(Θ)];
                yc = [yc;R*sin(Θ)];
                zc = [zc;Z];
                X = [X;(R+δre)*cos(Θ+δθn);
                       (R+δre)*cos(Θ+δθn);
                       (R-δrw)*cos(Θ+δθn);
                       (R-δrw)*cos(Θ+δθn);
                       (R+δre)*cos(Θ+δθn);NaN;
                       (R+δre)*cos(Θ-δθs);
                       (R+δre)*cos(Θ-δθs);
                       (R-δrw)*cos(Θ-δθs);
                       (R-δrw)*cos(Θ-δθs);
                       (R+δre)*cos(Θ-δθs);NaN;
                       (R-δrw)*cos(Θ+δθn);
                       (R-δrw)*cos(Θ-δθs);NaN;
                       (R+δre)*cos(Θ+δθn);
                       (R+δre)*cos(Θ-δθs);NaN;
                       (R+δre)*cos(Θ+δθn);
                       (R+δre)*cos(Θ-δθs);NaN;
                       (R-δrw)*cos(Θ+δθn);
                       (R-δrw)*cos(Θ-δθs);NaN]

               Y = [Y;(R+δre)*sin(Θ+δθn);
                      (R+δre)*sin(Θ+δθn);
                      (R-δrw)*sin(Θ+δθn);
                      (R-δrw)*sin(Θ+δθn);
                      (R+δre)*sin(Θ+δθn);NaN;
                      (R+δre)*sin(Θ-δθs);
                      (R+δre)*sin(Θ-δθs);
                      (R-δrw)*sin(Θ-δθs);
                      (R-δrw)*sin(Θ-δθs);
                      (R+δre)*sin(Θ-δθs);NaN;
                      (R-δrw)*sin(Θ+δθn);
                      (R-δrw)*sin(Θ-δθs);NaN;
                      (R+δre)*sin(Θ+δθn);
                      (R+δre)*sin(Θ-δθs);NaN;
                      (R+δre)*sin(Θ+δθn);
                      (R+δre)*sin(Θ-δθs);NaN;
                      (R-δrw)*sin(Θ+δθn);
                      (R-δrw)*sin(Θ-δθs);NaN]
              ZETA = [ZETA;(Z-δzb);
                     (Z+δzt);
                     (Z+δzt);
                     (Z-δzb);
                     (Z-δzb);NaN;
                     (Z-δzb);
                     (Z+δzt);
                     (Z+δzt);
                     (Z-δzb);
                     (Z-δzb);NaN;
                     (Z+δzt);
                     (Z+δzt);NaN;
                     (Z-δzb);
                     (Z-δzb);NaN;
                     (Z+δzt);
                     (Z+δzt);NaN;
                     (Z-δzb);
                     (Z-δzb);NaN]
            end
        end
    end
    return X,Y,ZETA,xc,yc,zc
end

X,Y,Z,xc,yc,zc=getXYZ()

scatter(xc,yc,zc,m=2.5,cam=(30,60))
plot(X,Y,Z,cam=(30,60))
