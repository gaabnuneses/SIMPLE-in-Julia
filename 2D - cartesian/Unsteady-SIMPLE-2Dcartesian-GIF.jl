using Plots

function meshgrid(x,y)
    X = zeros(length(x),length(y))
    Y = zeros(length(x),length(y))

    for i in 1:length(x), j in 1:length(y)
        X[i,j] = x[i]
        Y[i,j] = y[j]
    end
    return X,Y
end

function makeGIF(x,y,Uk,Vk,Pk)
    X = -5 .+ 10*rand(4000)#-3:.1:3
    Y = 2*rand(4000)#0:.05:1
    # X,Y = meshgrid(X,Y)
    anim = @animate for k in 1:(size(Uk,3)-1)

        plot()
        P′ = Pk[:,:,k]
        P′[10:20,7:13] = NaN*zeros(length(P′[10:20,7:13]))
        plot!(x[2:end-1],y[2:end-1],P′[2:end-1,2:end-1]',st=:heatmap,cam=(0,90),clim=(0,2.5),aspect_ratio=:equal,size=(1000,250),cbar=false)
        plot!(xlim=(minimum(x),maximum(x)),ylim=(minimum(y),maximum(y)))
        for i in 1:length(X)
            ix = findmin(abs.(X[i] .-x))[2]
            iy = findmin(abs.(Y[i] .-y))[2]
            Ux = Uk[ix,iy,k]
            Uy = Vk[ix,iy,k]
            X[i] = X[i] + 10*dt*Ux
            Y[i] = Y[i] + 10*dt*Uy
            if Ux==0 && Uy ==0
                X[i] = NaN
                Y[i] = NaN
            end
        end
        scatter!(X,Y,color=:cyan,leg=false)
    end

    gif(anim,"Flow.gif",fps=5)
end
