export fig

function fig(t)
    s=0

    if indexin(false,times.t.<=t)[1]==nothing
        index=length(times.t)
        else index=indexin(false,times.t.<=t)[1]-1
    end
    t=times.t[index]
    p1 = plot(title="temps t=$t, ϵ=$epsilon")
    
    for modele in range(1,size(Us)[end])
    
        (hhat,uhat)=(Us[modele][:,1,index],Us[modele][:,2,index])
        (h,u)=(real(ifft((Gamma.^s).*hhat)),real(ifft(uhat)))

        
        p1 = plot!(mesh.x,h)
    end
        
    p2 = plot()
    
    for modele in range(1,size(Us)[end])
    
        (hhat,uhat)=(Us[modele][:,1,index],Us[modele][:,2,index])
        (h,u)=(real(ifft((Gamma.^s).*hhat)),real(ifft(uhat)))

        p2 = plot!(fftshift(mesh.k),log10.(1e-18.+abs.(fftshift(hhat))))
        
    end
    p=plot(p1,p2,layout=(2,1),label=Labels)

    p
end


function fig(t, times, Gamma, Modeles::Dict, epsilon, mesh)
        
    Labels = keys(Modeles)
    s = 0
    if indexin(false,times.t.<=t)[1]==nothing
        index=length(times.t)
    else 
        index=indexin(false,times.t.<=t)[1]-1
    end
    t=times.t[index]
    
    p = plot(layout=(2,1))
    
    for label in Labels
        (hhat,uhat)=Modeles[label][index]
        (h,u)=(real(ifft((Gamma.^s).*hhat)),real(ifft(uhat)))
        plot!(p[1,1], mesh.x,h; label=string(label))
        plot!(p[2,1], fftshift(mesh.k),log10.(1e-18.+abs.(fftshift(hhat))); label=string(label))  
    end
    
    p
end


function fig(t, times::Times, models, mesh::Mesh)
        
    s = 0
    index = length(times.t)
    t = times.t[index]
    
    p = plot(layout=(2,1), title = "Deep Water Models")

    hr = zeros(Float64, mesh.N)
    ur = zeros(Float64, mesh.N)
    
    for model in models
        (hhat,uhat) = model.data[index]
        hr .= real(ifft((model.Gamma.^s).*hhat))
        ur .= real(ifft(uhat))
        plot!(p[1,1], mesh.x,hr; label=model.label)
        plot!(p[2,1], fftshift(model.mesh.k),log10.(1e-18.+abs.(fftshift(hhat))); label=model.label)  
    end

    p

end

function plot_model(times::Times, model::AbstractModel, mesh::Mesh)
        
    s = 0
    index = length(times.t)
    t = times.t[index]
    
    p = plot(layout=(2,1))
    
    (hhat,uhat) = model.data[index]
    (hr,ur)     = (real(ifft((model.Gamma.^s).*hhat)),real(ifft(uhat)))
    plot!(p[1,1], mesh.x,hr; label=model.label)
    plot!(p[2,1], fftshift(model.mesh.k),log10.(1e-18.+abs.(fftshift(hhat))); label=model.label)  
    
    p

end

function plot_model!(p, times::Times, model::AbstractModel, mesh::Mesh)
        
    s = 0
    index = length(times.t)
    t = times.t[index]
    
    (hhat,uhat) = model.data[index]
    (hr,ur)     = (real(ifft((model.Gamma.^s).*hhat)),real(ifft(uhat)))
    plot!(p[1,1], mesh.x,hr; label=model.label)
    plot!(p[2,1], fftshift(model.mesh.k),log10.(1e-18.+abs.(fftshift(hhat))); label=model.label)  
    
end

