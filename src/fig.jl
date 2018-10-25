export fig

function fig(t)
    s=0
    if indexin(false,times.t.<=t)[1]==nothing
        index=length(times.t)
        else index=indexin(false,times.t.<=t)
    end
    
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

        p2 = plot!(fftshift(freq.k),log10.(1e-18.+abs.(fftshift(hhat))))
        
    end
    p=plot(p1,p2,layout=(2,1),label=Labels)

    display(p)
end


function fig(t, times, Gamma, freq, Modeles::Dict, epsilon, mesh)
        
    Labels = keys(Modeles)
    s = 0
    if indexin(false,times.t.<=t)[1]==nothing
        index=length(times.t)
    else 
        index=indexin(false,times.t.<=t)
    end
    
    p = plot(layout=(2,1))
    
    for label in Labels
        (hhat,uhat)=Modeles[label][index]
        (h,u)=(real(ifft((Gamma.^s).*hhat)),real(ifft(uhat)))
        plot!(p[1,1], mesh.x,h; label=string(label))
        plot!(p[2,1], fftshift(freq.k),log10.(1e-18.+abs.(fftshift(hhat))); label=string(label))  
    end
    
    display(p)
end


function fig(t, times::Times, models, mesh::Mesh)
        
    s = 0
    if indexin(false,times.t.<=t)[1]==nothing
        index=length(times.t)
    else 
        index=indexin(false,times.t.<=t)
    end
    
    p = plot(layout=(2,1))
    
    for model in models
        (hhat,uhat)=model.data[index]
        (h,u)=(real(ifft((model.Gamma.^s).*hhat)),real(ifft(uhat)))
        plot!(p[1,1], mesh.x,h; label=model.label)
        plot!(p[2,1], fftshift(model.freq.k),log10.(1e-18.+abs.(fftshift(hhat))); 
            label=model.label)  
    end
    
    display(p)
end

