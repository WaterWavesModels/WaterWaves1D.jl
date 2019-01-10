export create_animation,fig_problem!

function create_animation( p::Problem )

    prog = Progress(p.times.Nt,1)

    anim = @animate for l in range(1,p.times.Nt-1)

        pl = plot(layout=(2,1))

		(hr,ur) = mapfro(p.model,p.data.U[l])

        plot!(pl[1,1], p.mesh.x, hr;
              ylims=(-0.6,1),
              title="physical space",
              label=p.model.label)

        plot!(pl[2,1], fftshift(p.mesh.k),
              log10.(1e-18.+abs.(fftshift(fft(hr))));
              title="frequency",
              label=p.model.label)

        next!(prog)

    end when mod(l, 200) == 0

    gif(anim, "anim.gif", fps=15); nothing

end

#----
#
#
#md # Plot results function
function fig_problem!( plt, p::Problem )

    (hr,ur) = mapfro(p.model,p.data.U[end])

    plot!(plt[1,1], p.mesh.x, hr;
		  title="physical space",
	          label=p.model.label)

    plot!(plt[2,1], fftshift(p.mesh.k),
                  log10.(1e-18.+abs.(fftshift(fft(hr))));
		  title="frequency",
    	          label=p.model.label)

end

function fig_problem!( plt, p::Problem, t::Real )

	t=max(t,0)
	t=min(t,p.times.tfin)
	index=1+round(Int,(t/p.times.tfin)*(p.times.Nt-1))


    (hr,ur) = mapfro(p.model,p.data.U[index])

    plot!(plt[1,1], p.mesh.x, hr;
		  title=string("physical space at t=",p.times.t[index]),
	          label=p.model.label)

    plot!(plt[2,1], fftshift(p.mesh.k),
                  log10.(1e-18.+abs.(fftshift(fft(hr))));
		  title="frequency",
    	          label=p.model.label)

end

function norm_problem!( plt, p::Problem, s::Real )
	N=[];
	Λ = sqrt.((p.mesh.k.^2).+1);
 	for index in range(1,stop=p.times.Nt)
    	(hr,ur) = mapfro(p.model,p.data.U[index])
		push!(N,norm(ifft(Λ.^s.*fft(hr))))
	end

    plot!(plt, p.times.t, N;
		  title=string("norm H^s avec s=",s),
	          label=p.model.label)

end


#
# function fig(t)
#     s=0
#
#     if indexin(false,times.t.<=t)[1]==nothing
#         index=length(times.t)
#         else index=indexin(false,times.t.<=t)[1]-1
#     end
#     t=times.t[index]
#     p1 = plot(title="temps t=$t, ϵ=$epsilon")
#
#     for modele in range(1,size(Us)[end])
#
#         (hhat,uhat)=(Us[modele][:,1,index],Us[modele][:,2,index])
#         (h,u)=(real(ifft((Gamma.^s).*hhat)),real(ifft(uhat)))
#
#
#         p1 = plot!(mesh.x,h)
#     end
#
#     p2 = plot()
#
#     for modele in range(1,size(Us)[end])
#
#         (hhat,uhat)=(Us[modele][:,1,index],Us[modele][:,2,index])
#         (h,u)=(real(ifft((Gamma.^s).*hhat)),real(ifft(uhat)))
#
#         p2 = plot!(fftshift(mesh.k),log10.(1e-18.+abs.(fftshift(hhat))))
#
#     end
#     p=plot(p1,p2,layout=(2,1),label=Labels)
#
#     p
# end
#
#
# function fig(t, times, Gamma, Modeles::Dict, epsilon, mesh)
#
#     Labels = keys(Modeles)
#     s = 0
#     if indexin(false,times.t.<=t)[1]==nothing
#         index=length(times.t)
#     else
#         index=indexin(false,times.t.<=t)[1]-1
#     end
#     t=times.t[index]
#
#     p = plot(layout=(2,1))
#
#     for label in Labels
#         (hhat,uhat)=Modeles[label][index]
#         (h,u)=(real(ifft((Gamma.^s).*hhat)),real(ifft(uhat)))
#         plot!(p[1,1], mesh.x,h; label=string(label))
#         plot!(p[2,1], fftshift(mesh.k),log10.(1e-18.+abs.(fftshift(hhat))); label=string(label))
#     end
#
#     p
# end
#
#
# function fig(t, times::Times, models, mesh::Mesh)
#
#     s = 0
#     index = length(times.t)
#     t = times.t[index]
#
#     p = plot(layout=(2,1), title = "Deep Water Models")
#
#     hr = zeros(Float64, mesh.N)
#     ur = zeros(Float64, mesh.N)
#
#     for model in models
#         (hhat,uhat) = model.data[index]
#         hr .= real(ifft((model.Gamma.^s).*hhat))
#         ur .= real(ifft(uhat))
#         plot!(p[1,1], mesh.x,hr; label=model.label)
#         plot!(p[2,1], fftshift(model.mesh.k),log10.(1e-18.+abs.(fftshift(hhat))); label=model.label)
#     end
#
#     p
#
# end
#
# function plot_model(times::Times, model::AbstractModel, mesh::Mesh)
#
#     s = 0
#     index = length(times.t)
#     t = times.t[index]
#
#     p = plot(layout=(2,1))
#
#     (hhat,uhat) = model.data[index]
#     (hr,ur)     = (real(ifft((model.Gamma.^s).*hhat)),real(ifft(uhat)))
#     plot!(p[1,1], mesh.x,hr; label=model.label)
#     plot!(p[2,1], fftshift(model.mesh.k),log10.(1e-18.+abs.(fftshift(hhat))); label=model.label)
#
#     p
#
# end
#
# function plot_model!(p, times::Times, model::AbstractModel, mesh::Mesh)
#
#     s = 0
#     index = length(times.t)
#     t = times.t[index]
#
#     (hhat,uhat) = model.data[index]
#     (hr,ur)     = (real(ifft((model.Gamma.^s).*hhat)),real(ifft(uhat)))
#     plot!(p[1,1], mesh.x,hr; label=model.label)
#     plot!(p[2,1], fftshift(model.mesh.k),log10.(1e-18.+abs.(fftshift(hhat))); label=model.label)
#
# end
