"""
Computes the solitary wave of the Whitham equation with prescribed velocity (more than 1, less than c=1.2290408)
"""
#using ShallowWaterModels
include("../src/dependencies.jl")

#----
function SWW(c)
 function solKdV(x,α)
   	 2*α*sech.(sqrt(3*2*α)/2*x).^2
 end
 if c<=1 || c>1.2290408
	 error("the velocity should be between `1` and `1.2290408`")
 elseif c<1.15
	param = ( μ  = 1,
 			ϵ  = 1,
         	N  = 2^10,
             L  = 10/sqrt(c-1),
 			)
 	mesh = Mesh(param)
	u = SolitaryWaveWhitham(mesh, merge(param,(c=c,)), solKdV(mesh.x,c-1); iterative = false)

 else
	param = ( μ  = 1,
			ϵ  = 1,
        	N  = 2^10,
            L  = 15,
			)
	mesh = Mesh(param)

	#
	# function solve!(u₀,u₁,c,dc)
	# 	while dc > 0.00001
	# 		new_dc = dc
	# 		@info string("c = ",c,"dc = ",dc,"\n")
	# 		u₀[1] .= u₁[1]
	# 		flag = 1
	# 		while flag==1
	# 			u₁ = SolitaryWaveWhitham(mesh, merge(param,(c=c+new_dc,)), u₁[1]+new_dc/dc*(u₁[1]-u₀[1]); iterative = true, verbose = false)
	# 			if u₁[2]==1
	# 				new_dc = new_dc/2
	# 			else
	# 				flag=0
	# 			end
	# 		end
	# 		c += new_dc
	# 		dc = new_dc
	# 	end
	# 	return (u₀,u₁,c,dc)
	# end

	#(u₀,u₁,c,dc)=solve!(u₀,u₁,c,dc)
	#c=1.101
	#dc = 0.01
	#
	# function solve0!(u₀,c,old_dc,new_dc)
	# 		u₀ = SolitaryWaveWhitham(mesh, merge(param,(c=c,)), u₀[1]; iterative = true, verbose = true)
	# end
	#
	# function solve1!(u₀,u₁,c,old_dc,new_dc)
	# 		tu = copy(u₁[1])
	# 		u₁ = SolitaryWaveWhitham(mesh, merge(param,(c=c,)), u₁[1]+new_dc/old_dc*(u₁[1]-u₀[1]); iterative = true, verbose = true)
	# 		u₀[1] .= tu
	# end

	# barycentric Lagrange inerpolation (4.2) in https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
	function P(X,dc,u₀,u₁,u₂)
		(u₀/(2*(X/dc+2))-u₁/(X/dc+1)+u₂/(2*X/dc))/(1/(2*(X/dc+2))-1/(X/dc+1)+1/(2*X/dc))
	end

	function Q(X,dc,u₀,u₁,u₂)
		(u₀./(2*(X/dc+2)) .-u₁./(X/dc+1) .+u₂./(2*X/dc))./(1/(2*(X/dc+2))-1/(X/dc+1)+1/(2*X/dc))
	end

	function solve2!(u₀,u₁,u₂,c,dc)
			tu₂ = copy(u₂)
			u₂ .= SolitaryWaveWhitham(mesh, merge(param,(c=c,)), P.(dc,dc,u₀,u₁,u₂); iterative = false, verbose = true, max_iter = 5)
			u₀ .= u₁
			u₁ .= tu₂
	end

	function solve2!(u₀,u₁,u₂,c,old_dc,new_dc)
			tu₀ = copy(u₀)
			tu₁ = copy(u₁)
			tu₂ = copy(u₂)
			u₀ .= tu₂
			u₁ .= SolitaryWaveWhitham(mesh, merge(param,(c=c+1*new_dc,)), P.(1*new_dc,old_dc,tu₀,tu₁,tu₂); iterative = false, verbose = true)
			u₂ .= SolitaryWaveWhitham(mesh, merge(param,(c=c+2*new_dc,)), P.(2*new_dc,old_dc,tu₀,tu₁,tu₂); iterative = false, verbose = true)
	end

	u₀ = SolitaryWaveWhitham(mesh, merge(param,(c=1.1,)), solKdV(mesh.x,0.1); iterative = false)
	u₁ = SolitaryWaveWhitham(mesh, merge(param,(c=1.11,)), solKdV(mesh.x,0.11); iterative = false)
	u₂ = SolitaryWaveWhitham(mesh, merge(param,(c=1.12,)), solKdV(mesh.x,0.12); iterative = false)

	if c < 1.22
		for cs = range(1.11; step = 0.01, stop = c)
			print(string("c = ",cs,"\n"))
			solve2!(u₀,u₁,u₂,cs,0.01)
		end
	else
		for cs = range(1.11; step = 0.01, stop = 1.22)
			print(string("c = ",cs,"\n"))
			solve2!(u₀,u₁,u₂,cs,0.01)
		end

		solve2!(u₀,u₁,u₂,1.221,0.01,1e-4)

		if c < 1.229
			for cs = range(1.221+3e-4; step = 1e-4, stop = c)
				print(string("c = ",cs,"\n"))
				solve2!(u₀,u₁,u₂,cs,1e-4)
			end
		else
			c=1.2290408
			for cs = range(1.221+3e-4; step = 1e-4, stop = 1.229)
				print(string("c = ",cs,"\n"))
				solve2!(u₀,u₁,u₂,cs,1e-4)
			end

			solve2!(u₀,u₁,u₂,1.229,1e-4,1e-6)


			for cs = range(1.229+3e-6; step = 1e-6, stop = 1.22904)
				print(string("c = ",cs,"\n"))
				solve2!(u₀,u₁,u₂,cs,1e-6)
			end

			solve2!(u₀,u₁,u₂,1.22904,1e-6,1e-8)

			for cs = range(1.22904+3e-8; step = 1e-8, stop = 1.2290407)
				print(string("c = ",cs,"\n"))
				solve2!(u₀,u₁,u₂,cs,1e-8)
			end
			for cs = range(1.2290407+1e-8; step = 1e-8, stop = 1.2290408)
				print(string("c = ",cs,"\n"))
				solve2!(u₀,u₁,u₂,cs,1e-8)
			end
		end
	end
	u = u₂
 end
 plt = plot(layout=(1,2))
 plot!(plt[1,1], mesh.x, [u solKdV(mesh.x,c-1)];
  title=string("c=",c),
  label=["Whitham" "KdV"])

 plot!(plt[1,2], fftshift(mesh.k),
  log10.(abs.(fftshift(fft(u))));
  title="frequency",
  label="Whitham")
end
