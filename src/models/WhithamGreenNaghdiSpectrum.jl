export Spectrum

"""
    Spectrum(model,η,u,c,ν;kwargs)

Computes the spectrum about a specified travelling solution of the
the linearized Green-Naghdi model proposed by V. Duchêne, S. Israwi and R. Talhouk.

# Argument
- `model` should be of type `WhithamGreenNaghdi`,
	constructed e.g. with `model = WhithamGreenNaghdi(param;SGN=true)` where
	e.g. `param = ( μ  = 0.1, ϵ  = 1, N  = 2^9, L  = 10*π)`

- `(η,u)` are real vectors of the aforementioned solution where
where

    - `η` is the surface deformation;
    - `u` corresponds to the layer-averaged velocity.
- 'c' is the velocity of the solution,
- 'ν' is a Floquet coefficient.

## Keywords
- `SGN :: Bool`: if `true` computes the Serre-Green-Naghdi (SGN) instead of Whitham-Green-Naghdi (WGN) system (default is `false`);

# Return values
an Eigen factorization object `L` which
contains the eigenvalues in `L.values` and the
eigenvectors in the columns of the matrix `L.vectors`.
(The `k`th eigenvector can be obtained from the slice `F.vectors[:, k]``.)

"""
function Spectrum(m::WhithamGreenNaghdi,η,u,c,ν)
	h = 1 .+ m.ϵ*η
	DxF(v) = ifft(m.Π⅔.* m.F₀ .* fft(v))
	v= u - 1/3 ./h .* (DxF(h.^3 .*DxF(u)))

	Dx = m.IFFT*Diagonal(m.Π⅔.*m.∂ₓ)*m.FFT + ν*m.Id#Diagonal(m.x./sqrt.(1 .+ m.x.^2))
	DxF₀ = m.IFFT*Diagonal(m.Π⅔.*m.F₀)*m.FFT + ν*m.Id#Diagonal(m.x./sqrt.(1 .+ m.x.^2))

	A11 = -Dx*Diagonal(m.ϵ*u)
	A12 = -Dx*Diagonal(h)
	A13 = zeros(size(A11))
	A21 = -Dx* ( m.Id - Diagonal( h.*(DxF(m.ϵ*u).^2) ) )
	A22 = -Dx*( m.ϵ*Diagonal( v - u ) - Diagonal( (h.^2).*DxF(m.ϵ*u)) *DxF₀ )
	A23 = -Dx*Diagonal(m.ϵ*u)

	B1 = -Diagonal(1 ./h)*DxF₀*Diagonal( (h.^2).*DxF(m.ϵ*u) ) + 1/3*Diagonal(DxF(h.^3 .* DxF(m.ϵ*u))./h.^2)
	B2 = m.Id - 1/3*Diagonal(1 ./h)*DxF₀*Diagonal(h.^3)*DxF₀

	C2 = inv(B2)
	C1 = -C2*B1

	M11 = A11+A12*C1 + c*Dx
	M12 = A13+A12*C2
	M21 = A21+A22*C1
	M22 = A23+A22*C2 + c*Dx

	#Up = Diagonal(exp.(1im*0*m.x))
	#Um = Diagonal(exp.(-1im*0*m.x))
	Up = m.Id; Um=m.Id;
	return eigen([ Um*M11*Up Um*M12*Up ; Um*M21*Up Um*M22*Up ])
end
