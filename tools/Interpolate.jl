export Interpolate
using ShallowWaterModels,FFTW

"""
    Interpolate(mesh,vector;k=2^4)

Interpolate a vector `vector` defined on a uniform collocation grid defined by `mesh`.

Returns `(new_mesh,new_vector)` a new uniform mesh with `n` times as many values, and the vector of values on collocation points.

"""
function Interpolate(mesh::Mesh,vector;n=2^3::Int)

    fourier=fft(vector)
    m=Int(mesh.N/2)
    new_fourier = [fourier[1:m] ;zeros((n-1)*2*m) ;fourier[m+1:end]]


    new_mesh = Mesh(mesh.xmin,mesh.xmax,n*mesh.N)
    new_vector=ifft(new_fourier)*n
    if all((vector[i] isa Real) for i in length(vector))
        new_vector=real.(new_vector)
    end
    return new_mesh,new_vector
end
"""
    Interpolate(mesh,vector;k=2^4)

Interpolate a vector `vector` defined on a uniform collocation grid defined by `mesh`, on collocation points given by `new_mesh`.

Returns `new_vector` the vector of values on collocation points.

"""
function Interpolate(mesh::Mesh,vector,new_mesh::Mesh)

    fourier=fft(vector)
    k = mesh.k
    x₀ = mesh.xmin

    new_vector=exp.(1im*(new_mesh.x.-x₀)*k')*fourier/length(k)
    if all((vector[i] isa Real) for i in length(vector))
        new_vector=real.(new_vector)
    end
    return new_vector
end
