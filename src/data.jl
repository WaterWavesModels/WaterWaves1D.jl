export Data

"""
    Data( mm :: Matrix )

Data structure to store the solution of an initial-value problem along time.

`data=Data(m)` is of parametric type and offers
- `data.U`, a 1-element vector with a copy of the matrix `m`;
- `(data.datalength,data.datasize)=size(m)`  where \
`datalength` is the number of computed modes, and `datasize` the number of involved equations, typically 2.
"""
mutable struct Data

    U :: Vector{Array{ComplexF64,2}}
    datasize :: Int
    datalength :: Int

    function Data( v :: AbstractArray )

        (datalength , datasize ) = size(v)

        U = [copy(v)]

        new(U, datasize, datalength)

    end

end

import Base.length

length( data :: Data ) = length( data.U )

import Base.size

size( data :: Data ) = size( first(data.U) )

function dump( h5file :: String, data :: Data )

    h5write(h5file * ".h5", "/data/size", collect(size(first(data.U))))
    h5write(h5file * ".h5", "/data/length", length(data.U))

    for (i,u) in enumerate(data.U)
       h5write(h5file * ".h5", "/data/U$i", u)
    end

end

 
 
function load_data( h5file :: String )

    nsteps = h5read(h5file * ".h5", "/data/length")
    np, nv = h5read(h5file * ".h5", "/data/size")

    v = zeros(ComplexF64, np, nv)

    data = Data(v)

    @assert length( data ) == 1
    @assert np == size(first(data.U))[1]
    @assert nv == size(first(data.U))[2]
     
    data.U[1] .= h5read(h5file * ".h5", "/data/U1")
    for i in 2:nsteps
        push!(data.U, h5read(h5file * ".h5", "/data/U$i"))
    end

    return data

end
