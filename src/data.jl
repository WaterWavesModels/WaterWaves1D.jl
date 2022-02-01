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


function dump( h5file :: String, data :: Data )

    h5write(h5file * ".h5", "/data/size", size(first(data)))
    h5write(h5file * ".h5", "/data/length", length(data))

    @show size(data.U)
    for (i,u) in enumerate(data.U)
       h5write(h5file * ".h5", "/data/U$i", u)
    end

end

function load!( data :: Data, h5file :: String )

    n = h5read(h5file * ".h5", "data/size")

    @assert n == length(data.U)

    for i in 1:n
       U = h5read(h5file * ".h5", "data/U$i")
       @show size(U)
    end

end

function load_data( h5file :: String )

    datasize = h5read(h5file * ".h5", "data/size")
    datalength = h5read(h5file * ".h5", "data/length")

    v = zeros(ComplexF64, datalength, datasize)

    data = Data(v)

    for i in 1:size
       push!(data.U, h5read(h5file * ".h5", "data/U$i"))
    end

    return data

end
