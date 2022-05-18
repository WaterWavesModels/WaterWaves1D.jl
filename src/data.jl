export Data

"""
    Data( m :: Matrix )

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

Base.:length( data :: Data ) = length( data.U )

Base.:size( data :: Data ) = size( first(data.U) )

Base.:(==)(d1::Data, d2::Data) =
    d1.datalength == d2.datalength &&
    d1.datasize == d2.datasize &&
    d1.U == d2.U