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

    U :: Vector{AbstractArray}
    datasize :: Int
    datalength :: Int

    function Data( v :: AbstractArray )

        dim = length(size(v))
        if dim == 1 # if v is a vector of data (η,vx,vy...)
            datasize = size(v,1)
            datalength = size(v[1],1)
        elseif dim == 2 # data are columns of v
            (datalength , datasize ) = size(v)
        elseif dim == 3 # data are in a three-dimensional array
            datasize = size(v,3)
            datalength = size(v,1)
        end
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