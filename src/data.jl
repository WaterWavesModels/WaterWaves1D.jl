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

    U
    datasize :: Int
    datalength :: Int

    function Data( v )

        (datalength , datasize ) = size(v)

        U = [copy(v)]

        new(U, datasize, datalength)

    end

end

import Base.length

length( data :: Data ) = length( data.U )
