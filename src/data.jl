export Data

"""
Data structure to store the solution of the problem along time
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
