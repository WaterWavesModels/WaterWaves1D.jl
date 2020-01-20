export Data

mutable struct Data

    U :: Array{Array{Complex{Float64},2}}
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

