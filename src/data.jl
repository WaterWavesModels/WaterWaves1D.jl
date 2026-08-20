export Data

"""
    Data( m :: Matrix )

Data structure to store the solution of an initial-value problem along time.

`data=Data(m)` is of parametric type and offers
- `data.U`, a 1-element vector with a copy of the matrix `m`;
- `(data.datalength,data.datasize)=size(m)`  where \
`datalength` is the number of computed modes, and `datasize` the number of involved equations, typically 2.
"""
struct Data{T, N}

    U::Vector{Vector{Array{T, N}}}
    datasize::Int
    datalength::Int

    function Data(v::Vector{Array{T, N}}) where {T, N}

        datasize = size(v, 1)
        datalength = size(v[1], 1)
        U = [copy(v)]

        return new{T, N}(U, datasize, datalength)

    end

    function Data(v::Array{T, N}) where {T, N}

        datalength, datasize = size(v)
        U = [[v[:, j] for j in axes(v, 2)]]

        return new{T, 1}(U, datasize, datalength)

    end

end


Base.:length(data::Data) = length(data.U)

Base.:size(data::Data) = size(first(data.U))

Base.:(==)(d1::Data, d2::Data) =
    d1.datalength == d2.datalength &&
    d1.datasize == d2.datasize &&
    d1.U == d2.U
