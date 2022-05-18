export InitialData

"""
Abstract type defining initial data from which initial-value problems can be built,
through `Problem( model :: AbstractModel, initial :: InitialData, param :: NamedTuple )`
"""
abstract type InitialData end

show(io::IO, i::InitialData) =
    try print(io,i.info)
    catch
        print(io,"Initial data: $(i.label)")
    end

function Base.:(==)(i1::InitialData, i2::InitialData)
    x=rand(2)
    i1.η(x) == i2.η(x) &&
    i1.v(x) == i2.v(x) &&
    i1.label == i2.label
end