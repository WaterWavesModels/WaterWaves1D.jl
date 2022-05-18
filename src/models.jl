export AbstractModel
"""
Abstract type whose subtypes are the models from which initial-value problems can be built,
through `Problem( model :: AbstractModel, initial :: InitialData, param :: NamedTuple )`
"""
abstract type AbstractModel end
show(io::IO, m::AbstractModel) =
    try print(io,m.info)
    catch
        print(io,"Model: $(m.label)")
    end

function Base.:(==)(m1::AbstractModel, m2::AbstractModel)
    x=rand(2)
    init=Random(x)
    U1=m1.mapto(init);U2=m2.mapto(init);
    data1=m1.mapfro(U1);data2=m2.mapfro(U2);
    m1.f!(U1);m2.f!(U2);

    U1 == U2 && data1 == data2 && m1.label == m2.label
end