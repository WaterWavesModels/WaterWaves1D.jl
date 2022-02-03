export TimeSolver

abstract type TimeSolver end

show(io::IO, s::TimeSolver) =
    try print(io,s.info)
    catch
        print(io,"Time solver: $(s.label)")
    end

function dump( h5file, ts::TimeSolver) 
    datasize = collect(size(ts.U1))
    h5write(joinpath( h5file * ".h5"), "/timesolver/label", ts.label)
    h5write(joinpath( h5file * ".h5"), "/timesolver/datasize", datasize)
end

function load_timesolver( h5file )  
    label = h5read(joinpath( h5file * ".h5"), "/timesolver/label")
    datasize = h5read(joinpath( h5file * ".h5"), "/timesolver/datasize")
    label == "Euler" && (return Euler(datasize))
    label == "RK4" && (return RK4(datasize...))
end

Base.:(==)(a::TimeSolver, b::TimeSolver) = begin

     datasize_a = collect(size(a.U1))
     datasize_b = collect(size(b.U1))

     a.label == b.label && datasize_a ≈ datasize_b

end

