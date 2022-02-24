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
    h5write(joinpath( h5file * ".h5"), "/timesolver/isreal", isreal(ts.U1))

end

function load_timesolver( h5file )
    label = h5read(joinpath( h5file * ".h5"), "/timesolver/label")
    datasize = h5read(joinpath( h5file * ".h5"), "/timesolver/datasize")
    isreal = h5read(joinpath( h5file * ".h5"), "/timesolver/isreal")
    label == "Euler" && (return Euler((N=datasize[1],),datasize[2];realdata=isreal))
    label == "RK4" && (return RK4((N=datasize[1],),datasize[2];realdata=isreal))
end

Base.:(==)(a::TimeSolver, b::TimeSolver) = begin

     datasize_a = collect(size(a.U1))
     datasize_b = collect(size(b.U1))

     a.label == b.label && datasize_a ≈ datasize_b

end
