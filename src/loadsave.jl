export dump, load_data, load_data!, load_init

import Base.dump

function dump( h5file :: String, data :: Data )
    h5open(joinpath(h5file * ".h5"), "w") do file  #"w" for write ("r" for read) 
        write(file,  "/data/size", collect(size(first(data.U))) )
        write(file, "/data/length", length(data.U) )

        for (i,u) in enumerate(data.U)
            write(file, "/data/U$i", u )
        end
    end

end
 
 
function load_data( h5file :: String )

    nsteps = h5read(h5file * ".h5", "/data/length")
    np, nv = h5read(h5file * ".h5", "/data/size")

    v = zeros(ComplexF64, np, nv)

    data = Data(v)

    @assert length( data ) == 1
    @assert np == size(first(data.U))[1]
    @assert nv == size(first(data.U))[2]
     
    data.U[1] .= h5read(h5file * ".h5", "/data/U1")
    for i in 2:nsteps
        push!(data.U, h5read(h5file * ".h5", "/data/U$i"))
    end

    return data

end

function dump( h5file, pb::Problem) 
    dump( h5file, pb.data )
end

function load_data!( pb::Problem, h5file :: String )
    data = load_data( h5file )
    pb.data = data
end


function dump( h5file, x::Vector, init::InitialData) 
    h5open(joinpath(h5file * ".h5"), "w") do file  #"w" for write ("r" for read) 

        write(file, "/init/x", x )
        write(file, "/init/η", init.η(x) )
        write(file, "/init/v", init.v(x) )
        write(file, "/init/label", init.label)

    end
end



function load_init( h5file ; fast = false) 

    x = h5read(joinpath( h5file * ".h5"), "/init/x")
    η = h5read(joinpath( h5file * ".h5"), "/init/η")
    v = h5read(joinpath( h5file * ".h5"), "/init/v")
    label = h5read(joinpath(h5file * ".h5"), "/init/label")


    return Init(x,η,v ; fast = fast, label = label)
end
