import Base.dump

export dump, load_data, load_data!, load_init

"""
    dump(file_name :: String, data :: Data)

Save `data` to the file with name `file_name` (and extension ".h5").

"""
function dump( h5file :: String, data :: Data )
    l=length(data.U)
    n=length(first(data.U))
    N=length(first(data.U)[1])
    h5open(joinpath(h5file * ".h5"), "cw") do file  #"w" for write ("r" for read) 
        write(file,  "/data/size", N) 
        write(file,  "/data/num", n) 
        write(file, "/data/length", l )

        for i in 1:l 
            for j in 1:n
                write(file, "/data/U$i$j", data.U[i][j] )
            end
        end
    end

end
 
"""
    load_data(file_name :: String)

Load data from the file with name `file_name` (and extension ".h5").

"""
function load_data( h5file :: String )

    nsteps = h5read(h5file * ".h5", "/data/length")
    np = h5read(h5file * ".h5", "/data/size")
    nv = h5read(h5file * ".h5", "/data/num")

    v = [zeros(ComplexF64, np) for _ in 1:nv]

    data = Data(v)

    @assert length( data ) == 1
    @assert np == length(first(data.U)[1])
    @assert nv == length(first(data.U))
     
    data.U[1] .= [h5read(h5file * ".h5", "/data/U1$j") for j in 1:nv]
    for i in 2:nsteps
        push!(data.U, [h5read(h5file * ".h5", "/data/U$i$j") for j in 1:nv])
    end

    return data

end

"""
    dump(file_name  :: String, problem :: Problem)

Save the data of `problem` to the file with name `file_name` (and extension ".h5").

"""
function dump( h5file :: String, pb :: Problem) 
    dump( h5file, pb.data )
end

"""
    load_data!(file_name :: String, problem :: Problem)

Fills `problem` with raw data extracted from the file with name `file_name` (and extension ".h5").

"""
function load_data!( h5file :: String, pb::Problem )
    data = load_data( h5file )
    pb.data = data
end

"""
    dump(file_name :: String, x::Vector, init::InitialData )

Save the values of the initial data `init` at collocation points `x` to the file with name `file_name` (and extension ".h5").

"""
function dump( h5file :: String, x::Vector, init::InitialData) 
    h5open(joinpath(h5file * ".h5"), "cw") do file  #"w" for write ("r" for read) 

        write(file, "/init/x", x )
        write(file, "/init/η", init.η(x) )
        write(file, "/init/v", init.v(x) )
        write(file, "/init/label", init.label)

    end
end


"""
    load_init(file_name :: String; fast = false)

Load initial data from  the file with name `file_name` (and extension ".h5").

Keyword argument `fast` is optional (default is `false`), and corresponds the the keyword argument of the function `Init`.

Return an object of type `InitialData`.

"""
function load_init( h5file :: String ; fast = false) 

    x = h5read(joinpath( h5file * ".h5"), "/init/x")
    η = h5read(joinpath( h5file * ".h5"), "/init/η")
    v = h5read(joinpath( h5file * ".h5"), "/init/v")
    label = h5read(joinpath(h5file * ".h5"), "/init/label")


    return Init(x,η,v ; fast = fast, label = label)
end
