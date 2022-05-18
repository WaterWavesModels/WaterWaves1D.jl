export TimeSolver

abstract type TimeSolver end

show(io::IO, s::TimeSolver) =
    try print(io,s.info)
    catch
        print(io,"Time solver: $(s.label)")
    end


Base.:(==)(a::TimeSolver, b::TimeSolver) = begin
  
    init = Random((N=10, L=10))
    
    Na = try size(a.U1,1)
    catch
        2
    end
    ma = Airy((μ=1, N=Na, L=10))
    Ua = ma.mapto(init)

    Nb = try size(b.U1,1)
    catch
        Complex.(rand(2,2))
    end
    mb = Airy((μ=1, N=Nb, L=10))
    Ub = mb.mapto(init)


    step!(a, ma, Ua, 1) 
    step!(b, mb, Ub, 1) 

    a.label == b.label && Ua == Ub

end
