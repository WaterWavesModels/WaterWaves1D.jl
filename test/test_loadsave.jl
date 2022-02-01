using Test
using WaterWaves1D
import WaterWaves1D: dump, load!, load_data

@testset "LoadSave" begin
 
    param  = ( ϵ  = 1/2, μ = 1, θ = 1)
    paramX = ( N  = 2^8, L  = 10)
    paramT = ( T  = 5, dt = 0.1)
    
    init  = BellCurve(param)
    model = WhithamGreenNaghdi(merge(param,paramX);ktol=1e-10)
    prob  = Problem( model, init, merge(paramT,paramX); solver = RK4(paramX) )

    solve!( prob )
    
    # model   :: AbstractModel
    # initial :: InitialData
    # param   :: NamedTuple
    # solver  :: TimeSolver
    # times   :: Times
    # mesh    :: Mesh
    # data    :: Data

    rm("testsave.h5", force=true)
    dump("testsave", prob.data)
    load!(prob.data, "testsave")

    data = load_data("testsave")

    @test data .≈ prob.data
 
end
