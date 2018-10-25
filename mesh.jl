struct Mesh
    
    nx   :: Int
    xmin :: Float64
    xmax :: Float64
    dx   :: Float64
    x    :: Vector{Float64}
    
    function Mesh( xmin, xmax, nx)
        dx = (xmax-xmin)/nx
        x = range(xmin, stop=xmax, length=nx+1)[1:end-1] 
        new( nx, xmin, xmax, dx, x)
    end
end
