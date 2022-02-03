function dump( h5file, model::AbstractModel) 
    h5write(joinpath( h5file * ".h5"), "/model/label", model.label)
end

function load_model( h5file )  
    h5read(joinpath( h5file * ".h5"), "/model/label")
end
