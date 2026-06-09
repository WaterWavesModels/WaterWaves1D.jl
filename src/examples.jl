export examples_dir, get_examples, default_example, convergence_test, data_dir
"""
    examples_dir()

Return the directory where the example files provided with WaterWaves1D.jl are located. 
If WaterWaves1D.jl is installed as a regular package (with `]add DispersiveShallowWater`), 
these files are read-only and should *not* be modified. 
To find out which files are available, use, e.g., `readdir`.

Copied from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).

# Example
```@example
readdir(examples_dir())
```
"""
examples_dir() = pkgdir(WaterWaves1D, "examples")::String

"""
    get_examples()

Return a list of all examples that are provided by DispersiveShallowWater.jl. See also
[`examples_dir`](@ref) and [`default_example`](@ref).
To run an example, use, e.g., `include`.


Copied from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).

# Example
```@example
include(get_examples()[1])
```
"""
function get_examples()
    examples = String[]
    for (root, dirs, files) in walkdir(examples_dir())
        for f in files
            if endswith(f, ".jl") && !startswith(f, "Study")
                push!(examples, joinpath(root, f))
            end
        end
    end

    return examples
end

"""
    default_example()

Return the path to an example that can be used to quickly see DispersiveShallowWater.jl in action.
See also [`examples_dir`](@ref) and [`get_examples`](@ref).

Copied from [Trixi.jl](https://github.com/trixi-framework/Trixi.jl).

To run an example, use, e.g., `include`.
# Example
```@example
include(default_example())
```
"""
function default_example()
    return joinpath(examples_dir(), "default.jl")
end