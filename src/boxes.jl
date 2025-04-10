export Bounds
struct Bounds{D,S}
    bounds::NTuple{D,Vector{S}}

    function Bounds{D,S}() where {D,S}
        bounds = Tuple(S[-Inf, +Inf] for d in 1:D)
        return new{D,S}(bounds)
    end
end

export add_point!
function add_point!(bounds::Bounds{D}, p::SVector{D}) where {D}
    for d in 1:D
        push!(bounds.bounds[d], p[d])
    end
    return nothing
end

export Boxes
struct Boxes{D,S,T}
    bounds::NTuple{D,Vector{S}}
    values::Array{T,D}

    function Boxes{D,S,T}(bounds::Bounds{D,S}) where {D,S,T}
        bounds = bounds.bounds
        for d in 1:D
            sort!(bounds[d])
            unique!(bounds[d])
        end
        values = zeros(T, length.(bounds) .- 1)
        return new{D,S,T}(bounds, values)
    end
end

export map_past!
function map_past!(f, boxes::Boxes{D}, p::SVector{D}) where {D}
    i = ntuple(D) do d
        r = searchsorted(boxes.bounds[d], p[d])
        @assert !isempty(r)
        first(r)
    end
    i = CartesianIndex(i .- 1)
    b = CartesianIndex(Tuple(1 for d in 1:D))
    for j in b:i
        boxes.values[j] = f(boxes.values[j])
    end
    return nothing
end

export map_future!
function map_future!(f, boxes::Boxes{D}, p::SVector{D}) where {D}
    i = ntuple(D) do d
        r = searchsorted(boxes.bounds[d], p[d])
        @assert !isempty(r)
        first(r)
    end
    i = CartesianIndex(i)
    e = CartesianIndex(size(boxes.values))
    for j in i:e
        boxes.values[j] = f(boxes.values[j])
    end
    return nothing
end

export set_past!
set_past!(boxes::Boxes{D}, p::SVector{D}, v) where {D} = map_past!(Returns(v), boxes, p)

export set_future!
set_future!(boxes::Boxes{D}, p::SVector{D}, v) where {D} = map_future!(Returns(v), boxes, p)
