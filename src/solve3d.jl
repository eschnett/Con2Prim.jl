export solve3d
function solve3d(f′, xlo::SVector{3,T}, xhi::SVector{3,T}) where {T}
    # Check inputs
    @assert all(isfinite.(xlo) .&& isfinite.(xhi))
    @assert all(xlo .<= xhi)

    # Validating function evaluation
    function f(x)
        @assert all(isfinite.(x))
        local y = f′(x)::SVector{2,T}
        all(isfinite.(y)) || @show x y
        @assert all(isfinite.(y))
        return y
    end

    # Ensure we can find a zero (assuming there is one)
    @assert all(f(xlo) .<= 0 .<= f(xhi))

    # Find the zero of f[1] on the lower boundary of the domain
    if f(SVector(xhi[1], xhi[2], xlo[3]))[1] >= 0
        # Check lower xy boundary
        xeqnene1 = let
            xq(q) = SVector(q[1], q[2], xlo[3])
            fq(q) = f(xq(q))[1]
            q = solve1d(fq, SVector(xlo[1], xlo[2]), SVector(xhi[1], xhi[2])).x
            xq(q)
        end
    elseif f(SVector(xhi[1], xlo[2], xlo[3]))[1] >= 0
        # Check lower xz boundary
        xeqnene1 = let
            xq(q) = SVector(q[1], xlo[2], q[2])
            fq(q) = f(xq(q))[1]
            q = solve1d(fq, SVector(xlo[1], xlo[3]), SVector(xhi[1], xhi[3])).x
            xq(q)
        end
    elseif f(SVector(xhi[1], xlo[2], xlo[3]))[1] >= 0
        # Check lower yz boundary
        xeqnene1 = let
            xq(q) = SVector(xlo[1], q[1], q[2])
            fq(q) = f(xq(q))[1]
            q = solve1d(fq, SVector(xlo[2], xlo[3]), SVector(xhi[2], xhi[3])).x
            xq(q)
        end
    else
        # impossible
        @assert false
    end

    # Find the zero of f[1] on the upper boundary of the domain
    if f(SVector(xhi[1], xhi[2], xlo[3]))[1] >= 0
        # Check lower xy boundary
        xeqnene1 = let
            xq(q) = SVector(q[1], q[2], xlo[3])
            fq(q) = f(xq(q))[1]
            q = solve1d(fq, SVector(xlo[1], xlo[2]), SVector(xhi[1], xhi[2])).x
            xq(q)
        end
    elseif f(SVector(xhi[1], xlo[2], xlo[3]))[1] >= 0
        # Check lower xz boundary
        xeqnene1 = let
            xq(q) = SVector(q[1], xlo[2], q[2])
            fq(q) = f(xq(q))[1]
            q = solve1d(fq, SVector(xlo[1], xlo[3]), SVector(xhi[1], xhi[3])).x
            xq(q)
        end
    elseif f(SVector(xhi[1], xlo[2], xlo[3]))[1] <= 0
        # Check lower yz boundary
        xeqnene1 = let
            xq(q) = SVector(xlo[1], q[1], q[2])
            fq(q) = f(xq(q))[1]
            q = solve1d(fq, SVector(xlo[2], xlo[3]), SVector(xhi[2], xhi[3])).x
            xq(q)
        end
    else
        # impossible
        @assert false
    end

    flolohi = f(SVector(xlo[1], xlo[2], xhi[3]))
    flohilo = f(SVector(xlo[1], xhi[2], xlo[3]))
    flohihi = f(SVector(xlo[1], xhi[2], xhi[3]))
    fhilolo = f(SVector(xhi[1], xlo[2], xlo[3]))
    fhilohi = f(SVector(xhi[1], xlo[2], xhi[3]))
    fhihilo = f(SVector(xhi[1], xhi[2], xlo[3]))

    @assert false
end
