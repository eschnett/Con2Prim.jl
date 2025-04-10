chop(x::Real) = ifelse(abs(x) < eps(typeof(x))^(3/4), zero(x), x)
round6(x::Real) = round(x; digits=6)

export solve2d

function solve2d(f′, xlo::SVector{2,T}, xhi::SVector{2,T}) where {T}
    D = 2

    # Check inputs
    @assert all(isfinite.(xlo) .&& isfinite.(xhi))
    @assert all(xlo .<= xhi)

    # Validating function evaluation
    function f(x::SVector{2,T})
        D = 2
        all(isfinite.(x)) || @show x
        @assert all(isfinite.(x))
        local y = f′(x)::SVector{D,T}
        all(isfinite.(y)) || @show x y
        @assert all(isfinite.(y))
        local h = 1.0e-6
        local y1 = f′(x + SVector(h, 0))::SVector{D,T}
        local y2 = f′(x + SVector(0, h))::SVector{D,T}
        all(y1 .> y) || @show x y 1 (y1-y)/h
        all(y2 .> y) || @show x y 2 (y2-y)/h
        @assert all(y1 .> y)
        @assert all(y2 .> y)
        return y::SVector{D,T}
    end

    # println("solve2d:")
    # println("    Inputs:")
    # println("        xlo=$(round6.(xlo))   f=$(round6.(f(xlo)))")
    # println("        xhi=$(round6.(xhi))   f=$(round6.(f(xhi)))")

    # Ensure we can find a zero (assuming there is one)
    @assert all(f(xlo) .<= 0 .<= f(xhi))

    # Find the zeros of f on the "left" and "right" boundaries of the domain
    # println("    Find bracketing struts:")
    # xzeros[ltrt][val]::SVector{D,T}
    xzeros = ntuple(2) do ltrt
        let
            xcorners = (xlo, ltrt==1 ? SVector(xlo[1], xhi[2]) : SVector(xhi[1], xlo[2]), xhi)
            fcorners = Tuple(f(x) for x in xcorners)
            xzeros = ntuple(D) do val
                @assert fcorners[1][val] <= 0
                @assert fcorners[3][val] >= 0
                if fcorners[2][val] >= 0
                    # Find zero on first edge
                    xmin, xmax = xcorners[1], xcorners[2]
                else
                    # Find zero on second edge
                    xmin, xmax = xcorners[2], xcorners[3]
                end
                xzero = let
                    xq(q) = (1-q) * xmin + q * xmax
                    fq(q) = f(xq(q))[val]
                    q = solve1d(fq, 0.0, 1.0).x
                    xq(q)
                end
                @assert abs(f(xzero)[val]) <= 10*eps(T)
                return xzero
            end
            # Find point closest to diagonal through centre
            xmid = (xlo + xhi) / 2
            val = findmin(x -> abs((x - xmid)[1] - (x - xmid)[2]), xzeros)[2]
            xval = xzeros[val]
            # Draw a line along the diagonal, find intersections with domain boundary
            λ = minimum(xval - xlo)
            xmin = xval - λ * SVector(1, 1)
            λ = minimum(xhi - xval)
            xmax = xval + λ * SVector(1, 1)
            # Find the other zero
            xval′ = let
                xq(q) = (1-q) * xmin + q * xmax
                fq(q) = f(xq(q))[3 - val]
                q = solve1d(fq, 0.0, 1.0).x
                xq(q)
            end
            val == 1 ? (xval, xval′) : (xval′, xval)
        end
    end
    xzeros = (lt=xzeros[1], rt=xzeros[2])
    for ltrt in 1:2
        Δx = xzeros[ltrt][2][1] - xzeros[ltrt][1][1]
        Δy = xzeros[ltrt][2][2] - xzeros[ltrt][1][2]
        @assert abs(Δy - Δx) <= 10*eps(T)
    end
    # # println("        xzeros.lt[1]=$(round6.(xzeros.lt[1]))   f=$(round6.(f(xzeros.lt[1])))")
    # # println("        xzeros.lt[2]=$(round6.(xzeros.lt[2]))   f=$(round6.(f(xzeros.lt[2])))")
    # # println("        xzeros.rt[1]=$(round6.(xzeros.rt[1]))   f=$(round6.(f(xzeros.rt[1])))")
    # # println("        xzeros.rt[2]=$(round6.(xzeros.rt[2]))   f=$(round6.(f(xzeros.rt[2])))")

    # xzeros[val][ltgt]::SVector{2,T}
    xzeros = ntuple(D) do val
        if f(xzeros.lt[val])[3 - val] <= 0 && f(xzeros.rt[val])[3 - val] >= 0
            (lt=xzeros.lt[val], gt=xzeros.rt[val])
        elseif f(xzeros.rt[val])[3 - val] <= 0 && f(xzeros.lt[val])[3 - val] >= 0
            (lt=xzeros.rt[val], gt=xzeros.lt[val])
        else
            # No proof we're bracketing a root
            @assert false
        end
    end
    for val in 1:D
        # println("        xzeros[$val].lt=$(round6.(xzeros[val].lt))   f=$(round6.(f(xzeros[val].lt)))")
        # println("        xzeros[$val].gt=$(round6.(xzeros[val].gt))   f=$(round6.(f(xzeros[val].gt)))")
    end

    return solve2d(f, xzeros)
end

function solve2d(f′, xzeros::NTuple{2,@NamedTuple{lt::SVector{2,T},gt::SVector{2,T}}}) where {T}
    D = 2

    # Validating function evaluation
    function f(x::SVector{2,T})
        D = 2
        all(isfinite.(x)) || @show x
        @assert all(isfinite.(x))
        local y = f′(x)::SVector{D,T}
        all(isfinite.(y)) || @show x y
        @assert all(isfinite.(y))
        local h = 1.0e-6
        local y1 = f′(x + SVector(h, 0))::SVector{D,T}
        local y2 = f′(x + SVector(0, h))::SVector{D,T}
        all(y1 .> y) || @show x y 1 (y1-y)/h
        all(y2 .> y) || @show x y 2 (y2-y)/h
        @assert all(y1 .> y)
        @assert all(y2 .> y)
        return y::SVector{D,T}
    end

    domainsize = Inf

    iter = 0
    while true
        iter += 1
        # println("    Iteration $iter")

        # Consistency check
        for val in 1:D
            for ltgt in 1:2
                @show val ltgt
                @show xzeros[val][ltgt]
                @show f(xzeros[val][ltgt])
                @assert abs(f(xzeros[val][ltgt])[val]) <= 10 * eps(T)
            end
            @assert f(xzeros[val].lt)[3 - val] <= 0
            @assert f(xzeros[val].gt)[3 - val] >= 0
        end

        # Find lower and upper end of bracket
        xlo = min.(xzeros[1].lt, xzeros[2].lt)
        xhi = max.(xzeros[1].gt, xzeros[2].gt)
        # Ensure we can find a zero (assuming there is one)
        @assert all(f(xlo) .<= 0 .<= f(xhi))

        # The domain size is the length of the struts times the distance between the two struts
        olddomainsize = domainsize
        domainsize = norm(xhi - xlo, Inf)
        # println("    Domain size: $(round6(domainsize))")
        @assert 0 <= domainsize < olddomainsize

        # Pick a new point in the centre of the bracket
        # println("    Choosing new point:")
        xnew = (xzeros[1].lt + xzeros[2].lt + xzeros[2].gt + xzeros[2].gt) / 4
        # println("        xnew=$(round6.(xnew))   f=$(round6.(f(xnew)))")

        # Check the domain size. If the domain is small enough then we are done.
        if domainsize <= 10*eps(T)
            sol = (x=xnew, xlo=xlo, xhi=xhi, iters=iter)
            # println("    Found solution (domain size small):")
            # println("        x=$(round6.(sol.x))   f=$(round6.(f(sol.x)))")
            # println("        xlo=$(round6.(sol.xlo))   f=$(round6.(f(sol.xlo)))")
            # println("        xhi=$(round6.(sol.xhi))   f=$(round6.(f(sol.xhi)))")
            # println("        iterations=$(sol.iters)")
            return sol
        end

        # Find the two zeros of f along the diagonal
        # xzeronew[val]
        xzeronew = ntuple(D) do val
            xzero = let
                xq(q) = (1-q) * xlo + q * xhi
                fq(q) = f(xq(q))[val]
                q = solve1d(fq, 0.0, 1.0).x
                xq(q)
            end
            @assert abs(f(xzero)[val]) <= 10*eps(T)
            return xzero
        end
        # println("    New strut:")
        for val in 1:D
            # println("        xzeronew[$val]=$(round6.(xzeronew[val]))   f=$(round6.(f(xzeronew[val])))")
        end

        for val in 1:2
            if abs(f(xzeronew[val])[3 - val]) <= 10 * eps(T)
                # Found solution
                sol = (x=xzeronew[val], xlo=xzeronew[val], xhi=xzeronew[val], iters=iter)
                # println("    Found solution (function value small):")
                # println("        x=$(round6.(sol.x))   f=$(round6.(f(sol.x)))")
                # println("        xlo=$(round6.(sol.xlo))   f=$(round6.(f(sol.xlo)))")
                # println("        xhi=$(round6.(sol.xhi))   f=$(round6.(f(sol.xhi)))")
                # println("        iterations=$(sol.iters)")
                return sol
            end
        end

        @assert (f(xzeronew[1])[2] <= 0 && f(xzeronew[2])[1] >= 0) || (f(xzeronew[1])[2] >= 0 && f(xzeronew[2])[1] <= 0)
        for val in 1:2
            if f(xzeronew[val])[3 - val] <= 0
                # # println("        updating xzeros[$val].lt")
                # xzeros[val].lt = xzeronew[val]
                xzeros = setindex(xzeros, setindex(xzeros[val], xzeronew[val], :lt), val)
            elseif f(xzeronew[val])[3 - val] >= 0
                # # println("        updating xzeros[$val].gt")
                # xzeros[val].gt = xzeronew[val]
                xzeros = setindex(xzeros, setindex(xzeros[val], xzeronew[val], :gt), val)
            else
                # impossible
                @assert false
            end
        end
    end
end
