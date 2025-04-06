chop(x::Real) = ifelse(abs(x) < eps(typeof(x))^(3/4), zero(x), x)
round6(x::Real) = round(x;digits=6)

export solve2d
function solve2d(f′, xlo::SVector{2,T}, xhi::SVector{2,T}) where {T}
    # Check inputs
    @assert all(isfinite.(xlo) .&& isfinite.(xhi))
    @assert all(xlo .<= xhi)

    # Validating function evaluation
    function f(x)
        all(isfinite.(x)) || @show x
        @assert all(isfinite.(x))
        local y = f′(x)::SVector{2,T}
        all(isfinite.(y)) || @show x y
        @assert all(isfinite.(y))
        return y
    end

    # Ensure we can find a zero (assuming there is one)
    @assert all(f(xlo) .<= 0 .<= f(xhi))

    xbnd = (xlo, xhi)

    # Find the zeros of f on the boundary of the domain
    # xzeros[val][ltgt]
    xzeros = Tuple(let
        xzeros = SVector{2,T}[]
        for dir in 1:2
            for face in 1:2
                xmin = setindex(xbnd[face], xbnd[1][dir], dir)
                xmax = setindex(xbnd[face], xbnd[2][dir], dir)
                fmin = f(xmin)[val]
                fmax = f(xmax)[val]
                if fmin <= 0 && fmax >= 0 || fmin >= 0 && fmax <= 0
                    xzero = let
                        xq(q) = setindex(xbnd[face], q, dir)
                        fq(q) = f(xq(q))[val]
                        q = solve1d(fq, xbnd[1][dir], xbnd[2][dir]).x
                        xq(q)
                    end
                    push!(xzeros, xzero)
                end
            end
        end
        # There should be 2 zeros, but we might have more if we double-counted a zero in a corner
        @assert length(xzeros) >= 2
        xzerolt = xzeros[findfirst(x -> f(x)[3 - val] <= 0, xzeros)]
        xzerogt = xzeros[findfirst(x -> f(x)[3 - val] >= 0, xzeros)]
        (lt=xzerolt, gt=xzerogt)
    end for val in 1:2)

    domainsize = SVector{2,T}(Inf,Inf)

    iter = 0
    while true

        # Consistency check
        for val in 1:2
            for ltgt in 1:2
                @assert abs(f(xzeros[val][ltgt])[val]) <= 10*eps(T)
            end
            @assert f(xzeros[val].lt)[3-val] <= 0
            @assert f(xzeros[val].gt)[3-val] >= 0
        end
        xlo = min.(xzeros[1].lt, xzeros[2].lt)::SVector{2,T}
        xhi = max.(xzeros[1].gt, xzeros[2].gt)::SVector{2,T}
        flo = f(xlo)
        fhi = f(xhi)
        @assert all(flo .<= 0) && all(fhi .>= 0)
        olddomainsize = domainsize
        domainsize = xhi - xlo

        @assert all(0 .<= domainsize .<= olddomainsize) && any(domainsize .< olddomainsize)


        # Check the domain size. If the domain is small enough then we are done.
        if norm(xhi - xlo, Inf) <= 10 * max(eps.(xlo)..., eps.(xhi)..., eps(T))
            return (x=(xlo + xhi) / 2, xlo=xlo, xhi=xhi, iters=iter)
        end

        iter += 1
        println("iter=$iter")
        println("    Current state:")
        println("        xzeros[1].lt=$(round6.(xzeros[1].lt))   f=$(round6.(f(xzeros[1].lt)))")
        println("        xzeros[1].gt=$(round6.(xzeros[1].gt))   f=$(round6.(f(xzeros[1].gt)))")
        println("        xzeros[2].lt=$(round6.(xzeros[2].lt))   f=$(round6.(f(xzeros[2].lt)))")
        println("        xzeros[2].gt=$(round6.(xzeros[2].gt))   f=$(round6.(f(xzeros[2].gt)))")
        println("        xlo=$(round6.(xlo))   f=$(round6.(f(xlo)))")
        println("        xhi=$(round6.(xhi))   f=$(round6.(f(xhi)))")

        # Pick a new point in the centre of the bracket
        x = (xzeros[1].lt + xzeros[1].gt + xzeros[2].lt + xzeros[2].gt) / 4
        fx = f(x)
        if all(fx .== 0)
            return (x=x, xlo=x, xhi=x, iters=iter)
        end

        # Find the two zeros of f along the diagonal
        # xzeronew[val]
        xzeronew = Tuple(let
            xq(q) = (1 - q) * xlo + q * xhi
            fq(q) = f(xq(q))[val]
            q = solve1d(fq, 0.0, 1.0).x
            xq(q)
        end for val in 1:2)
        println("    New point:")
        println("        xzeronew[1]=$(round6.(xzeronew[1]))   f=$(round6.(f(xzeronew[1])))")
        println("        xzeronew[2]=$(round6.(xzeronew[2]))   f=$(round6.(f(xzeronew[2])))")

        for val in 1:2
            if f(xzeronew[val])[3 - val] == 0
                # Found solution
                return (x=xzeronew[val], xlo=xzeronew[val], xhi=xzeronew[val], iters=iter)
            end
        end
        @assert (f(xzeronew[1])[2] <= 0 && f(xzeronew[2])[1] >= 0) || (f(xzeronew[1])[2] >= 0 && f(xzeronew[2])[1] <= 0)
        for val in 1:2
            if f(xzeronew[val])[3 - val] <= 0
                # xzeros[val].lt = xzeronew[val]
                println("        updating xzeros[$val].lt")
                xzeros = setindex(xzeros, setindex(xzeros[val], xzeronew[val], :lt), val)
            elseif f(xzeronew[val])[3 - val] >= 0
                # xzeros[val].gt = xzeronew[val]
                println("        updating xzeros[$val].gt")
                xzeros = setindex(xzeros, setindex(xzeros[val], xzeronew[val], :gt), val)
            else
                # impossible
                @assert false
            end
        end
    end
end
