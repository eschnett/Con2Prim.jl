function solve1d(newpos, choose, xlo::T, xhi::T) where {T}
    iter = 0
    while true
        iter += 1

        x::T = newpos(xlo, xhi)
        Δx = xhi - xlo
        abs(Δx) < 10 * max(eps(xlo), eps(xhi), eps(T)) && return (x=x, xlo=xlo, xhi=xhi, iters=iter)

        ch = choose(x)
        ch == 0 && return (x=x, xlo=x, xhi=x, iters=iter)

        if ch < 0
            xlo = x
        else
            xhi = x
        end
    end
end

export solve2d
function solve2d(f, xlo::SVector{2,T}, xhi::SVector{2,T}) where {T}
    perm = let
        @assert all(xlo .<= xhi)
        flolo = f(SVector(xlo[1], xlo[2]))
        flohi = f(SVector(xlo[1], xhi[2]))
        fhilo = f(SVector(xhi[1], xlo[2]))
        fhihi = f(SVector(xhi[1], xhi[2]))
        negposlohi = flohi[1] <= 0 && flohi[2] >= 0
        posneglohi = flohi[1] >= 0 && flohi[2] <= 0
        negposhilo = fhilo[1] <= 0 && fhilo[2] >= 0
        posneghilo = fhilo[1] >= 0 && fhilo[2] <= 0
        @assert all(flolo .<= flohi .<= fhihi)
        @assert all(flolo .<= fhilo .<= fhihi)
        @assert (negposlohi ⊻ posneglohi)
        @assert (negposhilo ⊻ posneghilo)
        @assert (negposlohi ⊻ negposhilo)
        if negposlohi
            SVector(1, 2)
        else
            SVector(2, 1)
        end
    end

    function invariant()
        local flolo = f(SVector(xlo[1], xlo[2]))[perm]
        local flohi = f(SVector(xlo[1], xhi[2]))[perm]
        local fhilo = f(SVector(xhi[1], xlo[2]))[perm]
        local fhihi = f(SVector(xhi[1], xhi[2]))[perm]
        return (
            all(xlo .<= xhi) &&
            all(flolo .<= flohi .<= fhihi) &&
            all(flolo .<= fhilo .<= fhihi) &&
            flohi[1] <= 0 &&
            flohi[2] >= 0 &&
            fhilo[1] >= 0 &&
            fhilo[2] <= 0
        )
    end

    @assert invariant()

    iter = 0
    while true
        iter += 1

        # Find a bisection point along the (xlolo, xhihi) diagonal
        x = let
            xavg = (xlo + xhi) / 2
            xrad = (xhi - xlo) / 2
            xq(q) = xavg + q * xrad
            function newpos(qlo, qhi)
                # TODO: Use secant method to improve guess
                return (qlo + qhi) / 2
            end
            function choose(q)
                fq = f(xq(q))[perm]
                if (fq[1] <= 0 && fq[2] >= 0) || (fq[1] >= 0 && fq[2] <= 0)
                    return 0
                elseif all(fq .<= 0)
                    return +1
                elseif all(fq .>= 0)
                    return -1
                else
                    @assert false
                end
            end
            sol = solve1d(newpos, choose, -1.0, 1.0)
            xq(sol.x)
        end

        Δx = xhi - xlo
        norm(Δx, Inf) < 10 * max(eps.(xlo)..., eps.(xhi)..., eps(T)) && return (x=x, xlo=xlo, xhi=xhi, iters=iter)

        # Choose a quadrant to keep
        let
            fx = f(x)[perm]
            (fx[1] == 0 && fx[2] == 0) && return (x=x, xlo=xlo, xhi=xhi, iters=iter)
            negposx = fx[1] <= 0 && fx[2] >= 0
            posnegx = fx[1] >= 0 && fx[2] <= 0
            @assert negposx ⊻ posnegx
            if negposx
                # Move lohi corner
                xlo, xhi = SVector(x[1], xlo[2]), SVector(xhi[1], x[2])
            else
                # Move hilo corner
                xlo, xhi = SVector(xlo[1], x[2]), SVector(x[1], xhi[2])
            end
        end

        @assert invariant()
    end
end
