export solve1d
function solve1d(f, xlo::T, xhi::T) where {T}
    @assert xlo <= xhi
    ylo::T = f(xlo)
    yhi::T = f(xhi)
    @assert ylo <= 0
    @assert yhi >= 0

    use_secant = true
    iter = 0
    while true
        iter += 1

        xavg = (xlo + xhi) / 2
        yavg = (ylo + yhi) / 2
        Δx = xhi - xlo
        Δy = yhi - ylo
        abs(Δx) < 10 * max(eps(xlo), eps(xhi), eps(T)) && return (x=xavg, y=yavg, xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, iters=iter)

        x::T = if use_secant
            # secant method
            xavg - Δx / Δy * yavg
        else
            # bisect
            xavg
        end
        x = clamp(x, xlo, xhi)

        y::T = f(x)
        y == 0 && return (x=x, y=y, xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, iters=iter)

        oldsize = abs(xhi - xlo)
        if y < 0
            xlo = x
            ylo = y
        else
            xhi = x
            yhi = y
        end
        newsize = abs(xhi - xlo)

        # Ensure progress
        if newsize >= oldsize
            if use_secant
                use_secant = false
            else
                # lack of progress
                @assert false
            end
        else
            use_secant = true
        end
    end
end

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
