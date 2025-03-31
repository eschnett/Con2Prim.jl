export solve1d
function solve1d(f, xlo::T, xhi::T) where {T}
    ylo::T = f(xlo)
    yhi::T = f(xhi)
    @assert ylo <= 0
    @assert yhi >= 0

    iter = 0
    while true
        iter += 1

        xavg = (xlo + xhi) / 2
        yavg = (ylo + yhi) / 2
        Δx = xhi - xlo
        Δy = yhi - ylo
        abs(Δx) < 10 * max(eps(xlo), eps(xhi), eps(T)) && return (x=xavg, y=yavg, xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, iters=iter)

        # bisect
        # x::T = xavg
        # secant method
        x::T = xavg - Δx / Δy * yavg

        y::T = f(x)
        y == 0 && return (x=x, y=y, xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, iters=iter)

        if y < 0
            xlo = x
            ylo = y
        else
            xhi = x
            yhi = y
        end
    end
end
