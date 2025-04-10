export solve3d

function triangle_area(xs::NTuple{3,SVector{3,T}}) where {T}
    a = norm(xs[1] - xs[2])
    b = norm(xs[2] - xs[3])
    c = norm(xs[3] - xs[1])
    s = (a+b+c)/2
    # Heron's formula
    A = sqrt(s * (s-a) * (s-b) * (s-c))
    return A
end

function cross(x::SVector{3}, y::SVector{3})
    return SVector(x[2] * y[3] - y[2] * x[3], x[3] * y[1] - y[3] * x[1], 3[1] * y[2] - y[2] * x[1])
end

function getindex2(x::SVector{3}, dir::Integer)
    # Skip direction `dir`
    dir==1 && return SVector(x[2], x[3])
    dir==2 && return SVector(x[1], x[3])
    dir==3 && return SVector(x[1], x[2])
    @assert false
end

function setindex2(x::SVector{3}, q::SVector{2}, dir::Integer)
    # Skip direction `dir`
    dir==1 && return SVector(x[1], q[1], q[2])
    dir==2 && return SVector(q[1], x[2], q[2])
    dir==3 && return SVector(q[1], q[2], x[3])
    @assert false
end

# z = A + B * x + C * y + D * x*y
function make_interp2(xyzs::NTuple{4,NTuple{3,T}}) where {T}
    mat = zeros(T, 4, 4)
    rhs = zeros(T, 4)
    for n in 1:4
        xyz = xyzs[n]
        x, y, z = xyz
        # z = A + B * x + C * y + D * x*y
        mat[n, :] .= [1, x, y, x*y]
        rhs[n] = z
    end
    ABCD = mat \ rhs
    ABCD = Tuple(ABCD)
    # Test
    for n in 1:4
        xyz = xyzs[n]
        x, y, z = xyz
        @assert abs(interp2(ABCD, SVector(x, y)) - z) <= 10*eps(T)
    end
    return ABCD
end
function interp2(ABCD::NTuple{4,T}, xy::SVector{2,T}) where {T}
    A, B, C, D = ABCD
    x, y = xy
    return A + B * x + C * y + D * x * y
end

function solve3d(f′, xlo::SVector{3,T}, xhi::SVector{3,T}) where {T}
    D = 3

    # Check inputs
    @assert all(isfinite.(xlo) .&& isfinite.(xhi))
    @assert all(xlo .<= xhi)

    # Validating function evaluation
    function f(x::SVector{3,T})
        D = 3
        @assert all(isfinite.(x))
        local y = f′(x)::SVector{D,T}
        all(isfinite.(y)) || @show x y
        @assert all(isfinite.(y))
        local h = 1.0e-6
        local y1 = f′(x + SVector(h, 0, 0))::SVector{D,T}
        local y2 = f′(x + SVector(0, h, 0))::SVector{D,T}
        local y3 = f′(x + SVector(0, 0, h))::SVector{D,T}
        all(y1 .> y) || @show x y 1 (y1-y)/h
        all(y2 .> y) || @show x y 2 (y2-y)/h
        all(y3 .> y) || @show x y 3 (y3-y)/h
        @assert all(y1 .> y)
        @assert all(y2 .> y)
        @assert all(y3 .> y)
        return y::SVector{D,T}
    end

    println("solve3d:")
    println("    Inputs:")
    println("        xlo=$(round6.(xlo))   f=$(round6.(f(xlo)))")
    println("        xhi=$(round6.(xhi))   f=$(round6.(f(xhi)))")

    # Ensure we can find a zero (assuming there is one)
    @assert all(f(xlo) .<= 0 .<= f(xhi))

    xbnd = (lo=xlo, hi=xhi)

    # Find six double-zeros on the domain faces
    println("    Finding double zeros of f on faces")
    # xzeros[ltgt][nzval]::SVector{D,T}
    xzeros = ntuple(2) do ltgt
        ntuple(D) do nzval
            for dir in 1:D
                for face in 1:2
                    xq(q) = setindex2(xbnd[face], q, dir)
                    fq(q) = getindex2(f(xq(q)), nzval)
                    qmin = getindex2(xlo, dir)
                    qmax = getindex2(xhi, dir)
                    # xmin = xq(qmin)
                    # xmax = xq(qmax)
                    fmin = fq(qmin)
                    fmax = fq(qmax)
                    if all(fmin .<= 0) && all(fmax .>= 0)
                        q = solve2d(fq, qmin, qmax).x
                        x = xq(q)
                        fx = f(x)
                        if ltgt==1 && fx[nzval]<=0 || ltgt==2 && fx[nzval]>=0
                            return x
                        end
                    end
                end
            end
        end
    end
    xzeros = (lt=xzeros[1], gt=xzeros[2])

    domainsize = Inf

    iter = 0
    while true
        println("    Iteration $iter:")

        println("        Bracketing points:")
        for nzval in 1:D
            for ltgt in (:lt, :gt)
                x = getindex(xzeros, ltgt)[nzval]
                fx = f(x)
                println("            xzeros.$ltgt[$nzval]=$(round6.(x))   f=$(round6.(fx))")
            end
        end

        xmid = sum(xzeros.lt[nzval] + xzeros.gt[nzval] for nzval in 1:D) / 6
        println("        Midpoint:")
        println("            xmid=$(round6.(xmid))   f=$(round6.(f(xmid)))")

        olddomainsize = domainsize
        domainsize = sum(norm(xzeros.gt[nzval] - xzeros.lt[nzval], Inf) for nzval in 1:D)
        println("        Domain size: $domainsize")
        @assert domainsize >= 0
        # Require progress
        @assert domainsize < olddomainsize

        xlo = min.(xzeros.lt[1], xzeros.lt[2], xzeros.lt[3])
        xhi = max.(xzeros.gt[1], xzeros.gt[2], xzeros.gt[3])
        println("        Box:")
        println("            xlo=$(round6.(xlo))   f=$(round6.(f(xlo)))")
        println("            xhi=$(round6.(xhi))   f=$(round6.(f(xhi)))")

        if domainsize <= 10*eps(T)
            sol = (x=xmid, xlo=xlo, xhi=xhi, iters=iter)
            println("    Found solution (domain size small):")
            println("        x=$(round6.(sol.x))    f=$(round6.(f(sol.x)))")
            println("        xlo=$(round6.(sol.xlo))    f=$(round6.(f(sol.xlo)))")
            println("        xhi=$(round6.(sol.xhi))    f=$(round6.(f(sol.xhi)))")
            println("        iters=$(sol.iters)")
            return sol
        end

        for nzval in 1:D
            for ltgt in (:lt, :gt)
                x = getindex(xzeros, ltgt)[nzval]
                fx = f(x)
                if all(abs.(fx) .<= 10*eps(T))
                    sol = (x=x, xlo=xlo, xhi=xhi, iters=iter)
                    println("    Found solution (function value small):")
                    println("        x=$(round6.(sol.x))    f=$(round6.(f(sol.x)))")
                    println("        xlo=$(round6.(sol.xlo))    f=$(round6.(f(sol.xlo)))")
                    println("        xhi=$(round6.(sol.xhi))    f=$(round6.(f(sol.xhi)))")
                    println("        iters=$(sol.iters)")
                    return sol
                end
            end
        end

        iter += 1

        println("        Choosing new bracketing point")
        nzval = argmax(norm(xzeros.gt[nzval] - xzeros.lt[nzval], Inf) for nzval in 1:D)
        println("            Choosing to bisect for f[$nzval]")
        xnew = let
            # For simplicity, choose as solving plane one of the coordinate planes
            xsep = xzeros.gt[nzval] - xzeros.lt[nzval]
            dir = argmax(abs.(xsep))
            dir1, dir2 = getindex2(SVector(1, 2, 3), dir)
            val1, val2 = getindex2(SVector(1, 2, 3), nzval)

            ABCD = make_interp2((
                (xzeros.lt[val1][dir1], xzeros.lt[val1][dir2], xzeros.lt[val1][dir]),
                (xzeros.lt[val2][dir1], xzeros.lt[val2][dir2], xzeros.lt[val2][dir]),
                (xzeros.gt[val1][dir1], xzeros.gt[val1][dir2], xzeros.gt[val1][dir]),
                (xzeros.gt[val2][dir1], xzeros.gt[val2][dir2], xzeros.gt[val2][dir]),
            ))

            qzeros = (
                lt=(getindex2(xzeros.lt[val1], dir), getindex2(xzeros.lt[val2], dir)),
                gt=(getindex2(xzeros.gt[val1], dir), getindex2(xzeros.gt[val2], dir)),
            )
            qzeros = ((lt=qzeros.lt[2], gt=qzeros.gt[2]), (lt=qzeros.lt[1], gt=qzeros.gt[1]))

            function xq(q)
                xdir = interp2(ABCD, q)
                return setindex2(SVector(xdir, xdir, xdir), q, dir)
            end
            fq(q) = getindex2(f(xq(q)), nzval)

            println("            Lower-dimensional control points:")
            println("                qzeros[1].lt=$(round6.(qzeros[1].lt))   f=$(round6.(fq(qzeros[1].lt)))")
            println("                qzeros[1].gt=$(round6.(qzeros[1].gt))   f=$(round6.(fq(qzeros[1].gt)))")
            println("                qzeros[2].lt=$(round6.(qzeros[2].lt))   f=$(round6.(fq(qzeros[2].lt)))")
            println("                qzeros[2].gt=$(round6.(qzeros[2].gt))   f=$(round6.(fq(qzeros[2].gt)))")
            println("                xzeros[1].lt=$(round6.(xq(qzeros[1].lt)))   f=$(round6.(f(xq(qzeros[1].lt))))")
            println("                xzeros[1].gt=$(round6.(xq(qzeros[1].gt)))   f=$(round6.(f(xq(qzeros[1].gt))))")
            println("                xzeros[2].lt=$(round6.(xq(qzeros[2].lt)))   f=$(round6.(f(xq(qzeros[2].lt))))")
            println("                xzeros[2].gt=$(round6.(xq(qzeros[2].gt)))   f=$(round6.(f(xq(qzeros[2].gt))))")

            q = solve2d(fq, qzeros).x
            xq(q)
        end
        fnew = f(xnew)
        println("            Found xnew=$(round6.(xnew))    f=$(round6.(fnew))")
        if fnew[nzval]<=0
            println("            Updating xzeros.lt[$nzval]")
            # xzeros.lt[nzval] = xnew
            xzeros = setindex(xzeros, setindex(xzeros.lt, xnew, nzval), :lt)
        else
            println("            Updating xzeros.gt[$nzval]")
            # xzeros.gt[nzval] = xnew
            xzeros = setindex(xzeros, setindex(xzeros.gt, xnew, nzval), :gt)
        end
    end
end
