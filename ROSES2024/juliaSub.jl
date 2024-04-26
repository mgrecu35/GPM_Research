function bisectm(xvec, nv, r)
    n1 = 1
    n2 = nv
    if r < xvec[1]
        return 1
    end
    if r >= xvec[n2]
        return nv
    end
    nmid = div(n1 + n2, 2)
    it = 0
    while !(r >= xvec[nmid] && r < xvec[nmid + 1]) && it < 7
        it += 1
        if r > xvec[nmid]
            n1 = nmid
        else
            n2 = nmid
        end
        nmid = div(n1 + n2, 2)
    end
    return nmid
end

function hb(zKum, alpha, beta, dr, srt_piaKu)
    q = 0.2 * log(10)
    zeta = q * beta * alpha * 10.0 .^(0.1 * zKum * beta) * dr
    zetamax = 1.0 - 10.0^(-srt_piaKu / 10 * beta)
    
    if sum(zeta) > zetamax
        eps = 0.9999 * zetamax / sum(zeta)
    else
        eps = 1.0
    end
    corrc = eps * cumsum(zeta)
    dpia= -10.0 ./ beta * log10.(1.0 .- corrc)
    zc = zKum .+ dpia
    return zc, eps, -10 / beta * log10(1 - corrc[end])
end