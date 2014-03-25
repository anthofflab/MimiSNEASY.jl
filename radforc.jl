module radforccomponent
using IAMF

type radforcpar
    deltat::Float64
    atmco2::Vector{Float64}
    other_forcing::Vector{Float64}

    function radforcpar(nsteps::Int)
        p = new()
        return p
    end
end

type radforcvar
    rf::Vector{Float64}

    function radforcvar(nsteps::Int)
        vars = new(zeros(nsteps))
        return vars
    end
end

type radforc <: ComponentState
    p::radforcpar
    v::radforcvar
    nsteps::Int

    function radforc(nsteps)
        s = new()
        s.nsteps = nsteps
        s.p = radforcpar(nsteps)
        s.v = radforcvar(nsteps)
        return s
    end
end

import IAMF.init
import IAMF.timestep

function init(s::radforc)    
end

function timestep(s::radforc, t::Int)
    v = s.v
    p = s.p
    v.rf[t] = 3.7 * log(p.atmco2[t]/p.atmco2[1])/log(2.0) + p.other_forcing[t]
end

end