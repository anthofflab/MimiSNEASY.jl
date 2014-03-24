module radforc
using Iamf

export radforcpar, radforcvar

type radforcpar <: ComponentParameters
    nsteps::Int
    deltat::Float64
    atmco2::Vector{Float64}
    other_forcing::Vector{Float64}

    function radforcpar(nsteps,deltat)
    	p = new()
    	p.nsteps=nsteps
    	p.deltat=deltat
    	return p
    end
end

type radforcvar <: ComponentVariables
	rf::Vector{Float64}

    function radforcvar(p::radforcpar)
        vars = new(zeros(p.nsteps))
        return vars
    end
end

end

importall radforc

function init(p::radforcpar, v::radforcvar)
    println("init radforc")
end

function timestep(p::radforcpar, v::radforcvar, t::Int)
    v.rf[t] = 3.7 * log(p.atmco2[t]/p.atmco2[1])/log(2.0) + p.other_forcing[t]
end