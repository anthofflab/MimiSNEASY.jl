module radforc

type radforcpar
    nsteps::Int
    atmco2::Vector{Float64}
    other_forcing::Vector{Float64}
end

type radforcvar
	rf::Vector{Float64}

    function ccmvar(p::radforcpar)
        vars = new(zeros(p.nsteps))
        return vars
    end
end

function init(p::radforcpar, v::radforcvar)
end

function timestep(p::ccmpar, v::ccmvar, t::Int)
	v.rf[t] = 3.7 * log(p.atmco2[t]/p.atmco2[1])/log(2.0) + p.other_forcing[t]
end

end