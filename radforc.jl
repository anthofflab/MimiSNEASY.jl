module radforccomponent
using IAMF

@defcomp radforc begin
    addParameter(deltat, Float64)
    addParameter(atmco2, Float64, index=[time])
    addParameter(other_forcing, Float64, index=[time])

    addVariable(rf, Float64, index=[time])
end

function init(s::radforc)    
end

function timestep(s::radforc, t::Int)
    v = s.Variables
    p = s.Parameters
    v.rf[t] = 3.7 * log(p.atmco2[t]/p.atmco2[1])/log(2.0) + p.other_forcing[t]
end

end