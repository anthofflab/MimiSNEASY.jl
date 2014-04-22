module radforccomponent
using IAMF

@defcomp radforc begin
    deltat = Parameter()
    atmco2 = Parameter(index=[time])
    other_forcing = Parameter(index=[time])

    rf = Variable(index=[time])
end

function init(s::radforc)    
end

function timestep(s::radforc, t::Int)
    v = s.Variables
    p = s.Parameters
    v.rf[t] = 3.7 * log(p.atmco2[t]/p.atmco2[1])/log(2.0) + p.other_forcing[t]
end

end