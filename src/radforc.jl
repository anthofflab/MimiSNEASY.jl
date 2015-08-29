module radforccomponent
using Mimi

@defcomp radforc begin
    deltat = Parameter()
    atmco2 = Parameter(index=[time])
    other_forcing = Parameter(index=[time])
    aerosol_forcing = Parameter(index=[time])
    alpha = Parameter()

    rf = Variable(index=[time])
end

function timestep(s::radforc, t::Int)
    v = s.Variables
    p = s.Parameters
    v.rf[t] = 3.7 * log(p.atmco2[t]/p.atmco2[1])/log(2.0) + p.other_forcing[t] + p.alpha * p.aerosol_forcing[t]
end

end
