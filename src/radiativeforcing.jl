module radiativeforcingcomponent
using Mimi

@defcomp radiativeforcing begin
    deltat = Parameter()
    rf_co2 = Parameter(index=[time])
    rf_aerosol = Parameter(index=[time])
    #rf_ch4 = Parameter(index=[time])
    rf_other = Parameter(index=[time])
    alpha = Parameter()
    rf = Variable(index=[time])
end

function run_timestep(s::radiativeforcing, t::Int)
    v = s.Variables
    p = s.Parameters

v.rf[t] = p.rf_co2[t]  + p.rf_other[t] + p.alpha * p.rf_aerosol[t]

end

end