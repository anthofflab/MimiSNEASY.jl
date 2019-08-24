@defcomp radiativeforcing begin
    deltat = Parameter()
    rf_co2 = Parameter(index=[time])
    rf_aerosol = Parameter(index=[time])
    #rf_ch4 = Parameter(index=[time])
    rf_other = Parameter(index=[time])
    alpha = Parameter()
    rf = Variable(index=[time])

    function run_timestep(p, v, d, t)
        v.rf[t] = p.rf_co2[t]  + p.rf_other[t] + p.alpha * p.rf_aerosol[t]    
    end
    
end
