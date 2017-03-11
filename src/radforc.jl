module radforccomponent
using Mimi

@defcomp radforc begin
    deltat = Parameter()
    atmco2 = Parameter(index=[time])
    #other_forcing = Parameter(index=[time])
    #aerosol_forcing = Parameter(index=[time])
    forcing_other = Parameter(index=[time])
    forcing_volcanic = Parameter(index=[time])
    forcing_solar = Parameter(index=[time])
    forcing_ghg_nonco2 = Parameter(index=[time])
    forcing_aerosol_direct = Parameter(index=[time])
    forcing_aerosol_indirect = Parameter(index=[time])
	alpha = Parameter()

    #forcing_aer = Variable(index=[time])
    forcing_nonco2  = Variable(index=[time])
    rf = Variable(index=[time])
end

function run_timestep(s::radforc, t::Int)
    v = s.Variables
    p = s.Parameters

    v.forcing_nonco2[t] = p.alpha * (p.forcing_aerosol_direct[t] + p.forcing_aerosol_indirect[t]) + p.forcing_ghg_nonco2[t] + p.forcing_solar[t] + p.forcing_volcanic[t] + p.forcing_other[t]

    v.rf[t] = 3.7 * log(p.atmco2[t]/p.atmco2[1])/log(2.0) + v.forcing_nonco2[t]


end

end