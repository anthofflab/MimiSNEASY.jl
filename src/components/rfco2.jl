@defcomp rfco2 begin
 	a₁      = Parameter(default=-2.4e-7)             # CO₂ forcing constant.
    b₁      = Parameter(default=7.2e-4)             # CO₂ forcing constant.
    c₁      = Parameter(default=-2.1e-4)             # CO₂ forcing constant.
    CO₂_0   = Parameter(default=278.05158)             # Initial (pre-industrial) CO₂ concentration (ppm).
	N₂O_0   = Parameter(default=272.95961)             # Initial (pre-industrial) N₂O concentration (ppb).
    N₂O     = Parameter(index=[time]) # N₂O concentration (ppb).

    scale_CO₂ = Parameter(default=1.) #Scaling factor to ensure "Tune the coefficient of CO2 forcing to acheive desired F2x, using # pre-industrial CO2 and N2O FOLLOWING FAIR APPROACH
    CO₂ = Parameter(index=[time])
    rf_co2 = Variable(index=[time])

    function run_timestep(p, v, d, t)
        #CO₂ Radiative Forcing
        CO₂_diff = p.CO₂[t]-p.CO₂_0
        N_hat = 0.5 * (p.N₂O[t] + p.N₂O_0)
    
        v.rf_co2[t] = ((p.a₁*CO₂_diff^2 + p.b₁*abs(CO₂_diff) + p.c₁*N_hat + 5.36) * log(p.CO₂[t] / p.CO₂_0)) * p.scale_CO₂
        
    end    
end
