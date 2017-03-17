close all
clear all

emis = load('emis_data_sep09.txt') ;

ES = load('earth_system.out') ;
yrs_emis = emis(:,1) ;
emissions = emis(:,2) ;
years = ES(:,1);    
MOC = ES(:,2);      
rad_forc = ES(:,3); 
co2 = ES(:,4);      
ocflux = ES(:,5);   
temp = ES(:,6); 
ocheat = ES(:,7); 

mytallfigure(1,2) ;

subplot(4,1,1)
plot(yrs_emis,emissions)
ylabel('CO_2 Emis. [GtC]')
title({'Cs = 3.4 deg.C, Q10 = 1.311, Beta = 0.502,';'Eta = 17.722, Hyd.sens. = 0.047 Sv/deg.C'})

subplot(4,1,2)
plot(years,MOC)
ylabel('MOC str. [Sv]')

subplot(4,1,3)
plot(years,rad_forc)
ylabel('Rad. forcing [W/m^2]')

subplot(4,1,4)
plot(years,co2)
ylabel('Atm. CO_2 [ppm]')
xlabel('year')

stampit

mytallfigure(2,2) ;

subplot(3,1,1)
plot(years,ocflux)
ylabel('Atm./ocean flux [GtC]')

subplot(3,1,2)
plot(years,temp)
ylabel('Global surf. temp. [K]')

subplot(3,1,3)
plot(years,ocheat)
ylabel('Global ocean heat [10^{22} J]')
xlabel('year')

stampit

!/bin/rm ES_figs.ps
plotallfigures(gcf,'ES_figs.ps') ;
