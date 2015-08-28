using DataFrames
# loadobs.R
#
# Nathan M. Urban (nurban@princeton.edu)
# Woodrow Wilson School, Princeton
#
# Loads observational data.
#
# Also defines index variables to align data and model output
# vectors.

# HADCRUT3 annual global mean surface temperature
# http://www.cru.uea.ac.uk/cru/data/temperature/hadcrut3gl.txt
function loaddata()
  dat = readtable("data/temp.txt", separator=' ')
  dat[:temperature] = dat[:temperature]-mean(dat[1:20,:temperature])
  rename!(dat, :temperature, :obs_temperature)
  df = dat

  # annual global ocean heat content (0-3000 m)
  # Gouretski & Koltermann, "How much is the ocean really warming?",
  # Geophysical Research Letters 34, L01610 (2007).
  dat = readtable("data/gouretski_ocean_heat_3000m.txt", header=false, separator=' ', names=[:year, :obs_ocheat, :obs_ocheatsigma], commentmark='%', allowcomments=true)
  df = join(df, dat, on=:year, kind=:outer)

  # Mauna Loa instrumental CO2
  dat = readtable("data/co2instobs.txt", separator=' ')
  rename!(dat, :co2, :obs_co2inst)
  df = join(df, dat, on=:year, kind=:outer)

  # Law Dome ice core CO2
  dat = readtable("data/co2iceobs.txt", separator=' ', header=false, names=[:year, :obs_co2ice, :unknown])
  dat = DataFrame(year=dat[:year], obs_co2ice=dat[:obs_co2ice])
  df = join(df, dat, on=:year, kind=:outer)
  # obs.co2ice.err = rep(4, length(obs.co2ice))

  # decadal ocean carbon fluxes, McNeil et al. (2003)
  dat = DataFrame(year=[1985, 1995], obs_ocflux=[-1.6, -2.0], obs_ocflux_err=[0.4, 0.4])
  df = join(df, dat, on=:year, kind=:outer)

  # MOC strength, Bryden et al. (2005)
  #dat = read.table("../data/bryden_moc_strength.txt")
  #obs.moc.bryden = dat[,2]
  #obs.moc.bryden.err = rep(6, rep(length(obs.moc.bryden)))
  #obs.moc.bryden.time = dat[,1]
  # MOC strength, Lumpkin & Speer (2007)
  #dat = read.table("../data/lumpkin_speer_moc_strength.txt")
  #obs.moc.ls = dat[,2]
  #obs.moc.ls.err = 2.5
  #obs.moc.ls.time = dat[,1]
  # MOC strength, Cunningham et al. (2007)
  #dat = read.table("../data/cunningham_moc_strength.txt")
  #obs.moc.cunn = dat[,2]
  #obs.moc.cunn.err = 2
  #obs.moc.cunn.time = dat[,1]

  #obs.moc = c(obs.moc.bryden, obs.moc.ls, obs.moc.cunn)
  #obs.moc.err = c(obs.moc.bryden.err, obs.moc.ls.err, obs.moc.cunn.err)
  #obs.moc.time = c(obs.moc.bryden.time, obs.moc.ls.time, obs.moc.cunn.time)

  dat = readtable("data/RCP85_EMISSIONS.csv")
  dat = DataFrame(year=dat[:YEARS], forcing_emis_co2=dat[:FossilCO2]+dat[:OtherCO2])
  df = join(df, dat, on=:year, kind=:outer)

  dat = readtable("data/forcing_rcp85.txt", separator=' ')
  dat = DataFrame(year=dat[:year], forcing_rf_nonco2=dat[:ghg_nonco2] + dat[:solar] + dat[:volcanic] + dat[:other], forcing_rf_aerolsols=dat[:aerosol_direct]+dat[:aerosol_indirect])
  df = join(df, dat, on=:year, kind=:outer)

  df = sort!(df, cols=[:year])
  return df
end
