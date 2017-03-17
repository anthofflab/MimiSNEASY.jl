SUBROUTINE init_ccm(tab_ni, tab_nj, oceanomtab)

    USE global
    USE CCM

    implicit none

    integer(i4b), intent(IN) :: tab_ni
    integer(i4b), intent(IN) :: tab_nj
    real(DP), dimension(tab_ni,tab_nj), intent(IN) :: oceanomtab

    ATsize = (/ tab_ni, tab_nj /)   ! Ocean anomaly table size

    call alloc_anomtab()

    anomtable = oceanomtab

    RETURN

END SUBROUTINE init_ccm

SUBROUTINE fin_ccm()

    USE global
    USE CCM

    implicit none

    call dealloc_anomtab()

    RETURN

END SUBROUTINE fin_ccm

SUBROUTINE run_ccm(ns, &
        temp_forcing, CO2_emis_forcing, &
        Climate_sens, Soil_respiration, Carbon_fertilization, Thermocline_diff, init_atm_CO2, &
        atmco2_out)

    USE global
    USE CCM

    implicit none

    integer(i4b), intent(IN) :: ns
    real(DP), intent(IN) :: Climate_sens
    real(DP), intent(IN) :: Soil_respiration
    real(DP), intent(IN) :: Carbon_fertilization
    real(DP), intent(IN) :: Thermocline_diff
    real(DP), intent(IN) :: init_atm_CO2
    real(DP), dimension(ns), intent(IN) :: temp_forcing
    real(DP), dimension(ns), intent(IN) :: CO2_emis_forcing
    real(DP), dimension(ns), intent(OUT) :: atmco2_out

    integer(i4b) :: i

! Assign global variables.
    nsteps = ns
    deltat = 1.0d0

    call init_CCM_arrays()

    call init_CCM_parameters(Climate_sens, Soil_respiration, &
                             Carbon_fertilization, Thermocline_diff)

    atmco2(1) = init_atm_CO2

    do i = 1,nsteps
        call CC_model(i, temp_forcing(i), CO2_emis_forcing(i))
    end do

    atmco2_out = atmco2

    call dealloc_CCM()

    RETURN

END SUBROUTINE run_ccm
