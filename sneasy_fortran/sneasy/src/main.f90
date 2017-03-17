!    SNEASY model:  Composite DOECLIM, Carbon Cycle, and MOC Boxmodel
!
!    Copyright (C) 2010  Klaus Keller, Nathan Urban, and Brian Tuttle
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!    Klaus Keller, klaus@psu.edu
!    Nathan Urban, nurban@psu.edu
!    Brian Tuttle, btuttle@psu.edu
!
!===========================================================================

PROGRAM main

    USE global
    USE sneasy
    USE CCM
    USE doeclim

    implicit none

    integer(i4b) :: stepsize    = 1
    real(DP) :: Clim_sens       = 2.0d0 !3.4d0     ! Cs [K]
    real(DP) :: Oc_vert_diff    = 1.1d0 !0.55d0    ! kappa2 [cm^2/s]
    real(DP) :: Soil_resp       = 1.311d0   ! Q10
    real(DP) :: CO2_fert        = 0.502d0   ! Beta
    real(DP) :: TC_diff         = 17.722d0  ! Eta
    real(DP) :: HS_NA_surf      = 0.047d0   ! hysens2 [Sv/deg.C]

    character(len=30) :: emis_fname = 'emis_data_sep09.txt'
    integer(i4b), parameter :: emisIOU = 10
    integer(i4b), parameter :: ndata   = 566
    character(len=30) :: anom_fname = 'anomtable.txt'
    integer(i4b), parameter :: anomIOU = 30
    integer(i4b), parameter :: nanom   = 16000
    character(len=30) :: out_fname = 'earth_system.out'
    integer(i4b), parameter :: outIOU  = 40

    real(DP), dimension(:), allocatable :: CO2_emissions
    real(DP), dimension(:), allocatable :: MOC_str              ! [Sv]
    real(DP), dimension(:), allocatable :: rad_forcing          ! [W/m^2]
    real(DP), dimension(:), allocatable :: atmos_CO2            ! [ppm]
    real(DP), dimension(:), allocatable :: atmos_oc_flux        ! [W/m^2]
    real(DP), dimension(:), allocatable :: global_surf_temp     ! [K]
    real(DP), dimension(:), allocatable :: global_oc_heat       ! [10^22 J]
!    real(SP), dimension(:,:), allocatable :: anomtable
    integer(i4b), dimension(:), allocatable :: year
    character(len=100) :: line
    integer(i4b) :: ntempobs, astat, i, j
    real :: readyr
    real(DP) :: readCO2

! Random number generator has not been seeded yet.
    seeded = .false.

    nsteps = ndata
    deltat = real(stepsize,DP)
    ATsize = (/ 100, nanom /)

    call init_SNES_arrays()

! Load ocean anomaly table.
    
    open(UNIT=anomIOU, FILE=trim(anom_fname), STATUS='OLD')
    do i = 1,nanom
        read(anomIOU,*), anomtable(1:100,i)
    end do
    close(anomIOU)

! Allocate arrays for the main program.
    call alloc_arrays(ndata)

! Load CO2 emissions data.
    CO2_emissions = 0.0d0
    open(UNIT=emisIOU, FILE=trim(emis_fname), STATUS='OLD')
    read(emisIOU,*), readyr, readCO2
    CO2_emissions(1) = readCO2
    year(1) = nint(readyr)
    do i = 2,ndata
      do j = 1,stepsize
        read(emisIOU,*), readyr, readCO2
      end do
      CO2_emissions(i) = readCO2
      year(i) = nint(readyr)
    end do
    close(emisIOU)

! Define the model parameters.
    call init_SNES_params(Clim_sens, Oc_vert_diff, Soil_resp, CO2_fert, &
                          TC_diff, HS_NA_surf)

! Run the model.
    do i = 1,nsteps
        call SNES_step(i, CO2_emissions(i), MOC_str(i), rad_forcing(i), &
                        global_surf_temp(i))
    end do

! Write results.
    open(UNIT=outIOU, FILE=trim(out_fname), STATUS='REPLACE', ACTION='WRITE')

    do i = 1,nsteps
        write(outIOU,"(i6,6f12.3)"), year(i), MOC_str(i), rad_forcing(i), &
         co2(i), ocflux(i), global_surf_temp(i), heat_interior(i)
    end do
    close(outIOU)

! Clean up.
    call dealloc_SNES()
    call dealloc_arrays()

CONTAINS
!-----------------------------------------------------------------------------
SUBROUTINE alloc_arrays(n)

    implicit none

    integer(i4b), intent(IN) :: n
    integer(i4b) :: astat = 0

    allocate(CO2_emissions(n), STAT=astat)
    allocate(MOC_str(n), STAT=astat)
    allocate(rad_forcing(n), STAT=astat)
    allocate(atmos_CO2(n), STAT=astat)
    allocate(atmos_oc_flux(n), STAT=astat)
    allocate(global_surf_temp(n), STAT=astat)
    allocate(global_oc_heat(n), STAT=astat)
    allocate(year(n), STAT=astat)

    RETURN

END SUBROUTINE alloc_arrays
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
SUBROUTINE dealloc_arrays()

    implicit none

    integer(i4b) :: astat = 0

    deallocate(CO2_emissions, STAT=astat)
    deallocate(MOC_str, STAT=astat)
    deallocate(rad_forcing, STAT=astat)
    deallocate(atmos_CO2, STAT=astat)
    deallocate(atmos_oc_flux, STAT=astat)
    deallocate(global_surf_temp, STAT=astat)
    deallocate(global_oc_heat, STAT=astat)
    deallocate(anomtable, STAT=astat)
    deallocate(year, STAT=astat)

    RETURN

END SUBROUTINE dealloc_arrays
!-----------------------------------------------------------------------------
END PROGRAM main
