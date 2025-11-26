module Atmospheres
    implicit none
    integer, parameter :: atm_settings_unit = 710
    integer, parameter :: atm_file_unit = 790
    character(len=7), parameter :: test_atm_files(4) = ['HAUS.00','VIRA.11','VIRA.14','VIRA.A6']
    ! ZA - height in km, PA - pressure in atm, WA - log(P) in mbar, ROCO2/ROH2O/ROSO2 - concentration in particles/(cm^2*km)
    real, allocatable :: ZA(:), PA(:), WA(:), TA(:), ROCO2(:), ROH2O(:), ROSO2(:)
    character(len=3) :: gas(3)
    integer :: N_ATM
    integer :: NG
    integer :: JAM ! number of atmospheric levels in an atmospheric file
    integer :: ICO2, IH2O, ISO2

contains

    subroutine initialize_atm_settings()  
        implicit none
        integer :: i

        ICO2=0; IH2O=0; ISO2=0
        open(atm_settings_unit, file='atm_settings.ini')
        read(atm_settings_unit, *) N_ATM
        read(atm_settings_unit, *) NG, JAM
        do i = 1, NG
            read(atm_settings_unit, '(A3)') gas(i)
            select case (gas(i))
                case ('CO2'); ICO2 = i
                case ('H2O'); IH2O = i
                case ('SO2'); ISO2 = i
            end select
        end do
        close(atm_settings_unit)
    end subroutine initialize_atm_settings

    subroutine read_atm_profile(atmo) ! legacy: NEW_ATM
        implicit none
        character(len=7), intent(in) :: atmo
        integer :: k, l ! loop variables
        real :: RORO(3)
        
        allocate(ZA(JAM), PA(JAM), WA(JAM), TA(JAM), ROCO2(JAM), ROH2O(JAM), ROSO2(JAM))
        ROCO2=0.; ROH2O=0.; ROSO2=0.
        ! TODO: estimate if it is proper to have atm_file_unit the same value for all atmospheric files
        open(atm_file_unit, file='data/atmospheres/'//atmo)

        do k = 1, JAM
            read(atm_file_unit, *) ZA(k), PA(k), TA(k), (RORO(l), l=1,NG) 
            WA(k) = alog(PA(k)*1013.16)
            ! TODO: remove these ifs
            IF(ICO2>0) ROCO2(k) = RORO(ICO2)
            IF(IH2O>0) ROH2O(k) = RORO(IH2O)
            IF(ISO2>0) ROSO2(k) = RORO(ISO2)
        end do
        close(atm_file_unit)

    end subroutine read_atm_profile

end module Atmospheres
