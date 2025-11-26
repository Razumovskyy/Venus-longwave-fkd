program main
    use Atmospheres
    use KD
    implicit none
    integer :: j
    character(len=7) :: atmo_name
    
    call initialize_atm_settings() ! reads N_ATM, NG, JAM, GAS(:)

    do j = 1, N_ATM
        
        atmo_name = test_atm_files(j)
        call read_atm_profile(atmo_name) ! allocates ZA, PA, TA, ROCO2, ROH2O, ROSO2
        call read_kd_data() ! reads reference log-P grid, and three effective cross-section files: CO2, H2O, SO2
        call interpolate_temperature() ! interpolation initial temperature profile onto reference log-P grid
        call temperature_correction() ! temperature correction for 7 effective cross-sections for CO2
        call kd_calc()
        call planck_calc()
        deallocate(ZA, PA, WA, TA, ROCO2, ROH2O, ROSO2)
        deallocate(VAC, ACO2, AH2O, ASO2)
    
    end do

end program main
