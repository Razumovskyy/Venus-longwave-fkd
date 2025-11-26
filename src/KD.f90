module KD
    use Atmospheres
! !       - Exposes: load_w_grid(data_dir, w), load_kd_tables(data_dir, sco2, sh2o, sso2),
!         open_planck_table(data_dir), read_planck_slice(unit, idx, a)
!       - Handles all F_LW files (W_STAND, O.CO2, O.H2O, O.SO2, KD-Planck, WS(T)).
!       - No dependencies.
    implicit none
    integer, parameter :: planck_fileunit = 910
    integer, parameter :: vac_output_fileunit = 520
    integer, parameter :: planck_output_fileunit = 530
    integer, parameter :: NKTR = 32 ! total number of k-terms
    integer, parameter :: MCO2=28, MH2O=6, MSO2=4 ! number of effective cross-sections per gas
    integer, parameter :: JM = 65 ! number of points in log p
    integer, parameter :: J_T=26 ! number of lines in the T-table
    integer, parameter :: I_M=7 ! number of k-terms to be corrected
    integer, parameter :: KT_SO2(32) = (/ &
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
            1,2,3,4, &
            0,0,0,0,0  &
            /)
    integer, parameter :: KT_H2O(32) = (/ &
            1,2,3,4, &
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
            5, &
            0,0,0, &
            6, &
            0 &
            /)
    integer, parameter :: KT_CO2(32) = (/ &
            0,0,0,0, &
            1,2,3,4,5,6,7,8,9,10, &
            11,12,13,14,15,16,17,18,19,20, &
            21,22,23,24,25,26,27,28 &
            /)
    integer, parameter :: w_stand_fileunit = 750
    integer, parameter :: co2_kd_fileunit = 810, h2o_kd_fileunit = 820, so2_kd_fileunit = 830, ws_t_fileunit = 840

    real :: S0(JM, I_M)
    real :: PL(NKTR)
    real :: WT(J_T)
    real :: COR_T(J_T,I_M)
    real, allocatable :: VAC(:,:), ACO2(:,:), AH2O(:,:), ASO2(:,:)
    real :: W(JM) ! array of reference log-P grid points (65 values)
    real :: SCO2(JM,MCO2),SH2O(JM,MH2O),SSO2(JM,MSO2)
    real :: SCO2T(JM,MCO2) ! T-corrected SCO2
    real :: TV(JM) ! temperature profile interpolated onto reference log-p grid

contains

    subroutine read_kd_data()
        ! TODO: consider of moving this subroutine
        implicit none
        integer :: ii, jj
    
        ! TODO: refactor this to avoid redundancies
        open(w_stand_fileunit, file='data/F_LW/W_STAND')
        do ii = 1, JM
            read(w_stand_fileunit, *) W(ii) ! reading the k-grid
        end do
        close(w_stand_fileunit)

        if (ICO2 > 0) then
            open(co2_kd_fileunit, file='data/F_LW/O.CO2')
        do jj = 1, JM
            read(co2_kd_fileunit, *) SCO2(jj,:)
        end do
        close(co2_kd_fileunit)
        end if
        
        if (IH2O > 0) then
            open(h2o_kd_fileunit, file='data/F_LW/O.H2O')
        do jj = 1, JM
            read(h2o_kd_fileunit, *) SH2O(jj,:)
        end do
        close(h2o_kd_fileunit)
        end if

        if (ISO2 > 0) then
            open(so2_kd_fileunit, file='data/F_LW/O.SO2')
        do jj = 1, JM
            read(so2_kd_fileunit, *) SSO2(jj,:)
        end do
        close(so2_kd_fileunit)
        end if

        open(ws_t_fileunit, file='data/F_LW/WS(T)')
        do jj = 1, J_T
            read(ws_t_fileunit, *) WT(jj), COR_T(jj,:)
        end do
        close(ws_t_fileunit)

    end subroutine read_kd_data

    subroutine interpolate_temperature()
        ! TODO: consider join this with temperature_correction subroutine
        implicit none
        integer :: kk, ll ! loop variable
        real :: c1, c2, w1, w2, wwwww
        
        ! interpolation of temperature profile onto reference log-p grid !
        TV = 0.
        do ll = 1, JM ! loop over reference log-p
            wwwww = W(ll)
            do kk = 2, JAM ! loop over atmospheric profile
                w1=WA(kk-1); w2=WA(kk)
                if (w2==0.) exit
                if (w1>=wwwww .and. w2<=wwwww) exit
            end do
            if (kk>JAM) exit ! TODO: remove this if
            ! TODO: consider fully remove c1, c2 (just simple linear interpolation)
            c2 = (w1-wwwww) / (w1-w2)
            c1 = 1. - c2
            TV(ll) = c1*TA(kk-1) + c2*TA(kk)
        end do 

    end subroutine interpolate_temperature

    subroutine temperature_correction()
        implicit none
        real :: S_(JM, I_M) ! corrected cross-section
        real :: RAT
        integer :: i, j
        real, parameter :: TH(JM) = [725.00,709.00,693.00,677.00,661.00,645.00,630.00,616.00,602.00,587.00,570.50,553.50,537.00,520.50,&
            503.50,486.50,469.50,453.50,438.00,423.00,409.50,397.00,385.00,372.50,357.50,340.50,322.00,301.50,&
            282.50,268.50,257.50,248.00,240.00,236.50,235.00,230.50,225.50,220.00,214.00,207.00,198.50,189.00,&
            180.00,173.50,169.00,165.50,163.00,162.00,162.50,166.50,167.50,162.50,157.50,152.50,147.50,142.50,&
            134.50,123.00,121.00,130.00,138.50,144.00,146.50,143.50,137.00 &
        ]
        do i = 1, JM
            S0(i,1:7)=SCO2(i,3:9) ! temperature correction is applied only on 3-9 columns of the reference CO2 k-terms
        end do
        
        S_ = S0
        do j = 36,61
            RAT = (TV(j)-TH(j)) / TH(j)
            do i = 1, I_M ! iteration over each concerned k-term
                S_(j,i) = S0(j,i) + COR_T(j-35,i)*RAT
            end do
        end do

        SCO2T = SCO2
        do J=36,61
            SCO2T(J,3:9) = S_(J,1:7)
        end do
    end subroutine temperature_correction

    subroutine kd_calc()
        implicit none
        integer :: j, jj, kk
        integer :: k1, k2, k3
        real :: wwwww, w1, w2, c1, c2

        allocate (VAC(JAM,NKTR), ACO2(JAM,MCO2), AH2O(JAM,MH2O), ASO2(JAM,MSO2))
        
        do j = 1, JAM
            wwwww = WA(j)
            do jj= 2, JM
                w1=W(JJ-1); w2=W(JJ)
                if (w2 == 0.) exit
                if (w2 <= wwwww) exit
            end do
            c1 = 1.; c2 = 0.
            
            ! TODO: this if seems to me redundant
            if (w1 >= wwwww) then
                c2=(w1-wwwww) / (w1-w2); c1=1.-c2
            end if
            
            if (JJ<=JM) then ! TODO: these ifs also look redundant
                if (ISO2>0) ASO2(j,:) = C1*SSO2(jj-1,:) + C2*SSO2(jj,:) ! *** SO2 *** !
                if (IH2O>0) AH2O(j,:) = C1*SH2O(jj-1,:) + C2*SH2O(jj,:) ! *** H2O *** !
                if (ICO2>0) ACO2(j,:) = C1*SCO2T(jj-1,:) + C2*SCO2T(jj,:) ! *** CO2-T-corrected *** !
            else
                ASO2(J,:)=0. ; AH2O(J,:)=0. ; ACO2(J,:)=0.
            end if
        end do
    
        ! TODO: join with previous loop
        ! account for concentration !
        do j = 1, JAM
            ACO2(j,:) = exp(ACO2(j,:)) * ROCO2(j)
            AH2O(j,:) = exp(AH2O(j,:)) * ROH2O(j)
            ASO2(j,:) = exp(ASO2(j,:)) * ROSO2(j)
        end do

        VAC = 0.
        do kk = 1, NKTR
            k1 = KT_CO2(kk); k2 = KT_H2O(kk); k3 = KT_SO2(kk)
            ! TODO: Probably these if-s can be removed
            if (K1>0) VAC(:,KK) = VAC(:,KK) + ACO2(:,K1)
            if (K2>0) VAC(:,KK) = VAC(:,KK) + AH2O(:,K2)
            if (K3>0) VAC(:,KK) = VAC(:,KK) + ASO2(:,K3)
        end do

        open(vac_output_fileunit, file='VAC') ! TODO: set new filename
        do j=1, JAM
            write(vac_output_fileunit, '(F10.3, F11.5, 32E12.4)') ZA(J), WA(J), VAC(J,:)
        end do
        close(vac_output_fileunit)
    end subroutine kd_calc

    subroutine planck_calc()
        implicit none
        real :: A1(NKTR), A2(NKTR)
        real :: TT
        real :: TPP, C
        real, parameter :: TS_TS=100.
        real, parameter :: HH_HH=0.1
        integer :: I_I
        integer :: l ! loop variable
        
        open(planck_fileunit, access='direct', form='unformatted', recl=128, file='data/F_LW/KD-Planck')
        open(planck_output_fileunit, file='PLANCK')
        do l = 1, JAM
            TT = TA(l)
            I_I = (TT-TS_TS) / HH_HH+1 
            TPP = TS_TS + HH_HH*(I_I-1)
            read(planck_fileunit, REC=I_I) A1
            I_I = I_I+1
            read(planck_fileunit, REC=I_I) A2
            C = (TT-TPP) / HH_HH
            PL = (1.-C)*A1 + C*A2
            write(planck_output_fileunit, '(F10.3, F8.2, 32E12.4)') ZA(l), TA(l), PL(:)
        end do
        close(planck_output_fileunit)
        close(planck_fileunit)
    end subroutine planck_calc

end module KD
