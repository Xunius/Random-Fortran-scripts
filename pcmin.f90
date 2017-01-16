module PCMIN_MOD
implicit none

private

integer, parameter :: ikind=selected_real_kind(p=15,r=10)

!------------Thermal dynamic constants------------
real(kind=ikind), parameter :: CPD = 1005.7       ! Specific heat of dry air, J/kg/K
real(kind=ikind), parameter :: CPV = 1870.0       ! Specific heat of vapor, J/kg/K
real(kind=ikind), parameter :: CL = 2500.0        ! Specific heat of liquid water (Reduced
                                                  ! from 4190 to 2500 to partially offset the
                                                  ! lack of latent heat of fusion, J/kg/K
real(kind=ikind), parameter :: CPVMCL = CPV-CL   
real(kind=ikind), parameter :: RV = 461.5         ! Gas constant of vapor, J/kg/K
real(kind=ikind), parameter :: RD = 287.04        ! Gas constant of dry air, J/kg/K
real(kind=ikind), parameter :: EPS = RD/RV 
real(kind=ikind), parameter :: ALV0 = 2.501E6     ! Latent heat at 0 C, J/kg/K

public :: PCMIN, CAPE

contains


SUBROUTINE PCMIN(SST,PSL,P,T,R,NA,N,PMIN,VMAX,IFL)

    !  Calculates the maximum wind speed and mimimum central pressure            
    !  achievable in tropical cyclones, given a sounding and a sea surface temperature.          
    !
    !  <SST>: real, sea surface temperature in K.
    !  <PSL>: real, sea level pressure in hPa.
    !  <P>, <T>, <R>: 1d real array of dimension <NA>, representing
    !                 pressure (hPa), temperature (K) and mixing ratio (kg/kg).
    !                 The arrays MUST be arranged so that the lowest index corresponds
    !                 to the lowest model level, with increasing index
    !                 corresponding to decreasing pressure. The temperature
    !                 sounding should extend to at least the tropopause and 
    !                 preferably to the lower stratosphere. 
    !  <NA>: integer, the dimension of <P>, <T>, and <R>.
    !  <N>: integer, the actual number of points in the sounding.
    !
    !  Return:  <PMIN>: real, the minimum central pressure, in hPa.
    !           <VMAX>: real, the maximum surface wind speed, in m/s.
    !                   (reduced to reflect surface drag)
    !           <IFL>: integer, return code, a value of 0 means OK; a value of 1
    !              indicates no convergence (hypercane); a value of 2
    !              means that the CAPE routine failed.


    integer, parameter :: ikind=selected_real_kind(p=15,r=10)
    real(kind=ikind) :: SST, PSL, PMIN, VMAX
    integer :: NA, N, IFL
    real(kind=ikind), dimension(NA) :: T, P, R

    !---------------Adjustable constants---------------
    
    ! Ratio of C_k to C_D
    real(kind=ikind), parameter :: CKCD = 0.9

    ! Buoyancy of displaced parcels, 0=Reversible ascent;  1=Pseudo-adiabatic ascent
    integer, parameter :: SIG = 0

    ! if IDISS = 0, no dissipative heating is allowed; otherwise, it is
    integer, parameter :: IDISS = 1

    ! Exponent, b, in assumed profile of azimuthal velocity in eye,
    ! V=V_m(r/r_m)^b. Used only in calculation of central pressure
    real(kind=ikind), parameter :: b = 2.0

    !  Factor to reduce gradient wind to 10 m wind
    real(kind=ikind), parameter :: VREDUC = 0.8

    !--------------Intermediate variables--------------
    real(kind=ikind) :: SSTC,ES0,PM,TOM,TOMS,RAT,TP,RP,PP,RS0,TV1,TVAV,TOA,TO
    real(kind=ikind) :: CAT,CAPEMS,CAPEM,CAPEA,PNEW,PMOLD,FAC,T0,R0
    integer :: IFLAG,NP

    !-----------Normalize certain quantities-----------
    SSTC=SST-273.15

    !------------------Default values------------------
    VMAX=0.0
    PMIN=PSL 
    TOMS=230.0
    IFL=0

    if (SSTC <= 5.0) then
        IFL=0
        return
    end if

    T0=T(1)
    R0=R(1)
    ES0=6.112*EXP(17.67*SSTC/(243.5+SSTC))
    TV1=T0*(1.+R0/EPS)/(1.+R0)

    !-------------Find environmental CAPE-------------
    TP=T0
    RP=R0
    PP=P(1)
    CALL CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEA,TOA,IFLAG)

    if (IFLAG /= 0) then
        IFL=2
    end if

    !-----Begin iteration to find minimum pressure-----
    NP=0
    PM=PSL*0.9
    PMOLD=PM
    PNEW=PM+2.

    do while (abs(PNEW-PMOLD) > 0.5)

        !-------Find CAPE at radius of maximum winds-------
        TP=T0
        PP=min(PM,1000.0)
        RP=EPS*R0*PSL/(PP*(EPS+R0)-R0*PSL)

        call CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEM,TOM,IFLAG)

        if (IFLAG /= 0) then
            IFL=2
        end if

        !-Find saturation CAPE at radius of maximum winds-
        TP=SST
        PP=min(PM,1000.0)
        RP=EPS*ES0/(PP-ES0)

        call CAPE(TP,RP,PP,T,R,P,NA,N,SIG,CAPEMS,TOMS,IFLAG)
        TO=TOMS

        if (IFLAG /= 0) then
            IFL=2
        end if

        RAT=SST/TOMS

        if (IDISS == 0) then
            RAT=1.0
        end if

        !-------Initial estimate of minimum pressure-------
        RS0=RP
        TVAV=0.5*(TV1+SST*(1.+RS0/EPS)/(1.+RS0))
        CAT=CAPEM-CAPEA+0.5*CKCD*RAT*(CAPEMS-CAPEM)
        CAT=max(CAT,0.0)
        PNEW=PSL*exp(-CAT/RD/TVAV)

        !---------------Test for convergence---------------
        PMOLD=PM
        PM=PNEW
        NP=NP+1

        if (NP > 200 .OR. PM < 400) then
            PMIN=PSL
            VMAX=0.
            IFL=2
            return
        end if
    end do

    !-----------------Find pmin, vmax-----------------
    CAT=CAPEM-CAPEA+CKCD*RAT*0.5*(1.+1./b)*(CAPEMS-CAPEM)
    CAT=max(CAT,0.0)
    PMIN=PSL*exp(-CAT/RD/TVAV)

    FAC=max(0.0,(CAPEMS-CAPEM))
    VMAX=VREDUC*sqrt(CKCD*RAT*FAC)

END SUBROUTINE PCMIN





SUBROUTINE CAPE(TP,RP,PP,T,R,P,ND,N,SIG,CAPED,TOB,IFLAG)

    ! Calculates the CAPE of a parcel given a sounding.

    ! <TP>: real, parcel temperature in K.
    ! <RP>: real, parcel mixing ratio in kg/kg.
    ! <PP>: real, parcel pressure in hPa.
    ! <T>: 1d real array, temperature profile in K.
    ! <R>: 1d real array, mixing ratio profile in kg/kg.
    ! <P>: 1d real array, pressure profile in hPa.
    ! <ND>: integer, dimension of arrays <T>, <R> and <P>.
    ! <N>: integer, actual number of points in the profile.
    ! <SIG>: integer, 0=Reversible ascent;  1=Pseudo-adiabatic ascent.
     
    ! Return <CAPED>: real, CAPE in m^2/s^2.
    !        <TOB>: real, temperature at the level of neutral Buoyancy, in K.
    !        <IFLAG>: integer, return code, 0 for successful run; 1 for improper
    !                 sounding; 2 for non-convergence.

    integer, parameter :: ikind=selected_real_kind(p=15,r=10)
    real(kind=ikind) :: TP,RP,PP,CAPED,TOB
    integer :: N, ND, SIG, IFLAG
    real(kind=ikind), dimension(ND) :: T, P, R, TVRDIF

    !--------------Intermediate variables--------------
    real(kind=ikind) :: TPC, EVP, RH, ALV, S, CHI, PLCL, ESP
    real(kind=ikind) :: TLVR, TG, RG, TJC, ES, TGNEW, TC, ENEW, AP, SG, SL
    real(kind=ikind) :: RMEAN, NA, PA, PINB, PMA, PFAC, PAT
    integer :: JMIN, J, NC, INB

    !------------------Default values------------------
    CAPED=0.0
    TOB=T(1)
    TVRDIF=0.
    IFLAG=0

    !---------Check that sounding is suitable---------
    if (RP < 1.0E-6 .OR. TP < 200.0) then
        IFLAG=1
        return
    end if

    !---------Define various parcel quantities---------
    TPC=TP-273.15_ikind
    ESP=6.112_ikind*EXP(17.67_ikind*TPC/(243.5_ikind+TPC))
    EVP=RP*PP/(EPS+RP)
    RH=RP*(PP-ESP)/EPS/ESP
    RH=MIN(RH,1.0)
    ALV=ALV0+CPVMCL*TPC
    S=(CPD+RP*CL)*LOG(TP)-RD*LOG(PP-EVP)+ALV*RP/TP-RP*RV*LOG(RH)            

    !------Find lifted condensation pressure PLCL------
    CHI=TP/(1669.0_ikind-122.0_ikind*RH-TP)
    PLCL=PP*(RH**CHI)

    !----------------Begin updraft loop----------------
    JMIN=1E3

    do J = 1,N
        !-----------Don't lift parcel above 60mb-----------
        if (P(J) < 50.0 .OR. P(J) >= PP) then
            cycle
        end if

        JMIN=min(JMIN,J)

        !-----------Parcel quantities below LCL-----------
        if (P(J) >= PLCL) then
            TG=TP*(P(J)/PP)**(RD/CPD)
            RG=RP

            !----------------Calculate buoyancy----------------
            TLVR=TG*(1.+RG/EPS)/(1.+RG)
            TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J))
        else
            !-----------Parcel quantities above LCL-----------
            TGNEW=T(J)
            TJC=T(J)-273.15_ikind
            ES=6.112_ikind*exp(17.67_ikind*TJC/(243.5_ikind+TJC))
            RG=EPS*ES/(P(J)-ES)

            ! Iteratively calculate lifted parcel temp and mixing ratio
            ! for reversible ascent
            NC=0
            TG=T(J)-2.

            do while (abs(TGNEW-TG) > 0.001)
                TG=TGNEW
                TC=TG-273.15_ikind
                ENEW=6.112_ikind*exp(17.67_ikind*TC/(243.5_ikind+TC))
                RG=EPS*ENEW/(P(J)-ENEW)

                !----Calculate estimates of the rates of change----
                ! of the entropy with T at constant pressure
                ALV=ALV0+CPVMCL*TC
                SL=(CPD+RP*CL+ALV*ALV*RG/(RV*TG*TG))/TG
                SG=(CPD+RP*CL)*log(TG)-RD*log(P(J)-ENEW)+ALV*RG/TG

                NC=NC+1
                AP=1.0

                TGNEW=TG+AP*(S-SG)/SL

                !--------Bail out if things get out of hand--------
                if (NC > 500 .OR. ENEW > (P(J) -1)) then
                    IFLAG=2
                    return
                end if
            end do

            !----------------Calculate buoyancy----------------
            RMEAN=SIG*RG+(1-SIG)*RP
            TLVR=TG*(1+RG/EPS)/(1+RMEAN)
            TVRDIF(J)=TLVR-T(J)*(1.+R(J)/EPS)/(1.+R(J))

        end if
    end do

    !-------Begin loop to find NA, PA, and CAPE-------
    NA=0.0
    PA=0.0

    !-----Find max level of positive buoyancy, INB-----
    INB=1

    do J = N,JMIN,-1
        if (TVRDIF(J) > 0) then
            INB=max(INB,J)
        end if
    end do
    if (INB == JMIN) then
        return
    end if

    !----Find positive and negative areas and CAPE----
    do J = (JMIN+1),INB
        PFAC=RD*(TVRDIF(J)+TVRDIF(J-1))*(P(J-1)-P(J))/(P(J)+P(J-1))
        PA=PA+max(PFAC,0.0)
        NA=NA-min(PFAC,0.0)
    end do

    !-Find area between parcel pressure and 1st level above it--
    PMA=(PP+P(JMIN))
    PFAC=RD*(PP-P(JMIN))/PMA
    PA=PA+PFAC*max(TVRDIF(JMIN),0.0)
    NA=NA-PFAC*min(TVRDIF(JMIN),0.0)

    !---Find residual positive area above INB and TO---
    if (INB < N) then
        PINB=(P(INB+1)*TVRDIF(INB)-P(INB)*TVRDIF(INB+1))/(TVRDIF(INB)-TVRDIF(INB+1))
        PAT=RD*TVRDIF(INB)*(P(INB)-PINB)/(P(INB)+PINB)
        TOB=(T(INB)*(PINB-P(INB+1))+T(INB+1)*(P(INB)-PINB))/(P(INB)-P(INB+1));
    else
        PAT=0.0
        TOB=T(INB)
    end if

    !--------------------Find CAPE--------------------
    CAPED=PA+PAT-NA
    CAPED=max(CAPED,0.0)


end SUBROUTINE CAPE


end module PCMIN_MOD
