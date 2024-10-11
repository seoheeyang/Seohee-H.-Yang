      MODULE SIMPLE_SOA_MOD
     ! this module is developed by seohee yang (shy, 2023) 
      USE CMN_SIZE_MOD                 ! Size parameters

      IMPLICIT NONE

      ! Make everything PRIVATE ...
      PRIVATE
      ! ... except these routines
      PUBLIC :: DO_SIMPLE_SOA
      PUBLIC :: INIT_SOA
      PUBLIC :: CLEANUP_SOA
      PUBLIC :: SOA_CALCULATE

      ! Arrays
      REAL*8,  ALLOCATABLE :: SOAS(:,:)
      REAL*8,  ALLOCATABLE :: SOAP(:,:)
      REAL*8,  ALLOCATABLE :: BROWN_SOP(:,:)
      REAL*8,  ALLOCATABLE :: A_CM2(:,:)
      REAL*8,  PARAMETER   :: SEC_IN_DAY  = 86400d0
      REAL*8,  PARAMETER   :: SEC_IN_2000 = SEC_IN_DAY * 366d0

      INTEGER             :: DRYSOAS, DRYSOAP, DRYBRSOA, DRYBRSOP

      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_SIMPLE_SOA (I,J,FAC1)

      USE LOGICAL_MOD
      USE TRACER_MOD,         ONLY : STT
      USE TRACERID_MOD,       ONLY : IDTSOAS,   IDTSOAP
      USE TRACERID_MOD,       ONLY : IDTBRSOP,  IDTBRSOA
      USE TRACERID_MOD,       ONLY : IDBCO,     IDBFCO,   IDTCO
      USE BPCH2_MOD,          ONLY : GET_TAU0
      USE TIME_MOD,           ONLY : GET_TS_EMIS,  GET_LOCALTIME
      USE MEGANUT_MOD,        ONLY : XLTMMP
      USE DAO_MOD,            ONLY : SUNCOS_MID, SUNCOS
      USE DAO_MOD,            ONLY : PARDF, PARDR
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE EDGAR_MOD_2,        ONLY : GET_EDGAR_CO, GET_EDGAR_NMVOC
      USE EDGAR_MOD_2,        ONLY : GET_EDGAR_NOx
      USE REAS_MOD,           ONLY : GET_REAS_NMVOC
      USE REAS_MOD,           ONLY : GET_REAS_CO
      USE REAS_MOD,           ONLY : GET_REAS_NOx
      USE TIME_MOD,           ONLY : GET_DAY_OF_WEEK, GET_MONTH
      USE EDGAR_MOD_2,        ONLY : GET_EDGAR_TODN
      USE TRACER_MOD,         ONLY : XNUMOL
      USE BIOFUEL_MOD,        ONLY : BIOFUEL_KG
      USE BIOMASS_MOD,        ONLY : BIOMASS
      USE MEGAN_MOD,          ONLY : GET_EMMONOT_MEGAN
      USE MEGAN_MOD,          ONLY : GET_EMISOP_MEGAN
      USE EPA_NEI_MOD,        ONLY : GET_USA_MASK, GET_EPA_BIOFUEL

      USE CMN_SIZE_MOD               ! Size parameters
      USE CMN_MONOT_MOD              ! Monoterpenes

      ! Local variables
      LOGICAL               :: WEEKDAY
      LOGICAL, SAVE         :: FIRST      = .TRUE.

      INTEGER               :: IJLOOP, MONTH
      INTEGER, INTENT(IN)   :: I       ! GEOS-Chem longitude index
      INTEGER, INTENT(IN)   :: J       ! GEOS-Chem latitude index
      REAL*8 , INTENT(IN)   :: FAC1    ! Factor for inhibition of isop
      REAL*8                :: BRRAT
      REAL*8                :: TMMP,  VALUE
      REAL*8                :: TODX, DTSRCE, DOW
      REAL*8                :: SC, PDF, PDR
      REAL*8                :: EMIS, EMMO
      REAL*8,  EXTERNAL     :: EMMONOT, EMISOP
      REAL*8,  PARAMETER    :: XNUMOL_CO   = 6.0225d23  / 28d-3
      REAL*8                :: EEE, BFF, BMS, ISO, MON, REAS
      REAL*8                :: AVOC, NOXRAT, REAS_NM


      DTSRCE   = GET_TS_EMIS() * 60d0
      DOW      = GET_DAY_OF_WEEK()
      WEEKDAY  = ( DOW > 0 .and. DOW < 6 )
      MONTH    = GET_MONTH()

       IF ( LNEI99 ) THEN
         IF ( GET_USA_MASK (I,J)>0d0 ) THEN

          BFF = ( GET_EPA_BIOFUEL ( I, J, IDTCO, WEEKDAY )*
     &    ( DTSRCE * A_CM2(I,J)) / XNUMOL_CO )
     
         ENDIF
       ELSE

        BFF = BIOFUEL_KG(IDBFCO,I,J)*DTSRCE / SEC_IN_2000
       ENDIF

      IF ( LBIOMASS ) THEN

        IF ( BIOMASS(I,J,IDBCO) > 0d0 ) THEN 
         BMS = BIOMASS(I,J,IDBCO) * (DTSRCE * A_CM2(I,J)) / XNUMOL_CO
        ELSE
         BMS=0d0
        ENDIF

      ENDIF

        ! 1-D loop index
        IJLOOP         = ( (J-1) * IIPAR ) + I
        ! Monoterpene

        ! Surface temperature [K]
        TMMP           = XLTMMP(I,J,IJLOOP)
            ! +++++++++++++++++++++++++++++++++++++++++
            ! ** MEGAN v2.1 **
            ! Cosine of solar zenith angle   (mpb,2009) 
            ! Diffuse and direct PAR         (mpb,2009)
            SC   = SUNCOS_MID(IJLOOP)
            PDR  = PARDR(I,J)
            PDF  = PARDF(I,J)
            ! +++++++++++++++++++++++++++++++++++++++++
            ! Get monoterpenes from MEGAN or GEIA [kg C/box]
            IF ( LMEGANMONO ) THEN
               ! +++++++++++++++++++++++++++++++++++++++++
               ! New formulation for MEGAN (mpb,2009)
               EMMO = GET_EMMONOT_MEGAN( I, J, SC, TMMP,
     &                                   PDR, PDF, 1d0 )
               ! +++++++++++++++++++++++++++++++++++++++++
               ! EMMO = GET_EMMONOT_MEGAN( I, J, TMMP, 1d0 )
            ELSE
               EMMO = EMMONOT( IJLOOP, TMMP, 1d0 )

            ENDIF

               MON= EMMO

        ! Isoprene

             IF ( LMEGAN ) THEN
                !------------------
                ! MEGAN Isoprene
                !------------------

                ! Cosine of solar zenith angle
!------------------------------------------------------------------------------
! Prior to 10/7/11:
! Now call MEGAN emissions with the cos(SZA) computed @ the midpoint of
! the chemistry timestep (bmy, 10/7/11)
!                  SC   = SUNCOS(IJLOOP)
!------------------------------------------------------------------------------
                SC   = SUNCOS_MID(IJLOOP)

                ! Diffuse and direct PAR
                PDR  = PARDR(I,J)
                PDF  = PARDF(I,J)

                 ! Isoprene         
                EMIS = GET_EMISOP_MEGAN(  I, J,     SC, TMMP,
     &                                      PDR, PDF, 1d0 )
             ELSE
                !------------------
                ! GEIA Isoprene 
                !------------------
                 ! Isoprene
                EMIS = EMISOP( I,J,IJLOOP,SUNCOS_MID,TMMP,1d0 )
             ENDIF
              ISO = EMIS * FAC1 ! Isoprene inhibition factor
!-------------------------------------------------------------------------------

        EEE = GET_EDGAR_CO ( I, J,MOLEC_CM2_S=.TRUE.) 
        EEE = EEE * (DTSRCE * A_CM2(I,J)) / XNUMOL_CO !kg/timesep CO

        IF ( LREAS )THEN

        REAS= GET_REAS_CO( I, J, MOLEC_CM2_S=.TRUE.)
        IF ( REAS .GT. 0 ) THEN
        EEE = REAS * (DTSRCE * A_CM2(I,J)) / XNUMOL_CO !kg/timesep CO
        ENDIF
        
        ENDIF

        SOAP(I,J) = (MON*0.05 )+(ISO*0.015)+(EEE*0.069)+
     &              (BFF*0.013)+(BMS*0.013)
        SOAS(I,J) = (MON*0.05 )+(ISO*0.015)


        IF ( LBRC ) THEN 

        IF ( MONTH .EQ. 1  ) BRRAT=37.00/100.
        IF ( MONTH .EQ. 2  ) BRRAT=37.30/100.
        IF ( MONTH .EQ. 3  ) BRRAT=40.62/100.
        IF ( MONTH .EQ. 4  ) BRRAT=44.91/100.
        IF ( MONTH .EQ. 5  ) BRRAT=50.67/100.
        IF ( MONTH .EQ. 6  ) BRRAT=54.84/100.
        IF ( MONTH .EQ. 7  ) BRRAT=54.31/100.
        IF ( MONTH .EQ. 8  ) BRRAT=50.88/100.
        IF ( MONTH .EQ. 9  ) BRRAT=46.41/100.
        IF ( MONTH .EQ. 10 ) BRRAT=43.28/100.
        IF ( MONTH .EQ. 11 ) BRRAT=38.44/100.
        IF ( MONTH .EQ. 12 ) BRRAT=34.92/100.
        
        BROWN_SOP(I,J) =((EEE*0.069))
     &                      * BRRAT !(VBS)
!                            2000-2010 YEAR
!                            ((ASOA1+ASOA2)/ASOA(1+2+3+4+5+6))AVG
!                            NOT ENSEMBLE, CODE 6 LONG_EV FILE
        SOAP(I,J)=SOAP(I,J)-BROWN_SOP(I,J)
        

        ENDIF
          
      END SUBROUTINE DO_SIMPLE_SOA

!------------------------------------------------------------------------------

      SUBROUTINE SOA_CALCULATE

      USE TRACER_MOD,      ONLY : STT
      USE TRACERID_MOD,    ONLY : IDTSOAS, IDTSOAP
      USE TRACERID_MOD,    ONLY : IDTBRSOA, IDTBRSOP
      USE DRYDEP_MOD,      ONLY : DEPNAME, NUMDEP
      USE TIME_MOD
 

      LOGICAL :: WEEKDAY
      LOGICAL, SAVE         :: FIRST      = .TRUE.
      REAL*8                :: NEWSOA, NEWBRC
      REAL*8                :: SOAP_LIFETIME,BRCP_LIFETIME
      REAL*8                :: CONVERT(NVEGTYPE)
      REAL*8                :: GMONOT(NVEGTYPE)
      REAL*8                :: TMMP, EMMO, VALUE
      REAL*8                :: DTSRCE, NOXRAT, AVOC
      INTEGER               :: I, J, L, N

      INTEGER             :: MONTH, YEAR
      REAL*8,  EXTERNAL   :: INHI_ISOP
      REAL*8              :: INHI


      DTSRCE        = GET_TS_EMIS() * 60d0
      SOAP_LIFETIME = 1.0d0 * 24.0d0 * 60.0d0 * 60.0d0
      BRCP_LIFETIME = 1.0d0 * 24.0d0 * 60.0d0 * 60.0d0
      IF ( FIRST ) THEN
       CALL INIT_SOA
       CALL RDLIGHT
       CALL RDISOPT( CONVERT )
       CALL RDMONOT( GMONOT  )
       CALL SETBASE( CONVERT, GMONOT )
 
       DO N = 1, NUMDEP
         SELECT CASE ( TRIM( DEPNAME(N) ) )
            CASE ( 'SOAP' )
               DRYSOAP = N
            CASE ( 'SOAS' )
               DRYSOAS = N
            CASE ( 'BRSOA' ) ! BRC
               DRYBRSOA = N
            CASE ( 'BRSOP' ) ! BRC
               DRYBRSOP = N
            CASE DEFAULT
               ! Nothing
         END SELECT
       ENDDO

       FIRST = .FALSE.

      ENDIF
      YEAR=GET_YEAR()
      MONTH=GET_MONTH()
      INHI = INHI_ISOP (YEAR, MONTH)
      
      IF (INHI .EQ. 0. .OR. INHI .EQ. 1. ) THEN
      PRINT*,'INHI ERROR',INHI
      STOP
      ENDIF
      IF (INHI .LE. 0. .OR. INHI .GE. 2. ) THEN
      PRINT*,'INHI ERROR',INHI
      STOP
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J=1, JJPAR
      DO I=1, IIPAR
        CALL DO_SIMPLE_SOA(I,J, INHI)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL EMITHIGH_SOA( SOAS, SOAP, BROWN_SOP)

      NEWSOA = DEXP( -DTSRCE / SOAP_LIFETIME ) !DTCHEM
      NEWBRC = DEXP( -DTSRCE / BRCP_LIFETIME ) !DTCHEM

      STT(:,:,:,IDTSOAS) = STT(:,:,:,IDTSOAS) +
     &                      (STT(:,:,:,IDTSOAP) * (1.0d0 - NEWSOA))
      STT(:,:,:,IDTSOAP) = STT(:,:,:,IDTSOAP) * NEWSOA

      STT(:,:,:,IDTBRSOA) = STT(:,:,:,IDTBRSOA) +
     &                      (STT(:,:,:,IDTBRSOP) * (1.0d0 - NEWBRC))
      STT(:,:,:,IDTBRSOP) = STT(:,:,:,IDTBRSOP) * NEWBRC

      ! YSH
      CALL SOA_DEPO( STT(:,:,:,IDTSOAP),   DRYSOAP,   IDTSOAP   ) 
      CALL SOA_DEPO( STT(:,:,:,IDTSOAS),   DRYSOAS,   IDTSOAS   )
      CALL SOA_DEPO( STT(:,:,:,IDTBRSOA),  DRYBRSOA,  IDTBRSOA  )
      CALL SOA_DEPO( STT(:,:,:,IDTBRSOP),  DRYBRSOP,  IDTBRSOP  )

      END SUBROUTINE SOA_CALCULATE

!------------------------------------------------------------------------------

      SUBROUTINE EMITHIGH_SOA(SPEC1, SPEC2, SPEC3, SPEC4)

      ! References to F90 modules
      USE PBL_MIX_MOD,  ONLY  : GET_FRAC_OF_PBL,  GET_PBL_MAX_L
      USE TRACER_MOD,   ONLY  : STT
      USE LOGICAL_MOD,  ONLY  : LBRC
      USE CMN_SIZE_MOD  ! Size parameters
      USE TRACERID_MOD, ONLY  : IDTSOAS, IDTSOAP, IDTBRSOP

      ! Arguments
      REAL*8, INTENT(IN) :: SPEC1(IIPAR,JJPAR)
      REAL*8, INTENT(IN) :: SPEC2(IIPAR,JJPAR)
      REAL*8, INTENT(IN), OPTIONAL :: SPEC3(IIPAR,JJPAR)
      REAL*8, INTENT(IN), OPTIONAL :: SPEC4(IIPAR,JJPAR)

      ! Local variables
      INTEGER            :: I,  J,  L,  PBL_MAX
      REAL*8             :: F_OF_PBL


      ! Maximum extent of PBL [model levels]
      PBL_MAX = GET_PBL_MAX_L()

      !=================================================================
      ! Partition emissions throughout the boundary layer
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_OF_PBL )
      DO L = 1, PBL_MAX
      DO J = 1, JJPAR
      DO I = 1, IIPAR

        F_OF_PBL = GET_FRAC_OF_PBL( I, J, L )

        STT(I,J,L,IDTSOAS) = STT(I,J,L,IDTSOAS) +
     &                        ( F_OF_PBL * SPEC1(I,J ) )

        STT(I,J,L,IDTSOAP) = STT(I,J,L,IDTSOAP) +
     &                        ( F_OF_PBL * SPEC2(I,J) )
        IF ( LBRC ) THEN

        STT(I,J,L,IDTBRSOP) = STT(I,J,L,IDTBRSOP) +
     &                        ( F_OF_PBL * SPEC3(I,J) )

        ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE EMITHIGH_SOA

!------------------------------------------------------------------------------

      SUBROUTINE SOA_DEPO( TC, DEPID, TRID )
!
!******************************************************************************
!  Subroutine SOA_DEPO computes dry-deposition of a particular SOA
!  species.
!  (rjp, bmy, 7/8/04, 10/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TC    (REAL*8 ) : Array of SOA tracer 
!  (2 ) DEPID (INTEGER) : Dry deposition ID # (from DEPVEL) 
!  (3 ) TRID  (INTEGER) : GEOS-CHEM tracer number 
! 
!  NOTES:
!  (1 ) Remove reference to CMN, it's obsolete (bmy, 7/20/04)
!  (2 ) Replace PBLFRAC from "drydep_mod.f" with  GET_FRAC_UNDER_PBLTOP
!  from 
!        "pbl_mix_mod.f" (bmy, 2/17/05)
!  (3 ) Bug fix: Add BL_FRAC to the PRIVATE list (mak, bmy, 10/3/05)
!  (4 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (5 ) Add non-local PBL scheme option for dry deposition (lin,
!  06/09/08)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : AD
      USE DIAG_MOD,     ONLY : AD44
      USE DRYDEP_MOD,   ONLY : DEPSAV
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACER_MOD,   ONLY : XNUMOL
      USE LOGICAL_MOD
      USE CMN_SIZE_MOD     ! Size parameters
      USE CMN_DIAG_MOD     ! ND44
      ! Arguments
      REAL*8,  INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR)
      INTEGER, INTENT(IN)    :: DEPID, TRID

      ! Local variable
      INTEGER                :: I, J, L
      REAL*8                 :: DTCHEM, FLUX, BL_FRAC
      REAL*8                 :: TC0, CNEW, FREQ, AREA_CM2
      REAL*8                 :: ND44_TMP(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! SOA_DEPO begins here!
      !=================================================================
      ! Return if tracer ID or tracer ID is undefined
      IF ( TRID == 0 .OR. DEPID == 0 ) RETURN

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

      ! Initialize for drydep diagnostic
      IF ( ND44 > 0 ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            ND44_TMP(I,J,L) = 0d0
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF


!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TC0, BL_FRAC, FREQ, CNEW, AREA_CM2, FLUX )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ! Initial SOA [kg]
         TC0 = TC(I,J,L)
         ! Fraction of box under the PBL top [unitless]
         BL_FRAC = GET_FRAC_UNDER_PBLTOP( I, J, L )
         ! Move drydep to vdiff_mod.f for non-local PBL mixing (Lin,
         ! 06/09/08) 
         IF (LNLPBL) BL_FRAC = 0.D0

         ! Only apply drydep to boxes w/in the PBL
         IF ( BL_FRAC > 0d0 ) THEN

            ! Drydep frequency [1/s]
            FREQ = DEPSAV(I,J,DEPID) * BL_FRAC
            ! Amount of SOA[G] left after drydep [kg]
            CNEW = TC0 * EXP( -FREQ * DTCHEM )
            !===========================================================
            ! ND44 diagnostic: drydep loss [atoms C/cm2/s]
            !===========================================================
            IF ( ND44 > 0 ) THEN

               ! Surface area [cm2]
               AREA_CM2 = GET_AREA_CM2(I, J )

               ! Convert drydep loss from [kg/timestep] to [atoms
               ! C/cm2/s]
               FLUX = ( TC0 - CNEW )
               FLUX = FLUX * XNUMOL(TRID) / ( AREA_CM2 * DTCHEM )

               ! Store in ND44_TMP as a placeholder
               ND44_TMP(I,J,L) = ND44_TMP(I,J,L) + FLUX
            ENDIF

         ELSE
            ! Otherwise, avoid doing the exponential
            ! to preserve precision and clock cycles
            CNEW = TC0

         ENDIF

         ! Prevent underflow condition
         IF ( CNEW < 1d-20 ) CNEW = 0d0
         ! Store modified concentration back in tracer array [kg]
         TC(I,J,L) = CNEW

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO  
      !=================================================================
      ! ND44: Sum drydep fluxes by level into the AD44 array in
      ! order to ensure that  we get the same results w/ sp or mp 
      !=================================================================
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR
            AD44(I,J,DEPID,1) = AD44(I,J,DEPID,1) + ND44_TMP(I,J,L)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE SOA_DEPO


!------------------------------------------------------------------------------

      FUNCTION GET_CO( I, J, L ) RESULT( CO_MOLEC_CM3 )
!
!******************************************************************************
!  Function GET_CO returns OH from SMVGEAR's CSPEC array (for coupled
!  runs)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical
!  level
!******************************************************************************

      ! References to F90 modules
      USE COMODE_MOD,    ONLY : CSPEC, JLOP
      USE DAO_MOD,       ONLY : SUNCOS
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE TIME_MOD,      ONLY : GET_TS_CHEM
      USE TRACER_MOD,    ONLY : ITS_A_FULLCHEM_SIM
      USE TRACERID_MOD,  ONLY : IDCO

      USE CMN_SIZE_MOD  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER               :: JLOOP
      REAL*8                :: CO_MOLEC_CM3
      REAL*8,  PARAMETER    :: XNUMOL_CO   = 6.0225d23  / 28d-3
      ! External functions
      REAL*8,  EXTERNAL     :: BOXVL
      !=================================================================
      ! GET_OH begins here!
      !=================================================================

         ! Take CO from the SMVGEAR array CSPEC
         ! CO is defined only in the troposphere

          CO_MOLEC_CM3 = CSPEC(I,J,L,IDCO)
         ! Convert from molec/cm3 to kg
          CO_MOLEC_CM3 = CO_MOLEC_CM3 * A_CM2(I,J)* BOXVL(I,J,L)
     &     / XNUMOL_CO

     
      ! Return to calling program
      END FUNCTION GET_CO

!------------------------------------------------------------------------------

      SUBROUTINE INIT_SOA
      ! Reference to F90 modules
      USE GRID_MOD,    ONLY : GET_AREA_CM2
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD!, ONLY : LEDGARSHIP

      USE CMN_SIZE_MOD  ! Size parameters

      ! Local Variables
      INTEGER            :: AS, I, J

      
      ALLOCATE( SOAS( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SOAS' )
      SOAS = 0d0
      
      ALLOCATE( SOAP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SOAP' )
      SOAP = 0d0

      ALLOCATE( BROWN_SOP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BROWN_SOP' )
      BROWN_SOP = 0d0

      ALLOCATE( A_CM2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'A_CM2' )

      DO J = 1, JJPAR
      DO I = 1, IIPAR
         A_CM2(I,J) = GET_AREA_CM2( I, J )
      ENDDO

      ENDDO
      ! Return to calling program
      END SUBROUTINE INIT_SOA

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_SOA
!
      IF ( ALLOCATED( A_CM2 ) ) DEALLOCATE( A_CM2 )
      IF ( ALLOCATED( SOAP  ) ) DEALLOCATE( SOAP  )
      IF ( ALLOCATED( SOAS  ) ) DEALLOCATE( SOAS  )
      IF ( ALLOCATED( BROWN_SOP ) ) DEALLOCATE( BROWN_SOP )
      END SUBROUTINE CLEANUP_SOA

!------------------------------------------------------------------------------

      ! End of module
      END MODULE SIMPLE_SOA_MOD
