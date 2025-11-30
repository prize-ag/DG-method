! @@@@@

MODULE MODAL_DG_1D_GLOBALS

  IMPLICIT NONE

  ! COUNTERS
  INTEGER                                                            :: ORD ! INDEX OF THE DEGREE FOR LEGENDRE POLYNOMIALS
  INTEGER                                                            :: I   ! SPATIAL (X-AXIS) INDEX
  INTEGER                                                            :: N   ! TIME INDEX
  INTEGER                                                            :: I0  ! GRID REFINEMENT INDEX

  ! REFINEMENT CONTROL
  INTEGER, PARAMETER                                                 :: NLEV      = 4 ! NUMBER OF GRID LEVELS

  ! REFERENCE GRID LEVEL (FOR EXACT / REFERENCE SOLUTION)
  ! INTEGER, PARAMETER                                                 :: IMAX_REF  = 12800 * 1

  ! ORDER OF THE DISCONTINUOUS GALERKIN METHOD AND POLYNOMIALS
  INTEGER, PARAMETER                                                 :: DG_ORD    = 3          ! ORDER OF DG SOLVER
  INTEGER, PARAMETER                                                 :: POLY_ORD  = DG_ORD - 1 ! ORDER OF POLYNOMIALS

  ! USEFUL CONSTANTS AND FUNCTIONS
  REAL( 8 ), PARAMETER                                               :: PI        = 4.0D0 * DATAN( 1.0D0 )
  REAL( 8 ), PARAMETER                                               :: CFL       = 0.1D0 / DBLE( 2 * POLY_ORD + 1 ) ! CFL NUMBER

  ! SPATIAL VARIABLES
  INTEGER                                                            :: IMAX                                 ! ACTIVE NUMBER OF SPATIAL CELLS
  INTEGER,   PARAMETER                                               :: I_BC      = 1                        ! GHOST BOUNDARY CELLS FOR WENO SCHEMES
  REAL( 8 ), PARAMETER                                               :: X_LEFT    = - 1.0D0                  ! LEFT BOUNDARY ON X-AXIS
  REAL( 8 ), PARAMETER                                               :: X_RIGHT   = + 1.0D0                  ! RIGHT BOUNDARY ON X-AXIS
  REAL( 8 ), PARAMETER                                               :: X_LENGTH  = DABS( X_RIGHT - X_LEFT ) ! LENGTH OF DOMAIN

  REAL( 8 ), ALLOCATABLE, DIMENSION( : )                             :: X     ! CELL CENTERS
  REAL( 8 ), ALLOCATABLE, DIMENSION( : )                             :: DX    ! CELL VOLUMES
  REAL( 8 )                                                          :: MAXDX ! MAXIMAL CELL VOLUME
  
  ! TIME VARIABLES
  INTEGER,   PARAMETER                                               :: RK_ORD   = 4          ! ORDER OF RUNGE-KUTTA SOLVER
  INTEGER,   PARAMETER                                               :: NMAX     = 999999999  ! NUMBER OF THE MARCHING STEP
  REAL( 8 ), PARAMETER                                               :: TIME_OUT = 0.3D0      ! FINAL TIME OF THE SOLUTION
  REAL( 8 )                                                          :: DT                    ! TIME STEP SIZE
  REAL( 8 )                                                          :: FLAG_0, FLAG_1        ! TIME FLAGS

  ! CPU TIME CHECKERS
  INTEGER                                                            :: RATE, TIME_0, TIME_1

  ! VARIABLES FOR QUADRATURES
  REAL( 8 ), DIMENSION( 1 : 4, 1 : 2 )                               :: Q 

  ! MASS COEFFICIENT MATRIX
  REAL( 8 ), DIMENSION( 0 : POLY_ORD, 0 : POLY_ORD )                 :: M

  ! CELL INTERFACE CONTRIBUTIONS
  REAL( 8 ), DIMENSION( 0 : POLY_ORD )                               :: L_B, R_B

  ! DEGREES OF FREEDOM IN THE MOMENTS BY THE L2 PROJECTION
  REAL( 8 ), ALLOCATABLE, DIMENSION( :, : )                          :: OLD_DEG 
  REAL( 8 ), ALLOCATABLE, DIMENSION( :, : )                          :: NEW_DEG 

  ! CONVEX SUMMATION OF DEGREES OF FREEDOM AND LEGENDRE POLYNOMIALS 
  REAL( 8 )                                                          :: TEMP
  REAL( 8 ), ALLOCATABLE, DIMENSION( : )                             :: U, NEW_U

  ! REFERENCE GRID STORAGE (ON COARSE GRID POINTS)
  !   REAL( 8 ), ALLOCATABLE, DIMENSION( : )                             :: U_REF
  !   REAL( 8 ), ALLOCATABLE, DIMENSION( : )                             :: X_REF

  ! REFERENCE GRID STORAGE (ON FIXED FINE GRID POINTS)
  !   REAL( 8 ), ALLOCATABLE, DIMENSION( : )                             :: U_FINE_REF
  !   REAL( 8 ), ALLOCATABLE, DIMENSION( : )                             :: X_FINE_REF

  ! ERROR STORAGE FOR CONVERGENCE STUDY
  REAL( 8 ), DIMENSION( 1 : NLEV + 1 )                               :: L2_ERR
    !   REAL( 8 ), DIMENSION( 1 : NLEV + 1 )                               :: LINF_ERR

  ! PLOTTING / STRING VARIABLES
  INTEGER, PARAMETER                                                 :: WIDTH = 36
  INTEGER                                                            :: PAD
  INTEGER                                                            :: LEN_ROW
  CHARACTER( LEN = 80 )                                              :: ROW

  INTEGER                                                            :: UNIT_ID
  CHARACTER( LEN = 40 )                                              :: FNAME

  INTEGER                                                            :: J_REF
  REAL( 8 )                                                          :: RATIO
  REAL( 8 )                                                          :: POS, ALPHA, DX_REF_LOC

  ! REFERENCE GRID FILE NAME
  CHARACTER( LEN = 40 ), PARAMETER                                   :: REF_FILE = 'DG1D_REFERENCE_GRID.TXT'

  ! LOCAL REFINEMENT RATIO BETWEEN COARSE GRID AND REFERENCE GRID (MUST BE ODD)
  INTEGER, PARAMETER                                                 :: REF_RATIO = 3

END MODULE MODAL_DG_1D_GLOBALS

! @@@@@

PROGRAM MODAL_DG_1D

  USE MODAL_DG_1D_GLOBALS

  IMPLICIT NONE

  !LOGICAL :: REF_EXISTS

  PRINT *, "================================================================================================="
  PRINT *, "                            SCALAR CONSERVATION LAWS SOLVER. (BUGERS')                           "
  PRINT *, "================================================================================================="
  PRINT *, "                                    CALCULATIONS HAVE STARTED.                                   "
  PRINT *, "================================================================================================="
  PRINT *, " "

  CALL GET_QUADRATURE( Q )

    !   PRINT *, "------------------------------------"
    !   PRINT *, "CFL NUMBER:", REAL( CFL )
    !   PRINT *, "------------------------------------"
    !   PRINT *, " "

  ! CHECK REFERENCE GRID FILE (FIXED FINE GRID)
  !INQUIRE( FILE = REF_FILE, EXIST = REF_EXISTS )
  
    !   PRINT *, "------------------------------------"
    !   IF ( REF_EXISTS ) THEN

    !       PRINT *, "REFERENCE GRID FILE FOUND."
    !       PRINT *, REF_FILE
    !       PRINT *, "------------------------------------"
    !       CALL LOAD_REFERENCE_FROM_FILE

    !   ELSE

    !       PRINT *, "REFERENCE GRID FILE NOT FOUND."
    !       PRINT *, "COMPUTING REFERENCE GRID ON FINE MESH ..."
    !       PRINT *, "------------------------------------"
    !       CALL COMPUTE_REFERENCE_AND_SAVE

!   END IF

  CALL SYSTEM_CLOCK( COUNT_RATE = RATE )
  CALL SYSTEM_CLOCK( TIME_0 )

  DO I0 = 1, NLEV + 1

    !   PRINT *, "------------------------------------"
    !   PRINT *, "GRID LEVEL: ", I0
      IMAX = 10 * ( 2 ** ( I0 + 1 ) )
    !   PRINT *, "NUMBER OF CELLS: ", IMAX

      FLAG_0 = 0.0D0
      FLAG_1 = FLAG_0

      ALLOCATE( X      ( - I_BC + 1 : IMAX + I_BC ) )
      ALLOCATE( DX     ( - I_BC + 1 : IMAX + I_BC ) )
      ALLOCATE( OLD_DEG( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC ) )
      ALLOCATE( NEW_DEG( 0 : POLY_ORD, 1 : IMAX ) )
      ALLOCATE( U      ( 1 : IMAX ) )
      ALLOCATE( NEW_U  ( 1 : IMAX ) )
    !   ALLOCATE( U_REF  ( 1 : IMAX ) )
    !   ALLOCATE( X_REF  ( 1 : IMAX ) )

      CALL GET_SPATIAL( I_BC, IMAX, X_LEFT, X_LENGTH, X, DX, MAXDX )
      CALL GET_MASS   ( POLY_ORD, MAXDX, Q, M )

      CALL GET_INITIAL( PI, POLY_ORD, Q, M, I_BC, IMAX, X, DX, MAXDX, OLD_DEG, L_B, R_B )

      N = 1

      DO WHILE ( ( FLAG_0 < TIME_OUT ) .AND. ( N <= NMAX ) )

          FLAG_1 = FLAG_0

          ! COMPUTE THE TIME STEP SIZE BY USING THE FIRST-DERIVATIVE OF FLUX
          CALL GET_TIMESTEP( I_BC, IMAX, DX, MAXDX, FLAG_1, TIME_OUT, CFL, OLD_DEG, DT )
          
          ! UPDATE THE VARIABLE BY USING TVD RK SOLVER
          CALL GET_RK      ( RK_ORD, POLY_ORD, Q, M, L_B, R_B, I_BC, IMAX, X, DX, MAXDX, DT, OLD_DEG, NEW_DEG )

          ! UPDATE THE REAL-TIME
          FLAG_1 = FLAG_1 + DT
          FLAG_0 = FLAG_1

          ! RECEIVE THE UPDATED DATA
          DO I = 1, IMAX

              DO ORD = 0, POLY_ORD

                  OLD_DEG( ORD, I ) = NEW_DEG( ORD, I )

              END DO

          END DO

          N = N + 1

      END DO

      ! BUILD NUMERICAL SOLUTION ON COARSE GRID
      DO I = 1, IMAX

          TEMP = 0.0D0

          DO ORD = 0, POLY_ORD

              TEMP = TEMP + OLD_DEG( ORD, I ) * POLY( ORD, 0.0D0, DX( I ) )

          END DO

          NEW_U( I ) = TEMP

      END DO    

      ! INTERPOLATE REFERENCE SOLUTION FROM FIXED FINE GRID TO COARSE GRID CENTERS
      !CALL INTERPOLATE_REFERENCE_ON_COARSE( IMAX, X, U_REF, X_REF )

      ! OPTIONAL GRID ALIGNMENT CHECK FOR THE FIRST GRID LEVEL
    !   IF ( I0 .EQ. 1 ) THEN

    !       PRINT *, "   I           X_COARSE          X_REF(I)"

    !       DO I = 1, IMAX

    !           WRITE( *, '(I4,2X,F16.8,2X,F16.8)' ) I, X( I ), X_REF( I )

    !       END DO

    !   END IF

      ! WRITE GRID DATA (X, U_REF, NEW_U) TO FILE FOR EXTERNAL PLOTTING
      WRITE( FNAME, '(A,I1,A)' ) 'DG1D_GRIDLEVEL_', I0, '.TXT'

      OPEN( NEWUNIT = UNIT_ID, FILE = FNAME, STATUS = 'REPLACE', &
            ACTION = 'WRITE', FORM = 'FORMATTED' )

      DO I = 1, IMAX

          WRITE( UNIT_ID, '(3(1X, ES24.16E3))' ) X( I ), NEW_U( I )

      END DO

      CLOSE( UNIT_ID )

    !   ! COMPUTE ERRORS
    !   PRINT *, "--------- ERRORS (FINAL) -----------"

    !   CALL GET_L2    ( I_BC, IMAX, DX, U_REF, NEW_U, L2_ERR  ( I0 ) )
    !   CALL GET_LINFTY(       IMAX,     U_REF, NEW_U, LINF_ERR( I0 ) )

    !   PRINT *, "L_2 ERROR:     ", REAL( L2_ERR( I0 ) )
    !   PRINT *, "L_INFTY ERROR: ", REAL( LINF_ERR( I0 ) )
    !   PRINT *, "------------------------------------"
    !   PRINT *, " "

      DEALLOCATE( X )
      DEALLOCATE( DX )
      DEALLOCATE( OLD_DEG )
      DEALLOCATE( NEW_DEG )
      DEALLOCATE( U )
      DEALLOCATE( NEW_U )
    !   DEALLOCATE( U_REF )
    !   DEALLOCATE( X_REF )

  END DO 

  CALL SYSTEM_CLOCK( TIME_1 )

    !   PRINT *, "---------- L^INFTY ORDERS ----------"
    !   PRINT '(A, 1X, F12.6)', "1ST STAGE ORDER: ", REAL( LOG2( LINF_ERR( 1 ) / LINF_ERR( 2 ) ) )
    !   PRINT '(A, 1X, F12.6)', "2ND STAGE ORDER: ", REAL( LOG2( LINF_ERR( 2 ) / LINF_ERR( 3 ) ) )
    !   PRINT '(A, 1X, F12.6)', "3RD STAGE ORDER: ", REAL( LOG2( LINF_ERR( 3 ) / LINF_ERR( 4 ) ) )
    !   PRINT '(A, 1X, F12.6)', "4TH STAGE ORDER: ", REAL( LOG2( LINF_ERR( 4 ) / LINF_ERR( 5 ) ) )
    ! !   PRINT '(A, 1X, F12.6)', "5TH STAGE ORDER: ", REAL( LOG2( LINF_ERR( 5 ) / LINF_ERR( 6 ) ) )
    ! !   PRINT '(A, 1X, F12.6)', "6TH STAGE ORDER: ", REAL( LOG2( LINF_ERR( 6 ) / LINF_ERR( 7 ) ) )

    !   PRINT *, "------------------------------------"

    !   PRINT *, "-------------- L^2 ORDERS ----------"
    !   PRINT '(A, 1X, F12.6)', "1ST STAGE ORDER: ", REAL( LOG2( L2_ERR( 1 ) / L2_ERR( 2 ) ) )
    !   PRINT '(A, 1X, F12.6)', "2ND STAGE ORDER: ", REAL( LOG2( L2_ERR( 2 ) / L2_ERR( 3 ) ) )
    !   PRINT '(A, 1X, F12.6)', "3RD STAGE ORDER: ", REAL( LOG2( L2_ERR( 3 ) / L2_ERR( 4 ) ) )
    !   PRINT '(A, 1X, F12.6)', "4TH STAGE ORDER: ", REAL( LOG2( L2_ERR( 4 ) / L2_ERR( 5 ) ) )
    ! !   PRINT '(A, 1X, F12.6)', "5TH STAGE ORDER: ", REAL( LOG2( L2_ERR( 5 ) / L2_ERR( 6 ) ) )
    ! !   PRINT '(A, 1X, F12.6)', "6TH STAGE ORDER: ", REAL( LOG2( L2_ERR( 6 ) / L2_ERR( 7 ) ) )

    !   PRINT *, "------------------------------------"

  PRINT *, "------------------------------------"
  PRINT '(A, F12.6, 1X, A)', "COMPUTATIONAL TIME : ", REAL( TIME_1 - TIME_0 ) / REAL( RATE ), "SEC"
  PRINT *, "------------------------------------"
  PRINT *, " "

  PRINT *, "================================================================================================="
  PRINT *, "                                   CALCULATIONS HAVE FINISHED.                                   "
  PRINT *, "================================================================================================="

  STOP

CONTAINS

  ! @@@@@

    !   SUBROUTINE COMPUTE_REFERENCE_AND_SAVE

    !     IMPLICIT NONE

    !     INTEGER                                                         :: IMAX_FINE
    !     REAL( 8 ), ALLOCATABLE, DIMENSION( : )                          :: X_FINE
    !     REAL( 8 ), ALLOCATABLE, DIMENSION( : )                          :: DX_FINE
    !     REAL( 8 )                                                       :: MAXDX_FINE
    !     REAL( 8 ), ALLOCATABLE, DIMENSION( :, : )                       :: OLD_DEG_FINE
    !     REAL( 8 ), ALLOCATABLE, DIMENSION( :, : )                       :: NEW_DEG_FINE
    !     REAL( 8 ), ALLOCATABLE, DIMENSION( : )                          :: U_FINE
    !     REAL( 8 )                                                       :: FLAG_FINE_0, FLAG_FINE_1
    !     REAL( 8 )                                                       :: DT_FINE
    !     INTEGER                                                         :: I_FINE, ORD_FINE
    !     INTEGER                                                         :: N_FINE
    !     !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    !     IMAX_FINE   = IMAX_REF
    !     FLAG_FINE_0 = 0.0D0
    !     FLAG_FINE_1 = FLAG_FINE_0

    !     ALLOCATE( X_FINE      ( - I_BC + 1 : IMAX_FINE + I_BC ) )
    !     ALLOCATE( DX_FINE     ( - I_BC + 1 : IMAX_FINE + I_BC ) )
    !     ALLOCATE( OLD_DEG_FINE( 0 : POLY_ORD, - I_BC + 1 : IMAX_FINE + I_BC ) )
    !     ALLOCATE( NEW_DEG_FINE( 0 : POLY_ORD, 1 : IMAX_FINE ) )
    !     ALLOCATE( U_FINE      ( 1 : IMAX_FINE ) )

    !     CALL GET_SPATIAL( I_BC, IMAX_FINE, X_LEFT, X_LENGTH, X_FINE, DX_FINE, MAXDX_FINE )
    !     CALL GET_MASS   ( POLY_ORD, MAXDX_FINE, Q, M )

    !     CALL GET_INITIAL( PI, POLY_ORD, Q, M, I_BC, IMAX_FINE, X_FINE, DX_FINE, MAXDX_FINE, OLD_DEG_FINE, L_B, R_B )

    !     N_FINE = 1

    !     DO WHILE ( ( FLAG_FINE_0 < TIME_OUT ) .AND. ( N_FINE <= NMAX ) )

    !         FLAG_FINE_1 = FLAG_FINE_0

    !         CALL GET_TIMESTEP( I_BC, IMAX_FINE, DX_FINE, MAXDX_FINE, FLAG_FINE_1, TIME_OUT, CFL, OLD_DEG_FINE, DT_FINE )
    !         CALL GET_RK      ( RK_ORD, POLY_ORD, Q, M, L_B, R_B, I_BC, IMAX_FINE, X_FINE, DX_FINE, &
    !          MAXDX_FINE, DT_FINE, OLD_DEG_FINE, NEW_DEG_FINE )

    !         FLAG_FINE_1 = FLAG_FINE_1 + DT_FINE
    !         FLAG_FINE_0 = FLAG_FINE_1

    !         DO I_FINE = 1, IMAX_FINE

    !             DO ORD_FINE = 0, POLY_ORD

    !                 OLD_DEG_FINE( ORD_FINE, I_FINE ) = NEW_DEG_FINE( ORD_FINE, I_FINE )

    !             END DO

    !         END DO

    !         N_FINE = N_FINE + 1

    !     END DO

    !     ALLOCATE( X_FINE_REF( 1 : IMAX_REF ) )
    !     ALLOCATE( U_FINE_REF( 1 : IMAX_REF ) )

    !     DO I_FINE = 1, IMAX_FINE

    !         TEMP = 0.0D0

    !         DO ORD_FINE = 0, POLY_ORD

    !             TEMP = TEMP + OLD_DEG_FINE( ORD_FINE, I_FINE ) * POLY( ORD_FINE, 0.0D0, DX_FINE( I_FINE ) )

    !         END DO

    !         U_FINE_REF( I_FINE ) = TEMP
    !         X_FINE_REF( I_FINE ) = X_FINE( I_FINE )

    !     END DO

    !     ! SAVE REFERENCE GRID TO FILE
    !     OPEN( NEWUNIT = UNIT_ID, FILE = REF_FILE, STATUS = 'REPLACE', &
    !           ACTION = 'WRITE', FORM = 'FORMATTED' )

    !     DO I_FINE = 1, IMAX_FINE

    !         WRITE( UNIT_ID, '(I8, 1X, ES24.16E3, 1X, ES24.16E3)' ) I_FINE, X_FINE_REF( I_FINE ), U_FINE_REF( I_FINE )

    !     END DO

    !     CLOSE( UNIT_ID )

    !     DEALLOCATE( X_FINE )
    !     DEALLOCATE( DX_FINE )
    !     DEALLOCATE( OLD_DEG_FINE )
    !     DEALLOCATE( NEW_DEG_FINE )
    !     DEALLOCATE( U_FINE )

    !     PRINT *, "REFERENCE GRID COMPUTED ON FINE MESH AND SAVED TO FILE:"
    !     PRINT *, REF_FILE

    !     RETURN
    !     !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
    !   END SUBROUTINE COMPUTE_REFERENCE_AND_SAVE

    !   SUBROUTINE LOAD_REFERENCE_FROM_FILE

    !     IMPLICIT NONE

    !     INTEGER   :: IOS
    !     REAL( 8 ) :: X_DUMMY
    !     REAL( 8 ) :: U_DUMMY
    !     INTEGER   :: I_DUMMY
    !     INTEGER   :: I_FINE
    !     !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    !     ALLOCATE( X_FINE_REF( 1 : IMAX_REF ) )
    !     ALLOCATE( U_FINE_REF( 1 : IMAX_REF ) )

    !     OPEN( NEWUNIT = UNIT_ID, FILE = REF_FILE, STATUS = 'OLD', &
    !           ACTION = 'READ', FORM = 'FORMATTED' )

    !     DO I_FINE = 1, IMAX_REF

    !         READ( UNIT_ID, *, IOSTAT = IOS ) I_DUMMY, X_DUMMY, U_DUMMY

    !         IF ( IOS /= 0 ) THEN

    !             PRINT *, "ERROR READING REFERENCE FILE AT LINE ", I_FINE
    !             STOP

    !         END IF

    !         X_FINE_REF( I_FINE ) = X_DUMMY
    !         U_FINE_REF( I_FINE ) = U_DUMMY

    !     END DO

    !     CLOSE( UNIT_ID )

    !     PRINT *, "REFERENCE GRID LOADED FROM FILE."
    !     PRINT *, REF_FILE
    !     PRINT *, "------------------------------------"

    !     RETURN
    !     !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
    !   END SUBROUTINE LOAD_REFERENCE_FROM_FILE

    !   SUBROUTINE INTERPOLATE_REFERENCE_ON_COARSE( IMAX_COARSE, X_COARSE, U_REF_COARSE, X_REF_COARSE )

    !     IMPLICIT NONE

    !     INTEGER, INTENT( IN )                                           :: IMAX_COARSE
    !     REAL( 8 ), DIMENSION( - I_BC + 1 : IMAX_COARSE + I_BC ), INTENT( IN )  :: X_COARSE
    !     REAL( 8 ), DIMENSION( 1 : IMAX_COARSE ), INTENT( OUT )          :: U_REF_COARSE
    !     REAL( 8 ), DIMENSION( 1 : IMAX_COARSE ), INTENT( OUT )          :: X_REF_COARSE

    !     INTEGER                                                         :: I_COARSE
    !     INTEGER                                                         :: J0, J1
    !     !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    !     DX_REF_LOC = X_FINE_REF( 2 ) - X_FINE_REF( 1 )

    !     DO I_COARSE = 1, IMAX_COARSE

    !         POS = ( X_COARSE( I_COARSE ) - X_LEFT ) / DX_REF_LOC + 0.5D0

    !         J0  = INT( POS )

    !         IF ( J0 < 1 ) J0 = 1
    !         IF ( J0 > IMAX_REF - 1 ) J0 = IMAX_REF - 1

    !         J1    = J0 + 1
    !         ALPHA = POS - DBLE( J0 )

    !         U_REF_COARSE( I_COARSE ) = ( 1.0D0 - ALPHA ) * U_FINE_REF( J0 ) + ALPHA * U_FINE_REF( J1 )
    !         X_REF_COARSE( I_COARSE ) = X_COARSE( I_COARSE )

    !     END DO

    !     RETURN
    !     !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
    !   END SUBROUTINE INTERPOLATE_REFERENCE_ON_COARSE

  ! @@@@@

  ! @@@@@

  REAL( 8 ) FUNCTION INI_U( X, PI )

    IMPLICIT NONE

    REAL( 8 ), INTENT( IN ) :: X, PI
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!
    
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!
    ! SET UP THE INITIAL FUNCTION
    ! INI_U = DSIN( PI * X ) ** 2
    INI_U = 1.0D0 / 4.0D0 + DSIN( PI * X ) / 2.0D0
    

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION INI_U

  REAL( 8 ) FUNCTION POLY( ORD, X, DX )

    IMPLICIT NONE

    INTEGER,   INTENT( IN ) :: ORD
    REAL( 8 ), INTENT( IN ) :: X, DX

    ! REAL( 8 )               :: LAMBDA
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!
   
    ! SET UP SCALED-ORTHONORMAL BASES
    ! LEGENDRE POLYNOMIALS
    IF ( ORD .EQ. 0 ) THEN

        POLY = X ** 0

    ELSEIF ( ORD .EQ. 1 ) THEN

        POLY = X ** 1

    ELSEIF ( ORD .EQ. 2 ) THEN

        POLY = X ** 2 - ( 1.0D0 / 12.0D0 )

    END IF

    ! DETERMINE SHAPE PARAMETER
    ! LAMBDA = 1.0D0 * 1.0D-2
  
    ! EXPONENTIAL POLYNOMIALS
    ! IF ( ORD .EQ. 0 ) THEN

    !     POLY = ( LAMBDA * X ) ** 0

    ! ELSEIF ( ORD .EQ. 1 ) THEN

    !     POLY = ( 1.0D0  * ( LAMBDA * X ) ** 5 ) / 1.0D0 &
    !          + ( 5.0D0  * ( LAMBDA * X ) ** 3 ) / 1.0D0 &
    !          + ( 15.0D0 * ( LAMBDA * X ) ** 1 ) / 2.0D0

    ! ELSEIF ( ORD .EQ. 2 ) THEN

    !     POLY = ( 1.0D0 * ( LAMBDA * X ) ** 4 ) / 1.0D0 &
    !          + ( 3.0D0 * ( LAMBDA * X ) ** 2 ) / 1.0D0 &
    !          - ( 1.0D0 * ( LAMBDA ** 4 ) ) / 80.0D0    &
    !          - ( 1.0D0 * ( LAMBDA ** 2 ) ) / 4.0D0

    ! END IF

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION POLY

  REAL( 8 ) FUNCTION D_POLY( ORD, X, DX )

    IMPLICIT NONE
  
    INTEGER,   INTENT( IN ) :: ORD
    REAL( 8 ), INTENT( IN ) :: X, DX

    ! REAL( 8 )               :: LAMBDA
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!
  
    ! SET UP DERIVATIVES OF SCALED-ORTHONORMAL BASES
    ! LEGENDRE POLYNOMIALS
    IF ( ORD .EQ. 0 ) THEN

        D_POLY = 0.0D0

    ELSEIF ( ORD .EQ. 1 ) THEN

        D_POLY = X ** 0

    ELSEIF ( ORD .EQ. 2 ) THEN

        D_POLY = 2.0D0 * X ** 1

    END IF

    ! DETERMINE SHAPE PARAMETER
    ! LAMBDA = 1.0D0 * 1.0D-2

    ! EXPONENTIAL POLYNOMIALS
    ! IF ( ORD .EQ. 0 ) THEN

    !     D_POLY = 0.0D0

    ! ELSEIF ( ORD .EQ. 1 ) THEN

    !     D_POLY = ( 5.0D0  * ( LAMBDA ** 5 ) * ( X ** 4 ) ) / 1.0D0 &
    !            + ( 15.0D0 * ( LAMBDA ** 3 ) * ( X ** 2 ) ) / 1.0D0 &
    !            + ( 15.0D0 * ( LAMBDA ** 1 ) * ( X ** 0 ) ) / 2.0D0

    ! ELSEIF ( ORD .EQ. 2 ) THEN

    !     D_POLY = ( 4.0D0 * ( LAMBDA ** 4 ) * ( X ** 3 ) ) / 1.0D0 &
    !            + ( 6.0D0 * ( LAMBDA ** 2 ) * ( X ** 1 ) ) / 1.0D0

    ! END IF

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION D_POLY

  REAL( 8 ) FUNCTION FLUX( U )

    IMPLICIT NONE

    REAL( 8 ), INTENT( IN ) :: U
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!
  
    ! SET UP THE FLUX
    FLUX = ( U ** 2 ) / 2.0D0

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION FLUX

  REAL( 8 ) FUNCTION LOG2( X )

    IMPLICIT NONE

    REAL( 8 ), INTENT( IN )                                      :: X
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    LOG2 = LOG( X ) / LOG( 2.0D0 )

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION LOG2

  SUBROUTINE GET_QUADRATURE( Q )

    IMPLICIT NONE

    REAL( 8 ), INTENT( OUT ), DIMENSION( 1 : 4, 1 : 2 ) :: Q
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! QUADRATURE WEIGHTS ARE SCALED FOR INTEGRATIONS ON [ -1.0D0 / 2.0D0, 1.0D0 / 2.0D0 ]
    ! JACOBIAN OF THE AFFINE MAPPING IS APPLIED TO EACH QUADRATURE WEIGHT

    ! NODES OF 6TH-ORDER GAUSS-LOBATTO QUADRATURE
    Q( 1, 1 ) =  1.0D0 / 2.0D0
    Q( 2, 1 ) = -Q( 1, 1 )
    Q( 3, 1 ) =  DSQRT( 5.0D0 ) / 10.0D0
    Q( 4, 1 ) = -Q( 3, 1 )

    ! WEIGHTS OF 6TH-ORDER GAUSS-LOBATTO QUADRATURE
    Q( 1, 2 ) = 1.0D0 / 12.0D0
    Q( 2, 2 ) = Q( 1, 2 )
    Q( 3, 2 ) = 5.0D0 / 12.0D0
    Q( 4, 2 ) = Q( 3, 2 )

    ! (OPTIONAL) 8TH-ORDER GAUSS-LEGENDRE (COMMENTED)
    ! Q( 1, 1 ) = DSQRT( 3.0D0 / 7.0D0 + ( 2.0D0 / 7.0D0 ) * DSQRT( 6.0D0 / 5.0D0 ) ) / 2.0D0
    ! Q( 2, 1 ) = DSQRT( 3.0D0 / 7.0D0 - ( 2.0D0 / 7.0D0 ) * DSQRT( 6.0D0 / 5.0D0 ) ) / 2.0D0
    ! Q( 3, 1 ) = -Q( 2, 1 )
    ! Q( 4, 1 ) = -Q( 1, 1 )

    ! Q( 1, 2 ) = ( 18.0D0 - DSQRT( 30.0D0 ) ) / 72.0D0
    ! Q( 2, 2 ) = ( 18.0D0 + DSQRT( 30.0D0 ) ) / 72.0D0
    ! Q( 3, 2 ) = Q( 2, 2 )
    ! Q( 4, 2 ) = Q( 1, 2 )
  
    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_QUADRATURE

  SUBROUTINE GET_MASS( POLY_ORD, MAXDX, Q, M )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )  :: POLY_ORD
    REAL( 8 ), INTENT( IN )  :: MAXDX
    REAL( 8 ), INTENT( IN )  :: Q( 1 : 4, 1 : 2 )

    REAL( 8 ), INTENT( OUT ) :: M( 0 : POLY_ORD, 0 : POLY_ORD )

    ! COUNTERS
    INTEGER                  :: ORD_0, ORD_1
    !----------------------------------------- CALCULATIONS HAVE STARTED -----------------------------------------!

    ! COMPUTE THE INVERSE OF THE MASS MATRIX
    DO ORD_0 = 0, POLY_ORD

        DO ORD_1 = 0, POLY_ORD

            M( ORD_0, ORD_1 ) = ( POLY( ORD_0, Q( 1, 1 ), MAXDX ) * POLY( ORD_1, Q( 1, 1 ), MAXDX ) ) * Q( 1, 2 ) &
                              + ( POLY( ORD_0, Q( 2, 1 ), MAXDX ) * POLY( ORD_1, Q( 2, 1 ), MAXDX ) ) * Q( 2, 2 ) &
                              + ( POLY( ORD_0, Q( 3, 1 ), MAXDX ) * POLY( ORD_1, Q( 3, 1 ), MAXDX ) ) * Q( 3, 2 ) &
                              + ( POLY( ORD_0, Q( 4, 1 ), MAXDX ) * POLY( ORD_1, Q( 4, 1 ), MAXDX ) ) * Q( 4, 2 )

            ! INVERSION 
            IF ( M( ORD_0, ORD_1 ) .GE. 1.0D-12 ) THEN

                M( ORD_0, ORD_1 ) = 1.0D0 / M( ORD_0, ORD_1 )

            ELSE

                M( ORD_0, ORD_1 ) = 0.0D0

            END IF

        END DO

    END DO

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ----------------------------------------!
  END SUBROUTINE GET_MASS

  SUBROUTINE GET_SPATIAL( I_BC, IMAX, X_LEFT, X_LENGTH, X, DX, MAXDX )

    IMPLICIT NONE 

    INTEGER,   INTENT( IN )          :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )          :: X_LEFT, X_LENGTH
  
    REAL( 8 ), INTENT( OUT )         :: X ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( OUT )         :: DX( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( OUT )         :: MAXDX

    ! COUNTER
    INTEGER                          :: I

    ! CALCULATION VARIABLE
    REAL( 8 ), DIMENSION( 0 : IMAX ) :: TEMP_X 
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!
  
    ! GENERATE PHYSICAL CELL CENTER LOCATIONS AND VOLUMES OF CELLS
    DO I = 0, IMAX

        TEMP_X( I ) = X_LEFT + DBLE( I ) * DBLE( X_LENGTH / IMAX )

    END DO
  
    DO I = 1, IMAX

        X ( I ) = ( TEMP_X( I ) + TEMP_X( I - 1 ) ) / 2.0D0
        DX( I ) = DABS( TEMP_X( I ) - TEMP_X( I - 1 ) )

    END DO
    

    MAXDX = MAXVAL( DX( 1 : IMAX ) )

    ! SET UP GHOST LOCATIONS AND VOLUMES OF CELLS
    DO I = 1, I_BC

        X ( - I + 1 ) = X( IMAX - I + 1 ) - X_LENGTH
        X ( IMAX + I ) = X( I )           + X_LENGTH

        DX( - I + 1 ) = DX( IMAX - I + 1 )
        DX( IMAX + I ) = DX( I )

    END DO
  
    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_SPATIAL

  SUBROUTINE GET_BOUNDARY( POLY_ORD, I_BC, IMAX, DEG )

    IMPLICIT NONE
 
    INTEGER, INTENT( IN )      :: POLY_ORD
    INTEGER, INTENT( IN )      :: I_BC, IMAX

    REAL( 8 ), INTENT( INOUT ) :: DEG( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )

    ! COUNTERS
    INTEGER                    :: ORD
    INTEGER                    :: I
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! APPLY PERIODIC CONDITIONS TO GET GHOST CELLS
    DO I = 1, I_BC

        DO ORD = 0, POLY_ORD 

            DEG( ORD, - I + 1 ) = DEG( ORD, IMAX - I + 1 )
            DEG( ORD, IMAX + I ) = DEG( ORD, I )

        END DO

    END DO

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_BOUNDARY

  SUBROUTINE GET_INITIAL( PI, POLY_ORD, Q, M, I_BC, IMAX, X, DX, MAXDX, U_0, L_B, R_B )

    IMPLICIT NONE

    REAL( 8 ), INTENT( IN )                 :: PI
    INTEGER,   INTENT( IN )                 :: POLY_ORD
    INTEGER,   INTENT( IN )                 :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )                 :: Q ( 1 : 4, 1 : 2 )
    REAL( 8 ), INTENT( IN )                 :: M ( 0 : POLY_ORD, 0 : POLY_ORD )
    REAL( 8 ), INTENT( IN )                 :: X ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )                 :: DX( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )                 :: MAXDX

    REAL( 8 ), INTENT( OUT )                :: U_0( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( OUT )                :: L_B( 0 : POLY_ORD )
    REAL( 8 ), INTENT( OUT )                :: R_B( 0 : POLY_ORD )
 
    ! COUNTERS
    INTEGER                                 :: ORD
    INTEGER                                 :: I

    ! CALCULATION VARIABLES
    REAL( 8 )                               :: TEMP
    REAL( 8 ), DIMENSION( 1 : IMAX )        :: SUM_U_0, EXACT_U_0
    REAL( 8 ), DIMENSION( 1 : IMAX, 1 : 2 ) :: E
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!
  
    ! COMPUTE L2-PROJECTION OF DATA INTO THE FINITE ELEMENT SPACE
    DO I = 1, IMAX

        ! COMPUTE DEGREES OF FREEDOM OF INITIAL DATA BY USING QUADRATURES
        DO ORD = 0, POLY_ORD

            U_0( ORD, I ) =  INI_U( X( I ) + Q( 1, 1 ) * DX( I ), PI ) * POLY( ORD, Q( 1, 1 ), DX( I ) ) * Q( 1, 2 ) &
                          +  INI_U( X( I ) + Q( 2, 1 ) * DX( I ), PI ) * POLY( ORD, Q( 2, 1 ), DX( I ) ) * Q( 2, 2 ) &
                          +  INI_U( X( I ) + Q( 3, 1 ) * DX( I ), PI ) * POLY( ORD, Q( 3, 1 ), DX( I ) ) * Q( 3, 2 ) &
                          +  INI_U( X( I ) + Q( 4, 1 ) * DX( I ), PI ) * POLY( ORD, Q( 4, 1 ), DX( I ) ) * Q( 4, 2 )

            ! PRODUCT WITH THE MASS MATRIX FOR ORTHONORMALIZATION                 
            U_0( ORD, I ) = M( ORD, ORD ) * U_0( ORD, I )

        END DO

    END DO
  
    ! COMPUTE THE CONVEX SUMMATION BY USING CALCULATED DEGREES OF FREEDOM
    ! (( X( I ) - X( I ) ) / DX( I ) = 0.0D0 : CELL CENTER )
    ! DO I = 1, IMAX

    !   TEMP = 0.0D0

    !   DO ORD = 0, POLY_ORD

    !       TEMP         = TEMP + U_0( ORD, I ) * POLY( ORD, 0.0D0, DX( I ) )
    !       SUM_U_0( I ) = TEMP

    !   END DO

    ! END DO

    ! SET THE EXACT SOLUTION OF THE INITIAL FUNCTION
    ! DO I = 1, IMAX

    !     EXACT_U_0( I ) = INI_U( X( I ), PI )

    ! END DO
  
    ! COMPUTE DIFFERENCES 
    ! DO I = 1, IMAX

    !     E( I, 1 ) = DABS( EXACT_U_0( I ) - SUM_U_0( I ) )
    !     E( I, 2 ) = ( EXACT_U_0( I ) - SUM_U_0( I ) ) ** 2
    !     E( I, 2 ) = E( I, 2 ) * DX( I )

    ! END DO

    ! COMPUTE THE SQUARE ROOT AVERAGE SUM
    ! PRINT *, "-------- ERRORS (INITIAL) ----------"
    ! PRINT *, "L_2 ERROR:     ", REAL( DSQRT( SUM( E( :, 2 ) ) ) )
    ! PRINT *, "L_INFTY ERROR: ", REAL( MAXVAL( E( :, 1 ) ) )
    ! PRINT *, "------------------------------------"

    ! COMPUTE CELL INTERFACE CONTRIBUTIONS
    ! (( X( I ) - DX( I ) / 2.0D0 ) - X( I ) ) / DX( I ) = - 1.0D0 / 2.0D0 : LEFT CELL INTERFACE
    ! (( X( I ) + DX( I ) / 2.0D0 ) - X( I ) ) / DX( I ) = + 1.0D0 / 2.0D0 : RIGHT CELL INTERFACE

    DO ORD = 0, POLY_ORD

      L_B( ORD ) = POLY( ORD, - 1.0D0 / 2.0D0, MAXDX )
      R_B( ORD ) = POLY( ORD, + 1.0D0 / 2.0D0, MAXDX )

    END DO
    
    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_INITIAL

  SUBROUTINE GET_RK( RK_ORD, POLY_ORD, Q, M, L_B, R_B, I_BC, IMAX, X, DX, MAXDX, DT, OLD_DEG, NEW_DEG )

    IMPLICIT NONE
  
    INTEGER,   INTENT( IN )                                               :: RK_ORD
    INTEGER,   INTENT( IN )                                               :: POLY_ORD
    INTEGER,   INTENT( IN )                                               :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )                                               :: Q      ( 1 : 4, 1 : 2 ) 
    REAL( 8 ), INTENT( IN )                                               :: M      ( 0 : POLY_ORD, 0 : POLY_ORD )
    REAL( 8 ), INTENT( IN )                                               :: L_B    ( 0 : POLY_ORD )
    REAL( 8 ), INTENT( IN )                                               :: R_B    ( 0 : POLY_ORD )
    REAL( 8 ), INTENT( IN )                                               :: X      ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )                                               :: DX     ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )                                               :: MAXDX
    REAL( 8 ), INTENT( IN )                                               :: DT

    REAL( 8 ), INTENT( INOUT )                                            :: OLD_DEG( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( OUT )                                              :: NEW_DEG( 0 : POLY_ORD, 1 : IMAX )

    ! COUNTERS
    INTEGER                                                               :: ORD
    INTEGER                                                               :: I

    ! CALCULATION VARIABLES
    REAL( 8 ), DIMENSION( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC, 1 : 3 ) :: TEMP_DEG
    REAL( 8 ), DIMENSION( 0 : POLY_ORD, 1 : IMAX )                        :: F_FLUX
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    IF ( RK_ORD .EQ. 3 ) THEN

        ! 3RD-ORDER TVD RK SOLVER
        !==========================================================================================================!
        ! STAGE 1
        !==========================================================================================================!

        ! COMPUTE THE RESIDUAL
        CALL GET_RESIDUAL( POLY_ORD, Q, M, L_B, R_B, I_BC, IMAX, X, DX, MAXDX, OLD_DEG, F_FLUX )

        DO I = 1, IMAX

            DO ORD = 0, POLY_ORD

              TEMP_DEG( ORD, I, 1 ) = + 1.0D0 * OLD_DEG( ORD, I ) &
                                      + 1.0D0 * DT * F_FLUX( ORD, I )

            END DO

        END DO

        !==========================================================================================================!
        ! STAGE 2
        !==========================================================================================================!

        ! COMPUTE THE RESIDUAL
        CALL GET_RESIDUAL( POLY_ORD, Q, M, L_B, R_B, I_BC, IMAX, X, DX, MAXDX, TEMP_DEG( :, :, 1 ), F_FLUX )

        DO I = 1, IMAX

            DO ORD = 0, POLY_ORD

              TEMP_DEG( ORD, I, 2 ) = + 3.0D0 * OLD_DEG( ORD, I ) &
                                      + 1.0D0 * TEMP_DEG( ORD, I, 1 ) &
                                      + 1.0D0 * DT * F_FLUX( ORD, I )

              TEMP_DEG( ORD, I, 2 ) = TEMP_DEG( ORD, I, 2 ) / 4.0D0

            END DO

        END DO

        !==========================================================================================================!
        ! STAGE 3
        !==========================================================================================================!

        ! COMPUTE THE RESIDUAL
        CALL GET_RESIDUAL( POLY_ORD, Q, M, L_B, R_B, I_BC, IMAX, X, DX, MAXDX, TEMP_DEG( :, :, 2 ), F_FLUX )

        DO I = 1, IMAX
        
            DO ORD = 0, POLY_ORD

                NEW_DEG( ORD, I ) = + 1.0D0 * OLD_DEG( ORD, I ) &
                                    + 2.0D0 * TEMP_DEG( ORD, I, 2 ) &
                                    + 2.0D0 * DT * F_FLUX( ORD, I )

                NEW_DEG( ORD, I ) = NEW_DEG( ORD, I ) / 3.0D0

            END DO

        END DO
    
    ELSEIF ( RK_ORD .EQ. 4 ) THEN

        ! 4TH-ORDER NON-TVD RK SOLVER
        !==========================================================================================================!
        ! STAGE 1
        !==========================================================================================================!
  
        ! COMPUTE THE RESIDUAL
        CALL GET_RESIDUAL( POLY_ORD, Q, M, L_B, R_B, I_BC, IMAX, X, DX, MAXDX, OLD_DEG, F_FLUX )
  
        DO I = 1, IMAX

            DO ORD = 0, POLY_ORD
  
                TEMP_DEG( ORD, I, 1 ) = + 1.0D0 * OLD_DEG( ORD, I ) &
                                        + 0.5D0 * DT * F_FLUX( ORD, I )

            END DO

        END DO

        !==========================================================================================================!
        ! STAGE 2
        !==========================================================================================================!
  
        ! COMPUTE THE RESIDUAL
        CALL GET_RESIDUAL( POLY_ORD, Q, M, L_B, R_B, I_BC, IMAX, X, DX, MAXDX, TEMP_DEG( :, :, 1 ), F_FLUX )

        DO I = 1, IMAX

            DO ORD = 0, POLY_ORD

                TEMP_DEG( ORD, I, 2 ) = + 1.0D0 * OLD_DEG( ORD, I ) &
                                        + ( DT / 2.0D0 ) * F_FLUX( ORD, I )

            END DO

        END DO

        !==========================================================================================================!
        ! STAGE 3
        !==========================================================================================================!

        ! COMPUTE THE RESIDUAL
        CALL GET_RESIDUAL( POLY_ORD, Q, M, L_B, R_B, I_BC, IMAX, X, DX, MAXDX, TEMP_DEG( :, :, 2 ), F_FLUX ) 

        DO I = 1, IMAX
      
            DO ORD = 0, POLY_ORD

                TEMP_DEG( ORD, I, 3 ) = + 1.0D0 * OLD_DEG( ORD, I ) &
                                        + 1.0D0 * DT * F_FLUX( ORD, I )

            END DO 

        END DO
  
        !==========================================================================================================!
        ! STAGE 4
        !==========================================================================================================!

        ! COMPUTE THE RESIDUAL
        CALL GET_RESIDUAL( POLY_ORD, Q, M, L_B, R_B, I_BC, IMAX, X, DX, MAXDX, TEMP_DEG( :, :, 3 ), F_FLUX )
 
        DO I = 1, IMAX

            DO ORD = 0, POLY_ORD

                NEW_DEG( ORD, I ) = - 1.0D0 * OLD_DEG( ORD, I ) &
                                    + 1.0D0 * TEMP_DEG( ORD, I, 1 ) &
                                    + 2.0D0 * TEMP_DEG( ORD, I, 2 ) &
                                    + 1.0D0 * TEMP_DEG( ORD, I, 3 ) &
                                    + ( DT / 2.0D0 ) * F_FLUX( ORD, I )
    
                NEW_DEG( ORD, I ) = NEW_DEG( ORD, I ) / 3.0D0

            END DO

        END DO

    END IF
  
    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_RK

  SUBROUTINE GET_L2( I_BC, IMAX, DX, U_EXACT, U_NUM, ERR_L2 )

    IMPLICIT NONE

    INTEGER,                                          INTENT( IN )  :: I_BC, IMAX
    REAL( 8 ), DIMENSION( - I_BC + 1 : IMAX + I_BC ), INTENT( IN )  :: DX
    REAL( 8 ), DIMENSION( 1 : IMAX ),                 INTENT( IN )  :: U_EXACT, U_NUM

    REAL( 8 ),                                        INTENT( OUT ) :: ERR_L2

    INTEGER                                                         :: I
    REAL( 8 )                                                       :: SUM
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    SUM = 0.0D0

    DO I = 1, IMAX
        SUM = SUM + ( ( U_EXACT( I ) - U_NUM( I ) ) ** 2 ) * DX( I )
    END DO

    ERR_L2 = DSQRT( SUM )

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_L2

  SUBROUTINE GET_LINFTY( IMAX, U_EXACT, U_NUM, ERR_LINF )

    IMPLICIT NONE

    INTEGER,  INTENT( IN )                          :: IMAX
    REAL( 8 ), DIMENSION( 1 : IMAX ), INTENT( IN )  :: U_EXACT, U_NUM

    REAL( 8 ),                        INTENT( OUT ) :: ERR_LINF

    INTEGER                                         :: I
    REAL( 8 )                                       :: DIFF
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ERR_LINF = 0.0D0

    DO I = 1, IMAX

        DIFF = DABS( U_EXACT( I ) - U_NUM( I ) )

        IF ( DIFF .GT. ERR_LINF ) THEN
 
            ERR_LINF = DIFF

        END IF

    END DO 
  
    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_LINFTY

  ! @@@@@

  ! @@@@@

  SUBROUTINE GET_TIMESTEP( I_BC, IMAX, DX, MAXDX, FLAG_1, TIME_OUT, CFL, DEG, DT )

    IMPLICIT NONE
  
    INTEGER,   INTENT( IN )          :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )          :: DX( - I_BC + 1 : IMAX + I_BC ), MAXDX
    REAL( 8 ), INTENT( IN )          :: FLAG_1, TIME_OUT, CFL
    REAL( 8 ), INTENT( IN )          :: DEG( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )

    REAL( 8 ), INTENT( OUT )         :: DT

    ! CALCULATION VARIABLES
    REAL( 8 )                        :: TEMP, F_MAX
    REAL( 8 ), DIMENSION( 1 : IMAX ) :: U
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!
    
    F_MAX = 0.0D0

    DO I = 1, IMAX

        TEMP = 0.0D0

        DO ORD = 0, POLY_ORD

            TEMP   = TEMP + DEG( ORD, I ) * POLY( ORD, 0.0D0, DX( I ) )
            U( I ) = TEMP

        END DO

    END DO

    F_MAX = MAXVAL( DABS( U ) )

    ! DETERMINE THE TIME STEP SIZE
    DT = CFL * MAXDX / F_MAX
  
    ! RESTRICTION
    IF ( ( FLAG_1 + DT ) .GE. TIME_OUT ) THEN

        DT = ( TIME_OUT - FLAG_1 )

    END IF

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_TIMESTEP

  SUBROUTINE GET_RESIDUAL( POLY_ORD, Q, M, L_B, R_B, I_BC, IMAX, X, DX, MAXDX, DEG, F_FLUX )

    IMPLICIT NONE

    INTEGER,   INTENT( IN )                        :: POLY_ORD
    INTEGER,   INTENT( IN )                        :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )                        :: Q     ( 1 : 4, 1 : 2 ) 
    REAL( 8 ), INTENT( IN )                        :: M     ( 0 : POLY_ORD, 0 : POLY_ORD )
    REAL( 8 ), INTENT( IN )                        :: L_B   ( 0 : POLY_ORD )
    REAL( 8 ), INTENT( IN )                        :: R_B   ( 0 : POLY_ORD )
    REAL( 8 ), INTENT( IN )                        :: X     ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )                        :: DX    ( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )                        :: MAXDX

    REAL( 8 ), INTENT( INOUT )                     :: DEG   ( 0 : POLY_ORD, - I_BC + 1 : IMAX + I_BC )

    REAL( 8 ), INTENT( OUT )                       :: F_FLUX( 0 : POLY_ORD, 1 : IMAX )

    ! COUNTERS
    INTEGER                                        :: ORD
    INTEGER                                        :: I

    ! VOLUME INTEGRATION VARIABLES
    REAL( 8 ), DIMENSION( 1 : 4 )                  :: QUAD_SUM
    REAL( 8 ), DIMENSION( 1 : 4, 1 : IMAX )        :: QUAD_U

    ! VOLUME INTEGRATION VARIABLE
    REAL( 8 ), DIMENSION( 0 : POLY_ORD, 1 : IMAX ) :: F_FLUX_V

    ! CELL INTERFACE VALUES
    REAL( 8 ), DIMENSION( 0 : 1 )                  :: TEMP
    REAL( 8 ), DIMENSION( 1 : IMAX )               :: U
    REAL( 8 ), DIMENSION( 0 : IMAX + 1 )           :: U_M
    REAL( 8 ), DIMENSION( -1 : IMAX )              :: U_P
   
    ! FLUX SPLITTING VARIABLE
    REAL( 8 )                                      :: F_MAX
    REAL( 8 ), DIMENSION( 0 : IMAX )               :: HAT_F

    ! BOUNDARY INTEGRATION VARIABLE
    REAL( 8 ), DIMENSION( 0 : POLY_ORD, 1 : IMAX ) :: F_FLUX_B
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! CALCULATE THE LOCAL RESIDUAL

    !============================================================================================================!
    ! VOLUME INTEGRAL PART
    !============================================================================================================!

    ! COMPUTE THE VOLUME INTEGRAL BY USING QUADRATURES
    DO I = 1, IMAX

        QUAD_SUM( : ) = 0.0D0

        DO ORD = 0, POLY_ORD

            QUAD_SUM( 1 ) = QUAD_SUM( 1 ) + DEG( ORD, I ) * POLY( ORD, Q( 1, 1 ), DX( I ) )
            QUAD_SUM( 2 ) = QUAD_SUM( 2 ) + DEG( ORD, I ) * POLY( ORD, Q( 2, 1 ), DX( I ) )
            QUAD_SUM( 3 ) = QUAD_SUM( 3 ) + DEG( ORD, I ) * POLY( ORD, Q( 3, 1 ), DX( I ) )
            QUAD_SUM( 4 ) = QUAD_SUM( 4 ) + DEG( ORD, I ) * POLY( ORD, Q( 4, 1 ), DX( I ) )

        END DO

        QUAD_U( 1, I ) = QUAD_SUM( 1 )
        QUAD_U( 2, I ) = QUAD_SUM( 2 )
        QUAD_U( 3, I ) = QUAD_SUM( 3 )
        QUAD_U( 4, I ) = QUAD_SUM( 4 )

    END DO
   
    DO I = 1, IMAX

        DO ORD = 0, POLY_ORD

            F_FLUX_V( ORD, I ) =   FLUX( QUAD_U( 1, I ) ) * D_POLY( ORD, Q( 1, 1 ), DX( I ) ) * Q( 1, 2 ) &
                                 + FLUX( QUAD_U( 2, I ) ) * D_POLY( ORD, Q( 2, 1 ), DX( I ) ) * Q( 2, 2 ) &
                                 + FLUX( QUAD_U( 3, I ) ) * D_POLY( ORD, Q( 3, 1 ), DX( I ) ) * Q( 3, 2 ) &
                                 + FLUX( QUAD_U( 4, I ) ) * D_POLY( ORD, Q( 4, 1 ), DX( I ) ) * Q( 4, 2 )

        END DO

    END DO

    !============================================================================================================!
    ! BOUNDARY INTEGRAL PART
    !============================================================================================================!
  
    CALL GET_BOUNDARY( POLY_ORD, I_BC, IMAX, DEG )

    ! APPLY LIMITING PROCESS (IF NEEDED)
  
    ! FOR GIVEN DEGREES OF FREEDOM, COMPUTE CONVEX SUMMATIONS AT CELL INTERFACES
    DO I = 0, IMAX + 1

        TEMP( : ) = 0.0D0

        DO ORD = 0, POLY_ORD

            ! LEFT BOUNDARY
            TEMP( 0 ) = TEMP( 0 ) + DEG( ORD, I ) * L_B( ORD )

            ! RIGHT BOUNDARY 
            TEMP( 1 ) = TEMP( 1 ) + DEG( ORD, I ) * R_B( ORD )

        END DO

        U_P( I - 1 ) = TEMP( 0 )
        U_M( I + 0 ) = TEMP( 1 )

    END DO

    ! TAKE THE APPROPRIATE NUMERICAL TRACE 
    ! AS THE AVERAGED MONOTONE FLUX, APPLY LAX-FRIEDRICHS (LF) FLUX
    
    F_MAX = 0.0D0

    DO I = 1, IMAX

        TEMP( 0 ) = 0.0D0

        DO ORD = 0, POLY_ORD

            TEMP( 0 ) = TEMP( 0 ) + DEG( ORD, I ) * POLY( ORD, 0.0D0, DX( I ) )
            U( I )    = TEMP( 0 )

        END DO

    END DO

    F_MAX = MAXVAL( DABS( U ) )

    DO I = 0, IMAX

        HAT_F( I ) = FLUX( U_M( I ) ) + FLUX( U_P( I ) ) + F_MAX * ( U_M( I ) -  U_P( I ) )
        HAT_F( I ) = HAT_F( I ) / 2.0D0

    END DO

    ! COMPUTE THE BOUNDARY INTEGRAL BY USING THE FUNDAMENTAL THEOREM OF CALCULUS
    DO I = 1, IMAX

        DO ORD = 0, POLY_ORD  

            F_FLUX_B( ORD, I ) = HAT_F( I - 1 ) * L_B( ORD ) - HAT_F( I + 0 ) * R_B( ORD )
 
        END DO
 
    END DO

    !============================================================================================================!
    ! GATHERING ALL INTEGRALS
    !============================================================================================================!

    ! COMPUTE LOCAL RESIDUALS BY USING CALCULATED INTEGRAL VALUES
    DO I = 1, IMAX

        DO ORD = 0, POLY_ORD

            F_FLUX( ORD, I ) = M( ORD, ORD ) * ( F_FLUX_B( ORD, I ) + F_FLUX_V( ORD, I ) ) / DX( I )

        END DO

    END DO
  
    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_RESIDUAL

  ! @@@@@

END PROGRAM MODAL_DG_1D