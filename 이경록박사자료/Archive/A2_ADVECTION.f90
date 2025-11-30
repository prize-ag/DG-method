! @@@@@

MODULE MODAL_DG_1D_GLOBALS

  IMPLICIT NONE

  ! COUNTERS
  INTEGER                                                            :: ORD ! INDEX OF THE DEGREE FOR LEGENDRE POLYNOMIALS
  INTEGER                                                            :: I   ! SPATIAL (X-AXIS) INDEX
  INTEGER                                                            :: N   ! TIME INDEX
  INTEGER                                                            :: I0  ! GRID REFINEMENT INDEX

  ! REFINEMENT CONTROL
  INTEGER, PARAMETER                                                 :: NLEV     = 4 ! NUMBER OF GRID LEVELS

  ! ORDER OF THE DISCONTINUOUS GALERKIN METHOD AND POLYNOMIALS
  INTEGER, PARAMETER                                                 :: DG_ORD   = 3          ! ORDER OF DG SOLVER
  INTEGER, PARAMETER                                                 :: POLY_ORD = DG_ORD - 1 ! ORDER OF POLYNOMIALS

  ! USEFUL CONSTANTS AND FUNCTIONS
  REAL( 8 ), PARAMETER                                               :: PI  = 4.0D0 * DATAN( 1.0D0 )
  REAL( 8 ), PARAMETER                                               :: CFL = 0.4D0 / DBLE( 2 * POLY_ORD + 1 ) ! CFL NUMBER

  ! SPATIAL VARIABLES
  INTEGER                                                            :: IMAX                                ! ACTIVE NUMBER OF SPATIAL CELLS
  INTEGER,   PARAMETER                                               :: I_BC     = 1                        ! GHOST BOUNDARY CELLS FOR WENO SCHEMES
  REAL( 8 ), PARAMETER                                               :: X_LEFT   = - 1.0D0                  ! LEFT BOUNDARY ON X-AXIS
  REAL( 8 ), PARAMETER                                               :: X_RIGHT  = + 1.0D0                  ! RIGHT BOUNDARY ON X-AXIS
  REAL( 8 ), PARAMETER                                               :: X_LENGTH = DABS( X_RIGHT - X_LEFT ) ! LENGTH OF DOMAIN

  INTEGER,   PARAMETER                                               :: IMAX_MAX = 40 * ( 2 ** NLEV )

  REAL( 8 ), DIMENSION( - I_BC + 1 : IMAX_MAX + I_BC )               :: X     ! CELL CENTERS
  REAL( 8 ), DIMENSION( - I_BC + 1 : IMAX_MAX + I_BC )               :: DX    ! CELL VOLUMES
  REAL( 8 )                                                          :: MAXDX ! MAXIMAL CELL VOLUME
  
  ! TIME VARIABLES
  INTEGER,   PARAMETER                                               :: RK_ORD    = 3     ! ORDER OF RUNGE-KUTTA SOLVER
  INTEGER,   PARAMETER                                               :: NMAX      = 99999 ! NUMBER OF THE MARCHING STEP
  REAL( 8 ), PARAMETER                                               :: TIME_OUT  = 2.0D0 ! FINAL TIME OF THE SOLUTION
  REAL( 8 )                                                          :: DT                ! TIME STEP SIZE
  REAL( 8 )                                                          :: FLAG_0, FLAG_1    ! TIME FLAGS

  ! CPU TIME CHECKERS
  INTEGER                                                            :: RATE, TIME_0, TIME_1

  ! VARIABLES FOR QUADRATURES
  REAL( 8 ), DIMENSION( 1 : 4, 1 : 2 )                               :: Q 

  ! MASS COEFFICIENT MATRIX
  REAL( 8 ), DIMENSION( 0 : POLY_ORD, 0 : POLY_ORD )                 :: M

  ! CELL INTERFACE CONTRIBUTIONS
  REAL( 8 ), DIMENSION( 0 : POLY_ORD )                               :: L_B, R_B

  ! DEGREES OF FREEDOM IN THE MOMENTS BY THE L2 PROJECTION
  REAL( 8 ), DIMENSION( 0 : POLY_ORD, - I_BC + 1 : IMAX_MAX + I_BC ) :: OLD_DEG 
  REAL( 8 ), DIMENSION( 0 : POLY_ORD, 1 : IMAX_MAX )                 :: NEW_DEG 

  ! CONVEX SUMMATION OF DEGREES OF FREEDOM AND LEGENDRE POLYNOMIALS 
  REAL( 8 )                                                          :: TEMP
  REAL( 8 ), DIMENSION( 1 : IMAX_MAX )                               :: U, NEW_U

  ! ERROR STORAGE FOR CONVERGENCE STUDY
  REAL( 8 ), DIMENSION( 0 : NLEV )                                   :: L2_ERR, LINF_ERR
  
  ! PLOTTING VARIABLES
  INTEGER, PARAMETER                                                 :: WIDTH = 36
  INTEGER                                                            :: PAD
  INTEGER                                                            :: LEN_ROW
  CHARACTER( LEN = 80 )                                              :: ROW

END MODULE MODAL_DG_1D_GLOBALS

! @@@@@

PROGRAM MODAL_DG_1D

  USE MODAL_DG_1D_GLOBALS
  USE OMP_LIB

  IMPLICIT NONE

  INTEGER :: NTHREADS

  PRINT *, "================================================================================================="
  PRINT *, "                            SCALAR CONSERVATION LAWS SOLVER. (ADVECTION)                         "
  PRINT *, "================================================================================================="
  PRINT *, "                                    CALCULATIONS HAVE STARTED.                                   "
  PRINT *, "================================================================================================="
  PRINT *, " "

  ! FIX NUMBER OF OPENMP THREADS
  CALL OMP_SET_NUM_THREADS( 1 )
  NTHREADS = OMP_GET_MAX_THREADS()

  PRINT *, "------------------------------------"
  PRINT *, "OMP NUM THREADS: ", NTHREADS
  PRINT *, "------------------------------------"
  PRINT *, " "

  CALL GET_QUADRATURE( Q )

  PRINT *, "------------------------------------"
  PRINT *, "CFL NUMBER:", REAL( CFL )
  PRINT *, "------------------------------------"
  PRINT *, " "

  CALL SYSTEM_CLOCK( COUNT_RATE = RATE )
  CALL SYSTEM_CLOCK( TIME_0 )

  DO I0 = 0, NLEV

      PRINT *, "------------------------------------"
      PRINT *, "GRID LEVEL: ", I0
      IMAX = 40 * ( 2 ** I0 )
      PRINT *, "NUMBER OF CELLS: ", IMAX

      FLAG_0 = 0.0D0
      FLAG_1 = FLAG_0

      CALL GET_SPATIAL( I_BC, IMAX, X_LEFT, X_LENGTH, X, DX, MAXDX )
      CALL GET_MASS   ( POLY_ORD, MAXDX, Q, M )

      ! @@@@@

      ! PRINT *, "------------------------------------"
      ! PRINT *, "------------ MASS MATRIX -----------"
      ! PRINT *, "------------------------------------"

      ! MASS MATRIX PRINT (OPTIONAL)
      ! DO ORD = 0, SIZE( M, 1 ) - 1
      ! 
      !     WRITE( ROW, '( 3I10 )' ) NINT( M( ORD, : ) )
      ! 
      !     LEN_ROW = LEN_TRIM( ROW )
      !     PAD     = ( WIDTH - LEN_ROW ) / 2
      !     PAD     = PAD - 3
      ! 
      !     IF ( PAD < 0 ) PAD = 0
      ! 
      !     WRITE( *, '( A, A )' ) REPEAT( ' ', PAD ), ROW( 1 : LEN_ROW )
      ! 
      ! END DO

      ! PRINT *, "------------------------------------"
      ! PRINT *, " "

      ! @@@@@

      CALL GET_INITIAL( PI, POLY_ORD, Q, M, I_BC, IMAX, X, DX, MAXDX, OLD_DEG, L_B, R_B )

      N = 1

    DO WHILE ( ( FLAG_0 < TIME_OUT ) .AND. ( N <= NMAX ) )

          FLAG_1 = FLAG_0

          ! COMPUTE THE TIME STEP SIZE BY USING THE FIRST-DERIVATIVE OF FLUX
          CALL GET_TIMESTEP( I_BC, IMAX, DX, FLAG_1, TIME_OUT, CFL, DT )
          
          ! UPDATE THE VARIABLE BY USING TVD RK SOLVER
          CALL GET_RK( RK_ORD, POLY_ORD, Q, M, L_B, R_B, I_BC, IMAX, X, DX, MAXDX, DT, OLD_DEG, NEW_DEG )

          ! UPDATE THE REAL-TIME
          FLAG_1 = FLAG_1 + DT
          FLAG_0 = FLAG_1

          ! RECEIVE THE UPDATED DATA
          !$OMP PARALLEL DO PRIVATE( ORD )
          DO I = 1, IMAX

              DO ORD = 0, POLY_ORD

                  OLD_DEG( ORD, I ) = NEW_DEG( ORD, I )

              END DO

          END DO
          !$OMP END PARALLEL DO

          N = N + 1

      END DO
          

      !$OMP PARALLEL DO PRIVATE( ORD, TEMP )
      DO I = 1, IMAX

          TEMP = 0.0D0

          DO ORD = 0, POLY_ORD

              TEMP = TEMP + OLD_DEG( ORD, I ) * POLY( ORD, 0.0D0, DX( I ) )

          END DO

          NEW_U( I ) = TEMP

      END DO
      !$OMP END PARALLEL DO


      ! SET UP THE EXACT SOLUTION AT TIME_OUT FOR THIS GRID
      !$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I )
      DO I = 1, IMAX

          U( I ) = EXT_U( X( I ), PI, TIME_OUT )

      END DO
      !$OMP END PARALLEL DO

      ! COMPUTE VARIOUS ERRORS AND STORE THEM FOR CONVERGENCE STUDY
      PRINT *, "--------- ERRORS (FINAL) -----------"

      CALL GET_L2    ( I_BC, IMAX, DX, U, NEW_U, L2_ERR  ( I0 ) )
      CALL GET_LINFTY(       IMAX,     U, NEW_U, LINF_ERR( I0 ) )

      PRINT *, "L_2 ERROR:     ", REAL( L2_ERR( I0 ) )
      PRINT *, "L_INFTY ERROR: ", REAL( LINF_ERR( I0 ) )
      PRINT *, "------------------------------------"
      PRINT *, " "

  END DO 

  CALL SYSTEM_CLOCK( TIME_1 )

  PRINT *, "---------- L^INFTY ORDERS ----------"
  PRINT '(A, 1X, F12.6)', "1ST STAGE ORDER: ", REAL( LOG2( LINF_ERR( 0 ) / LINF_ERR( 1 ) ) )
  PRINT '(A, 1X, F12.6)', "2ND STAGE ORDER: ", REAL( LOG2( LINF_ERR( 1 ) / LINF_ERR( 2 ) ) )
  PRINT '(A, 1X, F12.6)', "3RD STAGE ORDER: ", REAL( LOG2( LINF_ERR( 2 ) / LINF_ERR( 3 ) ) )
  PRINT '(A, 1X, F12.6)', "4TH STAGE ORDER: ", REAL( LOG2( LINF_ERR( 3 ) / LINF_ERR( 4 ) ) )
  PRINT *, "------------------------------------"

  PRINT *, "-------------- L^2 ORDERS ----------"
  PRINT '(A, 1X, F12.6)', "1ST STAGE ORDER: ", REAL( LOG2( L2_ERR( 0 ) / L2_ERR( 1 ) ) )
  PRINT '(A, 1X, F12.6)', "2ND STAGE ORDER: ", REAL( LOG2( L2_ERR( 1 ) / L2_ERR( 2 ) ) )
  PRINT '(A, 1X, F12.6)', "3RD STAGE ORDER: ", REAL( LOG2( L2_ERR( 2 ) / L2_ERR( 3 ) ) )
  PRINT '(A, 1X, F12.6)', "4TH STAGE ORDER: ", REAL( LOG2( L2_ERR( 3 ) / L2_ERR( 4 ) ) )
  PRINT *, "------------------------------------"

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

  REAL( 8 ) FUNCTION INI_U( X, PI )

    IMPLICIT NONE

    REAL( 8 ), INTENT( IN ) :: X, PI
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!
  
    ! SET UP THE INITIAL FUNCTION
    ! SMOOTH FUNCTIONS
    INI_U = DSIN( PI * X ) ** 1
    ! INI_U = DCOS( PI * X ) ** 1
    ! INI_U = DSIN( PI * X - DSIN( PI * X ) / PI )
  
    ! NONSMOOTH FUNCTION
    ! SQUARE WAVE
    ! IF ( ( X .GT. - 0.5D0 ) .AND. ( X .LT. 0.5D0 ) ) THEN

    !     INI_U = 1.0D0

    ! ELSE

    !     INI_U = 0.0D0

    ! END IF

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION INI_U

  REAL( 8 ) FUNCTION EXT_U( X, PI, TIME_OUT )
  
    IMPLICIT NONE
  
    REAL( 8 ), INTENT( IN ) :: X, PI
    REAL( 8 ), INTENT( IN ) :: TIME_OUT
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!
  
    ! SET UP THE EXACT SOLUTION
    ! SMOOTH FUNCTIONS
    EXT_U = DSIN( PI * ( X - TIME_OUT ) ) ** 1
    ! EXT_U = DCOS( PI * ( X - TIME_OUT ) ) ** 1
    ! EXT_U = DSIN( PI * ( X - TIME_OUT ) - DSIN( PI * ( X - TIME_OUT ) ) / PI )

    ! NONSMOOTH FUNCTION
    ! SQUARE WAVE
    ! IF ( ( X .GT. - 0.5D0 ) .AND. ( X .LT. 0.5D0 ) ) THEN

    !     EXT_U = 1.0D0

    ! ELSE

    !     EXT_U = 0.0D0

    ! END IF

    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END FUNCTION EXT_U

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
    FLUX = U

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
    !$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( ORD_1 )
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
    !$OMP END PARALLEL DO

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
    !$OMP PARALLEL DO DEFAULT( SHARED )
    DO I = 0, IMAX

        TEMP_X( I ) = X_LEFT + DBLE( I ) * DBLE( X_LENGTH / IMAX )

    END DO
    !$OMP END PARALLEL DO
  
    !$OMP PARALLEL DO DEFAULT( SHARED )
    DO I = 1, IMAX

        X ( I ) = ( TEMP_X( I ) + TEMP_X( I - 1 ) ) / 2.0D0
        DX( I ) = DABS( TEMP_X( I ) - TEMP_X( I - 1 ) )

    END DO
    !$OMP END PARALLEL DO

    MAXDX = MAXVAL( DX( 1 : IMAX ) )

    ! SET UP GHOST LOCATIONS AND VOLUMES OF CELLS (VERY SMALL, KEEP SERIAL)
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

    ! APPLY PERIODIC CONDITIONS TO GET GHOST CELLS (VERY SMALL, KEEP SERIAL)
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
    !$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( ORD )
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
    !$OMP END PARALLEL DO
  
    ! COMPUTE THE CONVEX SUMMATION BY USING CALCULATED DEGREES OF FREEDOM
    ! (( X( I ) - X( I ) ) / DX( I ) = 0.0D0 : CELL CENTER )
    !$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( ORD, TEMP )
    DO I = 1, IMAX

      TEMP = 0.0D0

      DO ORD = 0, POLY_ORD

          TEMP         = TEMP + U_0( ORD, I ) * POLY( ORD, 0.0D0, DX( I ) )
          SUM_U_0( I ) = TEMP

      END DO

    END DO
    !$OMP END PARALLEL DO

    ! SET THE EXACT SOLUTION OF THE INITIAL FUNCTION
    !$OMP PARALLEL DO DEFAULT( SHARED )
    DO I = 1, IMAX

        EXACT_U_0( I ) = INI_U( X( I ), PI )

    END DO
    !$OMP END PARALLEL DO
  
    ! COMPUTE DIFFERENCES 
    !$OMP PARALLEL DO DEFAULT( SHARED )
    DO I = 1, IMAX

        E( I, 1 ) = DABS( EXACT_U_0( I ) - SUM_U_0( I ) )
        E( I, 2 ) = ( EXACT_U_0( I ) - SUM_U_0( I ) ) ** 2
        E( I, 2 ) = E( I, 2 ) * DX( I )

    END DO
    !$OMP END PARALLEL DO

    ! COMPUTE THE SQUARE ROOT AVERAGE SUM (SERIAL FOR REPRODUCIBILITY)
    PRINT *, "-------- ERRORS (INITIAL) ----------"
    PRINT *, "L_2 ERROR:     ", REAL( DSQRT( SUM( E( :, 2 ) ) ) )
    PRINT *, "L_INFTY ERROR: ", REAL( MAXVAL( E( :, 1 ) ) )
    PRINT *, "------------------------------------"

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

  SUBROUTINE GET_TIMESTEP( I_BC, IMAX, DX, FLAG_1, TIME_OUT, CFL, DT )

    IMPLICIT NONE
  
    INTEGER,   INTENT( IN )  :: I_BC, IMAX
    REAL( 8 ), INTENT( IN )  :: DX( - I_BC + 1 : IMAX + I_BC )
    REAL( 8 ), INTENT( IN )  :: FLAG_1, TIME_OUT, CFL

    REAL( 8 ), INTENT( OUT ) :: DT
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!
  
    ! DETERMINE THE TIME STEP SIZE
    DT = CFL * MAXVAL( DX ) / 1.0D0
  
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
    REAL( 8 ), DIMENSION( 0 : IMAX + 1 )           :: U_M
    REAL( 8 ), DIMENSION( -1 : IMAX )              :: U_P
   
    ! FLUX SPLITTING VARIABLE
    REAL( 8 ), DIMENSION( 0 : IMAX )               :: HAT_F

    ! BOUNDARY INTEGRATION VARIABLE
    REAL( 8 ), DIMENSION( 0 : POLY_ORD, 1 : IMAX ) :: F_FLUX_B
    !----------------------------------------- CALCULATIONS HAVE STARTED ----------------------------------------!

    ! CALCULATE THE LOCAL RESIDUAL

    !============================================================================================================!
    ! VOLUME INTEGRAL PART
    !============================================================================================================!

    ! COMPUTE THE VOLUME INTEGRAL BY USING QUADRATURES
    !$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, ORD, QUAD_SUM )
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
    !$OMP END PARALLEL DO
   
    !$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, ORD )
    DO I = 1, IMAX

        DO ORD = 0, POLY_ORD

            F_FLUX_V( ORD, I ) =   FLUX( QUAD_U( 1, I ) ) * D_POLY( ORD, Q( 1, 1 ), DX( I ) ) * Q( 1, 2 ) &
                                 + FLUX( QUAD_U( 2, I ) ) * D_POLY( ORD, Q( 2, 1 ), DX( I ) ) * Q( 2, 2 ) &
                                 + FLUX( QUAD_U( 3, I ) ) * D_POLY( ORD, Q( 3, 1 ), DX( I ) ) * Q( 3, 2 ) &
                                 + FLUX( QUAD_U( 4, I ) ) * D_POLY( ORD, Q( 4, 1 ), DX( I ) ) * Q( 4, 2 )

        END DO

    END DO
    !$OMP END PARALLEL DO

    !============================================================================================================!
    ! BOUNDARY INTEGRAL PART
    !============================================================================================================!
  
    CALL GET_BOUNDARY( POLY_ORD, I_BC, IMAX, DEG )

    ! APPLY LIMITING PROCESS (IF NEEDED)
  
    ! FOR GIVEN DEGREES OF FREEDOM, COMPUTE CONVEX SUMMATIONS AT CELL INTERFACES
    !$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, ORD, TEMP )
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
    !$OMP END PARALLEL DO

    ! TAKE THE APPROPRIATE NUMERICAL TRACE 
    ! AS THE AVERAGED MONOTONE FLUX, APPLY LAX-FRIEDRICHS (LF) FLUX
    !$OMP PARALLEL DO DEFAULT( SHARED )
    DO I = 0, IMAX

        HAT_F( I ) = FLUX( U_M( I ) ) + FLUX( U_P( I ) ) + 1.0D0 * U_M( I ) - 1.0D0 * U_P( I )
        HAT_F( I ) = HAT_F( I ) / 2.0D0

    END DO
    !$OMP END PARALLEL DO

    ! COMPUTE THE BOUNDARY INTEGRAL BY USING THE FUNDAMENTAL THEOREM OF CALCULUS
    !$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, ORD )
    DO I = 1, IMAX

        DO ORD = 0, POLY_ORD  

            F_FLUX_B( ORD, I ) = HAT_F( I - 1 ) * L_B( ORD ) - HAT_F( I + 0 ) * R_B( ORD )
 
        END DO
 
    END DO
    !$OMP END PARALLEL DO

    !============================================================================================================!
    ! GATHERING ALL INTEGRALS
    !============================================================================================================!

    ! COMPUTE LOCAL RESIDUALS BY USING CALCULATED INTEGRAL VALUES
    !$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, ORD )
    DO I = 1, IMAX

        DO ORD = 0, POLY_ORD

            F_FLUX( ORD, I ) = M( ORD, ORD ) * ( F_FLUX_B( ORD, I ) + F_FLUX_V( ORD, I ) ) / DX( I )

        END DO

    END DO
    !$OMP END PARALLEL DO
  
    RETURN
    !----------------------------------------- CALCULATIONS HAVE FINISHED ---------------------------------------!
  END SUBROUTINE GET_RESIDUAL

  ! @@@@@
  
END PROGRAM MODAL_DG_1D
