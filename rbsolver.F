
!------------------------------------------------------------------------!
!  The Community Multiscale Air Quality (CMAQ) system software is in     !
!  continuous development by various groups and is based on information  !
!  from these groups: Federal Government employees, contractors working  !
!  within a United States Government contract, and non-Federal sources   !
!  including research institutions.  These groups give the Government    !
!  permission to use, prepare derivative works of, and distribute copies !
!  of their work in the CMAQ system to the public and to permit others   !
!  to do so.  The United States Environmental Protection Agency          !
!  therefore grants similar permission to use the CMAQ system software,  !
!  but users are requested to provide copies of derivative works or      !
!  products designed to operate in the CMAQ system to the United States  !
!  Government without restrictions as to use by others.  Software        !
!  that is used with the CMAQ system but distributed under the GNU       !
!  General Public License or the GNU Lesser General Public License is    !
!  subject to their copyright restrictions.                              !
!------------------------------------------------------------------------!

C RCS file, release, date & time of last delta, author, state, [and locker]
C $Header: /project/yoj/arc/CCTM/src/gas/ros3/rbsolver.F,v 1.5 2011/10/21 16:11:11 yoj Exp $ 

C what(1) key, module and SID; SCCS file; date and time of last delta:
C %W% %P% %G% %U%

       SUBROUTINE RBSOLVER( JDATE, JTIME, CHEMSTEP, NCSP,
     &                      LIRRFLAG, NIRRCLS, IRRCELL )

C**********************************************************************
C
C  Function: ROS3 solver 
C
C  Preconditions: None
C
C  Key Subroutines/Functions Called: 
C                                    DEGRADE_BLK
C
C  Revision History: Prototype created by Jerry Gipson, August, 2004
C                    Based on the solver described by Sandu et al
C                    ( Atm. Env., Vol. 31, No. 20, 1997 ) and included
C                    in the Kinetic PreProcessor-KPP ( see for example 
C                    Sandu et al., At, Env., Vol. 37, 5097-5114, 
C                    2003). This code also incorporates efficiency
C                    concepts originally developed by M. Jacobson
C                    for SMVGEAR (Atm. Env., Vol 28, No 2, 1994)
C                                          
C                    31 Jan 05 J.Young: get BLKSIZE from dyn alloc horizontal
C                    & vertical domain specifications module (GRID_CONF)
C                    29 Jul 05     WTH: Added IF blocks that call degrade 
C                                       routines if CALL_DEG is true, i.e.,
C                                       if MECHNAME contains 'TX' substring.
C                    28 Jun 10 J.Young: convert for Namelist redesign
C                    29 Mar 11 S.Roselle: Replaced I/O API include files
C                               with UTILIO_DEFN
C                    15 Jul 14 B.Hutzell: replaced mechanism include files with 
C                    RXNS_DATA module, updated explicit interfaces and revised
C                    block for FPP flag redebug is set in compilation
C                    21 Mar 22 G . Sarwar: changed BLKLAND to BLKSEAWATER
C**********************************************************************
!     USE GRID_CONF               ! horizontal & vertical domain specifications
      USE RXNS_DATA
      USE RBDATA                  ! ROS3 solver data
      USE CGRID_SPCS              ! CGRID mechanism species
      USE UTILIO_DEFN
      USE DEGRADE_ROUTINES, ONLY: DEGRADE_BLK
      USE PA_IRR_MODULE
      USE RXNS_FUNCTION
      
      IMPLICIT NONE 

C..Includes:
      INCLUDE '/global/scratch/projects/fc_aqmel/NSF_CDSandE/rocky8/compiled-programs/installations/openmpi/include/mpif.h'

C..Arguments:
      INTEGER,   INTENT( IN ) :: JDATE         ! Current date (YYYYDDD)
      INTEGER,   INTENT( IN ) :: JTIME         ! Current time (HHMMSS)
      REAL( 8 ), INTENT( IN ) :: CHEMSTEP      ! Chem integration interval (min)
      INTEGER,   INTENT( IN ) :: NCSP          ! Index of chem mech to use
                                                     ! 1=gas/day, 2=gas/night
      LOGICAL,   INTENT( IN ) :: LIRRFLAG      ! Flag for IRR calculations
      INTEGER,   INTENT( INOUT ) :: NIRRCLS    ! No. of cells in block for IRR
      INTEGER,   INTENT( IN ) :: IRRCELL( : )  ! Cell No. of an IRR cell

C..Parameters:

c..ROS3 solver parameters - from KPP
      REAL( 8 ), PARAMETER :: GAM =  0.43586652150845899941601945119356D+00
      REAL( 8 ), PARAMETER :: C21 = -0.10156171083877702091975600115545D+01
      REAL( 8 ), PARAMETER :: C31 =  0.40759956452537699824805835358067D+01
      REAL( 8 ), PARAMETER :: C32 =  0.92076794298330791242156818474003D+01
      REAL( 8 ), PARAMETER :: B1  =  0.10000000000000000000000000000000D+01
      REAL( 8 ), PARAMETER :: B2  =  0.61697947043828245592553615689730D+01
      REAL( 8 ), PARAMETER :: B3  = -0.42772256543218573326238373806514D+00
      REAL( 8 ), PARAMETER :: D1  =  0.50000000000000000000000000000000D+00
      REAL( 8 ), PARAMETER :: D2  = -0.29079558716805469821718236208017D+01
      REAL( 8 ), PARAMETER :: D3  =  0.22354069897811569627360909276199D+00
      REAL( 8 ), PARAMETER :: A21 =  1.0D+00
      REAL( 8 ), PARAMETER :: A31 =  1.0D+00
      REAL( 8 ), PARAMETER :: A32 =  0.0D+00
      REAL( 8 ), PARAMETER :: G1  =  0.43586652150845899941601945119356D+00
      REAL( 8 ), PARAMETER :: G2  =  0.24291996454816804366592249683314D+00
      REAL( 8 ), PARAMETER :: G3  =  0.21851380027664058511513169485832D+01
      REAL( 8 ), PARAMETER :: GROW = 1.0D+00 / 3.00D+00
      REAL( 8 ), PARAMETER :: RGAM = 1.0D+00 / GAM

      REAL( 8 ), PARAMETER :: DTSTART = 0.5D-01  ! Starting time step (min)
      REAL( 8 ), PARAMETER :: DTMIN   = 1.0D-08  ! Min time step
      REAL( 8 ), PARAMETER :: DTMAX   = 1.0D+01  ! Max time step
      REAL( 8 ), PARAMETER :: UROUND  = 1.0D-18  ! Roundoff parameter

      REAL( 8 ), PARAMETER :: FACMAX  = 1.0D+01  ! Max time step factor
      REAL( 8 ), PARAMETER :: FACMIN  = 1.0D-01  ! Min time step factor
      REAL( 8 ), PARAMETER :: FACONE  = 1.0D+00  ! Time step fac of 1.0

      REAL( 8 ), PARAMETER :: CONMIN  = 1.0D-30  ! Min conc

C..External FUNCTIONS:
 
C..Local Variables:
      CHARACTER( 16 ), SAVE :: PNAME = 'RBSOLVER'  ! Procedure name
      CHARACTER( 96 ) :: XMSG = ' '
      LOGICAL, SAVE :: LFIRST = .TRUE. ! Flag for first call

      INTEGER I, J, N, JSPC       ! Loop indices
      INTEGER IDIAGBEG            ! Index of diagonal start in Jac array 
      INTEGER OFFSET              ! Cell offset for blcock
      INTEGER NCELL               ! Cell loop index
      INTEGER NRX                 ! Loop index for reactions
      INTEGER ISPOLD              ! Species index for old array order
      INTEGER NCALL_DEGRADE       ! WTH
      INTEGER IOS                 ! status

      LOGICAL LPASS               ! Flag for convergence achieved

      REAL(8),    ALLOCATABLE, SAVE :: CIRR (  :,: )  ! Species concs for IRR analysis
      REAL(8),    ALLOCATABLE, SAVE :: RKIRR(  :,: )  ! Rate constants for IRR analysis
      INTEGER,    ALLOCATABLE, SAVE :: DUMMY(  : )    ! Dummy array for IRR call

      REAL( 8 ) :: D
      REAL( 8 ) :: DT             ! Time step
      REAL( 8 ) :: DTCELL         ! Time step for each cell for IRR
      REAL( 8 ) :: DTINV          ! Inverse of time step
      REAL( 8 ) :: DTFAC          ! Time step scale factor
      REAL( 8 ) :: GDTINV         ! Inverse of gamma x time step
      REAL( 8 ) :: TNOW           ! Elapsed time at start of integration step
      REAL( 8 ) :: TEND           ! Elapsed time at end of integration step
      REAL( 8 ) :: DT_DEGRADE     ! WTH: Time step for degradation routine
      REAL( 8 ) :: YTOL           ! Species tolerance
      REAL( 8 ), SAVE :: RNSPEC   ! Recipricol of # of species
      REAL( 8 ) :: X1, X2         ! Temp ROS3 variables
      REAL( 8 ) :: ERRYMAX        ! Cell/species stiffness estimate
      REAL( 8 ) :: YLOWEPS        ! Tolerance ratio used in stiffness calc
      REAL( 8 ) :: MAXERR         ! Max of cell error estimates
      REAL( 8 ) :: OLDERR         ! OLD max error 

      INTEGER    COL_ERR          ! column for max of cell error estimate
      INTEGER    ROW_ERR          ! row for max of cell error estimate
      INTEGER    LAY_ERR          ! layer for max of cell error estimate
      INTEGER    CELL_MAXERR      ! cell with maximum error estimate
      REAL( 8 ) :: MAX_SPC_ERR      ! species error in the cell
      
      REAL( 8 ), ALLOCATABLE, SAVE :: YDOT( :,: )   ! dc/dt array

c.....ROS3 intermediate variables
      REAL( 8 ), ALLOCATABLE, SAVE :: K1( :,: )    
      REAL( 8 ), ALLOCATABLE, SAVE :: K2( :,: )
      REAL( 8 ), ALLOCATABLE, SAVE :: K3( :,: )
      REAL( 8 ), ALLOCATABLE, SAVE :: K4( :,: )
      REAL( 8 ), ALLOCATABLE, SAVE :: YP( :,: )     ! Predicted conc
      REAL( 8 ), ALLOCATABLE, SAVE :: ERR( : )      ! Error est for each cell


#ifdef rbdebug

      INTEGER COL, CD                       ! Column for debug output 
      INTEGER ROW, RD                       ! Row for debug output
      INTEGER LEV, LD                       ! Level for debug output
      INTEGER DBGOUT                        ! Output unit for debug output

      LOGICAL LDEBUG                       ! Debug output flag
      LOGICAL, SAVE :: LOPEN = .FALSE.     ! Flag for debug file opened

#endif

      INTERFACE
!         SUBROUTINE RBFEVAL( NCSP, YIN, YDOT )
!            INTEGER,   INTENT(    IN ) :: NCSP          ! Index of mech to use: 1=gas/day, 2=gas/night
!            REAL( 8 ), INTENT( INOUT ) :: YIN(  :, : )  ! Species concs, ppm
!            REAL( 8 ), INTENT(   OUT ) :: YDOT( :, : )  ! Species rates of change, ppm/min
!         END SUBROUTINE RBFEVAL
!         SUBROUTINE RBSOLVE( NCSP, RHS )
!           INTEGER,   INTENT( IN ) :: NCSP        ! Index of chem mech to use: 1=gas/day, 2=gas/night
!           REAL( 8 ), INTENT( INOUT ) :: RHS( :,: )  ! Right hand side = {b}
!         END SUBROUTINE RBSOLVE
!         SUBROUTINE RBJACOB( NCSP, YIN )
!           INTEGER,   INTENT( IN ) :: NCSP        ! Index of chem mech to use; 1=gas/day, 2=gas/night
!           REAL( 8 ), INTENT( IN ) :: YIN( :,: )  ! Species concs, ppm
!         END SUBROUTINE RBJACOB
!         SUBROUTINE RBDECOMP( NCSP )
!           INTEGER, INTENT( IN ) :: NCSP  ! Index of chem mech to use: 1=gas/day, 2=gas/night
!         END SUBROUTINE
          SUBROUTINE KERNELDRIVER( BLKSIZE, Y, RKI, DTSTART,
     &                   CHEMSTEP, NCSP, LORDERING, ISCHAN,   
     &                   IDIAGBEG, DT, TNOW, LPASS, OFFSET,   
     &                   ATOL, RTOL, ERRMX2, RGAM, IARRAY,    
     &                   NUMCELLS, N_SPEC, C21, C31, 
     &                   C32, B1, B2, B3, D1, D2, D3, UROUND, 
     &                   GROW, CONMIN, FACMIN, FACMAX, 
     &                   FACONE, DTMAX, DTMIN, MXCELLS, NCS2,  
     &                   NUMB_MECH_SPC, RNSPEC, NUSERAT,      
     &                   NKUSERAT, NREACT, IRM2, NPRDCT, SC,  
     &                   NRXNS, MXRCT, MXPRD, MXRR, MXRP, NDERIVL, 
     &                   JARRL, JLIAL, NDERIVP, JARRP, JPIAL, ICOEFF, JZLO, IDEC1LO, 
     &                   IDEC1HI, MXCOUNT2, IJDECA, IJDECB, 
     &                   IKDECA, IKDECB, KJDECA, KJDECB, JARRAYPT, 
     &                   JHIZ1, JHIZ2, MXARRAY, JZEROA, JZEROB, 
     &                   KZHI0, KZERO, KZLO1, KZHI1, KZLO2, 
     &                   KZILCH, MZHI0, MZLO1, MZHI1, MZLO2, MZILCH, NCS, MPI_COMM_WORLD,
     &                   BLKTEMP, BLKDENS, NSPECIAL_RXN )
            INTEGER, INTENT( IN ) :: BLKSIZE
            REAL( 8 ), INTENT( INOUT ) :: Y( BLKSIZE,NUMB_MECH_SPC )
            REAL( 8 ), INTENT( IN ) :: RKI( BLKSIZE,NRXNS )
            REAL( 8 ), INTENT( IN ) :: DTSTART
            REAL( 8 ), INTENT( IN ) :: CHEMSTEP
            INTEGER, INTENT( IN ) :: NCSP
            LOGICAL, INTENT( IN ) :: LORDERING
            INTEGER, INTENT( IN ) :: ISCHAN
            INTEGER, INTENT( IN ) :: IDIAGBEG
            REAL( 8 ), INTENT( IN ) :: DT
            REAL( 8 ), INTENT( IN ) :: TNOW
            LOGICAL, INTENT( IN ) :: LPASS
            INTEGER, INTENT( IN ) :: OFFSET
            REAL( 8 ), INTENT( IN ) :: ATOL( NUMB_MECH_SPC )
            REAL( 8 ), INTENT( IN ) :: RTOL( NUMB_MECH_SPC )
            REAL( 8 ), INTENT( INOUT ) :: ERRMX2( MXCELLS )
            REAL( 8 ), INTENT( IN ) :: RGAM
            INTEGER, INTENT( IN ) :: IARRAY( NCS2 )
            !INTEGER, INTENT( IN ) :: BLKSIZE
            INTEGER, INTENT( IN ) :: NUMCELLS
            INTEGER, INTENT( IN ) :: N_SPEC
            REAL( 8 ), INTENT( IN ) :: C21
            REAL( 8 ), INTENT( IN ) :: C31
            REAL( 8 ), INTENT( IN ) :: C32
            REAL( 8 ), INTENT( IN ) :: B1
            REAL( 8 ), INTENT( IN ) :: B2
            REAL( 8 ), INTENT( IN ) :: B3
            REAL( 8 ), INTENT( IN ) :: D1
            REAL( 8 ), INTENT( IN ) :: D2
            REAL( 8 ), INTENT( IN ) :: D3
            REAL( 8 ), INTENT( IN ) :: UROUND
            REAL( 8 ), INTENT( IN ) :: GROW
            REAL( 8 ), INTENT( IN ) :: CONMIN
            REAL( 8 ), INTENT( IN ) :: FACMIN
            REAL( 8 ), INTENT( IN ) :: FACMAX
            REAL( 8 ), INTENT( IN ) :: FACONE
            REAL( 8 ), INTENT( IN ) :: DTMAX
            REAL( 8 ), INTENT( IN ) :: DTMIN
            INTEGER, INTENT( IN ) :: MXCELLS
            INTEGER, INTENT( IN ) :: NCS2
            INTEGER, INTENT( IN ) :: NUMB_MECH_SPC
            REAL( 8 ), INTENT( IN ) :: RNSPEC
            INTEGER, INTENT( IN ) :: NUSERAT( NCS2 )
            INTEGER, INTENT( IN ) :: NKUSERAT( NRXNS,NCS2 )
            INTEGER, INTENT( IN ) :: NREACT( NRXNS )
            INTEGER, INTENT( IN ) :: IRM2( NRXNS,MXRCT+MXPRD,NCS2 )
            INTEGER, INTENT( IN ) :: NPRDCT( NRXNS )
            REAL( 8 ), INTENT( IN ) :: SC( NRXNS,MXPRD )
            INTEGER, INTENT( IN ) :: NRXNS
            INTEGER, INTENT( IN ) :: MXRCT
            INTEGER, INTENT( IN ) :: MXPRD
            INTEGER, INTENT( IN ) :: MXRR
            INTEGER, INTENT( IN ) :: MXRP
            INTEGER, INTENT( IN ) :: NDERIVL( NRXNS,NCS2 )
            INTEGER, INTENT( IN ) :: JARRL( NRXNS,MXRR,NCS2 )
            INTEGER, INTENT( IN ) :: JLIAL( NRXNS,MXRR,NCS2 )
            INTEGER, INTENT( IN ) :: NDERIVP( NRXNS,NCS2 )
            INTEGER, INTENT( IN ) :: JARRP( NRXNS,MXRP,NCS2 )
            INTEGER, INTENT( IN ) :: JPIAL( NRXNS,MXRP,NCS2 )
            INTEGER, INTENT( IN ) :: ICOEFF( NRXNS,MXRP,NCS2 )
            INTEGER, INTENT( IN ) :: JZLO( NCS2 )
            INTEGER, INTENT( IN ) :: IDEC1LO( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: IDEC1HI( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: MXCOUNT2
            INTEGER, INTENT( IN ) :: IJDECA( MXCOUNT2 )
            INTEGER, INTENT( IN ) :: IJDECB( MXCOUNT2 )
            INTEGER, INTENT( IN ) :: IKDECA( MXCOUNT2 )
            INTEGER, INTENT( IN ) :: IKDECB( MXCOUNT2 )
            INTEGER, INTENT( IN ) :: KJDECA( MXCOUNT2 )
            INTEGER, INTENT( IN ) :: KJDECB( MXCOUNT2 )
            INTEGER, INTENT( IN ) :: JARRAYPT( NUMB_MECH_SPC,NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: JHIZ1( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: JHIZ2( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: MXARRAY
            INTEGER, INTENT( IN ) :: JZEROA( MXARRAY )
            INTEGER, INTENT( IN ) :: JZEROB( MXARRAY )
            INTEGER, INTENT( IN ) :: KZHI0( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: KZERO( MXARRAY,NCS2 )
            INTEGER, INTENT( IN ) :: KZLO1( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: KZHI1( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: KZLO2( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: KZILCH( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: MZHI0( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: MZLO1( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: MZHI1( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: MZLO2( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: MZILCH( NUMB_MECH_SPC,NCS2 )
            INTEGER, INTENT( IN ) :: NCS
            INTEGER, INTENT( IN ) :: MPI_COMM_WORLD
            REAL( 8 ), INTENT( IN ) :: BLKTEMP( BLKSIZE )
            REAL( 8 ), INTENT( IN ) :: BLKDENS( BLKSIZE )
            INTEGER, INTENT( IN ) :: NSPECIAL_RXN  
          END SUBROUTINE KERNELDRIVER
      END INTERFACE

C**********************************************************************

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C The s-stage Rosenbrock method solves the following equations
C                        s
C         Y(n+1) = Yn + SUM( Bi * Ki) 
C                       i=1
C     where
C                          i-1                          i
C         Ki = H * F[ Yn + SUM( Aij * Kj ) ] + H * J * SUM( GAMMAij * Kj ) 
C                          j=1                         j=1
C
C See Sandu et al. for details and the values of Bi, GAMMAij, Aij, etc.
C
C For computational efficiency, the equations are re-arranged as  
C follows (e.g., see Press, Numerical Recipes, Chap 16 on ODEs )
C
C     Gi = SUM( GAMMAij * Kj ) + GAMMA * Ki  i=1,...s
C
C    [ I / GAMMA H - J ] G1 = F[ Yn ]
C    [ I / GAMMA H - J ] G2 = F[ Yn + A21 * G1 ] + ( C21 * G1 ) / H
C    [ I / GAMMA H - J ] G3 = F[ Yn + A31 * G1 + A32 * G2 ] + 
C                              ( C31 * G1 + C32 G2 ) / H   
C 
C The code below sequentially calculates the Gi, and then computes
C Y(n+1) via the first formula.  Note that the parameter values will be
C different from those shown in Sandu et al. because the code computes
C Gi instead of Ki
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF ( LFIRST ) THEN
         LFIRST = .FALSE.
         ALLOCATE ( YDOT( BLKSIZE,N_SPEC ),
     &              K1  ( BLKSIZE,N_SPEC ),
     &              K2  ( BLKSIZE,N_SPEC ),
     &              K3  ( BLKSIZE,N_SPEC ),
     &              K4  ( BLKSIZE,N_SPEC ),
     &              YP  ( BLKSIZE,N_SPEC ),
     &              ERR ( BLKSIZE ), STAT = IOS )
         IF ( IOS .NE. 0 ) THEN
            XMSG = '*** Memory Allocation Error'
            CALL M3EXIT( PNAME, 0, 0, XMSG, XSTAT1 )
         END IF

!        IF ( LIRRFLAG ) THEN   !*** This works only if in irr subdomain window
            ALLOCATE ( CIRR  ( BLKSIZE,N_SPEC ),
     &                 RKIRR ( BLKSIZE,N_RXNS ), 
     &                 DUMMY ( BLKSIZE ),         STAT = IOS )
            IF ( IOS .NE. 0 ) THEN
               XMSG = '*** Memory Allocation Error'
               CALL M3EXIT( PNAME, 0, 0, XMSG, XSTAT1 )
            END IF
!        END IF

         RNSPEC = 1.0D0 / REAL( N_SPEC,8 )

      END IF   ! LFIRST

c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++
#ifdef rbdebug 

! CD, RD, LD set by modifying code
      CD = 0
      RD = 0
      LD = 0
      DO NCELL = 1, NUMCELLS
         COL = CCOL( NCELL )
         ROW = CROW( NCELL )
         LEV = CLEV( NCELL )
         IF ( CD .EQ. COL .AND. RD .EQ. ROW .AND. LD .EQ. LEV ) THEN
!        IF ( JTIME .EQ. 160000 ) THEN
              LDEBUG = .TRUE.
              EXIT
         ELSE
              LDEBUG = .FALSE.
         END IF
       END DO

         IF ( LDEBUG ) THEN
              IF ( .NOT. LOPEN ) THEN
                 DBGOUT = JUNIT()
                 OPEN( UNIT = DBGOUT, FILE = 'debug.out' )
                 LOPEN = .TRUE.
              END IF

              WRITE( DBGOUT, '( A, 2I4, I3, 1X, I7, 1X, I6 ) ' )
     &              'Debug output for col/row/lev/date/time:', 
     &               C, R, L, JDATE, JTIME
              WRITE( DBGOUT, '( A, F7.2) ' )
     &              'CHEMTMSTEP = ', CHEMSTEP
              WRITE( DBGOUT, '( A )' ) 'Starting concs and rate constants'
              DO N = 1, N_SPEC
                 WRITE( DBGOUT,  '( A, I3, 1X, A, 1X, 1PE13.5 )' )
     &                           'SP ',N, CHEMISTRY_SPC( N ), Y( NCELL, N )
              END DO
              DO N = 1, N_RXNS
                 WRITE( DBGOUT, '( A, I3, 1X, 1PE13.5 )' )
     &                          'RKI ', N, RKI( NCELL, N )
              END DO
        END IF

#endif
c++++++++++++++++++++++++Debug section++++++++++++++++++++++++++++++++++


      ISCHAN = ISCHANG( NCS )

      IDIAGBEG = IARRAY( NCSP ) - ISCHAN + 1

      DT = MAX( DTMIN, DTSTART )

      TNOW = 0.0D+00

      NCALL_DEGRADE = 0

      LPASS = .FALSE.

      OFFSET = BLKCNO( BLKID )

      IF ( LIRRFLAG ) THEN

         DO NRX = 1, N_RXNS
            DO NCELL = 1, NIRRCLS
               RKIRR( NCELL,NRX ) = RKI( IRRCELL( NCELL ),NRX )
            END DO
         END DO

         DTCELL = 0.0D0
         DO I = 1, ISCHAN
            ISPOLD = INEW2OLD( I, NCS )
            DO NCELL = 1, NIRRCLS
               CIRR( NCELL, ISPOLD ) = Y( IRRCELL( NCELL ),I )
            END DO
         END DO
         CALL PA_IRR ( .TRUE., .FALSE., RKIRR, CIRR, DTCELL, NIRRCLS, DUMMY )

      END IF


#ifdef sens
       YAVE = 0.0D0
#endif

      !PRINT *, "NKUSERAT ", SIZEOF( NKUSERAT )
      !PRINT *, "NCS2 ", NCS2

      CALL KERNELDRIVER( BLKSIZE, Y, RKI, DTSTART,
     &                   CHEMSTEP, NCSP, LORDERING, ISCHAN,   
     &                   IDIAGBEG, DT, TNOW, LPASS, OFFSET,   
     &                   ATOL, RTOL, ERRMX2, RGAM, IARRAY,    
     &                   NUMCELLS, N_SPEC, C21, C31, 
     &                   C32, B1, B2, B3, D1, D2, D3, UROUND, 
     &                   GROW, CONMIN, FACMIN, FACMAX, 
     &                   FACONE, DTMAX, DTMIN, MXCELLS, NCS2,  
     &                   NUMB_MECH_SPC, RNSPEC, NUSERAT,      
     &                   NKUSERAT, NREACT, IRM2, NPRDCT, SC,  
     &                   NRXNS, MXRCT, MXPRD, MXRR, MXRP, NDERIVL, 
     &                   JARRL, JLIAL, NDERIVP, JARRP, JPIAL, ICOEFF, JZLO, IDEC1LO, 
     &                   IDEC1HI, MXCOUNT2, IJDECA, IJDECB, 
     &                   IKDECA, IKDECB, KJDECA, KJDECB, JARRAYPT, 
     &                   JHIZ1, JHIZ2, MXARRAY, JZEROA, JZEROB, 
     &                   KZHI0, KZERO, KZLO1, KZHI1, KZLO2, 
     &                   KZILCH, MZHI0, MZLO1, MZHI1, MZLO2, MZILCH, NCS, MPI_COMM_WORLD,
     &                   BLKTEMP, BLKDENS, NSPECIAL_RXN ) 


#ifdef sens
!   complete calculation for YAVE
       YAVE(1:NUMCELLS,1:NUMB_MECH_SPC) = YAVE(1:NUMCELLS,1:NUMB_MECH_SPC) / CHEMSTEP
!        YAVE = 0.5D0*(YAVE+Y)
#endif
      RETURN

92100 FORMAT( '      Convergence failure ', 
     &       '  JDATE = ', I7, '  JTIME = ' , I6 )

92110 FORMAT( ' Convergence failure in Gas Chemistry Solver ', 
     &       '  JDATE = ', I7, '  JTIME = ' , I6,
     &       ' at COL = ', I4, ' ROW = ', I4, ' LAY = ', I4 ) 

92113 FORMAT( 'Number of Species solved = ', I4, ' out of ',
     &         I4, ' Total GC Species ')

92114 FORMAT('CHEMISTRY_SPC( ', I4, ' ) = ', A16,' Error Contribution =', 
     &        ES12.4, ' Initial Conc = ', ES12.4 , 
     &        ' Predicted Conc = ', ES12.4 )
     
92115 FORMAT(' Rates used in Photolysis Reactions ')

92116 FORMAT('Reaction: ', A16, ' uses PHOTAB ', A16, ' = ', ES12.4)

92117 FORMAT('Reaction: ', A16, ' uses Reaction ', A16, ' = ', ES12.4)

92216 FORMAT(/ 'MAXERR = ', ES12.4
     &       / 'SQRT( ERR( CELL_MAXERR ) * RNSPEC ) = ', ES12.4)

92214 FORMAT( I3, A16,' = ', E12.4)

92215 FORMAT(/ 'Cell Properties '
     &       / 'Temp = ', ES12.4, ' K '
     &       / 'Press = ', ES12.4,' Pa '
     &       / 'Dens = ', ES12.4, ' Kg/m3 ' 
     &       / 'H2O Vapor = ', ES12.4, ' ppm '
     &       / 'SEAWATER( - ) = ', ES12.4 )
      
      END

