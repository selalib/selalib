!-----------------------------------------------------------------------
!**BEGIN PROLOGUE Input
!**DATE WRITTEN   1999/06/01  (year/month/day)
!**AUTHORS
!
!  Peter Kravanja, Marc Van Barel
!
!    Department of Computer Science
!    Katholieke Universiteit Leuven
!    Celestijnenlaan 200 A
!    B-3001 Heverlee
!    BELGIUM
!    E-mail: Peter.Kravanja@na-net.ornl.gov
!            Marc.VanBarel@cs.kuleuven.ac.be 
!
!  Omiros Ragos, Michael N. Vrahatis, and Filareti A. Zafiropoulos
!
!    Department of Mathematics
!    University of Patras
!    GR-261.10 Patras
!    GREECE
!    E-mail: {ragos,vrahatis,phikapa}@math.upatras.gr
!
!**DESCRIPTION
!
!  The user can specify the values of the input parameters LV and H,
!  which determine the geometry of the rectangular region that is to be
!  considered, the values of M, ICON and NR, which determine the type
!  of calculation that ZEAL will perform, and the values of VERBOSE,
!  FILES and IFAIL. 
!
!**END PROLOGUE Input
!-----------------------------------------------------------------------

MODULE Zeal_Input_Module

     USE Precision_Module

     IMPLICIT NONE

!-----------------------------------------------------------------------
!**ACCESSIBILITY
!
     PUBLIC

!-----------------------------------------------------------------------
!  Specify the geometry of the rectangular region. 
!  The lower left vertex has coordinates LV(1) and LV(2). The stepsizes 
!  along the X- and Y-direction are H(1) and H(2), respectively.
!

     REAL(KIND=DP), DIMENSION(2) :: LV = (/-2.0_DP,-2.0_DP/)
     REAL(KIND=DP), DIMENSION(2) :: H  = (/ 4.0_DP, 4.0_DP/)


!     REAL(KIND=DP), DIMENSION(2), PARAMETER :: LV = (/-2.0_DP,-2.0_DP/)
!     REAL(KIND=DP), DIMENSION(2), PARAMETER :: H  = (/ 4.0_DP, 5.0_DP/)

!
!  Other examples
!
!    REAL(KIND=DP), DIMENSION(2), PARAMETER :: LV = (/-0.5_DP,-0.5_DP/)
!    REAL(KIND=DP), DIMENSION(2), PARAMETER :: H  = (/ 6.0_DP, 2.0_DP/)

!    REAL(KIND=DP), DIMENSION(2), PARAMETER :: LV = (/-1.0_DP,-1.0_DP/)
!    REAL(KIND=DP), DIMENSION(2), PARAMETER :: H  = (/ 4.0_DP, 2.0_DP/) 
     
!-----------------------------------------------------------------------
!  The value of M determines the maximum number of zeros 
!  (counting multiplicities) that are considered within a subregion.
!  M has to be larger than the maximum of the multiplicities of the
!  zeros. A recommended value is 5.
!  --------
!  There are conditionning problems for plasma dispersion relations
!  when M=5, M=3 seems to be a robust choice (Eric Sonnendrucker)
     INTEGER, PARAMETER :: M = 3

!-----------------------------------------------------------------------
!  Specify the values of ICON and NR.
!  The parameter ICON may take the values 1, 2, 3 and 4.   
!
!    1  =>  calculation of the total number of zeros, only.
!
!    2  =>  calculation of the total number of zeros and isolation 
!           of subregions that contain at most M zeros.
!
!    3  =>  calculation of the total number of zeros, isolation of 
!           subregions that contain at most M zeros, and computation 
!           of all the zeros and their multiplicity.
!
!    4  =>  calculation of the total number of zeros, isolation of 
!           subregions that contain at most M zeros, and computation 
!           of NR mutually distinct zeros and their multiplicity.
!
     INTEGER, PARAMETER :: ICON = 3 
     INTEGER            :: NR   = 2 

!-----------------------------------------------------------------------
!  Specify the value of VERBOSE.
!  ZEAL is allowed to print information (concerning the user's input 
!  and the computed results) if and only if VERBOSE is equal to .TRUE. 
!
     LOGICAL, PARAMETER :: VERBOSE = .FALSE.

!-----------------------------------------------------------------------
!  Specify the value of FILES.
!  If FILES is set equal to .TRUE. then ZEAL generates the files 
!  "zeros.dat" and "mult.dat". They contain the computed approximations 
!  for the zeros as well as their respective multiplicities. ZEAL also 
!  writes the file "fzeros.dat", which contains the values that the 
!  function takes at the computed approximations for the zeros.
!
     LOGICAL :: FILES = .FALSE.

!-----------------------------------------------------------------------
!  Specify the value of IFAIL.
!  This variable determines how errors are to be handled. 
!  We follow the NAG convention:
!
!    1  =>  "soft silent error"
!           Control returned to calling program.
!
!   -1  =>  "soft noisy error"
!           Error message is printed.
!           Control returned to calling program.
!
!    0  =>  "hard noisy error"
!           Error message is printed and program is stopped.
!
     INTEGER :: IFAIL = 1 


   END MODULE Zeal_Input_Module
