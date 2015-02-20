!File: Module Zone
!Module commun contenant notammment 
!les types derivees et les donnees physiques du probleme.
module zone

implicit none

type mesh_fields
   real(8), dimension(:,:), pointer :: e, b, j
   real(8), dimension(:,:,:), pointer :: j_dg   ! ajout Martin
   real(8), dimension(:,:), pointer :: r_dg     ! ajout Pierre
end type mesh_fields


integer, private :: i, j

!Constant: mitmax
!Nombre maximum d'iterations 
integer :: mitmax = 1000000     

!Constant: nesmx
!nombre max de types de particules
integer, parameter :: nesmx = 5 
!Constant: pi 
!constante mathematique
real(8), parameter :: pi = 3.14159265358979323846264338327950d0
!Constants: petitx, grandx
!petitx - petite valeur absolue    
!grandx - grande valeur absolue   
real(8), parameter :: petitx = 1.e-20, grandx = 1.e+20
!Constant: xnbdig
!plus petit nb en flottant sans exposant
real(8), parameter :: xnbdig = 1.e-08
!Constant: c
!Vitesse de la lumiere et vitesse de la lumiere au carre
real(8) :: c    = 2.9979e+8, csq
!Constant qelec
!Charge de l'electron
real(8), parameter :: qelec  = 1.6021d-19
!Constant: eps0     
!Permitivite electrique du vide
real(8) :: eps0 = 8.8542e-12
!Constant: xmu0     
!Permeabilite magnetique du vide
real(8) :: xmu0 = 4.e-7 * pi
!logical: lmodte, lmodtm
!resolution du systeme (e1,e2,b3) ou (b1,b2,e3)
logical :: lmodte  = .true., lmodtm = .false.

!integer: nschme
!Type de schema temporel utilise en mode TE
!0 - saute-mouton
!1 - Newmark      
!2 - schema a 3 pas
integer :: nschme  = 0
!integer: nschmm
!Type de schema temporel utilise en mode TM
!0 - saute-mouton
!1 - Newmark      
!2 - schema a 3 pas
integer :: nschmm = 0
!real: epsche, epschm
!coefficient pour le schema temporel mode TE, TM
real(8) :: epsche  = 0.0, epschm = 0.0

!logical: lcorrp
!correction de Poisson 
logical :: lcorrp  = .true.

logical :: lchexn  = .false., lchexp = .false.

integer :: nstep
integer :: nesp = 0     !nombre d'especes de particules
integer :: iout

!real: dt
!pas de temps
real(8) :: dt, time = 0.0

!character: titre
!Nom du directory devant contenir les sorties
character(len=40)  :: titre

!Donnees communes pour le calcul avec le schema de Yee
type tm_mesh_fields
   real(8), dimension(:,:), pointer :: ex, ey, bz
   real(8), dimension(:,:), pointer :: r0, r1
   real(8), dimension(:,:), pointer :: jx, jy
end type tm_mesh_fields

!*** an electron
real(8) :: rho_0,E_0
character(len=7) :: solver
character(len=6) :: nomcas, nomcor 
character(len=6) :: bcname, jname
integer          :: nx, ny, icrea, idiag, nbpm
! pour toutes les corrections
real(8) :: chi, cfl, omega
! Champs externes pour l'utilisation du solveur Yee
integer :: shapeB
real(8) :: S, eta
real(8) :: exext, eyext, bzext
real(8) :: kx, ky,  poids, alpha, q_sur_m
logical :: relativ

! parametres pour le cas statio
integer :: mm=4, nn=2

logical :: ldtfrc = .false. !permet de forcer le calcul si dt trop grand 

! minter : nb d'iterations dans la boucle interne
! mcentr : nb de boucles internes par boucle centrale
! mexter : nb de boucles centrales par boucle externe 

integer :: minter, mexter, mcentr
integer :: nwrpar        !Niveau d'impression des diagnostics 

!----------------------------------------------------------------------------!
logical :: lctrac = .true.   !Activateur general des traces en ligne 
logical :: ldebug = .false.  !Active les sorties pour le debuggage
logical :: landau = .false.  !Amortissement Landau
logical :: lmtest = .false.  !Cas Test Maxwell 
logical :: lcopic = .false.  !Methode conservative du calcul des courants

real(8) :: gsq = 0d0 !(gamma au carre)  !

real(8) :: xO, yO   !Rayon du faisceau
real(8) :: vth  	!Vitesse thermique
real(8) :: vbeam  !Vitesse du faisceau

end module zone
