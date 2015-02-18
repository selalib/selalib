!File: Module Zone
!Module commun contenant notammment 
!les types derivees et les donnees physiques du probleme.
module zone

implicit none

integer, private :: i, j

!Type: particle
!Structure pour une espece de particule 
!nbpart - Nombre de particules
!nlpa   - Index de l'element ou se trouve la particule
!pos    - Position
!vit    - Vitesse
!poid   - Poids
!date   - Date d'emission
!epx    - Composante suivant x du champ electrique
!epy    - Composante suivant y du champ electrique
!epz    - Composante suivant z du champ electrique
!bpx    - Composante suivant x du champ magnetique
!bpy    - Composante suivant y du champ magnetique
!bpz    - Composante suivant z du champ magnetique
!xlm    - Coordonnees dans l'element

!Type: objet_fonction
!Type derive d'une forme temporelle
!ifor - Nature de la fonction
!tdb  - Instant de debut
!tmt  - Duree de la montee
!tpl  - Duree du plateau
!tds  - Duree de la descente
!freq - Frequence
!phas - Phase
!sig  - coeffiscient
!tau  - Autre coeffiscient
!eta  - Autre coeffiscient

type particle
   integer :: nbpart
   integer, pointer :: nlpa(:)
   real(4), pointer :: fnlpa(:)
   real(8), pointer :: pos(:,:), vit(:,:), poid(:), date(:)
   real(8), pointer :: epx(:), epy(:), epz(:)
   real(8), pointer :: bpx(:), bpy(:), bpz(:)
   real(8), pointer :: xlm(:,:), xlb(:,:)
end type particle

type mesh_fields
   real(8), dimension(:,:), pointer :: e, b, j
   real(8), dimension(:,:,:), pointer :: j_dg   ! ajout Martin
   real(8), dimension(:,:), pointer :: r_dg     ! ajout Pierre
end type mesh_fields

type objet_fonction
   integer :: ifor
   real(8)    :: tmt, tds, tdb, tpl
   real(8)    :: freq, phas, sig, tau, eta
   integer, dimension(:),   pointer :: nbpdon, ippdon
   real(8), dimension(:,:), pointer :: ttpdon, vtpdon 
end type objet_fonction

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

!integer: nbpam 
!nombre max de particules par espece
integer :: nbpam(nesmx) 
!
!integer: nbpmax 
!nombre max de particules pour toutes especes
integer :: nbpmax

!real: pcharg
!charge elementaire d'une particule
real(8) :: pcharg(nesmx)

!real: pmasse
!masse  elementaire d'une particule
real(8) :: pmasse(nesmx)

!real: rapbol
real(8) :: rapbol(nesmx)
!rapport de Boltzmann

!real: alprjt
!Coefficient de compensation de charge d'espace
!par particule a appliquer a la projection
real(8) :: alprjt(nesmx)

!logical: lpousp
!Activateur general du pousseur
logical :: lpousp

!logical: lpoust
!Activateur du pousseur par espece
logical :: lpoust(nesmx)

!logical: lchapa
!Activateur general de l'interpolation
logical :: lchapa

!logical: lchapt
!Activateur de l'interpolation par espece
logical :: lchapt(nesmx)

!logical: lprojp
!Activateur general de l'assignation 
logical :: lprojp

!logical: lprojt
!Activateur de l'assignation par espece
logical :: lprojt(nesmx)

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
integer :: Mx, My
real(8) :: dx, dy, dimx, dimy, tfinal
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
logical :: lctrac = .true. 	!Activateur general des traces en ligne 
logical :: ldebug = .false.	!Active les sorties pour le debuggage
logical :: landau = .false.	!Amortissement Landau
logical :: lmtest = .false.	!Cas Test Maxwell 
logical :: lcopic = .false.	!Methode conservative du calcul des courants

real(8) :: gsq = 0d0 !(gamma au carre)	!

real(8) :: xO, yO 	!Rayon du faisceau
real(8) :: vth		!Vitesse thermique
real(8) :: vbeam	!Vitesse du faisceau

contains

subroutine init_chem(nsolve, tm, tm1, inpfil)

implicit none

integer, intent(inout) :: nsolve
type(tm_mesh_fields) :: tm, tm1, sol
!character(len=72) :: s
integer :: ios
character(len=*) :: inpfil

namelist/nlchem/ dimx,   &  !largeur
                  dimy,  &  !longueur
                    nx,  &  !nbre de points suivant x
                    ny,  &  !nbre de points suivant y
                  nbpm,  &  !nbre de particules par maille
                   cfl,  &  !nbre de Courant
                tfinal,  &  !duree maxi
                nomcas,  &  !nom du cas ce calcul
                solver,  &  !poisson ou maxwell
                nomcor,  &  !type de correction
                 jname,  &  !calcul de j   
         	icrea,   &  !frequence d'emission des particules
         	idiag,   &  !frequence des diagnostics
         	bcname,  &  !type de conditions limites
		shapeB,  &  !Type de focalisation 
		S, eta,  &  !Type de focalisation 
		vbeam,   &  !Vitesse de deplacement du faisceau
         	exext,   &  !champ electrique exterieur
         	eyext,   &  !champ electrique exterieur
         	bzext,   &  !champ magnetique exterieur
         	poids,   &  !poids total
                relativ, &  !calcul relativiste de la vitesse
		xO, yO, vth !Quantites pour les faisceaux

rho_0   = 1.d0
E_0     = 1.d0
relativ = .true.
nbpm    = 50
csq     = c * c
xO	= 0.0
yO	= 0.0
vth	= 0.0
vbeam   = 0.0
shapeB  = 0; S = 0.0; eta = 0.0
exext   = 0.0; eyext = 0.0; bzext = 0.0

open(10, file=inpfil, status="old")
write(*,*)"Lecture de la namelist $nlchem"
read(10,nlchem)
close(10)

if (nomcas == "envelo") then
   write(*,*) " xO, yO = ", xO, yO
   if (xO == 0. .or. yO == 0.) then
      stop ' Rayon du faisceau nul ! '
   end if
   if (vth == 0.) then
      stop ' Vitesse thermique du faisceau nulle ! '
   end if
end if
if (nomcas == "plasma") then
   alpha = 0.1
   kx = 0.5; ky = 0.
   dimx = 2*pi/kx
   poids = dimx * dimy ! car int(f0) = dimx*dimy
endif

if (nomcas == "weibel") poids = dimx * dimy

dx = dimx / (nx-1)
dy = dimy / (ny-1)
if (solver == 'maxwell') &
dt = cfl  / sqrt (1./(dx*dx)+1./(dy*dy)) / c

if ( nstep*dt > tfinal ) nstep = floor(tfinal/dt)

!if (nomcor == 'marder') dt = 0.5*cfl /(1./(dx*dx)+1./(dy*dy)) /csq
if (nomcor == 'hyperb') then
   chi = 1.2
   dt  = dt/chi
endif

if (solver=='poisson') then 
   nsolve = 21
else if (solver=='maxwell') then
   nsolve = 22
end if

write(*,"( 10x,'Poids des particules = ',g15.3)") poids 
write(*,"( 10x,'Nombre de Courant    = ',g15.5)") cfl
write(*,"( 10x,'Nom du cas           = ',a)") nomcas
write(*,"(/10x,'Largeur  dimx        = ',g15.3)") dimx
write(*,"( 10x,'Longueur dimy        = ',g15.3)") dimy
write(*,"( 10x,'Nombre nx            = ',i3)"   ) nx
write(*,"( 10x,'Nombre ny            = ',i3)"   ) ny
write(*,"(/10x,'Solveur de             ',a)")  solver
write(*,"( 10x,'Type de correction   = ',a)") nomcor
write(*,"( 10x,'Calcul de j          = ',a)")  jname
write(*,"( 10x,'Frequence d''emission des particules',i3)") icrea
write(*,"( 10x,'Frequence des diagnostics',i3)") idiag  
write(*,"( 10x,'Type de conditions limites',a)") bcname 
write(*,"( 10x,'Champ electrique exterieur',g15.3)") exext
write(*,"( 10x,'Champ electrique exterieur',g15.3)") eyext
write(*,"( 10x,'Champ magnetique exterieur',g15.3)") bzext
write(*,"(/10x,'dx = ',g15.3)") dx
write(*,"(/10x,'dy = ',g15.3)") dy
write(*,"(/10x,'Pas de temps dt    = ',g15.3)") dt
write(*,"( 10x,'Temps              = ',g15.3) ") tfinal
write(*,"( 10x,'Nombre de particules par maille = ',g15.3) ") nbpm
write(*,"( 10x,'Nombre de Courant  = ',g15.3)") cfl
write(*,"( 10x,'Charge de la particule',g15.3)") pcharg(1)
write(*,"( 10x,'Masse de la particule ',g15.3)") pmasse(1)
q_sur_m = pcharg(1) / pmasse(1)
write(*,"( 10x,'q_sur_m =',g15.3) ") q_sur_m    
write(*,"( 10x,'Cas relativiste = ',g15.3)") relativ

end subroutine init_chem


!subroutine: readin
!Lecture des donnees systeme de la namelist "nlsyst" 
subroutine readin (prefix, dirpr, nsolve)

character(len=*), intent(in)   :: prefix
character(len=72), intent(out) :: dirpr
character(len=72) :: inpfil, expfil
integer, intent(out) :: nsolve
integer :: ios
!Solveur utilise
!0 - on n'avance pas les champs
!1 - solveur electrostatique
!2 - Maxwell volumes finis
!7 - Maxwell Poisson couples

!Anciennes variables utilisees pour M2V
logical :: ldemar = .false.     !Demarrage d'un cas uniquement 
logical :: lpospr = .false.     !Utilisation en post-processeur
logical :: lprote = .false.     !Protection apres integration 
logical :: lrepri = .false.     !Reprise apres protection 
integer :: nmaill = 0       !numero de mailleur ayant servi a creer le maillage 
logical :: lstmai = .false. !Stockage des tableaux du maillage 
real    :: pasdt  = 0.0     !Pas de temps
integer :: metpoi = 0       !methode pour Poisson
real    :: errpoi = 0.0     !max du carre de l'erreur pour POISSON
integer :: metamp = 0       !methode pour Poisson
real    :: erramp = 0.0     !max du carre de l'erreur pour POISSON
real    :: pivamp = 0.      !seuil de precision pour AMPERE 
integer :: nitgc  = 0       !nombre max d'iterations du gradient conjugue
real    :: pivpoi = 0.      !seuil de precision pour POISSON 
logical :: lsecur = .false. !Protection de securite
integer :: numprt = 0       !Numero de la protection a lire
real    :: tlmcpu = 1.e10   !Temps limite CPU associe a lsecur
logical :: lrztrk = .false. !remise a zero des traceurs
logical :: lrzsmp = .false.     !remise a zero des echantillons 
logical :: lrzhst = .false. !remise a zero des diag de  m2v.hst.nomj
logical :: lcputm = .false. !Calcul des temps CPU par theme 
logical :: lmainv = .false. !inversion du maillage (r,z->r,r)
integer :: iter = 1
real :: temps = 0.0

integer :: nwraldy= 0       !Niveau d'impression de l'alloc dynamique 
integer :: nwrmai = 0       !Niveau d'impression des caracteristiques  
integer :: nwrint = 0       !Niveau d'impression des tableaux
logical :: lerr

character(len=72) :: s1, s2, s3, s4, s5, s6
real              :: x1, x2, x3, x4
integer           :: n1, n2, n3, n4

namelist/nlsyst/ titre,  &      !Titre
                 dirpr,  &      !Repertoire d'ecriture des fichiers sortie
                 iout,   &      !Sortie = 0 (X11), 1 (fichier.out)
                 dt,     &      !pas de temps
                 nstep,  &      !nbre d'iterations maxi
                 lmodte, &      !Mode TE
                 lmodtm, &      !Mode TM
                 lcorrp, &      !Activateur de la correction de Poisson
                 c,      &      !vitesse de la lumiere
                 eps0,   &      !permittivite electrique
                 nschme, &      !Numero du schema, mode TE
                 nschmm, &      !Numero du schema, mode TM
                 epsche, &      !Coefficient pour le schema mode TE
                 epschm, &      !Coefficient pour le schema mode TM
                 nsolve, &  !Type de solveur
		 gsq,    &  !gamma au carre pour la correction hyperbolique
                 ldemar, lpospr, lprote, lrepri, nmaill,    &
                 lstmai, ldtfrc, pasdt, minter, mcentr, mexter, &
                 mitmax, metpoi, errpoi, nitgc, pivpoi, lctrac,	&
		 nwrpar, ldebug, landau, lmtest, lcopic
    
!Valeurs par defaut
dt    = 5.0e-13
dirpr = './'

inpfil = trim(prefix)//".inp"
open(10,file=inpfil,status='old')
write(*,*)"Lecture de la namelist $nlsyst"
read(10,nlsyst,err=100) 
goto 200
100 continue
write(*,*)"Erreur de lecture de la namelist $nlsyst"
write(*,*)"Valeurs prises par defaut"
write(*,nlsyst)
stop
200 continue
close(10)

if (iout == 0 ) then 
   iout = 6
else
   iout = 7
   open(7, file=trim(titre)//".out")
end if

if (minter>0 .and. mcentr>0 .and. mexter>0) then
   nstep = minter*mcentr*mexter
   write(iout,900)
   write(iout,"(/10x,a)")      " Les valeurs de minter, mcentr, mexter et mcentr "
   write(iout,"( 10x,a)")      " ne sont pas compatibles avec nstep "
   write(iout,"(/10x,a,i6,/)") " nstep = ", nstep
end if

if (pasdt > 0) then
   dt = dble(pasdt)
else
   pasdt = dt
end if
csq     = c * c

write(iout,1100)titre, dt, eps0, nstep, " Poisson ", " 2D "

!write(iout,937) crtupd
write(iout,915) nmaill,nsolve,metpoi,metamp  &
                ,lmodte,lmodtm,lcorrp   &
                ,nschme,nschmm,epsche,epschm
write(iout,911) ldemar,lpospr,lrepri,numprt,lprote,dirpr    &
                ,lsecur,tlmcpu,lrztrk,lrzsmp,lrzhst
write(iout,912) pasdt,minter,mcentr,mexter,mitmax,iter,temps    &
                ,ldtfrc
write(iout,913) lctrac,lcputm,lmainv,lstmai &
                ,nwraldy,nwrmai,nwrint,nwrpar
write(iout,921) errpoi,erramp,nitgc,pivpoi,pivamp

write(iout,901) 

write(iout,902) 

if(nmaill == 0) then
   write(iout,933) 
else if(nmaill == 1) then
   write(iout,905) 
else if(nmaill == 2) then
   write(iout,936) 
else if(nmaill == 3) then
   write(iout,906) 
else
   write(iout,907) nmaill
   goto 800
end if


if(nsolve == 0) then
   write(iout,908) 
else if(nsolve == 01) then
   write(iout,920) 
else if(nsolve == 02) then
   write(iout,922) 
else if(nsolve == 03) then
   write(iout,923) 
else if(nsolve == 10) then
   write(iout,910) 
else if(nsolve == 11) then
   write(iout,911) 
else if(nsolve == 21) then
   write(iout,911) 
else
   write(iout,909) nsolve
   goto 800
end if

if (nschme < 0.OR.nschme > 2) then
   write(iout,9091) nschme
   goto 800
end if

if (nschmm < 0.OR.nschmm > 2) then
   write(iout,9092) nschmm
   goto 800
end if

if (epsche < 0..OR.epsche > 1.) then
   write(iout,9093) epsche
   goto 800
end if

if (epschm < 0..OR.epschm > 1.) then
   write(iout,9094) epschm
   goto 800
end if

if (nsolve == 1.OR.nsolve == 3) then
   if(metpoi == 0) then
      write(iout,927) 
   elseif(metpoi == 01) then
      write(iout,928) 
   else
      write(iout,929) metpoi
      goto 800
   end if
end if

if (nsolve == 2.OR.nsolve == 3) then
   if(metamp == 0) then
      write(iout,930) 
   elseif(metamp == 01) then
      write(iout,931) 
   else
      write(iout,932) metamp
      goto 800
   end if
end if

if(lmodte) then
   write(iout,916)
   if(lcorrp) then
      write(iout,919)
   end if
end if

if(lmodtm) write(iout,917)

!write(iout,924) minter*mcentr*mexter
!if(nstep > mitmax) then
!   write(iout,925) mitmax
!   !goto 800
!end if

!lerr=.FALSE.
!if(pasdt.LE.0)  lerr=.TRUE.
!if(minter < 0) lerr=.TRUE.
!if(mcentr < 0) lerr=.TRUE.
!if(mexter < 0) lerr=.TRUE.
!if(mitmax < 0) lerr=.TRUE.
!if(errpoi.LE.0.OR.errpoi.GE.1) lerr=.TRUE.
!if(erramp.LE.0.OR.erramp.GE.1) lerr=.TRUE.
!if(nitgc.LE.0) lerr=.TRUE.

!if(lerr) then
!   write(iout,926)
!   goto 800
!end if
      
! --- 7.0 --- Parametres du code ---------------------------------------

!write(iout,935) ndaux,nesmx,nfrmx,ninjv,ntrkv,ninjv,nphmx  &
!                ,nplmx,nsdmx,nanmx,nspmx,ncpmx         &
!                ,nbuff,ncxmx,nemimx                &
!                ,nfcrmx,ndommx
!
! --- 9.0 --- Formats --------------------------------------------------
800 continue
 
 900  format(//10x,     &
     & '++++++++  Erreur dans  nlsyst   ++++++++'///)
 901  format(//10x,'======== CONDITIONS DU CALCUL ========'/)
 902  format(//10x,'Geometrie cartesienne X,Y'/)
 903  format(//10x,'Geometrie axisymetrique Z,R'/)
 905  format(//10x,'Fichier de maillage MODULEF'/)
 906  format(//10x,'Fichier de maillage FLUX2D '/)
 907  format(//10x,'Type de mailleur :    nmaill =',i3/     &
     &         10x,'Option non implementee'/)
 908  format(//10x,'Pas de solveur de champs'/)
 909  format(//10x,'Numero de solveur  :    nsolve =',i3/       &
     &         10x,'Option non implementee'/)
 9091 format(//10x,'Type de schema temporel en mode TE : nschme = ',I10/&
     &         10x,'Option non implementee'/)
 9092 format(//10x,'Type de schema temporel en mode TM : nschmm = ',I10/&
     &         10x,'Option non implementee'/)
 9093 format(//10x,'Coef. du schema temp. en mode TE : epsche = ',E12.3/&
     &         10x,'Choix impossible'/)
 9094 format(//10x,'Coef. du schema temp. en mode TM : epschm = ',E12.3/&
     &         10x,'Choix impossible'/)
 910  format(//10x,'Solveur de Maxwell instationnaire'/&
     &         10x,'Methode duale des volumes finis'/) 
 911  format( 10x,'Demarrage                  ldemar : ',L10/   &
     &        10x,'Post-processeur            lpospr : ',L10/   &
     &        10x,'Reprise                    lrepri : ',L10/   &
     &        10x,'Numero de la prot. lue     numprt : ',I10/   &
     &        10x,'Protection                 lprote : ',L10/   &
     &        10x,'Directory des protections  dirpr  : '/,A80/  &
     &        10x,'Protection de securite     lsecur : ',L10/   &
     &        10x,'Temps limite associe       tlmcpu : ',E10.1/ &
     &        10x,'RAZ traceurs               lrztrk : ',L10/   &
     &        10x,'RAZ echantillons           lrzsmp : ',L10/   &
     &        10x,'RAZ de m2v.hst.nomj        lrzhst : ',L10/)
 912  format(/10x,'Pas de temps               pasdt  : ',E12.4/ &
     &        10x,'Boucles internes           minter : ',I10/   &
     &        10x,'Boucles centrales          mcentr : ',I10/   &
     &        10x,'Boucles externes           mexter : ',I10/   &
     &        10x,'Nombre max iterations      mitmax : ',I10/   &
     &        10x,'Numero initial iteration   iter   : ',I10/   &
     &        10x,'Temps courant initial      temps  : ',E12.3/ &
     &        10x,'On force le calcul         ldtfrc : ',L10)
 913  format(/10x,'Activation des traces      lctrac : ',L10/   &
     &        10x,'Calcul des temps CPU       lcputm : ',L10/   &
     &        10x,'Inversion du maillage      lmainv : ',L10/   &
     &        10x,'Stockage du maillage       lstmai : ',L10/   &
     &        10x,'Niv. imp. alloc. dyn.      nwraldy: ',I10/   &
     &        10x,'Niv. imp. du maillage      nwrmai : ',I10/   &
     &        10x,'Niv. imp. methode E.F.     nwrint : ',I10/   &
     &        10x,'Niv. imp. diagno. part.    nwrpar : ',I10)
 915  format(/10x,'Fichier de maillage        nmaill : ',I10/   &
     &        10x,'Numero de solveur          nsolve : ',I10/   &
     &        10x,'Methode pour POISSON       metpoi : ',I10/   &
     &        10x,'Methode pour AMPERE        metamp : ',I10/   &
     &        10x,'Resolution systeme 1       lmodte : ',L10/   &
     &        10x,'Resolution systeme 2       lmodtm : ',L10/   &
     &        10x,'Correction de Poisson      lcorrp : ',L10/   &
     &        10x,'Schema temporel mode TE    nschme : ',I10/   &
     &        10x,'Schema temporel mode TM    nschmm : ',I10/   &
     &        10x,'Coef. schema temp. mode TE epsche : ',E12.3/ &
     &        10x,'Coef. schema temp. mode TM epschm : ',E12.3)
 916  format(/10x,'Resolution du systeme E1,E2,B3'/)
 917  format(/10x,'Resolution du systeme B1,B2,E3'/)
 919  format( 10x,'Avec la correction de POISSON'/)
 920  format(//10x,'Solveur electrostatique '/  &
     &         10x,'Methode des elements finis '/) 
 921  format( 10x,'Max erreur pour Poisson    errpoi : ',E12.3/ &
     &        10x,'Max erreur pour Ampere     erramp : ',E12.3/ &
     &        10x,'Nb max iter grad conj      nitgc  : ',I10/   &
     &        10x,'Pivot Cholesky Poisson     pivpoi : ',E12.3/ &
     &        10x,'Pivot Cholesky Ampere      pivamp : ',E12.3)
 922  format(//10x,'Solveur magnetostatique '/  &
     &         10x,'Methode des elements finis '/) 
 923  format(//10x,'Solveur electrostatique et magnetostatique '/   &
     &         10x,'Methode des elements finis '/) 
 924  format( /10x,'Nombre total d''iteration en temps ',I8)
 925  format(  10x,'Ce nombre est superieur a la limite  mitmax =',I8)
 926  format( /10x,'Erreur flagrante dans les donnees')
 927  format( /10x,'Resolution de POISSON par un gradient conjugue')
 928  format( /10x,'Resolution de POISSON par Cholesky')
 929  format( /10x,'metpoi = ',I5,'     mauvais choix')
 930  format( /10x,'Resolution de AMPERE  par un gradient conjugue')
 931  format( /10x,'Resolution de AMPERE  par Cholesky')
 932  format( /10x,'metamp = ',I5,'     mauvais choix')
 933  format(//10x,'Fichier de maillage ASCII'/)
 934  format(//10x,'***** REPRISE APRES PROTECTION *****'//)
 935  format(//10x,'PARAMETRES du CODE DEGAS2D'             &
     &       //10x,'ndaux - dimension de 2 tableaux de travail ',I8 &
     &        /10x,'nesmx - nombre max de types de particules  ',I8 &
     &        /10x,'nfrmx - nb max de frontieres referencees   ',I8 &
     &        /10x,'ninjv - nb max de volumes d injection      ',I8 &
     &        /10x,'ntrkv - nb max de vol de select de traceurs',I8 &
     &        /10x,'ninjv - nb max de front d injection        ',I8 &
     &        /10x,'nphmx - nb max d espace de phase           ',I8 &
     &        /10x,'nplmx - nb max de fenetres graphiques      ',I8 &
     &        /10x,'nsdmx - nombre maximum de sous-domaines    ',I8 &
     &        /10x,'nanmx - nombre maximum d antennes          ',I8 &
     &        /10x,'nspmx - nb max de points d echantillonnage ',I8 &
     &        /10x,'ncpmx - nb max de coupes                   ',I8 &
     &        /10x,'nbuff - taille du buffer de transfert      ',I8 &
     &        /10x,'ncxmx - nb max de champs externes par mode ',I8 &
     &        /10x,'nemimx- nb max de surf de diag d emittance ',I8 &
     &        /10x,'nfcrmx- nb max de frontieres diagnostiquees',I8 &
     &        /10x,'ndommx- nb max de domaines independants    ',I8/)
 936  format(//10x,'Fichier de maillage IDEAS'/)
 937  format(//10x,'Titre du calcul :'/10x,a80)

1100 format ( /5x,' name      :',  a,     &
     &       //5x,' time step = ', g15.3, &
     &        /5x,' eps0      = ', g15.3, &
     &        /5x,' nstep     = ', i10,   &
     &        /5x,' solveur   = ', 2a )       

end subroutine readin

!************************************************************

!Function: init_particules
!Initialisation des donnees concernant les particules

subroutine init_particules(inpfil)

character(len=*), intent(in) :: inpfil

namelist/nlpart/nesp, nbpam, pcharg, pmasse, alprjt, &
        lpousp, lpoust, lchapa, lchapt, lprojp, lprojt

nesp=0
lpousp = .true.; lpoust = .true.
lchapa = .true.; lchapt = .true.
lprojp = .true.; lprojt = .true.

do i=1,nesmx
   alprjt(i)=0.
   nbpam( i) =5000
   pcharg(i)=-1.6021e-19
   pmasse(i)= 9.1091e-31
end do

open(10, file = inpfil)
write(*,"(//10x,'*Lecture de la namelist $nlpart*',/)")
read(10, nlpart,err=100)
goto 200
100 continue
write(*,*)"Erreur de lecture de la namelist $nlpart"
write(*,*)"Valeurs prises par defaut"
write(*,nlpart)
stop
200 close(10)

if(nesp <= 0) then
   lpousp = .false.; lpoust = .false.
   lchapa = .false.; lchapt = .false.
   lprojp = .false.; lprojt = .false.
   lcorrp = .false.
   write(iout,901)
else if (nesp > nesmx) then
   write(iout,902) nesp,nesmx
   call errout(6,"F","particules.f90"," ")
else
   write(iout,900) nesp,lpousp,lchapa,lprojp
   do i=1,nesp
      write(iout,903) i,nbpam(i) ,pcharg(i),pmasse(i),alprjt(i) &
                     ,lpoust(i),lchapt(i),lprojt(i)
   end do
end if

nbpmax=1
if(nesp > 0) then
   do i=1,nesp
      nbpmax=max(nbpmax,nbpam(i))
      rapbol(i)=pcharg(i)/pmasse(i)
   end do
end if

 900 format(/10x,'Nombre d''especes de particules         nesp ',I10    &
            /10x,'Activa. general du pousseur           lpousp ',L10    &
            /10x,'Activa. general de l''interpolation   lchapa ',L10    &
            /10x,'Activa. general de l''assignation     lprojp ',L10)
 901  format(//10x,'nesp=0 implique : === PAS DE PARTICULES ==='//)
 902  format(//10x,'nesp=',I10,'  > nesmx=',I3,' ...IMPOSSIBLE'     &
     &       /10x,'prendre un nombre d''especes de part. plus petit ou' &
     &       /10x,'demander une compilation du code avec nesmx + grand')
 903  format(/10x,'Espece de particules No ',I3             &
     &       /10x,'Nombre maximum de particules          nbpam  ',I10   &
     &       /10x,'Charge elementaire                    pcharg ',E12.5 &
     &       /10x,'Masse  elementaire                    pmasse ',E12.5 &
     &       /10x,'coef de compensation de charge        alprjt ',E12.5 &
     &       /10x,'Activation du pousseur                lpoust ',L10   &
     &       /10x,'Activation de l''interpolation        lchapt ',L10   &
     &       /10x,'Activation de l''assignation          lprojt ',L10)

end subroutine init_particules


!!$subroutine uouvre
!!$ 
!!$!     But:                                                             *
!!$!     Date, heure, carte UPDATE                                        *
!!$!     Ouverture de tous les fichiers utilises par le code              *
!!$!     Attribution des etiquettes logiques                              *
!!$!                                                                      *
!!$ 
!!$LOGICAL :: lexist
!!$CHARACTER*13 ::nomfich
!!$
!!$
!!$! --- 1.0 --- Extraire le nom du job de l'environnement ----------------
!!$ 
!!$!CALL utntrv(cnomjb)
!!$
!!$! --- 1.5 --- Fabrication de la carte UPDATE ---------------------------
!!$! --- cette sequence est eliminee car on utilise la carte UPDATE pour 
!!$! --- definir un titre de cas de calcul
!!$!CALL utdate(cdate)
!!$!CALL uttime(ctime)
!!$!CALL utcput(tutotl)
!!$! --- cette sequence est eliminee car on utilise la carte UPDATE pour 
!!$! --- definir un titre de cas de calcul
!!$!     CALL utcupd(cnomjb,cdate,ctime,crtupd)
!!$
!!$! --- 2.0 --- Nombres de reference -------------------------------------
!!$ 
!!$!      CALL utrefl(petitx,grandx,xnbdig)
!!$ 
!!$! --- 3.0 --- Ouverture des fichiers -----------------------------------
!!$ 
!!$! ... Flags d'ouverture de fichiers ....................................
!!$
!!$lfdat=.FALSE.
!!$lflst=.FALSE.
!!$lfmai=.FALSE.
!!$lfmaf=.FALSE.
!!$lfmas=.FALSE.
!!$lfmaa=.FALSE.
!!$lftrk=.FALSE.
!!$lfsmp=.FALSE.
!!$lfxcr=.FALSE.
!!$lfcha=.FALSE.
!!$lfpar=.FALSE.
!!$lfhst=.FALSE.
!!$lfinj=.FALSE.
!!$lfprz=.FALSE.
!!$
!!$! ... Unite nflst ........ Fichier de sortie (listing) .................
!!$ 
!!$jerr = 1
!!$nflst = 6
!!$nomfich = 'm2v.out.'//cnomjb(1:5)
!!$ 
!!$CALL utunit(nflst,iuerr)
!!$if(iuerr.ne.0) goto 800
!!$
!!$INQUIRE(file=nomfich,exist=lexist,err=800)    
!!$
!!$if(.NOT.lexist) then
!!$   OPEN(nflst,file=nomfich,status='new',form='formatted',    &
!!$                           access='sequential',err=800)
!!$else
!!$   OPEN(nflst,file=nomfich,status='old',form='formatted',    &
!!$                           access='sequential',err=800)
!!$   REWIND nflst
!!$endif
!!$                                                           
!!$write(nflst,900)
!!$write(nflst,908) cdate,ctime,cnomjb
!!$!write(nflst,907) crtupd
!!$write(nflst,909)
!!$
!!$write(nflst,903) nflst,nomfich
!!$lflst=.TRUE.
!!$ 
!!$! ... Unite nfdat ........ Fichier de donnees ..........................
!!$ 
!!$jerr = 2
!!$nomfich = 'm2v.inp.'//cnomjb(1:5)
!!$ 
!!$CALL utunit(nfdat,iuerr)
!!$if(iuerr.NE.0) goto 800
!!$INQUIRE(file=nomfich,exist=lexist,err=800)    
!!$ 
!!$if(lexist) then
!!$   OPEN(nfdat,file=nomfich,status='old', &
!!$              form='formatted',access='sequential',err=800)
!!$   write(nflst,903) nfdat,nomfich
!!$   REWIND nfdat
!!$else
!!$   write(nflst,905)
!!$   goto 800
!!$endif
!!$lfdat=.TRUE.
!!$ 
!!$! --- 8.0 --- Problemes rencontres ==> arret du job --------------------
!!$ 
!!$      goto 850
!!$ 800  CONTINUE
!!$      write(nflst,902) jerr
!!$      !CALL utabrt( nom_routine )
!!$      stop
!!$ 850  CONTINUE
!!$ 
!!$! --- 9.0 --- Formats --------------------------------------------------
!!$ 
!!$900 FORMAT(//10x,'-------------------------------------------------'     &
!!$   &        /10x,'                      DEGAS2D                    '     &
!!$   &        /10x,'            Code MAXWELL couple VLASOV           '     &
!!$   &        /10x,'                       2D1/2                     '     &
!!$   &        /10x,'                    Version 1.7                  '     &
!!$   &       //10x,'    Service Mathematiques et Codes Numeriques    '     &
!!$   &        /10x,'     Departement de Mathematiques Appliquees     '     &
!!$   &        /10x,'       Centre d''Etudes de Limeil-Valenton       '     &
!!$   &        /10x,'     F-94195 Villeneuve Saint Georges Cedex      ' &
!!$   &        /10x,'-------------------------------------------------'/)
!!$902 FORMAT(/10x,'Erreur a l''ouverture des fichiers dans UOUVRE'/    &
!!$   &        20x,' jerr = ',i5/)
!!$903 FORMAT(10x,'Unite logique',i3,'    fichier  ',a14)
!!$905 FORMAT(/10x,'Le fichier de donnees n''est pas present'//)
!!$907 FORMAT(/10x,'Carte UPDATE :'/1x,a80/)                         
!!$908 FORMAT(/10x,'Date  : ',a8/10x,'Heure : ',a8/10x,'Job   : ',a8/)
!!$909 FORMAT(//10x,'>>>>>>  Ouverture des fichiers de base  <<<<<<'    &
!!$   &        /10x,'(les autres fichiers sont ouverts si necessaire)'/)
!!$ 
!!$end subroutine uouvre

end module zone
