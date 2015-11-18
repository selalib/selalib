!File: Module Poisson
!Solveur de Poisson sur un maillage non structure
!
! Traduit en Fortran 90 a partir de M2V ou DEGAS2D
!
module sll_m_tri_poisson
#include "sll_working_precision.h"
#include "sll_memory.h"
#include "sll_errors.h"
use sll_m_triangular_meshes
use sll_m_choleski

implicit none

private

interface sll_create
  module procedure initialize_poisson_solver
  module procedure initialize_poisson_solver_from_file
end interface sll_create

interface sll_delete
  module procedure delete_tri_poisson
end interface sll_delete

public :: new_triangular_poisson_2d
public :: sll_create
public :: sll_compute_phi_from_rho
public :: sll_compute_e_from_rho
public :: sll_compute_e_from_phi 
public :: sll_delete 

!    Caracteristiques du maillage triangulaire:         
!
!        nbs     - nombre de noeuds du maillage         
!        nbt     - nombre de triangles du maillage     
!        ndiric  - nombre de noeuds verifiant Dirichlet              
!        ndirb3  - nombre de noeuds verifiant Dirichlet pour Ampere 
!        nbfrax  - nombre de noeuds a l'intersection axe/frontiere  
!        nmxfr   - nombre de frontieres referencees max            
!        nmxsd   - nombre de sous-domaines references max         
!        nefro   - nombre d'elements ayant au moins 1 sommet frontiere
!        nelfr   - nombre d'elements frontieres                      
!        nelin   - nombre d'elements internes                       
!        irefdir - references de frontieres Dirichlet non homogenes  
!        nnoeuf  - nombre de noeuds frontieres Dirichlet non homogenes
!
!    Dimensions de reference:
!
!        xlml   - limite inf x du domaine                          
!        xlmu   - limite sup x du domaine                         
!        ylml   - limite inf y du domaine                        
!        ylmu   - limite sup y du domaine                       
!        petitl - petite longueur de reference                 
!        grandl - grande longueur de reference                
!        imxref - nombre entier eleve de reference           
!                                                           
!    Caracteristiques du maillage de Delaunay-Voronoi:
!
!        nbcoti - nombre de cotes internes Delaunay       
!        nbtcot - nombre total de cotes Delaunay         
!        nbcfli - nombre de cotes internes ne satisfaisant pas CFL
!        ncotcu - nombre de cotes cumules par type de frontiere  
!        nnref  - nombre de references de frontieres Dirichlet non
!                 homogenes                                      
!                                                               
!    Caracteristiques des cotes frontieres:
!
!        nctfrt - nombre total de cotes frontieres            
!        nctfro - nombre de cotes frontieres par reference   
!        nctfrp - pointeur des tableaux de cotes frontieres 
!                                                          
!    Caracteristiques des frontieres internes: 
!
!        nnofnt - nombre de noeuds sur les frontieres internes     
!        ntrfnt - nombre total de triangles (ceux de droite)      
!        ntrfrn - nombre de triangles par frontiere              
!        ntrfrc - nombre cumule de triangles                    
!        nndfnt - noeuds Dirichlet sur les frontieres internes 
!        ncdfnt - cotes  Dirichlet sur les frontieres internes (VF)
!                                                                 
!    Tableau des classes d'elements:
!
!        nmxcol - nombre de classes des elements                
!        nclcol - nombre d'elements dans chaque classe         
!
sll_real64, parameter :: grandx = 1.d+20

!>@brief
!> Derived type for Poisson solver on unstructured mesh with triangles.
!> @details
!> We using here P1-conformin finite element method.
!> This program is derived from F. Assous and J. Segre code M2V.
type, public :: sll_triangular_poisson_2d

  private
  sll_real64, dimension(:), allocatable :: vnx    !< normal vector on node
  sll_real64, dimension(:), allocatable :: vny    !< normal vector on node
  logical,    dimension(:), allocatable :: naux   !< true if the node is on boundary
  sll_real64, dimension(:), allocatable :: gradx  !< x-grad matrix
  sll_real64, dimension(:), allocatable :: grady  !< y-grad matrix
  sll_real64, dimension(:), allocatable :: grgr   !< grad-grad matrix
  sll_int32,  dimension(:), allocatable :: mors1  !< array to store matrix
  sll_int32,  dimension(:), allocatable :: mors2  !< array to store matrix
  sll_int32,  dimension(:), allocatable :: iprof  !< array to store matrix
  sll_int32,  dimension(:), allocatable :: ifron
  sll_real64, dimension(:), allocatable :: amass
  sll_real64, dimension(:), allocatable :: vtantx
  sll_real64, dimension(:), allocatable :: vtanty
  sll_real64, dimension(:), allocatable :: sv1
  sll_real64, dimension(:), allocatable :: sv2

  sll_int32  :: ndiric
  sll_int32  :: ntypfr(5)
  sll_real64 :: potfr(5)
  sll_real64 :: eps0 = 1.0_f64 !8.8542e-12

  type(sll_triangular_mesh_2d), pointer :: mesh => null()

end type sll_triangular_poisson_2d

contains

function new_triangular_poisson_2d(mesh, ntypfr, potfr) result (solver)

type(sll_triangular_poisson_2d), pointer          :: solver
type(sll_triangular_mesh_2d),  intent(in), target :: mesh
sll_int32,  dimension(:),      intent(in)         :: ntypfr 
sll_real64, dimension(:),      intent(in)         :: potfr 

sll_int32 :: ierr

SLL_ALLOCATE(solver,ierr)
call initialize_poisson_solver(solver, mesh, ntypfr, potfr)

end function new_triangular_poisson_2d

!> Delete the solver derived type
subroutine delete_tri_poisson( this )
type(sll_triangular_poisson_2d), pointer :: this

#ifdef DEBUG
print*, 'delete poisson solver'
#endif
nullify(this)

end subroutine delete_tri_poisson

!> Compute electric field from charge density
!> @param[in]  this solver derived type
!> @param[in]  rho  charge density array on nodes
!> @param[out] phi  electric potential on nodes
!> @param[out] ex   electric field x component on nodes
!> @param[out] ey   electric field y component on nodes
subroutine sll_compute_e_from_rho( this, rho, phi, ex, ey)
type(sll_triangular_poisson_2d) :: this
sll_real64, intent(in)          :: rho(:)
sll_real64, intent(out)         :: phi(:)
sll_real64, intent(out)         :: ex(:)
sll_real64, intent(out)         :: ey(:)

call poissn(this, rho, phi, ex, ey)
call poifrc(this, ex, ey)

end subroutine sll_compute_e_from_rho

subroutine sll_compute_phi_from_rho( this, rho, phi)
type(sll_triangular_poisson_2d) :: this
sll_real64, intent(in)          :: rho(:)
sll_real64, intent(out)         :: phi(:)

call poissn(this, rho, phi)

end subroutine sll_compute_phi_from_rho

subroutine sll_compute_e_from_phi( this, phi, ex, ey)
type(sll_triangular_poisson_2d) :: this
sll_real64, intent(in)          :: phi(:)
sll_real64, intent(out)         :: ex(:)
sll_real64, intent(out)         :: ey(:)

call poliss(this, phi, ex, ey)
call poifrc(this, ex, ey)

end subroutine sll_compute_e_from_phi

subroutine read_data_solver(ntypfr, potfr)

sll_int32,  parameter       :: nfrmx = 5
sll_int32,  intent(out)     :: ntypfr(nfrmx)
sll_real64, intent(out)     :: potfr(nfrmx)

sll_int32                   :: ityp
sll_int32                   :: ifr
sll_int32                   :: i
character(len=72)           :: argv
character(len=132)          :: inpfil
logical                     :: lask
character(len=*), parameter :: this_sub_name = 'read_data_solver'

NAMELIST/nlcham/ntypfr,potfr

call get_command_argument( 1, argv); write(*,'(1x, a)') argv

!------------------------------------------------------------!
!     Reads in default parameters from input file (.inp)        !
!     If LASK=t, ask user if file is to be read.                !
!------------------------------------------------------------!

lask = .true.

inpfil = trim(argv)//'.inp' 

write(*,"(/10x,'Fichier d''entree : ',a)") inpfil

open(10,file=inpfil,status='OLD',err=80)
write(*,1050,advance='no') trim(inpfil)
lask = .false.

80 continue

if (lask) then
   write(*,1900) inpfil
   write(*,1800,advance='no')
   read(*,"(a)") inpfil
   write(*,1700) trim(inpfil)
end if

! ----------- Valeurs par defaut et initialisations -------------------t
write(6,900)

ntypfr = 1  
potfr  = 0.0_f64

!--- 2.0 --- Lecture des donnees --------------------------------------

open(10, file = inpfil)
write(*,*)"Lecture de la namelist nlcham"
read(10,nlcham,err=100)
goto 200
100 continue
write(*,*)"Erreur de lecture de la namelist nlcham"
write(*,*)"Valeurs par defaut"
write(*,nlcham)
stop
200 close(10)

!--- 2.5 --- Ecriture des donnees -------------------------------------

write(6,921)

do i=1,nfrmx
   write(6,922) i,ntypfr(i),potfr(i)
end do

!--- 3.0 --- Controle des donnees -------------------------------------

write(6,932)
do ifr=1,nfrmx
  ityp=ntypfr(ifr)
  if(ityp == 1) then
    write(6,902) ifr
  else if(ityp == 3) then
    write(6,904) ifr
  else
    write(6,910) ifr,ityp
    SLL_ERROR( this_sub_name, "Unspecified error.")
  end if
end do

900  format(//10x,'Conditions aux limites sur les frontieres')
902  format(/10x,'Reference :',i3,5x,'Dirichlet ')
904  format(/10x,'Reference :',i3,5x,'Neumann')
910  format(/10x,'Option non disponible'/  &
     &       10x,'Reference :',i3,5x,'Type :',i3/)
921  format(//5x,'Frontiere',5x,'ntypfr',7x,'potfr')
922  format(5x,I3,4x,I10,2E12.3,2I10,E12.3)
932  format(//10x,'CONDITIONS AUX LIMITES'/)
1050 format(/' Read settings from file  ', A, ' ?  Y')
1700 format(/' New parameters write to file  ', A, /)
1800 format(/' Settings may have been changed - New title :')
1900 format(/' Input file  ', A,'  not found')

end subroutine read_data_solver

subroutine initialize_poisson_solver_from_file(this, mesh)

type(sll_triangular_poisson_2d),  intent(out)     :: this
type(sll_triangular_mesh_2d),  intent(in), target :: mesh
sll_int32                                         :: ntypfr(5)
sll_real64                                        :: potfr(5)

call read_data_solver(ntypfr, potfr)  

call initialize_poisson_solver(this, mesh, ntypfr, potfr)

end subroutine initialize_poisson_solver_from_file

!Subroutine: init_solveur_poisson   
! Reservation de tableaux et calcul des 
! matrices necessaires au solveur de poisson   
!
! amass - matrice de masse diagonale     
! mors1 - matrice morse2       
! mors2 - elements des matrices morses
! prof  - matrice profil    
! grgr  - matrice grad-grad
! gradx - matrice gradx   
! grady - matrice grady  
! fron  - noeuds Dirichlet
! xaux  - tableaux auxiliaires reels          
! iaux  - tableaux auxiliaires entiers       

! Allocation des tableaux permettant de stocker des matrices sous forme morse
! Tableau donnant le numero du dernier terme de chaque ligne (mors1)
subroutine initialize_poisson_solver( this, mesh, ntypfr, potfr )

type(sll_triangular_poisson_2d), intent(out)         :: this
type(sll_triangular_mesh_2d)   , intent(in ), target :: mesh
sll_int32                      , intent(in )         :: ntypfr(:)
sll_real64                     , intent(in )         :: potfr(:)

sll_real64, allocatable     :: tmp1(:)
character(len=*), parameter :: this_sub_name = 'read_data_solver'
character(len=128)          :: err_msg

sll_int32 :: nref, nn, ndir
sll_int32 :: i
sll_int32 :: ierr

if ( mesh%analyzed) then
  this%mesh => mesh
else
  err_msg = "Call analyze_triangular_mesh before initialize_poisson_solver."
  SLL_ERROR( this_sub_name, err_msg )
endif

this%ntypfr = ntypfr
this%potfr  = potfr

allocate(this%sv1(mesh%nbtcot));    this%sv1    = 0.0_f64
allocate(this%sv2(mesh%nbtcot));    this%sv2    = 0.0_f64
allocate(this%vtantx(mesh%nbtcot)); this%vtantx = 0.0_f64
allocate(this%vtanty(mesh%nbtcot)); this%vtanty = 0.0_f64

allocate(this%mors1(mesh%num_nodes+1)); this%mors1 = 0

! Tableau contenant le numero des termes de chaque ligne (mors2)
! on choisit une taille a priori superieure a la taille reelle
! de ce tableau.

allocate(this%mors2(12*mesh%num_nodes)); this%mors2 = 0
 
! Calcul de mors1,mors2.

call morse(mesh%npoel1,        &
           mesh%npoel2,        &
           mesh%nodes,         &
           mesh%num_triangles, &
           mesh%num_nodes,     &
           this%mors1,         &
           this%mors2          )
 
! Ajustement de la taille de mors2.
! pas sur que ca fonctionne a tous les coups
! deallocate(mors2); allocate(mors2(mors1(mesh%num_nodes+1)))

! Adressage des tableaux permettant de stocker des matrices sous forme profil

allocate(this%iprof(mesh%num_nodes+1)); this%iprof = 0

call profil(mesh%nodes,     &
            mesh%num_nodes, &
            mesh%npoel1,    &
            mesh%npoel2,    &
            this%iprof      )

!======================================================================
!--- 2.0 --- POISSON par une methode d'elements finis -----------------
!======================================================================
 
!matrice de masse diagonalisee
allocate(this%amass(mesh%num_nodes)); this%amass = 0.0_f64 
 
!matrice "grad-grad" stockee sous forme profil.
allocate(this%grgr(this%iprof(mesh%num_nodes+1))); this%grgr = 0.0_f64
 
!gradx et grady 
allocate(this%gradx(this%mors1(mesh%num_nodes+1))); this%gradx = 0.0_f64
allocate(this%grady(this%mors1(mesh%num_nodes+1))); this%grady = 0.0_f64

!--- Tableau relatif aux frontieres Dirichlet -------------------------

ndir=0
do nn=1,mesh%num_nodes
   nref= mesh%refs(nn)
   if (nref > 0) then 
      if (this%ntypfr(nref)==1)then
         ndir=ndir+1
      end if 
   end if 
end do

allocate(this%ifron(ndir)); this%ifron = 0

ndir=0
do nn=1,mesh%num_nodes
   nref= mesh%refs(nn)
   if (nref > 0) then 
      if (this%ntypfr(nref)==1)then
         ndir=ndir+1
         this%ifron(ndir)=nn
      end if 
   end if 
end do

this%ndiric=ndir

!--- Calcul des matrices ----------------------------------------------

call poismc(this)

!Calcul de la matrice B tel que B*Bt = A dans le cas Cholesky

allocate(tmp1(this%iprof(mesh%num_nodes+1))); tmp1 = 0.0_f64

!write(*,"(//5x,a)")" *** Appel Choleski pour Poisson ***  "
call choles(this%iprof,this%grgr,tmp1)

do i=1,this%iprof(mesh%num_nodes+1)
   this%grgr(i)=tmp1(i)
end do

deallocate(tmp1)

!Tableaux pour les conditions limites
!Initialisation des normales aux noeuds et du tableau indiquant 
!les noeuds appartenant a la frontiere consideree

SLL_CLEAR_ALLOCATE(this%vnx(1:this%mesh%num_nodes),ierr)
SLL_CLEAR_ALLOCATE(this%vny(1:this%mesh%num_nodes),ierr)
SLL_ALLOCATE(this%naux(1:this%mesh%num_nodes),ierr)
this%naux = .false.

! --- 8.5 --- Ecriture des tableaux ------------------------------------

!if (ldebug) then
!   write(6,900) 
!   do i=1,mesh%num_nodes
!      write(6,901) i,(this%mors2(j), j=this%mors1(i)+1,this%mors2(i+1))
!   end do
!   write(6,902) 
!   write(6,903) (this%ifron(i),i=1,this%ndiric)
!end if
!
!! ======================================================================
!! --- 9.0 --- Formats --------------------------------------------------
! 
!900 format(//10x,'Tableau MORS2 pointe par le tableau MORS1'/  &
!              3x,'No de noeud           Noeuds associes'/)
!901 format(2x,I8,3x,12I8)
!902 format(//10x,'Noeuds frontiere du type DIRICLET pour POISSON'/)
!903 format(32000(2x,7I9/)/)

end subroutine initialize_poisson_solver

!======================================================================


!Function: morse
! Calcul du tableau contenant l'adresse du dernier terme de   
! chaque ligne des matrices "morse" et le tableau contenant les
! numeros des termes de ces matrices. Le terme diagonal est    
! range en derniere place de chaque ligne.                      
!          
! npoel1 - emplacement dans npoel2 du dernier element relatif a chaque
!          noeud avec npoel1(1)=0 et npoel1(i+1) relatif au noeud i   
! npoel2 - tableau des numeros des elements ayant un sommet en commun  
! ntri   - numeros des sommets des triangles                 
! nbt    - nombre de triangles du maillage                    
! nbs    - nombre de noeuds                                    
!                                                               
! mors1  - tableau des adresses des derniers termes de chaque ligne 
!          de la matrice avec la convention:                         
!          mors1(1)=0 et mors1(i+1) adresse de aii                    
! mors2  - tableau des numeros des termes des matrices "morse"         
subroutine morse(npoel1, npoel2, ntri, nbt, nbs, mors1, mors2)

sll_int32,                   intent(in)  :: nbt
sll_int32,                   intent(in)  :: nbs
sll_int32, dimension(:),     intent(in)  :: npoel1
sll_int32, dimension(:),     intent(in)  :: npoel2
sll_int32, dimension(3,nbt), intent(in)  :: ntri
sll_int32, dimension(:),     intent(out) :: mors1
sll_int32, dimension(:),     intent(out) :: mors2
sll_int32, dimension(20)                 :: ilign

sll_int32 :: l, itest1, itest2, js1, js2, is1, is2, is3, numel
sll_int32 :: iel, nlign, nel, is, im, k

im = 0
k  = 0

mors1(1)=0
 
do is=1,nbs  !Boucle sur les noeuds

  !Nombre d'elements ayant is comme sommet

  nel=npoel1(is+1)-npoel1(is)
  nlign=0

  !Boucle sur ces elements

  do iel=1,nel

     k=k+1

     !numero de l'element

     numel=npoel2(k)
     is1=ntri(1,numel); is2=ntri(2,numel); is3=ntri(3,numel)
     if (is1==is) then 
        js1=is2; js2=is3
     end if 
     if (is2==is) then 
        js1=is3; js2=is1
     end if 
     if (is3==is) then 
        js1=is1; js2=is2
     end if 

     !On regarde si les 2 noeuds autres que is de l'element courant ont 
     !deja ete pris en compte dans l'ensemble des noeuds interagissant 
     !avec is.

     itest1=0; itest2=0
     if (nlign.NE.0) then 
        do l=1,nlign
           if (js1==ilign(l)) then 
              itest1=1
           end if 
           if (js2==ilign(l)) then 
              itest2=1
           end if 
        end do
     end if 

     if (itest1==0) then 
        nlign=nlign+1
        ilign(nlign)=js1
     end if 
     if (itest2==0) then 
        nlign=nlign+1
        ilign(nlign)=js2
     end if 

  end do

  !Definition de l'adresse du dernier terme de la ligne

  mors1(is+1)=mors1(is)+nlign+1

  !Remplissage du tableau mors2 avec les numeros des termes
  !de la ligne is.

  if (nlign.NE.0) then 
     do l=1,nlign
        im=im+1
        mors2(im)=ilign(l)
     end do
  end if 
  im=im+1
  mors2(im)=is
end do

end subroutine morse


! ======================================================================


!Function: poismc
!                                                 
! But: 
!   Calcul de toutes les matrices du systeme   
!       pour POISSON cartesien                      
!                                                    
! Parametres d'entree:                                
!
! m      - super-tableau                                        
! coor   - coordonnees des noeuds                                
! refs   - indice permettant de savoir si un noeud appartient     
!          a une frontiere referencee                              
! ifron  - tableau des noeuds Dirichlet                             
! mors1  - tableau du nombre de termes par ligne de la matrice morse 
!          symetrique                                                 
! mors2  - tableau des numeros des termes de chaque ligne de la matrice
!          morse symetrique                              
! ntri   - numeros des sommets des triangles              
! aire   - aire des triangles                              
! iprof  - profil de la matrice grad-grad                   
!                                                            
! noefnt - noeuds internes Dirichlet                          
!                                                              
! Parametres resultats:
! 
! agrgr  - matrice de l'operateur "grad-grad" sous forme "morse"
! aliss  - matrice de lissage                                   
! amass  - matrice de masse sous forme diagonalisee              
! amclt  - matrice de masse diagonalisee relative a la condition  
!          aux limites absorbante                             
! aroro  - matrice de l'operateur "rot-rot" sous forme "morse" 
! d1dx,  - derives des fonctions de base dans chaque triangle   
!..d3dy                                                          
! gradx  - matrice de l'operateur gradx                           
! grady  - matrice de l'operateur grady                            
! rotx   - matrice de l'operateur rotx                              
! roty   - matrice de l'operateur roty                               
! aire   - surface de chaque element                                  
!Auteur:
!   612-POISMC     
!
! Puertolas - Version 1.0  Octobre  1992  
subroutine poismc(this)

type(sll_triangular_poisson_2d),  intent(inout) :: this
sll_real64 :: amloc(3),aggloc(9),grxloc(9),gryloc(9)
sll_real64 :: dntx1, dntx2, dntx3, dnty1, dnty2, dnty3 
sll_real64 :: x1t, x2t, x3t, y1t, y2t, y3t, coef
sll_int32 :: is1t, is2t, is3t, iel
sll_int32 :: is, j
 
!Boucle sur les elements.

do iel=1,this%mesh%num_triangles

  !Calcul des coefficients dependant de la geometrie du triangle.

  is1t = this%mesh%nodes(1,iel)
  is2t = this%mesh%nodes(2,iel)
  is3t = this%mesh%nodes(3,iel)

  x1t  = this%mesh%coord(1,is1t)
  x2t  = this%mesh%coord(1,is2t)
  x3t  = this%mesh%coord(1,is3t)

  y1t  = this%mesh%coord(2,is1t)
  y2t  = this%mesh%coord(2,is2t)
  y3t  = this%mesh%coord(2,is3t)

  dntx1 = y2t-y3t
  dntx2 = y3t-y1t
  dntx3 = y1t-y2t

  dnty1 = x3t-x2t
  dnty2 = x1t-x3t
  dnty3 = x2t-x1t

  !Contribution a la matrice de masse

  amloc(1) = this%mesh%aire(iel)/3.
  amloc(2) = this%mesh%aire(iel)/3.
  amloc(3) = this%mesh%aire(iel)/3.

  !Assemblage

  call asbld(amloc,is1t,is2t,is3t,this%amass)

  !Contribution a la matrice grad-grad

  coef=1./(4.*this%mesh%aire(iel))

  aggloc(1)=(dntx1**2   +dnty1**2   )*coef
  aggloc(2)=(dntx1*dntx2+dnty1*dnty2)*coef
  aggloc(3)=(dntx1*dntx3+dnty1*dnty3)*coef
  aggloc(4)=(dntx2**2   +dnty2**2   )*coef
  aggloc(5)=(dntx2*dntx3+dnty2*dnty3)*coef
  aggloc(6)=(dntx3**2   +dnty3**2   )*coef

  !Contribution aux matrices gradx et grady:
  !Calcul de matrices locales.

  grxloc(1)=-dntx1/6.; gryloc(1)=-dnty1/6.
  grxloc(2)=-dntx2/6.; gryloc(2)=-dnty2/6.
  grxloc(3)=-dntx3/6.; gryloc(3)=-dnty3/6.
  grxloc(4)=-dntx1/6.; gryloc(4)=-dnty1/6.
  grxloc(5)=-dntx2/6.; gryloc(5)=-dnty2/6.
  grxloc(6)=-dntx3/6.; gryloc(6)=-dnty3/6.
  grxloc(7)=-dntx1/6.; gryloc(7)=-dnty1/6.
  grxloc(8)=-dntx2/6.; gryloc(8)=-dnty2/6.
  grxloc(9)=-dntx3/6.; gryloc(9)=-dnty3/6.

  !Assemblage 

  call asblp(this%iprof, aggloc,is1t,is2t,is3t,this%grgr)

  call asblm2(grxloc,      &
              gryloc,      &
              this%mors1,  &
              this%mors2,  &
              is1t,        &
              is2t,        &
              is3t,        &
              this%gradx,  &
              this%grady)

end do

! ======================================================================
! ... Prise en compte des conditions aux limites Dirichlet

do j=1,this%ndiric
   is=this%ifron(j)             
   this%grgr(this%iprof(is+1))=grandx
end do

! ================================================================
! ... Ecriture des matrices mass, grgr, gradx et grady

!if (ldebug) then
!  write(6,900) 
!  do is=1,this%mesh%num_nodes
!     write(6,901) is,this%amass(is)
!  end do
!
!  write(6,907) 
!  write(6,902) 
!  do is=1,this%mesh%num_nodes
!     nis=this%iprof(is+1)-this%iprof(is)
!     write(6,903) is,(this%grgr(this%iprof(is)+il),il=1,nis)
!  end do
!
!  write(6,904) 
!  do is=1,this%mesh%num_nodes
!     nis=this%mors1(is+1)-this%mors1(is)
!     write(6,903) is,(this%gradx(this%mors1(is)+il),il=1,nis)
!  end do
!
!  write(6,905) 
!  do is=1,this%mesh%num_nodes
!     nis=this%mors1(is+1)-this%mors1(is)
!     write(6,903) is,(this%grady(this%mors1(is)+il),il=1,nis)
!  end do
!end if

! ================================================================

! 900 format(//10x,'Matrice de masse'/               &
!              10x,'No de noeud     Terme diagonal'/)
! 901 format(  10x,I10,E15.4)
! 902 format(//10x,'Matrice grad-grad'/              &
!               10x,'No de noeud     Termes de la ligne'/)
! 903 format(  I10,5E12.3/10(10x,5E12.3/)/)
! 904 format(//10x,'Matrice gradx'/              &
!              10x,'No de noeud     Termes de la ligne'/)
! 905 format(//10x,'Matrice grady'/              &
!               10x,'No de noeud     Termes de la ligne'/)
! 907 format(//10x,'Resolution du systeme par Choleski')

end subroutine poismc

!Function: poissn
!calculer potentiel et champ electrique de l'equation de
!Poisson en cartesien
! Variables en argument :                                
!
! ex     - e1 projete sur les noeuds        
! ey     - e2 projete sur les noeuds       
! rho    - densite de charge aux noeuds du maillage 
! phi    - potentiel aux noeuds du maillage        
! ifron  - tableau des noeuds frontaliers verifiant Dirichlet
! noefnt - noeuds sur les frontieres internes Dirichlet     
! irffnt - numero de reference de ces noeuds               
! grgr   - matrice de l'operateur "grad-grad"             
! amass  - matrice de masse sous forme diagonalisee      
! gra1   - matrice de l'operateur "grad1"               
! gra2   - matrice de l'operateur "grad2"              
! mors1  - tableau descriptif des matrices "morse"  
! mors2  - tableau descriptif des matrices "morse"
! iprof  - matrice profil                           
!
! Tableaux auxilliaires :                            
!
! grgrdd - produit de la matrice grgr par le vecteur 
!          "direction de descente" p (gradient conjugue) 
! dird   - direction de descente dans la methode de gradient conjugue 
! res    - residu dans la methode de gradient conjugue               
! sdmb   - second membre du systeme a resoudre pour traiter         
!          le potentiel                                            
! sdmb12 - second membre du systeme a resoudre pour traiter       
!          le champ electrique                                   
! precd  - tableau local utilise dans la resolution de la correction 
!          par un gradient conjugue preconditionne                  
! nbs    - nombre de noeuds du maillage                            
! nmxfr  - nombre de frontieres referencees                       
! errpoi - precision (max du carre de l'erreur relative)         
! nitgc  - nombre max d'iterations du gradient conjugue         
!                                                              
subroutine poissn(this,rho,phi,ex,ey)

type (sll_triangular_poisson_2d)                :: this
sll_real64, dimension(:), intent(in)            :: rho
sll_real64, dimension(:), intent(out)           :: phi
sll_real64, dimension(:), intent(out), optional :: ex
sll_real64, dimension(:), intent(out), optional :: ey

sll_real64, dimension(:), allocatable           :: sdmb
sll_real64, dimension(:), allocatable           :: sdmb12

sll_int32 :: i, is, nref, ierr

! ----------- CALCUL DU POTENTIEL  -------------------------------------
!
!... Calcul du second membre complet ...

SLL_ALLOCATE(sdmb(this%mesh%num_nodes),ierr)

!!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(is)

!!$OMP DO SCHEDULE(RUNTIME)
do is=1,this%mesh%num_nodes
   sdmb(is)=this%amass(is)*rho(is)/this%eps0
end do
!!$OMP END DO 
!!$OMP END PARALLEL

!... Condition aux limites Dirichlet homogene 

do is=1,this%ndiric
   nref=this%mesh%refs(this%ifron(is))
   sdmb(this%ifron(is))=this%potfr(nref)*grandx
end do

call desrem(this%iprof, this%grgr,sdmb,this%mesh%num_nodes,phi)

!*** CALCUL DES CHAMPS Ex et Ey:

if (present(ex) .and. present(ey)) then

  SLL_ALLOCATE(sdmb12(this%mesh%num_nodes),ierr)
  !*** Second membre pour la composante x:

  call m1p(this%gradx,this%mors1,this%mors2,phi,this%mesh%num_nodes,sdmb12)

  do i=1,this%mesh%num_nodes
     ex(i)=sdmb12(i)/this%amass(i)
  end do

  !*** Second membre pour la composante y:

  call m1p(this%grady,this%mors1,this%mors2,phi,this%mesh%num_nodes,sdmb12)

  do i=1,this%mesh%num_nodes
     ey(i)=sdmb12(i)/this%amass(i)
  end do

end if

end subroutine poissn


!Function: poifrc
!Correction des champs sur les frontieres, en particulier
!
! version qui ne tient pas compte de l'ordre des cotes de la frontiere
!                                                                    
!  - E.tau = 0 sur les frontieres Dirichlets
!  - E.nu =  0 sur les frontieres Neumann  
!                                                                
!Variables en argument:                                        
!
!  ksofro  - numeros des 2 sommets extremite du cote 
!  krefro  - numero de reference du cote            
!  vnofro  - composantes du vecteur normal (vers l'interieur)
!  vnx,vny - tableaux donnant les composantes d'une normale aux noeuds
!                                                                   
!Tableaux auxilliaires:                                           
!
! naux   - tableau permettant de reperer les noeuds d'une frontiere
!
subroutine poifrc(this, ex, ey)

type(sll_triangular_poisson_2d)   :: this
sll_real64 :: ex(:), ey(:)
sll_real64 :: pscal, xnor
sll_int32 :: is1, is2, ict
sll_int32 :: i
      
!!$ ======================================================================
!!$ ... On force E.tau = 0 sur toutes les frontieres Dirichlet 


!!$ ... Boucle sur les cotes frontieres pour construire les normales aux
!!$ ... noeuds "Dirichlet"
this%vnx=0.0_f64; this%vny=0.0_f64; this%naux= .false.

do ict=1,this%mesh%nctfrt

  if (this%ntypfr(this%mesh%krefro(ict))==1) then 

    is1=this%mesh%ksofro(1,ict)
    is2=this%mesh%ksofro(2,ict)

    this%vnx(is1)=this%vnx(is1)+this%mesh%vnofro(1,ict)
    this%vny(is1)=this%vny(is1)+this%mesh%vnofro(2,ict)
    this%vnx(is2)=this%vnx(is2)+this%mesh%vnofro(1,ict)
    this%vny(is2)=this%vny(is2)+this%mesh%vnofro(2,ict)

    this%naux(is1)=.true.
    this%naux(is2)=.true.

  end if 

end do 

!... on impose la condition E.tau=0

do i=1,this%mesh%num_nodes

  if (this%naux(i)) then 

     xnor=SQRT(this%vnx(i)**2+this%vny(i)**2)
     if (xnor>this%mesh%petitl) then 
        this%vnx(i)=this%vnx(i)/xnor
        this%vny(i)=this%vny(i)/xnor

        pscal=this%vnx(i)*ex(i)+this%vny(i)*ey(i)
        ex(i)=this%vnx(i)*pscal
        ey(i)=this%vny(i)*pscal
     end if 

  end if 

end do 

!======================================================================
! ... On force E.nu = 0 sur toutes les frontieres Neumann
!
! ... Initialisation des normales aux noeuds et du tableau indiquant 
! ... les noeuds appartenant a la frontiere consideree
!
this%vnx=0.0_f64; this%vny=0.0_f64; this%naux= .false.

! ... Boucle sur les cotes frontieres pour construire les normales aux
! ... noeuds "Neumann"

do ict=1,this%mesh%nctfrt

  if (this%ntypfr(this%mesh%krefro(ict))==3) then 

    is1=this%mesh%ksofro(1,ict)
    is2=this%mesh%ksofro(2,ict)

    this%vnx(is1)=this%vnx(is1)+this%mesh%vnofro(1,ict)
    this%vny(is1)=this%vny(is1)+this%mesh%vnofro(2,ict)
    this%vnx(is2)=this%vnx(is2)+this%mesh%vnofro(1,ict)
    this%vny(is2)=this%vny(is2)+this%mesh%vnofro(2,ict)

    this%naux(is1)=.true.
    this%naux(is2)=.true.

  end if 

end do 

!... on impose la condition E.nu=0

do i=1,this%mesh%num_nodes

  if (this%naux(i)) then 

    xnor=SQRT(this%vnx(i)**2+this%vny(i)**2)
    if (xnor>this%mesh%petitl) then 
      this%vnx(i)=this%vnx(i)/xnor
      this%vny(i)=this%vny(i)/xnor
      pscal=this%vnx(i)*ex(i)+this%vny(i)*ey(i)
      ex(i)=ex(i)-this%vnx(i)*pscal
      ey(i)=ey(i)-this%vny(i)*pscal
    end if 

  end if 

end do 

end subroutine poifrc


!Function: profil
!determination des numeros des elements situes sur la diagonale
!de la matrice profil associee a la methode d'elements finis.
!
!Parametres d'entree:
! npoel1 - emplacement dans npoel2 du dernier element relatif a chaque
!          noeud avec npoel1(1)=0 et npoel1(i+1) relatif au noeud i 
! npoel2 - tableau des numeros des elements ayant un sommet en commun
! ntri   - numeros des sommets des triangles                          
! nbt    - nombre de triangles         
! noe    - nombre de noeuds             
!                                                                    
!Parametre resultat:                                              
! iprof  - tableau des numeros des termes situes sur la diagonale avec
!        - iprof(i+1) numero du terme de la diagonale de la ligne i et
!        - iprof(1)=0 
!
!Auteur:
! J. Segre - Juillet 89
subroutine profil(ntri,nbs, npoel1, npoel2, iprof)

sll_int32, dimension(:) :: iprof
sll_int32, dimension(:), intent(in) :: npoel1, npoel2
sll_int32, dimension(:,:), intent(in) :: ntri
sll_int32 :: in, k, is1, is2, is3, numel, ind, iel, nel
sll_int32, intent(in) :: nbs      !Nombre de sommets

!************************************************************************
!*
!* determination de la longueur de chaque ligne de la matrice profil :
!* pour chaque noeud i on cherche les noeuds j qui lui sont lies;
!* si le noeud j n'a jamais ete rencontre, la longueur de la ligne j est
!* j-i+1
!
!* boucle sur les noeuds

k=0
do in=1,nbs
      
  !* nombre d'elements ayant in comme sommet
  nel=npoel1(in+1)-npoel1(in)

  !* boucle sur ces elements

  do iel=1,nel

    k=k+1
    !* numero de l'element
    numel=npoel2(k)
    is1=ntri(1,numel)
    ind=is1+1
    if (iprof(ind).eq.0)then
       iprof(ind)=is1-in+1
    end if 
    is2=ntri(2,numel)
    ind=is2+1
    if (iprof(ind).eq.0)then
       iprof(ind)=is2-in+1
    end if 
    is3=ntri(3,numel)
    ind=is3+1
    if (iprof(ind).eq.0)then
       iprof(ind)=is3-in+1
    end if 

  end do

end do

!* determination de la position des termes diagonaux de la matrice
!* profil (il suffit de sommer les nombres d'elements des lignes
!* precedentes et de la ligne courante).

do ind=3,nbs+1
   iprof(ind)=iprof(ind-1)+iprof(ind)
end do

end subroutine profil

!Function: asbld
!  Assembler une matrice elementaire    
!  dans une matrice globale dans          
!  le cas ou elle est diagonale            
!                                                                 
!Parametres d'entree: 
!                                 
!  aele         -    Matrice elementaire diagonale           
!                    (3 termes pour un element triangulaire)
!  i1,i2,i3     -    numeros des sommets de l'element      
!                                                                 
!Parametre resultat:                                     
!                                                              
!  xmass    -   matrice globale diagonalisee          
!
!Auteur:
!      J. Segre - Version 1.0  Juillet 1989
subroutine asbld(aele,i1,i2,i3,xmass)
 
sll_real64, dimension(:) :: aele, xmass
sll_int32 :: i1, i2, i3

xmass(i1)=xmass(i1)+aele(1)
xmass(i2)=xmass(i2)+aele(2)
xmass(i3)=xmass(i3)+aele(3)

end subroutine

!Function: asblm2
!   
!          assembler 3 matrices elementaires
!          dans 3 matrices globales stockees
!          sous forme "morse" non symetriques
!                                                          
!      Entrees:
!          aele1,aele2    -  Matrices elementaires          
!          mors1          -  Tableau du nombre de termes par 
!                            ligne de la matrice morse        
!          mors2          -  Tableau des numeros des termes    
!                            de chaque ligne de la matrice morse
!          i1,i2,i3       -  Numeros des sommets de l'element
!                                                         
!      Sorties:
!          a1,a2          -  Matrices globales stockees     
!                            sous forme "morse"            
!
!Auteur:
!J. Segre - Version 1.0  Aout 1989  
subroutine asblm2(aele1,aele2,mors1,mors2,i1,i2,i3,a1,a2)

sll_int32, intent(in) :: i1, i2, i3
sll_int32, dimension(:), intent(in) :: mors1, mors2
sll_real64,    dimension(:), intent(in)  :: aele1, aele2
sll_real64,    dimension(:), intent(out) :: a1, a2
sll_int32 :: j, j1, j2, ind1, ind2, ind3

! --- 1.1 --- Rangement des termes diagonaux ---------------------------

ind1=mors1(i1+1)
a1(ind1)=a1(ind1)+aele1(1)
a2(ind1)=a2(ind1)+aele2(1)

ind2=mors1(i2+1)
a1(ind2)=a1(ind2)+aele1(5)
a2(ind2)=a2(ind2)+aele2(5)

ind3=mors1(i3+1)
a1(ind3)=a1(ind3)+aele1(9)
a2(ind3)=a2(ind3)+aele2(9)
      
! --- 1.2 --- Rangement des autres termes ------------------------------

j2=ind1-1
j1=mors1(i1)+1
 
if (j2>=j1) then 
   do  j=j1,j2
      if (i2==mors2(j)) then 
         a1(j)=a1(j)+aele1(2)
         a2(j)=a2(j)+aele2(2)
      end if 
      if (i3==mors2(j)) then 
         a1(j)=a1(j)+aele1(3)
         a2(j)=a2(j)+aele2(3)
      end if 
   end do
end if 
 
j2=ind2-1
j1=mors1(i2)+1

if (j2>=j1) then 
   do  j=j1,j2
      if (i1==mors2(j)) then 
         a1(j)=a1(j)+aele1(4)
         a2(j)=a2(j)+aele2(4)
      end if 
      if (i3==mors2(j)) then 
         a1(j)=a1(j)+aele1(6)
         a2(j)=a2(j)+aele2(6)
      end if 
   end do
end if 

j2=ind3-1
j1=mors1(i3)+1

if (j2>=j1) then 
   do  j=j1,j2
      if (i1==mors2(j)) then 
        a1(j)=a1(j)+aele1(7)
        a2(j)=a2(j)+aele2(7)
      end if 
      if (i2==mors2(j)) then 
         a1(j)=a1(j)+aele1(8)
         a2(j)=a2(j)+aele2(8)
      end if 
   end do
end if 

end subroutine

!----------------------------------------------------------------------


!Function: asblp
!                                                      
!      assembler une matrice elementaire             
!      dans une matrice globale symetrique          
!      stockee sous forme "profil"                 
!                                                 
!Entrees:
!                                                                 
!          aele         -  Matrice elementaire (6 ou 10 termes   
!                          selon que l'element est un triangle  
!                          ou un quadrilatere)                 
!          iprof        -  Description du profil de la matrice 
!                          c'est-a-dire liste des numeros d'ordre 
!                          des termes diagonaux dans le tableau de
!                          stockage de la matrice                
!          i1,i2,i3     -  Numeros des sommets de l'element     
!                                                              
!Sorties:
!                                                             
!          xmass        -  matrice globale stockee sous forme "profil"
!                                                                    
!Auteur:
!       J. Segre - Version 1.0  Aout 1989
subroutine asblp(iprof, aele,i1,i2,i3,xmass)

sll_int32, dimension(:) :: iprof
sll_real64,    dimension(:) :: aele(*), xmass(*)
sll_int32 :: i1, i2, i3, idiag1, idiag2, idiag3, ind

!--- 1.1 --- Rangement des termes diagonaux ---------------------------
 
idiag1 = iprof(i1+1)
idiag2 = iprof(i2+1)
idiag3 = iprof(i3+1)

xmass(idiag1) = xmass(idiag1) + aele(1)
xmass(idiag2) = xmass(idiag2) + aele(4)
xmass(idiag3) = xmass(idiag3) + aele(6)

!--- 1.2 --- Rangement des autres termes ------------------------------
!           (on ne stocke que les aij tels que i>j)
 
if (i1 < i2) then 
   ind=idiag2+i1-i2
else
   ind=idiag1+i2-i1
end if 
 
xmass(ind)=xmass(ind)+aele(2)

if (i1 < i3) then 
   ind=idiag3+i1-i3
else
   ind=idiag1+i3-i1
end if 
 
xmass(ind)=xmass(ind)+aele(3)

if (i2 < i3) then 
   ind=idiag3+i2-i3
else
   ind=idiag2+i3-i2
end if 

xmass(ind)=xmass(ind)+aele(5)

!----------------------------------------------------------------------

end subroutine asblp

 
!Function: m1p
!                                             
!     faire l'operation : yvect =  xmors.xvect 
!     ou  "xmors" est une matrice "morse"       
!     xvect,yvect  des vecteurs                  
!
!Entrees:
!          xmors        -     matrice  "morse"          
!          xvect        -     vecteur operande           
!          mors1,mors2  -     tableaux descriptifs de     
!                             la matrice "morse"           
!          nlign        -     nombre de lignes des matrices 
!                                                      
!Sorties:
!          yvect    -   vecteur resultat            
!                                                  
!Auteur:
!     J. Segre - Version 1.0  Decembre 1989
subroutine m1p(xmors,mors1,mors2,xvect,nlign,yvect)

sll_real64    :: xmors(*),xvect(*),yvect(*)
sll_int32 :: mors1(*),mors2(*)
sll_int32 :: il, nlign, noeui

do il=1,nlign

   noeui=mors1(il+1)-mors1(il)
  
   select case (noeui)
   case(6)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+   &
                xmors(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+   &
                xmors(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))
  
   case(5)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+   &
                xmors(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))

   case(7)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+   &
                xmors(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+   &
                xmors(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+   &
                xmors(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))

   case(4)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))
 
   case(8)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+   &
                xmors(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+   &
                xmors(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+   &
                xmors(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))+   &
                xmors(mors1(il+1)-7)*xvect(mors2(mors1(il+1)-7))

   case(3)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))

   case(9)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+   &
                xmors(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+   &
                xmors(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+   &
                xmors(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))+   &
                xmors(mors1(il+1)-7)*xvect(mors2(mors1(il+1)-7))+   &
                xmors(mors1(il+1)-8)*xvect(mors2(mors1(il+1)-8))
 
   case(10)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+   &
                xmors(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+   &
                xmors(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+   &
                xmors(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))+   &
                xmors(mors1(il+1)-7)*xvect(mors2(mors1(il+1)-7))+   &
                xmors(mors1(il+1)-8)*xvect(mors2(mors1(il+1)-8))+   &
                xmors(mors1(il+1)-9)*xvect(mors2(mors1(il+1)-9))
  
   case(11)
      yvect(il)=xmors(mors1(il+1)   )*xvect(mors2(mors1(il+1)   ))+ &
                xmors(mors1(il+1)- 1)*xvect(mors2(mors1(il+1)- 1))+ &
                xmors(mors1(il+1)- 2)*xvect(mors2(mors1(il+1)- 2))+ &
                xmors(mors1(il+1)- 3)*xvect(mors2(mors1(il+1)- 3))+ &
                xmors(mors1(il+1)- 4)*xvect(mors2(mors1(il+1)- 4))+ &
                xmors(mors1(il+1)- 5)*xvect(mors2(mors1(il+1)- 5))+ &
                xmors(mors1(il+1)- 6)*xvect(mors2(mors1(il+1)- 6))+ &
                xmors(mors1(il+1)- 7)*xvect(mors2(mors1(il+1)- 7))+ &
                xmors(mors1(il+1)- 8)*xvect(mors2(mors1(il+1)- 8))+ &
                xmors(mors1(il+1)- 9)*xvect(mors2(mors1(il+1)- 9))+ &
                xmors(mors1(il+1)-10)*xvect(mors2(mors1(il+1)-10))

   case(12)
      yvect(il)=xmors(mors1(il+1)   )*xvect(mors2(mors1(il+1)   ))+ &
                xmors(mors1(il+1)- 1)*xvect(mors2(mors1(il+1)- 1))+ &
                xmors(mors1(il+1)- 2)*xvect(mors2(mors1(il+1)- 2))+ &
                xmors(mors1(il+1)- 3)*xvect(mors2(mors1(il+1)- 3))+ &
                xmors(mors1(il+1)- 4)*xvect(mors2(mors1(il+1)- 4))+ &
                xmors(mors1(il+1)- 5)*xvect(mors2(mors1(il+1)- 5))+ &
                xmors(mors1(il+1)- 6)*xvect(mors2(mors1(il+1)- 6))+ &
                xmors(mors1(il+1)- 7)*xvect(mors2(mors1(il+1)- 7))+ &
                xmors(mors1(il+1)- 8)*xvect(mors2(mors1(il+1)- 8))+ &
                xmors(mors1(il+1)- 9)*xvect(mors2(mors1(il+1)- 9))+ &
                xmors(mors1(il+1)-10)*xvect(mors2(mors1(il+1)-10))+ &
                xmors(mors1(il+1)-11)*xvect(mors2(mors1(il+1)-11))

   end select 

end do

end subroutine m1p



!  Effectuer le calcul des composantes Ex et Ey du    
!  ---    champ electrique a partir du potentiel au noeud.
!         On appelle ca un lissage.                        
subroutine poliss(this, phi, ex, ey)

type(sll_triangular_poisson_2d), intent(inout) :: this
sll_real64, dimension(:),        intent(in)    :: phi
sll_real64, dimension(:),        intent(out)   :: ex
sll_real64, dimension(:),        intent(out)   :: ey

sll_int32 :: is, ic, nbc, iac
LOGICAL :: lerr

lerr=.FALSE.

! --- 1.0 --- Calcul des termes individuels des seconds membres --------
  
do ic=1,this%mesh%nbtcot
   this%vtanty(ic)=(phi(this%mesh%nuvac(1,ic))-phi(this%mesh%nuvac(2,ic)))/this%mesh%xlcod(ic)
end do

do ic=1,this%mesh%nbtcot
   this%vtantx(ic)=this%vtanty(ic)*this%mesh%vtaux(ic)
end do

do ic=1,this%mesh%nbtcot
   this%vtanty(ic)=this%vtanty(ic)*this%mesh%vtauy(ic)
end do

! --- 2.0 --- Calcul des seconds membres -------------------------------
  
do is=1,this%mesh%num_nodes

   iac=this%mesh%nbcov(is)+1
   nbc=this%mesh%nbcov(is+1)-this%mesh%nbcov(is)

   if (nbc == 6) then

      this%sv1(is) = this%vtantx(this%mesh%nugcv(iac  ))  &
                   + this%vtantx(this%mesh%nugcv(iac+1))  &
                   + this%vtantx(this%mesh%nugcv(iac+2))  &
                   + this%vtantx(this%mesh%nugcv(iac+3))  &
                   + this%vtantx(this%mesh%nugcv(iac+4))  &
                   + this%vtantx(this%mesh%nugcv(iac+5))  

      this%sv2(is) = this%vtanty(this%mesh%nugcv(iac  ))  &
                   + this%vtanty(this%mesh%nugcv(iac+1))  &
                   + this%vtanty(this%mesh%nugcv(iac+2))  &
                   + this%vtanty(this%mesh%nugcv(iac+3))  &
                   + this%vtanty(this%mesh%nugcv(iac+4))  &
                   + this%vtanty(this%mesh%nugcv(iac+5))  

   else if (nbc == 5) then

      this%sv1(is) = this%vtantx(this%mesh%nugcv(iac  ))  &
                   + this%vtantx(this%mesh%nugcv(iac+1))  &
                   + this%vtantx(this%mesh%nugcv(iac+2))  &
                   + this%vtantx(this%mesh%nugcv(iac+3))  &
                   + this%vtantx(this%mesh%nugcv(iac+4))

      this%sv2(is) = this%vtanty(this%mesh%nugcv(iac  ))  &
                   + this%vtanty(this%mesh%nugcv(iac+1))  &
                   + this%vtanty(this%mesh%nugcv(iac+2))  &
                   + this%vtanty(this%mesh%nugcv(iac+3))  &
                   + this%vtanty(this%mesh%nugcv(iac+4))

   else if (nbc == 7) then

      this%sv1(is) = this%vtantx(this%mesh%nugcv(iac  ))  &
                   + this%vtantx(this%mesh%nugcv(iac+1))  &
                   + this%vtantx(this%mesh%nugcv(iac+2))  &
                   + this%vtantx(this%mesh%nugcv(iac+3))  &
                   + this%vtantx(this%mesh%nugcv(iac+4))  &
                   + this%vtantx(this%mesh%nugcv(iac+5))  &
                   + this%vtantx(this%mesh%nugcv(iac+6))

      this%sv2(is) = this%vtanty(this%mesh%nugcv(iac  ))  &
                   + this%vtanty(this%mesh%nugcv(iac+1))  &
                   + this%vtanty(this%mesh%nugcv(iac+2))  &
                   + this%vtanty(this%mesh%nugcv(iac+3))  &
                   + this%vtanty(this%mesh%nugcv(iac+4))  &
                   + this%vtanty(this%mesh%nugcv(iac+5))  &
                   + this%vtanty(this%mesh%nugcv(iac+6))

   else if (nbc == 4) then

      this%sv1(is) = this%vtantx(this%mesh%nugcv(iac  ))  &
                   + this%vtantx(this%mesh%nugcv(iac+1))  &
                   + this%vtantx(this%mesh%nugcv(iac+2))  &
                   + this%vtantx(this%mesh%nugcv(iac+3))

      this%sv2(is) = this%vtanty(this%mesh%nugcv(iac  ))  &
                   + this%vtanty(this%mesh%nugcv(iac+1))  &
                   + this%vtanty(this%mesh%nugcv(iac+2))  &
                   + this%vtanty(this%mesh%nugcv(iac+3))

   else if (nbc == 8) then

      this%sv1(is) = this%vtantx(this%mesh%nugcv(iac  ))  &
                   + this%vtantx(this%mesh%nugcv(iac+1))  &
                   + this%vtantx(this%mesh%nugcv(iac+2))  &
                   + this%vtantx(this%mesh%nugcv(iac+3))  &
                   + this%vtantx(this%mesh%nugcv(iac+4))  &
                   + this%vtantx(this%mesh%nugcv(iac+5))  &
                   + this%vtantx(this%mesh%nugcv(iac+6))  &
                   + this%vtantx(this%mesh%nugcv(iac+7))  

      this%sv2(is) = this%vtanty(this%mesh%nugcv(iac  ))  &
                   + this%vtanty(this%mesh%nugcv(iac+1))  &
                   + this%vtanty(this%mesh%nugcv(iac+2))  &
                   + this%vtanty(this%mesh%nugcv(iac+3))  &
                   + this%vtanty(this%mesh%nugcv(iac+4))  &
                   + this%vtanty(this%mesh%nugcv(iac+5))  &
                   + this%vtanty(this%mesh%nugcv(iac+6))  &
                   + this%vtanty(this%mesh%nugcv(iac+7))  

   else if (nbc == 3) then

      this%sv1(is) = this%vtantx(this%mesh%nugcv(iac  ))  &
                   + this%vtantx(this%mesh%nugcv(iac+1))  &
                   + this%vtantx(this%mesh%nugcv(iac+2))

      this%sv2(is) = this%vtanty(this%mesh%nugcv(iac  ))  &
                   + this%vtanty(this%mesh%nugcv(iac+1))  &
                   + this%vtanty(this%mesh%nugcv(iac+2))

   else if (nbc == 9) then

      this%sv1(is) = this%vtantx(this%mesh%nugcv(iac  ))  &
                   + this%vtantx(this%mesh%nugcv(iac+1))  &
                   + this%vtantx(this%mesh%nugcv(iac+2))  &
                   + this%vtantx(this%mesh%nugcv(iac+3))  &
                   + this%vtantx(this%mesh%nugcv(iac+4))  &
                   + this%vtantx(this%mesh%nugcv(iac+5))  &
                   + this%vtantx(this%mesh%nugcv(iac+6))  &
                   + this%vtantx(this%mesh%nugcv(iac+7))  &
                   + this%vtantx(this%mesh%nugcv(iac+8))

      this%sv2(is) = this%vtanty(this%mesh%nugcv(iac  ))  &
                   + this%vtanty(this%mesh%nugcv(iac+1))  &
                   + this%vtanty(this%mesh%nugcv(iac+2))  &
                   + this%vtanty(this%mesh%nugcv(iac+3))  &
                   + this%vtanty(this%mesh%nugcv(iac+4))  &
                   + this%vtanty(this%mesh%nugcv(iac+5))  &
                   + this%vtanty(this%mesh%nugcv(iac+6))  &
                   + this%vtanty(this%mesh%nugcv(iac+7))  &
                   + this%vtanty(this%mesh%nugcv(iac+8))

   else if (nbc == 2) then

      this%sv1(is) = this%vtantx(this%mesh%nugcv(iac  ))  &
                   + this%vtantx(this%mesh%nugcv(iac+1))

      this%sv2(is) = this%vtanty(this%mesh%nugcv(iac  ))  &
                   + this%vtanty(this%mesh%nugcv(iac+1))

   else if (nbc == 10) then

      this%sv1(is) = this%vtantx(this%mesh%nugcv(iac  ))  &
                   + this%vtantx(this%mesh%nugcv(iac+1))  &
                   + this%vtantx(this%mesh%nugcv(iac+2))  &
                   + this%vtantx(this%mesh%nugcv(iac+3))  &
                   + this%vtantx(this%mesh%nugcv(iac+4))  &
                   + this%vtantx(this%mesh%nugcv(iac+5))  &
                   + this%vtantx(this%mesh%nugcv(iac+6))  &
                   + this%vtantx(this%mesh%nugcv(iac+7))  &
                   + this%vtantx(this%mesh%nugcv(iac+8))  &
                   + this%vtantx(this%mesh%nugcv(iac+9))  

      this%sv2(is) = this%vtanty(this%mesh%nugcv(iac  ))  &
                   + this%vtanty(this%mesh%nugcv(iac+1))  &
                   + this%vtanty(this%mesh%nugcv(iac+2))  &
                   + this%vtanty(this%mesh%nugcv(iac+3))  &
                   + this%vtanty(this%mesh%nugcv(iac+4))  &
                   + this%vtanty(this%mesh%nugcv(iac+5))  &
                   + this%vtanty(this%mesh%nugcv(iac+6))  &
                   + this%vtanty(this%mesh%nugcv(iac+7))  &
                   + this%vtanty(this%mesh%nugcv(iac+8))  &
                   + this%vtanty(this%mesh%nugcv(iac+9)) 

   else if (nbc == 11) then

      this%sv1(is) = this%vtantx(this%mesh%nugcv(iac  ))  &
                   + this%vtantx(this%mesh%nugcv(iac+1))  &
                   + this%vtantx(this%mesh%nugcv(iac+2))  &
                   + this%vtantx(this%mesh%nugcv(iac+3))  &
                   + this%vtantx(this%mesh%nugcv(iac+4))  &
                   + this%vtantx(this%mesh%nugcv(iac+5))  &
                   + this%vtantx(this%mesh%nugcv(iac+6))  &
                   + this%vtantx(this%mesh%nugcv(iac+7))  &
                   + this%vtantx(this%mesh%nugcv(iac+8))  &
                   + this%vtantx(this%mesh%nugcv(iac+9))  &
                   + this%vtantx(this%mesh%nugcv(iac+10))

      this%sv2(is) = this%vtanty(this%mesh%nugcv(iac  ))  &
                   + this%vtanty(this%mesh%nugcv(iac+1))  &
                   + this%vtanty(this%mesh%nugcv(iac+2))  &
                   + this%vtanty(this%mesh%nugcv(iac+3))  &
                   + this%vtanty(this%mesh%nugcv(iac+4))  &
                   + this%vtanty(this%mesh%nugcv(iac+5))  &
                   + this%vtanty(this%mesh%nugcv(iac+6))  &
                   + this%vtanty(this%mesh%nugcv(iac+7))  &
                   + this%vtanty(this%mesh%nugcv(iac+8))  &
                   + this%vtanty(this%mesh%nugcv(iac+9))  &
                   + this%vtanty(this%mesh%nugcv(iac+10))

   else if (nbc == 12) then

      this%sv1(is) = this%vtantx(this%mesh%nugcv(iac  ))  &
                   + this%vtantx(this%mesh%nugcv(iac+1))  &
                   + this%vtantx(this%mesh%nugcv(iac+2))  &
                   + this%vtantx(this%mesh%nugcv(iac+3))  &
                   + this%vtantx(this%mesh%nugcv(iac+4))  &
                   + this%vtantx(this%mesh%nugcv(iac+5))  &
                   + this%vtantx(this%mesh%nugcv(iac+6))  &
                   + this%vtantx(this%mesh%nugcv(iac+7))  &
                   + this%vtantx(this%mesh%nugcv(iac+8))  &
                   + this%vtantx(this%mesh%nugcv(iac+9))  &
                   + this%vtantx(this%mesh%nugcv(iac+10)) &
                   + this%vtantx(this%mesh%nugcv(iac+11))

      this%sv2(is) = this%vtanty(this%mesh%nugcv(iac  ))  &
                   + this%vtanty(this%mesh%nugcv(iac+1))  &
                   + this%vtanty(this%mesh%nugcv(iac+2))  &
                   + this%vtanty(this%mesh%nugcv(iac+3))  &
                   + this%vtanty(this%mesh%nugcv(iac+4))  &
                   + this%vtanty(this%mesh%nugcv(iac+5))  &
                   + this%vtanty(this%mesh%nugcv(iac+6))  &
                   + this%vtanty(this%mesh%nugcv(iac+7))  &
                   + this%vtanty(this%mesh%nugcv(iac+8))  &
                   + this%vtanty(this%mesh%nugcv(iac+9))  &
                   + this%vtanty(this%mesh%nugcv(iac+10)) &
                   + this%vtanty(this%mesh%nugcv(iac+11))
   else

      lerr=.TRUE.

   end if

end do

if (lerr) then
   write(6,900)
   stop
end if

! --- 3.0 --- Resolution des systemes lineaires 2*2 --------------------

do is=1,this%mesh%num_nodes
   ex(is)=this%mesh%xmal2(is)*this%sv1(is)-this%mesh%xmal3(is)*this%sv2(is)
   ey(is)=this%mesh%xmal1(is)*this%sv2(is)-this%mesh%xmal3(is)*this%sv1(is)
end do


! --- 9.0 --- Formats --------------------------------------------------
    
900 format(//10x,'Erreur dans POLISS'   &
            /10x,'On a trouve plus de 12 cotes')

end subroutine poliss

end module sll_m_tri_poisson
