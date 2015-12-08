! module pour gérer et résoudre le laplacien sur un domaine
! difféomorphe à un carré
! la solution est le potentiel: phi
! le terme source est la charge: rho
! conditions aux limites de Dirichlet données par la fonction
! "potexact"
module sll_m_lobalap
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  use sll_m_map_function, only: &
    map

  implicit none

  public :: &
    assemb, &
    assemb_rhs, &
    compute_electric_field, &
    compute_phi, &
    computelu, &
    init, &
    plotgmsh, &
    release

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ordre de l'interpolation élément fini
  integer :: order
  ! nombre de noeuds locaux dans chaque élément
  integer :: nloc

  ! tableaux pour intégration de Gauss-Lobatto
  ! points et poids
  sll_real64,dimension(:),allocatable :: xpg,wpg
  ! dérivées des polynômes de Lagrange aux points de gauss
  sll_real64,dimension(:,:),allocatable :: dlag

  ! variable pour vérifier que les initialisations ont eu lieu
  logical :: is_init=.false.

  ! nombre d'éléments en x et en y
  integer :: nx,ny
  ! nombre d'éléments
  integer :: nel
  ! nombre d'inconnues en x, en y, en tout
  integer :: neqx,neqy,neq
  ! nombre de conditions de bord
  integer :: nbc

  ! noeuds du maillage
  sll_real64,dimension(:,:),allocatable :: node
  ! connectivité
  integer,dimension(:,:),allocatable :: connec

  ! tableaux pour le stockage morse de la matrice
  ! indices de ligne et de colonne des valeurs non nulles
  !integer,dimension(:),allocatable :: icol,irow
  ! valeurs correspondantes de la matrice
  !sll_real64,dimension(:),allocatable :: matval

  ! tableaux pour le stockage skyline de la matrice
  ! profil, début des colonnes
  integer,dimension(:),allocatable :: prof,kld
  sll_real64,dimension(:),allocatable :: vdiag,vsup,vinf
  ! taille de ces tableaux
  integer :: nsky
  ! grand pivot pour les cl
  sll_real64,parameter :: big=1.e20_f64

  ! vecteur pour le second membre et la solution
  sll_real64,dimension(:),allocatable :: phi,rho

  ! liste des noeuds du bord
  integer,dimension(:),allocatable :: indexbc

contains

  ! fonction donnant le potentiel exact (debug) et/ou les conditions aux limites
  function potexact(x,y)
    implicit none
    sll_real64,intent(in) :: x,y
    sll_real64 :: potexact
    !potexact=x*x+y*y
    potexact=0.0_f64+x-x+y-y
  end function potexact

  ! fonction donnant le terme source
  function source(x,y)
    implicit none
    sll_real64,intent(in) :: x,y
    sll_real64 :: source
    source=-4.0_f64+x-x+y-y
  end function source

  ! fonction qui envoie le carré [0,1]x[0,1] sur le vrai domaine de calcul
  ! variables de référence: (u,v)
  ! variables physiques: (x,y)
  ! autres données calculées:
  ! jac, invjac, det: jacobienne, son inverse et déterminant de la jacobienne
!  ! subroutine map(u,v,x,y,jac,invjac,det)
!  ! pour l'instant on n'utilise pas la jacobienne
!  subroutine map(u,v,x,y)
!    implicit none
!    sll_real64,intent(in) :: u,v
!    sll_real64,intent(out) :: x,y
!    sll_real64 :: jac(2,2),invjac(2,2),det
!    sll_real64,parameter :: pi=4*atan(1._f64)
!
!    x=(1+u)*(1+v)*cos(pi*v)
!    y=(1+u)*sin(pi*v)
!
!    x=u
!    y=v
!
!    !x=u
!    !y=v
!    ! non utilisé
!    jac=0
!    jac(1,1)=1
!    jac(2,2)=1
!    invjac=jac
!    det=1
!
!  end subroutine map

  ! remplissage des tableaux de pg
  ! et de polynômes de Lagrange
  ! dlag(i,j): dérivé polynôme i au point j
  subroutine init_gauss()
    implicit none


    write(*,*) 'Init tableaux points de Gauss...'

    allocate(xpg(order+1))
    allocate(wpg(order+1))
    allocate(dlag(order+1,order+1))

    select case(order)
    case(1)
       xpg(1)=0.0_f64
       xpg(2)=1.0_f64
       wpg(1)=0.5_f64
       wpg(2)=0.5_f64
       dlag(1,1) = -1.0_f64
       dlag(1,2) = -1.0_f64
       dlag(2,1) = 1.0_f64
       dlag(2,2) = 1.0_f64

    case(2)
       xpg(1)=0.0_f64
       xpg(2)=0.5_f64
       xpg(3)=1.0_f64
       wpg(1)=1.0_f64/6.0_f64
       wpg(2)=4.0_f64/6.0_f64
       wpg(3)=1.0_f64/6.0_f64
       dlag(1,1) = -0.300000000000000000000000000000e1_f64
       dlag(1,2) = -0.100000000000000000000000000000e1_f64
       dlag(1,3) = 0.100000000000000000000000000000e1_f64
       dlag(2,1) = 0.400000000000000000000000000000e1_f64
       dlag(2,3) = -0.400000000000000000000000000000e1_f64
       dlag(3,1) = -0.100000000000000000000000000000e1_f64
       dlag(3,2) = 0.100000000000000000000000000000e1_f64
       dlag(3,3) = 0.300000000000000000000000000000e1_f64
    case(3)
       xpg(1)=0.0_f64
       xpg(2)=(1.0_f64-sqrt(1.0_f64/5.0_f64))/2.0_f64
       xpg(3)=(1.0_f64+sqrt(1.0_f64/5.0_f64))/2.0_f64
       xpg(4)=1.0_f64
       wpg(1)=1.0_f64/12.0_f64
       wpg(2)=5.0_f64/12.0_f64
       wpg(3)=5.0_f64/12.0_f64
       wpg(4)=1.0_f64/12.0_f64
       dlag(1,1) = -0.599999999999999999999999999998e1_f64
       dlag(1,2) = -0.161803398874989484820458683436e1_f64
       dlag(1,3) = 0.618033988749894848204586834362e0_f64
       dlag(1,4) = -0.999999999999999999999999999994e0_f64
       dlag(2,1) = 0.809016994374947424102293417177e1_f64
       dlag(2,2) = -0.1D-28
       dlag(2,3) = -0.223606797749978969640917366872e1_f64
       dlag(2,4) = 0.309016994374947424102293417184e1_f64
       dlag(3,1) = -0.309016994374947424102293417182e1_f64
       dlag(3,2) = 0.223606797749978969640917366872e1_f64
       dlag(3,3) = 0.1D-28
       dlag(3,4) = -0.809016994374947424102293417177e1_f64
       dlag(4,1) = 0.999999999999999999999999999994e0_f64
       dlag(4,2) = -0.618033988749894848204586834362e0_f64
       dlag(4,3) = 0.161803398874989484820458683436e1_f64
       dlag(4,4) = 0.599999999999999999999999999998e1_f64
    case(4)
       xpg(1)=0.0_f64
       xpg(2)=(1.0_f64-sqrt(3.0_f64/7.0_f64))/2.0_f64
       xpg(3)=0.5_f64
       xpg(4)=(1.0_f64+sqrt(3.0_f64/7.0_f64))/2.0_f64
       xpg(5)=1.0_f64
       wpg(1)=1.0_f64/20.0_f64
       wpg(2)=49.0_f64/180.0_f64
       wpg(3)=32.0_f64/90.0_f64
       wpg(4)=49.0_f64/180.0_f64
       wpg(5)=1.0_f64/20.0_f64
       dlag(1,1) = -0.100000000000000000000000000000e2_f64
       dlag(1,2) = -0.248198050606196571569743868439e1_f64
       dlag(1,3) = 0.750000000000000000000000000000e0_f64
       dlag(1,4) = -0.518019493938034284302561315633e0_f64
       dlag(1,5) = 0.100000000000000000000000000002e1_f64
       dlag(2,1) = 0.135130049774484800076860550594e2_f64
       dlag(2,2) = -0.2e-28_f64
       dlag(2,3) = -0.267316915539090667050969419638e1_f64
       dlag(2,4) = 0.152752523165194666886268239794e1_f64
       dlag(2,5) = -0.282032835588485332564727827395e1_f64
       dlag(3,1) = -0.533333333333333333333333333328e1_f64
       dlag(3,2) = 0.349148624377587810025755976661e1_f64
       dlag(3,3) = 0.2e-28_f64
       dlag(3,4) = -0.349148624377587810025755976659e1_f64
       dlag(3,5) = 0.533333333333333333333333333332e1_f64
       dlag(4,1) = 0.282032835588485332564727827398e1_f64
       dlag(4,2) = -0.152752523165194666886268239791e1_f64
       dlag(4,3) = 0.267316915539090667050969419635e1_f64
       dlag(4,4) = -0.4e-28_f64
       dlag(4,5) = -0.135130049774484800076860550595e2_f64
       dlag(5,1) = -0.100000000000000000000000000000e1_f64
       dlag(5,2) = 0.518019493938034284302561315635e0_f64
       dlag(5,3) = -0.750000000000000000000000000006e0_f64
       dlag(5,4) = 0.248198050606196571569743868439e1_f64
       dlag(5,5) = 0.100000000000000000000000000000e2_f64
    case default
       write(*,*) 'pas prévu...'
       stop
    end select



  end subroutine init_gauss


  ! tracé de la solution avec gmsh
  subroutine plotgmsh()
    implicit none
    integer,parameter :: nlocmax=25
    integer :: typelem
    integer,dimension(4) :: invpermut1=(/1,2,4,3/)
    integer,dimension(9) :: invpermut2=(/1,5,2,8,9,6,4,7,3/)
    integer,dimension(16) :: invpermut3=(/1,5,6,2,12,13,14,7,11,16,15,8,4,10,9,3/)
    integer :: invpermut(nloc),permut(nloc),i,ii

    write(*,*) 'Tracé gmsh...'

    if (order.gt.3) then
       write(*,*) 'ordre non prévu'
       stop
    end if

    select case(order)
    case(1)
       invpermut(1:nloc)=invpermut1
       typelem=3
    case(2)
       invpermut(1:nloc)=invpermut2
       typelem=10
    case(3)
       invpermut(1:nloc)=invpermut3
       typelem=36
    case default
       write(*,*) 'ordre non prévu'
       stop
    end select

    do i=1,nloc
       permut(invpermut(i))=i
    end do
       

    open(101,file='sll_m_lobalap.msh')
    write(101,'(A)') '$MeshFormat'
    write(101,'(A)') '2 0 8'
    write(101,'(A)') '$EndMeshFormat'
    write(101,'(A)') '$Nodes'
    write(101,*) neq
    do i=1,neq
       write(101,*) i,node(1,i),node(2,i),0._f64
    end do
    write(101,'(A)') '$EndNodes'
    write(101,'(A)') '$Elements'
    write(101,*) nel
    do i=1,nel
       write(101,*) i,typelem,0,(connec(permut(ii),i),ii=1,nloc)
    end do
    write(101,'(A)') '$EndElements'
    write(101,'(A)') '$NodeData'
    write(101,*) 1
    write(101,*) 'PHI'
    write(101,*) 1
    write(101,*) 0
    write(101,*) 3
    write(101,*) 0
    write(101,*) 1
    write(101,*) neq
    do i=1,neq
       write(101,*) i,phi(i)-potexact(node(1,i),node(2,i))
    end do
    write(101,'(A)') '$EndNodeData'

    close(101)

  end subroutine plotgmsh

  ! construction du maillage
  subroutine build_mesh()
    implicit none
    integer    :: iix,ix,iiy,iy,inoloc,ino,iel
    sll_real64 :: du,dv,u,v,x,y,eps

    ! construit les tableaux de pg
    call init_gauss()

    write(*,*) 'Construction du maillage...'
    write(*,*) 'Elements:',nel,' Noeuds:',neq
    allocate(node(2,neq))
    allocate(connec(nloc,neq))

    du=1._f64/nx
    dv=1._f64/ny

    ! compteur pour les noeuds du bord
    nbc=0
    ! tolérance
    eps=1.0e-12_f64

    ! construit la connectivité éléments --> noeuds
    ! boucle sur les éléments
    do ix=0,nx-1
       do iy=0,ny-1
          iel=1+ix+iy*nx
          !boucle sur les noeuds de l'éléments
          do iix=0,order
             do iiy=0,order
                ! numéro local dans l'élément
                inoloc=1+iix+(order+1)*iiy
                ! numéro global dans le maillage
                ino=1+iix+order*ix+(order*nx+1)*(iiy+order*iy)
                connec(inoloc,iel)=ino
                ! coordonnées
                u=ix*du+du*xpg(iix+1)
                v=iy*dv+dv*xpg(iiy+1)
                node(1,ino)=u
                node(2,ino)=v
             end do
          end do
       end do
    end do
    write(*,*) " ino = ", ino

    write(*,*) 'Repérage noeuds du bord'

    ! repère les noeuds du bord
    ! première passe pour compter
    do ino=1,neq
       u=node(1,ino)
       v=node(2,ino)
       if (abs(u*(1-u)*v*(1-v)).lt.eps) nbc=nbc+1
    end do
    
    allocate(indexbc(nbc))

    ! deuxième passe pour remplir
    ! et calculer les coords physiques
    ! des noeuds
    nbc=0
    do ino=1,neq
       u=node(1,ino)
       v=node(2,ino)
       if (abs(u*(1-u)*v*(1-v)).lt.eps) then
          nbc=nbc+1
          indexbc(nbc)=ino
       end if
       call map(u,v,x,y)
       node(1,ino)=x
       node(2,ino)=y
    end do

!!$    write(*,*) 'Points de bord: ',nbc
!!$   do ino=1,nbc
!!$      write(*,*) indexbc(ino)
!!$   end do


    

  end subroutine build_mesh

  ! alloue et prépare les données pour le calcul
  subroutine init(nx0,ny0,order0)
    implicit none
    integer,intent(in) :: nx0,ny0,order0
    integer :: ino,iel,i,ii,j,jj
    nx=nx0
    ny=ny0
    order=order0
    nloc=(order+1)*(order+1)
    is_init=.true.

    ! nombre d'inconnues
    neqx=nx*order+1
    neqy=ny*order+1
    neq=neqx*neqy
    nel=nx*ny

    ! construction du maillage
    ! et des listes de conditions aux limites
    call build_mesh()

    write(*,*) 'Allocation solution et source...'

    ! allocation solution et source
    allocate(phi(neq))
    allocate(rho(neq))

    ! une solution bidon pour tester la visu
    do ino=1,neq
       phi(ino)=potexact(node(1,ino),node(2,ino))
    end do

    ! initialisation du second pour l'assemblage
    rho=0.0_f64

    write(*,*) 'Calcul structure matrice creuse...'

    ! préparation de la matrice skyline
    allocate(prof(neq))
    allocate(kld(neq+1))

    ! construction du profil
    prof=0
    do iel=1,nel
       do ii=1,nloc
          do jj=1,nloc
             i=connec(ii,iel)
             j=connec(jj,iel)
             prof(j)=max(prof(j),j-i)
             prof(i)=max(prof(i),i-j)
          enddo
       enddo
    enddo
    ! tableau des débuts de colonnes
    kld(1)=1
    do i=1,neq
       kld(i+1)=kld(i)+prof(i)
    enddo
    
    nsky=kld(neq+1)-1

    allocate(vdiag(neq))
    allocate(vinf(nsky))
    allocate(vsup(nsky))

    vdiag=0.0_f64
    vinf=0.0_f64
    vsup=0.0_f64
    

  end subroutine init

  ! symbole de kronecker delta
  function delta(i,j)
    implicit none
    integer    :: i,j
    sll_real64 :: delta
    if (i.eq.j) then
       delta=1.0_f64
    else
       delta=0.0_f64
    end if
  end function delta

  ! gradient de la fonction de base iloc au point de gauss ipg
  ! sur l'élément de référence dans les variables de référence
  ! on utilise le fait que les fonctions de base sont des produits 
  ! tensoriels de polynômes de Lagrange associés aux 
  ! point de Guass-Lobatto
  subroutine gradpg(iloc,ipg,grad,poids)
    implicit none
    integer,    intent(in)  :: iloc,ipg
    sll_real64, intent(out) :: grad(2),poids
    integer :: ilocx,ilocy,ipgx,ipgy
    ! indices locaux du pg
    ! dans les directions de référence 
    ipgx=mod(ipg-1,order+1)+1
    ipgy=(ipg-1)/(order+1)+1

    ! indices de la fonction de base
    ! suivant u et v dans l'élément de référence
    ilocx=mod(iloc-1,order+1)+1
    ilocy=(iloc-1)/(order+1)+1
 
    ! gradient de la fonction de base au pg 
    ! dans les variables de référence
    grad(1)=dlag(ilocx,ipgx)*delta(ilocy,ipgy)
    grad(2)=dlag(ilocy,ipgy)*delta(ilocx,ipgx)

    poids=wpg(ipgx)*wpg(ipgy)

  end subroutine gradpg


  ! calcul du champ électrique
  subroutine compute_electric_field(dg_ex,dg_ey)
    implicit none
    ! matrice locale
    sll_real64 :: jac(2,2),cojac(2,2),det
    sll_real64 :: gradref(2),xg,yg
    sll_real64 :: grad(2),dxy(2),poids
    integer :: iel,ipg,ii,jj,ig,ib,iib,ielx,iely
    sll_real64, dimension(order+1,order+1,nx,ny) :: dg_ex,dg_ey
    
    ! assemblage de la matrice de rigidité
    ! et du second membre

    write(*,*) 'Calcul champ éléectrique...'

    dg_ex=0.0_f64
    dg_ey=0.0_f64

    ! boucle sur les éléments
    do iel=1,nel
       ! boucle sur les pg pour le calcul du champ
       ielx=mod(iel-1,nx)+1
       iely=(iel-1)/nx+1
       do ipg=1,nloc
          ! numéros locaux du pg
          ii=mod(ipg-1,order+1)+1
          jj=(ipg-1)/(order+1)+1
          ! numéro global du pg
          ig=connec(ipg,iel)
          xg=node(1,ig)
          yg=node(2,ig)
          ! calcul de la jacobienne au pg
          ! on pourra le faire directement avec map plus tard
          jac=0.0_f64
          ! boucle sur les fonctions de base d'interpolation
          ! pour construire la jacobienne de la transformation
          ! géométrique
          do iib=1,nloc
             ! gradient de la fonction de base au pg 
             ! dans les variables de référence
             ! calcul du poids de gauss
             call gradpg(iib,ipg,dxy,poids)
             ib=connec(iib,iel)
             ! contribution à la jacobienne
             jac(1,1)=jac(1,1)+node(1,ib)*dxy(1)
             jac(1,2)=jac(1,2)+node(1,ib)*dxy(2)
             jac(2,1)=jac(2,1)+node(2,ib)*dxy(1)
             jac(2,2)=jac(2,2)+node(2,ib)*dxy(2)
          end do
          !write(*,*) ipg,jac
          ! déterminant et transposée de l'inverse
          det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
          cojac(1,1)=jac(2,2)/det
          cojac(2,2)=jac(1,1)/det
          cojac(1,2)=-jac(2,1)/det
          cojac(2,1)=-jac(1,2)/det

          ! nouvelle boucle sur les noeuds locaux de l'élément iel
          ! pour calculer le gradient du potentiel
          do iib=1,nloc
             ! gradients dans les variables de référence
             call gradpg(iib,ipg,gradref,poids)
             ib=connec(iib,iel)
             grad=matmul(cojac,gradref) 
             dg_ex(ii,jj,ielx,iely)=dg_ex(ii,jj,ielx,iely)+phi(ib)*grad(1)
             dg_ey(ii,jj,ielx,iely)=dg_ey(ii,jj,ielx,iely)+phi(ib)*grad(2)
          end do
       end do
    end do

       
  end subroutine compute_electric_field





  ! assemblage de la matrice élément fini et des conditions aux limites
  subroutine assemb()
    implicit none
    ! matrice locale
    sll_real64 :: jac(2,2),cojac(2,2),det
    sll_real64 :: gradref_i(2),gradref_j(2),xg,yg
    sll_real64 :: grad_i(2),grad_j(2),dxy(2),v,poids,vf
    integer :: iel,ipg,i,ii,j,jj,ig,ib,iib
    
    ! assemblage de la matrice de rigidité
    ! et du second membre

    write(*,*) 'Assemblage matrice...'

    ! boucle sur les éléments
    do iel=1,nel
       ! boucle sur les pg pour intégrer
       do ipg=1,nloc
          ! numéro global du pg
          ig=connec(ipg,iel)
          ! coordonnées du pg
          xg=node(1,ig)
          yg=node(2,ig)
          ! on calcule la charge au point de Gauss
          vf=source(xg,yg)
          ! calcul de la jacobienne au pg
          ! on pourra le faire directement avec map plus tard
          jac=0.0_f64
          ! boucle sur les fonctions de base d'interpolation
          ! pour construire la jacobienne de la transformation
          ! géométrique
          do iib=1,nloc
             ! gradient de la fonction de base au pg 
             ! dans les variables de référence
             ! calcul du poids de gauss
             call gradpg(iib,ipg,dxy,poids)
             ib=connec(iib,iel)
             ! contribution à la jacobienne
             jac(1,1)=jac(1,1)+node(1,ib)*dxy(1)
             jac(1,2)=jac(1,2)+node(1,ib)*dxy(2)
             jac(2,1)=jac(2,1)+node(2,ib)*dxy(1)
             jac(2,2)=jac(2,2)+node(2,ib)*dxy(2)
          end do
          !write(*,*) ipg,jac
          ! déterminant et transposée de l'inverse
          det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
          cojac(1,1)=jac(2,2)/det
          cojac(2,2)=jac(1,1)/det
          cojac(1,2)=-jac(2,1)/det
          cojac(2,1)=-jac(1,2)/det

          ! assemblage du second membre
          !rho(ig)=rho(ig)+vf*det*poids

          ! nouvelle boucle sur les noeuds locaux de l'élément iel
          ! pour calculer la matrice élémentaire
          do ii=1,nloc
             ! gradients dans les variables de référence
             call gradpg(ii,ipg,gradref_i,poids)  
             ! gradient dans les variables physiques
             grad_i=matmul(cojac,gradref_i)
             ! gradient 
             ! indice global
             i=connec(ii,iel)
             do jj=1,nloc
                call gradpg(jj,ipg,gradref_j,poids)
                j=connec(jj,iel)
                ! gradient dans les variables physiques
                grad_j=matmul(cojac,gradref_j)
                ! produit scalaire gradi . gradj
                v=(grad_i(1)*grad_j(1)+grad_i(2)*grad_j(2))*det*poids
                call add(v,i,j)
             end do
          end do
       end do
    end do

    write(*,*) 'Assemblage conditions aux limites matrice...'

    ! assemblage des conditions aux limites de dirichlet
    ! par pénalisation
    do ii=1,nbc
       i=indexbc(ii)
       call add(big,i,i)
    end do


       
  end subroutine assemb
    
  ! assemblage du second membre
  subroutine assemb_rhs(dg_rho)
    implicit none
    ! matrice locale
    sll_real64 :: jac(2,2),det
    sll_real64 :: dxy(2),poids,vf
    sll_real64, dimension(order+1,order+1,nx,ny) :: dg_rho
    integer :: iel,ipg,i,ii,jj,ig,ib,iib,ielx,iely
    
    ! assemblage de la matrice de rigidité
    ! et du second membre

    write(*,*) 'Assemblage second membre...'

    rho=0.0_f64

    ! boucle sur les éléments
    do iel=1,nel
       ielx=mod(iel-1,nx)+1
       iely=(iel-1)/nx+1
       ! boucle sur les pg pour intégrer
       do ipg=1,nloc
          ! numéro global du pg
          ig=connec(ipg,iel)

          ! on récupère la charge au point de Gauss
          ii=mod(ipg-1,order+1)+1
          jj=(ipg-1)/(order+1)+1

          vf=dg_rho(ii,jj,ielx,iely)
          !write(*,*) 'vf=',vf

          ! calcul de la jacobienne au pg
          ! on pourra le faire directement avec map plus tard
          jac=0.0_f64
          ! boucle sur les fonctions de base d'interpolation
          ! pour construire la jacobienne de la transformation
          ! géométrique
          do iib=1,nloc
             ! gradient de la fonction de base au pg 
             ! dans les variables de référence
             ! calcul du poids de gauss
             call gradpg(iib,ipg,dxy,poids)
             ib=connec(iib,iel)
             ! contribution à la jacobienne
             jac(1,1)=jac(1,1)+node(1,ib)*dxy(1)
             jac(1,2)=jac(1,2)+node(1,ib)*dxy(2)
             jac(2,1)=jac(2,1)+node(2,ib)*dxy(1)
             jac(2,2)=jac(2,2)+node(2,ib)*dxy(2)
          end do
          !write(*,*) ipg,jac
          ! déterminant et transposée de l'inverse
          det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
!!$          cojac(1,1)=jac(2,2)/det
!!$          cojac(2,2)=jac(1,1)/det
!!$          cojac(1,2)=-jac(2,1)/det
!!$          cojac(2,1)=-jac(1,2)/det

          ! assemblage du second membre
          rho(ig)=rho(ig)+vf*det*poids
       end do
    end do
       
    write(*,*) 'Assemblage conditions aux limites rhs...'

    ! assemblage des conditions aux limites de dirichlet
    ! par pénalisation
    do ii=1,nbc
       i=indexbc(ii)
       rho(i)=potexact(node(1,i),node(2,i))*big
       !rho(i)=0
    end do




  end subroutine assemb_rhs




  ! compute (in place) the LU decomposition
  subroutine computeLU()
    implicit none

    integer :: nsym=1,mp=6,ifac=1,isol=0,ier
    sll_real64 :: energ,vu,vfg

    write(*,*) 'Factorisation LU...'

    call sol(vsup,vdiag,vinf,   &
         vfg,kld,vu,neq,mp,ifac,isol, &
         nsym,energ,ier,nsky)  
  end subroutine computeLU


  ! résout le système linéaire
  subroutine compute_phi()
    implicit none

    integer :: nsym=1,mp=6,ifac=0,isol=1,ier
    sll_real64 :: energ !,solution(neq),rhs(neq)

    write(*,*) 'Résolution...'
    
    call sol(vsup,vdiag,vinf,   &
         rho,kld,phi,neq,mp,ifac,isol, &
         nsym,energ,ier,nsky)  

  end subroutine compute_phi


  ! ajoute la valeur val  à la position (i,j)
  ! dans la matrice ligne de ciel
  subroutine add(val,i,j)
    implicit none
    integer,intent(in) :: i,j
    integer :: ll
    sll_real64, intent(in) :: val
 
    if (i-j.gt.prof(i).or.j-i.gt.prof(j)) then
       write(*,*) '(',i,',',j,') out of matrix profile'
       write(*,*) prof(i)
       write(*,*) prof(j)
       stop
    end if
    if (i.eq.j) then
       vdiag(i)=vdiag(i)+val
    else if (j.gt.i) then
       ll=kld(j+1)-j+i
       vsup(ll)=vsup(ll)+val
    else 
       ll=kld(i+1)-i+j
       vinf(ll)=vinf(ll)+val
    end if


  end subroutine add

  ! libère la mémoire
  subroutine release()
    implicit none

    write(*,*) 'Libération mémoire...'

    deallocate(xpg)
    deallocate(wpg)
    deallocate(dlag)
    deallocate(node)
    deallocate(connec)
    deallocate(indexbc)
    deallocate(prof)
    deallocate(kld)
    deallocate(vdiag)
    deallocate(vsup)
    deallocate(vinf)
    deallocate(phi)
    deallocate(rho)

  end subroutine release



end module sll_m_lobalap

