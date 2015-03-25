!> This fortran module is dedicated to advection on
!> triangular mesh.
module sll_advection_2d_tri_mesh
#include "sll_working_precision.h"
#include "sll_constants.h"
#include "sll_memory.h"
#include "sll_utilities.h"
#include "sll_assert.h"
#include "sll_boundary_condition_descriptors.h"

use sll_meshes_base
use sll_triangular_meshes

implicit none


!> @brief 2d advection on triangular mesh
type :: sll_advection_tri_mesh

  type(sll_triangular_mesh_2d), pointer :: mesh  !< mesh
  sll_int32,  dimension(:,:),   pointer :: nvoiv
  sll_int32,  dimension(:),     pointer :: nlpa
  sll_int32,  dimension(:),     pointer :: nbpama
  sll_int32,  dimension(:),     pointer :: iad1
  sll_int32,  dimension(:),     pointer :: indice
  sll_int32,  dimension(:),     pointer :: itabor
  sll_int32,  dimension(:),     pointer :: itest
  sll_int32,  dimension(:),     pointer :: nlmloc
  sll_real64, dimension(:),     pointer :: xbas
  sll_real64, dimension(:,:),   pointer :: xlm
  sll_real64, dimension(:,:),   pointer :: coef
  sll_real64, dimension(:),     pointer :: xp
  sll_real64, dimension(:),     pointer :: yp
  sll_real64, dimension(:),     pointer :: f_out
  logical,    dimension(:),     pointer :: inzone
  sll_int32, dimension(:),      pointer :: numres

end type sll_advection_tri_mesh

contains

!> @brief allocates the memory space for a new 2D advection
!> on triangular mesh on the heap,
!> initializes it with the given triangular mesh and returns a pointer to the
!> object.
!> @param mesh triangular mesh
!> @return a pointer to the newly allocated object.
function new_advection_2d_tri_mesh( mesh ) result(adv)

  type(sll_triangular_mesh_2d), intent(in), target  :: mesh
  type(sll_advection_tri_mesh),             pointer :: adv

  sll_int32  :: ierr
  sll_int32  :: i
  sll_int32  :: j
  sll_int32  :: is1, is2, is3

  SLL_ALLOCATE(adv, ierr)

  adv%mesh => mesh

    
  ! --- Remplissage de "nvoiv" -----------------------------------
  !     Identique a "nvois" sauf pour les aretes appartenant a 
  !     une frontiere. Le chiffre correspond ici a un code pour 
  !     le traitement des conditions aux limites sur les 
  !     particules, alors que dans "nvois" ce chiffre est 
  !     l'oppose du numero de reference de la frontiere concernee
  
  allocate(adv%nvoiv(3,mesh%num_triangles))
  do i = 1,mesh%num_triangles
    do j = 1, 3
      if (mesh%nvois(j,i)>0) then
        adv%nvoiv(j,i) = mesh%nvois(j,i)
      else
        adv%nvoiv(j,i) =  0
      end if
    end do    
  end do
  
  allocate(adv%nlpa(mesh%num_nodes))
  do i = 1,mesh%num_nodes
    adv%nlpa(i) = mesh%npoel2(mesh%npoel1(i)+1)
  end do

  SLL_ALLOCATE(adv%itest(mesh%num_nodes),ierr)
  SLL_ALLOCATE(adv%coef(4,mesh%num_nodes),ierr)
  SLL_ALLOCATE(adv%xlm(3,mesh%num_nodes),ierr)
  SLL_ALLOCATE(adv%nlmloc(mesh%num_nodes),ierr)
  SLL_ALLOCATE(adv%xp(mesh%num_nodes),ierr)
  SLL_ALLOCATE(adv%yp(mesh%num_nodes),ierr)
  SLL_ALLOCATE(adv%inzone(mesh%num_nodes),ierr)
  allocate(adv%numres(mesh%num_nodes))
  allocate(adv%iad1(mesh%num_triangles))
  allocate(adv%indice(mesh%num_triangles))
  allocate(adv%itabor(mesh%num_nodes))
  allocate(adv%f_out(mesh%num_nodes))
  allocate(adv%nbpama(mesh%num_triangles))
  allocate(adv%xbas(mesh%num_nodes))
  
  do i = 1, mesh%num_triangles
  
     is1 = mesh%nodes(1,i) 
     is2 = mesh%nodes(2,i) 
     is3 = mesh%nodes(3,i) 
  
     adv%xbas(is1) = adv%xbas(is1) + mesh%aire(i)/3.
     adv%xbas(is2) = adv%xbas(is2) + mesh%aire(i)/3.
     adv%xbas(is3) = adv%xbas(is3) + mesh%aire(i)/3.
  end do


end function new_advection_2d_tri_mesh


!> @brief 
!> Compute characterisitic origin in triangular mesh
!> @details
!!    xlm1   - 1ere coordonnee barycentrique                   
!!    xlm2   - 2eme coordonnee barycentrique                   
!!    xlm3   - 3eme coordonnee barycentrique                   
!!                                                              
!!    coord  - coordonnees des noeuds                          
!!    nodes  - numero des sommets des triangles                
!!    nvois  - numero des voisins des elements                 
!!    aire   - aire de chaque element                          
!!                                                              
!!    coef1  - tableau temporaire des determinants             
!!    coef2  - tableau temporaire des determinants             
!!    coef3  - tableau temporaire des determinants             
!!    coef4  - tableau temporaire                              
!!                                                              
!!    numpt  - tableau auxiliaire contenant les numeros des    
!!              particules a eliminer                           
!!    nelet  - tableau auxiliaire contenant les numeros des   
!!              elements qui contenaient ces particules         
!!    nlmloc - tableau auxiliaire contenant les numeros des    
!!              elements ou l'on cherche les particules         
!!    numres - tableau auxiliaire contenant les numeros des    
!!              particules non encore localisees                
!!                                                              
!!    itest  - tableau auxiliaire pour preciser le comportement:            
!!    - si  itest=0 la particules reste dans son triangle       
!!    - si  itest=1,2,3 la particule traverse le cote 1,2ou3   
!!                       mais ne traverse pas de frontiere      
!!    - si  itest=11,12,13 la particule est absorbee par le     
!!                         cote frontiere 1,2,3                
!!    - si  itest=21,22,23 la particule est reflechie par       
!!                         le cote frontiere 1,2,3             
!!  Dans tous les cas le chiffre des unites de itest designe    
!!  le numero du cote traverse par la particule.                
!!  Le type de frontiere est defini dans le tableau nvoiv       
!!    - si nvoiv(i,n) > 0  le cote i n'est pas une frontiere    
!!    - si nvoiv(i,n) = 0  le cote i absorbe les particules     
!!    - si nvoiv(i,n) =-1  le cote i reflechit les particules  
!!                                                              
!!    nbpert - nombre de particules a eliminer                 
!!    rho    - densite de charge                        
!!    ad1    - tableau temporaire (adresse de la 1ere    
!!             particule de chaque maille dans le tableau 
!!             ordonne des particules)                     
!!    indice - tableau temporaire (incrementation du nombre 
!!             de particules dja reperees)               
!!    itabor - tableau temporaire (numeros des particules 
!!             ordonnes suivant les numeros des mailles)
!!    nbpama - tableau temporaire (nombre de particules  
!!             par maille)                                 
!!                                                   
!!    petitl - petite longueur de reference             
!<
subroutine advection_2d(adv, f_in, ex, ey, dt )

type(sll_advection_tri_mesh), intent(inout) :: adv  !< mesh
sll_real64, dimension(:),     intent(inout) :: f_in !< distribution function on nodes
sll_real64, dimension(:),     intent(in)    :: ex   !< electric field on x1
sll_real64, dimension(:),     intent(in)    :: ey   !< electric field on x2
sll_real64,                   intent(in)    :: dt   !< time step

sll_real64 :: eps
sll_int32  :: nbpert
sll_int32  :: nbpres
sll_int32  :: num
sll_int32  :: nrest
sll_int32  :: nfin
sll_int32  :: nbpr

sll_int32 :: ip
sll_int32 :: jp

sll_real64 :: pa1x, pa1y, pa2x, pa2y, pa3x, pa3y

sll_real64 :: phi1, phi2, phi3
sll_real64 :: xm11, xm12, xm13
sll_int32  :: nprest
sll_int32  :: mpa, inum, ip1, ip2, ks, ind
sll_int32  :: it, is1, is2, is3
sll_real64 :: x1, x2, x3, y1, y2, y3, det

eps = -adv%mesh%petitl**2

!Set the position of the characteristic origin
!Cell is arbitrary set in one of the cell close to the ip vertex.
do ip = 1, adv%mesh%num_nodes
  adv%numres(ip) = ip
  adv%xp(ip)     = adv%mesh%coord(1,ip) + ex(ip) * dt
  adv%yp(ip)     = adv%mesh%coord(2,ip) + ey(ip) * dt
  adv%nlmloc(ip) = adv%nlpa(ip)
end do
     
nbpres = adv%mesh%num_nodes
nbpert = 0
num    = 0

adv%inzone = .false.

do while( nbpres > 0 )

  !*** Premiere boucle dans le meme element

  nfin  = 0
  nrest = 0
  num   = num + 1

  if ( num == 100 ) then

    write(*,"(//2x,'Arret dans POSITIONS:')",advance='no')
    write(*,"('on n''arrive pas a positionner ')",advance='no')
    write(*,"(i4)",advance='no')nbpres
    write(*,"('  particules')",advance='no')
    write(*,"(/2x,'Particules non positionnees :')",advance='no')
    write(*,"(5x,'numero - position - vitesse - element'/)") 
    do ip=1,min(50,nbpres)
       write(*,"(i10,2x,4e12.4,i10)") adv%numres(ip),      &
                                      adv%xp(ip), adv%yp(ip),  &
                                      ex(ip), ey(ip), ip
    end do
    write(*,"(/5x,a)") & 
    'Il y a certainement un probleme avec le solveur de champ'
    write(*,"(/5x,'(en general la CFL n''est pas verifiee)'/)")

    stop

  end if

  do ip = 1, nbpres    

    jp = adv%numres(ip)
    it = adv%nlmloc(ip)

    x1 = adv%mesh%coord(1,adv%mesh%nodes(1,it))
    y1 = adv%mesh%coord(2,adv%mesh%nodes(1,it))
    x2 = adv%mesh%coord(1,adv%mesh%nodes(2,it))
    y2 = adv%mesh%coord(2,adv%mesh%nodes(2,it))
    x3 = adv%mesh%coord(1,adv%mesh%nodes(3,it))
    y3 = adv%mesh%coord(2,adv%mesh%nodes(3,it))

    pa1x = x1 - adv%xp(jp)
    pa1y = y1 - adv%yp(jp)
    pa2x = x2 - adv%xp(jp)
    pa2y = y2 - adv%yp(jp)
    pa3x = x3 - adv%xp(jp)
    pa3y = y3 - adv%yp(jp)

    adv%coef(1,ip) = pa1x*pa2y - pa1y*pa2x
    adv%coef(2,ip) = pa2x*pa3y - pa2y*pa3x
    adv%coef(3,ip) = pa3x*pa1y - pa3y*pa1x

    if(      adv%coef(1,ip) >= eps    &
       .and. adv%coef(2,ip) >= eps    &
       .and. adv%coef(3,ip) >= eps ) then

       nfin = nfin + 1

       det = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)

       adv%xlm(1,jp) = adv%coef(1,ip) / det
       adv%xlm(2,jp) = adv%coef(2,ip) / det
       adv%xlm(3,jp) = adv%coef(3,ip) / det

       adv%nlpa(jp)   = adv%nlmloc(ip)
       adv%itest(ip)  = 0
       adv%inzone(jp) = .true.

    end if

  end do

  !*** Deuxieme boucle pour celles qui sont sorties

  nbpr = nbpres - nfin

  if( nbpr .ne. 0 ) then

    do ip = 1, nbpres

      jp = adv%numres(ip)
      it = adv%nlmloc(ip)

      if (       adv%coef(1,ip) <  eps   &
           .and. adv%coef(2,ip) >= eps   &
           .and. adv%coef(3,ip) >= eps   ) then

         adv%itest(ip) = 1 + 10*(1-min(1,adv%nvoiv(1,it)))

      end if

      if (       adv%coef(1,ip) >= eps   &
           .and. adv%coef(2,ip) <  eps   &
           .and. adv%coef(3,ip) >= eps   ) then
 
         adv%itest(ip) = 2 + 10*(1-min(1,adv%nvoiv(2,it)))

      end if
 
      if (       adv%coef(1,ip) >= eps   &
           .and. adv%coef(2,ip) >= eps   &
           .and. adv%coef(3,ip) <  eps   ) then
 
         adv%itest(ip) = 3 + 10*(1-min(1,adv%nvoiv(3,it)))

      end if
    
      if (   adv%coef(1,ip) < eps    &
       .and. adv%coef(2,ip) < eps )  then

         pa2x = adv%mesh%coord(1,adv%mesh%nodes(2,it))-adv%xp(jp)
         pa2y = adv%mesh%coord(2,adv%mesh%nodes(2,it))-adv%yp(jp)
     
         adv%coef(4,ip) = pa2x*ey(jp) - pa2y*ex(jp)

         adv%itest(ip) = 1 + max(0,nint(sign(1d0,adv%coef(4,ip))))
         adv%itest(ip) = adv%itest(ip)+10*(1-min(1,adv%nvoiv(adv%itest(ip),it)))

      end if

      if (       adv%coef(2,ip) < eps     &
           .and. adv%coef(3,ip) < eps )  then

         pa3x = adv%mesh%coord(1,adv%mesh%nodes(3,it))-adv%xp(jp) 
         pa3y = adv%mesh%coord(2,adv%mesh%nodes(3,it))-adv%yp(jp)
 
         adv%coef(4,ip) = pa3x*ey(jp) - pa3y*ex(jp)

         adv%itest(ip) = 2 + max(0,nint(sign(1d0,adv%coef(4,ip))))
         adv%itest(ip) = adv%itest(ip)+10*(1-min(1,adv%nvoiv(adv%itest(ip),it)))

      end if

      if (       adv%coef(3,ip) < eps    &
           .and. adv%coef(1,ip) < eps )  then

         pa1x = adv%mesh%coord(1,adv%mesh%nodes(1,it))-adv%xp(jp) 
         pa1y = adv%mesh%coord(2,adv%mesh%nodes(1,it))-adv%yp(jp)

         adv%coef(4,ip) = pa1x*ey(jp) - pa1y*ex(jp)

         adv%itest(ip) = 1 +mod(2+max(0,nint(sign(1d0,adv%coef(4,ip)))),3)
         adv%itest(ip) = adv%itest(ip)+10*(1-min(1,adv%nvoiv(adv%itest(ip),it)))

      end if

    end do

    !*** Particules absorbees par une frontiere -------------------
    
    do ip=1,nbpres
  
      if( adv%itest(ip) > 10 .and. adv%itest(ip) < 14 )  then
         nbpert = nbpert + 1
      end if
  
    end do
  
    !*** Reorganisation des tableaux temporaires ------------------
    !    Particules traversant un cote interne ou reflechie
  
    do ip=1,nbpres
  
      if( (adv%itest(ip) >  0 .and. adv%itest(ip) <  4) .or.     &
          (adv%itest(ip) > 20 .and. adv%itest(ip) < 24)) then
        nrest = nrest + 1
        adv%numres(nrest) = adv%numres(ip)
        
        if (adv%itest(ip)> 20) then
          adv%nlmloc(nrest) = adv%mesh%nvois(adv%itest(ip)-20,adv%nlmloc(ip))  &
                          * max(0,sign(1,10-adv%itest(ip)))        &
                          + adv%nlmloc(ip)* max(0,sign(1,adv%itest(ip)-20))    
        else
          adv%nlmloc(nrest) = adv%mesh%nvois(adv%itest(ip),adv%nlmloc(ip)) &
                          * max(0,sign(1,10-adv%itest(ip)))    &
                          + adv%nlmloc(ip)* max(0,sign(1,adv%itest(ip)-20))    
        end if

      end if
  
    end do

  end if

  nbpres = nrest

end do

!Recherche du nombre de particules de chaque maille -------

adv%nbpama = 0
do ip = 1 , adv%mesh%num_nodes
  if (adv%inzone(ip)) then
    mpa             = adv%nlpa(ip)
    adv%nbpama(mpa) = adv%nbpama(mpa) + 1
  end if
end do

!--- Determination de l'adresse de la premiere particule 
!    de chaque maille dans le tableau ordonne


adv%iad1 = 0
ks = 1
do it = 1 , adv%mesh%num_triangles
   if ( adv%nbpama(it) .ne. 0 )  then
     adv%iad1(it) = ks
     ks           = ks + adv%nbpama(it)
   end if
end do

!--- Construction du tableau ordonne des particules -----------


adv%indice = 0
adv%itabor = 0

do ip = 1, adv%mesh%num_nodes

   if (adv%inzone(ip)) then
     mpa             = adv%nlpa(ip)
     ind             = adv%iad1(mpa) + adv%indice(mpa)
     adv%itabor(ind) = ip
     adv%indice(mpa) = adv%indice(mpa) + 1
   end if

end do

nprest = adv%mesh%num_nodes

adv%f_out = 0.0_f64

do it = 1 , adv%mesh%num_triangles
      
   xm11 = 0.; xm12 = 0.; xm13 = 0.

   !nbpama(it)  !Nombre de particules dans la maille numero it
            
   if (adv%nbpama(it) .ne. 0) then

      ip1 = adv%iad1(it)
      ip2 = ip1 + adv%nbpama(it) - 1

      !Boucle sur les particules situees dans la maille courante

      do ip = ip1 , ip2 
            
         inum = adv%itabor(ip)
       
         phi1 = adv%xlm(2,inum) * f_in(inum) 
         phi2 = adv%xlm(3,inum) * f_in(inum)
         phi3 = adv%xlm(1,inum) * f_in(inum)
               
         xm11 = xm11 + phi1
         xm12 = xm12 + phi2
         xm13 = xm13 + phi3

      end do

      is1 = adv%mesh%nodes(1,it)
      is2 = adv%mesh%nodes(2,it)
      is3 = adv%mesh%nodes(3,it)

      adv%f_out(is1) = adv%f_out(is1) + xm11 
      adv%f_out(is2) = adv%f_out(is2) + xm12 
      adv%f_out(is3) = adv%f_out(is3) + xm13 
   
      if (adv%nbpama(it) == nprest ) exit

      nprest = nprest - adv%nbpama(it)

   end if

end do

f_in = adv%f_out

end subroutine advection_2d

end module sll_advection_2d_tri_mesh
