module sll_generate_tri_mesh

implicit none
private

integer, dimension(:,:), allocatable :: nvois
integer, dimension(:,:), allocatable :: nvoif

integer, dimension(:), allocatable   :: refs
integer, dimension(:), allocatable   :: npoel1
integer, dimension(:), allocatable   :: npoel2

integer :: i, j, k, ifr

integer :: imxref=99999999 ! entier designant la reference par defaut
integer, parameter :: iout = 6

integer :: nmxfr   ! nb total de frontieres referencees
integer :: nmxsd   ! nb de sous domaines references
integer :: nelin   ! nb de triangles internes
integer :: nelfr   ! nb de triangles sur une frontiere
integer :: nefro   ! nb de triangles ayant 1 noeud sur une frontiere

integer, dimension(:), allocatable :: nar

logical, public :: PERIODIC_BC = .false.

public :: plaqx, plaqy

CONTAINS

subroutine plaqy( x_min, x_max, nbox, &
                  y_min, y_max, nboy, &
                  coor, ntri, nbs, nbt)

!*-----------------------------------------------------------------------
!*  Generation du maillage pour une plaque rectangulaire dans le plan xOy.
!*  Maillage regulier de type TRIANGLES P1 - Q4T.
!*-----------------------------------------------------------------------

real(8), dimension(:,:), pointer :: coor
integer, dimension(:,:), pointer :: ntri
integer :: nbs     ! nb de sommets
integer :: nbt     ! nb de triangles
integer :: nd1, nd2, nd3
integer :: nelt
integer :: l, l1, ll1, ll2
integer :: in, iin, n, noeud, nsom, nxm, nym, neltot, ndd
integer :: nbox, nboy

real(8) :: xx1, xx2, pasx0, pasx1, pasy0, pasy1
real(8) :: alx, aly
real(8) :: x_min, x_max, y_min, y_max

ndd = max0(nbox,nboy)
alx = x_max-x_min
aly = y_max-y_min

!*---- Determination de quelques utilitaires

nxm    = nbox - 1      !nombre de quadrangles sur l'axe Ox du maillage
nym    = nboy - 1      !nombre de quadrangles sur l'axe Oy du maillage
neltot = 2 * nxm * nym !nombre total de triangles du maillage 
nsom   = nbox * nboy   !nombre de "sommets" du maillage quadrangulaire
noeud  = nsom          !nombre total de noeuds

nbt = neltot; nbs = noeud
allocate(coor(1:2,nbs))
allocate(ntri(1:3,1:nbt))

!*---- Ecriture des coordonnees des sommets ----*!
                  
pasx0 = alx / nxm                   !pas le long de OX
pasx1 = pasx0 * 0.5                 !demi - pas       
pasy0 = aly / nym                   !pas le long de OY
pasy1 = pasy0 * 0.5                 !demi - pas       

do n = 0 , nym       !Coordonnees des noeuds
   in = nbox * n
   do i = 1 , nbox
      iin = in + i                  !numero du noeud
      xx1 = x_min + alx - (i - 1) * pasx0   !coordonnees du noeud
      xx2 = y_min + n * pasy0 
      coor(1,iin) = xx1
      coor(2,iin) = xx2
   end do
end do

!*---- Table de connectivite ----*!

nelt = 0                        !initialisation
do l = 1 , nym                  !boucle sur les lignes

   l1  = l - 1
   ll1 = l1 * nbox
   ll2 = l  * nbox 

   do i = 1 , nxm               !boucle sur les colonnes

      nd1 = ll1 + i             
      nd2 = ll2 + i             
      nd3 = ll2 + i + 1         

      nelt = nelt + 1           !premier e.f. du "quadrangle"

      ntri(1,nelt) = nd1       !premier  noeud
      ntri(2,nelt) = nd2       !deuxieme noeud
      ntri(3,nelt) = nd1 + 1   !troisieme noeud

      nelt = nelt + 1           !second e.f. du "quadrangle"
       
      ntri(1,nelt) = nd1 + 1       !premier noeud        
      ntri(2,nelt) = nd2       !deuxieme noeud
      ntri(3,nelt) = nd3       !troisieme noeud

   end do

end do

call write_mesh(coor, ntri, nbs, nbt, 'y')

end subroutine plaqy

subroutine plaqx( x_min, x_max, nbox, &
                  y_min, y_max, nboy, &
                  coor, ntri, nbs, nbt)

!*-----------------------------------------------------------------------
!*  Generation du maillage pour une plaque rectangulaire dans le plan xOy.
!*  Maillage regulier de type TRIANGLES P1 - Q4T.
!*-----------------------------------------------------------------------

real(8), dimension(:,:), pointer :: coor
integer, dimension(:,:), pointer :: ntri
integer :: nbs     ! nb de sommets
integer :: nbt     ! nb de triangles
integer :: nd1, nd2, nd3
integer :: nelt
integer :: l, l1, ll1, ll2
integer :: in, iin, n, noeud, nsom, nxm, nym, neltot, ndd
integer :: nbox, nboy, cent1, cent2

real(8) :: xx1, xx2, pasx0, pasx1, pasy0, pasy1
real(8) :: alx, aly
real(8) :: x_min, x_max, y_min, y_max

ndd  = max(nbox,nboy)

!*---- Determination de quelques utilitaires

nxm    = nbox - 1           !nombre de quadrangles sur l'axe Ox du maillage
nym    = nboy - 1           !nombre de quadrangles sur l'axe Oy du maillage
neltot = 4 * nxm * nym      !nombre total de triangles du maillage 
nsom   = nbox * nboy        !nombre de "sommets" du maillage quadrangulaire
noeud  = nsom + nxm * nym       !nombre total de noeuds
nmxsd  = 1

nbt = neltot; nbs = noeud
allocate(coor(1:2,nbs))
allocate(ntri(1:3,1:nbt))

!*---- Ecriture des coordonnees des sommets ----*!

alx = x_max - x_min
aly = y_max - y_min
                  
pasx0 = alx / nxm                   !pas le long de OX
pasx1 = pasx0 * 0.5                 !demi - pas       
pasy0 = aly / nym                   !pas le long de OY
pasy1 = pasy0 * 0.5                 !demi - pas       

!Coordonnees des noeuds

do n = 0 , nym       
   in = nbox * n
   do i = 1 , nbox
      iin = in + i                          !numero du noeud
      xx1 = x_min + alx - (i - 1) * pasx0           !coordonnees du noeud
      xx2 = y_min + n * pasy0 
      coor(1,iin) = xx1
      coor(2,iin) = xx2
   end do
end do

do n = 0 , nym-1
   in = nxm * n
   do i = 1 , nxm
      iin = in + i + nsom                   !numero du noeud
      xx1 = x_min + alx - (i - 1) * pasx0 - pasx1       !coordonnees du noeud
      xx2 = y_min + n * pasy0 + pasy1
      coor(1,iin) = xx1
      coor(2,iin) = xx2
   end do
end do

!*---- Table de connectivite ----*!

nelt = 0                        !initialisation
do l = 1 , nym                  !boucle sur les lignes

   l1  = l - 1
   ll1 = l1 * nbox
   ll2 = l  * nbox 
   cent1 = nsom + (l1*nxm)

   do i = 1 , nxm               !boucle sur les colonnes

      cent2 = cent1 + i

      nd1 = ll1 + i             
      nd2 = ll2 + i             
      nd3 = cent2

      nelt = nelt + 1    

      ntri(1,nelt) = nd1
      ntri(2,nelt) = nd2
      ntri(3,nelt) = nd3

      nelt = nelt + 1    
       
      ntri(1,nelt) = nd2 
      ntri(2,nelt) = nd2+1
      ntri(3,nelt) = cent2

      nelt = nelt + 1    
       
      ntri(1,nelt) = nd2+1 
      ntri(2,nelt) = nd1+1
      ntri(3,nelt) = cent2

      nelt = nelt + 1    
       
      ntri(1,nelt) = nd1+1
      ntri(2,nelt) = nd1
      ntri(3,nelt) = cent2

   end do

end do

call write_mesh(coor, ntri, nbs, nbt, 'x')

end subroutine plaqx

subroutine write_mesh(coor, ntri, nbs, nbt, choix)
real(8), dimension(:,:), pointer :: coor
integer, dimension(:,:), pointer :: ntri
integer :: nbs     ! nb de sommets
integer :: nbt     ! nb de triangles
character(len=*) :: choix
real(8) :: x1, y1

!Trace du maillage au format plotmtv
open( 10, file="plaqu"//choix//".mtv")

!--- Trace du maillage ---

write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Maillage' "

do i = 1, nbt
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)coor(1,ntri(2,i)),coor(2,ntri(2,i)),0.
   write(10,*)coor(1,ntri(3,i)),coor(2,ntri(3,i)),0.
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)
end do

!--- Numeros des noeuds et des triangles

write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Numeros des noeuds et des triangles' "

do i = 1, nbt
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)coor(1,ntri(2,i)),coor(2,ntri(2,i)),0.
   write(10,*)coor(1,ntri(3,i)),coor(2,ntri(3,i)),0.
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)
end do

do i = 1, nbt
   x1 = (  coor(1,ntri(1,i))  &
         + coor(1,ntri(2,i))  &
         + coor(1,ntri(3,i))    )/3.
   y1 = (  coor(2,ntri(1,i))  &
         + coor(2,ntri(2,i))  &
         + coor(2,ntri(3,i))    )/3.
   write(10,"(a)"   , advance="no")"@text x1="
   write(10,"(f8.5)", advance="no") x1
   write(10,"(a)"   , advance="no")" y1="
   write(10,"(f8.5)", advance="no") y1
   write(10,"(a)"   , advance="no")" z1=0. lc=4 ll='"
   write(10,"(i4)"  , advance="no") i
   write(10,"(a)")"'"
end do

do i = 1, nbs
   x1 = coor(1,i)
   y1 = coor(2,i)
   write(10,"(a)"   , advance="no")"@text x1="
   write(10,"(g15.3)", advance="no") x1
   write(10,"(a)"   , advance="no")" y1="
   write(10,"(g15.3)", advance="no") y1
   write(10,"(a)"   , advance="no")" z1=0. lc=5 ll='"
   write(10,"(i4)"  , advance="no") i
   write(10,"(a)")"'"
end do

!--- Numeros des noeuds 

write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Numeros des noeuds' "

do i = 1, nbt
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)coor(1,ntri(2,i)),coor(2,ntri(2,i)),0.
   write(10,*)coor(1,ntri(3,i)),coor(2,ntri(3,i)),0.
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)
end do

do i = 1, nbs
   x1 = coor(1,i)
   y1 = coor(2,i)
   write(10,"(a)"   , advance="no")"@text x1="
   write(10,"(g15.3)", advance="no") x1
   write(10,"(a)"   , advance="no")" y1="
   write(10,"(g15.3)", advance="no") y1
   write(10,"(a)"   , advance="no")" z1=0. lc=5 ll='"
   write(10,"(i4)"  , advance="no") i
   write(10,"(a)")"'"
end do

!--- Numeros des triangles

write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Numeros des triangles' "

do i = 1, nbt
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)coor(1,ntri(2,i)),coor(2,ntri(2,i)),0.
   write(10,*)coor(1,ntri(3,i)),coor(2,ntri(3,i)),0.
   write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
   write(10,*)
end do

do i = 1, nbt
   x1 = (  coor(1,ntri(1,i))  &
         + coor(1,ntri(2,i))  &
         + coor(1,ntri(3,i))    )/3.
   y1 = (  coor(2,ntri(1,i))  &
         + coor(2,ntri(2,i))  &
         + coor(2,ntri(3,i))    )/3.
   write(10,"(a)"   , advance="no")"@text x1="
   write(10,"(g15.3)", advance="no") x1
   write(10,"(a)"   , advance="no")" y1="
   write(10,"(g15.3)", advance="no") y1
   write(10,"(a)"   , advance="no")" z1=0. lc=4 ll='"
   write(10,"(i4)"  , advance="no") i
   write(10,"(a)")"'"
end do

write(10,*)"$END"
close(10)
   
end subroutine write_mesh

end module sll_generate_tri_mesh
