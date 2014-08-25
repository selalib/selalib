module sll_generate_tri_mesh

implicit none

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

logical, public :: PERIODIC_BC = .true.

CONTAINS

subroutine plaqy( x_min, x_max, nbox, &
                  y_min, y_max, nboy, &
                  coor, ntri, nbs, nbt)

!*-----------------------------------------------------------------------
!*  Generation du maillage pour une plaque rectangulaire dans le plan xOy.
!*  Maillage regulier de type TRIANGLES P1 - Q4T.
!*-----------------------------------------------------------------------

real(8), dimension(:,:), allocatable :: coor
integer, dimension(:,:), allocatable :: ntri
integer :: nbs     ! nb de sommets
integer :: nbt     ! nb de triangles
integer :: nd1, nd2, nd3
integer :: nelt, nlp, nhp
integer :: l, l1, ll1, ll2
integer :: in, iin, n, noeud, nsom, nxm, nym, neltot, ndd
integer :: nbox, nboy, iel, iev

real(8) :: xx1, xx2, pasx0, pasx1, pasy0, pasy1
real(8) :: alx, aly, y
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
allocate(refs(nbs)); refs = 0
allocate(nvois(3,nbt)) ; nvois = imxref !Numeros des voisins
allocate(nvoif(3,nbt))

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

!*** Recherche des elements sur la frontiere
nefro = 4*(nbox-2) + 4*(nboy-2)
nelfr = 2*(nbox-1) + 2*(nboy-1) - 2
nelin = nbt - nelfr 

allocate(nar(2*(nbox+nboy)))

nmxsd = 1

!---- Description des frontieres ----

nhp = nbox + 1
do i = 1 , nbox
   nar(i) = nhp - i 
   refs(nar(i)) = 2
end do        

nhp = nbox * nym
do i = 1 , nbox
   nar(i) = nhp + i 
   refs(nar(i)) = 3
end do        

nlp = nboy + 1
do i = 1 , nboy
   nar(i) = nbox * (nlp - i)
   y = coor(2,nar(i))
   if ( abs(y) > 0.01875) then
      refs(nar(i)) = 1
   else 
      refs(nar(i)) = 5
   end if
end do        

do i = 1 , nboy
   nar(i) = i * nbox - nxm
   refs(nar(i)) = 4
end do        

nmxfr = 0
!*** Initialisation des voisins ***

if (PERIODIC_BC) then

   do i = 1, nym
      iel = (i-1) * 2 * nym + 1
      iev =  iel + 2 * nym - 1 
      nvois(1,iel) = iev
      nvois(3,iev) = iel
      nvoif(1,iel) = 3
      nvoif(3,iev) = 1
   end do
   
   do i = 1, nxm
      iel = (i-1) * 2 + 1
      iev =  iel + 2 * nxm * (nym-1) + 1 
      nvois(3,iel) = iev
      nvois(2,iev) = iel
      nvoif(3,iel) = 2
      nvoif(2,iev) = 3
   end do

   nelfr = 0
   nefro = 0

else

   !*** Voisins en bas et en haut
   
   do i = 1, 2*(nbox-1), 2
      nvois(3,i      ) = -2
      nvois(2,nbt-i+1) = -3
   end do
   
   !*** Voisins de chaque cote lateral
   
   nvois(1,1) = -4
   nhp = 2*(nbox-1)
   do i = 1, nboy - 2
      j = i * nhp
      nvois(3,j  ) = -1
      nvois(1,j+1) = -4 
   end do
   nvois(3,nbt) = -1

   do i = 1, nbs
      if (refs(i) >= nmxfr) nmxfr = nmxfr+1
   end do

end if

call write_mesh(coor, ntri, nbs, nbt, nbox , nboy, alx, aly, 'y')

end subroutine plaqy

subroutine plaqx( x_min, x_max, nbox, &
                  y_min, y_max, nboy, &
                  coor, ntri, nbs, nbt)

!*-----------------------------------------------------------------------
!*  Generation du maillage pour une plaque rectangulaire dans le plan xOy.
!*  Maillage regulier de type TRIANGLES P1 - Q4T.
!*-----------------------------------------------------------------------

real(8), dimension(:,:), allocatable :: coor
integer, dimension(:,:), allocatable :: ntri
integer :: nbs     ! nb de sommets
integer :: nbt     ! nb de triangles
integer :: nd1, nd2, nd3
integer :: nelt, nlp, nhp, iel
integer :: l, l1, ll1, ll2
integer :: in, iin, n, noeud, nsom, nxm, nym, neltot, ndd
integer :: nbox, nboy, n1, n2, cent1, cent2

real(8) :: xx1, xx2, pasx0, pasx1, pasy0, pasy1
real(8) :: alx, aly
real(8) :: x_min, x_max, y_min, y_max

integer :: iev

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
allocate(refs(nbs)); refs = 0
allocate(nvois(3,nbt)) ; nvois = imxref !Numeros des voisins
allocate(nvoif(3,nbt)) ; nvoif = imxref !Numeros des voisins)
nbt = nbt
nbs = nbs

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

!*** Recherche des elements sur la frontiere
nefro = 6*(nbox-3) + 6*(nboy-3) + 16
nelfr = 2*(nbox-1) + 2*(nboy-1) 
nelin = nbt - nelfr 
write(*,*)" Nbre d'elements ayant au moins un noeud sur la frontiere ", nefro
write(*,*)" Nbre d'elements internes", nelin
write(*,*)" Nbre d'elements situes sur la frontiere ", nelfr

allocate(nar(2*(nbox+nboy)))

!---- Description des frontieres ----

!*** Frontiere Sud ***
nhp = nbox + 1 
do i = 1 , nbox
   nar(i) = nhp - i 
   refs(nar(i)) = 1
end do        

!*** Frontiere Nord ***
nhp = nbox * nym
do i = 1 , nbox
   nar(i) = nhp + i 
   refs(nar(i)) = 3
end do        

!*** Frontiere Est ***

do i = 1 , nboy
   nar(i) = i * nbox - nxm
   refs(nar(i)) = 2
end do        
refs(nar(1)) = 1
refs(nar(nboy)) = 2

!*** Frontiere Ouest ***

nlp = nboy + 1 
do i = 1 , nboy
   nar(i) = nbox * (nlp - i)
   refs(nar(i)) = 4
end do        
refs(nar(1)) = 4
refs(nar(nboy)) = 1


if(PERIODIC_BC) then


   do i = 1, nym
      iel = 1 + nxm*4*(i-1)
      iev = iel+nxm*4-2
      nvois(1,iel) = iev
      nvois(1,iev) = iel
      nvoif(1,iel) = 1
      nvoif(1,iev) = 1
      print*, iel, iev
   end do
   
   do i = 1, nxm
      iel = i * 4
      iev = iel + 4 * (nxm-1) + (nym-2)*4*nxm+2
      nvois(1,iel) = iev
      nvois(1,iev) = iel
      nvoif(1,iel) = 1
      nvoif(1,iev) = 1
      print*, iel, iev
   end do
   stop

else

   do iel = 4, 4*(nbox-1), 4
      n1 = refs(ntri(1,iel))
      n2 = refs(ntri(2,iel))
      nvois(1,iel) = - max(n1,n2)
   end do

   do i = 4, 4*(nbox-1), 4
      nvois(1,nbt-i+2) = -3
   end do
   
   nhp = 4*(nbox-1)
   do i = 1, nboy - 1
      nvois(1,(i-1)*nhp+1) = - 2
   end do
   
   do i = 1, nboy-1
      iel = i* nhp -1
      n1 = refs(ntri(1,iel))
      n2 = refs(ntri(2,iel))
      nvois(1,iel) = - max(n1,n2)
   end do

end if

nmxfr = 0
do i = 1, nbs
   if (refs(i) >= nmxfr) nmxfr = nmxfr+1
end do

call write_mesh(coor, ntri, nbs, nbt, nbox , nboy, alx, aly, 'x')

end subroutine plaqx

subroutine write_mesh(coor, ntri, nbs, nbt, nbox, nboy, alx, aly, choix)
real(8), dimension(:,:), allocatable :: coor
integer, dimension(:,:), allocatable :: ntri
integer :: nbs     ! nb de sommets
integer :: nbt     ! nb de triangles
character(len=*) :: choix
real(8) :: x1, y1, alx, aly
integer :: nbox, nboy, iel, nhp, nlp, nxm, nym

nxm    = nbox - 1           !nombre de quadrangles sur l'axe Ox du maillage
nym    = nboy - 1           !nombre de quadrangles sur l'axe Oy du maillage

!Trace du maillage au format plotmtv
open( 10, file="plaque.mtv")

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

!Reference des noeuds

write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='References des noeuds' "

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
   write(10,"(i4)"  , advance="no") refs(i)
   write(10,"(a)")"'"
end do

write(10,*)

write(6,"(//'Nb de noeuds                : ',i10/   &
        &   'Nb de triangles             : ',i10/   &
        &   'Nb max de front referencees : ',i10/   &
        &   'Nb max de SD references     : ',i10/   &
        &  /'Nb de triangles ayant au moins 1 sommet'&
        & ,' sur une frontiere : ',i10/)")      &
        &    nbs,nbt,nmxfr,nmxsd,nefro

write(6,"(//'Nb d''elements internes     : ',i10/   &
        &   'Nb d''elements frontieres   : ',i10/)") nelin,nelfr


if (choix == 'x') then
   !---- Description des frontieres ----
   write(10,*)"$DATA=CURVE3D"
   write(10,*)"%equalscale=T"
   write(10,*)"%toplabel='Description des frontieres' "
   write(10,*)"%xmin=",-0.1*alx, " xmax = ", 1.1*alx
   write(10,*)"%ymin=",-0.1*aly, " ymax = ", 1.1*aly
   
   !*** Noeuds du bas de gauche a droite
   nhp = nbox + 1 
   do i = 1 , nbox
      nar(i) = nhp - i 
      x1 = coor(1,nar(i))
      y1 = coor(2,nar(i))
      write(10,"(a)"   ,advance="no")"@text x1="
      write(10,"(g15.3)",advance="no") x1
      write(10,"(a)"   ,advance="no")" y1="
      write(10,"(g15.3)",advance="no") y1
      write(10,"(a)"   ,advance="no")" z1=0. ll="
      write(10,"(i5)") refs(nar(i))
   end do        
   
   !*** Noeuds du haut de droite a gauche
   nhp = nbox * nym
   do i = 1 , nbox
      nar(i) = nhp + i 
      x1 = coor(1,nar(i))
      y1 = coor(2,nar(i))
      write(10,"(a)"   ,advance="no")"@text x1="
      write(10,"(g15.3)",advance="no") x1
      write(10,"(a)"   ,advance="no")" y1="
      write(10,"(g15.3)",advance="no") y1
      write(10,"(a)"   ,advance="no")" z1=0. ll="
      write(10,"(i5)") refs(nar(i))
   end do        
   
   !*** Noeuds a droite de bas en haut
   do i = 1 , nboy
      nar(i) = i * nbox - nxm
      x1 = coor(1,nar(i))
      y1 = coor(2,nar(i))
      write(10,"(a)"   ,advance="no")"@text x1="
      write(10,"(g15.3)",advance="no") x1
      write(10,"(a)"   ,advance="no")" y1="
      write(10,"(g15.3)",advance="no") y1
      write(10,"(a)"   ,advance="no")" z1=0. ll="
      write(10,"(i5)") refs(nar(i))
   end do        
   
   !*** Noeuds a gauche de haut en bas
   nlp = nboy + 1 
   do i = 1 , nboy
      nar(i) = nbox * (nlp - i)
      x1 = coor(1,nar(i))
      y1 = coor(2,nar(i))
      write(10,"(a)"   ,advance="no")"@text x1="
      write(10,"(g15.3)",advance="no") x1
      write(10,"(a)"   ,advance="no")" y1="
      write(10,"(g15.3)",advance="no") y1
      write(10,"(a)"   ,advance="no")" z1=0. ll="
      write(10,"(i5)") refs(nar(i))
      write(10,*)
   end do        
   
   !*** Voisins en bas de droite a gauche
   
   do i = 4, 4*(nbox-1), 4
      write(10,*)"%linecolor=",-nvois(1,i)
      write(10,*)coor(1,ntri(1,i)),coor(2,ntri(1,i)),0.
      write(10,*)coor(1,ntri(2,i)),coor(2,ntri(2,i)),0.
      write(10,*)
      x1 = 0.5*(coor(1,ntri(1,i))+coor(1,ntri(2,i)))
      y1 = 0.5*(coor(2,ntri(1,i))+coor(2,ntri(2,i)))
      write(10,"(a)"   ,advance="no")"@text x1="
      write(10,"(g15.3)",advance="no") x1
      write(10,"(a)"   ,advance="no")" y1="
      write(10,"(g15.3)",advance="no") y1
      write(10,"(a)"   ,advance="no")" z1=0. ll="
      write(10,"(i4)"  ,advance="no") -nvois(1,i)
      write(10,"(a)") " "
      write(10,*)
   end do
   
   !*** Voisins en haut de gauche a droite
   
   do i = 4, 4*(nbox-1), 4
      iel = nbt-i+2
      write(10,*)"%linecolor=",-nvois(1,iel)
      write(10,*)coor(1,ntri(1,iel)),coor(2,ntri(1,iel)),0.
      write(10,*)coor(1,ntri(2,iel)),coor(2,ntri(2,iel)),0.
      write(10,*)
      x1 = 0.5*(coor(1,ntri(1,iel))+coor(1,ntri(2,iel)))
      y1 = 0.5*(coor(2,ntri(1,iel))+coor(2,ntri(2,iel)))
      write(10,"(a)"   ,advance="no")"@text x1="
      write(10,"(g15.3)",advance="no") x1
      write(10,"(a)"   ,advance="no")" y1="
      write(10,"(g15.3)",advance="no") y1
      write(10,"(a)"   ,advance="no")" z1=0. ll="
      write(10,"(i4)"  ,advance="no") -nvois(1,iel)
      write(10,"(a)") " "
      write(10,*)
   end do
   
   !*** Voisins a droite de bas en haut ***
   
   nhp = 4*(nbox-1)
   do i = 1, nboy - 1
      iel = (i-1)*nhp+1
      write(10,*)"%linecolor=",-nvois(1,iel)
      write(10,*)coor(1,ntri(1,iel)),coor(2,ntri(1,iel)),0.
      write(10,*)coor(1,ntri(2,iel)),coor(2,ntri(2,iel)),0.
      write(10,*)
      x1 = 0.5*(coor(1,ntri(1,iel))+coor(1,ntri(2,iel)))
      y1 = 0.5*(coor(2,ntri(1,iel))+coor(2,ntri(2,iel)))
      write(10,"(a)"   ,advance="no")"@text x1="
      write(10,"(g15.3)",advance="no") x1
      write(10,"(a)"   ,advance="no")" y1="
      write(10,"(g15.3)",advance="no") y1
      write(10,"(a)"   ,advance="no")" z1=0. ll="
      write(10,"(i4)"  ,advance="no") -nvois(1,iel)
      write(10,"(a)") " "
      write(10,*)
   end do
   
   !*** Voisins a gauche de bas en haut ***
   
   do i = 1, nboy - 1
      iel = i* nhp -1
      write(10,*)"%linecolor=",-nvois(1,iel)
      write(10,*)coor(1,ntri(1,iel)),coor(2,ntri(1,iel)),0.
      write(10,*)coor(1,ntri(2,iel)),coor(2,ntri(2,iel)),0.
      write(10,*)
      x1 = 0.5*(coor(1,ntri(1,iel))+coor(1,ntri(2,iel)))
      y1 = 0.5*(coor(2,ntri(1,iel))+coor(2,ntri(2,iel)))
      write(10,"(a)"   ,advance="no")"@text x1="
      write(10,"(g15.3)",advance="no") x1
      write(10,"(a)"   ,advance="no")" y1="
      write(10,"(g15.3)",advance="no") y1
      write(10,"(a)"   ,advance="no")" z1=0. ll="
      write(10,"(i5)") -nvois(1,iel)
      write(10,*)
   end do
end if

write(10,*)"$END"
close(10)
   
end subroutine write_mesh


subroutine calmai(ntri, nbs, nbt)

integer, dimension(:,:), allocatable :: ntri
integer :: nbs     ! nb de sommets
integer :: nbt     ! nb de triangles
integer, dimension (:), allocatable :: indc

integer :: it, iel
integer :: i1, i2, i3, is, is1, is2, is3
integer :: jel1, jel2, jel3, nel1, nel2, nel3
logical :: ldefr, ltemp, lcalsu, lmainv
integer :: itab(3,3)

!-----
!---- itab(isomloc1,isomloc2) est le numero local de l'arete qui a pour
!---- sommets locaux isomloc1 et isomloc2
!-----

itab = reshape(source=(/0,1,3,1,0,2,3,2,0/), shape=(/3,3/))

ldefr  = .false.
ltemp  = .true.
lcalsu = .false.
lmainv = .false.

write(iout,"(/////10x,'>>> Calcul des quantites liees au maillage <<<'/)")

! --- Gestion des triangles ayant un noeud en commun -----------
write(iout,*)"*** Gestion des triangles ayant un noeud en commun ***"
 
! ... recherche des elements ayant un sommet commun
!     creation du tableau npoel1(i+1)  contenant le nombre de 
!     triangles ayant le noeud i en commun

allocate(npoel1(nbs+1))

npoel1 = 0
do i=1,nbt
   is1 = ntri(1,i)
   is2 = ntri(2,i)
   is3 = ntri(3,i)
   npoel1(is1+1) = npoel1(is1+1)+1
   npoel1(is2+1) = npoel1(is2+1)+1
   npoel1(is3+1) = npoel1(is3+1)+1
end do

! ... le tableau npoel1 devient le tableau donnant l'adresse 
!     dans npoel2 du dernier element dans la suite des triangles
!     communs a un noeud

npoel1(1)=0
do i=3,nbs+1
   npoel1(i)=npoel1(i-1)+npoel1(i)
end do

! ... creation du tableau npoel2 contenant sequentiellement les 
!     numeros des triangles ayant un noeud en commun      
!     le premier triangle s'appuyant sur le noeud i est
!     adresse par "npoel1(i)+1" 
!     le nombre de triangles ayant le noeud i en commun est
!     "npoel1(i+1)-npoel1(i)"


allocate(npoel2(npoel1(nbs+1)))
allocate(indc(nbs))

indc   = 1  !Le tableau temporaire indc doit etre initialise a 1

do it = 1,nbt
   do k = 1,3
      is = ntri(k,it)
      npoel2(npoel1(is)+indc(is)) = it
      indc(is) = indc(is)+1
   end do
end do
 
! --- Recherche des numeros des triangles voisins d'un triangle 

write(iout,*)"*** Recherche des numeros des triangles voisins d'un triangle ***"

do iel=1,nbt

   ! ... numeros des 3 sommets du triangle

   is1=ntri(1,iel)
   is2=ntri(2,iel)
   is3=ntri(3,iel)

   ! ... boucles imbriquees sur les elements pointant vers
   !     les 2 noeuds extremites de l'arete consideree
   !     Le voisin est le triangle commun (hormis iel)

   ! ... premiere arete (entre le sommet is1 et is2)

   nel1=npoel1(is1+1)-npoel1(is1) !nb de triangles communs a is1
   nel2=npoel1(is2+1)-npoel1(is2) !nb de triangles communs a is2

   loop1:do i1=1,nel1
            jel1=npoel2(npoel1(is1)+i1) !premier triangle is1
            if(jel1.ne.iel) then
               do i2=1,nel2
                  jel2=npoel2(npoel1(is2)+i2)
                  if(jel2 == jel1) then
                     nvois(1,iel)  = jel1
                     exit loop1
              end if
           end do
            end if
     end do loop1

   ! ... deuxieme arete (entre le sommet is2 et is3)

   nel2=npoel1(is2+1)-npoel1(is2)
   nel3=npoel1(is3+1)-npoel1(is3)

   loop2:do i2=1,nel2
            jel2=npoel2(npoel1(is2)+i2)
            if(jel2 /= iel) then
               do i3=1,nel3
                  jel3=npoel2(npoel1(is3)+i3)
                  if(jel3 == jel2) then
                     nvois(2,iel)=jel2
             exit loop2
              end if
           end do
            end if
     end do loop2

   ! ... troisieme arete (entre le sommet is3 et is1)

   nel3=npoel1(is3+1)-npoel1(is3)
   nel1=npoel1(is1+1)-npoel1(is1)

   loop3:do i3=1,nel3
            jel3=npoel2(npoel1(is3)+i3)
            if(jel3 /= iel) then
               do i1=1,nel1
                  jel1=npoel2(npoel1(is1)+i1)
                  if(jel1 == jel3) then
                     nvois(3,iel)=jel3
             exit loop3
              end if
           end do
            end if
     end do loop3
end do

!Calcul de nvoif : Numero local dans le voisin de la face j commune a l'elt i

nvoif(1:3,:) = nvois(1:3,:)

do i = 1, nbt
   do j = 1,3
      jel1 = nvois(j,i)
      if (jel1 > 0) then
         do k = 1,3
            jel2 = nvois(k,jel1)
            if (jel2 == i) then 
               nvoif(j,i) = k
               exit
            end if
         end do
      end if
   end do
end do

end subroutine calmai

end module sll_generate_tri_mesh
