!!File: gnuplot_module.f90
module gnuplot_module

use zone, only: lmodte, lmodtm, mesh_fields, iout,  &
        lcorrp, titre, time
use maillage, only: sll_triangular_mesh_2d, nsomare, ngnoeelt, nnoeselt

implicit none

integer, private :: i, j, k

contains

subroutine champs_gnuplot(ebj, rho, phi, mesh, iplot, dir,  &
              ltre1, ltre2, ltre3, ltrb1, ltrb2, ltrb3, &
              ltrj1, ltrj2, ltrj3, ltrpo, ltrrh)

type(mesh_fields) :: ebj
double precision, dimension(:) :: phi, rho
type(sll_triangular_mesh_2d) :: mesh
integer :: iplot
character(len=*) :: dir
logical, intent(in) :: ltre1, ltre2, ltre3, ltrb1, ltrb2, ltrb3
logical, intent(in) :: ltrj1, ltrj2, ltrj3
logical, intent(in) :: ltrpo, ltrrh

integer :: kk0, kk1, kk2, kk3, kk4
character(len=4) :: fin
character(len=1) :: aa,bb,cc,dd
character(len=40):: nom

kk0 = iplot
kk1 = kk0/1000
aa  = char(kk1 + 48)
kk2 = (kk0 - kk1*1000)/100
bb  = char(kk2 + 48)
kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
cc  = char(kk3 + 48)
kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
dd  = char(kk4 + 48)
fin = aa//bb//cc//dd
   
nom = trim(dir)//trim(titre)//fin

if (lmodte) then

   if (ltre1) call chloop(mesh,trim(nom)//'ex',ebj%e(1,:),trim(titre)//"_plotex.gnu",iplot)
   if (ltre2) call chloop(mesh,trim(nom)//'ey',ebj%e(2,:),trim(titre)//"_plotey.gnu",iplot)
   if (ltrb3) call chloop(mesh,trim(nom)//'bz',ebj%b(3,:),trim(titre)//"_plotbz.gnu",iplot)
   !if (ltrb1) call chloop(mesh,trim(nom)//'bx',ebj%b(1,:),trim(titre)//"_plotbx.gnu",iplot)
   !if (ltrb2) call chloop(mesh,trim(nom)//'by',ebj%b(2,:),trim(titre)//"_plotby.gnu",iplot)
   !if (ltre3) call chloop(mesh,trim(nom)//'ez',ebj%e(3,:),trim(titre)//"_plotez.gnu",iplot)
   if (ltrj1) call chloop(mesh,trim(nom)//'jx',ebj%j(1,:),trim(titre)//"_plotjx.gnu",iplot)
   if (ltrj2) call chloop(mesh,trim(nom)//'jy',ebj%j(2,:),trim(titre)//"_plotjy.gnu",iplot)

end if

if (lmodtm) then

      if (ltrb1) call chloop(mesh,trim(nom)//'bx',ebj%b(1,:),trim(titre)//"_plotbx.gnu",iplot)
      if (ltrb2) call chloop(mesh,trim(nom)//'by',ebj%b(2,:),trim(titre)//"_plotby.gnu",iplot)
      if (ltre3) call chloop(mesh,trim(nom)//'ez',ebj%e(3,:),trim(titre)//"_plotez.gnu",iplot)
      if (ltrj3) call chloop(mesh,trim(nom)//'jz',ebj%j(3,:),trim(titre)//"_plotjz.gnu",iplot)

end if

if (ltrrh) call chloop(mesh,trim(nom)//'r0',rho,trim(titre)//"_plotrho.gnu",iplot)
if (ltrpo) call chloop(mesh,trim(nom)//'fi',phi,trim(titre)//"_plotphi.gnu",iplot)

end subroutine champs_gnuplot 

!************************************************************************

subroutine chloop(mesh, nomchamp, champ, f_gnu, iplot)

type(sll_triangular_mesh_2d) :: mesh
character(len=*) :: nomchamp, f_gnu
double precision, dimension(:) :: champ
double precision :: xs1, xs2, xs3, ys1, ys2, ys3
integer :: iplot

write(*,"(10x, 'Fichier de sortie GNUplot ', a)") nomchamp

open(83, file = f_gnu, position="append")
if ( iplot == 1 ) rewind(83)
write(83,*)"set title'Field "//f_gnu//"'"
write(83,*)"splot '"//nomchamp//"' title 'Time = ",time,"' w l"
close(83)

open(94, file = nomchamp)

do i = 1, mesh%num_cells

   xs1 = mesh%coord(1,mesh%nodes(1,i)); ys1 = mesh%coord(2,mesh%nodes(1,i))
   xs2 = mesh%coord(1,mesh%nodes(2,i)); ys2 = mesh%coord(2,mesh%nodes(2,i))
   xs3 = mesh%coord(1,mesh%nodes(3,i)); ys3 = mesh%coord(2,mesh%nodes(3,i))

   write(94,*)xs1, ys1, champ(mesh%nodes(1,i))
   write(94,*)xs2, ys2, champ(mesh%nodes(2,i))
   write(94,*)xs3, ys3, champ(mesh%nodes(3,i))
   write(94,*)xs1, ys1, champ(mesh%nodes(1,i))
   write(94,*)
   write(94,*)

end do

close(94)

end subroutine chloop

subroutine gnu_output( ebj, rho, phi, mesh, ltre1, ltre2, ltre3     &
             , ltrb1, ltrb2, ltrb3, ltrj1, ltrj2, ltrj3,    &
               ltrpo, ltrrh )

type (mesh_fields) :: ebj
double precision, dimension(:) :: phi, rho
type (sll_triangular_mesh_2d) :: mesh
logical, intent(in) :: ltre1, ltre2, ltre3, ltrb1, ltrb2, ltrb3
logical, intent(in) :: ltrj1, ltrj2, ltrj3
logical, intent(in) :: ltrpo, ltrrh

if (lmodte) then

   if (ltre1) call gnu_loop(mesh, "ex", ebj%e(1,:))
   if (ltre2) call gnu_loop(mesh, "ey", ebj%e(2,:))
   if (ltrb3) call gnu_loop(mesh, "bz", ebj%b(3,:))
   if (ltrj1) call gnu_loop(mesh, "jx", ebj%j(1,:))
   if (ltrj2) call gnu_loop(mesh, "jy", ebj%j(2,:))

end if

if (lmodtm) then

   if (ltrb1) call gnu_loop(mesh, "bx", ebj%b(1,:))
   if (ltrb2) call gnu_loop(mesh, "by", ebj%b(2,:))
   if (ltre3) call gnu_loop(mesh, "ez", ebj%b(3,:))
   if (ltrj3) call gnu_loop(mesh, "jz", ebj%j(3,:))

end if

if (ltrrh) call gnu_loop(mesh, "rho", rho)
if (ltrpo) call gnu_loop(mesh, "phi", phi)

end subroutine gnu_output

subroutine gnu_loop(mesh, nomchamp, champ)

type(sll_triangular_mesh_2d) :: mesh
character(len=*) :: nomchamp
double precision, dimension(:) :: champ
integer :: is1, is2, is3

write(*,"(/10x,a)")     &
& " Sortie pour GNUplot de "//nomchamp//" dans fichier "//trim(titre)//nomchamp//".dat"

open(27, file = trim(titre)//nomchamp//".dat")

do i = 1, mesh%num_cells

    is1 = mesh%nodes(1,i)
    is2 = mesh%nodes(2,i)
    is3 = mesh%nodes(3,i)

    write(27,*)sngl(mesh%coord(1,is1)),sngl(mesh%coord(2,is1)),sngl(champ(is1))
    write(27,*)sngl(mesh%coord(1,is2)),sngl(mesh%coord(2,is2)),sngl(champ(is2))
    write(27,*)sngl(mesh%coord(1,is3)),sngl(mesh%coord(2,is3)),sngl(champ(is3))
    write(27,*)sngl(mesh%coord(1,is1)),sngl(mesh%coord(2,is1)),sngl(champ(is1))
    write(27,*)
    write(27,*)

end do

close(27)
   
end subroutine gnu_loop

!Function: vigie
!Ecriture de fichier pour visualisation avec VIGIE (INRIA)
!
!Parametres:
! eb    - Champs electromagnetiques
! cr    - Courants
! mesh  - Maillage
! rho   - Densite de charge
! phi   - Potentiel
! iplot - Numero du trace
!
subroutine vigie( nom, mesh, phi, iplot )

character(len=*) :: nom
type(sll_triangular_mesh_2d), intent(in) :: mesh
double precision, dimension(:), intent(in) :: phi
integer, intent(in) :: iplot

integer :: iunit0 = 18, iunit1 = 19

real, dimension (:,:), allocatable :: q, xx
integer, dimension (6) :: nnperelement(6)
integer :: l, nb, nblock, nu, nnu, nnperel, ielemtype, ihmg, itype

integer :: kk0, kk1, kk2, kk3, kk4
character(len=04):: fin
character(len=01):: aa,bb,cc,dd

if ( minval(phi(:)) /= maxval(phi(:)) ) then


   kk0 = iplot
   kk1 = kk0/1000
   aa  = char(kk1 + 48)
   kk2 = (kk0 - kk1*1000)/100
   bb  = char(kk2 + 48)
   kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
   cc  = char(kk3 + 48)
   kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
   dd  = char(kk4 + 48)
   fin = aa//bb//cc//dd

   write(*,"(/5x,a,/)")" --> sortie "//fin//" vigie  pour "//nom

   !  unstructured homogeneous mesh
   iunit1 = 1
   nblock = 1      !  total number of blocks
   itype  = 0      !  itype = 0: Unstructured mesh
               !  itype = 1: Structured mesh
   ihmg   = 0      !  ihmg  = 0: homogeneous mesh
               !  ihmg  = 1: non-homogeneous mesh
   ielemtype = 3       !  ielemtype: type of surface elements
               !            3: T3, 3-nodes triangular cell
               !            4: Q4, 4-nodes quadrilateral cell
   nnperelement(3) = 3     !  nnperel: number of nodes per element.
   nnperelement(4) = 4
   nnperel = nnperelement(ielemtype)

   !  nbs: number of nodes
   !  nbt: number of surface elements
   !  mesh%coor: nodes coordinates
   !  mesh%nodes: connectivity table
   !  q: unknowns array

   nnu  = 1            !  nnu: number of unknowns stored
   allocate(q(mesh%num_nodes,nnu))
   allocate(xx(2,mesh%num_nodes))
   do i = 1, mesh%num_nodes
      xx(1,i) = sngl(mesh%coord(1,i))
      xx(2,i) = sngl(mesh%coord(2,i))
      q(i,1)  = sngl(phi(i))
   end do

   open(iunit1, file=trim(titre)//'_'//nom//"_vigie"//fin//".dat", form="unformatted")
   write (iunit1) nblock
   do nb = 1,nblock
      write(iunit1) itype, nnu
      write(iunit1) ihmg, mesh%num_nodes
      write(iunit1) mesh%num_cells, ielemtype
      write(iunit1) ((xx(l,i),l=1,2),i=1,mesh%num_nodes)
      write(iunit1) ((q(i,nu),i=1,mesh%num_nodes),nu=1,nnu)
      write(iunit1) ((mesh%nodes(l,i),l=1,nnperel),i=1,mesh%num_cells)
   end do

   !write description file

   open(iunit0, file=trim(titre)//'_'//nom//"_vigie"//fin//".desc" )
   write(iunit0,*) "hfep"
   write(iunit0,*) "./"//trim(titre)//'_'//nom//"_vigie"//fin//".dat"
   write(iunit0,"(a)") nom
   close(iunit0); close(iunit1)
   
   write(iout,"(10x,a,/)")" Sortie pour VIGIE fichier: vigie"//fin

   !Sortie Vigie - Fichier ASCII

   open(iunit1,file=trim(titre)//'_'//nom//'_mesh.ascii2d_'//fin//'.dat')

   write(iunit1,*) 'points',mesh%num_nodes   !ecriture des noeuds
   do i = 1, mesh%num_nodes
      write(iunit1,*) sngl(mesh%coord(1:2,i))
   end do

   write(iunit1,*) 'triangles', mesh%num_cells   !ecriture des triangles
   do i = 1 , mesh%num_cells
      write(iunit1,"(3i6)") mesh%nodes(1:3,i)
   enddo

   write(iunit1,*) 'scalars '//nom     !ecriture des variables
   do i = 1 , mesh%num_nodes
      write(iunit1,*) sngl(phi(i))
   enddo
   close(iunit1)

   open(iunit0, file=trim(titre)//'_'//nom//"_vigie_ascii2d_"//fin//'.desc' )
   rewind(iunit0)
   write(iunit0,"(a)")'ascii2d'
   write(iunit0,"(a)")trim(titre)//'_'//nom//'_mesh.ascii2d_'//fin//'.dat'
   write(iunit0,"(a)") nom
   close(iunit0)
   deallocate(q)
   deallocate(xx)

else

    write(*,"(/5x,a,/)")" pas de sortie vigie  pour "//nom//" (champ constant)"

end if

end subroutine vigie

subroutine vigie2( nom, mesh, phi, iplot )

character(len=*) :: nom
type(sll_triangular_mesh_2d), intent(in) :: mesh
double precision, dimension(:), intent(in) :: phi
integer, intent(in) :: iplot

integer :: iunit0 = 18, iunit1 = 19

real, dimension (:,:), allocatable :: xx

integer :: i1, i2
integer :: kk0, kk1, kk2, kk3, kk4
character(len=04):: fin
character(len=01):: aa,bb,cc,dd

if ( minval(phi(:)) /= maxval(phi(:)) ) then

   kk0 = iplot
   kk1 = kk0/1000
   aa  = char(kk1 + 48)
   kk2 = (kk0 - kk1*1000)/100
   bb  = char(kk2 + 48)
   kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
   cc  = char(kk3 + 48)
   kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
   dd  = char(kk4 + 48)
   fin = aa//bb//cc//dd

   iunit1 = 1
   allocate(xx(2,mesh%num_nodes+mesh%nbtcot))
   k = 0
   do i = 1, mesh%num_nodes
      k = k+1
      xx(1,k) = sngl(mesh%coord(1,i))
      xx(2,k) = sngl(mesh%coord(2,i))
   end do
   do i = 1, mesh%nbtcot
      k = k+1
      xx(1,k) = sngl(0.5*(mesh%coord(1,nsomare(1,i))+mesh%coord(1,nsomare(2,i))))
      xx(2,k) = sngl(0.5*(mesh%coord(2,nsomare(1,i))+mesh%coord(2,nsomare(2,i))))
   end do

   !Sortie Vigie - Fichier ASCII

   open(iunit1,file=trim(titre)//'_'//nom//'_mesh.ascii2d_'//fin//'.dat')

   write(iunit1,*) 'points',mesh%num_nodes+mesh%nbtcot   !ecriture des noeuds
   do i = 1, mesh%num_nodes+mesh%nbtcot
      write(iunit1,*) xx(1:2,i)
   end do

   write(iunit1,*) 'triangles', 4*mesh%num_cells   !ecriture des triangles
   do i = 1 , mesh%num_cells
      do j = 1, 4
         write(iunit1,"(3i6)") (ngnoeelt(nnoeselt(k,j),i),k=1,3)
      end do
   enddo

   write(iunit1,*) 'scalars '//nom     !ecriture des variables
   do i = 1 , mesh%num_nodes+mesh%nbtcot
      write(iunit1,*) sngl(phi(i))
   enddo
   close(iunit1)

   open(iunit0, file=trim(titre)//'_'//nom//"_vigie_ascii2d_"//fin//'.desc' )
   rewind(iunit0)
   write(iunit0,"(a)")'ascii2d'
   write(iunit0,"(a)")trim(titre)//'_'//nom//'_mesh.ascii2d_'//fin//'.dat'
   write(iunit0,"(a)") nom
   close(iunit0)
   deallocate(xx)

else

    write(*,*) " pas de sortie vigie  pour "//nom//" (champ constant)"

end if

end subroutine vigie2

end module gnuplot_module
