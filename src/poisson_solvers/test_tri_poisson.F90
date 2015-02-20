!   Pierre NAVARO
!
!   Institut de Recherche en Mathematique Avancee
!
!   IRMA   UMR 7501 CNRS/Unistra
!
!   7 rue René Descartes F-67084 Strasbourg Cedex, FRANCE.
!
!   http://www-irma.u-strasbg.fr/~navaro
!
!   email: navaro@math.cnrs.fr
!
!   tel: 03 68 85 01 73
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program test_tri_poisson
                        
use tri_poisson
use sll_triangular_meshes
use sll_mesh_calculus_2d_module
use sll_gnuplot

!----------------------------------------------------------------------

implicit none

real(8) :: tcpu
real(8), dimension(:), allocatable :: rho, phi
real(8), dimension(:), allocatable :: ex, ey

character(len=72) :: argv

type(sll_triangular_poisson_2d)       :: solver
type(sll_triangular_mesh_2d), pointer :: mesh
character(len=132):: inpfil
character(len=132):: maafil
logical :: lask

call getarg( 1, argv); write(*,'(1x, a)') argv

!------------------------------------------------------------!
!     Reads in default parameters from input file (.inp)        !
!     If LASK=t, ask user if file is to be read.                !
!------------------------------------------------------------!

lask = .true.

inpfil = trim(argv)//'.inp' 
maafil = trim(argv)//'.maa'

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

mesh => new_triangular_mesh_2d(maafil) 
call write_triangular_mesh_mtv(mesh, "picsou.mtv")

call analyze_triangular_mesh(mesh) 

allocate(ex(mesh%num_nodes)); ex=0.0; 
allocate(ey(mesh%num_nodes)); ey=0.0; 
allocate(rho(mesh%num_nodes)); rho = 0.0
allocate(phi(mesh%num_nodes)); phi = 0.0

call cpu_time(tcpu)  !Initialisation du temps CPU

call sll_create(solver, mesh, inpfil)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Equation de POISSON - elements finis !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call poissn(solver, ex, ey, rho, phi)
call poliss(solver, phi, ex, ey)
call poifrc(solver, ex, ey)

call sll_gnuplot_2d( phi, "phi", mesh%coord, mesh%nodes, 1)

1050 format(/' Read settings from file  ', A, ' ?  Y')
1700 format(/' New parameters write to file  ', A, /)
1800 format(/' Settings may have been changed - New title :')
1900 format(/' Input file  ', A,'  not found')

end program test_tri_poisson
