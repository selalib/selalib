!File: Program Main
!*Programme principal*
!   
!   PICSOU 
!
!   Pierre NAVARO
!
!   Institut de Recherche en Mathematique Avancee
!
!   IRMA   UMR 7501 CNRS/ULP
!
!   7 rue René Descartes F-67084 Strasbourg Cedex, FRANCE.
!
!   http://www-irma.u-strasbg.fr/~navaro
!
!   email: navaro@math.u-strasbg.fr
!
!   tel: 03 90 24 01 73
!
!Use:
! - <Module Zone>
! - <Module Maillage>
! - <Module Sorties>
! - <Module Champs_Externes>
! - <Module Solveurs_Module>
! - <Module Injections>
! - <Module Source_Module>
! - <Module Particules>
! - <Module Projections>
! - <Module Interpolations>
! - <Module Save_Module>
! - <Module Mpi_Module>
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program PICSOU
                        
use zone,            only: readin, lmodte, iout, mesh_fields, objet_fonction
use solveurs_module, only: lecture_donnees_solveur
use poisson,         only: poissn, poifrc, init_solveur_poisson, poliss

use sll_triangular_meshes
use sll_mesh_calculus_2d_module
use sll_gnuplot

!----------------------------------------------------------------------

implicit none

integer :: i
integer :: n, lu, nsolve
integer :: istep

real(8) :: tcpu
real(8), dimension(:), allocatable :: rho, phi

character(len=72) :: argv, dirpr

type(mesh_fields)  :: mxw
type(sll_triangular_mesh_2d)    :: mesh
character(len=132):: inpfil
character(len=132):: maafil
logical :: lask

integer :: ntypfr(5)
real(8) :: potfr(5)

n = 1 !iargc()
do i = 1, n
   call getarg( i, argv); write(*,'(i2, 1x, a)') i, argv
end do

!------------------------------------------------------------!
!     Reads in default parameters from input file (.inp)        !
!     If LASK=t, ask user if file is to be read.                !
!------------------------------------------------------------!

lask = .true.; lu = 10

inpfil = trim(argv)//'.inp' 
maafil = trim(argv)//'.maa'

write(*,"(/10x,'Fichier d''entree : ',a)") inpfil

open(lu,file=inpfil,status='OLD',err=80)
write(*,1050,advance='no') trim(inpfil)
lask = .false.

80 continue

if (lask) then
   write(*,1900) inpfil
   write(*,1800,advance='no')
   read(*,"(a)") inpfil
   write(*,1700) trim(inpfil)
end if

call readin(trim(argv), dirpr, nsolve)           !Donnees generales

call read_from_file(mesh, maafil)                    !Lecture maillage
call write_triangular_mesh_mtv(mesh, "picsou.mtv")
call lecture_donnees_solveur(inpfil,ntypfr, potfr)  

call analyze_triangular_mesh(mesh, ntypfr) 

allocate(mxw%e(3,mesh%num_nodes)); mxw%e=0.0; 
allocate(mxw%b(3,mesh%num_nodes)); mxw%b=0.0; 
allocate(mxw%j(3,mesh%num_nodes)); mxw%j=0.0; 
allocate(rho(mesh%num_nodes)); rho = 0.0
allocate(phi(mesh%num_nodes)); phi = 0.0

call cpu_time(tcpu)  !Initialisation du temps CPU

call init_solveur_poisson(mesh, ntypfr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Equation de POISSON - elements finis !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
istep = 1

call poissn(potfr, mxw, rho, phi, mesh, istep)
call poliss(phi, mxw, mesh)
call poifrc(mxw, mesh, ntypfr)

call sll_gnuplot_2d( phi, "phi", mesh%coord, mesh%nodes, 1)

1050 format(/' Read settings from file  ', A, ' ?  Y')
1700 format(/' New parameters write to file  ', A, /)
1800 format(/' Settings may have been changed - New title :')
1900 format(/' Input file  ', A,'  not found')

end program PICSOU
