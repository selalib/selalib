!File: Module Solveurs_Module
!Definition des solveurs de Maxwell et de Poisson
module solveurs_module
#include "sll_utilities.h"

!----------------------------------------------------------------------
implicit none


!Variable: nfrmx 
!nombre maximum de frontieres referencees
!Variable: ntypfr 
!Variables: Conditions aux limites
!            potfr  - potentiel applique sur les frontieres Dirichlet  
!            1 - DIRICHLET 
!            3 - NEUMAN    

contains

!Subroutine: lecture_donnees_solveur
!  Lecture des donnees relatives aux champs                 
!   Namelist "nlcham"                                        
!
!Traduction en Fortran 90 de leccha
!
!  33-LECCHA  -   Version 1.0  Novembre 1991   A. Adolf  
subroutine lecture_donnees_solveur(nomfich, ntypfr, potfr)

Character(len=*), intent(in) :: nomfich
integer                      :: ityp
integer                      :: ifr
integer                      :: i
integer                      :: nfrmx
integer, dimension(:)        :: ntypfr
double precision             :: potfr(:)

NAMELIST/nlcham/ntypfr,potfr

! ----------- Valeurs par defaut et initialisations -------------------t
write(6,900)

nfrmx  = size(ntypfr)
ntypfr = 1  
potfr  = 0.

!--- 2.0 --- Lecture des donnees --------------------------------------

open(10, file = nomfich)
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
    SLL_ERROR(" ")
  end if
end do

900  format(//10x,'Conditions aux limites sur les frontieres')
902  format(/10x,'Reference :',i3,5x,'Dirichlet ')
904  format(/10x,'Reference :',i3,5x,'Neumann')
910  format(/10x,'Option non disponible'/  &
    &        10x,'Reference :',i3,5x,'Type :',i3/)
921  format(//5x,'Frontiere',5x,'ntypfr',7x,'potfr')
922  format(5x,I3,4x,I10,2E12.3,2I10,E12.3)
932  format(//10x,'CONDITIONS AUX LIMITES'/)

end subroutine lecture_donnees_solveur

end module solveurs_module
