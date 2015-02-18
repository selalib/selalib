!File: Module Sorties
!Module contenant les subroutines de diagnostics
module sorties

use zone, only: mesh_fields, lmodte, lmodtm, nstep,     &
        time, titre,        &
        iout,  dt,        &
        pi, c, xmu0, petitx, grandx, &
        objet_fonction,mcentr,minter,mexter

use maillage, only: mesh_data, voronoi, nctfrp, nctfro

use gnuplot_module, only: gnu_output, champs_gnuplot
use utlib, only: utunit
use solveurs_module, only: nfrmx, nspmx

implicit none

integer, private :: i, j, k, isf, iz
integer :: iplot1 = 0, iplot2 = 0

!* ================================================================= *
!*                                                    *
!*   Caracteristiques des images                                     *
!*                                                                   *
!*                                                                   *
logical :: ldiach ! Activateur general de diag  des champs           *
logical :: gnuplot! Activateur general de diag  GNUplot           *
!*                                                    *
!*   Champs                                            *
!*                                                      *
integer :: itrcha ! Iteration pour le debut des traces               *
integer :: jtrcha ! Frequence des traces                             *
logical :: ltrb1  ! graphe des isovaleurs de B1                      *
logical :: ltrb2  ! graphe des isovaleurs de B2                      *
logical :: ltrb3  ! graphe des isovaleurs de B1                      *
logical :: ltre1  ! graphe des isovaleurs de E1                      *
logical :: ltre2  ! graphe des isovaleurs de E2                      *
logical :: ltre3  ! graphe des isovaleurs de E3                      *
logical :: ltrj1  ! graphe des isovaleurs de J1                      *
logical :: ltrj2  ! graphe des isovaleurs de J2                      *
logical :: ltrj3  ! graphe des isovaleurs de J3                      *
logical :: ltrpo  ! graphe des isovaleurs de V                       *
logical :: ltrrh  ! graphe des isovaleurs de rho                     *
!*                                                    *
!*                                                    *
!* ================================================================= *

contains

!Function: donnees_diag
!                                                                   
!But:  
!Lecture des donnees 'diagnostics champs et particules' 
!          namelist "nldiag", pour le controle                     
!                - des traces d'isovaleurs et de champs de vecteurs
!                - des positions des particules                   
!                - des espaces de phase                          
!                - des sorties de fichiers sur listing        
!                - des ecritures sur disques pour le post-processeur
subroutine donnees_diag(inpfil)

logical :: lerr
integer :: ispl
logical :: ltrc1,ltrc2,ltrbt,ltret,ltrct,ltrev,ltrcv
integer :: nbrfen
logical :: ltrmai
logical :: ltreph(4,4,4), ltrpfn(4,4), lstock,lstchb,lwrite
integer :: istcha, jstcha
character(len=16) :: cdesti = 'sel ldummy;ex'
character(len=*), intent(in) :: inpfil

namelist/nldiag/ ldiach,itrcha,jtrcha               &
                ,ltrb1,ltrb2,ltrb3,ltre1,ltre2,ltre3            &
                ,ltrj1,ltrj2,ltrj3,ltrpo,ltrrh,           &
                ltrc1,ltrc2,ltrbt,ltret,ltrct,ltrev,ltrcv,  &
                nbrfen,ltrmai,   &
                gnuplot

!Valeurs pour la compatibilite avec les anciens fichiers M2V (degas2d)  
ltrc1=.false.;ltrc2=.false.
ltrbt=.false.;ltret=.false.;ltrct=.false.;ltrev=.false.;ltrcv=.false.
nbrfen=0;ltrmai=.false.
ltreph=.false.;ltrpfn=.false.;lstock=.false.
istcha=0;jstcha=0;lstchb=.false.;lwrite=.false.
! --- 1.0 --- Valeurs par defaut ---------------------------------------

ldiach = .false.

itrcha = 0; jtrcha = 1
ltrb1  = .false.; ltrb2  = .false.; ltrb3  = .false.
ltre1  = .false.; ltre2  = .false.; ltre3  = .false.
ltrj1  = .false.; ltrj2  = .false.; ltrj3  = .false.
ltrpo = .false.; ltrrh = .false.

gnuplot = .false.

open(10, file=inpfil, status="old")
write(*,"(//10x,'Lecture de la namelist $nldiag',/)")
read(10, nldiag,err=100)
goto 200
100 continue
write(*,*)"Erreur de lecture de la namelist $nldiag"
write(*,*)"Valeurs prises par defaut"
write(*,nldiag)
stop
200 close(10)

if (ltrc1) ltrj1  = .true.
if (ltrc2) ltrj2  = .true.

! --- 3.0 --- Ecriture des donnees -------------------------------------
   
write(iout,808) ldiach
write(iout,8091) gnuplot
write(iout,810) itrcha,jtrcha
write(iout,815) ltrb1,ltrb2,ltrb3
write(iout,816) ltre1,ltre2,ltre3
write(iout,817) ltrj1,ltrj2,ltrj3
write(iout,818) ltrrh,ltrpo

lerr =.false.
if (jtrcha <= 0) lerr=.true.
if (lerr) then 
   write(iout,855)
   goto 800
end if 

! --- 8.0 --- Point de sortie-------------------------------------------
 
      goto 850
 800  continue
      call errout(6,"F","Erreur lecture donnees","module diagnostics")
 850  continue

! --- 9.0 --- Formats --------------------------------------------------
 
 808  format( /10x,'      TRACES EN LIGNE'              &
     &       //10x,'Active les diag des champs         ldiach : ',L10   )
 8091 format( /10x,'Active les sorties GNUPLOT         gnuplot: ',L10   )
 810  format( /10x,'Iteration pour le debut des traces itrcha : ',I10   &
     &        /10x,'Frequence des traces               jtrcha : ',I10) 
 815  format( 10x,'Isovaleurs B            ltrb1,ltrb2,ltrb3 : ',3L4)
 816  format( 10x,'Isovaleurs E            ltre1,ltre2,ltre3 : ',3L4)
 817  format( 10x,'Isovaleurs J            ltrj1,ltrj2,ltrj3 : ',3L4)
 818  format( 10x,'Isovaleurs rho et V         ltrrh,ltrpo : ',2L4)
 855  format( 10x,'Erreur flagrante dans les donnees'/)
 

end subroutine donnees_diag

!************************************************************************

!Function: diagcha
!Ecriture des fichiers de sortie pour les champs electromagnetiques
!
!Parametres:
!istep  - Numero de l'iteration
!eb     - Champs electromagnetiques
!cr     - Courants
!rho    - Densite de Charge
!phi    - Potentiel
!mesh   - Maillage
!dirpr  - Repertoire de sortie
subroutine diagcha( istep, ebj, rho, phi, mesh, dirpr )

  integer, save :: iplot

  integer,        intent(in) :: istep
  type(mesh_fields),  intent(in) :: ebj
  type(mesh_data),    intent(in) :: mesh
  double precision, dimension(:),     intent(in) :: phi, rho
  character(len=*),   intent(in) :: dirpr

  !Ecriture d'un fichier dx pour la visu du maillage + particule

  if (ldiach) then

     if (istep >= itrcha .and. mod(istep, jtrcha) == 0.) then

        iplot = iplot + 1

        if (gnuplot) then
           call champs_gnuplot(ebj, rho, phi, mesh, iplot, trim(dirpr),   &
                ltre1, ltre2, ltre3, ltrb1, ltrb2, ltrb3,  &
                ltrj1, ltrj2, ltrj3, ltrpo, ltrrh)
        end if

     end if

  end if

end subroutine diagcha

end module sorties
