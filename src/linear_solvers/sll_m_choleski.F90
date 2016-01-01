module sll_m_choleski
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "sll_working_precision.h"

  implicit none

  public :: &
    sll_s_choles, &
    sll_s_desrem

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

contains

!
!factoriser 'cholesky' la matrice profil "ae" dans "as"
!
!in:
!
!  mudl,ae      - pointeur,matrice(symetrique definie >0)
!  eps          - test pour les pivots
!
!out:
!
!  as           - matrice factorisee
!  ntest        - = 0 si aucun pivot < eps
!               - = 1 sinon
!
!  programmeur : 
!  f hecht - juin 84 inria
!
subroutine sll_s_choles(mudl,ae,as)

!**********************************************************************

sll_int32,  dimension(:), intent(in)    :: mudl
sll_real64, dimension(:), intent(in)    :: ae
sll_real64, dimension(:), intent(inout) :: as

sll_int32  :: kj, jid, jmi, ij, jj, id, imi, ii, ntest
sll_int32  :: i, j
sll_real64 :: s, xii

sll_real64, parameter  :: eps = 1.0e-10_f64

ntest = 0
as(1) = sqrt(ae(1))
ii    = mudl(1)

!--- 1.0 --- Matrice symetrique ---------------------------------------

do  i=2,size(mudl)

  xii = 0.0_f64
  imi = ii+1
  ii  = mudl(i)
  id  = i - ii
  jj  = mudl(imi+id-1)

  do ij = imi,ii-1
    j   = ij+id
    jmi = jj+1
    jj  = mudl(j)
    jid = j - jj -id
    s   = 0.0_f64

    do  kj = max( jmi , imi-jid )  ,  jj-1
      s = s + as( kj ) * as( kj + jid )
    end do

    as(ij) = ( ae(ij) - s ) / as(jj)
    xii    = xii + as(ij)*as(ij)
  end do 

  xii = ae(ii)  - xii
  if ( xii  <  eps*abs(ae(ii))) then
    write(6,900) i,eps
    write(6,901)xii,eps,ae(ii)
    ntest = 1
  end if

  as(ii) = sqrt ( xii )

end do

if(ntest==1) then 
  write(6,902) 
  stop "sll_m_choleski"
end if

!--- 9.0 --- Formats ---------------------------------------------------

 900 format(/10x,'resultats a verifier : le coefficient diagonal du dl'  &
             ,i6,' est inferieur au seuil de precision',e14.7)
 901 format(/10x,'xii:',e14.7,e14.7,e14.7)
 902 format(/10x,'************  Erreur dans sll_s_choles *************'/)

end subroutine sll_s_choles


!Function: sll_s_desrem
!descente et/ou remontee d'un systeme lineaire ( cholesky )
!
! in:
!  mudl,a,ntdl - pointeur,matrice et son ordre
!  be          - le second membre
!
! out:
!  bs          - le resultat
!
!  programmeur: 
!f hecht  - juin 84 inria , f. hermeline aout 89 cel/v
subroutine sll_s_desrem(mudl,a,be,ntdl,bs)

sll_int32,  intent(in)    :: mudl(0:*)
sll_real64, intent(in)    :: a(*)
sll_real64, intent(in)    :: be(*)
sll_int32,  intent(in)    :: ntdl
sll_real64, intent(inout) :: bs(*)

sll_int32   :: ii, ij, kj, il, i, j
sll_real64  :: y

ii = mudl(1)

!**********************************************************************
! matrice symetrique
!**********************************************************************
! descentes

do i=1,ntdl
   ij = ii + 1
   ii = mudl(i)
   kj = i  - ii
   y  = 0.0_f64
   do  il=ij,ii-1
      y = y +  a(il) * bs(il+kj)
   end do
   bs(i) = ( be(i) - y ) / a(ii)
end do

! remontees

ij = mudl(ntdl) + 1

do i = ntdl , 1 , -1
   ii = ij - 1
   ij = mudl(i-1) + 1
   kj =  ii - i
   bs(i) =  bs(i) / a(ii)
   do j = ij-kj , i - 1
      bs(j) = bs(j) - a(j+kj) * bs(i)
   end do
end do

end subroutine sll_s_desrem


end module sll_m_choleski
