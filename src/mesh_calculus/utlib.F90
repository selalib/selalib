!!File: utlib.f90
module utlib
use zone, only: iout
implicit double precision (a-h,o-z)

contains

double precision function utrand()

!************************************************************************
!*    But :    Genere un nombre aleatoire entre 0 et 1                  *
!************************************************************************

call random_number(utrand)

end function utrand

!-----------------------------------------------------------------

double precision function utfact(temps,pasdt,t0,t1,t2,t3,iforme)

integer :: npuiss, iforme
double precision :: t0, t1, t2, t3, coef
double precision :: temps, pasdt, dts3, dtm3, arg, del, pis2
double precision :: sigma, tau, omeg, xlg2, xlog2, ts, xpuiss
!* ------------------------------------------------------------- *
!*                                                               *
!* Definir une forme temporelle                                  *
!* ---                                                           *
!* temps  :  Temps courant                                       *
!* pasdt  :  Pas de temps                                        *
!* iforme :  Forme de l'impulsion                                *
!*          (0=pas de dependance                                 *
!*           1=trapeze                                           *
!*           2=sinusoides                                        *
!*           3=polynomes degres 3                                *
!*           4=super-gaussienne                                  *
!*           5=gaussienne                                        *
!*           6=fonction a spectre passe-bas de la forme:         *
!*                                                               *
!*          2                                                2   *
!*     f=sin (omeg(t-tau))/(omeg(t-tau)).exp(-((t-tau)/sigma) )  *
!*                                                               *
!*           7=creneau         )                          *
!*                                                               *
!*     Si 0 < iforme < 5 :                                       *
!*             t0     :  Temps initial (debut de la rampe)       *
!*             t1     :  Debut du plateau                        *
!*             t2     :  Fin du plateau                          *
!*             t3     :  Fin de l'impulsion                      *
!*     Si iforme = 5 :                                           *
!*             t0     :  centre de la gaussienne                 *
!*             t1     :  demi-largeur a mi-hauteur               *
!*     Si iforme = 6 :                                           *
!*             t0     :  omeg                                    *
!*             t1     :  tau                                     *
!*             t2     :  sigma                                   *
!*     Si iforme = 7 :                                           *
!*             t0     :  period                                    *
!*             t1     :  tau                                     *
!*             t2     :  sigma                                   *
!* --------------------------------------------------------------*

select case(iforme)

case(0)

   utfact=1.

case(1)

   del=pasdt/5.
   if ((temps>t0+del).and.(temps<t3-del)) then 
      if (temps<t1-del) then 
         utfact=(temps-t0)/(t1-t0)
      else if(temps>t2+del) then 
         utfact=(temps-t3)/(t2-t3)
      else
         utfact=1.
      end if 
   else
      utfact=0.
   end if 

case(2)

   del=pasdt/5.
   if ((temps>t0+del).and.(temps<t3-del)) then 
      pis2=asin(1.)
      if (temps<t1-del) then 
         utfact=sin(pis2*(temps-t0)/(t1-t0))**2
      else if(temps>t2+del) then 
         utfact=sin(pis2*(temps-t3)/(t2-t3))**2
      else
         utfact=1.
      end if 
   else
      utfact=0.
   end if 

case(3)

   del=pasdt/5.
   if ((temps>t0+del).and.(temps<t3-del)) then 
      if (temps<t1-del) then 
         dtm3=(t1-t0)*(t1-t0)*(t1-t0)
         utfact=-(temps-t0)*(temps-t0)*(2*temps+t0-3*t1)/dtm3
      else if(temps>t2+del) then 
         dts3=(t3-t2)*(t3-t2)*(t3-t2)
         utfact= (temps-t3)*(temps-t3)*(2*temps+t3-3*t2)/dts3
      else
         utfact=1.
      end if 
   else
      utfact=0.
   end if 

case(4)

   if ((t1-t0)<pasdt) then 
      if ((temps<t0).or.(temps>t2)) then 
         utfact=0.
      else
         utfact=1.
      end if 
   else
      npuiss=MAX(1,1+NINT((t2-t1)/(t1-t0)))
      xpuiss=npuiss
      xlog2=LOG(2.)
      if (npuiss>100) then 
         coef=xlog2/xpuiss
      else
         coef=2.**(1./xpuiss) - 1.
      end if 
      ts=temps-0.5*MAX(0.d0,t2-t1)-(t1-t0)
      arg=xlog2*(coef*((2.*(ts-t0)/(t1-t0))-1.))**(2*npuiss)
      utfact=exp(-arg)
   end if 

case(5)

   xlg2=LOG(2.)
   arg=xlg2*((temps-t0)/t1)**2
   arg=((temps-t0)/t1)**2/2
        
   if (arg>100) then 
      utfact=0.
   else
      utfact=exp(-arg)
   end if 

case(6)

   omeg=t0
   tau=t1
   sigma=t2
   del=1.e-2*pasdt

   if (abs(temps-tau)>del) then 
      utfact=(sin(omeg*(temps-tau)))**2/(omeg*(temps-tau))
   else
      utfact=0
   end if 

   utfact=utfact*exp(-((temps-tau)/sigma)**2)

case(7)    !correspond a shapeB = 4,5

   T     = 1./t0
   eta   = t1
   z0    = temps - T*floor(temps/T)
   if ( (z0 >= eta*T*0.5) .and. (z0 <=  T*(1-eta*0.5)) ) then
      utfact = 0.0  
   else 
      utfact = 1.0
   end if

case(8)    !correspond a shapeB = 6

   T     = 1./t0
   eta   = t1
   z0    = temps - T*floor(temps/T)
   if ( (z0 >= eta*T*0.5) .and. (z0 <=  T*(1-eta*0.5)) ) then
      utfact = -1.0  
   else 
      utfact = +1.0
   end if

case(9)    !correspond a shapeB = 7

   T   =1./t0
   eta =t1
   z  = temps - T*floor(temps/T)
   if ((z < eta*T*.25) .or. (z >=  T*(1-eta*.25))) then
      utfact = -1.0
   else  if ((z>=(1.-eta*.5)*T*.5) .and. &
             (z< (1.+eta*.5)*T*.5)) then
      utfact = 1.0
   else
      utfact = 0.0
   end if

case default

   utfact=1.

end select
                  
!                  select case(shapeB)
!                  case (1) !Continuous
!                      e1ext = 0.0
!                      e2ext = 0.0
!                      b3ext = Bmax;  
!                  case(2)  !Periodic focusing with positive  sinusoid  
!                      e1ext = 0.0
!                      e2ext = 0.0
!                      b3ext = Bmax*0.5*(1+cos(2*pi*time*freq));  
!                  case(3)  !Periodic focusing with  alternated sinusoid  
!                      e1ext = 0.0
!                      e2ext = 0.0
!                      b3ext = Bmax*cos(2*pi*time*freq);  
!                  case(4)  !periodic focusing with  positive, symetric lattice  
!                      e1ext = 0.0
!                      e2ext = 0.0
!                      if ( (pz >= S/4) .and. (pz <=  3*S/4) ) then
!                    b3ext = 0.0  
!                      else 
!                    b3ext = Bmax  
!                      end if
!                  case(5)  !periodic with positive  disymetric lattice  
!                      e1ext = 0.0
!                      e2ext = 0.0
!                      if ( (pz >= eta*S*0.5) .and. &
!               (pz <=  S*(1-eta*0.5)) ) then
!                          b3ext = 0.0  
!                      else 
!                          b3ext = Bmax      
!                      end if 
!                  case(6)  !periodic focusing with alternated lattice  
!                      e1ext = 0.0
!                      e2ext = 0.0
!                      if ( (pz >= S/4) .and. (pz <= 3*S/4) ) then
!                         b3ext = -Bmax
!                                  else 
!                         b3ext = Bmax
!                      end if
!                  case(7) 
!                     if ((pz < eta*S*0.25) .or. (pz >=  S*(1-eta*0.25))) then
!                        Bprime = -Bmax
!                     else  if ((pz>=(1.-eta*.5)*S*.5) .and. &
!                   (pz<(1.+eta*0.5)*S*0.5)) then
!                        Bprime =  Bmax
!                     end if
!                     e1ext =  -0.5*Bprime*vbeam
!                     e2ext =  +0.5*Bprime*vbeam
!                     b3ext = 0
!                  case default
!                     e1ext = xextep(1,i)
!                     e2ext = xextep(2,i)
!                     b3ext = xextep(3,i)
!                  end select
!

end function utfact

!-----------------------------------------------------------------------

subroutine utvppp(tt,val,ttab,vtab,ntab,ipp)
 
!************************************************************************
!*                                                                      *
!*         Version 1.0   Janvier 94  A.Adolf                            *
!*                                                                      *
!*         But :     Interpolation lineaire a partir de donnees         *
!*         ---       point par point.                                   *
!*                   Condition : la valeur de l'abscisse est            *
!*                   croissante d'un appel a l'autre.                   *
!*                   Application : formes temporelles manuelles.        *
!*                                                                      *
!*         Parametres d'entree  :                                       *
!*         --------------------                                         *
!*                     tt     :  temps courant                          *
!*                     ttab   :  tableau des abscisses                  *
!*                     vtab   :  tableau des ordonnees                  *
!*                     ntab   :  nombre de points                       *
!*                     ipp    :  numero du point utilise a l'appel      *
!*                               precedent                              *
!*                                                                      *
!*         Parametres de sortie :                                       *
!*         --------------------                                         *
!*                     val    :  valeur de la fonction en tt            *
!*                                                                      *
!************************************************************************
 
integer :: ntab, ipp, iq
double precision :: tt, val, t1, t2, v1, v2
double precision ::  ttab(*),vtab(*)
 
!--- 1.0 --- Tests ----------------------------------------------------

if ( (ipp<1).OR.(ipp>ntab) ) then 
   write(iout,901) ipp,ntab
   write(iout,900)
   call errout(iout,"F","uvppp"," ")
end if 

!--- 2.0 --- tt reste inferieur a la valeur pointee precedemment ------

if (tt<=ttab(1)) then 
   val=vtab(1)
else if(tt<=ttab(ipp)) then 
   t1=ttab(ipp-1)
   t2=ttab(ipp  )
   v1=vtab(ipp-1)
   v2=vtab(ipp  )
   val=v1+(tt-t1)*((v2-v1)/(t2-t1))

!--- 3.0 --- tt devient superieur a la valeur pointee precedemment ----
         
else if(tt>=ttab(ntab)) then 
   val=vtab(ntab)
else

   do iq=ipp,ntab
      if (ttab(iq)>=tt) then 
         ipp=iq
         t1=ttab(ipp-1)
         t2=ttab(ipp  )
         v1=vtab(ipp-1)
         v2=vtab(ipp  )
         val=v1+(tt-t1)*((v2-v1)/(t2-t1))
         exit
      end if 
   end do

end if 

! --- 4.0 --- Format ---------------------------------------------------

 900  FORMAT(//10x,'Erreur dans les interpolations des formes'  &
     &        /10x,'temporelles dans UTVPPP'            &
     &        /10x,'Appelez l''assistance technique...'///)
 901  FORMAT(//10x,'Mauvaise valeur du pointeur  ipp=',I10      &
     &        /10x,'Valeur min = 1   ;   Valeur max = ',I10/)

end subroutine utvppp



!subroutine: utranx
!Genere une distribution gaussienne de deviation 1   
!(de largeur a mi-hauteur egale a 1/e )               
!(la valeur moyenne est de 0.55 sur une demi-gaussienne)
!                                                      
! trng  - tableau recevant les nombres aleatoires 
! ntrng - dimension de trng                      
! trav  - tableau de travail                    
! ntrav - dimension de trav                    
!
!Version 1.0   
!
!Octobre 1994  A. Adolf  
!
subroutine utranx(trng,ntrng,trav,ntrav)

integer :: ntrng, ntrav, nt10, i, j, k
double precision :: trng(ntrng),trav(ntrav)
 
trng = 0.0; trav = 0.0

IF(ntrav >= ntrng*10) THEN
   nt10=ntrng*10
   CALL utrant(trav,nt10)
   DO i=1,10
      DO j=1,ntrng
         k=ntrng*(i-1)+j
         trng(j)=trng(j)+trav(k)-0.5
      enddo
   enddo
   DO j=1,ntrng
      trng(j)=trng(j)*0.755
   enddo
ELSE
   DO j=1,ntrng
      trng(j)=utrang()
   enddo
ENDIF

end subroutine utranx

      double precision function utrang()

!* -------------------------------------------------------------------- *
!*                                                                      *
!*    Version 1.0   Mai  1991   A. Adolf                                *
!*                                                                      *
!*    But :    Genere une distribution gaussienne de deviation 1        *
!*    ---      (de largeur a mi-hauteur egale a 1/e )                   *
!*             (la valeur moyenne est de 0.55 sur une demi-gaussienne)  *
!*                                                                      *
!* -------------------------------------------------------------------- *

      integer :: i, n
      double precision :: xi
      double precision :: rnd(10)
 
      n=10
      call utrant(rnd,n)

      xi = 0.0
      do i = 1, 10
         xi = xi + rnd(i) - 0.5
      end do

!     utrang = 1.12*xi  (erreur?)
      utrang =  0.755*xi
 
end function utrang

subroutine utrant(tab,n)

!************************************************************************
!*                                                                      *
!*    Version 1.0   Avril 1991  A. Adolf                                *
!*                                                                      *
!*    But :    Genere un tableau de n nombres aleatoires entre 0 et 1   *
!*    ---                                                               *
!*                                                                      *
!************************************************************************

integer :: i
integer, intent(in) :: n
double precision, dimension(:) :: tab
             
do i = 1, n
   call random_number(tab(i))
end do
                
end subroutine utrant

subroutine utunit(kunit,ierr)

!*********************************************************************
!*                                                                   *
!*      Version 1.0   Aout 1990   A.Adolf                            *
!*                                                                   *
!*      But :     trouver un numero d'unite libre                    *
!*      ---       connecte a aucun fichier.                          *
!*                                                                   *
!*      Parametres de sortie :                                       *
!*      --------------------                                         *
!*                                                                   *
!*                  kunit  :  numero d'unite du fichier              *
!*                  ierr   :  code d`erreur                          *
!*                                                                   *
!*********************************************************************

integer  :: kunit, ierr, iunit
logical  :: lopen

ierr=0

l1: do iunit=20,99

   INQUIRE(unit=iunit,opened=lopen)
   IF(lopen) cycle

   !on verifie que le fichier n'est pas preconnecte ..................

   open (iunit,status = 'scratch',err=732)
   close(iunit,status = 'delete' ,err=732)

   kunit=iunit
   exit 
   ierr=1
   732 continue
end do l1


end subroutine utunit

 
!subroutine: utexpc
! -------------------------------------------------------------------- *
!                                                                      *
!     Version 1.0 Fevrier  1994  A. Adolf                              *
!                                                                      *
!     Calcul des champs Ez-Er ou Bz-Br par expansion paraxiale         *
!     a partir des valeurs sur l'axe.                                  *
!                                                                      *
!                Bz(r,z) = Bz0(z) - 0.25*r*r*(d2 Bz0(z)/dz2)           *
!                Br(r,z) =        - 0.5*r*(d Bz0(z)/dz)                *
!                                                                      *
!     Arguments :      p1(i)   Position z des particules               *
!     -----------      p2(i)   Position r des particules               *
!                      b1p(i)  Champ Bz vu par les particules          *
!                      b2p(i)  Champ Br vu par les particules          *
!                      nbp     nombre de points p1-p2                  *
!                      np      nombre de points B sur l'axe            *
!                      nder    ordre de derivation                     *
!                      xmul    coefficient multiplicateur              *
!                      zin     tableau des abscisses                   *
!                      bz0     tableau des valeurs de B sur l'axe      *
!                                                                      *
! -------------------------------------------------------------------- *
subroutine utexpc(p1,p2,b1p,b2p,nbp,np,nder,xmul,zin,bz0)
 
double precision, intent(in) :: xmul
double precision :: p1(*),p2(*),b1p(*),b2p(*),zin(*),bz0(*)
double precision :: dz0, dzn
 
!--- 1.0 --- Quelques coefficients utiles ------------------------------
 
IF(nder.GE.2) THEN
   xder2=1.
ELSE
   xder2=0.
ENDIF
 
dz0 = zin(2 )-zin(1)
dzn = zin(np)-zin(np-1)
 
! --- 0.5 --- Points compris entre -~ et Zin(1)-dz0 --------------------

DO ip=1,nbp

   zz = p1(ip)
   rr = p2(ip)
   rr2= rr*rr

   IF(zz.LE.(zin(1)-dz0)) THEN
    b1p(ip)=b1p(ip)+xmul*bz0(1)
   ENDIF

end do

! --- 1.0 --- Points compris entre zin(1)-dz0 et zin(1) ----------------

DO ip=1,nbp
   zz = p1(ip)
   rr = p2(ip)
   rr2= rr*rr

   IF( (zz.GT.(zin(1)-dz0)).AND.(zz.LE.zin(1)) ) THEN

      dz3=zin(2)-zin(1)
      dz1=dz3
      dz2=dz3
      d2z1=(dz1+dz2)*0.5
      d2z2=(dz2+dz3)*0.5
      frac = (zz-(zin(1)-dz1))/dz1

      bz01=bz0(1)
      bz02=bz0(1)

      dbdz1=0.
      dbdz2=0.
      dbdz3=(bz0(2)-bz0(1))/dz3

      d2bdz21=(dbdz2-dbdz1)/d2z1
      d2bdz22=(dbdz3-dbdz2)/d2z2
 
      b1p(ip) = b1p(ip) + xmul * (bz01 - frac*0.25*rr2*d2bdz22*xder2)

      IF(nder.GE.1) THEN
         b2p(ip) = b2p(ip) - xmul*0.25*rr*frac*dbdz3
      ENDIF

   ENDIF
end do

! --- 2.0 --- Points compris entre zin(1) et zin(2) --------------------

      DO 20 ip=1,nbp
         zz = p1(ip)
         rr = p2(ip)
         rr2= rr*rr

         IF( (zz.GT.zin(1)).AND.(zz.LE.zin(2)) ) THEN

        dz2=zin(2)-zin(1)
        dz3=zin(3)-zin(2)
            dz1=dz2
        frac = (zz-zin(1))/dz2

            d2z1=(dz1+dz2)*0.5
            d2z2=(dz2+dz3)*0.5

            bz01=bz0(1)
        bz02=bz0(2)

            dbdz1=0.
            dbdz2=(bz0(2)-bz0(1))/dz2
            dbdz3=(bz0(3)-bz0(2))/dz3

            d2bdz21=(dbdz2-dbdz1)/d2z1
            d2bdz22=(dbdz3-dbdz2)/d2z2

            b1p(ip) = b1p(ip) + xmul *                  &
     &               ( (1-frac) * (bz01 - 0.25*rr2*d2bdz21*xder2)   &
     &                   + frac * (bz02 - 0.25*rr2*d2bdz22*xder2) )


            IF(nder.GE.1) THEN
               b2p(ip) = b2p(ip) - xmul*0.5*rr*     &
     &             (0.5*dbdz2 + 0.5*frac*dbdz3)
            ENDIF

         ENDIF
 20   CONTINUE

! --- 3.0 --- Points compris entre zin(2) et zin(np-1) -----------------

      DO 31 k=3,np-1
       DO 30 ip=1,nbp
         zz = p1(ip)
         rr = p2(ip)
         rr2= rr*rr

         IF( (zz.GT.zin(k-1)).AND.(zz.LE.zin(k)) ) THEN

            dz1=zin(k-1)-zin(k-2)
        dz2=zin(k  )-zin(k-1)
        dz3=zin(k+1)-zin(k  )
        frac = (zz-zin(k-1))/dz2

            d2z1=(dz1+dz2)*0.5
            d2z2=(dz2+dz3)*0.5

            bz01=bz0(k-1)
        bz02=bz0(k)

            dbdz1=(bz0(k-1)-bz0(k-2))/dz1
            dbdz2=(bz0(k  )-bz0(k-1))/dz2
            dbdz3=(bz0(k+1)-bz0(k  ))/dz3

            d2bdz21=(dbdz2-dbdz1)/d2z1
            d2bdz22=(dbdz3-dbdz2)/d2z2

            b1p(ip) = b1p(ip) + xmul *                  &
     &               ( (1-frac) * (bz01 - 0.25*rr2*d2bdz21*xder2)   &
     &                   + frac * (bz02 - 0.25*rr2*d2bdz22*xder2) )


            IF(nder.GE.1) THEN
               b2p(ip) = b2p(ip) - xmul*0.5*rr*     &
     &             (0.5*dbdz2 + 0.5*(1-frac)*dbdz1 + 0.5*frac*dbdz3)
            ENDIF

         ENDIF
 30    CONTINUE
 31   CONTINUE

! --- 4.0 --- Points compris entre zin(np-1) et zin(np) ----------------

      DO 40 ip=1,nbp
         zz = p1(ip)
         rr = p2(ip)
         rr2= rr*rr

         IF( (zz.GT.zin(np-1)).AND.(zz.LE.zin(np)) ) THEN

            dz1=zin(np-1)-zin(np-2)
        dz2=zin(np  )-zin(np-1)
        dz3=dz2
        frac = (zz-zin(np-1))/dz2

            d2z1=(dz1+dz2)*0.5
            d2z2=(dz2+dz3)*0.5

            bz01=bz0(np-1)
        bz02=bz0(np)

            dbdz1=(bz0(np-1)-bz0(np-2))/dz1
            dbdz2=(bz0(np  )-bz0(np-1))/dz2
            dbdz3= 0.

            d2bdz21=(dbdz2-dbdz1)/d2z1
            d2bdz22=(dbdz3-dbdz2)/d2z2

            b1p(ip) = b1p(ip) + xmul *                  &
     &               ( (1-frac) * (bz01 - 0.25*rr2*d2bdz21*xder2)   &
     &                   + frac * (bz02 - 0.25*rr2*d2bdz22*xder2) )


            IF(nder.GE.1) THEN
               b2p(ip) = b2p(ip) - xmul*0.5*rr*     &
     &             (0.5*dbdz2 + 0.5*(1-frac)*dbdz1)
            ENDIF

         ENDIF
 40   CONTINUE

! --- 5.0 --- Points compris entre zin(np) et zin(np)+dzn --------------

      DO 50 ip=1,nbp
         zz = p1(ip)
         rr = p2(ip)
         rr2= rr*rr

         IF( (zz.GT.zin(np)).AND.(zz.LT.(zin(np)+dzn)) ) THEN

            dz1=zin(np)-zin(np-1)
        dz2=dz1
        dz3=dz1
        frac = (zz-zin(np))/dz2

            d2z1=(dz1+dz2)*0.5
            d2z2=(dz2+dz3)*0.5

            bz01=bz0(np)

            dbdz1=(bz0(np)-bz0(np-1))/dz1
            dbdz2=0.

            d2bdz21=(dbdz2-dbdz1)/d2z1

            b1p(ip) = b1p(ip) + xmul *  (bz01 -     &
     &               (1-frac)*0.25*rr2*d2bdz21*xder2 )


            IF(nder.GE.1) THEN
               b2p(ip) = b2p(ip) - xmul*0.25*rr*(1-frac)*dbdz1
            ENDIF

         ENDIF
 50   CONTINUE

! --- 6.0 --- Points compris entre zin(np)+dzn et +~ -------------------

      DO 60 ip=1,nbp
         zz = p1(ip)
         rr = p2(ip)
         rr2= rr*rr

         IF(zz.GE.(zin(np)+dzn)) THEN
        b1p(ip)=b1p(ip)+xmul*bz0(np)
         ENDIF
 60   CONTINUE

end subroutine utexpc


double precision function utelrd(x,y,z)

!* -------------------------------------------------------------------- *
!*                                                                      *
!*    Version 1.0   Juillet 1999 Numerical Recipes p. 257               *
!*                                                                      *
!*    But :    Calculer l'integrale elliptique de Carlson de 2eme espece*
!*    ---                                                               *
!*    Arguments :      x     | Parametres tels que x, y et z ne soient  *
!*    -----------      y     | pas negatifs. Au plus 1 parametre parmi  *
!*                     z     | ces 3 peut etre nul                      *
!*                                                                      *
!* -------------------------------------------------------------------- *
double precision ::  x,y,z
double precision, parameter :: errtol=.05,c1=3./14,c2=1./6.,c3=9./22.
double precision, parameter :: c4=3./26.,c5=.25*c3,c6=1.5*c4
double precision :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac
double precision :: sqrtx,sqrty,sqrtz,somme,xt,yt,zt

xt=x
yt=y
zt=z
somme=0.
fac=1.
delx=1.e20
dely=0.
delz=0.
DO WHILE (MAX(ABS(delx),ABS(dely),ABS(delz)).GT.errtol)
   sqrtx=SQRT(xt)
   sqrty=SQRT(yt)
   sqrtz=SQRT(zt)
   alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
   somme=somme+fac/(sqrtz*(zt+alamb))
   fac=0.25*fac
   xt=0.25*(xt+alamb)
   yt=0.25*(yt+alamb)
   zt=0.25*(zt+alamb)
   ave=0.2*(xt+yt+3.*zt)
   delx=(ave-xt)/ave
   dely=(ave-yt)/ave
   delz=(ave-zt)/ave
enddo

ea=delx*dely
eb=delz*delz
ec=ea-eb
ed=ea-6.*eb
ee=ed+ec+ec
utelrd=3.*somme+fac*(1.+ed*(-c1+c5*ed-c6*delz*delz*ee)  &
      +delz*(c2*ee+delz*(-c3*ec+delz*c4*ea)))/(ave*SQRT(ave))

end function utelrd

double precision function utelrf(x,y,z)

!* -------------------------------------------------------------------- *
!*                                                                      *
!*    Version 1.0   Juillet 1999 Numerical Recipes p. 257               *
!*                                                                      *
!*    But :    Calculer l'integrale elliptique de Carlson de 1ere espece*
!*    ---                                                               *
!*    Arguments :      x     | Parametres tels que x, y et z ne soient  *
!*    -----------      y     | pas negatifs. Au plus 1 parametre parmi  *
!*                     z     | ces 3 peut etre nul                      *
!*                                                                      *
!* -------------------------------------------------------------------- *
double precision :: x,y,z
double precision, parameter  :: errtol=.08,third=1./3.,c1=1./24,c2=.1,c3=3./44.,c4=1./14.
double precision :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt

xt=x
yt=y
zt=z
delx=1.e20
dely=0.
delz=0.
DO WHILE (MAX(ABS(delx),ABS(dely),ABS(delz)).GT.errtol)
   sqrtx=SQRT(xt)
   sqrty=SQRT(yt)
   sqrtz=SQRT(zt)
   alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
   xt=0.25*(xt+alamb)
   yt=0.25*(yt+alamb)
   zt=0.25*(zt+alamb)
   ave=third*(xt+yt+zt)
   delx=(ave-xt)/ave
   dely=(ave-yt)/ave
   delz=(ave-zt)/ave
enddo

e2=delx*dely-delz**2
e3=delx*dely*delz
utelrf=(1.+(c1*e2-c2-c3*e3)*e2+c4*e3)/SQRT(ave)

end function utelrf


double precision function utelrc(x,y)

!* -------------------------------------------------------------------- *
!*                                                                      *
!*    Version 1.0   Juillet 1999 Numerical Recipes p. 258               *
!*                                                                      *
!*    But :    Calculer l'integrale elliptique de Carlson de 3eme espece*
!*    ---                                                               *
!*    Arguments :      x     | Parametres tels que x n'est pas negatif  *
!*    -----------      y     | et y non nul. Si y<0 on calcule la valeur*
!*                             principale de Cauchy                     *
!*                                                                      *
!* -------------------------------------------------------------------- *
double precision :: x,y
double precision, parameter :: errtol=.04,third=1./3.,c1=.3,c2=1./7.,c3=.375,c4=9./22
double precision :: alamb,ave,s,w,xt,yt

IF(y.GT.0.) THEN
   xt=x
   yt=y
   w=1
ELSE
   xt=x-y
   yt=-y
   w=sqrt(x)/sqrt(xt)
ENDIF

s=1.e20
DO WHILE (ABS(s).GT.errtol)
   alamb=2.*SQRT(xt)*SQRT(yt)+yt
   xt=0.25*(xt+alamb)
   yt=0.25*(yt+alamb)
   ave=third*(xt+yt+yt)
   s=(yt-ave)/ave
enddo

utelrc=w*(1.+s*s*(c1+s*(c2+s*(c3+s*c4))))/SQRT(ave)

end function utelrc


double precision function utelrj(x,y,z,p)
!* -------------------------------------------------------------------- *
!*                                                                      *
!*    Version 1.0   Juillet 1999 Numerical Recipes p. 258               *
!*                                                                      *
!*    But :    Calculer l'integrale elliptique de Carlson de 3eme espece*
!*    ---                                                               *
!*    Arguments :      x     | Parametres tels que x, y et z ne soient  *
!*    -----------      y     | pas negatifs. Au plus 1 parametre parmi  *
!*                     z     | ces 3 peut etre nul                      *
!*                     p       parametre tel que p<0 ou p>0             *
!*                             si p<0 on calcule la valeur principale   *
!*                             de Cauchy                                *
!*                                                                      *
!* -------------------------------------------------------------------- *
double precision :: p,x,y,z
double precision, parameter :: errtol=.05,c1=3./14,c2=1./3.,c3=3./22.
double precision, parameter :: c4=3./26.,c5=.75*c3,c6=1.5*c4,c7=.5*c2,c8=c3+c3
double precision :: a,alamb,alpha,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee
double precision :: fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,somme,tau,xt,yt,zt

somme=0.
fac=1.
IF(p.GT.0.) THEN
   xt=x
   yt=y
   zt=z
   pt=p
ELSE
   xt=MIN(x,y,z)
   zt=MAX(x,y,z)
   yt=x+y+z-xt-zt
   a=1./(yt-p)
   b=a*(zt-yt)*(yt-xt)
   pt=yt+b
   rho=xt*zt/yt
   tau=p*pt/yt
   rcx=utelrc(rho,tau)
ENDIF

delx=1.e20
dely=0.
delz=0.
delp=0.
DO WHILE (MAX(ABS(delx),ABS(dely),ABS(delz),ABS(delp)).GT.errtol)
   sqrtx=SQRT(xt)
   sqrty=SQRT(yt)
   sqrtz=SQRT(zt)
   alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
   alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
   beta=pt*(pt+alamb)**2
   somme=somme+fac*utelrc(alpha,beta)
   fac=0.25*fac
   xt=0.25*(xt+alamb)
   yt=0.25*(yt+alamb)
   zt=0.25*(zt+alamb)
   pt=0.25*(pt+alamb)
   ave=0.2*(xt+yt+zt+pt+pt)
   delx=(ave-xt)/ave
   dely=(ave-yt)/ave
   delz=(ave-zt)/ave
   delp=(ave-pt)/ave
enddo

ea=delx*(dely+delz)+dely*delz
eb=delx*dely*delz
ec=delp**2
ed=ea-3.*ec
ee=eb+2.*delp*(ea-ec)
utelrj=3.*somme+fac*(1.+ed*(-c1+c5*ed-c6*ee)  &
      &     +eb*(c7+delp*(-c8+delp*c4))    &
      &     +delp*ea*(c2-delp*c3)-c2*delp*ec)/(ave*SQRT(ave))
IF(p.LE.0.) THEN
   utelrj=a*(b*utelrj+3.*(rcx-utelrf(xt,yt,zt)))
ENDIF

end function utelrj

end module utlib
