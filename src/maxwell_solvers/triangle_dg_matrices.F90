!File: Module MatriceElem
! Calcule les matrices elementaires pour des elements finis 2D
! sur l'element de reference (0,0), (1,0), (0,1)
! en utilisant une base orthonormale de polynomes calculee a partir
! des polynomes de Jacobi
!
!  Contact : Pierre Navaro http://wwww-irma.u-strasbg.fr/~navaro
!
module triangle_dg_matrices 

implicit none
private


type, public :: ElementRef
   integer :: degre   ! degre du polynome associe
   real(8), dimension(:,:), pointer :: MassMat ! matrice de masse
   real(8), dimension(:,:), pointer :: MassMatInv ! inverse de la matrice de masse
   real(8), dimension(:,:), pointer :: DxMat,DyMat ! matrices derivees 2D
   real(8), dimension(:,:), pointer :: BndMassMat ! matrice de masse sur le bord
end type ElementRef

public :: AssMatElem

contains

!-------------------------------
! function JacobiP(N,alpha,beta,x)
! Evalue le polynome de jacobi de degre N de type (alpha,beta)
! renvoie la valeur du polynome P aux points x
!    
! Auteur : Eric Sonnendrucker 29/12/2005
  function JacobiP(N,alpha,beta,x)

    real(8), dimension(:) :: x       ! points d'evaluation du polynome
    real(8) :: alpha,beta            ! parametres definissant le polynome
    real(8), dimension(size(x,1)):: JacobiP ! valeur du polynomes aux points x
    ! Variables locales
    real(8), dimension(:,:), allocatable :: PL
    real(8) :: ipa, ipb, bb, dd, apbp2i, a2mb2
    integer :: i, n
 
    allocate(PL(N+1,size(x,1)))
    ! Initial values P_0(x) and P_1(x)
    PL(1,:) = 1.0
    if (N==0) then
       JacobiP=PL(1,:)
       return
    end if
    PL(2,:) = (alpha+beta+2)*x/2 + (alpha-beta)/2
    if (N==1) then
       JacobiP=PL(2,:)
       return
    end if
    ! utiliser la formule de recurrence pour les autres
    apbp2i = alpha+beta
    a2mb2 = alpha**2-beta**2
    ipa = alpha
    ipb = beta
    do i=1,N-1
       apbp2i = apbp2i + 2
       dd = 2*(i+1)*(apbp2i-i+1)*apbp2i
       ipa = ipa + 1
       ipb = ipb + 1
       bb = 2*ipa*ipb*(apbp2i+2)
       PL(i+2,:) = ((apbp2i+1)*(a2mb2 + apbp2i*(apbp2i + 2)*x) *PL(i+1,:)  &
            -bb*PL(i,:))/dd
   end do
    JacobiP = PL(N+1,:)
  end function JacobiP

  subroutine AssMatElem(Elem,k,x,y)
   ! -----------------
   ! Assemble les matrices sur l'element de referece pour un degre k donne
   ! Utilise les matrices de  Vandermonde avec des polynomes orthogonaux
   ! et des degres de libertes nodaux.
   ! Utilise les pol. orthonormaux de Proriol sur le triangle de reference
   ! Ces polynomes sont definis a partir des pol de Jacobi par
   ! P_m,n(x,y)=sqrt((2*m+1)*2*(n+m+1))*JacobiP(m,0,0,2*x/(1-y)-1)*(1-y)**m 
   !                    *JacobiP(n,2*m+1,0,2*y-1)
   !                          pour 0 <= x,y <= 1 et x+y <= 1
   ! Utilise les polynomes de Legendre (Jacobi(n,0,0)) normalises sur le bord
   !------------------
   type (ElementRef) :: Elem ! Element de reference  
   integer, intent(in) :: k
   ! variables locales
   real(8), dimension((k+1)*(k+2)/2, (k+1)*(k+2)/2) :: vdm, dxvdm, dyvdm, work1
   real(8), dimension(k+1,k+1) :: vdm1d
   real(8), dimension((k+1)*(k+2)/2) :: x,y   ! coordonnees des points
   real(8), dimension((k+1)*(k+2)/2) :: xi, eta, dum
   real(8), dimension((k+1)*(k+2)/2) :: jacp1, jacp2, djacp1, djacp2, djacp1x, jacp2y
   real(8), dimension((k+1)*(k+2)/2) :: work
   integer, dimension((k+1)*(k+2)/2) :: ipivot
   integer :: i,n,m,nn,Nk,info,ldwork,l
   real(8) :: coef

   ! On conserve le degre de l'espace de polynomes
   Elem%degre = k
   Nk = (k+1)*(k+2)/2

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Allocation des matrices  (FAIT DANS LE MODULE GALERKINE)
   !allocate(Elem%MassMat(Nk,Nk),stat=ierr)
   !if (ierr.ne.0) print*, 'error in allocation of Elem%MassMat'
   !allocate(Elem%MassMatInv(Nk,Nk),stat=ierr)
   !if (ierr.ne.0) print*, 'error in allocation of Elem%MassMatInv'
   !allocate(Elem%DxMat(Nk,Nk), stat=ierr)
   !if (ierr.ne.0) print*, 'error in allocation of Elem%DxMat'
   !allocate(Elem%DyMat(Nk,Nk), stat=ierr)
   !if (ierr.ne.0) print*, 'error in allocation of Elem%DyMat'
   !allocate(Elem%BndMassMat(k+1,k+1), stat=ierr)
   !if (ierr.ne.0) print*, 'error in allocation of Elem%BndMassMat'
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! On definit les points d'interp de Lagrange sur le triangle de reference
   nn = 0
   ! On commence par les points sur l'arete 1-2 (y=0)
   do n=0,k
      nn = nn +1
      x(nn) = dble(n)/k
      y(nn) = 0.0d0
   end do

   ! On continue avec les points sur l'arete 2-3 (x+y=1)
   do n=1,k
      nn = nn +1
      y(nn) = dble(n)/k
      x(nn) = 1.0d0-y(nn)
   end do

   ! puis avec les points sur l'arete 3-1 (x=0)
   do n=1,k-1
      nn = nn +1
      y(nn) = 1.0d0-dble(n)/k
      x(nn) = 0.0d0
   end do

   ! et on finit par les points interieurs
   do n=1,k-1
      do m=1,k-1-n
         nn = nn +1
         x(nn) = dble(n)/k
         y(nn) = dble(m)/k
      end do
   end do

   print"(//,a,/)", 'Positions des DDL sur le triangle de reference'
   do n=1,Nk
      print*,n,x(n),y(n)
   end do

   ! coordonnees pour le passage du carre au triangle
   do n = 1, Nk
      xi(n) = 2.d0*x(n)/(1.d0-y(n)+epsilon(y(n)))-1.0d0
      eta(n) = 2.d0*y(n)-1.0d0
   end do

   ! On calcule les matrices de Vandermonde des P et des derivees DxP et DyP
   do n = 0,k
      do m = 0,k-n
         nn = 1 + (n+m)*(n+m+1)/2 + n 
         jacp1 = JacobiP(m,0.d0,0.d0,xi)
         jacp2 = JacobiP(n,dble(2*m+1),0.0d0,eta)
         if (m > 0) then
            dum(:) = 1. !calcul de (1.0d0-y)**m-1
            do i = 1,m-1
               dum(:) = dum(:) * (1.0d0-y(:))
            end do
            djacp1  = JacobiP(m-1,1.d0,1.d0,xi)
            djacp1x = djacp1*dum(:) !(1.d0-y)**(m-1)
            jacp2y  = jacp2*dum(:)  !(1.0d0-y)**(m-1)
         else 
            djacp1  = 0.0
            djacp1x = 0.0
            do l = 1, Nk
               jacp2y(l) = jacp2(l)/(1.0d0-y(l)+epsilon(y(l)))
            end do
         end if
         if (n > 0) then
            djacp2 = JacobiP(n-1,real(2*m+2,8),1.d0,eta)
         else
            djacp2 = 0.0d0
         end if
         coef = dsqrt(dble((2*m+1)*2*(n+m+1)))
         dum(:) = 1. !calcul de (1.0d0-y)**m
         do i = 1,m
            dum(:) = dum(:) * (1.0d0-y(:))
         end do
         vdm(:,nn) = jacp1(:) * jacp2(:) * coef * dum(:)
         dxvdm(:,nn) = coef*(m+1)*djacp1x*jacp2
         dyvdm(:,nn) = coef*(((xi+1.d0)*0.5*(m+1)*djacp1-m*jacp1)*jacp2y &
                 + jacp1*dum(:)*(2*m+n+2)*djacp2)
      end do
   end do

   print*,'dxvdm'
   do n=1,(k+1)*(k+2)/2
      write(*,'(6e15.3)') (dxvdm(n,m) , m=1,(k+1)*(k+2)/2)
   end do
    
   print*,'dyvdm'
   do n=1,(k+1)*(k+2)/2
      write(*,'(6e15.3)') (dyvdm(n,m) , m=1,(k+1)*(k+2)/2)
   end do

   print*,'vdm'
   do n=1,(k+1)*(k+2)/2
      write(*,'(6e15.3)') (vdm(n,m) , m=1,(k+1)*(k+2)/2)
   end do

   ! calcule l'inverse de la matrice de masse qui vaut vdm*vdm^T
   Elem%MassMatInv(:,:) = 0.0
   do m=1,Nk
      do n=1,Nk
         do nn=1,Nk
            Elem%MassMatInv(m,n) = Elem%MassMatInv(m,n) + vdm(m,nn)*vdm(n,nn)
         end do
      end do
   end do
   ! inverse la matrice de Vandermonde en utilisant LAPACK
   ! l'inverse sera stockee dans vdm
   Nk =(k+1)*(k+2)/2
   ldwork = Nk
   call DGETRF(Nk, Nk, vdm, Nk, ipivot, info)
   call DGETRI(Nk, vdm, Nk, ipivot, work, ldwork, info)
   ! Calcule la matrice de masse = vdm^T*vdm ou vdm est maintenant 
   ! l'inverse de la matrice de Vandermonde
   Elem%MassMat(:,:) = 0.0
   do m=1,Nk
      do n=1,Nk
         do nn=1,Nk
            Elem%MassMat(m,n) = Elem%MassMat(m,n) + vdm(nn,m)*vdm(nn,n)
         end do
      end do
   end do

   ! On en deduit les matrices derivees : Dx = M^-1*(dxvdm*vdminv)^T*M
   work1 = transpose(matmul(dxvdm,vdm))
   Elem%DxMat = matmul(matmul(Elem%MassMatInv,work1),Elem%MassMat)
 
   work1 = transpose(matmul(dyvdm,vdm))
   Elem%DyMat = matmul(matmul(Elem%MassMatInv,work1),Elem%MassMat)

   ! Assemblage de la matrice de masse sur le bord
   !----------------------------------------------
   ! On commence a calculer la matrice de VanderMonde 1d

   ! Les points du bord sont les k+1 premiers, il faut les ramener
   ! de l'intervalle [0,1] à [-1,1] 
   do n=0,k
      vdm1d(:,n+1)=JacobiP(n,0.0d0,0.0d0,2.0d0*x(1:k+1)-1.0d0)
      vdm1d(:,n+1)=vdm1d(:,n+1)*sqrt(real(2*n+1,8))  ! normalisation
   end do

   ! inverse la matrice de Vandermonde sur le bord en utilisant LAPACK
   Nk =k+1
   ldwork = Nk
   call DGETRF(Nk, Nk, vdm1d, Nk, ipivot, info)
   call DGETRI(Nk, vdm1d, Nk, ipivot, work, ldwork, info)
   ! On en deduit la matrice de masse sur le bord
   Elem%BndMassMat(:,:)=matmul(transpose(vdm1d),vdm1d)

   print*,'Inverse of Mass Matrix'
   do n=1,(k+1)*(k+2)/2
      write(*,'(6f12.3)') (Elem%MassMatInv(n,m) , m=1,(k+1)*(k+2)/2)
   end do
   print*,'Mass Matrix'
   do n=1,(k+1)*(k+2)/2
      write(*,'(6f12.5)') (Elem%MassMat(n,m) , m=1,(k+1)*(k+2)/2)
   end do

   print*,'Product Matrix = Identity'
   if (k==1) then
      write(*,"(3f12.5)") matmul(Elem%MassMatInv(:,:),Elem%MassMat(:,:))
   else if (k==2) then
      write(*,"(6f12.5)") matmul(Elem%MassMatInv(:,:),Elem%MassMat(:,:))
   end if

   print*,'Dx Matrix'
   do n=1,(k+1)*(k+2)/2
      write(*,'(6f12.5)') (Elem%DxMat(n,m) , m=1,(k+1)*(k+2)/2)
   end do
   print*,'Dy Matrix'
   do n=1,(k+1)*(k+2)/2
      write(*,'(6f12.5)') (Elem%DyMat(n,m) , m=1,(k+1)*(k+2)/2)
   end do

   print*,'Boundary Mass Matrix'
   do n=1,k+1
      write(*,'(4g15.3)') (Elem%BndMassMat(n,m) , m=1,k+1)
   end do
   
   return
end subroutine AssMatElem

end module triangle_dg_matrices
