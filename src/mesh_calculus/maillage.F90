!File: Module Maillage
module maillage

use zone, only: nesp, iout, c, dt, pi, nesmx, lmodtm, titre, &
        ldtfrc, Mx, My

use solveurs_module, only: mesh_bound, degrek
use sll_triangular_meshes

implicit none

integer, private :: i, j, k, ifr

logical :: ldebug = .false.

!Variable: imxref
! entier designant la reference par defaut
integer :: imxref=99999999 

!Variable: nelin
!nombre de triangles internes
integer :: nelin   
!Variable: nefro
!nombre de triangles ayant 1 noeud sur une frontiere
integer :: nefro   
integer :: nelmatf
integer, dimension(:),   allocatable :: ipoint

!Variables:
! Caracteristiques des cotes situes sur les frontieres
! kelfro - element auquel appartient un cote frontiere
! kctfro - numero local de cote (1,2,3)                
! krefro - numero de reference du cote frontiere        
! ksofro - numeros des 2 sommets extremite du cote       
! vnofro - composantes du vecteur normal a la frontiere (vers l'interieur)


!Variables:
! Caracteristiques des frontieres internes                   
! nnofnt - nombre de noeuds sur les frontieres internes   
! ntrfnt - nombre total de triangles (ceux de droite)      
! ntrfrn - nombre de triangles par frontiere                
! ntrfrc - nombre cumule de triangles                        
! ncdfnt - cotes  Dirichlet sur les frontieres internes (VF)    

integer, dimension(:), allocatable :: ntrfrn, ntrfrc
integer :: nnofnt
integer, private, dimension(:),   allocatable :: itrfnt,ictfnt
integer, private, dimension(:,:), allocatable :: isofnt

!Variables:
!  Vecteurs tangeants
!  vtaux  - composante x des vecteurs tangeants         
!  vtauy  - composante y des vecteurs tangeants        
real(8), dimension(:),   allocatable :: vtaux, vtauy

!Variables:
!  Quantites liees au maillage    
!  xlml   - limite inferieure x du domaine           
!  xlmu   - limite superieure x du domaine          
!  ylml   - limite inferieure y du domaine         
!  ylmu   - limite superieure y du domaine        
real(8) :: xlml, xlmu, ylml, ylmu

integer, dimension (:), allocatable :: nctfro, nctfrp

!Variables: pour le solveur Galerkine Discontinue
integer, dimension(:,:), allocatable :: ngareelt
integer, dimension(:,:), allocatable :: ngnoeelt
integer, dimension(:,:), allocatable :: nsomare
real(8), dimension(:), allocatable :: rsf
integer, dimension (3,4) :: nnoeselt

real(8) :: cfl   !(reel = c*dt)

integer :: nmaill	!Entier caracteristique du maillage

CONTAINS

!========================================================================


!Subroutine: calmai
!
! Calculer toutes les quantites liees au maillage    
!                                                           
! Variables d'entree:                                         
!
!             iout - etiquette logique du fichier "listing" 
!                                                             
!             nbs  - nombre de noeuds                          
!             nbt  - nombre de triangles du maillage            
!             coor - coordonnees des noeuds           
!             refn - numeros de references des noeuds 
!             ntri - numeros des sommets              
!             vois - numeros des voisins des triangles
!             voiv - numeros des voisins des triangles
!             aire - aires des triangles              
!             base - integrales des fonctions de base  
!             nusd - numeros de sous-domaine            
!
!             npoel1 - pointeur du tableau "npoel2"       
!             npoel2 - numero des elements ayant un sommet en commun                       
!                                                       
!             nucfl - num de cotes ne satisfaisant pas    
!                      la condition CFL                           
!                                                                
!             petitl - petite longueur de reference         
!             grandl - grande longueur de reference       
!             imxref - grand nombre entier de reference  
!                                                       
!             nbtcot - nombre total de cotes           
!             nbcoti - nombre de cotes internes       
!             nbcfli - nbre de cotes internes ne satisfaisant pas CFL 
!             ncotcu - nombre de cotes internes et frontieres        
!             ncotcu - nombre de cotes cumules :                    
!             - 1 = nombre de cotes internes effectifs             
!             - 2 = 1 + nombre de cotes internes fictifs          
!             - 3 = 2 + nombre de cotes frontieres references 1  
!             - 4 = 3 + nombre de cotes frontieres references 2 , etc... 
!
!             nelmatf - nombre d'elements de la matrice profil Laplacien
!                                                                      
!             c   - vitesse de la lumiere dans le vide              
!             dt  - pas de temps                                   
!             ldtfrc - permet de forcer le calcul si dt trop grand
!                                                                
!Auteur:
!
!A. Adolf / L. Arnaud - Version 1.0  Octobre 1991 
subroutine calmai(mesh, vmsh, bcnd)

type(mesh_bound), intent(in)   :: bcnd
type(sll_triangular_mesh_2d), intent(inout) :: mesh
type(voronoi), intent(out)     :: vmsh

integer, dimension (:), allocatable :: indc
integer, dimension (:), allocatable :: nuctfr
integer, dimension (:), allocatable :: nucfl
integer, dimension (4) :: ieltmp

integer :: it, ntmp, id1, nct, iel, ind, iel1, iel2, nel
integer :: i1, i2, i3, is, is1, is2, is3, iv1, iv2, iv3
integer :: jel1, jel2, jel3, nel1, nel2, nel3
integer :: jelt, iare, ielt1, jelt1, ielt2, jelt2
integer :: isom1, isom2, isomloc1, isomloc2
integer :: nucti, nc, nbti, icot, jcot
integer :: nm1, nm2, n1, n2, num1, num2, ivois
integer :: ncv, nuctf, indv1, indv2, indn1, indn2
integer :: ic, ic1, ic2, itmp, numc1, numc2, icfl
integer :: nusop, ncndel, nd1, nd2, nv3, nv4
integer :: ntrobt, ind1, ind2, indv, ntiter, nbcflj, nbcflf
integer :: ict, ictcl, iref, jref, ntrfnt, nmxfrp1
integer :: keltmp, kcttmp, kretmp, ks1tmp, ks2tmp
integer :: nbcot, ict1, ict2, ictcl1, ictcl2, ie, ntrcf
integer :: iac, nbc, ncot1, ncot2, ncot3, nref
integer :: itm1, itm2, itm3, iesp, isl1, isl2
integer :: nuct, ktrtmp, nelprm, itrg1, itrg2, itr, ntr

real(8) :: lx1, lx2, ly1, ly2, x1, y1, x2, y2
real(8) :: det, syca, syba, xa, ya, xb, yb, xc, yc
real(8) :: dd, xd, yd, xo, yo, xm, ym, xtmp
real(8) :: prvc, xt1, xt2, yt1, yt2, xs1, xs2, ys1, ys2
real(8) :: xcc, ycc, somme, rc, ab, bc
real(8) :: sina, sinb, sinc, cosa, cosb, cosc, sin2a, sin2b, sin2c
real(8) :: xlvnmn, xlvnmx, xldnmn, xldnmx
real(8) :: ca, cb, cc, pabc, airemn, airemx, airtot, surf
real(8) :: xl, xr, xtm1, xtm2, yc1, yc2, yc3, yn1, yn2, yn3
real(8) :: x1s, y1s, x2s, y2s, x3s, y3s, xxs, yys
real(8) :: xfr1, xfr2, yfr1, yfr2

logical :: ldefr, ltemp, lcalsu, lmainv
logical :: lflag0,lflag1,lflag2,lflag3,lerr

integer :: itab(3,3), inoe1, inoe2, inoe3, inoe4, inoe5, inoe6
real(8) :: coef

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

!*** Calcul des longueurs de reference
if (ldebug) write(iout,*)"*** Calcul des longueurs de reference ***"

xlml = minval(mesh%coord(1,:))
xlmu = maxval(mesh%coord(1,:))
ylml = minval(mesh%coord(2,:))
ylmu = maxval(mesh%coord(2,:))

mesh%petitl = 1.e-04 * min(xlmu-xlml,ylmu-ylml)/sqrt(float(mesh%num_nodes))
mesh%grandl = 1.e+04 * max(xlmu-xlml,ylmu-ylmu)

!*** Correction des erreurs de precision ----------------------
!    pour les noeuds sur l'axe en axisymetrique.
      
!*** Calcul des aires des triangles
if (ldebug) write(iout,*)"*** Calcul des aires des triangles ***"

allocate(mesh%aire(mesh%num_cells)); mesh%aire=0.0

airtot = 0.

do it = 1, mesh%num_cells

   lx1 = mesh%coord(1,mesh%nodes(2,it))-mesh%coord(1,mesh%nodes(1,it))
   ly1 = mesh%coord(2,mesh%nodes(3,it))-mesh%coord(2,mesh%nodes(1,it))
   lx2 = mesh%coord(1,mesh%nodes(3,it))-mesh%coord(1,mesh%nodes(1,it))
   ly2 = mesh%coord(2,mesh%nodes(2,it))-mesh%coord(2,mesh%nodes(1,it))

   mesh%aire(it) = 0.5 * abs(lx1*ly1 - lx2*ly2)

   if( mesh%aire(it) <= 0. ) then
     write(iout,*) " Triangle : ", it
     write(iout,*) mesh%nodes(1,it), ":",mesh%coord(1:2,mesh%nodes(1,it))
     write(iout,*) mesh%nodes(2,it), ":",mesh%coord(1:2,mesh%nodes(2,it))
     write(iout,*) mesh%nodes(3,it), ":",mesh%coord(1:2,mesh%nodes(3,it))
     !call errout(iout,"F","maillage.f90","Aire de triangle negative")
   end if

   airtot = airtot + mesh%aire(it)

end do

!Calcul des integrales des fonctions de base
if (ldebug) write(iout,*)"*** Calcul des integrales des fonctions de base ***"

allocate(mesh%xbas(mesh%num_nodes)); mesh%xbas=0.0

xtmp = 2. * pi / 6.

do it = 1, mesh%num_cells

   is1 = mesh%nodes(1,it) 
   is2 = mesh%nodes(2,it) 
   is3 = mesh%nodes(3,it) 

   mesh%xbas(is1) = mesh%xbas(is1) + mesh%aire(it)/3.
   mesh%xbas(is2) = mesh%xbas(is2) + mesh%aire(it)/3.
   mesh%xbas(is3) = mesh%xbas(is3) + mesh%aire(it)/3.

end do

write(iout,"(/10x,'Longueurs de reference :',2E15.5/)") mesh%petitl,mesh%grandl
write(iout,"(/10x,'Limites x du domaine   :',2E15.5/    &
  &        10x,'Limites y du domaine   :',2E15.5/   &
  &        10x,'Aire des triangles     :', E15.5/)") xlml,xlmu,ylml,ylmu,airtot


! --- Gestion des triangles ayant un noeud en commun -----------
if (ldebug) &
write(iout,*)"*** Gestion des triangles ayant un noeud en commun ***"
 
! ... recherche des elements ayant un sommet commun
!     creation du tableau npoel1(i+1)  contenant le nombre de 
!     triangles ayant le noeud i en commun

allocate(mesh%npoel1(mesh%num_nodes+1))

mesh%npoel1 = 0
do i=1,mesh%num_cells
   is1 = mesh%nodes(1,i)
   is2 = mesh%nodes(2,i)
   is3 = mesh%nodes(3,i)
   mesh%npoel1(is1+1) = mesh%npoel1(is1+1)+1
   mesh%npoel1(is2+1) = mesh%npoel1(is2+1)+1
   mesh%npoel1(is3+1) = mesh%npoel1(is3+1)+1
end do

! ... le tableau npoel1 devient le tableau donnant l'adresse 
!     dans npoel2 du dernier element dans la suite des triangles
!     communs a un noeud

mesh%npoel1(1)=0
do i=3,mesh%num_nodes+1
   mesh%npoel1(i)=mesh%npoel1(i-1)+mesh%npoel1(i)
end do

! ... creation du tableau npoel2 contenant sequentiellement les 
!     numeros des triangles ayant un noeud en commun      
!     le premier triangle s'appuyant sur le noeud i est
!     adresse par "npoel1(i)+1" 
!     le nombre de triangles ayant le noeud i en commun est
!     "npoel1(i+1)-npoel1(i)"


allocate(mesh%npoel2(mesh%npoel1(mesh%num_nodes+1)))
allocate(indc(mesh%num_nodes))

indc   = 1  !Le tableau temporaire indc doit etre initialise a 1

do it = 1,mesh%num_cells
   do k = 1,3
      is = mesh%nodes(k,it)
      mesh%npoel2(mesh%npoel1(is)+indc(is)) = it
      indc(is) = indc(is)+1
   end do
end do
 
! --- Recherche des numeros des triangles voisins d'un triangle 

if (ldebug) then
   write(iout,*)"*** Recherche des numeros des triangles voisins d'un triangle ***"
   do i = 1, mesh%num_cells
      write(iout,*) " Triangle ", i, " Voisins :", mesh%nvois(1:3,i)
   end do
end if


do iel=1,mesh%num_cells

   ! ... numeros des 3 sommets du triangle

   is1=mesh%nodes(1,iel)
   is2=mesh%nodes(2,iel)
   is3=mesh%nodes(3,iel)

   ! ... boucles imbriquees sur les elements pointant vers
   !     les 2 noeuds extremites de l'arete consideree
   !     Le voisin est le triangle commun (hormis iel)

   ! ... premiere arete (entre le sommet is1 et is2)

   nel1=mesh%npoel1(is1+1)-mesh%npoel1(is1) !nb de triangles communs a is1
   nel2=mesh%npoel1(is2+1)-mesh%npoel1(is2) !nb de triangles communs a is2

   loop1:do i1=1,nel1
            jel1=mesh%npoel2(mesh%npoel1(is1)+i1) !premier triangle is1
            if(jel1.ne.iel) then
               do i2=1,nel2
                  jel2=mesh%npoel2(mesh%npoel1(is2)+i2)
                  if(jel2 == jel1) then
                     mesh%nvois(1,iel)  = jel1
                     exit loop1
              end if
           end do
            end if
     end do loop1

   ! ... deuxieme arete (entre le sommet is2 et is3)

   nel2=mesh%npoel1(is2+1)-mesh%npoel1(is2)
   nel3=mesh%npoel1(is3+1)-mesh%npoel1(is3)

   loop2:do i2=1,nel2
            jel2=mesh%npoel2(mesh%npoel1(is2)+i2)
            if(jel2 /= iel) then
               do i3=1,nel3
                  jel3=mesh%npoel2(mesh%npoel1(is3)+i3)
                  if(jel3 == jel2) then
                     mesh%nvois(2,iel)=jel2
             exit loop2
              end if
           end do
            end if
     end do loop2

   ! ... troisieme arete (entre le sommet is3 et is1)

   nel3=mesh%npoel1(is3+1)-mesh%npoel1(is3)
   nel1=mesh%npoel1(is1+1)-mesh%npoel1(is1)

   loop3:do i3=1,nel3
            jel3=mesh%npoel2(mesh%npoel1(is3)+i3)
            if(jel3 /= iel) then
               do i1=1,nel1
                  jel1=mesh%npoel2(mesh%npoel1(is1)+i1)
                  if(jel1 == jel3) then
                     mesh%nvois(3,iel)=jel3
             exit loop3
              end if
           end do
            end if
     end do loop3
end do

!Calcul de nvoif : Numero local dans le voisin de la face j commune a l'elt i

allocate(mesh%nvoif(3,mesh%num_cells))
mesh%nvoif(1:3,:) = mesh%nvois(1:3,:)

do i = 1, mesh%num_cells
   !write(*,"(/,4i7)") i, mesh%nodes(1:3,i)
   do j = 1,3
      jel1 = mesh%nvois(j,i)
      if (jel1 > 0) then
         do k = 1,3
            jel2 = mesh%nvois(k,jel1)
            if (jel2 == i) then 
               mesh%nvoif(j,i) = k
               exit
            end if
         end do
         !write(*,"(5i7)") jel1, mesh%nodes(1:3,jel1), mesh%nvoif(j,i)
      end if
   end do
end do

! --- Definition de nctfrt: le nombre de cotes frontieres

mesh%nctfrt=0
do i=1,mesh%num_cells
   if (mesh%nvois(1,i) < 0) mesh%nctfrt=mesh%nctfrt+1
   if (mesh%nvois(2,i) < 0) mesh%nctfrt=mesh%nctfrt+1
   if (mesh%nvois(3,i) < 0) mesh%nctfrt=mesh%nctfrt+1
end do

! --- Rangement de npoel2 dans l'ordre trigonometrique ---------

do is=1,mesh%num_nodes

   nel =mesh%npoel1(is+1)-mesh%npoel1(is)

   if ( nel > 1 ) then

      !*** Noeuds internes (Numero de reference nul) ***

      if( mesh%refs(is) == 0) then

        ind =1
        iel1=mesh%npoel2(mesh%npoel1(is)+1)

   loop4:do iel=2,nel-1
            do j=1,3
               if(mesh%nodes(j,iel1) == is) nct=mod(j+1,3)+1
            end do

            iel1=mesh%nvois(nct,iel1)
            do id1=ind+1,nel
               if(iel1 == mesh%npoel2(mesh%npoel1(is)+id1)) then
                  ind=ind+1
                  ntmp=mesh%npoel2(mesh%npoel1(is)+ind)
                  mesh%npoel2(mesh%npoel1(is)+ind)=iel1
                  mesh%npoel2(mesh%npoel1(is)+id1)=ntmp
              cycle loop4
               end if
            end do
         end do loop4

      ! Noeuds frontieres

      else 

      ! --> Recherche du premier triangle dans l'ordre trigonometrique

   loop5:do id1=1,nel
            iel1=mesh%npoel2(mesh%npoel1(is)+id1)
            do j=1,3
               if(mesh%nvois(j,iel1).le.0 .and. mesh%nodes(j,iel1) == is) then
                  ntmp=mesh%npoel2(mesh%npoel1(is)+1)
                  mesh%npoel2(mesh%npoel1(is)+1)=iel1
                  mesh%npoel2(mesh%npoel1(is)+id1)=ntmp
              exit loop5
               end if
            end do
         end do loop5
            
      ! --> Rangement des autres triangles dans l'ordre trigonometrique
      !     (s'il y en a plus que 2) 

         if(nel  > 2) then
            ind =1
            iel1=mesh%npoel2(mesh%npoel1(is)+1)
   
      loop6:do iel=2,nel-1
               do j=1,3
              if(mesh%nodes(j,iel1)==is) then
                 nct=mod(j+1,3)+1
              end if
               end do

               iel1=mesh%nvois(nct,iel1)
   
               do id1=ind+1,nel
                  if(iel1 == mesh%npoel2(mesh%npoel1(is)+id1)) then
                     ind=ind+1
                     ntmp=mesh%npoel2(mesh%npoel1(is)+ind)
                     mesh%npoel2(mesh%npoel1(is)+ind)=iel1
                     mesh%npoel2(mesh%npoel1(is)+id1)=ntmp
                 cycle loop6
              end if
               end do

            end do loop6

         end if

     end if

  end if

end do

!======================================================================
!----------- Nombre de noeuds sur les frontieres internes -------------
!======================================================================

nnofnt=0
if(bcnd%nbfrnt > 0) then
   do ifr=1,bcnd%nbfrnt
      xfr1=bcnd%z1frnt(ifr)-mesh%petitl
      xfr2=bcnd%z2frnt(ifr)+mesh%petitl
      yfr1=bcnd%r1frnt(ifr)-mesh%petitl
      yfr2=bcnd%r2frnt(ifr)+mesh%petitl
      do is=1,mesh%num_nodes
         xxs=mesh%coord(1,is)
     yys=mesh%coord(2,is)
     if( (xxs  > xfr1).and.(xxs.lt.xfr2).and.   &
             (yys  > yfr1).and.(yys.lt.yfr2) ) then
        nnofnt=nnofnt+1
         end if
      end do
   end do

   if(nnofnt < 2) then
      write(iout,"(//10x,'On ne trouve aucune arete ')")
      write(iout,"('sur les frontieres internes')")
      write(iout,"(10x,'(Verifier le namelist nlcham ou refaire')")
      write(iout,"(' le maillage)'/)")
      !call errout(iout,"F","maillage.f90","nnofnt < 2")
   end if

end if

do iel = 1, mesh%num_cells
   do j = 1, 3
   if( mesh%nvois(j,iel) == imxref ) then
      write(iout,*) " Triangle ", iel, " Voisins :", (mesh%nvois(i,iel),i=1,3)
      write(iout,*) " Coordonnees x =", mesh%coord(1,mesh%nodes(1,iel))
      write(iout,*) " Coordonnees y =", mesh%coord(2,mesh%nodes(1,iel))
      stop
   end if
   end do
end do

!Calcul des voisins pour les particules
if (ldebug) write(iout,*)"*** Calcul des voisins pour les particules ***"

! --- Remplissage de "nvoiv" -----------------------------------
!     Identique a "nvois" sauf pour les aretes appartenant a 
!     une frontiere. Le chiffre correspond ici a un code pour 
!     le traitement des conditions aux limites sur les 
!     particules, alors que dans "nvois" ce chiffre est 
!     l'oppose du numero de reference de la frontiere concernee

allocate(mesh%nvoiv(3,mesh%num_cells))

do i = 1,mesh%num_cells

   ! ... Cotes internes
   do j = 1, 3
      if (mesh%nvois(j,i)>0) then
         mesh%nvoiv(j,i) = mesh%nvois(j,i)
      end if
   end do

   ! ... Cotes frontieres
   do j = 1, 3
      if (mesh%nvois(j,i)<0) then

         select case (bcnd%ntypfr(-mesh%nvois(j,i)))

         case(0)
            mesh%nvoiv(j,i) = -1
         case(1)
            mesh%nvoiv(j,i) =  0
         case(2)
            mesh%nvoiv(j,i) = -1
         case(3)
            mesh%nvoiv(j,i) =  0
         case(4)
            mesh%nvoiv(j,i) =  0
         case(5)
            mesh%nvoiv(j,i) = -1
         case(6)
            mesh%nvoiv(j,i) = -1
         case(7) 
            mesh%nvoiv(j,i) =  0
         end select

      end if
   end do    

end do

! ... Verification approchee de la condition CFL sur les triangles

cfl = c*dt

!  nbtcot : nombre total de cotes                           *
!  nbcoti : nombre de cotes internes                        *
!  ncotcu : nombre de cotes internes et frontieres          *
!  ncotcu : nombre de cotes cumules :                       *
!  1 = nombre de cotes internes effectifs                   *
!  2 = 1 + nombre de cotes internes fictifs                 *
!  3 = 2 + nombre de cotes frontieres references 1          *
!  4 = 3 + nombre de cotes frontieres references 2 , etc... *

write(iout,"(/10x,a,i3)") 'Nombre maximum de frontieres referencees ', mesh%nmxfr

!*** Calcul du nb de cotes internes et total (nbcoti,nbtcot)

mesh%nbtcot = (3*mesh%num_cells+mesh%nctfrt)/2

write(iout,"( 10x,a,i6)") 'Nombre total de cotes =', mesh%nbtcot

mesh%nbcoti = mesh%nbtcot - mesh%nctfrt

write(iout,"( 10x,a,i6)") 'Nombre de cotes internes =', mesh%nbcoti
write(iout,"( 10x,a,i6/)")'Nombre de cotes frontieres =', mesh%nctfrt
 
!  Calcul du nb de cotes cumules par type de traitement (ncotcu) ....

allocate(vmsh%ncotcu(mesh%nbcoti+mesh%nmxfr*mesh%nctfrt))

vmsh%ncotcu=0
vmsh%ncotcu(1)=mesh%nbcoti
vmsh%ncotcu(2)=mesh%nbcoti

do ifr=1,mesh%nmxfr
   do i=1,mesh%num_cells
      do j = 1,3
         if(mesh%nvois(j,i) == -ifr) then
            vmsh%ncotcu(ifr+2) = vmsh%ncotcu(ifr+2) + 1
         end if
      end do
   end do
end do

do ifr=1,mesh%nmxfr
   vmsh%ncotcu(ifr+2)=vmsh%ncotcu(ifr+1)+vmsh%ncotcu(ifr+2)
end do

if (ldebug) then
   write(iout,"(10x,'Nombre de cotes cumules a partir de nbcoti :')")
   do ifr=1,mesh%nmxfr
      write(iout,"(20x,'ref = ',i3,'   ncotcu(i) = ',i6)") ifr,vmsh%ncotcu(ifr+2)
   end do
end if

allocate(vtaux(mesh%nbtcot),vtauy(mesh%nbtcot))
vtaux = 0.; vtauy = 0.

!==============================================================!
!--- Calcul complementaires relatifs au lissage des champs ----!
!    dans le cas de la resolution de POISSON               !
!==============================================================!

!tableau des numeros des noeuds associes a un cote 
allocate(vmsh%nuvac(2,mesh%nbtcot)); vmsh%nuvac = 0

!tableau des longueurs des cotes 
allocate(vmsh%xlcod(mesh%nbtcot)); vmsh%xlcod = 0.0

!tableaux de lissage et des cotes tangeants 
allocate(mesh%xmal1(mesh%num_nodes),mesh%xmal2(mesh%num_nodes),mesh%xmal3(mesh%num_nodes))
mesh%xmal1=0.;mesh%xmal2=0.;mesh%xmal3=0.

allocate(vmsh%nbcov(mesh%num_nodes+1))  !pointeur des cotes pointant sur le meme noeud
allocate(vmsh%nugcv(10*mesh%num_nodes)) !tableau contenant les numeros de ces cotes
allocate(nuctfr(mesh%nmxfr)) !tableau temporaire                              

!Creation du maillage de Voronoi et quantites associees

call poclis(mesh, vmsh, nuctfr, vtaux, vtauy)

!On tue le tableau temporaire         I                      

deallocate(nuctfr) 


!Event: Numerotation des cotes frontieres
!
!    Numerotation des cotes frontieres dans le sens 
!    trigonometrique, par numero de reference, pour
!    effectuer tous les diagnostics relatifs aux  
!    frontieres.                                 
!                                               
!    kelfro - numero d'element auquel appartient le cote 
!    kctfro - numero local (1,2,3) de cote              
!    krefro - numero de reference du cote              
!    ksofro - numeros des 2 sommets extremite du cote 
!    vnofro - composantes du vecteur normal (vers l'interieur)
!                                                            
!    nctfrt - Nb total de cotes frontieres                  
!    nctfro - Nb de cotes frontieres par reference         
!    nctfrp - Pointeur de tableau (nb de cote par ref)    
!                                                        

! --- Initialisation (numerotation quelconque) -----------------

allocate(nctfrp(0:mesh%nmxfr))
allocate(nctfro(mesh%nmxfr))
allocate(mesh%kctfro(mesh%nctfrt))
allocate(mesh%kelfro(mesh%nctfrt))
allocate(mesh%krefro(mesh%nctfrt))
allocate(mesh%ksofro(2,mesh%nctfrt))
allocate(mesh%vnofro(2,mesh%nctfrt))

nctfro = 0

!*** Premiere numerotation des cotes frontieres ***

ifr=0
do ict=1,3
   do iel=1,mesh%num_cells
      if ( mesh%nvois(ict,iel) < 0 ) then 

         iref = -mesh%nvois(ict,iel)

         ifr  = ifr+1
         mesh%kelfro(ifr) = iel  !element auquel appartient le cote
         mesh%kctfro(ifr) = ict  !numero local du cote dans cet element
         mesh%krefro(ifr) = iref !reference de ce cote
         mesh%ksofro(1,ifr) = mesh%nodes(ict,iel) !numeros des 2 sommets de ce cote
         mesh%ksofro(2,ifr) = mesh%nodes(mod(ict,3)+1,iel) 

         nctfro(iref)  = nctfro(iref)+1  !nombre de cote par reference

      end if 
   end do
end do

!... Pointeur de tableau .........................................

nctfrp(0)=0
do ifr=1,mesh%nmxfr
   nctfrp(ifr)=nctfrp(ifr-1)+nctfro(ifr)
end do

!--- Balayage sur les numeros de reference --------------------
!    Premier classement par numero de reference        

ictcl=1
do iref=1,mesh%nmxfr

   nbcot = nctfro(iref)
   if (nbcot > 0) then 

      ict1 = ictcl

      do ict=ict1,mesh%nctfrt
         jref = mesh%krefro(ict)

         if (jref == iref) then 

            keltmp = mesh%kelfro(ict)
            kcttmp = mesh%kctfro(ict)
            kretmp = mesh%krefro(ict)
            ks1tmp = mesh%ksofro(1,ict)
            ks2tmp = mesh%ksofro(2,ict)

            mesh%kelfro(ict) = mesh%kelfro(ictcl)
            mesh%kctfro(ict) = mesh%kctfro(ictcl)
            mesh%krefro(ict) = mesh%krefro(ictcl)
            mesh%ksofro(1,ict) = mesh%ksofro(1,ictcl)
            mesh%ksofro(2,ict) = mesh%ksofro(2,ictcl)

            mesh%kelfro(ictcl) = keltmp
            mesh%kctfro(ictcl) = kcttmp
            mesh%krefro(ictcl) = kretmp
            mesh%ksofro(1,ictcl) = ks1tmp
            mesh%ksofro(2,ictcl) = ks2tmp

            ictcl=ictcl+1

         end if 
      end do
   end if 
end do

! --- Rangement dans l'ordre trigonometrique -------------------

do iref=1,mesh%nmxfr

   nbcot=nctfro(iref)

   !Une seule reference pour toute la frontiere ................

   if(nbcot  ==  mesh%nctfrt) then 
      ict1=1
 35   continue
      is2=mesh%ksofro(2,ict1)
      do ict2=ict1,mesh%nctfrt
         if (ict1 /= ict2) then 
            is1=mesh%ksofro(1,ict2)
            if (is1==is2) then 
            
               keltmp = mesh%kelfro(ict1+1)
               kcttmp = mesh%kctfro(ict1+1)
               kretmp = mesh%krefro(ict1+1)
               ks1tmp = mesh%ksofro(1,ict1+1)
               ks2tmp = mesh%ksofro(2,ict1+1)

               mesh%kelfro(ict1+1) = mesh%kelfro(ict2)
               mesh%kctfro(ict1+1) = mesh%kctfro(ict2)
               mesh%krefro(ict1+1) = mesh%krefro(ict2)
               mesh%ksofro(1,ict1+1) = mesh%ksofro(1,ict2)
               mesh%ksofro(2,ict1+1) = mesh%ksofro(2,ict2)

               mesh%kelfro(ict2) = keltmp
               mesh%kctfro(ict2) = kcttmp
               mesh%krefro(ict2) = kretmp
               mesh%ksofro(1,ict2) = ks1tmp
               mesh%ksofro(2,ict2) = ks2tmp

               if (ict1<mesh%nctfrt) then 
              ict1=ict1+1
          GOTO 35
           end if

       end if
        end if
     end do

     ! ... Plusieurs references sur la frontiere ...........................

  else if (nbcot>1) then 

     ictcl1=nctfrp(iref-1)+1
     ictcl2=nctfrp(iref)

     ! ... Premier cote d'un segment ........................................

     ict1=ictcl1
 31  continue
     is1=mesh%ksofro(1,ict1)
     do ict2=ictcl1,ictcl2
        if (ict1 /= ict2) then 
           is2=mesh%ksofro(2,ict2)
           if (is1==is2) then 
              ict1=ict2
              if (ict1==ictcl1) then !Test si on tourne en rond ............
                 GOTO 37
              end if
              GOTO 31
           end if
        end if
     end do

! ... Permutation ......................................................
               
     keltmp=mesh%kelfro(ict1)
     kcttmp=mesh%kctfro(ict1)
     kretmp=mesh%krefro(ict1)
     ks1tmp=mesh%ksofro(1,ict1)
     ks2tmp=mesh%ksofro(2,ict1)

     mesh%kelfro(ict1)=mesh%kelfro(ictcl1)
     mesh%kctfro(ict1)=mesh%kctfro(ictcl1)
     mesh%krefro(ict1)=mesh%krefro(ictcl1)

     mesh%ksofro(1,ict1)=mesh%ksofro(1,ictcl1)
     mesh%ksofro(2,ict1)=mesh%ksofro(2,ictcl1)

     mesh%kelfro(ictcl1)=keltmp
     mesh%kctfro(ictcl1)=kcttmp
     mesh%krefro(ictcl1)=kretmp

     mesh%ksofro(1,ictcl1)=ks1tmp
     mesh%ksofro(2,ictcl1)=ks2tmp
        
! ... Classement des cotes suivants ....................................

 37  continue
     ict1=ictcl1
 33  continue
     is2=mesh%ksofro(2,ict1)
     do ict2=ict1,ictcl2
        if (ict1 /= ict2) then 
           is1=mesh%ksofro(1,ict2)
           if (is1==is2) then 
        
              keltmp=mesh%kelfro(ict1+1)
              kcttmp=mesh%kctfro(ict1+1)
              kretmp=mesh%krefro(ict1+1)
              ks1tmp=mesh%ksofro(1,ict1+1)
              ks2tmp=mesh%ksofro(2,ict1+1)

              mesh%kelfro(ict1+1)=mesh%kelfro(ict2)
              mesh%kctfro(ict1+1)=mesh%kctfro(ict2)
              mesh%krefro(ict1+1)=mesh%krefro(ict2)
              mesh%ksofro(1,ict1+1)=mesh%ksofro(1,ict2)
              mesh%ksofro(2,ict1+1)=mesh%ksofro(2,ict2)

              mesh%kelfro(ict2)=keltmp
              mesh%kctfro(ict2)=kcttmp
              mesh%krefro(ict2)=kretmp
              mesh%ksofro(1,ict2)=ks1tmp
              mesh%ksofro(2,ict2)=ks2tmp

              ict1=ict1+1
              GOTO 33
           end if
        end if
     end do

! ...Test s'il y a d'autres segments de la meme reference .............

     if (ict1<ictcl2) then 
        ictcl1=ict1+1
        ict1=ictcl1
        GOTO 31
     end if
   end if
end do

!  Modification de nvois ------------------------------------
! (Numero de reference change par le numero de cote)
 
do ict=1,mesh%nctfrt
   ie=mesh%kelfro(ict)
   ic=mesh%kctfro(ict)
   mesh%nvois(ic,ie)=-ict
end do

!--- Calcul des composantes du vecteur normal -----------------

do ict=1,mesh%nctfrt

   is1 = mesh%ksofro(1,ict)
   is2 = mesh%ksofro(2,ict)

   x1 = mesh%coord(1,is1)
   y1 = mesh%coord(2,is1)
   x2 = mesh%coord(1,is2)
   y2 = mesh%coord(2,is2)

   mesh%vnofro(1,ict) = -y2+y1
   mesh%vnofro(2,ict) =  x2-x1

end do

deallocate(indc)

!Event: Frontieres Internes
!
!             Mise en oeuvre des frontieres internes. 
!                - identification des noeuds et cotes          
!                - modifications des tableaux                   
!                - traitement des conditions sur les champs      
!                                                                 
!             nctfnt - numeros globaux des cotes sur ces frontieres
!                                                                   
!             itrfnt - tableau temporaire                            
!                      (triangle s'appuyant sur une frontiere interne)
!             ictfnt - tableau temporaire                            
!                      (numreo local du cote sur la frontiere interne)
!             isofnt - tableau temporaire                             
!                      (numreo des 2 sommets de chaque triangle)       
!

if (bcnd%nbfrnt > 0 ) then

   write(iout,911) nnofnt
 
   allocate(mesh%noefnt(nnofnt))
   allocate(itrfnt(nnofnt))
   allocate(mesh%nctfnt(nnofnt))
   allocate(mesh%irffnt(nnofnt))
   allocate(ictfnt(nnofnt))
   allocate(isofnt(2,nnofnt))
   allocate(ntrfrn(bcnd%nbfrnt),ntrfrc(0:bcnd%nbfrnt))
   
   !----------- Tests preliminaires --------------------------------------
   
   do ifr=1,bcnd%nbfrnt
      if(bcnd%irefnt(ifr) <= mesh%nmxfr .or. bcnd%irefnt(ifr) > mesh%nmxfr) then
         nmxfrp1=mesh%nmxfr+1
         write(iout,901) ifr,bcnd%irefnt(ifr),nmxfrp1,mesh%nmxfr
         !call errout(iout,"F",'maillage.f90',' ligne 2276 ')
      end if
   end do
   
   !----------- Determination des triangles sur une frontiere interne ----
   !            (on ne prend que ceux de droite ou du bas)
   
   ntrfnt=0
   
   do ifr=1,bcnd%nbfrnt
   
      xfr1=bcnd%z1frnt(ifr)-mesh%petitl
      xfr2=bcnd%z2frnt(ifr)+mesh%petitl
      yfr1=bcnd%r1frnt(ifr)-mesh%petitl
      yfr2=bcnd%r2frnt(ifr)+mesh%petitl
   
      ntrfrn(ifr)=0
   
      do iel=1,mesh%num_cells
   
         !... numeros des somnmets
   
         is1=mesh%nodes(1,iel)
         is2=mesh%nodes(2,iel)
         is3=mesh%nodes(3,iel)
   
         !... coordonnees des sommets des triangles
   
         x1s=mesh%coord(1,is1)
         y1s=mesh%coord(2,is1)
         x2s=mesh%coord(1,is2)
         y2s=mesh%coord(2,is2)
         x3s=mesh%coord(1,is3)
         y3s=mesh%coord(2,is3)
   
         lflag0=.false.
         lflag1=.false.
         lflag2=.false.
         lflag3=.false.
   
         ! noeud 1 sur la frontiere
   
         if( (x1s  > xfr1).and.(x1s <  xfr2).and.   &
             (y1s  > yfr1).and.(y1s <  yfr2) ) then
            lflag1=.true.
         end if
   
         ! noeud 2 sur la frontiere
   
         if( (x2s  > xfr1).and.(x2s <  xfr2).and.   &
             (y2s  > yfr1).and.(y2s <  yfr2) ) then
            lflag2=.true.
         end if
   
         ! noeud 3 sur la frontiere
   
         if( (x3s  > xfr1).and.(x3s <  xfr2).and.   &
             (y3s  > yfr1).and.(y3s <  yfr2) ) then
               lflag3=.true.
         end if
   
         ! tests si le triangle est sur la frontiere interne
         ! et reperage des sommets et cotes
   
         if(lflag1.and.lflag2) then
            isl1=is1
        isl2=is2
        ict=1
        lflag0=.true.
         else if(lflag2.and.lflag3) then
            isl1=is2
            isl2=is3
            ict=2
            lflag0=.true.
         else if(lflag3.and.lflag1) then
            isl1=is3
            isl2=is1
            ict=3
            lflag0=.true.
         end if
   
         ! test si le triangle n'est pas sur une frontiere externe
         ! et s'il est a droite ou en bas de la frontiere interne
   
         if (lflag0) then
            if((mesh%nvois(ict,iel) > 0).and. (             &
                 (mesh%coord(2,isl1) > mesh%coord(2,isl2)+mesh%petitl).or.    &
                 (mesh%coord(1,isl1) > mesh%coord(1,isl2)+mesh%petitl) ) ) then
               ntrfnt=ntrfnt+1
               ntrfrn(ifr)=ntrfrn(ifr)+1
               itrfnt(ntrfnt)=iel
               ictfnt(ntrfnt)=ict
               isofnt(1,ntrfnt)=isl1
               isofnt(2,ntrfnt)=isl2
            end if
         end if
      end do
   end do
        
   !Pointeur (nb de triangles par frontiere)
   
   ntrfrc(0)=0
   do ifr=1,bcnd%nbfrnt
      ntrfrc(ifr)=ntrfrc(ifr-1)+ntrfrn(ifr)
   end do
   
   write(iout,905) ntrfnt
   
   if(ntrfnt == 0) then
      write(iout,909) 
      !call errout(iout,"F",'maillage.f90',' ')
   end if
   
   do i=1,bcnd%nbfrnt
   
      if(ntrfrn(i) == 0) then
         write(iout,910) 
         !call errout(iout,"F",'maillage.f90',' ')
      end if
   
      write(iout,906) i,ntrfrn(i),ntrfrc(i),bcnd%itfrnt(i),bcnd%irefnt(i)

   end do
   
   !----------- On ordonne le tableau dans le sens decroissant -----------

   do ifr=1,bcnd%nbfrnt
      ntr=ntrfrn(ifr)
      if (ntr > 1) then
   
         !numeros locaux des triangles de la frontiere ifr
   
         itrg1=ntrfrc(ifr-1)+1
         itrg2=ntrfrc(ifr)
   
         !Recherche du premier triangle pour chaque frontiere
         !et test de discontinuite
   
         nelprm=0
   
         do iel1=itrg1,itrg2
   
            !sommet du haut (ou de droite) d'un 1er triangle
   
        is1=isofnt(1,iel1)
        lflag0=.true.
   
            do iel2=itrg1,itrg2
   
           if(iel1 /= iel2) then
   
               !sommet du bas (ou de gauche) d'un 2eme triangle
   
              is2=isofnt(2,iel2)
   
                  !les deux sommets sont identiques donc le triangle 1
                  !n'est pas le plus haut ou le plus a gauche
   
              if(is1 == is2) then
                 lflag0=.false.
              end if
   
               end if
            end do
   
            !Si on a trouve deux fois un triangle extremite 
            !alors la frontiere est discontinue (c'est interdit)

        if((lflag0).and.(nelprm /= 0)) then
           write(iout,907) ifr,bcnd%z1frnt(ifr),bcnd%r1frnt(ifr),bcnd%irefnt(ifr)
               !call errout(iout,"F",'maillage.f90',' ')
            end if

            !le triangle 1 est bien le plus haut

        if((lflag0).and.(nelprm == 0)) then
           nelprm=iel1
            end if

         end do

         !Permutation pour mettre le triangle 1 en 1ere position
   
         ktrtmp=itrfnt(nelprm)
         kcttmp=ictfnt(nelprm)
         ks1tmp=isofnt(1,nelprm)
         ks2tmp=isofnt(2,nelprm)

         itrfnt(nelprm)=itrfnt(itrg1)
         ictfnt(nelprm)=ictfnt(itrg1)
         isofnt(1,nelprm)=isofnt(1,itrg1)
         isofnt(2,nelprm)=isofnt(2,itrg1)

         itrfnt(itrg1)=ktrtmp
         ictfnt(itrg1)=kcttmp
         isofnt(1,itrg1)=ks1tmp
         isofnt(2,itrg1)=ks2tmp

         !Classement des cotes suivants

         iel1=itrg1

         135  continue

         is2=isofnt(2,iel1)
         do iel2=iel1,itrg2
            if(iel1 /= iel2) then
               is1=isofnt(1,iel2)
               if(is1 == is2) then

                  ktrtmp=itrfnt(iel1+1)
                  kcttmp=ictfnt(iel1+1)
                  ks1tmp=isofnt(1,iel1+1)
                  ks2tmp=isofnt(2,iel1+1)
   
                  itrfnt(iel1+1)=itrfnt(iel2)
                  ictfnt(iel1+1)=ictfnt(iel2)
                  isofnt(1,iel1+1)=isofnt(1,iel2)
                  isofnt(2,iel1+1)=isofnt(2,iel2)
   
                  itrfnt(iel2)=ktrtmp
                  ictfnt(iel2)=kcttmp
                  isofnt(1,iel2)=ks1tmp
                  isofnt(2,iel2)=ks2tmp
   
                  iel1=iel1+1
              GOTO 135
               end if
        end if
         end do

      end if
   end do

   ! ----------- Numeros des noeuds et references ( TYPE 1)        --------
   !             pour les champs
   lerr   = .false.
   mesh%nndfnt = 0
   do ifr=1,bcnd%nbfrnt

      ntr=ntrfrn(ifr)
   
      if (ntr  > 0.and.(bcnd%itfrnt(ifr) == 1)) then
         itrg1=ntrfrc(ifr-1)+1
         itrg2=ntrfrc(ifr)
         do iel=itrg1,itrg2 
            mesh%nndfnt=mesh%nndfnt+1
            mesh%noefnt(mesh%nndfnt)=isofnt(1,iel)
            mesh%irffnt(mesh%nndfnt)=bcnd%irefnt(ifr)

            !Numero global des cotes sur une frontiere interne (VF)
            !avec test CFL pour eviter les couplages entre 2 domaines
               
         enddo
         mesh%nndfnt=mesh%nndfnt+1
         mesh%noefnt(mesh%nndfnt)=isofnt(2,itrg2)
         mesh%irffnt(mesh%nndfnt)=bcnd%irefnt(ifr)
      endif
   enddo

   if (lerr) then
      !call errout(iout,"F",'maillage.f90',' ')
   endif

   ! ======================================================================
   ! ----------- Ecriture des resultats -----------------------------------
      
   write(iout,902)
   if(mesh%nndfnt  > 0) then
      do i=1,mesh%nndfnt
         write(iout,903) i,mesh%noefnt(i),mesh%coord(2,mesh%noefnt(i)),mesh%irffnt(i)
      end do
   else
      write(iout,904)
   end if

   !Event: Frontieres Internes (suite)
   !
   !             Mise en oeuvre des frontieres internes.      
   !             traitement des conditions sur les particules
   !                                                        
   !             itrfnt - tableau temporaire               
   !                      (triangle s'appuyant sur une frontiere interne) 
   !             ictfnt - tableau temporaire                             
   !                      (numreo local du cote sur la frontiere interne)
   !             isofnt - tableau temporaire                            
   !                      (numreo des 2 sommets de chaque triangle)    
   !             nctfrt - Nb total de cotes frontieres                
   !             nctfro - Nb de cotes frontieres par reference       
   !             nctfrp - Pointeur de tableau (nb de cote par ref)  
   !                                                               
   !
      
   do iesp = 1, nesp

      !----------- Tests preliminaires --------------------------------------
      do ifr=1,bcnd%nbfrnt
         if (bcnd%ipfrnt(ifr,iesp).lt.1.and.bcnd%ipfrnt(ifr,iesp)  > 2) then
             write(iout,912) iesp,ifr,bcnd%ipfrnt(ifr,iesp),1,2
             !call errout(iout,"F","maillage.f90"," ")
         endif
      enddo

      !----------- Modifications de nvoiv pour les conditions de conducteur
      !            parfait pour les champs et des particules absorbees
      !          
      do ifr=1,bcnd%nbfrnt
         ntr=ntrfrn(ifr)
         if (ntr  > 0.and.bcnd%ipfrnt(ifr,iesp) == 2) then
            itrg1=ntrfrc(ifr-1)+1
            itrg2=ntrfrc(ifr)
            do itr=itrg1,itrg2 
      
               iel1=itrfnt(itr)
               ict1=ictfnt(itr)
               iel2=mesh%nvois(ict1,iel1)
      
               is=isofnt(1,itr)
      
               is1=mesh%nodes(1,iel2)
               is2=mesh%nodes(2,iel2)
               is3=mesh%nodes(3,iel2)
      
               if (is1 == is) then
                  ict2=3
               elseif(is2 == is) then
                  ict2=1
               elseif(is3 == is) then
                  ict2=2
               else
                  write(iout,913)
                  !call errout(iout,"F","maillage.f90",' ')
               endif
   
               mesh%nvoiv(ict1,iel1)=0
               mesh%nvoiv(ict2,iel2)=0
   
            enddo
         endif
      end do
   end do

end if  !Fin du traitement lie au frontieres internes

901 format(//10x,'Mauvaise reference sur une frontiere interne'     &
            /10x,'(voir namelist  nlcham)'              &
            /10x,'Frontiere interne N0 ',I3             &
            /10x,'Reference            ',I3             &
            /10x,'Limites autorisees  ',2I4/)
902 format( /10x,'Noeuds Dirichlet des frontieres de type 1'/)
903 format( /10x,I3,5X,'N0 ',I7,5X,'Ordonnee',E12.3,5x,'Ref ',I3)
904 format( /20x,'Aucun noeud de ce type'/)
905 format( /10x,'Nombre total de triangles sur les frontieres'     &
            /10x,'internes :',I6)
906 format( /10x,'N0 de frontiere',I3                   &
            /15x,'Nb de triangles ',I6,5x,'Nb cumule ',I6       &
            /15x,'Type            ',I6,5x,'Reference ',I6)
907 format(//10x,'La frontiere interne',I4,'  est discontinue'      &
            /10x,'z1frnt :',E12.3,5x,                   &
            /10x,'r1frnt :',E12.3,5x,                   &
            /10x,'irefnt :',I4,                     &
            /10x,'Revoir les donnees ou le maillage'//)
908 format(//10x,'Detection d''un probleme de voisins'/)
909 format(//10x,'Aucun triangle ne repose sur les'         &
            /10x,'frontieres internes'/)
910 format(//10x,'Aucun triangle ne repose sur la'          &
            /10x,'frontiere interne No :',I5/)
911 format( /10x,'*****  FRONTIERES INTERNES  *****'            &
           //10x,'Nombre max de noeuds sur ces frontieres :',I6/)
921 format(  10x,'Pb de CFL sur cote',I8,'  Triangle',I8        &
            /10x,'Frontiere interne No :',I5)
912 format(//10x,'Mauvaise reference sur une frontiere interne'     &
            /10x,'pour l espece de particules NO ',I3           &
            /10x,'(voir namelist  nlcham)'              &
            /10x,'Frontiere interne N0 ',I3             &
            /10x,'Reference            ',I3             &
            /10x,'Limites autorisees  ',2I4/)
913 format(//10x,'Detection d''un probleme de voisins'/)


end subroutine calmai

!**************************************************************


!Subroutine: poclis
!                                               
!   Calcul des matrices de lissage associees a chaque     
!   noeud du maillage.                                  
!   Necessaires au calcul des composantes de E1,E2     
!   a partir des composantes tangeantielles de E      
!   connues sur les cotes des triangles.             
!   Calcul des composantes des vecteurs unitaires tangeants.                                     
!                                                              
!   Variables d'entree:                                    
!  
!      nuvac  - numeros des PV associes aux cotes          
!      coor   - coordonnees des noeuds Delaunay           
!      xlcod  - longueur des cotes Delaunay              
!      npoel1 - pointeur du tableau npoel2              
!      npoel2 - numeros des triangles entourant un noeud       
!      xmal1  - somme des taux*taux entourant un noeud (/det) 
!      xmal2  - somme des tauy*tauy entourant un noeud (/det)
!      xmal3  - somme des taux*tauy entourant un noeud (/det)
!      vtaux  - composante x des vecteurs tangeants         
!      vtauy  - composante y des vecteurs tangeants        
!                                                               
!Auteur:
! A. Adolf - Version 1.0   Septembre 1994  
subroutine poclis(mesh, vmsh, nuctfr, vtaux, vtauy)

type(sll_triangular_mesh_2d) :: mesh
type(voronoi)   :: vmsh
double precision, dimension(:) :: vtaux, vtauy
integer, dimension(:) :: nuctfr
logical :: lerr
double precision :: det, s1, s2, s3, x21, y21, xa, ya, xb, yb
integer :: ic, nuctf, ind, nm1, nm2, nel1, nel2, indv1, indv2
integer :: indn1, indn2, num1, num2, n1, n2, ivois, nc, iel, nucti
integer :: nbti, is
     
!======================================================================
! --- 1.0 --- Pointeur des numeros de cotes pointant vers un noeud -----

do is=1,mesh%num_nodes
   nbti=mesh%npoel1(is+1)-mesh%npoel1(is)
   vmsh%nbcov(is+1)=nbti
   if(mesh%refs(is) /= 0) then
      vmsh%nbcov(is+1)=nbti+1
   end if
end do

vmsh%nbcov(1)=0
do is=1,mesh%num_nodes
   vmsh%nbcov(is+1)=vmsh%nbcov(is)+vmsh%nbcov(is+1)
end do

! --- 1.5 --- Tableau temporaire (cumul des cotes frontieres) ----------
nucti=0
do i=1,mesh%nmxfr
   nuctfr(i)=0
end do

! ======================================================================
! --- 2.0 --- Numerotation des cotes -----------------------------------

do iel=1,mesh%num_cells
   do nc=1,3
      ivois=mesh%nvois(nc,iel)
      n1=nc
      n2=mod(nc,3)+1

      num1 =mesh%nodes(n1,iel)
      num2 =mesh%nodes(n2,iel)
      indn1=mesh%npoel1(num1)
      indn2=mesh%npoel1(num2)
      indv1=vmsh%nbcov(num1)
      indv2=vmsh%nbcov(num2)
      nel1 =mesh%npoel1(num1+1)-mesh%npoel1(num1)
      nel2 =mesh%npoel1(num2+1)-mesh%npoel1(num2)

      !Cas des cotes internes ...........................................

      if(ivois >  iel) then

         nucti=nucti+1

         !Numeros globaux de cotes pointant vers le meme noeud 

         do nm1=1,nel1
            if(mesh%npoel2(indn1+nm1) == ivois) then
               vmsh%nugcv(indv1+nm1)=nucti
            end if
         end do

         do nm2=1,nel2
            if(mesh%npoel2(indn2+nm2) == iel) then
               vmsh%nugcv(indv2+nm2)=nucti
            end if
         end do

         !Numeros des triangles ou polygones associes

         vmsh%nuvac(1,nucti)=num1
         vmsh%nuvac(2,nucti)=num2

         !Cas des cotes frontaliers ........................................

      else if(ivois < 0) then

         ind=-ivois
         nuctfr(ind)=nuctfr(ind)+1
         nuctf=vmsh%ncotcu(ind+1)+nuctfr(ind)

         vmsh%nugcv(indv1+nel1+1)=nuctf
         vmsh%nugcv(indv2+nel2  )=nuctf
         vmsh%nuvac(1,nuctf)=num1
         vmsh%nuvac(2,nuctf)=num2

      end if

   end do

end do

!======================================================================
!----------- Longueurs des cotes des triangles ------------------------

do ic=1,mesh%nbtcot
   xa=mesh%coord(1,vmsh%nuvac(1,ic))
   ya=mesh%coord(2,vmsh%nuvac(1,ic))
   xb=mesh%coord(1,vmsh%nuvac(2,ic))
   yb=mesh%coord(2,vmsh%nuvac(2,ic))
   vmsh%xlcod(ic)=sqrt((xa-xb)*(xa-xb)+(ya-yb)*(ya-yb))
end do

!======================================================================
!--- 4.0 --- Calcul des matrices de lissage ---------------------------

do ic=1,mesh%nbtcot
   n1 = vmsh%nuvac(1,ic)
   n2 = vmsh%nuvac(2,ic)
   x21= (mesh%coord(1,n2)-mesh%coord(1,n1))/vmsh%xlcod(ic)
   y21= (mesh%coord(2,n2)-mesh%coord(2,n1))/vmsh%xlcod(ic)
   s1 = x21*x21
   s2 = y21*y21
   s3 = x21*y21
   vtaux(ic)=x21
   vtauy(ic)=y21
   mesh%xmal1(n1)=mesh%xmal1(n1)+s1
   mesh%xmal2(n1)=mesh%xmal2(n1)+s2
   mesh%xmal3(n1)=mesh%xmal3(n1)+s3
   mesh%xmal1(n2)=mesh%xmal1(n2)+s1
   mesh%xmal2(n2)=mesh%xmal2(n2)+s2
   mesh%xmal3(n2)=mesh%xmal3(n2)+s3
end do

!--- 4.5 --- Normalisation par rapport aux determinants ---------------
lerr = .false.

do is=1,mesh%num_nodes
   det=mesh%xmal1(is)*mesh%xmal2(is)-mesh%xmal3(is)*mesh%xmal3(is)
   if(det /= 0) then
      mesh%xmal1(is)=mesh%xmal1(is)/det
      mesh%xmal2(is)=mesh%xmal2(is)/det
      mesh%xmal3(is)=mesh%xmal3(is)/det
   else
      lerr=.TRUE.
   end if
end do
 
if(lerr) then
   write(iout,900)
   !call errout(iout,"F","poclis","maillage.f90" )
end if

! ======================================================================
! --- 5.0 --- Impression des tableaux ----------------------------------
 
!write(iout,902)
!do is=1,mesh%num_nodes
!   write(iout,903) mesh%xmal1(is),mesh%xmal2(is),mesh%xmal3(is)
!end do

!write(iout,904)
!do ic=1,mesh%nbtcot
!   write(iout,905) vtaux(ic),vtauy(ic)
!end do

if(lerr) then
   write(iout,901)
   !call errout(iout,"F","poclis","maillage.f90" )
end if

 900  format(//10x,'Determinant des coefficients des matrices'  &
              /10x,'de lissage nul'             &
              /10x,'Il faut modifier le maillage'//)
 901  format(//10x,'Le nombre de triangles communs a un noeud'  &
              /10x,'est superieur a 12'             &   
              /10x,'Modifiez legerement le maillage'//)
 902  format(//10x,'Matrices de lissage'            &
              /10x,'xmal1',7x,'xmal2',7x,'xmal3')
 903  format(  10x,3e12.3)
 904  format(//10x,'Composantes des vecteurs tangeants'     &
              /10x,'vtaux',7x,'vtauy')
 905  format(  10x,2e12.3)

end subroutine poclis

!-------------------------------------------------------------------------

!Subroutine: write_mesh
!Ecriture du maillage pour visualisation au format GNUPLOT et PLOTMTV
subroutine write_mesh(mesh)
type(sll_triangular_mesh_2d) :: mesh
double precision :: x1, y1
integer :: iev, is1, is2, iel

!Trace du maillage

!format gnuplot
if (titre == "") titre = "dummy"
open(10,file=trim(titre)//"_mesh.dat")
do i = 1, mesh%num_cells
   write(10,*) mesh%coord(1:2,mesh%nodes(1,i)),0.0
   write(10,*) mesh%coord(1:2,mesh%nodes(2,i)),0.0
   write(10,*) mesh%coord(1:2,mesh%nodes(3,i)),0.0
   write(10,*) mesh%coord(1:2,mesh%nodes(1,i)),0.0
   write(10,*) 
   write(10,*) 
end do
close(10)
!format plotmtv


open( 10, file=trim(titre)//"_mesh.mtv")

!--- Trace du maillage ---

write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Maillage' "

do i = 1, mesh%num_cells
   write(10,*)mesh%coord(1:2,mesh%nodes(1,i)),0.0
   write(10,*)mesh%coord(1:2,mesh%nodes(2,i)),0.0
   write(10,*)mesh%coord(1:2,mesh%nodes(3,i)),0.0
   write(10,*)mesh%coord(1:2,mesh%nodes(1,i)),0.0
   write(10,*)
end do

!--- Numeros des noeuds et des triangles

write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='Numeros des noeuds et des triangles' "

do i = 1, mesh%num_cells
   write(10,*)mesh%coord(1:2,mesh%nodes(1,i)),0.0
   write(10,*)mesh%coord(1:2,mesh%nodes(2,i)),0.0
   write(10,*)mesh%coord(1:2,mesh%nodes(3,i)),0.0
   write(10,*)mesh%coord(1:2,mesh%nodes(1,i)),0.0
   write(10,*)
end do

do i = 1, mesh%num_cells
   x1 = (  mesh%coord(1,mesh%nodes(1,i))  &
         + mesh%coord(1,mesh%nodes(2,i))  &
     + mesh%coord(1,mesh%nodes(3,i))    )/3.
   y1 = (  mesh%coord(2,mesh%nodes(1,i))  &
         + mesh%coord(2,mesh%nodes(2,i))  &
     + mesh%coord(2,mesh%nodes(3,i))    )/3.
   write(10,"(a)"   , advance="no")"@text x1="
   write(10,"(f8.5)", advance="no") x1
   write(10,"(a)"   , advance="no")" y1="
   write(10,"(f8.5)", advance="no") y1
   write(10,"(a)"   , advance="no")" z1=0. lc=4 ll='"
   write(10,"(i4)"  , advance="no") i
   write(10,"(a)")"'"
end do

do i = 1, mesh%num_nodes
   x1 = mesh%coord(1,i)
   y1 = mesh%coord(2,i)
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

do i = 1, mesh%num_cells
   write(10,*)mesh%coord(1:2,mesh%nodes(1,i)),0.
   write(10,*)mesh%coord(1:2,mesh%nodes(2,i)),0.
   write(10,*)mesh%coord(1:2,mesh%nodes(3,i)),0.
   write(10,*)mesh%coord(1:2,mesh%nodes(1,i)),0.
   write(10,*)
end do

do i = 1, mesh%num_nodes
   x1 = mesh%coord(1,i)
   y1 = mesh%coord(2,i)
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

do i = 1, mesh%num_cells
   write(10,*)mesh%coord(1:2,mesh%nodes(1,i)),0.
   write(10,*)mesh%coord(1:2,mesh%nodes(2,i)),0.
   write(10,*)mesh%coord(1:2,mesh%nodes(3,i)),0.
   write(10,*)mesh%coord(1:2,mesh%nodes(1,i)),0.
   write(10,*)
end do

do i = 1, mesh%num_cells
   x1 = (  mesh%coord(1,mesh%nodes(1,i))  &
         + mesh%coord(1,mesh%nodes(2,i))  &
     + mesh%coord(1,mesh%nodes(3,i))    )/3.
   y1 = (  mesh%coord(2,mesh%nodes(1,i))  &
         + mesh%coord(2,mesh%nodes(2,i))  &
     + mesh%coord(2,mesh%nodes(3,i))    )/3.
   write(10,"(a)"   , advance="no")"@text x1="
   write(10,"(g15.3)", advance="no") x1
   write(10,"(a)"   , advance="no")" y1="
   write(10,"(g15.3)", advance="no") y1
   write(10,"(a)"   , advance="no")" z1=0. lc=4 ll='"
   write(10,"(i4)"  , advance="no") i
   write(10,"(a)")"'"
end do

!Reference des triangles

write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='References des triangles' "

do i = 1, mesh%num_cells
   write(10,*)mesh%coord(1:2,mesh%nodes(1,i)),0.
   write(10,*)mesh%coord(1:2,mesh%nodes(2,i)),0.
   write(10,*)mesh%coord(1:2,mesh%nodes(3,i)),0.
   write(10,*)mesh%coord(1:2,mesh%nodes(1,i)),0.
   write(10,*)
end do

do i = 1, mesh%num_cells
   x1 = (  mesh%coord(1,mesh%nodes(1,i))  &
         + mesh%coord(1,mesh%nodes(2,i))  &
     + mesh%coord(1,mesh%nodes(3,i))    )/3.
   y1 = (  mesh%coord(2,mesh%nodes(1,i))  &
         + mesh%coord(2,mesh%nodes(2,i))  &
     + mesh%coord(2,mesh%nodes(3,i))    )/3.
   write(10,"(a)"   , advance="no")"@text x1="
   write(10,"(g15.3)", advance="no") x1
   write(10,"(a)"   , advance="no")" y1="
   write(10,"(g15.3)", advance="no") y1
   write(10,"(a)"   , advance="no")" z1=0. lc=4 ll='"
   write(10,"(i4)"  , advance="no") mesh%reft(i)
   write(10,"(a)")"'"
end do

!Reference des noeuds 
write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='References des noeuds ' "

do i = 1, mesh%num_cells
   write(10,*)mesh%coord(1:2,mesh%nodes(1,i)),0.
   write(10,*)mesh%coord(1:2,mesh%nodes(2,i)),0.
   write(10,*)mesh%coord(1:2,mesh%nodes(3,i)),0.
   write(10,*)mesh%coord(1:2,mesh%nodes(1,i)),0.
   write(10,*)
end do

do i = 1, mesh%num_nodes
   x1 = mesh%coord(1,i)
   y1 = mesh%coord(2,i)
   write(10,"(a)"   , advance="no")"@text x1="
   write(10,"(f8.4)", advance="no") x1
   write(10,"(a)"   , advance="no")" y1="
   write(10,"(f8.4)", advance="no") y1
   write(10,"(a)"   , advance="no")" z1=0. lc="
   write(10,"(i1)"  , advance="no") mesh%refs(i)
   write(10,"(a)"   , advance="no")" ll='"
   write(10,"(i4)"  , advance="no") mesh%refs(i)
   write(10,"(a)")"'"
end do

!Reference des frontieres 
write(10,*)"$DATA=CURVE3D"
write(10,*)"%equalscale=T"
write(10,*)"%toplabel='References des frontieres ' "

do iel = 1, mesh%num_cells      !Boucle sur les elements
   do iev = 1, 3     !Boucle sur les voisins

      if (mesh%nvois(iev,iel)<0) then
         select case(iev)
         case(1)
            is1 = 1; is2 = 2
         case(2)
            is1 = 2; is2 = 3
         case(3)
            is1 = 3; is2 = 1
         end select
         write(10,*)"%linecolor=",-mesh%nvois(iev,iel)
         write(10,*)mesh%coord(1,mesh%nodes(is1,iel)),&
                    mesh%coord(2,mesh%nodes(is1,iel)),0.
         write(10,*)mesh%coord(1,mesh%nodes(is2,iel)), &
                    mesh%coord(2,mesh%nodes(is2,iel)),0.
         write(10,*)
         x1 = 0.5*(  mesh%coord(1,mesh%nodes(is1,iel)) &
                   + mesh%coord(1,mesh%nodes(is2,iel)))
         y1 = 0.5*(  mesh%coord(2,mesh%nodes(is1,iel)) &
                   + mesh%coord(2,mesh%nodes(is2,iel)))
         write(10,"(a)"   ,advance="no")"@text x1="
         write(10,"(g15.3)",advance="no") x1
         write(10,"(a)"   ,advance="no")" y1="
         write(10,"(g15.3)",advance="no") y1
         write(10,"(a)"   ,advance="no")" z1=0. ll="
         write(10,"(i5)") -mesh%nvois(iev,iel)
         write(10,*)
      end if

   end do
end do



write(10,*)"$end"
close(10)

end subroutine write_mesh


end module maillage

