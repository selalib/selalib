!File: Module Poisson
!Solveur de Poisson sur un maillage non structure
!
! Traduit en Fortran 90 a partir de M2V ou DEGAS2D
!
module poisson

use zone, only: iout,                &
        grandx, petitx, time, dt, eps0,     &
        pcharg, mesh_fields, pi, xmu0

use maillage, only: sll_triangular_mesh_2d, voronoi, vtaux, vtauy

use solveurs_module, only: mesh_bound

use utlib, only: utfact

implicit none

integer, public  :: niem0, niemp0

integer, private :: i, j

double precision, dimension(:), allocatable :: vnx, vny
integer, dimension(:), allocatable :: naux
double precision, parameter  :: pivpoi = 1.0d-10
double precision, parameter  :: pivamp = 1.0d-10

!    Caracteristiques du maillage triangulaire:         
!
!        nbs     - nombre de noeuds du maillage         
!        nbt     - nombre de triangles du maillage     
!        ndiric  - nombre de noeuds verifiant Dirichlet              
!        ndirb3  - nombre de noeuds verifiant Dirichlet pour Ampere 
!        niem0   - nombre de noeuds de la frontiere emettrice      
!        niemp0  - nombre de noeuds voisins de la frontiere emettrice
!        nbfrax  - nombre de noeuds a l'intersection axe/frontiere  
!        nmxfr   - nombre de frontieres referencees max            
!        nmxsd   - nombre de sous-domaines references max         
!        nefro   - nombre d'elements ayant au moins 1 sommet frontiere
!        nelfr   - nombre d'elements frontieres                      
!        nelin   - nombre d'elements internes                       
!        nelmatf - nombre d'elements de la matrice profil du Laplacien
!        irefdir - references de frontieres Dirichlet non homogenes  
!        nnoeuf  - nombre de noeuds frontieres Dirichlet non homogenes
!
!    Dimensions de reference:
!
!        xlml   - limite inf x du domaine                          
!        xlmu   - limite sup x du domaine                         
!        ylml   - limite inf y du domaine                        
!        ylmu   - limite sup y du domaine                       
!        petitl - petite longueur de reference                 
!        grandl - grande longueur de reference                
!        imxref - nombre entier eleve de reference           
!                                                           
!    Caracteristiques du maillage de Delaunay-Voronoi:
!
!        nbcoti - nombre de cotes internes Delaunay       
!        nbtcot - nombre total de cotes Delaunay         
!        nbcfli - nombre de cotes internes ne satisfaisant pas CFL
!        ncotcu - nombre de cotes cumules par type de frontiere  
!        nnref  - nombre de references de frontieres Dirichlet non
!                 homogenes                                      
!                                                               
!    Caracteristiques des cotes frontieres:
!
!        nctfrt - nombre total de cotes frontieres            
!        nctfro - nombre de cotes frontieres par reference   
!        nctfrp - pointeur des tableaux de cotes frontieres 
!                                                          
!    Caracteristiques des frontieres internes: 
!
!        nnofnt - nombre de noeuds sur les frontieres internes     
!        ntrfnt - nombre total de triangles (ceux de droite)      
!        ntrfrn - nombre de triangles par frontiere              
!        ntrfrc - nombre cumule de triangles                    
!        nndfnt - noeuds Dirichlet sur les frontieres internes 
!        ncdfnt - cotes  Dirichlet sur les frontieres internes (VF)
!                                                                 
!    Tableau des classes d'elements:
!
!        nmxcol - nombre de classes des elements                
!        nclcol - nombre d'elements dans chaque classe         
!

integer :: ndiric, nnref, ndirb3, nnoeuf

double precision, private, dimension(:), allocatable :: gradx, grady
double precision, private, dimension(:), allocatable :: grgr

integer, dimension(:), allocatable :: iem0, iemp0
integer, dimension(:), allocatable :: mors1, iprof
integer, dimension(:), allocatable :: ifron, ifrb3, irefdir
integer, dimension(:), allocatable :: mors2
double precision,    dimension(:), allocatable :: rho2, dcr1, dcr2
double precision,    dimension(:), allocatable :: amass

integer, dimension(:), allocatable :: inoeuf

logical :: ldebug = .false.

double precision, dimension(:), private, allocatable :: vtantx, vtanty
double precision, dimension(:), private, allocatable :: sv1,   sv2

contains

!Subroutine: init_solveur_poisson   
! Reservation de tableaux et calcul des 
! matrices necessaires au solveur electrostatique   
! et magnetostatique.                               
!
! amass - matrice de masse diagonale     
! mors1 - matrice morse2       
! mors2 - elements des matrices morses
! prof  - matrice profil    
! grgr  - matrice grad-grad
! gradx - matrice gradx   
! grady - matrice grady  
! fron  - noeuds Dirichlet
! xaux  - tableaux auxiliaires reels          
! iaux  - tableaux auxiliaires entiers       
! ggp0  - vecteur  "grad-grad * potentiel" sur la frontiere emissive
! iem0  - tableau des noeuds situes sur la frontiere emissive     
! iemp0 - tableau des noeuds voisins de la frontiere emissive    
! rho2  - tableau densite de charge aux noeuds de la frontiere 
!         emettrice due aux particules crees au cours du pas de tps
! gg0p0 - matrice de l'operateur "grad-grad" limitee aux noeuds de
!         la frontiere emettrice et de leurs voisins             

! Allocation des tableaux permettant de stocker des matrices sous forme morse
! Tableau donnant le numero du dernier terme de chaque ligne (mors1)
subroutine init_solveur_poisson(mesh, bcnd)

type(sll_triangular_mesh_2d),  intent(in) :: mesh
type(mesh_bound), intent(in) :: bcnd
 
double precision, dimension(:), allocatable :: tmp1

integer :: nref, nn, ndir

allocate(sv1(mesh%nbtcot),sv2(mesh%nbtcot)); sv1 = 0.; sv2 = 0.
allocate(vtantx(mesh%nbtcot),vtanty(mesh%nbtcot));vtantx=0.;vtanty=0.

allocate(mors1(mesh%num_nodes+1)); mors1 = 0

! Tableau contenant le numero des termes de chaque ligne (mors2)
! on choisit une taille a priori superieure a la taille reelle
! de ce tableau.

allocate(mors2(12*mesh%num_nodes)); mors2 = 0
 
! Calcul de mors1,mors2.

call morse(mesh%npoel1,mesh%npoel2,mesh%nodes,mesh%num_cells,mesh%num_nodes,mors1,mors2)
 
! Ajustement de la taille de mors2.
! pas sur que ca fonctionne a tous les coups
! deallocate(mors2); allocate(mors2(mors1(mesh%num_nodes+1)))

! Adressage des tableaux permettant de stocker des matrices sous forme profil

allocate(iprof(mesh%num_nodes+1)); iprof = 0
call profil(mesh%nodes,mesh%num_nodes,mesh%npoel1,mesh%npoel2)

!======================================================================
!--- 2.0 --- POISSON par une methode d'elements finis -----------------
!======================================================================
 
!matrice de masse diagonalisee
allocate(amass(mesh%num_nodes)); amass = 0.0 
 
!matrice "grad-grad" stockee sous forme profil.
allocate(grgr(iprof(mesh%num_nodes+1))); grgr = 0.0
 
!gradx et grady 

allocate(gradx(mors1(mesh%num_nodes+1))); gradx = 0.0
allocate(grady(mors1(mesh%num_nodes+1))); grady = 0.0

niem0  = 0
niemp0 = 0

!--- Tableau relatif aux frontieres Dirichlet -------------------------

ndir=0
allocate(ifron(mesh%num_nodes)); ifron = 0

do nn=1,mesh%num_nodes
   nref= mesh%refs(nn)
   if (nref > 0) then 
      if (bcnd%ntypfr(nref)==1 .or. bcnd%ntypfr(nref) == 5)then
         ndir=ndir+1
         ifron(ndir)=nn
      end if 
   end if 
end do

ndiric=ndir
!Ajustement ....
!deallocate(ifron); allocate(ifron(ndiric))

!--- Calcul des matrices ----------------------------------------------

call poismc(mesh,bcnd)

!Calcul de la matrice B tel que B*Bt = A dans le cas Cholesky

allocate(tmp1(iprof(mesh%num_nodes+1))); tmp1 = 0.0

write(*,"(//5x,a)")" *** Appel Choleski pour Poisson ***  "
call choles(iprof,grgr,pivpoi,tmp1)

do i=1,iprof(mesh%num_nodes+1)
   grgr(i)=tmp1(i)
end do

deallocate(tmp1)

! --- 8.5 --- Ecriture des tableaux ------------------------------------

if (ldebug) then
   write(iout,900) 
   do i=1,mesh%num_nodes
      write(iout,901) i,(mors2(j), j=mors1(i)+1,mors2(i+1))
   end do
   write(iout,902) 
   write(iout,903) (ifron(i),i=1,ndiric)
end if

! ======================================================================
! --- 9.0 --- Formats --------------------------------------------------
 
 900 format(//10x,'Tableau MORS2 pointe par le tableau MORS1'/  &
               3x,'No de noeud           Noeuds associes'/)
 901 format(2x,I8,3x,12I8)
 902 format(//10x,'Noeuds frontiere du type DIRICLET pour POISSON'/)
 903 format(32000(2x,7I9/)/)

end subroutine init_solveur_poisson

!======================================================================


!Function: morse
! Calcul du tableau contenant l'adresse du dernier terme de   
! chaque ligne des matrices "morse" et le tableau contenant les
! numeros des termes de ces matrices. Le terme diagonal est    
! range en derniere place de chaque ligne.                      
!          
! npoel1 - emplacement dans npoel2 du dernier element relatif a chaque
!          noeud avec npoel1(1)=0 et npoel1(i+1) relatif au noeud i   
! npoel2 - tableau des numeros des elements ayant un sommet en commun  
! ntri   - numeros des sommets des triangles                 
! nbt    - nombre de triangles du maillage                    
! nbs    - nombre de noeuds                                    
!                                                               
! mors1  - tableau des adresses des derniers termes de chaque ligne 
!          de la matrice avec la convention:                         
!          mors1(1)=0 et mors1(i+1) adresse de aii                    
! mors2  - tableau des numeros des termes des matrices "morse"         
subroutine morse(npoel1, npoel2, ntri, nbt, nbs, mors1, mors2)

integer, intent(in) :: nbs, nbt
integer, dimension(:), intent(in) :: npoel1, npoel2
integer, dimension(3,nbt), intent(in) :: ntri
integer, dimension(:),   intent(out):: mors1, mors2
integer, dimension(20)  :: ilign
integer :: l, itest1, itest2, js1, js2, is1, is2, is3, numel
integer :: iel, nlign, nel, is, im = 0, k = 0

mors1(1)=0
 
do is=1,nbs  !Boucle sur les noeuds

   !Nombre d'elements ayant is comme sommet

   nel=npoel1(is+1)-npoel1(is)
   nlign=0

   !Boucle sur ces elements

   do  iel=1,nel

      k=k+1

      !numero de l'element

      numel=npoel2(k)
      is1=ntri(1,numel); is2=ntri(2,numel); is3=ntri(3,numel)
      if (is1==is) then 
         js1=is2; js2=is3
      end if 
      if (is2==is) then 
         js1=is3; js2=is1
      end if 
      if (is3==is) then 
         js1=is1; js2=is2
      end if 

      !On regarde si les 2 noeuds autres que is de l'element courant ont 
      !deja ete pris en compte dans l'ensemble des noeuds interagissant 
      !avec is.

      itest1=0; itest2=0
      if (nlign.NE.0) then 
         do l=1,nlign
            if (js1==ilign(l)) then 
               itest1=1
            end if 
            if (js2==ilign(l)) then 
               itest2=1
            end if 
         end do
      end if 

      if (itest1==0) then 
         nlign=nlign+1
         ilign(nlign)=js1
      end if 
      if (itest2==0) then 
         nlign=nlign+1
         ilign(nlign)=js2
      end if 

   end do
 
   !Definition de l'adresse du dernier terme de la ligne

   mors1(is+1)=mors1(is)+nlign+1
 
   !Remplissage du tableau mors2 avec les numeros des termes
   !de la ligne is.

   if (nlign.NE.0) then 
      do l=1,nlign
         im=im+1
         mors2(im)=ilign(l)
      end do
   end if 
   im=im+1
   mors2(im)=is
end do

end subroutine morse


! ======================================================================


!Function: poismc
!                                                 
! But: 
!   Calcul de toutes les matrices du systeme   
!       pour POISSON cartesien                      
!                                                    
! Parametres d'entree:                                
!
! m      - super-tableau                                        
! coor   - coordonnees des noeuds                                
! refs   - indice permettant de savoir si un noeud appartient     
!          a une frontiere referencee                              
! ifron  - tableau des noeuds Dirichlet                             
! mors1  - tableau du nombre de termes par ligne de la matrice morse 
!          symetrique                                                 
! mors2  - tableau des numeros des termes de chaque ligne de la matrice
!          morse symetrique                              
! ntri   - numeros des sommets des triangles              
! aire   - aire des triangles                              
! iprof  - profil de la matrice grad-grad                   
!                                                            
! noefnt - noeuds internes Dirichlet                          
!                                                              
! Parametres resultats:
! 
! agrgr  - matrice de l'operateur "grad-grad" sous forme "morse"
! aliss  - matrice de lissage                                   
! amass  - matrice de masse sous forme diagonalisee              
! amclt  - matrice de masse diagonalisee relative a la condition  
!          aux limites absorbante                             
! aroro  - matrice de l'operateur "rot-rot" sous forme "morse" 
! d1dx,  - derives des fonctions de base dans chaque triangle   
!..d3dy                                                          
! gradx  - matrice de l'operateur gradx                           
! grady  - matrice de l'operateur grady                            
! rotx   - matrice de l'operateur rotx                              
! roty   - matrice de l'operateur roty                               
! aire   - surface de chaque element                                  
!Auteur:
!   612-POISMC     
!
! Puertolas - Version 1.0  Octobre  1992  
subroutine poismc(mesh, bcnd)

type(sll_triangular_mesh_2d),  intent(in) :: mesh
type(mesh_bound), intent(in) :: bcnd
double precision, dimension(size(mesh%refs)) :: vectmp
double precision :: amloc(3),aggloc(9),grxloc(9),gryloc(9)
double precision :: dntx1, dntx2, dntx3, dnty1, dnty2, dnty3 
double precision :: x1t, x2t, x3t, y1t, y2t, y3t, coef
integer :: is1t, is2t, is3t, iel, nis
integer :: is, il
 
!Boucle sur les elements.

do iel=1,mesh%num_cells

   !Calcul des coefficients dependant de la geometrie du triangle.

   is1t = mesh%nodes(1,iel)
   is2t = mesh%nodes(2,iel)
   is3t = mesh%nodes(3,iel)

   x1t  = mesh%coord(1,is1t)
   x2t  = mesh%coord(1,is2t)
   x3t  = mesh%coord(1,is3t)

   y1t  = mesh%coord(2,is1t)
   y2t  = mesh%coord(2,is2t)
   y3t  = mesh%coord(2,is3t)

   dntx1 = y2t-y3t
   dntx2 = y3t-y1t
   dntx3 = y1t-y2t

   dnty1 = x3t-x2t
   dnty2 = x1t-x3t
   dnty3 = x2t-x1t

   !Contribution a la matrice de masse

   amloc(1) = mesh%aire(iel)/3.
   amloc(2) = mesh%aire(iel)/3.
   amloc(3) = mesh%aire(iel)/3.

   !Assemblage

   call asbld(amloc,is1t,is2t,is3t,amass)

   !Contribution a la matrice grad-grad

   coef=1./(4.*mesh%aire(iel))

   aggloc(1)=(dntx1**2   +dnty1**2   )*coef
   aggloc(2)=(dntx1*dntx2+dnty1*dnty2)*coef
   aggloc(3)=(dntx1*dntx3+dnty1*dnty3)*coef
   aggloc(4)=(dntx2**2   +dnty2**2   )*coef
   aggloc(5)=(dntx2*dntx3+dnty2*dnty3)*coef
   aggloc(6)=(dntx3**2   +dnty3**2   )*coef

   !Contribution aux matrices gradx et grady:
   !Calcul de matrices locales.

   grxloc(1)=-dntx1/6.; gryloc(1)=-dnty1/6.
   grxloc(2)=-dntx2/6.; gryloc(2)=-dnty2/6.
   grxloc(3)=-dntx3/6.; gryloc(3)=-dnty3/6.
   grxloc(4)=-dntx1/6.; gryloc(4)=-dnty1/6.
   grxloc(5)=-dntx2/6.; gryloc(5)=-dnty2/6.
   grxloc(6)=-dntx3/6.; gryloc(6)=-dnty3/6.
   grxloc(7)=-dntx1/6.; gryloc(7)=-dnty1/6.
   grxloc(8)=-dntx2/6.; gryloc(8)=-dnty2/6.
   grxloc(9)=-dntx3/6.; gryloc(9)=-dnty3/6.

   !Assemblage 

   call asblp(aggloc,is1t,is2t,is3t,grgr)
   call asblm2(grxloc,gryloc,mors1,mors2,is1t,is2t,is3t,gradx,grady)

end do

! ======================================================================
! ... Prise en compte des conditions aux limites Dirichlet

do j=1,ndiric
   is=ifron(j)             
   grgr(iprof(is+1))=grandx
end do

! ======================================================================
! ... Frontieres internes Dirichlet                           

if (bcnd%nbfrnt>0 .and. mesh%nndfnt>0) then 

   do j=1,mesh%nndfnt
      is=mesh%noefnt(j)             
      grgr(iprof(is+1))=grandx
   end do

   ! ======================================================================
   ! ... Stockage d'une matrice extraite de "grad-grad": ggop0
   ! ... Calcul de grgrp0
   ! ... Transformation de "grad-grad"

 
end if 

! ================================================================
! ... Ecriture des matrices mass, grgr, gradx et grady

if (ldebug) then
   write(iout,900) 
   do is=1,mesh%num_nodes
      write(iout,901) is,amass(is)
   end do

   write(iout,907) 
   write(iout,902) 
   do is=1,mesh%num_nodes
      nis=iprof(is+1)-iprof(is)
      write(iout,903) is,(grgr(iprof(is)+il),il=1,nis)
   end do

   write(iout,904) 
   do is=1,mesh%num_nodes
      nis=mors1(is+1)-mors1(is)
      write(iout,903) is,(gradx(mors1(is)+il),il=1,nis)
   end do

   write(iout,905) 
   do is=1,mesh%num_nodes
      nis=mors1(is+1)-mors1(is)
      write(iout,903) is,(grady(mors1(is)+il),il=1,nis)
   end do
end if

! ================================================================

 900 format(//10x,'Matrice de masse'/               &
              10x,'No de noeud     Terme diagonal'/)
 901 format(  10x,I10,E15.4)
 902 format(//10x,'Matrice grad-grad'/              &
               10x,'No de noeud     Termes de la ligne'/)
 903 format(  I10,5E12.3/10(10x,5E12.3/)/)
 904 format(//10x,'Matrice gradx'/              &
              10x,'No de noeud     Termes de la ligne'/)
 905 format(//10x,'Matrice grady'/              &
               10x,'No de noeud     Termes de la ligne'/)
 907 format(//10x,'Resolution du systeme par Choleski')

end subroutine poismc

!Function: poissn
!calculer potentiel et champ electrique de l'equation de
!Poisson en cartesien ou cylindrique
!
!
! Variables en argument :                                
!
! eb%f1  - e1 projete sur les noeuds        
! eb%f2  - e2 projete sur les noeuds       
! rho    - densite de charge aux noeuds du maillage 
! phi    - potentiel aux noeuds du maillage        
! ifron  - tableau des noeuds frontaliers verifiant Dirichlet
! noefnt - noeuds sur les frontieres internes Dirichlet     
! irffnt - numero de reference de ces noeuds               
! grgr   - matrice de l'operateur "grad-grad"             
! amass  - matrice de masse sous forme diagonalisee      
! gra1   - matrice de l'operateur "grad1"               
! gra2   - matrice de l'operateur "grad2"              
! mors1  - tableau descriptif des matrices "morse"  
! mors2  - tableau descriptif des matrices "morse"
! iprof  - matrice profil                           
!
! Tableaux auxilliaires :                            
!
! grgrdd - produit de la matrice grgr par le vecteur 
!          "direction de descente" p (gradient conjugue) 
! dird   - direction de descente dans la methode de gradient conjugue 
! res    - residu dans la methode de gradient conjugue               
! sdmb   - second membre du systeme a resoudre pour traiter         
!          le potentiel                                            
! sdmb12 - second membre du systeme a resoudre pour traiter       
!          le champ electrique                                   
! precd  - tableau local utilise dans la resolution de la correction 
!          par un gradient conjugue preconditionne                  
! nbs    - nombre de noeuds du maillage                            
! nmxfr  - nombre de frontieres referencees                       
! errpoi - precision (max du carre de l'erreur relative)         
! nitgc  - nombre max d'iterations du gradient conjugue         
!                                                              
! Auteur:
!  E. Puertolas - Version 1.0  Decembre 1991  
subroutine poissn(bcnd,ebj,rho,phi,mesh,istep)

type (mesh_bound)  :: bcnd
type (sll_triangular_mesh_2d)   :: mesh
type (mesh_fields) :: ebj
double precision, dimension(:) :: phi, rho
double precision :: factt, tt0, tt1, tt2, tt3
double precision :: sdmb(size(rho)), sdmb12(size(rho))
double precision :: fonc, e0, phas, freq, omeg

integer :: is, istep, nref, iform, iond

! ----------- CALCUL DU POTENTIEL  -------------------------------------
!
! ... Calcul du facteur temporel pour la frontiere emettrice 
!     dans le cas Child Langmuir.

!... Calcul du second membre complet ...

do is=1,mesh%num_nodes
   sdmb(is)=amass(is)*rho(is)/eps0
end do

!... Condition aux limites Dirichlet homogene avec forme temporelle

do is=1,ndiric

   nref=mesh%refs(ifron(is))
   iform=bcnd%ifopot(nref)

   if (iform>0 .and. iform<=10) then 
      tt0=bcnd%tdbpot(nref)
      tt1=tt0+bcnd%tmtpot(nref)
      tt2=tt1+bcnd%tplpot(nref)
      tt3=tt2+bcnd%tdspot(nref)
      factt=utfact(time,dt,tt0,tt1,tt2,tt3,iform)
   else
      factt=1.
   end if 

   iond = bcnd%numond(nref)
   e0   = bcnd%eleond(iond)
   freq = bcnd%freond(iond)
   phas = bcnd%phaond(iond)
   omeg = 2.*pi*freq
   fonc = e0 * cos(omeg*time+phas)

   if ( fonc /= 0 ) factt = factt * fonc

   sdmb(ifron(is))=bcnd%potfr(nref)*factt

end do

!*** Frontiere internes (Noeuds Dirichlet) ***

if (bcnd%nbfrnt>0.and.mesh%nndfnt > 0) then 

   do is=1,mesh%nndfnt

      nref=mesh%irffnt(is)
      iform=bcnd%ifopot(nref)

      if ( (iform>0).and.(iform<=10) ) then 
         tt0=bcnd%tdbpot(nref)
         tt1=tt0+bcnd%tmtpot(nref)
         tt2=tt1+bcnd%tplpot(nref)
         tt3=tt2+bcnd%tdspot(nref)
         factt=utfact(time,dt,tt0,tt1,tt2,tt3,iform)
      else
     factt=1.
      end if 

      sdmb(mesh%noefnt(is))=bcnd%potfr(nref)*factt

   end do

end if 

if (istep==1) then 
   do is=1,niem0
      phi(iem0(is))  = 0.
      sdmb(iem0(is)) = 0.
   end do
end if 

do is=1,ndiric
   sdmb(ifron(is))=sdmb(ifron(is))*grandx
end do

if (bcnd%nbfrnt>0.and.mesh%nndfnt > 0) then 
   do is=1,mesh%nndfnt
      sdmb(mesh%noefnt(is))=sdmb(mesh%noefnt(is))*grandx
   end do
end if 

call desrem(iprof, grgr,sdmb,mesh%num_nodes,phi)

!*** CALCUL DES CHAMPS E1N et E2N:

!*** Second membre pour la composante 1:

call m1p(gradx,mors1,mors2,phi,mesh%num_nodes,sdmb12)

do i=1,mesh%num_nodes
   ebj%e(1,i)=sdmb12(i)/amass(i)
end do

!*** Second membre pour la composante 2:

call m1p(grady,mors1,mors2,phi,mesh%num_nodes,sdmb12)
do i=1,mesh%num_nodes
   ebj%e(2,i)=sdmb12(i)/amass(i)
end do

end subroutine poissn


!Function: poifrc
!Correction des champs sur les frontieres, en particulier
!
! version qui ne tient pas compte de l'ordre des cotes de la frontiere
!                                                                    
!         - E2n=0 sur l'axe                            
!         - E.tau = 0 sur les frontieres Dirichlets
!         - E.nu =  0 sur les frontieres Neumann  
!                                                                
!Variables en argument:                                        
!
!            ksofro  - numeros des 2 sommets extremite du cote 
!            krefro  - numero de reference du cote            
!            vnofro  - composantes du vecteur normal (vers l'interieur)
!        vnx,vny - tableaux temporaires donnant les composantes d'une normale
!                      aux noeuds                                                
!                                                                   
!Tableaux auxilliaires:                                           
!
! naux   - tableau auxiliaire permettant de reperer les noeuds d'une
!          frontiere                                                
!Auteurs:
!  7711-POifRC   
!
! A. Adolf - Version 1.0  Octobre  1994
! J. Segre - Version 1.1  Avril    1998 
subroutine poifrc(ebj,mesh, bcnd)

type(sll_triangular_mesh_2d)   :: mesh
type(mesh_bound)  :: bcnd
type(mesh_fields) :: ebj
double precision :: pscal, xnor
integer :: is1, is2, ict
      
!!$ ======================================================================
!!$ ... On force E.tau = 0 sur toutes les frontieres Dirichlet 

!!$ ... Initialisation des normales aux noeuds et du tableau indiquant 
!!$ ... les noeuds appartenant a la frontiere consideree

if (.not. allocated(vnx)) then 
   allocate(vnx(mesh%num_nodes), vny(mesh%num_nodes), naux(mesh%num_nodes))
end if 

vnx=0.; vny=0.; naux=0

!!$ ... Boucle sur les cotes frontieres pour construire les normales aux
!!$ ... noeuds "Dirichlet"

do  ict=1,mesh%nctfrt

   if (bcnd%ntypfr(mesh%krefro(ict))==1 .or.  bcnd%ntypfr(mesh%krefro(ict))==5) then 

     is1=mesh%ksofro(1,ict)
     is2=mesh%ksofro(2,ict)

     vnx(is1)=vnx(is1)+mesh%vnofro(1,ict)
     vny(is1)=vny(is1)+mesh%vnofro(2,ict)
     vnx(is2)=vnx(is2)+mesh%vnofro(1,ict)
     vny(is2)=vny(is2)+mesh%vnofro(2,ict)

     naux(is1)=1
     naux(is2)=1

  end if 

end do 

!... on impose la condition E.tau=0

do  i=1,mesh%num_nodes

   if (naux(i)==1) then 

      xnor=SQRT(vnx(i)**2+vny(i)**2)
      if (xnor>mesh%petitl) then 
         vnx(i)=vnx(i)/xnor
         vny(i)=vny(i)/xnor

         pscal=vnx(i)*ebj%e(1,i)+vny(i)*ebj%e(2,i)
         ebj%e(1,i)=vnx(i)*pscal
         ebj%e(2,i)=vny(i)*pscal
      end if 

   end if 

end do 

!======================================================================
! ... On force E.nu = 0 sur toutes les frontieres Neumann
!
! ... Initialisation des normales aux noeuds et du tableau indiquant 
! ... les noeuds appartenant a la frontiere consideree
!
vnx=0.; vny=0.; naux=0

! ... Boucle sur les cotes frontieres pour construire les normales aux
! ... noeuds "Neumann"

do  ict=1,mesh%nctfrt

   if (bcnd%ntypfr(mesh%krefro(ict))==3 .or. bcnd%ntypfr(mesh%krefro(ict))==6) then 

     is1=mesh%ksofro(1,ict)
     is2=mesh%ksofro(2,ict)

     vnx(is1)=vnx(is1)+mesh%vnofro(1,ict)
     vny(is1)=vny(is1)+mesh%vnofro(2,ict)
     vnx(is2)=vnx(is2)+mesh%vnofro(1,ict)
     vny(is2)=vny(is2)+mesh%vnofro(2,ict)

     naux(is1)=1
     naux(is2)=1

   end if 

end do 

!!... on impose la condition E.nu=0

do i=1,mesh%num_nodes

   if (naux(i)==1) then 

      if (xnor>mesh%petitl) then 
        xnor=SQRT(vnx(i)**2+vny(i)**2)
        vnx(i)=vnx(i)/xnor
        vny(i)=vny(i)/xnor

        pscal=vnx(i)*ebj%e(1,i)+vny(i)*ebj%e(2,i)
        ebj%e(1,i)=ebj%e(1,i)-vnx(i)*pscal
        ebj%e(2,i)=ebj%e(2,i)-vny(i)*pscal
      end if 

   end if 

end do 

end subroutine poifrc

!Function: choles
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
subroutine choles(mudl,ae,eps,as)

!**********************************************************************

integer, dimension(:), intent(in) :: mudl
double precision, dimension(:), intent(in)  :: ae
double precision, dimension(:), intent(inout) :: as
integer :: kj, jid, jmi, ij, jj, id, imi, ii, ntest
double precision    :: s, xii, eps

ntest = 0
as(1) = sqrt(ae(1))
ii    = mudl(1)

!--- 1.0 --- Matrice symetrique ---------------------------------------

do  i=2,size(mudl)

   xii = 0
   imi = ii+1
   ii  = mudl(i)
   id  = i - ii
   jj  = mudl(imi+id-1)

   do  ij =imi,ii-1
      j   = ij+id
      jmi = jj+1
      jj  = mudl(j)
      jid = j - jj -id
      s = 0

      do  kj = max( jmi , imi-jid )  ,  jj-1
         s = s + as( kj ) * as( kj + jid )
      end do

      as(ij) = ( ae(ij) - s ) / as(jj)
      xii   = xii + as(ij)*as(ij)
   end do 

   xii = ae(ii)  - xii
   if ( xii  <  eps*abs(ae(ii))) then
      write(iout,900) i,eps
      write(iout,901)xii,eps,ae(ii)
      ntest = 1
   end if

   as(ii) = sqrt ( xii )

end do

if(ntest==1) then 
   write(iout,902) 
   call errout(iout,"F","choles","poisson.f90")
end if

!--- 9.0 --- Formats ---------------------------------------------------

 900 format(/10x,'resultats a verifier : le coefficient diagonal du dl'  &
             ,i6,' est inferieur au seuil de precision',e14.7)
 901 format(/10x,'xii:',e14.7,e14.7,e14.7)
 902 format(/10x,'************  Erreur dans CHOLES *************'/)

end subroutine choles


!Function: desrem
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
subroutine desrem(mudl,a,be,ntdl,bs)

integer :: ntdl
integer :: mudl(0:*)
double precision    :: a(*),be(*),bs(*)
integer :: ii, ij, kj, il
double precision :: y

ii = mudl(1)

!**********************************************************************
! matrice symetrique
!**********************************************************************
! descentes

do i=1,ntdl
   ij = ii + 1
   ii = mudl(i)
   kj = i  - ii
   y  = 0
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


end subroutine desrem

!Function: profil
!determination des numeros des elements situes sur la diagonale
!de la matrice profil associee a la methode d'elements finis.
!
!Parametres d'entree:
! npoel1 - emplacement dans npoel2 du dernier element relatif a chaque
!          noeud avec npoel1(1)=0 et npoel1(i+1) relatif au noeud i 
! npoel2 - tableau des numeros des elements ayant un sommet en commun
! ntri   - numeros des sommets des triangles                          
! nbt    - nombre de triangles         
! noe    - nombre de noeuds             
!                                                                    
!Parametre resultat:                                              
! iprof  - tableau des numeros des termes situes sur la diagonale avec
!        - iprof(i+1) numero du terme de la diagonale de la ligne i et
!        - iprof(1)=0 
!
!Auteur:
! J. Segre - Juillet 89
subroutine profil(ntri,nbs, npoel1, npoel2)

integer, dimension(:), intent(in) :: npoel1, npoel2
integer, dimension(:,:), intent(in) :: ntri
integer :: in, k, is1, is2, is3, numel, ind, iel, nel
integer, intent(in) :: nbs      !Nombre de sommets

!************************************************************************
!*
!* determination de la longueur de chaque ligne de la matrice profil :
!* pour chaque noeud i on cherche les noeuds j qui lui sont lies;
!* si le noeud j n'a jamais ete rencontre, la longueur de la ligne j est
!* j-i+1
!
!* boucle sur les noeuds

k=0
do in=1,nbs
      
   !* nombre d'elements ayant in comme sommet
   nel=npoel1(in+1)-npoel1(in)

   !* boucle sur ces elements

   do iel=1,nel

      k=k+1
      !* numero de l'element
      numel=npoel2(k)
      is1=ntri(1,numel)
      ind=is1+1
      if (iprof(ind).eq.0)then
         iprof(ind)=is1-in+1
      end if 
      is2=ntri(2,numel)
      ind=is2+1
      if (iprof(ind).eq.0)then
         iprof(ind)=is2-in+1
      end if 
      is3=ntri(3,numel)
      ind=is3+1
      if (iprof(ind).eq.0)then
         iprof(ind)=is3-in+1
      end if 

   end do

end do

!* determination de la position des termes diagonaux de la matrice
!* profil (il suffit de sommer les nombres d'elements des lignes
!* precedentes et de la ligne courante).

do ind=3,nbs+1
   iprof(ind)=iprof(ind-1)+iprof(ind)
end do

end subroutine profil

!Function: asbld
!     Assembler une matrice elementaire    
!     dans une matrice globale dans          
!     le cas ou elle est diagonale            
!                                                                 
!Parametres d'entree: 
!                                 
!          aele         -    Matrice elementaire diagonale           
!                           (3 termes pour un element triangulaire)
!          i1,i2,i3     -    numeros des sommets de l'element      
!                                                                 
!Parametre resultat:                                     
!                                                              
!          xmass    -   matrice globale diagonalisee          
!
!Auteur:
!      J. Segre - Version 1.0  Juillet 1989
subroutine asbld(aele,i1,i2,i3,xmass)
 
double precision, dimension(:) :: aele, xmass
integer :: i1, i2, i3

xmass(i1)=xmass(i1)+aele(1)
xmass(i2)=xmass(i2)+aele(2)
xmass(i3)=xmass(i3)+aele(3)

end subroutine


 
!Function: asbldc
!Assembler une matrice de masse elementaire 
!relative a une condition limite dans une  
!matrice globale dans le cas ou elle est diagonale
!                                                               
!  Entrees:
!  aele         -    matrice elementaire diagonale  
!                   (3 ou 4 termes selon que l'element  
!                    est un triangle ou un quadrilatere)
!  indl         -    tableau des numeros des noeuds    
!                    concernes par la condition limite
!  nbl          -    nombre de noeuds concernes par la
!                    condition limite                
!  i1,i2,i3     -    numeros des sommets de l'element
!                                                          
!  Sorties:
!  xmass    -   matrice globale diagonalisee    
!
!  Auteur:
!  J. Segre - Version 1.0  Juillet 1989
!                                                      
subroutine asbldc(aele,indl,i1,i2,i3,nbl,xmass)

integer :: i1, i2, i3, nbl
double precision,    dimension(:) :: aele, xmass
integer, dimension(:) :: indl
 
! ----------------------------------------------------------------------

do i=1,nbl
   if (indl(i)==i1) then 
      xmass(i)=aele(1)
   end if 
   if (indl(i)==i2) then 
      xmass(i)=aele(2)
   end if 
   if (indl(i)==i3) then 
      xmass(i)=aele(3)
   end if 
end do

end subroutine


!Function: asblm2
!   
!          assembler 3 matrices elementaires
!          dans 3 matrices globales stockees
!          sous forme "morse" non symetriques
!                                                          
!      Entrees:
!          aele1,aele2    -  Matrices elementaires          
!          mors1          -  Tableau du nombre de termes par 
!                            ligne de la matrice morse        
!          mors2          -  Tableau des numeros des termes    
!                            de chaque ligne de la matrice morse
!          i1,i2,i3       -  Numeros des sommets de l'element
!                                                         
!      Sorties:
!          a1,a2          -  Matrices globales stockees     
!                            sous forme "morse"            
!
!Auteur:
!J. Segre - Version 1.0  Aout 1989  
subroutine asblm2(aele1,aele2,mors1,mors2,i1,i2,i3,a1,a2)

integer, intent(in) :: i1, i2, i3
integer, dimension(:), intent(in) :: mors1, mors2
double precision,    dimension(:), intent(in)  :: aele1, aele2
double precision,    dimension(:), intent(out) :: a1, a2
integer :: j1, j2, ind1, ind2, ind3

! --- 1.1 --- Rangement des termes diagonaux ---------------------------

ind1=mors1(i1+1)
a1(ind1)=a1(ind1)+aele1(1)
a2(ind1)=a2(ind1)+aele2(1)

ind2=mors1(i2+1)
a1(ind2)=a1(ind2)+aele1(5)
a2(ind2)=a2(ind2)+aele2(5)

ind3=mors1(i3+1)
a1(ind3)=a1(ind3)+aele1(9)
a2(ind3)=a2(ind3)+aele2(9)
      
! --- 1.2 --- Rangement des autres termes ------------------------------

j2=ind1-1
j1=mors1(i1)+1
 
if (j2>=j1) then 
   do  j=j1,j2
      if (i2==mors2(j)) then 
         a1(j)=a1(j)+aele1(2)
         a2(j)=a2(j)+aele2(2)
      end if 
      if (i3==mors2(j)) then 
         a1(j)=a1(j)+aele1(3)
         a2(j)=a2(j)+aele2(3)
      end if 
   end do
end if 
 
j2=ind2-1
j1=mors1(i2)+1

if (j2>=j1) then 
   do  j=j1,j2
      if (i1==mors2(j)) then 
         a1(j)=a1(j)+aele1(4)
         a2(j)=a2(j)+aele2(4)
      end if 
      if (i3==mors2(j)) then 
         a1(j)=a1(j)+aele1(6)
         a2(j)=a2(j)+aele2(6)
      end if 
   end do
end if 

j2=ind3-1
j1=mors1(i3)+1

if (j2>=j1) then 
   do  j=j1,j2
      if (i1==mors2(j)) then 
        a1(j)=a1(j)+aele1(7)
        a2(j)=a2(j)+aele2(7)
      end if 
      if (i2==mors2(j)) then 
         a1(j)=a1(j)+aele1(8)
         a2(j)=a2(j)+aele2(8)
      end if 
   end do
end if 

end subroutine

!----------------------------------------------------------------------


!Function: asblp
!                                                      
!      assembler une matrice elementaire             
!      dans une matrice globale symetrique          
!      stockee sous forme "profil"                 
!                                                 
!Entrees:
!                                                                 
!          aele         -  Matrice elementaire (6 ou 10 termes   
!                          selon que l'element est un triangle  
!                          ou un quadrilatere)                 
!          iprof        -  Description du profil de la matrice 
!                          c'est-a-dire liste des numeros d'ordre 
!                          des termes diagonaux dans le tableau de
!                          stockage de la matrice                
!          i1,i2,i3     -  Numeros des sommets de l'element     
!                                                              
!Sorties:
!                                                             
!          xmass        -  matrice globale stockee sous forme "profil"
!                                                                    
!Auteur:
!       J. Segre - Version 1.0  Aout 1989
subroutine asblp(aele,i1,i2,i3,xmass)

double precision,    dimension(:) :: aele(*), xmass(*)
integer :: i1, i2, i3, idiag1, idiag2, idiag3, ind

!--- 1.1 --- Rangement des termes diagonaux ---------------------------
 
idiag1 = iprof(i1+1)
idiag2 = iprof(i2+1)
idiag3 = iprof(i3+1)

xmass(idiag1) = xmass(idiag1) + aele(1)
xmass(idiag2) = xmass(idiag2) + aele(4)
xmass(idiag3) = xmass(idiag3) + aele(6)

!--- 1.2 --- Rangement des autres termes ------------------------------
!           (on ne stocke que les aij tels que i>j)
 
if (i1 < i2) then 
   ind=idiag2+i1-i2
else
   ind=idiag1+i2-i1
end if 
 
xmass(ind)=xmass(ind)+aele(2)

if (i1 < i3) then 
   ind=idiag3+i1-i3
else
   ind=idiag1+i3-i1
end if 
 
xmass(ind)=xmass(ind)+aele(3)

if (i2 < i3) then 
   ind=idiag3+i2-i3
else
   ind=idiag2+i3-i2
end if 

xmass(ind)=xmass(ind)+aele(5)

!----------------------------------------------------------------------

end subroutine asblp

 
!Function: m1p
!                                             
!     faire l'operation : yvect =  xmors.xvect 
!     ou  "xmors" est une matrice "morse"       
!     xvect,yvect  des vecteurs                  
!
!Entrees:
!          xmors        -     matrice  "morse"          
!          xvect        -     vecteur operande           
!          mors1,mors2  -     tableaux descriptifs de     
!                             la matrice "morse"           
!          nlign        -     nombre de lignes des matrices 
!                                                      
!Sorties:
!          yvect    -   vecteur resultat            
!                                                  
!Auteur:
!     J. Segre - Version 1.0  Decembre 1989
subroutine m1p(xmors,mors1,mors2,xvect,nlign,yvect)

double precision    :: xmors(*),xvect(*),yvect(*)
integer :: mors1(*),mors2(*)
integer :: il, nlign, noeui

do il=1,nlign

   noeui=mors1(il+1)-mors1(il)
  
   select case (noeui)
   case(6)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+   &
                xmors(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+   &
                xmors(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))
  
   case(5)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+   &
                xmors(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))

   case(7)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+   &
                xmors(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+   &
                xmors(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+   &
                xmors(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))

   case(4)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))
 
   case(8)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+   &
                xmors(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+   &
                xmors(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+   &
                xmors(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))+   &
                xmors(mors1(il+1)-7)*xvect(mors2(mors1(il+1)-7))

   case(3)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))

   case(9)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+   &
                xmors(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+   &
                xmors(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+   &
                xmors(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))+   &
                xmors(mors1(il+1)-7)*xvect(mors2(mors1(il+1)-7))+   &
                xmors(mors1(il+1)-8)*xvect(mors2(mors1(il+1)-8))
 
   case(10)
      yvect(il)=xmors(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+   &
                xmors(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+   &
                xmors(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+   &
                xmors(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+   &
                xmors(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+   &
                xmors(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+   &
                xmors(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))+   &
                xmors(mors1(il+1)-7)*xvect(mors2(mors1(il+1)-7))+   &
                xmors(mors1(il+1)-8)*xvect(mors2(mors1(il+1)-8))+   &
                xmors(mors1(il+1)-9)*xvect(mors2(mors1(il+1)-9))
  
   case(11)
      yvect(il)=xmors(mors1(il+1)   )*xvect(mors2(mors1(il+1)   ))+ &
                xmors(mors1(il+1)- 1)*xvect(mors2(mors1(il+1)- 1))+ &
                xmors(mors1(il+1)- 2)*xvect(mors2(mors1(il+1)- 2))+ &
                xmors(mors1(il+1)- 3)*xvect(mors2(mors1(il+1)- 3))+ &
                xmors(mors1(il+1)- 4)*xvect(mors2(mors1(il+1)- 4))+ &
                xmors(mors1(il+1)- 5)*xvect(mors2(mors1(il+1)- 5))+ &
                xmors(mors1(il+1)- 6)*xvect(mors2(mors1(il+1)- 6))+ &
                xmors(mors1(il+1)- 7)*xvect(mors2(mors1(il+1)- 7))+ &
                xmors(mors1(il+1)- 8)*xvect(mors2(mors1(il+1)- 8))+ &
                xmors(mors1(il+1)- 9)*xvect(mors2(mors1(il+1)- 9))+ &
                xmors(mors1(il+1)-10)*xvect(mors2(mors1(il+1)-10))

   case(12)
      yvect(il)=xmors(mors1(il+1)   )*xvect(mors2(mors1(il+1)   ))+ &
                xmors(mors1(il+1)- 1)*xvect(mors2(mors1(il+1)- 1))+ &
                xmors(mors1(il+1)- 2)*xvect(mors2(mors1(il+1)- 2))+ &
                xmors(mors1(il+1)- 3)*xvect(mors2(mors1(il+1)- 3))+ &
                xmors(mors1(il+1)- 4)*xvect(mors2(mors1(il+1)- 4))+ &
                xmors(mors1(il+1)- 5)*xvect(mors2(mors1(il+1)- 5))+ &
                xmors(mors1(il+1)- 6)*xvect(mors2(mors1(il+1)- 6))+ &
                xmors(mors1(il+1)- 7)*xvect(mors2(mors1(il+1)- 7))+ &
                xmors(mors1(il+1)- 8)*xvect(mors2(mors1(il+1)- 8))+ &
                xmors(mors1(il+1)- 9)*xvect(mors2(mors1(il+1)- 9))+ &
                xmors(mors1(il+1)-10)*xvect(mors2(mors1(il+1)-10))+ &
                xmors(mors1(il+1)-11)*xvect(mors2(mors1(il+1)-11))

   end select 

end do

end subroutine m1p


!subroutine: m2p1a
!
! Faire l'operation: zvect=xmors1.xvect+xmors2.yvect    
! ou  "xmors1" et "xmors2" sont deux matrices "morse"  
! xvect,yvect,zvect des vecteurs.                     
!
!          PARAMETRES D'ENTREE :                    
!
!          xmors1,xmors2 -    deux matrices  "morse" 
!          xvect ,yvect  -    deux vecteurs         
!          mors1 ,mors2  -    tableaux descriptifs des matrices "morse"
!          nlign         -    nombre de lignes des matrices
!
!          PARAMETRE RESULTAT :                         
!
!          zvect    -   vecteur resultat              
! 
!Auteur:
!J. Segre - Version 1.0  Aout 1989 
!                                                                      
subroutine m2p1a(xmors1,xmors2,mors1,mors2,xvect,yvect,nlign,zvect)

integer :: mors1(*),mors2(*)
double precision :: xmors1(*),xmors2(*),xvect (*),yvect(*),zvect(*)

integer :: il, nlign, noeui
 
! ----------------------------------------------------------------------

do il=1,nlign

   noeui=mors1(il+1)-mors1(il)

   if (noeui == 6) then
      zvect(il)=xmors1(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+  &
                xmors1(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+  &
                xmors1(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+  &
                xmors1(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+  &
                xmors1(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+  &
                xmors1(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+  &
                xmors2(mors1(il+1)  )*yvect(mors2(mors1(il+1)  ))+  &
                xmors2(mors1(il+1)-1)*yvect(mors2(mors1(il+1)-1))+  &
                xmors2(mors1(il+1)-2)*yvect(mors2(mors1(il+1)-2))+  &
                xmors2(mors1(il+1)-3)*yvect(mors2(mors1(il+1)-3))+  &
                xmors2(mors1(il+1)-4)*yvect(mors2(mors1(il+1)-4))+  &
                xmors2(mors1(il+1)-5)*yvect(mors2(mors1(il+1)-5))

   else if (noeui == 5) then
      zvect(il)=xmors1(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+  &
                xmors1(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+  &
                xmors1(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+  &
                xmors1(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+  &
                xmors1(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+  &
                xmors2(mors1(il+1)  )*yvect(mors2(mors1(il+1)  ))+  &
                xmors2(mors1(il+1)-1)*yvect(mors2(mors1(il+1)-1))+  &
                xmors2(mors1(il+1)-2)*yvect(mors2(mors1(il+1)-2))+  &
                xmors2(mors1(il+1)-3)*yvect(mors2(mors1(il+1)-3))+  &
                xmors2(mors1(il+1)-4)*yvect(mors2(mors1(il+1)-4))

   else if (noeui == 7) then
      zvect(il)=xmors1(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+  &
                xmors1(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+  &
                xmors1(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+  &
                xmors1(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+  &
                xmors1(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+  &
                xmors1(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+  &
                xmors1(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))+  &
                xmors2(mors1(il+1)  )*yvect(mors2(mors1(il+1)  ))+  &
                xmors2(mors1(il+1)-1)*yvect(mors2(mors1(il+1)-1))+  &
                xmors2(mors1(il+1)-2)*yvect(mors2(mors1(il+1)-2))+  &
                xmors2(mors1(il+1)-3)*yvect(mors2(mors1(il+1)-3))+  &
                xmors2(mors1(il+1)-4)*yvect(mors2(mors1(il+1)-4))+  &
                xmors2(mors1(il+1)-5)*yvect(mors2(mors1(il+1)-5))+  &
                xmors2(mors1(il+1)-6)*yvect(mors2(mors1(il+1)-6))
 
   else if (noeui == 4) then
      zvect(il)=xmors1(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+  &
                xmors1(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+  &
                xmors1(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+  &
                xmors1(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+  &
                xmors2(mors1(il+1)  )*yvect(mors2(mors1(il+1)  ))+  &
                xmors2(mors1(il+1)-1)*yvect(mors2(mors1(il+1)-1))+  &
                xmors2(mors1(il+1)-2)*yvect(mors2(mors1(il+1)-2))+  &
                xmors2(mors1(il+1)-3)*yvect(mors2(mors1(il+1)-3))

   else if (noeui == 8) then
      zvect(il)=xmors1(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+  &
                xmors1(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+  &
                xmors1(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+  &
                xmors1(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+  &
                xmors1(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+  &
                xmors1(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+  &
                xmors1(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))+  &
                xmors1(mors1(il+1)-7)*xvect(mors2(mors1(il+1)-7))+  &
                xmors2(mors1(il+1)  )*yvect(mors2(mors1(il+1)  ))+  &
                xmors2(mors1(il+1)-1)*yvect(mors2(mors1(il+1)-1))+  &
                xmors2(mors1(il+1)-2)*yvect(mors2(mors1(il+1)-2))+  &
                xmors2(mors1(il+1)-3)*yvect(mors2(mors1(il+1)-3))+  &
                xmors2(mors1(il+1)-4)*yvect(mors2(mors1(il+1)-4))+  &
                xmors2(mors1(il+1)-5)*yvect(mors2(mors1(il+1)-5))+  &
                xmors2(mors1(il+1)-6)*yvect(mors2(mors1(il+1)-6))+  &
                xmors2(mors1(il+1)-7)*yvect(mors2(mors1(il+1)-7))

   else if (noeui == 3) then
      zvect(il)=xmors1(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+  &
                xmors1(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+  &
                xmors1(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+  &
                xmors2(mors1(il+1)  )*yvect(mors2(mors1(il+1)  ))+  &
                xmors2(mors1(il+1)-1)*yvect(mors2(mors1(il+1)-1))+  &
                xmors2(mors1(il+1)-2)*yvect(mors2(mors1(il+1)-2))

   else if (noeui == 9) then
      zvect(il)=xmors1(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+  &
                xmors1(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+  &
                xmors1(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+  &
                xmors1(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+  &
                xmors1(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+  &
                xmors1(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+  &
                xmors1(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))+  &
                xmors1(mors1(il+1)-7)*xvect(mors2(mors1(il+1)-7))+  &
                xmors1(mors1(il+1)-8)*xvect(mors2(mors1(il+1)-8))+  &
                xmors2(mors1(il+1)  )*yvect(mors2(mors1(il+1)  ))+  &
                xmors2(mors1(il+1)-1)*yvect(mors2(mors1(il+1)-1))+  &
                xmors2(mors1(il+1)-2)*yvect(mors2(mors1(il+1)-2))+  &
                xmors2(mors1(il+1)-3)*yvect(mors2(mors1(il+1)-3))+  &
                xmors2(mors1(il+1)-4)*yvect(mors2(mors1(il+1)-4))+  &
                xmors2(mors1(il+1)-5)*yvect(mors2(mors1(il+1)-5))+  &
                xmors2(mors1(il+1)-6)*yvect(mors2(mors1(il+1)-6))+  &
                xmors2(mors1(il+1)-7)*yvect(mors2(mors1(il+1)-7))+  &
                xmors2(mors1(il+1)-8)*yvect(mors2(mors1(il+1)-8))

   else if (noeui == 10) then
      zvect(il)=xmors1(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+  &
                xmors1(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+  &
                xmors1(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+  &
                xmors1(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+  &
                xmors1(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+  &
                xmors1(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+  &
                xmors1(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))+  &
                xmors1(mors1(il+1)-7)*xvect(mors2(mors1(il+1)-7))+  &
                xmors1(mors1(il+1)-8)*xvect(mors2(mors1(il+1)-8))+  &
                xmors1(mors1(il+1)-9)*xvect(mors2(mors1(il+1)-9))+  &
                xmors2(mors1(il+1)  )*yvect(mors2(mors1(il+1)  ))+  &
                xmors2(mors1(il+1)-1)*yvect(mors2(mors1(il+1)-1))+  &
                xmors2(mors1(il+1)-2)*yvect(mors2(mors1(il+1)-2))+  &
                xmors2(mors1(il+1)-3)*yvect(mors2(mors1(il+1)-3))+  &
                xmors2(mors1(il+1)-4)*yvect(mors2(mors1(il+1)-4))+  &
                xmors2(mors1(il+1)-5)*yvect(mors2(mors1(il+1)-5))+  &
                xmors2(mors1(il+1)-6)*yvect(mors2(mors1(il+1)-6))+  &
                xmors2(mors1(il+1)-7)*yvect(mors2(mors1(il+1)-7))+  &
                xmors2(mors1(il+1)-8)*yvect(mors2(mors1(il+1)-8))+  &
                xmors2(mors1(il+1)-9)*yvect(mors2(mors1(il+1)-9))

   else if (noeui == 11) then
      zvect(il)=xmors1(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+  &
              xmors1(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+    &
              xmors1(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+    &
              xmors1(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+    &
              xmors1(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+    &
              xmors1(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+    &
              xmors1(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))+    &
              xmors1(mors1(il+1)-7)*xvect(mors2(mors1(il+1)-7))+    &
              xmors1(mors1(il+1)-8)*xvect(mors2(mors1(il+1)-8))+    &
              xmors1(mors1(il+1)-9)*xvect(mors2(mors1(il+1)-9))+    &
              xmors1(mors1(il+1)-10)*xvect(mors2(mors1(il+1)-10))
     zvect(il)=zvect(il)+                       &
              xmors2(mors1(il+1)  )*yvect(mors2(mors1(il+1)  ))+    &
              xmors2(mors1(il+1)-1)*yvect(mors2(mors1(il+1)-1))+    &
              xmors2(mors1(il+1)-2)*yvect(mors2(mors1(il+1)-2))+    &
              xmors2(mors1(il+1)-3)*yvect(mors2(mors1(il+1)-3))+    &
              xmors2(mors1(il+1)-4)*yvect(mors2(mors1(il+1)-4))+    &
              xmors2(mors1(il+1)-5)*yvect(mors2(mors1(il+1)-5))+    &
              xmors2(mors1(il+1)-6)*yvect(mors2(mors1(il+1)-6))+    &
              xmors2(mors1(il+1)-7)*yvect(mors2(mors1(il+1)-7))+    &
              xmors2(mors1(il+1)-8)*yvect(mors2(mors1(il+1)-8))+    &
              xmors2(mors1(il+1)-9)*yvect(mors2(mors1(il+1)-9))+    &
              xmors2(mors1(il+1)-10)*yvect(mors2(mors1(il+1)-10))

   else if (noeui == 12) then
      zvect(il)=xmors1(mors1(il+1)  )*xvect(mors2(mors1(il+1)  ))+  &
              xmors1(mors1(il+1)-1)*xvect(mors2(mors1(il+1)-1))+    &
              xmors1(mors1(il+1)-2)*xvect(mors2(mors1(il+1)-2))+    &
              xmors1(mors1(il+1)-3)*xvect(mors2(mors1(il+1)-3))+    &
              xmors1(mors1(il+1)-4)*xvect(mors2(mors1(il+1)-4))+    &
              xmors1(mors1(il+1)-5)*xvect(mors2(mors1(il+1)-5))+    &
              xmors1(mors1(il+1)-6)*xvect(mors2(mors1(il+1)-6))+    &
              xmors1(mors1(il+1)-7)*xvect(mors2(mors1(il+1)-7))+    &
              xmors1(mors1(il+1)-8)*xvect(mors2(mors1(il+1)-8))+    &
              xmors1(mors1(il+1)-9)*xvect(mors2(mors1(il+1)-9))+    &
              xmors1(mors1(il+1)-10)*xvect(mors2(mors1(il+1)-10))+  &
              xmors1(mors1(il+1)-11)*xvect(mors2(mors1(il+1)-11))
     zvect(il)=zvect(il)+                       &
              xmors2(mors1(il+1)  )*yvect(mors2(mors1(il+1)  ))+    &
              xmors2(mors1(il+1)-1)*yvect(mors2(mors1(il+1)-1))+    &
              xmors2(mors1(il+1)-2)*yvect(mors2(mors1(il+1)-2))+    &
              xmors2(mors1(il+1)-3)*yvect(mors2(mors1(il+1)-3))+    &
              xmors2(mors1(il+1)-4)*yvect(mors2(mors1(il+1)-4))+    &
              xmors2(mors1(il+1)-5)*yvect(mors2(mors1(il+1)-5))+    &
              xmors2(mors1(il+1)-6)*yvect(mors2(mors1(il+1)-6))+    &
              xmors2(mors1(il+1)-7)*yvect(mors2(mors1(il+1)-7))+    &
              xmors2(mors1(il+1)-8)*yvect(mors2(mors1(il+1)-8))+    &
              xmors2(mors1(il+1)-9)*yvect(mors2(mors1(il+1)-9))+    &
              xmors2(mors1(il+1)-10)*yvect(mors2(mors1(il+1)-10))+  &
              xmors2(mors1(il+1)-11)*yvect(mors2(mors1(il+1)-11))

   end if

end do 

end subroutine m2p1a


subroutine poliss(phi,ebj,mesh,vmsh)

!  Effectuer le calcul des composantes Ex et Ey du    
!  ---    champ electrique a partir du potentiel au noeud.
!         On appelle ca un lissage.                        

type(mesh_fields) :: ebj
type(sll_triangular_mesh_2d) :: mesh
type(voronoi)   :: vmsh
double precision, dimension(:) :: phi

integer :: is, ic, nbc, iac
LOGICAL :: lerr

lerr=.FALSE.

! --- 1.0 --- Calcul des termes individuels des seconds membres --------
  
do ic=1,mesh%nbtcot
   vtanty(ic)=(phi(vmsh%nuvac(1,ic))-phi(vmsh%nuvac(2,ic)))/vmsh%xlcod(ic)
end do

do ic=1,mesh%nbtcot
   vtantx(ic)=vtanty(ic)*vtaux(ic)
end do

do ic=1,mesh%nbtcot
   vtanty(ic)=vtanty(ic)*vtauy(ic)
end do

! --- 2.0 --- Calcul des seconds membres -------------------------------
  
do is=1,mesh%num_nodes

   iac=vmsh%nbcov(is)+1
   nbc=vmsh%nbcov(is+1)-vmsh%nbcov(is)

   if (nbc == 6) then
      sv1(is) = vtantx(vmsh%nugcv(iac  ))+vtantx(vmsh%nugcv(iac+1)) &
              + vtantx(vmsh%nugcv(iac+2))+vtantx(vmsh%nugcv(iac+3)) &
              + vtantx(vmsh%nugcv(iac+4))+vtantx(vmsh%nugcv(iac+5)) 
      sv2(is) = vtanty(vmsh%nugcv(iac  ))+vtanty(vmsh%nugcv(iac+1)) &
              + vtanty(vmsh%nugcv(iac+2))+vtanty(vmsh%nugcv(iac+3)) &
              + vtanty(vmsh%nugcv(iac+4))+vtanty(vmsh%nugcv(iac+5)) 
   else if (nbc == 5) then
      sv1(is) = vtantx(vmsh%nugcv(iac  ))+vtantx(vmsh%nugcv(iac+1)) &
              + vtantx(vmsh%nugcv(iac+2))+vtantx(vmsh%nugcv(iac+3)) &
              + vtantx(vmsh%nugcv(iac+4))
      sv2(is) = vtanty(vmsh%nugcv(iac  ))+vtanty(vmsh%nugcv(iac+1)) &
              + vtanty(vmsh%nugcv(iac+2))+vtanty(vmsh%nugcv(iac+3)) &
              + vtanty(vmsh%nugcv(iac+4))
   else if (nbc == 7) then
      sv1(is) = vtantx(vmsh%nugcv(iac  ))+vtantx(vmsh%nugcv(iac+1)) &
              + vtantx(vmsh%nugcv(iac+2))+vtantx(vmsh%nugcv(iac+3)) &
              + vtantx(vmsh%nugcv(iac+4))+vtantx(vmsh%nugcv(iac+5))     &
              + vtantx(vmsh%nugcv(iac+6))
      sv2(is) = vtanty(vmsh%nugcv(iac  ))+vtanty(vmsh%nugcv(iac+1)) &
              + vtanty(vmsh%nugcv(iac+2))+vtanty(vmsh%nugcv(iac+3)) &
              + vtanty(vmsh%nugcv(iac+4))+vtanty(vmsh%nugcv(iac+5))     &
              + vtanty(vmsh%nugcv(iac+6))
   else if (nbc == 4) then
      sv1(is) = vtantx(vmsh%nugcv(iac  ))+vtantx(vmsh%nugcv(iac+1)) &
              + vtantx(vmsh%nugcv(iac+2))+vtantx(vmsh%nugcv(iac+3))
      sv2(is) = vtanty(vmsh%nugcv(iac  ))+vtanty(vmsh%nugcv(iac+1)) &
              + vtanty(vmsh%nugcv(iac+2))+vtanty(vmsh%nugcv(iac+3))
   else if (nbc == 8) then
      sv1(is) = vtantx(vmsh%nugcv(iac  ))+vtantx(vmsh%nugcv(iac+1)) &
              + vtantx(vmsh%nugcv(iac+2))+vtantx(vmsh%nugcv(iac+3)) &
              + vtantx(vmsh%nugcv(iac+4))+vtantx(vmsh%nugcv(iac+5))     &
              + vtantx(vmsh%nugcv(iac+6))+vtantx(vmsh%nugcv(iac+7)) 
      sv2(is) = vtanty(vmsh%nugcv(iac  ))+vtanty(vmsh%nugcv(iac+1)) &
              + vtanty(vmsh%nugcv(iac+2))+vtanty(vmsh%nugcv(iac+3)) &
              + vtanty(vmsh%nugcv(iac+4))+vtanty(vmsh%nugcv(iac+5))     &
              + vtanty(vmsh%nugcv(iac+6))+vtanty(vmsh%nugcv(iac+7)) 
   else if (nbc == 3) then
      sv1(is) = vtantx(vmsh%nugcv(iac  ))+vtantx(vmsh%nugcv(iac+1)) &
              + vtantx(vmsh%nugcv(iac+2))
      sv2(is) = vtanty(vmsh%nugcv(iac  ))+vtanty(vmsh%nugcv(iac+1)) &
              + vtanty(vmsh%nugcv(iac+2))
   else if (nbc == 9) then
      sv1(is) = vtantx(vmsh%nugcv(iac  ))+vtantx(vmsh%nugcv(iac+1)) &
              + vtantx(vmsh%nugcv(iac+2))+vtantx(vmsh%nugcv(iac+3)) &
              + vtantx(vmsh%nugcv(iac+4))+vtantx(vmsh%nugcv(iac+5))     &
              + vtantx(vmsh%nugcv(iac+6))+vtantx(vmsh%nugcv(iac+7))     &
              + vtantx(vmsh%nugcv(iac+8))
      sv2(is) = vtanty(vmsh%nugcv(iac  ))+vtanty(vmsh%nugcv(iac+1)) &
              + vtanty(vmsh%nugcv(iac+2))+vtanty(vmsh%nugcv(iac+3)) &
              + vtanty(vmsh%nugcv(iac+4))+vtanty(vmsh%nugcv(iac+5))     &
              + vtanty(vmsh%nugcv(iac+6))+vtanty(vmsh%nugcv(iac+7))     &
              + vtanty(vmsh%nugcv(iac+8))
   else if (nbc == 2) then
      sv1(is) = vtantx(vmsh%nugcv(iac  ))+vtantx(vmsh%nugcv(iac+1))
      sv2(is) = vtanty(vmsh%nugcv(iac  ))+vtanty(vmsh%nugcv(iac+1))
   else if (nbc == 10) then
      sv1(is) = vtantx(vmsh%nugcv(iac  ))+vtantx(vmsh%nugcv(iac+1)) &
              + vtantx(vmsh%nugcv(iac+2))+vtantx(vmsh%nugcv(iac+3)) &
              + vtantx(vmsh%nugcv(iac+4))+vtantx(vmsh%nugcv(iac+5))     &
              + vtantx(vmsh%nugcv(iac+6))+vtantx(vmsh%nugcv(iac+7))     &
              + vtantx(vmsh%nugcv(iac+8))+vtantx(vmsh%nugcv(iac+9)) 
      sv2(is) = vtanty(vmsh%nugcv(iac  ))+vtanty(vmsh%nugcv(iac+1)) &
              + vtanty(vmsh%nugcv(iac+2))+vtanty(vmsh%nugcv(iac+3)) &
              + vtanty(vmsh%nugcv(iac+4))+vtanty(vmsh%nugcv(iac+5))     &
              + vtanty(vmsh%nugcv(iac+6))+vtanty(vmsh%nugcv(iac+7))     &
              + vtanty(vmsh%nugcv(iac+8))+vtanty(vmsh%nugcv(iac+9)) 
   else if (nbc == 11) then
      sv1(is) = vtantx(vmsh%nugcv(iac  ))+vtantx(vmsh%nugcv(iac+1)) &
              + vtantx(vmsh%nugcv(iac+2))+vtantx(vmsh%nugcv(iac+3)) &
              + vtantx(vmsh%nugcv(iac+4))+vtantx(vmsh%nugcv(iac+5))     &
              + vtantx(vmsh%nugcv(iac+6))+vtantx(vmsh%nugcv(iac+7))     &
              + vtantx(vmsh%nugcv(iac+8))+vtantx(vmsh%nugcv(iac+9))     &
              + vtantx(vmsh%nugcv(iac+10))
      sv2(is) = vtanty(vmsh%nugcv(iac  ))+vtanty(vmsh%nugcv(iac+1)) &
              + vtanty(vmsh%nugcv(iac+2))+vtanty(vmsh%nugcv(iac+3)) &
              + vtanty(vmsh%nugcv(iac+4))+vtanty(vmsh%nugcv(iac+5))     &
              + vtanty(vmsh%nugcv(iac+6))+vtanty(vmsh%nugcv(iac+7))     &
              + vtanty(vmsh%nugcv(iac+8))+vtanty(vmsh%nugcv(iac+9))     &
              + vtanty(vmsh%nugcv(iac+10))
   else if (nbc == 12) then
      sv1(is) = vtantx(vmsh%nugcv(iac  ))+vtantx(vmsh%nugcv(iac+1)) &
              + vtantx(vmsh%nugcv(iac+2))+vtantx(vmsh%nugcv(iac+3)) &
              + vtantx(vmsh%nugcv(iac+4))+vtantx(vmsh%nugcv(iac+5))     &
              + vtantx(vmsh%nugcv(iac+6))+vtantx(vmsh%nugcv(iac+7))     &
              + vtantx(vmsh%nugcv(iac+8))+vtantx(vmsh%nugcv(iac+9))     &
              + vtantx(vmsh%nugcv(iac+10))+vtantx(vmsh%nugcv(iac+11))
      sv2(is) = vtanty(vmsh%nugcv(iac  ))+vtanty(vmsh%nugcv(iac+1)) &
              + vtanty(vmsh%nugcv(iac+2))+vtanty(vmsh%nugcv(iac+3)) &
              + vtanty(vmsh%nugcv(iac+4))+vtanty(vmsh%nugcv(iac+5))     &
              + vtanty(vmsh%nugcv(iac+6))+vtanty(vmsh%nugcv(iac+7))     &
              + vtanty(vmsh%nugcv(iac+8))+vtanty(vmsh%nugcv(iac+9))     &
              + vtanty(vmsh%nugcv(iac+10))+vtanty(vmsh%nugcv(iac+11))

   else
      lerr=.TRUE.
   end if

end do

if (lerr) then
   write(iout,900)
   stop
end if

! --- 3.0 --- Resolution des systemes lineaires 2*2 --------------------

do is=1,mesh%num_nodes
   ebj%e(1,is)=mesh%xmal2(is)*sv1(is)-mesh%xmal3(is)*sv2(is)
   ebj%e(2,is)=mesh%xmal1(is)*sv2(is)-mesh%xmal3(is)*sv1(is)
end do

! --- 9.0 --- Formats --------------------------------------------------
    
900 format(//10x,'Erreur dans POLISS'   &
            /10x,'On a trouve plus de 12 cotes')

end subroutine poliss


 
! ==================================================================== *
!  616-AMPLIM        Version 1.0     Juillet  1992   E.Puertolas       *
! ==================================================================== *
!                                                                      *
!         But :     Ordonner dans le sens trigonometrique les noeuds   *
!                   Dirichlet non homogenes.                           *
!                                                                      *
!      Variables en arguments :                                        *
!      ------------------------                                        *
!                                                                      *
!             ntri   : numero des sommets des triangles                *
!             nvois  : numeros des voisins des elements                *
!             krefro : numeros de reference des cotes frontieres       *
!             inoeuf : tableau noeuds frontieres                       *
!                                                                      *
!             ieldir : tableau temporaire (elements frontieres)        *
!             ictdir : tableau temporaire (cotes correspondants)       *

!             iout   : unite logique d'impression                      *
!             nnref  : nombre de references de frontieres Dirichlet    *
!                      non homogenes                                   *
!             irefdir: references de frontieres Dirichlet non          *
!                      homogenes                                       *
!                                                                      *
!             nnoeuf : nombre de noeuds Dirichlet non homogenes        *
!                                                                      *
! ==================================================================== *
!      Appelant : 61-INIC00                                          *
! ==================================================================== *
subroutine amplim(mesh,inoeuf,ieldir,ictdir)
type(sll_triangular_mesh_2d) :: mesh

integer :: iel, iref, neldir
integer :: is1, in1, jn1, ie, je, jel, numel1
integer :: is2, js2, in2, jn2, itmp1, itmp2
integer :: ieldir(*),ictdir(*), inoeuf(*)
logical :: lcomm

! ... Recherche des elements concernes
 
neldir=0
DO iref=1,nnref
   DO iel=1,mesh%num_cells
      IF(mesh%nvois(1,iel).LT.0) THEN
         IF(mesh%krefro(-mesh%nvois(1,iel)).EQ.irefdir(iref)) THEN
            neldir=neldir+1
            ieldir(neldir)=iel
            ictdir(neldir)=1
         ENDIF
      ENDIF
      IF(mesh%nvois(2,iel).LT.0) THEN
         IF(mesh%krefro(-mesh%nvois(2,iel)).EQ.irefdir(iref)) THEN
            neldir=neldir+1
            ieldir(neldir)=iel
            ictdir(neldir)=2
         ENDIF
      ENDIF
      IF(mesh%nvois(3,iel).LT.0) THEN
         IF(mesh%krefro(-mesh%nvois(3,iel)).EQ.irefdir(iref)) THEN
            neldir=neldir+1
            ieldir(neldir)=iel
            ictdir(neldir)=3
         ENDIF
      ENDIF
   end do
end do

! ... Recherche du 1er element dans le sens trigo

     numel1=0
         DO 20 ie=1,neldir
        iel=ieldir(ie)
        in1=ictdir(ie)
        is1=mesh%nodes(in1,iel)
        lcomm=.FALSE.
        DO 30 je=1,neldir
           IF(je.NE.ie) THEN
              jel=ieldir(je)
              jn1=ictdir(je)
              jn2=mod(jn1+3,3)+1
              js2=mesh%nodes(jn2,jel)
                  IF(is1.EQ.js2) THEN
             lcomm=.TRUE.
              ENDIF
           ENDIF
 30         CONTINUE
        
        IF((.NOT.lcomm).AND.(numel1.NE.0)) THEN
           WRITE(iout,901) 
               GOTO 800
        ENDIF
        IF((.NOT.lcomm).AND.(numel1.EQ.0)) THEN
           numel1=ie
        ENDIF
 20      CONTINUE


! ... Classement des elements et des noeuds dans l'ordre trigo

         itmp1=ieldir(1)
         itmp2=ictdir(1)
         ieldir(1)=ieldir(numel1)
         ictdir(1)=ictdir(numel1)
         ieldir(numel1)=itmp1
         ictdir(numel1)=itmp2
     in1=ictdir(1)
     in2=mod(in1+3,3)+1
     inoeuf(1)=mesh%nodes(in1,ieldir(1))
     inoeuf(2)=mesh%nodes(in2,ieldir(1))
     nnoeuf=2
     DO 40 ie=1,neldir-1
        iel=ieldir(ie)
        is2=inoeuf(nnoeuf)

            DO 50 je=ie+1,neldir
           jel=ieldir(je)
           jn1=ictdir(je)
           is1=mesh%nodes(jn1,jel)
           IF(is1.EQ.is2) THEN
                  itmp1=ieldir(ie+1)
                  itmp2=ictdir(ie+1)
                  ieldir(ie+1)=ieldir(je)
                  ictdir(ie+1)=ictdir(je)
                  ieldir(je)=itmp1
                  ictdir(je)=itmp2
              nnoeuf=nnoeuf+1
          jn2=mod(jn1+3,3)+1
              inoeuf(nnoeuf)=mesh%nodes(jn2,ieldir(ie+1))
           ENDIF
 50         CONTINUE
 40       CONTINUE


! --- 8.0 --- Point de sortie du sous-programme ------------------------

      GOTO 850
 800  CONTINUE
      WRITE(iout,900) 
!     CALL utabrt
 850  CONTINUE

! --- 9.0 --- Format ---------------------------------------------------
 
 900  FORMAT(//10x,'Erreur dans AMPLIM  '///)
 901  FORMAT(//10x,'La frontiere est constituee de plusieurs'   &
     &            ,' segments disjoints'/)

!     ___
end subroutine amplim

!                                                                      i
!     asblm3 : Assembler 3 matrices elementaires dans 3 matrices       i
!              globales stockees sous forme "morse" non symetriques.   i
!                                                                      i
subroutine asblm3(aele1,aele2,aele3,mors1,mors2,i1,i2,i3,a1,a2,a3)

! ==================================================================== i
!                                                                      i
!          Version 1.0  Aout 1989  J. Segre                            i
!                                                                      i
!                                                                      i
!          BUT :         assembler 3 matrices elementaires             i
!          -----         dans 3 matrices globales stockees             i
!                        sous forme "morse" non symetriques            i
!                                                                      i
!          PARAMETRES D'ENTREE :                                       i
!          ---------------------                                       i
!                                                                      i
!          aele1,..,aele3 :  Matrices elementaires                     i
!          mors1          :  Tableau du nombre de termes par           i
!                            ligne de la matrice morse                 i
!          mors2          :  Tableau des numeros des termes            i
!                            de chaque ligne de la matrice morse       i
!          i1,i2,i3       :  Numeros des sommets de l'element          i
!                                                                      i
!          PARAMETRES RESULTATS :                                      i
!          ---------------------                                       i
!                                                                      i
!          a1,..a3        :  Matrices globales stockees                i
!                            sous forme "morse"                        i
!                                                                      i
! ==================================================================== i
 
integer :: i1, i2, i3, j1, j2, ind1, ind2, ind3
double precision, dimension(:) :: aele1,aele2,aele3,a1,a2,a3
integer, dimension(:) ::  mors1, mors2

! ----------------------------------------------------------------------
! --- 1.1 --- Rangement des termes diagonaux ---------------------------

ind1=mors1(i1+1)
a1(ind1)=a1(ind1)+aele1(1)
a2(ind1)=a2(ind1)+aele2(1)
a3(ind1)=a3(ind1)+aele3(1)
 
ind2=mors1(i2+1)
a1(ind2)=a1(ind2)+aele1(5)
a2(ind2)=a2(ind2)+aele2(5)
a3(ind2)=a3(ind2)+aele3(5)
 
ind3=mors1(i3+1)
a1(ind3)=a1(ind3)+aele1(9)
a2(ind3)=a2(ind3)+aele2(9)
a3(ind3)=a3(ind3)+aele3(9)

! --- 1.2 --- Rangement des autres termes ------------------------------

j2=ind1-1
j1=mors1(i1)+1
 
IF(j2.GE.j1) THEN
   DO j=j1,j2
      IF(i2.EQ.mors2(j)) THEN
         a1(j)=a1(j)+aele1(2)
         a2(j)=a2(j)+aele2(2)
         a3(j)=a3(j)+aele3(2)
      ENDIF
      IF(i3.EQ.mors2(j)) THEN
         a1(j)=a1(j)+aele1(3)
         a2(j)=a2(j)+aele2(3)
         a3(j)=a3(j)+aele3(3)
      ENDIF
   end do
ENDIF
 
j2=ind2-1
j1=mors1(i2)+1
 
IF(j2.GE.j1) THEN
   DO j=j1,j2
      IF(i1.EQ.mors2(j)) THEN
         a1(j)=a1(j)+aele1(4)
         a2(j)=a2(j)+aele2(4)
         a3(j)=a3(j)+aele3(4)
      ENDIF
      IF(i3.EQ.mors2(j)) THEN
         a1(j)=a1(j)+aele1(6)
         a2(j)=a2(j)+aele2(6)
         a3(j)=a3(j)+aele3(6)
      ENDIF
   end do
ENDIF

j2=ind3-1
j1=mors1(i3)+1
 
IF(j2.GE.j1) THEN
   DO j=j1,j2
      IF(i1.EQ.mors2(j)) THEN
         a1(j)=a1(j)+aele1(7)
         a2(j)=a2(j)+aele2(7)
         a3(j)=a3(j)+aele3(7)
      ENDIF
      IF(i2.EQ.mors2(j)) THEN
         a1(j)=a1(j)+aele1(8)
         a2(j)=a2(j)+aele2(8)
         a3(j)=a3(j)+aele3(8)
      ENDIF
   end do
ENDIF

end subroutine asblm3   

end module poisson
