!File: Module Solveurs_Module
!Definition des solveurs de Maxwell et de Poisson
module solveurs_module

use zone,          only: mesh_fields, lmodte, lmodtm,   &
                 lcorrp, nesmx, dt, iout, nesp,     &
                 objet_fonction

!----------------------------------------------------------------------
implicit none

integer, private :: i, j, k

!Variable: ncxmx 
!nombre maximum de champs externes par mode TE ou TM
integer, parameter :: ncxmx = 10
!Variable: nanmx 
!nombre maximum d'antennes 
integer, parameter :: nanmx = 5
!Variable: nspmx 
!nombre maximum de points d'echantillonnage des champs
integer, parameter :: nspmx = 20
!Variable: nsdmx 
!nombre maximum de sous-domaines              
integer, parameter :: nsdmx = 20

!Variable: nfrmx 
!nombre maximum de frontieres referencees
integer, parameter :: nfrmx =  20
!Variable: ntypfr 
!Type de condition aux limites par frontiere      
!                     referencee sur les champs et les particules :    
!               0 - axe            (reflechissant pour les particules) 
!               1 - conducteur parfait (absorbant pour les particules) 
!               2 - mur magnetique (reflechissant pour les particules) 
!               3 - front absorbante   (absorbant pour les particules) 
!               4 - onde entrante      (absorbant pour les particules) 
!               5 - conducteur parfait (reflexion pour les particules)
!               6 - front absorbante   (reflexion pour les particules) 
!               7 - front ddp imposee  (absorbant pour les particules) 

                                                                      
!Variables: Conditions aux limites
!            potfr  - potentiel applique sur les frontieres Dirichlet  
!            ddpfr  - valeur de la difference de potentiel le long des 
!                     frontieres ou l'on impose une condition de type  
!                     ddp imposee (ntypfr = 7)                         
!            ddpft  - valeur de la ddp a chaque instant:               
!                     ddpfr*enveloppe temporelle definie par numddp    
!            numddp - numero de la dependance temporelle (idem que     
!                     numond) a appliquer a la ddp                     
!                                                                      
!          Dans le cas de POISSON ces frontieres sont                  
!
!               0 - axe            (reflechissant pour les particules) 
!               1 - DIRICHLET      (absorbant pour les particules)     
!               3 - NEUMAN         (absorbant pour les particules)     
!               5 - DIRICHLET      (reflexion pour les particules)     
!               6 - NEUMAN         (reflexion pour les particules)     
!
!                                                                      
!            ntypfb - type de frontiere pour Ampere (uniquement)  
!            b3fr   - champ magnetique sur une frontiere Dirichlet
!                      pour Ampere


!Variables: Caracteristiques electromagnetiques des sous-domaines       
!                                                                      
!             epssd - permittivite relative ou constante dielectrique 
!                     par sous-domaine                               
!                                                                   
!             xmusd - permeabilite relative                        
!                     par sous-domaine                            
!                                                                
!             sigsd - conductivite electrique (0=infini)        
!                     par sous-domaine                         
!                                                             

! -------------------------------------------------------------------- 
!                                                                      
!Variables: Definition de l'onde entrante (ntypfr=4 pour Maxwell)        
!         de la forme creneau*E0*cos(2*pi*freq*temps+phase)            
!                                                                      
!                     Mode 1 =  E1,E2,B3                               
!                                                                      
!            numond - numero d'onde par reference de frontiere         
!            nu2ond - 2eme onde par reference de frontiere             
!            nu3ond - 3eme onde par reference de frontiere             
!            nu4ond - 4eme onde par reference de frontiere             
!            nu5ond - 5eme onde par reference de frontiere             
!            nu6ond - 6eme onde par reference de frontiere             
!            ifoond - code forme de l'enveloppe                        
!            tdbond - temps de debut d'injection                       
!            tmtond - temps de montee                                  
!            tplond - duree du plateau                                 
!            tdsond - temps de descente                                
!            freond - frequence de l'onde                              
!            phaond - phase                                            
!            eleond - amplitude du champ electrique par arete          
!            tauond - grandeur "tau" intervenant quand ifoond=6        
!            sigond - grandeur "sigma" intervenant quand ifoond=6      
!                                                                      
!                     Mode 2 =  B1,B2,E3                               
!                                                                      
!            numon2 - numero d'onde par reference de frontiere         
!            ifoon2 - code forme de l'enveloppe                        
!            tdbon2 - temps de debut d'injection                       
!            tmton2 - temps de montee                                  
!            tplon2 - duree du plateau                                 
!            tdson2 - temps de descente                                
!            freon2 - frequence de l'onde                              
!            phaon2 - phase                                            
!            eleon2 - amplitude du champ electrique par arete          
!            tauon2 - grandeur "tau" intervenant quand ifoon2=6        
!            sigon2 - grandeur "sigma" intervenant quand ifoon2=6      
!                                                                      
! --------------------------------------------------------------------


!
!
!Variables: Frontieres superposees transparentes 
!pour les particules 
!         correspondant a des frontieres de 2 domaines colles        
!                                                                     
!            nbfrtp - nombre de frontieres transparentes pour les part 
!                     (donnees par le mailleur)                    
!            ireftp - numero de reference de chaque frontiere     
!                     (donnees par le mailleur)                  
!            ltrftp - activateur de "transparence" par particule et par
!                     frontiere double (domaines colles)              
!            nnftp  - nombre de noeuds sur ces frontieres            
!                                                                   
!

integer :: nbfrtp, ireftp(nfrmx), nnftp
logical :: ltrftp(nfrmx,nesmx) 

!
!Variables: Variation temporelle du potentiel 
!applique sur un Dirichlet  
!
!            tdbpot - debut de la rampe (s)                            
!            tmtpot - temps de montee   (s)                           
!            tplpot - duree du plateau  (s)                          
!            tdspot - temps de descente (s)                         
!            ifopot - forme de l'impulsion : 0=pas de dependance en t 
!                                            1=trapeze               
!                                            2=sinusoides aux 2 bouts 
!                                            3=  ...                 
!

!Variables: Definition de courants externes 
!(antennes)  
!
!
!            lanten - Activateur de cette option                      
!            nanten - nombre d'antennes                              
!            nopant - option                                        
!            curant - densite de courant                           
!            rayant - rayon ou largeur de l'antenne (option=1)    
!            xyantn - definition du volume                       
!            freant - frequence                                 
!            phaant - phase                                    
!            ifoant - code forme de l'impulsion               
!            tmtant - temps de montee                        
!            tdbant - debut de la rampe de montee           
!            tplant - duree du plateau                     
!            tdsant - temps de descente                   
!            tauant - grandeur "tau" intervenant quand ifoant=6     
!            sigant - grandeur "sigma" intervenant quand ifoant=6    
!            noeant - nombre de noeuds impliques                   
!                                                                 
!

logical :: lanten 
integer :: nanten
integer :: nopant(nanmx)
double precision    :: freant(nanmx), phaant(nanmx)
integer :: ifoant(nanmx) 
double precision    :: tdbant(nanmx), tmtant(nanmx), tplant(nanmx), tdsant(nanmx)
double precision    :: curant(nanmx), rayant(nanmx), xyantn(7,nanmx)
double precision    :: tauant(nanmx), sigant(nanmx)
integer :: noeant(nanmx) 

!
!
!Variables: Formes temporelles definies point par point  
!
!        Dedendances temporelles point par point (de 1 a ncxmx)     
!        (Nombre maximum de points = 100)                            
!
!        Appel avec iforme=11,12, .... ==>  forme=iforme-10           
!                                                                  
!             nbpdon - nombre de points                           
!             ttpdon - temps definis point par point               
!             vtpdon - valeurs de la fonction aux points precedents 
!             ippdon - pointeur (utilisation interne)                
!                                                                      


!      Definition des frontieres internes    
!
!         correspondant a des lignes de maillage                
!         (plan x=x0 en cartesien ; plan z=z0 en axisymetrique)  
!                                                        
!            nbfrnt - nombre de frontieres internes 
!            z1frnt - position z minimum                  
!            z2frnt - position z maximum                   
!            r1frnt - rayon ou ordonnee minimum             
!            r2frnt - rayon ou ordonnee maximum              
!            itfrnt - type de frontiere pour les champs       
!                       1 - conducteur parfait pour les champs 
!                       2 - inexistant  pour les champs         
!            ipfrnt - condition pour les particules              
!                       1 - transparent pour les particules       
!                       2 - absorbant pour les particules          
!            irefnt - numero de reference de la frontiere           
!                     (ce numero doit etre superieur au plus grand   
!                      numero de reference fourni par le mailleur)    
!                                                                      
!         Frontieres superposees transparentes pour les particules     
!         correspondant a des frontieres de 2 domaines colles         
!                                                                    
!            nbfrtp - nombre de frontieres transparentes pour les part
!                     (donnees par le mailleur)                      
!            ireftp - numero de reference de chaque frontiere       
!                     (donnees par le mailleur)                    
!            ltrftp - activateur de "transparence" par particule et par
!                     frontiere double (domaines colles)              
!            nnftp  - nombre de noeuds sur ces frontieres            
!                                                                   


!Type: mesh_bound
type mesh_bound
   integer :: nbfrnt
   integer, dimension(:),   pointer :: ntypfr, ntypfb, ntypfp
   integer, dimension(:,:), pointer :: ipfrnt
   integer, dimension(:),   pointer :: itfrnt, irefnt
   double precision,    dimension(:),   pointer :: r1frnt, r2frnt, z1frnt,z2frnt
   logical :: laonde, laond2
   integer, dimension(:),   pointer :: numond, nu2ond, nu3ond
   integer, dimension(:),   pointer :: nu4ond, nu5ond, nu6ond
   integer, dimension(:),   pointer :: ifoond
   double precision   , dimension(:),   pointer :: tdbond, tmtond, tplond, tdsond
   double precision   , dimension(:),   pointer :: freond, phaond, eleond
   double precision   , dimension(:),   pointer :: tauond, sigond
   integer, dimension(:),   pointer :: numon2, ifoon2
   double precision,    dimension(:),   pointer :: tdbon2, tmton2, tplon2, tdson2
   double precision,    dimension(:),   pointer :: freon2, phaon2, eleon2
   double precision,    dimension(:),   pointer :: tauon2, sigon2
   double precision,    dimension(:),   pointer :: tdbpot, tmtpot, tplpot, tdspot
   integer, dimension(:),   pointer :: ifopot
   logical :: laddp
   integer, dimension(:),   pointer :: numddp
   double precision,    dimension(:),   pointer :: potfr, ddpfr, ddpft, b3fr
   logical :: lepssd, lxmusd, lsigsd
   double precision,    dimension(:),   pointer :: epssd, xmusd, sigsd
end type mesh_bound

integer :: degrek=1 , Nk

contains

!Subroutine: lecture_donnees_solveur
!  Lecture des donnees relatives aux champs                 
!   Namelist "nlcham"                                        
!
!Traduction en Fortran 90 de leccha
!
!  33-LECCHA  -   Version 1.0  Novembre 1991   A. Adolf  
subroutine lecture_donnees_solveur(nomfich, nsolve, ntypfr, potfr, nmxfr, nmxsd)


type (mesh_bound) :: bcnd
integer, intent(in) :: nmxfr, nmxsd
integer, intent(in) :: nsolve
character(len=*), intent(in) :: nomfich
integer :: ip, jp, ityp, ifr, isd, iesp
integer :: nbfrnt

integer, dimension(nfrmx) :: ntypfr, ntypfb, ntypfp
integer, dimension(nfrmx,nesmx) :: ipfrnt
integer, dimension(nfrmx) :: itfrnt, irefnt
double precision,    dimension(nfrmx) :: r1frnt, r2frnt, z1frnt, z2frnt

integer :: numond(nfrmx), nu2ond(nfrmx), nu3ond(nfrmx)
integer :: nu4ond(nfrmx), nu5ond(nfrmx), nu6ond(nfrmx)
integer :: ifoond(nfrmx)

double precision    :: tdbond(nfrmx), tmtond(nfrmx), tplond(nfrmx), tdsond(nfrmx)
double precision    :: freond(nfrmx), phaond(nfrmx), eleond(nfrmx)
double precision    :: tauond(nfrmx), sigond(nfrmx)

integer :: numon2(nfrmx), ifoon2(nfrmx)
double precision    :: tdbon2(nfrmx), tmton2(nfrmx), tplon2(nfrmx), tdson2(nfrmx)
double precision    :: freon2(nfrmx), phaon2(nfrmx), eleon2(nfrmx)
double precision    :: tauon2(nfrmx), sigon2(nfrmx)

double precision    :: tdbpot(nfrmx), tmtpot(nfrmx), tplpot(nfrmx), tdspot(nfrmx)
integer :: ifopot(nfrmx)
integer :: numddp(nfrmx)
double precision    :: potfr(nfrmx), ddpfr(nfrmx), ddpft(nfrmx), b3fr(nfrmx)

logical :: lepssd=.false., lxmusd=.false., lsigsd=.false.
double precision    :: epssd(0:nsdmx), xmusd(0:nsdmx), sigsd(nsdmx)

integer :: nbpdon(ncxmx), ippdon(ncxmx)
double precision    :: ttpdon(100,ncxmx), vtpdon(100,ncxmx) 


NAMELIST/nlcham/ntypfr,potfr,ddpfr,numddp,ntypfb,b3fr       &
               ,numond,nu2ond,nu3ond,nu4ond,nu5ond,nu6ond   &
               ,ifoond,tdbond,tmtond,tplond,tdsond          &
               ,freond,phaond,eleond,tauond,sigond          &
               ,numon2,ifoon2,tdbon2,tmton2,tplon2,tdson2   &
               ,freon2,phaon2,eleon2,tauon2,sigon2          &
               ,nbfrnt,z1frnt,z2frnt,itfrnt,r1frnt,r2frnt   &
               ,irefnt,ipfrnt,nbfrtp,ireftp,ltrftp          &
               ,epssd ,xmusd , sigsd                        &
               ,tdbpot,tmtpot,tplpot,tdspot,ifopot          &
               ,lanten,nanten,nopant,ifoant,freant          &
               ,phaant,rayant,tmtant,tdbant,tplant          &
               ,tdsant,curant,xyantn,tauant,sigant          &
               ,nbpdon,ttpdon,vtpdon                        

! ----------- Valeurs par defaut et initialisations -------------------t
write(iout,900)

ntypfr=1;  potfr =0.; ddpfr =0.; numddp= 1; ntypfb=1; ntypfp=1;
tdbpot=0.; tmtpot=0.; tplpot=0.; tdspot=0.; ifopot=0
b3fr = 0.0

numond=1
nu2ond=0 ; nu3ond=0 ; nu4ond=0;  nu5ond=0 ; nu6ond=0
ifoond=0 ; tdbond=0.; tmtond=0.; tplond=1.; tdsond=0.
freond=0.; phaond=0.; eleond=0.; tauond=0.; sigond=0.

numon2=1
ifoon2=0 ; tdbon2=0.; tmton2=0.; tplon2=1.
tdson2=0.; freon2=0.; phaon2=0.; eleon2=1.
tauon2=0.; sigon2=0.

z1frnt=0.; z2frnt=0.; r1frnt=0.; r2frnt=0.; itfrnt=1
irefnt=0 ; ipfrnt=1 ; ltrftp=.false.
ireftp=0 ; nbfrnt=0 ; nbfrtp=0

epssd=1.; xmusd=1.; sigsd=0.

!... Sources de courants externes

lanten = .FALSE.; nanten = 1

nopant=0 ; ifoant=0 ; freant=0.  ; phaant=0.
tdbant=0.; tmtant=0.; tplant=1.  ; tdsant=0.
curant=0.; rayant=0.; xyantn=0.; tauant=0.
sigant=0.

!... Definition des formes temporelles point par point

nbpdon=2; ippdon=1; ttpdon=0.; vtpdon=0.; ttpdon=1.

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

!open(10, file = nomfich)
!!* Proprietes des materiaux (LIGNE 5) *
!read(10,*)s1, n1
!if ( s1 == '*MATERIAUX' ) then
!   do i = 1, n1
!      read(10,*) n2, s2, s3, s4, x1, s5, x2
!      select case(s2)
!      case('VIDE') 
!      case('AIR') 
!      case('DIELECTRIQUE') 
!         if (s3 == 'ISOTROPE') then
!            if ( s4 == 'EPSILON' ) eps0 = x1           
!       if ( s5 == 'MU'      ) xmu0 = x2
!         else
!            write(*,*) " Probleme de lecture LIGNE 5 "
!         end if
!      case('PLASMA') 
!      end select
!   end do
!else
!   write(*,*) " Probleme de lecture LIGNE 5 "
!end if
!close(10)

!--- 2.5 --- Ecriture des donnees -------------------------------------

write(iout,921)

do i=1,nfrmx

   write(iout,922) i,ntypfr(i),potfr(i),ddpfr(i),numddp(i),     &
             ntypfb(i), b3fr(i)

   write(iout,9223) i,ntypfp(i)

   if(nsolve /= 10 .and. nsolve /= 4 .and. ntypfr(i) == 7) then
      write(iout,9221) i,ntypfr(i),nsolve
      call errout(iout,"F","init_solveur_maxwell"," ")
   end if

   if(.NOT.lmodte  .and. ntypfr(i) == 7) then
      write(iout,9222) i,ntypfr(i),lmodte
      call errout(iout,"F","init_solveur_maxwell"," ")
   end if

end do

if(lmodte) then

   write(iout,945)

   do i=1,nfrmx
      write(iout,946) i,numond(i),nu2ond(i),nu3ond(i)   &
                       ,nu4ond(i),nu5ond(i),nu6ond(i)
   end do

   write(iout,938)

   do i=1,nfrmx

      write(iout,939) i,freond(i),phaond(i),eleond(i)   &
                       ,tdbond(i),tmtond(i),tplond(i)   &
               ,tdsond(i),ifoond(i)

      if(ifoond(i) == 6) then
         write(iout,9391) tauond(i),sigond(i)
         if(ABS(sigond(i)) <  1.e-20*dt) then
            write(iout,914)
            call errout(iout,"F","initialisation_solveurs"," ")
         end if
      end if

   end do
end if

if(lmodtm) then
   write(iout,947)
   do i=1,nfrmx
      write(iout,948) i,numond(i)
   end do
   write(iout,949)
   do i=1,nfrmx
      write(iout,939) i,freon2(i),phaon2(i),eleon2(i)   &
                       ,tdbon2(i),tmton2(i),tplon2(i)   &
               ,tdson2(i),ifoon2(i)
      if(ifoon2(i) == 6) then
         write(iout,9392) tauon2(i),sigon2(i)
         if(ABS(sigon2(i)) <  1.e-20*dt) then
            write(iout,915)
            call errout(iout,"F","initialisation_solveurs"," ")
         end if
      end if
   end do
end if

write(iout,934) nbfrnt
if(nbfrnt >  0) then
   do i=1,nbfrnt
      write(iout,935) i,z1frnt(i),z2frnt(i),r1frnt(i)   &
                   ,r2frnt(i),itfrnt(i),irefnt(i)
      do iesp=1,nesmx
         write(iout,9351) iesp,ipfrnt(i,iesp)
      enddo
   end do
end if

write(iout,936) nbfrtp
if(nbfrtp >  0) then
   do i=1,nbfrtp
      write(iout,937) i,ireftp(i)
      do iesp=1,nesmx
         write(iout,9371) iesp,ltrftp(i,iesp)
      end do
   end do
end if

write(iout,923)
do i=1,nfrmx
   write(iout,924) i,tdbpot(i),tmtpot(i),tplpot(i)  &
                    ,tdspot(i),ifopot(i)
end do

write(iout,906)
do isd=1,nsdmx
   write(iout,907) isd,epssd(isd),xmusd(isd),sigsd(isd)
end do

write(iout,906)
do ip=1,ncxmx
   if(nbpdon(ip) >  100) then
      write(iout,942) ip,nbpdon(ip)
      call errout(iout,"F","initialisation_solveurs"," ")
   end if
   if(nbpdon(ip) >  2) then
      write(iout,943) ip,nbpdon(ip)
      write(iout,944) (i,ttpdon(i,ip),vtpdon(i,ip),i=1,nbpdon(ip))
   end if
end do

!--- 3.0 --- Controle des donnees -------------------------------------

write(iout,932)
do ifr=1,nfrmx
   ityp=ntypfr(ifr)
   if(ityp == 0) then
      write(iout,901) ifr
   else if(ityp == 1) then
      write(iout,902) ifr
   else if(ityp == 2) then
      write(iout,903) ifr
      call errout(iout,"F","initialisation_solveurs"," ")
   else if(ityp == 3) then
      write(iout,904) ifr
   else if(ityp == 4) then
      write(iout,905) ifr
   else if(ityp == 5) then
      write(iout,952) ifr
   else if(ityp == 6) then
      write(iout,953) ifr
   else if(ityp == 7) then
      write(iout,954) ifr
   else
      write(iout,910) ifr,ityp
      call errout(iout,"F","initialisation_solveurs"," ")
   end if
   if (ntypfp(ifr) == 0) then
      write(iout,"(/10x,'Reference :',i3,5x,'P nul'/)") ntypfp(ifr)
   else
      write(iout,"(/10x,'Reference :',i3,5x,'P = E.n'/)") ntypfp(ifr)
   end if
end do

do ifr=1,nfrmx

   if(      (ntypfr(ifr) == 1.OR.ntypfr(ifr) == 5)  &
      .and. (ifopot(ifr) >  0)  &
      .and. (ifopot(ifr) <= 10) ) then

      if(tmtpot(ifr)+tplpot(ifr)+tdspot(ifr) == 0) then
         write(iout,913) ifr
         call errout(iout,"F","initialisation_solveurs"," ")
      end if

      if((tmtpot(ifr) <  0).OR. &
         (tplpot(ifr) <  0).OR. &
         (tdspot(ifr) <  0))     then
         write(iout,913) ifr
         call errout(iout,"F","initialisation_solveurs"," ")
      end if
   end if

   if((ntypfr(ifr) == 1.OR.ntypfr(ifr) == 5)    &
                     .and.(ifopot(ifr) == 6) ) then
      write(iout,917) ifr,ifopot(ifr)
      call errout(iout,"F","initialisation_solveurs"," ")
   end if

end do

do i=1,nsdmx
   if((epssd(i) <= 0).OR.(xmusd(i) <= 0).OR.(sigsd(i) <  0)) then
      write(iout,933)
      call errout(iout,"F","initialisation_solveurs"," ")
   end if
end do

!... Test de presence de materiaux eps,xmu,sig 

lepssd=.FALSE.
lxmusd=.FALSE.
lsigsd=.FALSE.

do isd=1,nsdmx
   if(epssd(isd).ne.1.) then
      lepssd=.true.
   end if
   if(xmusd(isd).ne.1.) then
      lxmusd=.true.
   end if
   if(sigsd(isd).ne.0.) then
      lsigsd=.true.
   end if
end do

do ip=1,ncxmx
   if(nbpdon(ip) >  2) then
      do jp=2,nbpdon(ip)
         if(ttpdon(jp-1,ip) >= ttpdon(jp,ip)) then
            write(iout,940)
            call errout(iout,"F","initialisation_solveurs"," ")
         end if
      end do
   else if(nbpdon(ip) <  2) then
      write(iout,941)
      call errout(iout,"F","initialisation_solveurs"," ")
   end if
end do

if(nbfrnt >  0) then
   do i=1,nbfrnt
      if((z1frnt(i).ne.z2frnt(i)).and.(r1frnt(i).ne.r2frnt(i))) then
         write(iout,950)
         write(iout,935) i,z1frnt(i),z2frnt(i),r1frnt(i),r2frnt(i)
         call errout(iout,"F","initialisation_solveurs"," ")
      end if
      if((z1frnt(i) >  z2frnt(i)).OR.(r1frnt(i) >  r2frnt(i)))  then
         write(iout,951)
         write(iout,935) i,z1frnt(i),z2frnt(i),r1frnt(i),r2frnt(i)
         call errout(iout,"F","initialisation_solveurs"," ")
      end if
   end do
end if

!*** Initialisation de la structure des conditions limites
allocate(bcnd%ntypfr(nmxfr),bcnd%ntypfb(nmxfr),bcnd%ntypfp(nmxfr))
allocate(bcnd%ipfrnt(nmxfr,nesp))
allocate(bcnd%itfrnt(nmxfr), bcnd%irefnt(nmxfr))
allocate(bcnd%r1frnt(nmxfr), bcnd%r2frnt(nmxfr), bcnd%z1frnt(nmxfr),bcnd%z2frnt(nmxfr))
allocate(bcnd%numond(nmxfr), bcnd%nu2ond(nmxfr), bcnd%nu3ond(nmxfr))
allocate(bcnd%nu4ond(nmxfr), bcnd%nu5ond(nmxfr), bcnd%nu6ond(nmxfr))
allocate(bcnd%ifoond(nmxfr))
allocate(bcnd%tdbond(nmxfr), bcnd%tmtond(nmxfr), bcnd%tplond(nmxfr), bcnd%tdsond(nmxfr))
allocate(bcnd%freond(nmxfr), bcnd%phaond(nmxfr), bcnd%eleond(nmxfr))
allocate(bcnd%tauond(nmxfr), bcnd%sigond(nmxfr))
allocate(bcnd%numon2(nmxfr), bcnd%ifoon2(nmxfr))
allocate(bcnd%tdbon2(nmxfr), bcnd%tmton2(nmxfr), bcnd%tplon2(nmxfr), bcnd%tdson2(nmxfr))
allocate(bcnd%freon2(nmxfr), bcnd%phaon2(nmxfr), bcnd%eleon2(nmxfr))
allocate(bcnd%tauon2(nmxfr), bcnd%sigon2(nmxfr))
allocate(bcnd%tdbpot(nmxfr), bcnd%tmtpot(nmxfr), bcnd%tplpot(nmxfr), bcnd%tdspot(nmxfr))
allocate(bcnd%ifopot(nmxfr))
allocate(bcnd%numddp(nmxfr))
allocate(bcnd%potfr(nmxfr), bcnd%ddpfr(nmxfr), bcnd%ddpft(nmxfr),bcnd%b3fr(nmxfr))
allocate(bcnd%epssd(0:nmxsd),bcnd%xmusd(0:nmxsd),bcnd%sigsd(nmxsd))

bcnd%ntypfr = ntypfr(1:nmxfr)
bcnd%ntypfb = ntypfb(1:nmxfr)
bcnd%ntypfp = ntypfp(1:nmxfr)

bcnd%nbfrnt = nbfrnt
bcnd%ipfrnt = ipfrnt(1:nmxfr,1:nesp)
bcnd%itfrnt = itfrnt(1:nmxfr); bcnd%irefnt = irefnt(1:nmxfr)
bcnd%r1frnt = r1frnt(1:nmxfr); bcnd%r2frnt = r2frnt(1:nmxfr)
bcnd%z1frnt = z1frnt(1:nmxfr); bcnd%z2frnt = z2frnt(1:nmxfr)

bcnd%laonde = .false.; 
bcnd%numond = numond(1:nmxfr)
bcnd%nu2ond = nu2ond(1:nmxfr); bcnd%nu3ond = nu3ond(1:nmxfr)
bcnd%nu4ond = nu4ond(1:nmxfr); bcnd%nu5ond = nu5ond(1:nmxfr)
bcnd%nu6ond = nu6ond(1:nmxfr)
bcnd%ifoond = ifoond(1:nmxfr)
bcnd%tdbond = tdbond(1:nmxfr); 
bcnd%tmtond = tmtond(1:nmxfr); bcnd%tplond = tplond(1:nmxfr)
bcnd%tdsond = tdsond(1:nmxfr); bcnd%freond = freond(1:nmxfr); 
bcnd%phaond = phaond(1:nmxfr)
bcnd%eleond = eleond(1:nmxfr); bcnd%tauond = tauond(1:nmxfr); 
bcnd%sigond = sigond(1:nmxfr)

bcnd%laond2 = .false.
bcnd%numon2 = numon2(1:nmxfr); bcnd%ifoon2 = ifoon2(1:nmxfr)
bcnd%tdbon2 = tdbon2(1:nmxfr); bcnd%tmton2 = tmton2(1:nmxfr)
bcnd%tplon2 = tplon2(1:nmxfr); bcnd%tdson2 = tdson2(1:nmxfr)
bcnd%freon2 = freon2(1:nmxfr); bcnd%phaon2 = phaon2(1:nmxfr)
bcnd%eleon2 = eleon2(1:nmxfr); bcnd%tauon2 = tauon2(1:nmxfr)
bcnd%sigon2 = sigon2(1:nmxfr); bcnd%tdbpot = tdbpot(1:nmxfr)

bcnd%tmtpot = tmtpot(1:nmxfr); bcnd%tplpot = tplpot(1:nmxfr)
bcnd%tdspot = tdspot(1:nmxfr); bcnd%ifopot = ifopot(1:nmxfr)
bcnd%numddp = numddp(1:nmxfr); bcnd%potfr  = potfr (1:nmxfr)
bcnd%ddpfr  = ddpfr (1:nmxfr); bcnd%ddpft  = ddpft (1:nmxfr)
bcnd%b3fr   = b3fr  (1:nmxfr)

bcnd%epssd  = epssd(0:nmxsd); bcnd%xmusd = xmusd(0:nmxsd)
bcnd%sigsd  = sigsd(1:nmxsd)

 900  format(//10x,'Conditions aux limites sur les frontieres'/)
 901  format(/10x,'Reference :',i3,5x,'Axe de symetrie'/    &
     &        10x,'           ',3x,5x,'Reflexion des particules')
 902  format(/10x,'Reference :',i3,5x,'Conducteur parfait'/ &
     &        10x,'           ',3x,5x,'Absorption des particules')
 903  format(/10x,'Reference :',i3,5x,'Mur magnetique'/ &
     &        10x,'           ',3x,5x,'Reflexion des particules')
 904  format(/10x,'Reference :',i3,5x,'Onde sortante'/  &
     &        10x,'           ',3x,5x,'Absorption des particules')
 905  format(/10x,'Reference :',i3,5x,'Onde entrante'/  &
              10x,'           ',3x,5x,'Absorption des particules')
 906  format(//5x,'Caracteristiques de chaque sous-domaine '    &
     &       //5x,'Sous-domaine',5x,'epssd',10x,'xmusd',10x,'sigsd')
 907  format(5x,I6,F17.5,F15.5,E17.5)
 910  format(/10x,'Option non disponible'/  &
     &        10x,'Reference :',i3,5x,'Type :',i3/)
 913  format(/10x,'ifr =',i5/10x,'Erreur dans les coefs temporels'/)
 914  format(/10x,'Mauvaise valeur de sigond')
 915  format(/10x,'Mauvaise valeur de sigon2')
 916  format(/10x,'Mauvaise valeur de sigant')
 917  format(/10x,'Frontiere ',I4,' ifopot = ',I4,' est impossible')
 921  format(//5x,'Frontiere',5x,'ntypfr',7x,'potfr',7x,'ddpfr' &
     &        ,6x,'ntypfb',6x,'b3fr')
 922  format(5x,I3,4x,I10,2E12.3,2I10,E12.3)
 9221 format(//5x,'Frontiere ',I4,' ntypfr = ',I2,' est impossible' &
     &           ,' avec nsolve = ',I10)
 9222 format(//5x,'Frontiere ',I4,' ntypfr = ',I2,' est impossible' &
     &           ,' avec lmodte = ',L10)
 9223 format(//5x,'Frontiere',5x,'ntypfp', 'C.L. Correction Galerkine', &
              /5x,I3,4x,I10)
 923  format(//5x,'Variation temporelle du potentiel applique'      &
     &       //1x,'Frontiere',4x,'tdbpot',6x,'tmtpot',6x,'tplpot'   &
     &                      ,6x,'tdspot',6x,'ifopot')
 924  format(1x,I5,4x,4E12.3,I10)
 932  format(//10x,'CONDITIONS AUX LIMITES'/)
 933  format(/10x,'Erreur dans les caracteristiques des materiaux'/)
 934  format(//10x,'FRONTIERES INTERNES   nbfrnt : ',I3//)
 935  format(/5x,'No',I3,5x,'z1frnt :',E12.3,5x,'z2frnt :',E12.3    &
     &      /15x,'r1frnt :',E12.3,5x,'r2frnt :',E12.3           &
     &      /15x,'itfrnt :',I3,5x,'irefnt :',I3)
 9351 format(15x,'Espece No ',I3,5x,'ipfrnt :',I3)
 936  format(//10x,'FRONTIERES INTERNES   nbfrtp : ',I3//)
 937  format(/5x,'No',I3,5x,'ireftp :',I3)
 9371 format(/10x,'Espece No ',I3,5x,'Transparence : ltrftp = ',L10)
 938  format(//1x,'   ',3x,'freond',3x,'phaond'             &
     &        ,3x,'eleond',3x,'tdbond',3x,'tmtond',3x,'tplond'      &
     &        ,3x,'tdsond',2x,'ifoond'/)
 939  format(1x,I3,7E9.2,I5)
 9391 format(56x,'tauond=',E9.2/56x,'sigond=',E9.2)
 9392 format(56x,'tauon2=',E9.2/56x,'sigon2=',E9.2)
 940  format(//10x,'La progression des valeurs de temps pour la'    &
     &        /10x,'definition des formes temporelles doit etre'    &
     &        /10x,'strictement croissante'//)
 941  format(//10x,'Le nombre de points pour la'            &
     &        /10x,'definition des formes temporelles doit etre'    &
     &        /10x,'superieur a 1 '//)
 942  format( /10x,'Numero de forme temporelle : ',I2           &
     &        /10x,'Nombre de points : ',I5             &
     &        /10x,'La limite est fixee a 100'//)
 943  format( /10x,'Numero de forme temporelle : ',I2           &
     &        /10x,'Nombre de points : ',I3             &
     &       //10x,'    n            t                f'/)
 944  format(  10x,I10,2E15.4)
 945  format(//5x,'Definition des ONDES ENTRANTES Mode (E1,E2,B3)'  &
     &       //1x,'Numero de ref',5x,'Numeros d''onde'/)
 946  format(1X,I7,5X,6I4)
 947  format(//5x,'Definition des ONDES ENTRANTES Mode (B1,B2,E3)'  &
     &       //1x,'Numero de ref',5x,'Numero d''onde'/)
 948  format(1x,I7,I12)
 949  format(//1x,'   ',3x,'freon2',3x,'phaon2'             &
     &        ,3x,'eleon2',3x,'tdbon2',3x,'tmton2',3x,'tplon2'      &
     &        ,3x,'tdson2',2x,'ifoon2'/)
 950  format(//10x,'Erreur de frontiere interne'            &
     &        /10x,'Elle doit etre horizontale ou verticale')
 951  format(//10x,'Erreur sur les coordonnees des extremites'      &
     &        /10x,'des frontieres internes')
 952  format(/10x,'Reference :',i3,5x,'Conducteur parfait'/ &
     &        10x,'           ',3x,5x,'Reflexion des particules')
 953  format(/10x,'Reference :',i3,5x,'Onde sortante ou Neumann'/   &
     &        10x,'           ',3x,5x,'Reflexion des particules')
 954  format(/10x,'Reference :',i3,5x,'DDP imposee'/            &
     &        10x,'           ',3x,5x,'Absorption des particules')
 955  format(//5x,'Degre des polynome pour le solveur Galerkine ', i2)

end subroutine lecture_donnees_solveur

end module solveurs_module
