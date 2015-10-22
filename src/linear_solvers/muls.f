      SUBROUTINE MULKU(VKGS,VKGD,VKGI,KLD,VFG,NEQ,NSYM,VRES,nsky)       !MULK   1
C=======================================================================MULK   2
C     CE SOUS-PROGRAMME AJOUTE AU VECTEUR RES LE PRODUIT DE LA          MULK   3
C     MATRICE KG PAR LE VECTEUR FG                                      MULK   4
C       ENTREES                                                         MULK   5
C          VKGS,VKGD,VKGI  MATRICE KG STOCKEE PAR LIGNE DE CIEL (SYM.   MULK   6
C                          OU NON SYM.)                                 MULK   7
C          KLD     TABLE DES POINTEURS DES HAUTS DE COLONNES DE KG      MULK   8
C          VFG     VECTEUR FG                                           MULK   9
C          NEQ     DIMENSION DES VECTEURS FG ET RES                     MULK  10
C          NSYM    .EQ.1 SI LE PROBLEME N'EST PAS SYMETRIQUE            MULK  11
C          VRES    VECTEUR RES                                          MULK  12
C       SORTIE                                                          MULK  13
C          VRES    VECTEUR RES                                          MULK  14
C=======================================================================MULK  15
      IMPLICIT REAL(8)(A-H,O-Z)                                         !MULK  16
      DIMENSION VKGS(nsky),VKGD(neq),VKGI(nsky),
     &  KLD(neq+1),VFG(neq),VRES(neq)                                   !MULK  17
C-----------------------------------------------------------------------MULK  18
C-------  POUR CHAQUE COLONNE DE LA MATRICE KG                          MULK  19
      DO 20 IK=1,NEQ                                                    !MULK  20
      JHK=KLD(IK)                                                       !MULK  21
      JHK1=KLD(IK+1)                                                    !MULK  22
      LHK=JHK1-JHK                                                      !MULK  23
C-------  TERME DIAGONAL                                                !MULK  24
      C=VKGD(IK)*VFG(IK)                                                !MULK  25
      IF(LHK.LE.0) GO TO 20                                             !MULK  26
      I0=IK-LHK                                                         !MULK  27
C-------  TERMES DE LIGNE                                               MULK  28
      IF(NSYM.NE.1) C=C+SCAL(VKGS(JHK),VFG(I0),LHK)                     !MULK  29
      IF(NSYM.EQ.1) C=C+SCAL(VKGI(JHK),VFG(I0),LHK)                     !MULK  30
C-------  TERMES DE COLONNE                                             MULK  31
      J=JHK                                                             !MULK  32
      I1=IK-1                                                           !MULK  33
      DO 10 IJ=I0,I1                                                    !MULK  34
      VRES(IJ)=VRES(IJ)+VKGS(J)*VFG(IK)                                 !MULK  35
10    J=J+1                                                             !MULK  36
20    VRES(IK)=VRES(IK)+C                                               !MULK  37
      RETURN                                                            !MULK  38
      END                                                               !MULK  39

