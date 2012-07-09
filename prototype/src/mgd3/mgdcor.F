      subroutine mgdcor(sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf,phif,
     1                  sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc,phic,
     2                  sx1,ex1,sy1,ey1,sz1,ez1,bd,phibc,IOUT)
# include "compdir.inc"
      include "mpif.h"
      integer sxf,exf,syf,eyf,szf,ezf,nxf,nyf,nzf
      integer sxc,exc,syc,eyc,szc,ezc,nxc,nyc,nzc
      integer sx1,ex1,sy1,ey1,sz1,ez1,bd(26),IOUT
      REALN phif(sxf-1:exf+1,syf-1:eyf+1,szf-1:ezf+1)
      REALN phic(sxc-1:exc+1,syc-1:eyc+1,szc-1:ezc+1),phibc(6)
c------------------------------------------------------------------------
c Add correction from coarse grid level to fine grid level. Uses
c bilinear interpolation for the old version of the multigrid code,
c and volume weighting for its new version.
c
c Code      : mgd3, 3-D parallel multigrid solver
c Author    : Bernard Bunner (bunner@engin.umich.edu), January 1998
c Called in : mgdkcyc
c Calls     : mgdbdry
c------------------------------------------------------------------------
      integer i,j,k,ic,jc,kc,i1,i2,j1,j2,k1,k2
# if cdebug
      double precision tinitial
      tinitial=MPI_WTIME()
# endif
c
# if WMGD
c------------------------------------------------------------------------
c new version: the correction is the weighted average of either two 
c or four points at the coarser grid level depending on whether 
c coarsifying takes place in all directions or not
c
c first, general case where coarsifying takes place in all directions;
c take the averages from the 8 surrounding points
c
      if ((nxc.lt.nxf).and.(nyc.lt.nyf).and.(nzc.lt.nzf)) then
        do kc=szc-1,ezc
          k=2*kc-1
          do jc=syc-1,eyc
            j=2*jc-1
            do ic=sxc-1,exc
              i=2*ic-1
              phif(i,j,k)=phif(i,j,k)+
     1          (27.0d0*phic(ic,jc,kc)+9.0d0*phic(ic,jc+1,kc)
     2          +9.0d0*phic(ic,jc,kc+1)+3.0d0*phic(ic,jc+1,kc+1)
     3          +9.0d0*phic(ic+1,jc,kc)+3.0d0*phic(ic+1,jc+1,kc)
     4          +3.0d0*phic(ic+1,jc,kc+1)+phic(ic+1,jc+1,kc+1))
     5          /64.0d0
              phif(i+1,j,k)=phif(i+1,j,k)+
     1          (9.0d0*phic(ic,jc,kc)+3.0d0*phic(ic,jc+1,kc)
     2          +3.0d0*phic(ic,jc,kc+1)+phic(ic,jc+1,kc+1)
     3          +27.0d0*phic(ic+1,jc,kc)+9.0d0*phic(ic+1,jc+1,kc)
     4          +9.0d0*phic(ic+1,jc,kc+1)+3.0d0*phic(ic+1,jc+1,kc+1))
     5          /64.0d0
              phif(i,j+1,k)=phif(i,j+1,k)+
     1          (9.0d0*phic(ic,jc,kc)+27.0d0*phic(ic,jc+1,kc)
     2          +3.0d0*phic(ic,jc,kc+1)+9.0d0*phic(ic,jc+1,kc+1)
     3          +3.0d0*phic(ic+1,jc,kc)+9.0d0*phic(ic+1,jc+1,kc)
     4          +phic(ic+1,jc,kc+1)+3.0d0*phic(ic+1,jc+1,kc+1))
     5          /64.0d0
              phif(i,j,k+1)=phif(i,j,k+1)+
     1          (9.0d0*phic(ic,jc,kc)+3.0d0*phic(ic,jc+1,kc)
     2          +27.0d0*phic(ic,jc,kc+1)+9.0d0*phic(ic,jc+1,kc+1)
     3          +3.0d0*phic(ic+1,jc,kc)+phic(ic+1,jc+1,kc)
     4          +9.0d0*phic(ic+1,jc,kc+1)+3.0d0*phic(ic+1,jc+1,kc+1))
     5          /64.0d0
              phif(i+1,j+1,k)=phif(i+1,j+1,k)+
     1          (3.0d0*phic(ic,jc,kc)+9.0d0*phic(ic,jc+1,kc)
     2          +phic(ic,jc,kc+1)+3.0d0*phic(ic,jc+1,kc+1)
     3          +9.0d0*phic(ic+1,jc,kc)+27.0d0*phic(ic+1,jc+1,kc)
     4          +3.0d0*phic(ic+1,jc,kc+1)+9.0d0*phic(ic+1,jc+1,kc+1))
     5          /64.0d0
              phif(i+1,j,k+1)=phif(i+1,j,k+1)+
     1          (3.0d0*phic(ic,jc,kc)+phic(ic,jc+1,kc)
     2          +9.0d0*phic(ic,jc,kc+1)+3.0d0*phic(ic,jc+1,kc+1)
     3          +9.0d0*phic(ic+1,jc,kc)+3.0d0*phic(ic+1,jc+1,kc)
     4          +27.0d0*phic(ic+1,jc,kc+1)+9.0d0*phic(ic+1,jc+1,kc+1))
     5          /64.0d0
              phif(i,j+1,k+1)=phif(i,j+1,k+1)+
     1          (3.0d0*phic(ic,jc,kc)+9.0d0*phic(ic,jc+1,kc)
     2          +9.0d0*phic(ic,jc,kc+1)+27.0d0*phic(ic,jc+1,kc+1)
     3          +phic(ic+1,jc,kc)+3.0d0*phic(ic+1,jc+1,kc)
     4          +3.0d0*phic(ic+1,jc,kc+1)+9.0d0*phic(ic+1,jc+1,kc+1))
     5          /64.0d0
              phif(i+1,j+1,k+1)=phif(i+1,j+1,k+1)+
     1          (phic(ic,jc,kc)+3.0d0*phic(ic,jc+1,kc)
     2          +3.0d0*phic(ic,jc,kc+1)+9.0d0*phic(ic,jc+1,kc+1)
     3          +3.0d0*phic(ic+1,jc,kc)+9.0d0*phic(ic+1,jc+1,kc)
     4          +9.0d0*phic(ic+1,jc,kc+1)+27.0d0*phic(ic+1,jc+1,kc+1))
     5          /64.0d0
            end do
          end do
        end do
      else
c
c no coarsifying in two directions; take the averages from the 2 
c surrounding points
c
        if ((nxf.eq.nxc).and.(nyf.eq.nyc)) then
          do kc=szc-1,ezc
            k=2*kc-1
            do j=syf-1,eyf+1
              jc=j
              do i=sxf-1,exf+1
                ic=i
                phif(i,j,k)=phif(i,j,k)+
     1            (3.0d0*phic(ic,jc,kc)+phic(ic,jc,kc+1))/4.0d0
                phif(i,j,k+1)=phif(i,j,k+1)+
     1            (phic(ic,jc,kc)+3.0d0*phic(ic,jc,kc+1))/4.0d0
              end do
            end do
          end do
        else if ((nxf.eq.nxc).and.(nzf.eq.nzc)) then
          do k=szf-1,ezf+1
            kc=k
            do jc=syc-1,eyc
              j=2*jc-1
              do i=sxf-1,exf+1
                ic=i
                phif(i,j,k)=phif(i,j,k)+
     1            (3.0d0*phic(ic,jc,kc)+phic(ic,jc+1,kc))/4.0d0
                phif(i,j+1,k)=phif(i,j+1,k)+
     1            (phic(ic,jc,kc)+3.0d0*phic(ic,jc+1,kc))/4.0d0
              end do
            end do
          end do
        else if ((nyf.eq.nyc).and.(nzf.eq.nzc)) then
          do k=szf-1,ezf+1
            kc=k
            do j=syf-1,eyf+1
              jc=j
              do ic=sxc-1,exc
                i=2*ic-1
                phif(i,j,k)=phif(i,j,k)+
     1            (3.0d0*phic(ic,jc,kc)+phic(ic+1,jc,kc))/4.0d0
                phif(i+1,j,k)=phif(i+1,j,k)+
     1            (phic(ic,jc,kc)+3.0d0*phic(ic+1,jc,kc))/4.0d0
              end do
            end do
          end do
c
c no coarsifying in one direction; take the averages from the 4
c surrounding points
c
        else if (nxf.eq.nxc) then
          do kc=szc-1,ezc
            k=2*kc-1
            do jc=syc-1,eyc
              j=2*jc-1
              do i=sxf-1,exf+1
                ic=i
                phif(i,j,k)=phif(i,j,k)+
     1            (9.0d0*phic(ic,jc,kc)+3.0d0*phic(ic,jc+1,kc)
     2            +3.0d0*phic(ic,jc,kc+1)+phic(ic,jc+1,kc+1))
     3            /16.0d0
                phif(i,j+1,k)=phif(i,j+1,k)+
     1            (3.0d0*phic(ic,jc,kc)+9.0d0*phic(ic,jc+1,kc)
     2            +phic(ic,jc,kc+1)+3.0d0*phic(ic,jc+1,kc+1))
     3            /16.0d0
                phif(i,j,k+1)=phif(i,j,k+1)+
     1            (3.0d0*phic(ic,jc,kc)+phic(ic,jc+1,kc)
     2            +9.0d0*phic(ic,jc,kc+1)+3.0d0*phic(ic,jc+1,kc+1))
     3            /16.0d0
                phif(i,j+1,k+1)=phif(i,j+1,k+1)+
     1            (phic(ic,jc,kc)+3.0d0*phic(ic,jc+1,kc)
     2            +3.0d0*phic(ic,jc,kc+1)+9.0d0*phic(ic,jc+1,kc+1))
     3            /16.0d0
              end do
            end do
          end do
        else if (nyf.eq.nyc) then
          do kc=szc-1,ezc
            k=2*kc-1
            do j=syf-1,eyf+1
              jc=j
              do ic=sxc-1,exc
                i=2*ic-1
                phif(i,j,k)=phif(i,j,k)+
     1            (9.0d0*phic(ic,jc,kc)+3.0d0*phic(ic+1,jc,kc)
     2            +3.0d0*phic(ic,jc,kc+1)+phic(ic+1,jc,kc+1))
     3            /16.0d0
                phif(i+1,j,k)=phif(i+1,j,k)+
     1            (3.0d0*phic(ic,jc,kc)+9.0d0*phic(ic+1,jc,kc)
     2            +phic(ic,jc,kc+1)+3.0d0*phic(ic+1,jc,kc+1))
     3            /16.0d0
                phif(i,j,k+1)=phif(i,j,k+1)+
     1            (3.0d0*phic(ic,jc,kc)+phic(ic+1,jc,kc)
     2            +9.0d0*phic(ic,jc,kc+1)+3.0d0*phic(ic+1,jc,kc+1))
     3            /16.0d0
                phif(i+1,j,k+1)=phif(i+1,j,k+1)+
     1            (phic(ic,jc,kc)+3.0d0*phic(ic+1,jc,kc)
     2            +3.0d0*phic(ic,jc,kc+1)+9.0d0*phic(ic+1,jc,kc+1))
     3            /16.0d0
              end do
            end do
          end do
        else if (nzf.eq.nzc) then
          do k=szf-1,ezf+1
            kc=k
            do jc=syc-1,eyc
              j=2*jc-1
              do ic=sxc-1,exc
                i=2*ic-1
                phif(i,j,k)=phif(i,j,k)+
     1            (9.0d0*phic(ic,jc,kc)+3.0d0*phic(ic+1,jc,kc)
     2            +3.0d0*phic(ic,jc+1,kc)+phic(ic+1,jc+1,kc))
     3            /16.0d0
                phif(i+1,j,k)=phif(i+1,j,k)+
     1            (3.0d0*phic(ic,jc,kc)+9.0d0*phic(ic+1,jc,kc)
     2            +phic(ic,jc+1,kc)+3.0d0*phic(ic+1,jc+1,kc))
     3            /16.0d0
                phif(i,j+1,k)=phif(i,j+1,k)+
     1            (3.0d0*phic(ic,jc,kc)+phic(ic+1,jc,kc)
     2            +9.0d0*phic(ic,jc+1,kc)+3.0d0*phic(ic+1,jc+1,kc))
     3            /16.0d0
                phif(i+1,j+1,k)=phif(i+1,j+1,k)+
     1            (phic(ic,jc,kc)+3.0d0*phic(ic+1,jc,kc)
     2            +3.0d0*phic(ic,jc+1,kc)+9.0d0*phic(ic+1,jc+1,kc))
     3            /16.0d0
              end do
            end do
          end do
        end if
      end if
c
c impose Neumann and Dirichlet boundary conditions
C TEMP: periodicity is not enforced to save one call to gxch1lin;
c check whether it has an impact or not...
c
      call mgdbdry(sxf,exf,syf,eyf,szf,ezf,phif,bd,phibc,IOUT)
# else
c------------------------------------------------------------------------
c old version
c
      if (nxc.lt.nxf) then
        i1=1
        i2=0
      else
        i1=0
        i2=1
      end if
      if (nyc.lt.nyf) then
        j1=1
        j2=0
      else
        j1=0
        j2=1
      end if
      if (nzc.lt.nzf) then
        k1=1
        k2=0
      else
        k1=0
        k2=1
      end if
c
c identity at the points of the fine grid which have odd indices 
c in i, j, and k
c
      do kc=sz1,ez1
        k=k1*(2*kc-1)+k2*kc
        do jc=sy1,ey1
          j=j1*(2*jc-1)+j2*jc
          do ic=sx1,ex1
            i=i1*(2*ic-1)+i2*ic
            phif(i,j,k)=phif(i,j,k)+phic(ic,jc,kc)
          end do
        end do
      end do
c
c interpolation of the two neighboring values for the points of the 
c fine grid with even index for i and odd indices for j and k
c
      if (nxc.lt.nxf) then
        do kc=sz1,ez1
          k=k1*(2*kc-1)+k2*kc
          do jc=sy1,ey1
            j=j1*(2*jc-1)+j2*jc
            do ic=sxc-1,exc
              i=2*ic
              phif(i,j,k)=phif(i,j,k)
     1                   +0.5d0*(phic(ic,jc,kc)+phic(ic+1,jc,kc))
            end do
          end do
        end do
      end if
c
c interpolation of the two neighboring values for the points of the 
c fine grid with even index for j and odd indices for i and k
c
      if (nyc.lt.nyf) then
        do kc=sz1,ez1
          k=k1*(2*kc-1)+k2*kc
          do jc=syc-1,eyc
            j=2*jc
            do ic=sx1,ex1
              i=i1*(2*ic-1)+i2*ic
              phif(i,j,k)=phif(i,j,k)
     1                   +0.5d0*(phic(ic,jc,kc)+phic(ic,jc+1,kc))
            end do
          end do
        end do
      end if
c
c interpolation of the two neighboring values for the points of the
c fine grid with even index for k and odd indices for i and j
c
      if (nzc.lt.nzf) then
        do kc=szc-1,ezc
          k=2*kc
          do jc=sy1,ey1
            j=j1*(2*jc-1)+j2*jc
            do ic=sx1,ex1
              i=i1*(2*ic-1)+i2*ic
              phif(i,j,k)=phif(i,j,k)
     1                   +0.5d0*(phic(ic,jc,kc)+phic(ic,jc,kc+1))
            end do
          end do
        end do
      end if
c
c interpolation of the four neighboring values for the points of the
c fine grid with even indices for i and j and odd index for k
c
      if (nxc.lt.nxf.and.nyc.lt.nyf) then
        do kc=sz1,ez1
          k=k1*(2*kc-1)+k2*kc
          do jc=syc-1,eyc
            j=2*jc
            do ic=sxc-1,exc
              i=2*ic
              phif(i,j,k)=phif(i,j,k)
     1                   +0.25d0*(phic(ic,jc,kc)+phic(ic+1,jc,kc)
     2                           +phic(ic,jc+1,kc)+phic(ic+1,jc+1,kc))
            end do
          end do
        end do
      end if
c
c interpolation of the four neighboring values for the points of the
c fine grid with even indices for i and k and odd index for j
c
      if (nxc.lt.nxf.and.nzc.lt.nzf) then
        do kc=szc-1,ezc
          k=2*kc
          do jc=sy1,ey1
            j=j1*(2*jc-1)+j2*jc
            do ic=sxc-1,exc
              i=2*ic
              phif(i,j,k)=phif(i,j,k)
     1                   +0.25d0*(phic(ic,jc,kc)+phic(ic+1,jc,kc)
     2                           +phic(ic,jc,kc+1)+phic(ic+1,jc,kc+1))
            end do
          end do
        end do
      end if
c
c interpolation of the four neighboring values for the points of the
c fine grid with even indices for j and k and odd index for i
c
      if (nyc.lt.nyf.and.nzc.lt.nzf) then
        do kc=szc-1,ezc
          k=2*kc
          do jc=syc-1,eyc
            j=2*jc
            do ic=sx1,ex1
              i=i1*(2*ic-1)+i2*ic
              phif(i,j,k)=phif(i,j,k)
     1                   +0.25d0*(phic(ic,jc,kc)+phic(ic,jc+1,kc)
     2                           +phic(ic,jc,kc+1)+phic(ic,jc+1,kc+1))
            end do
          end do
        end do
      end if
c
c interpolation of the eight neighboring values for the points of the
c fine grid with even indices for i, j, and k
c
      if (nxc.lt.nxf.and.nyc.lt.nyf.and.nzc.lt.nzf) then
        do kc=szc-1,ezc
          k=2*kc
          do jc=syc-1,eyc
            j=2*jc
            do ic=sxc-1,exc
              i=2*ic
              phif(i,j,k)=phif(i,j,k)
     1               +0.125d0*(phic(ic,jc,kc)+phic(ic+1,jc,kc)
     2                        +phic(ic,jc+1,kc)+phic(ic,jc,kc+1)
     3                        +phic(ic+1,jc+1,kc)+phic(ic+1,jc,kc+1)
     4                        +phic(ic,jc+1,kc+1)+phic(ic+1,jc+1,kc+1))
            end do
          end do
        end do
      end if
# endif
c
# if cdebug
      timing(92)=timing(92)+MPI_WTIME()-tinitial
# endif
      return
      end
