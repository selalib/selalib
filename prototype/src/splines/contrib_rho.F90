
module contrib_rho_module
#include "sll_working_precision.h"
#include "sll_assert.h"
#include "sll_memory.h"

  use numeric_constants
  !use utils
  implicit none
contains  




function compute_contrib_rho(integration_points,N_int1,N_int2,size_contrib,index_contrib,value_contrib,Ntot,N1,N2)
  sll_int :: compute_contrib_rho
  sll_int,intent(in)         ::  Ntot,N1,N2,N_int1,N_int2
  sll_real64,dimension(2,N_int2+1,N_int1+1)  ::  integration_points
  sll_int,dimension(-1:N2+1,N_int1+1)    ::  size_contrib
  sll_int,dimension(Ntot)         ::  index_contrib
  sll_real64,dimension(Ntot)         :: value_contrib
  sll_int                    ::  i,j,ii,jj,i_loc,j_loc
  sll_real64                  ::x1,x2,w(1:2,-1:2),xx
  sll_int,dimension(:,:,:),pointer :: point_flag
  sll_int,dimension(:,:),pointer :: line_flag
  sll_real64,dimension(:,:),pointer ::  point_val
  sll_int :: N_size,s,err
  SLL_ALLOCATE(point_flag(1:2,-1:N1+1,-1:N2+1),err)
  SLL_ALLOCATE(line_flag(1:2,-1:N2+1),err)
  SLL_ALLOCATE(point_val(-1:N1+1,-1:N2+1),err)
  N_size = 0
  line_flag=0
  point_flag=0
  s=0
  size_contrib=0
  do i=1,N_int1+1
    !buf=0
    do j=1,N_int2+1
      x1 = integration_points(1,j,i)
      x2 = integration_points(2,j,i)
      !treatment of boundary conditions
      ii=floor(x1*real(N1,f64))
      if(x1>=1._f64)then
        x1=1._f64
        ii=N1-1
      endif
      if(x1<=0._f64)then
        x1=0._f64
        ii=0
      endif
      jj=floor(x2*real(N2,f64))
      if(x2>=1._f64)then
        x2=1._f64
        jj=N2-1
      endif
      if(x2<=0._f64)then
        x2=0._f64
        jj=0
      endif
      
      xx=(x1*real(N1,f64))-real(ii,f64)
      
      if((xx>1_f64).or.(xx<0._f64))then
        print *,'x1=',xx
        stop
      endif
      
      w(1,-1)=(1.0_f64-xx)*(1.0_f64-xx)*(1.0_f64-xx)/6._f64    
      w(1,0)=2._f64/3._f64+xx*xx*xx/2._f64-xx*xx
      w(1,1)=xx*xx/2._f64+xx/2._f64+1._f64/6._f64-xx*xx*xx/2._f64
      w(1,2)=xx*xx*xx/6._f64
      xx=(x2*real(N2,f64))-real(jj,f64)      

      if((xx>1_f64).or.(xx<0._f64))then
        print *,'x2=',xx
        stop
      endif


      w(2,-1)=(1.0_f64-xx)*(1.0_f64-xx)*(1.0_f64-xx)/6._f64    
      w(2,0)=2._f64/3._f64+xx*xx*xx/2._f64-xx*xx
      w(2,1)=xx*xx/2._f64+xx/2._f64+1._f64/6._f64-xx*xx*xx/2._f64
      w(2,2)=xx*xx*xx/6._f64
      
      do j_loc=-1,2
        do i_loc=-1,2
          if(line_flag(1,jj+j_loc)/=i)then
            line_flag(1,jj+j_loc)=i
            line_flag(2,jj+j_loc)=-1
         endif
          if(point_flag(1,ii+i_loc,jj+j_loc)/=i)then
            size_contrib(jj+j_loc,i)=size_contrib(jj+j_loc,i)+1
            N_size  = N_size +1
            point_flag(1,ii+i_loc,jj+j_loc)=i
            point_val(ii+i_loc,jj+j_loc) = 0._f64
            point_flag(2,ii+i_loc,jj+j_loc) = line_flag(2,jj+j_loc)
            line_flag(2,jj+j_loc) = ii+i_loc
          endif          
          point_val(ii+i_loc,jj+j_loc)=point_val(ii+i_loc,jj+j_loc)+w(1,i_loc)*w(2,j_loc)
        enddo
      enddo
      
    enddo
    
    !fill the contrib arrays
    !do j=-1,N2+1
    !  !print *,i,j,size_contrib(j,i),line_flag(2,j)
    !  ii=line_flag(2,j)
    !  do jj=1,size_contrib(j,i)
    !    !print *,"##",ii,point_val(ii,j)
    !    ii=point_flag(2,ii,j) 
    !  enddo      
    !enddo
    !stop
    
    do j = -1,N2+1
      if(size_contrib(j,i)>=1)then
        ii=line_flag(2,j)
      endif
      do jj=1,size_contrib(j,i)
        s=s+1
        if(s<=Ntot)then          
          index_contrib(s) = ii
          value_contrib(s) = point_val(ii,j)
        endif
        ii=point_flag(2,ii,j) 
      enddo
    enddo
    
  enddo
  
  SLL_DEALLOCATE(point_flag,err)
  SLL_DEALLOCATE(line_flag,err)
  SLL_DEALLOCATE(point_val,err)
  !deallocate(point_flag,line_flag,point_val)
  
  if(s/=N_size)then
    print *,'s=',s,'N_size=',N_size
    !stop
  endif
  
  compute_contrib_rho = s
  
end function compute_contrib_rho




end module contrib_rho_module

