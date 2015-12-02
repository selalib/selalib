#ifndef _particle_representations_h_
#define _particle_representations_h_

! ------   ATTENTION   ------   !
!  the 'real' offset (i.e. x - x_j where x_j is the closest node, at left of x)
!  of a particle is  offset*dx.
#define COMPUTE_CELL_OFFSET(x,xmin,rdelta,icell,offset,tmp) \
  do; \
    tmp = (x - xmin)*rdelta; \
    icell  = int(tmp); \
    offset = real(tmp - icell,f32); \
    exit; \
 end do

! 'extended' cell offset returns arbitrary cell index, allows to track particles not within domain
#define COMPUTE_CELL_OFFSET_EXTENDED(x,xmin,rdelta,icell,offset,tmp) \
  do; \
    tmp = (x - xmin)*rdelta; \
    icell  = int(floor(tmp)); \
    offset = tmp - real(icell,f64); \
    exit; \
 end do


! x, y = particle positions (in variable)
! p = particle (out variable)
! xmin, ymin, ncx, rdx, rdy = in variables
! tmp1, tmp2, off_x, off_y  = temporary real variables (tmp are real64, off are real32)  
! ic_x, ic_y                = temporary integer variables 
#define SET_PARTICLE_POSITION(p,xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2) \
 do; \
    COMPUTE_CELL_OFFSET(x,xmin,rdx,ic_x,off_x,tmp1); \
    COMPUTE_CELL_OFFSET(y,ymin,rdy,ic_y,off_y,tmp2); \
    p%ic = ic_x+1+ic_y*ncx; \
    p%dx = off_x; \
    p%dy = off_y; \
    exit; \
 end do

! similar to SET_PARTICLE_POSITION but 'extended' in the sense that cell index is arbitrary, to track particles not within domain
! also, now the (LTP) particles store one cell index per dimension: ic_x and ic_y
! x, y = particle positions (in variable)
! p = particle (out variable)
! xmin, ymin, ncx, rdx, rdy = in variables
! tmp1, tmp2, off_x, off_y  = temporary real variables (tmp are real64, off are real32)
! ic_x, ic_y                = temporary integer variables
#define SET_PARTICLE_POSITION_EXTENDED(p,xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2) \
 do; \
    COMPUTE_CELL_OFFSET_EXTENDED(x,xmin,rdx,ic_x,off_x,tmp1); \
    COMPUTE_CELL_OFFSET_EXTENDED(y,ymin,rdy,ic_y,off_y,tmp2); \
    p%ic_x = ic_x; \
    p%ic_y = ic_y; \
    p%dx = off_x; \
    p%dy = off_y; \
    exit; \
 end do

#define SET_PARTICLE_VELOCITY(p,vx,vy) \
 do; \
    p%vx = vx; \
    p%vy = vy; \
    exit; \
 end do

#define SET_PARTICLE_VALUES(p,x,y,vx,vy,qq,xmin,ymin,ncx,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2) \
 do; \
    SET_PARTICLE_POSITION(p,xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2); \
    SET_PARTICLE_VELOCITY(p,vx,vy); \
    p%q  = qq; \
    exit; \
 enddo


#define GET_PARTICLE_POSITION(p,m2d,x,y) \
  do; \
     x = m2d%delta_eta1*(real(p%dx,f64) + real( mod(p%ic-1,m2d%num_cells1), f64) ); \
     y = m2d%delta_eta2*(real(p%dy,f64) + real( int( (p%ic-1)/m2d%num_cells1 ), f64)); \
    exit; \
 end do

#define GET_PARTICLE_POSITION_EXTENDED(p,m2d,x,y) \
  do; \
     x = m2d%eta1_min + m2d%delta_eta1*(real(p%dx,f64) + real( modulo(p%ic_x-1,m2d%num_cells1), f64) ); \
     y = m2d%eta2_min + m2d%delta_eta2*(real(p%dy,f64) + real( modulo(p%ic_y-1,m2d%num_cells2), f64) ); \
    exit; \
 end do

#define SET_2DPARTICLE_VALUES(p,x,y,qq,xmin,ymin,ncx,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2) \
 do; \
    SET_PARTICLE_POSITION(p,xmin,ymin,ncx,x,y,ic_x,ic_y,off_x,off_y,rdx,rdy,tmp1,tmp2); \
    p%q  = qq; \
    exit; \
 enddo

#endif
