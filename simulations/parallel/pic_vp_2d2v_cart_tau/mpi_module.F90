!File: Module Mpi_Module
!Module contenant les routines de parallelisation
module mpi_module 

use mpi
use zone

implicit none

integer, dimension(MPI_STATUS_SIZE) :: stat

integer :: nproc=1, rang=0
integer :: code, tag=1111

contains

!*********************************************************************

!Subroutine: init_mpi
!Initialisation de MPI
subroutine init_mpi( prank , psize)
integer :: prank
integer :: psize

call MPI_INIT(code)
call MPI_COMM_RANK(MPI_COMM_WORLD,prank,code)
call MPI_COMM_SIZE(MPI_COMM_WORLD,psize,code)
print*, ' Hello from mpi proc number ',prank, ' of ', psize
call MPI_BARRIER(MPI_COMM_WORLD,code)

end subroutine init_mpi

!*********************************************************************

!Subroutine: finish_mpi
!Arret de MPI
subroutine finish_mpi()

call MPI_BARRIER(MPI_COMM_WORLD,code)
call MPI_FINALIZE(code)

end subroutine finish_mpi

!*********************************************************************

subroutine work_mpi(nx, ny)
integer :: nx
integer :: ny


!integer, intent(in)   :: istep, nstep, nesp
!integer, dimension(:) :: nplus
!type (particle), dimension(:) :: spc
!type (mesh_fields), intent(inout) :: mxw
!type (ext_fields), intent(in) ::  en, bn
!double precision, dimension(:), intent(inout) :: rho
!type(mesh_data),    intent(inout) :: mesh
!integer              :: type_field, valeur, somme, ltab
!integer,dimension(3) :: types,lbloc,adr,dep
!integer :: type_ligne, ideb, ifin, ip, nbp, iesp, itmp, ifro, ict
!real(8) :: cour, past, cour_tot
!type(objet_fonction), dimension(:) :: chx
!character(len=*), intent(in) :: dir

call MPI_BCAST(nx,1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
call MPI_BCAST(ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,code)

!if (rang == 0) then
!
!      do ip = 1, nproc-1
!
!         ideb = (ip-1)*nplus(iesp)/(nproc-1)+1
!         ifin =  ip   *nplus(iesp)/(nproc-1)
!         nbp  = ifin - ideb + 1
!
!         !write(*,2000) rang, nbp, ideb, ifin
!         call MPI_SEND(nbp,1,MPI_INTEGER,ip,tag,MPI_COMM_WORLD,code)
!
!         if (nbp > 0 ) then
!
!            call MPI_TYPE_VECTOR(nbp,3,3,MPI_REAL8,type_ligne,code)
!            call MPI_TYPE_COMMIT(type_ligne,code)
!
!            call MPI_SEND(spc(iesp)%pos(1,ideb), 1,type_ligne,  &
!                          ip,tag,MPI_COMM_WORLD,code)
!            call MPI_SEND(spc(iesp)%vit(1,ideb), 1,type_ligne,  &
!                          ip,tag,MPI_COMM_WORLD,code)
!            call MPI_SEND(spc(iesp)%xlm(1,ideb), 1,type_ligne,  &
!                          ip,tag,MPI_COMM_WORLD,code)
!            call MPI_SEND(spc(iesp)%poid(ideb),nbp,MPI_REAL8,   &
!                          ip,tag,MPI_COMM_WORLD,code)
!            call MPI_SEND(spc(iesp)%nlpa(ideb),nbp,MPI_INTEGER, &
!                          ip,tag,MPI_COMM_WORLD,code)
!
!            call MPI_TYPE_FREE(type_ligne,code)
!
!         end if
!
!   end do
!
!else 
!
!      ideb  = spc(iesp)%nbpart + 1
!      call MPI_RECV(nbp,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,stat,code)
!      spc(iesp)%nbpart = spc(iesp)%nbpart + nbp
!
!      if (spc(iesp)%nbpart > nbpam(iesp)) then
!         call MPI_FINALIZE(code)
!         stop "mpi_module - nombre de particules trop grand"
!      end if
!
!      if (nbp > 0 ) then
!
!         call MPI_TYPE_VECTOR(nbp,3,3,MPI_REAL8,type_ligne,code)
!         call MPI_TYPE_COMMIT(type_ligne,code)
!
!         call MPI_RECV(spc(iesp)%pos(1,ideb), 1,type_ligne, &
!                       0,tag,MPI_COMM_WORLD,stat,code)
!         call MPI_RECV(spc(iesp)%vit(1,ideb), 1,type_ligne, &
!                       0,tag,MPI_COMM_WORLD,stat,code)
!         call MPI_RECV(spc(iesp)%xlm(1,ideb), 1,type_ligne, &
!                   0,tag,MPI_COMM_WORLD,stat,code)
!         call MPI_RECV(spc(iesp)%poid(ideb),nbp,MPI_REAL8,  &
!                   0,tag,MPI_COMM_WORLD,stat,code)
!         call MPI_RECV(spc(iesp)%nlpa(ideb),nbp,MPI_INTEGER,    &
!                   0,tag,MPI_COMM_WORLD,stat,code)
!
!         call MPI_TYPE_FREE(type_ligne,code)
!         nplus(iesp) = nbp
!
!      end if
!
!endif
!
!call MPI_BCAST(mxw%e(1,:),mesh%nbs,MPI_REAL8,0,MPI_COMM_WORLD,code)
!call MPI_BCAST(mxw%e(2,:),mesh%nbs,MPI_REAL8,0,MPI_COMM_WORLD,code)
!call MPI_BCAST(mxw%b(3,:),mesh%nbs,MPI_REAL8,0,MPI_COMM_WORLD,code)
!
!!*** Calcul du nombre de particules dans le cas parallele
!
!if (rang == 0 ) then
!   do iesp = 1, nesp
!      spc(iesp)%nbpart = 0
!      !write(iout,"(/5x,' Espece ',i3)")iesp
!      do ip = 1, nproc-1
!         call MPI_RECV(nbp,1,MPI_INTEGER,ip,tag,MPI_COMM_WORLD,stat,code)
!	 !write(iout,"(5x,i7,' particules sur le processeur ',i3)")nbp,ip
!	 nbpproc(ip+1) = nbp
!         spc(iesp)%nbpart = spc(iesp)%nbpart + nbp
!      end do
!   end do
!else
!   do iesp = 1, nesp
!      nbp = spc(iesp)%nbpart
!      call MPI_SEND(nbp,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,code)
!   end do
!end if 
!
!call MPI_BCAST(nbpproc, nproc,MPI_INTEGER,0,MPI_COMM_WORLD,code)
!
!if( rang == 0 ) then
!
!   mxw%j(:,:) = 0.0; rho(:) = 0.0
!      
!   do ip = 1, nproc-1
!
!      call MPI_RECV(tmp, nnoeuds,MPI_REAL8,ip,tag,MPI_COMM_WORLD,stat,code)
!      rho(:) = rho(:) + tmp(:)
!
!   end do
!
!else
!
!   call MPI_SEND(rho,   nnoeuds,MPI_REAL8,0,tag,MPI_COMM_WORLD,code)
!
!end if

end subroutine work_mpi

!*********************************************************************

subroutine mpi_global_master()

call MPI_BCAST(nx,    1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
call MPI_BCAST(ny,    1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
call MPI_BCAST(nstep, 1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
call MPI_BCAST(nbpart,1,MPI_INTEGER,0,MPI_COMM_WORLD,code)
call MPI_BCAST(dimx,  1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(dimy,  1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(dx,    1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(dy,    1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(alpha, 1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(kx,    1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(ky,    1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(poids, 1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)
call MPI_BCAST(dt,    1,MPI_REAL8  ,0,MPI_COMM_WORLD,code)

end subroutine mpi_global_master

end module mpi_module
