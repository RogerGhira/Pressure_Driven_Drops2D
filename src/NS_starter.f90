subroutine NS_starter()
use Lsys !It's necessary to allow reading of CGtol
implicit none

integer :: iostatus  !I/O status: 0 for sucessful
integer :: i,j !loop counters
real(8) :: xlocal,ylocal
character(len=filn_len) :: filename !I/O file name

namelist /grid_dim/ xmin, ymin, Lx, Ly, nx, ny !namelist to easy data reading from file
namelist /flow_data/ Re
namelist /num_data/ dt, CGtol, Tmax
namelist /control_setup/ verb,bck_sv_frq, plt_sv_frq, bck_read_flag, bck_file, single_bck_file, &
&                        first_it, disp_exec_time
namelist /boundary_setup/ boundary_kind, initial_kind

!Reading data from file++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
filename='Drops2D_2020.dat'
open(unit=1, file=filename, status='old', iostat = iostatus)

if(iostatus.eq.0) then
    read(unit=1, nml=grid_dim) !reading directly from the namelist. Like magic!!!
    read(unit=1, nml=flow_data)
    read(unit=1, nml=num_data)
    read(unit=1, nml=control_setup)
    read(unit=1, nml=boundary_setup)
    close(unit=1,status='keep')
else 
    !Impossible to open file. Tell user.
    write(unit=*,fmt=9001) filename,iostatus
    stop
    9001 format('Imposible to open "',a60,'". Status=',i7,/,'PROGRAM STOPPED!')
end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!Setting all sizes+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if     (boundary_kind.eq.1) then
    call set_dirichlet_sizes()
elseif ((boundary_kind.eq.2).or.(boundary_kind.eq.3)) then
    call set_x_periodic_sizes()
else
    !Impossible to decide for a BC. Tell user.
    write(unit=*,fmt=9002)
    stop
9002 format('Set an appropriate boundary condition in the data file and run again!')
end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!Memory allocation+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
call alloc_memory_for_LSys(max(niu,niv,nip),max(nju,njv,njp))

!grid
allocate(x(-3:nx+3,-3:ny+3))
allocate(y(-3:nx+3,-3:ny+3))

!flow variables
allocate(u(-2:niu+3,-2:nju+3)) !NOTE: Wider vectors are necessary to allow straight ENO2 calculations with periodic condictions
allocate(v(-2:niv+3,-2:njv+3)) !NOTE: Wider vectors are necessary to allow straight ENO2 calculations with periodic condictions
allocate(ugradu(0:niu+1,0:nju+1,0:1)) !0 means n-1; 1 means n;
allocate(ugradv(0:niv+1,0:njv+1,0:1)) !0 means n-1; 1 means n;
allocate(press(0:nip+1,0:njp+1))
allocate(divustar(0:nip+1,0:njp+1))
allocate(divu(0:nip+1,0:njp+1))

!projection method variables
allocate(uinitial(-2:niu+3,-2:nju+3))
allocate(vinitial(-2:niv+3,-2:njv+3))
allocate(ucenter(-2:nip+3,-2:njp+3))
allocate(vcenter(-2:nip+3,-2:njp+3))
allocate(pressinitial(0:nip+1,0:njp+1))
allocate(ustar(0:niu+1,0:nju+1))
allocate(vstar(0:niv+1,0:njv+1))
allocate(chi(0:nip+1,0:njp+1,0:1)) !0 means n-1; 1 means n;

!ENO and WENO variables
allocate(D1x(-1:max(niu,niv,nip)+1,1:max(nju,njv,njp))) !NOTE: D1x(1,j) means D1x(3/2,j+1/2); The same for D1y. D1 is a forward difference
allocate(D1y(1:max(niu,niv,nip),-1:max(nju,njv,njp)+1))
allocate(D2x(0:max(niu,niv,nip)+1,1:max(nju,njv,njp)))
allocate(D2y(1:max(niu,niv,nip),0:max(nju,njv,njp)+1))

!Finite differences variables --> five points stencil
allocate(au(1:niu,1:nju),av(1:niv,1:njv),achi(1:nip,1:njp)) 
allocate(bu(1:niu,1:nju),bv(1:niv,1:njv),bchi(1:nip,1:njp)) 
allocate(cu(1:niu,1:nju),cv(1:niv,1:njv),cchi(1:nip,1:njp)) 
allocate(du(1:niu,1:nju),dv(1:niv,1:njv),dchi(1:nip,1:njp)) 
allocate(eu(1:niu,1:nju),ev(1:niv,1:njv),echi(1:nip,1:njp)) 
allocate(fu(0:niu+1,0:nju+1),fv(0:niv+1,0:njv+1),fchi(0:nip+1,0:njp+1)) 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!Mesh creation+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dx=Lx/nx
dy=Ly/ny

do j=-3,ny+3
    do i=-3,nx+3
        x(i,j)=xmin+i*dx
        y(i,j)=ymin+j*dy
    end do
end do

!time step setup
if(dt.eq.0.0d0) dt=dx/4.0d0
nit=ceiling(Tmax/dt)
time=(first_it-1)*dt
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!safety initialization+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
u=0.0d0; v=0.0d0; ugradu=0.0d0; ugradv=0.0d0; press=0.0d0; divu=0.0d0;
chi=0.0d0; divustar=0.0d0; ustar=0.0d0; vstar=0.0d0;
uinitial=0.0d0; vinitial=0.0d0; pressinitial=0.0d0; ucenter=0.0d0; vcenter=0.0d0;
au=0.0d0; bu=0.0d0; cu=0.0d0; du=0.0d0; eu=0.0d0; fu=0.0d0;
av=0.0d0; bv=0.0d0; cv=0.0d0; dv=0.0d0; ev=0.0d0; fv=0.0d0;
achi=0.0d0; bchi=0.0d0; cchi=0.0d0; dchi=0.0d0; echi=0.0d0; fchi=0.0d0;
D1x=0.0d0; D1y=0.0d0; D2x=0.0d0; D2y=0.0d0;
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!starting +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(not(bck_read_flag)) then
    if(initial_kind.eq.1) then
        call set_initial_flow_field()
        u=uinitial
        v=vinitial
        press=pressinitial
    end if
end if 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!restarting from backup++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(bck_read_flag) then
    write(unit=*,fmt=*)'Resuming simulation from backup file....'
    open(unit=1, file=bck_file, status='old', iostat = iostatus, form='unformatted')
    if(iostatus.eq.0) then
        read(unit=1)u
        read(unit=1)v
        read(unit=1)press
        close(unit=1,status='keep')
    else 
        !Impossible to open file. Tell user.
        write(unit=*,fmt=9001) bck_file,iostatus
        stop
    end if
end if 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

return
contains
!==================================================================================================================================
subroutine set_dirichlet_sizes()
implicit none

niu=nx-1
nju=ny
niv=nx
njv=ny-1
nip=nx
njp=ny

return
end subroutine set_dirichlet_sizes
!==================================================================================================================================

!==================================================================================================================================
subroutine set_x_periodic_sizes()
implicit none

niu=nx
nju=ny
niv=nx
njv=ny-1
nip=nx
njp=ny

return
end subroutine set_x_periodic_sizes
!==================================================================================================================================

end subroutine NS_starter