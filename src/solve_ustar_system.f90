subroutine solve_ustar_system()
use Lsys    
implicit none
integer, dimension(8) :: tini,tfin

call DATE_AND_TIME(values=tini)

!solving ustar linear system
if     (boundary_kind.eq.1) then
    !In this block:
    !1 - Lid shear-driven cavity flow
    !2 - Enclosures
    call PCG55_SSOR_v3(niu,nju,au(1:niu,1:nju),bu(1:niu,1:nju),cu(1:niu,1:nju),du(1:niu,1:nju),eu(1:niu,1:nju),fu(0:niu+1,0:nju+1),ustar(0:niu+1,0:nju+1),tol=CGtol,h=dx)
    call PCG55_SSOR_v3(niv,njv,av(1:niv,1:niv),bv(1:niv,1:njv),cv(1:niv,1:njv),dv(1:niv,1:njv),ev(1:niv,1:njv),fv(0:niv+1,0:njv+1),vstar(0:niv+1,0:njv+1),tol=CGtol,h=dx)
elseif ((boundary_kind.eq.2).or.(boundary_kind.eq.3)) then
    !In this block:
    !1 - Simple shear flow
    !2 - Periodic pressure-driven flow
    call PCG55_SSOR_x_per_v3(niu,nju,au(1:niu,1:nju),bu(1:niu,1:nju),cu(1:niu,1:nju),du(1:niu,1:nju),eu(1:niu,1:nju),fu(0:niu+1,0:nju+1),ustar(0:niu+1,0:nju+1),tol=CGtol,h=dx)
    call PCG55_SSOR_x_per_v3(niv,njv,av(1:niv,1:niv),bv(1:niv,1:njv),cv(1:niv,1:njv),dv(1:niv,1:njv),ev(1:niv,1:njv),fv(0:niv+1,0:njv+1),vstar(0:niv+1,0:njv+1),tol=CGtol,h=dx)
end if

call DATE_AND_TIME(values=tfin)
if(disp_exec_time) write(unit=*,fmt=1000) 'Elapsed time to solve ustar and vstar systems in seconds:',time_difference(tfin,tini)
1000 format(a60,f8.3)
return
end subroutine solve_ustar_system