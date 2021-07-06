subroutine solve_chi_system(number_of_chi_it)
use Lsys    
implicit none
integer, optional :: number_of_chi_it
integer, dimension(8) :: tini,tfin

call DATE_AND_TIME(values=tini)

!solving chi linear system
if     (boundary_kind.eq.1) then
    !In this block:
    !1 - Lid shear-driven cavity flow
    !2 - Enclosures
    call PCG55_SSOR_v3(nip,njp,achi(1:nip,1:njp),bchi(1:nip,1:njp),cchi(1:nip,1:njp),dchi(1:nip,1:njp),echi(1:nip,1:njp),fchi(0:nip+1,0:njp+1),chi(0:nip+1,0:njp+1,1),tol=CGtol,h=dx,number_it=number_of_chi_it)
elseif ((boundary_kind.eq.2).or.(boundary_kind.eq.3)) then
    !In this block:
    !1 - Simple shear flow
    !2 - Periodic pressure-driven flow
    call PCG55_SSOR_x_per_v3(nip,njp,achi(1:nip,1:njp),bchi(1:nip,1:njp),cchi(1:nip,1:njp),dchi(1:nip,1:njp),echi(1:nip,1:njp),fchi(0:nip+1,0:njp+1),chi(0:nip+1,0:njp+1,1),tol=CGtol,h=dx,number_it=number_of_chi_it)
end if

call DATE_AND_TIME(values=tfin)
if(disp_exec_time)  write(unit=*,fmt=1000) '             Elapsed time to solve chi system in seconds:',time_difference(tfin,tini)
1000 format(a60,f8.3)
     
return
end subroutine solve_chi_system