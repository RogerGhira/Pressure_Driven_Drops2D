!This subroutine runs every case so far
!1 - Lid shear-driven flow;
!2 - Simple shear flow;
!3 - Presure-driven flow;
subroutine fill_f_for_chi()
implicit none

if     ((boundary_kind.eq.1).or.(boundary_kind.eq.2)) then
    !In this block
    !Lid shear-driven cavity flow (monophasic flow)
    !Simple shear flow (monophasic flow)
    call fill_f_for_chi_mono_no_body_forces()
elseif (boundary_kind.eq.3) then
    !In this block
    !Pressure-driven flow
    call fill_f_for_chi_pressure_driven()
end if


contains
!==================================================================================================================================
subroutine fill_f_for_chi_mono_no_body_forces()
implicit none

integer :: i,j

do j=1,njp
    do i=1,nip
        divustar(i,j)=(ustar(i,j)-ustar(i-1,j))/dx+(vstar(i,j)-vstar(i,j-1))/dy
        fchi(i,j)=divustar(i,j)/dt
    end do
end do

return
end subroutine fill_f_for_chi_mono_no_body_forces
!==================================================================================================================================

!==================================================================================================================================
subroutine fill_f_for_chi_pressure_driven()
implicit none

integer :: i,j

do j=1,njp
    do i=1,nip
        divustar(i,j)=(ustar(i,j)-ustar(i-1,j))/dx+(vstar(i,j)-vstar(i,j-1))/dy
        fchi(i,j)=Re*divustar(i,j)/(12*dt)
    end do
end do

return
end subroutine fill_f_for_chi_pressure_driven
!==================================================================================================================================
end subroutine fill_f_for_chi
