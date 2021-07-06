!Main NSProj_2D boundary condition subroutine.
subroutine fill_f_for_uv()
implicit none

if     ((boundary_kind.eq.1).or.(boundary_kind.eq.2)) then
    !In this block
    !Lid shear-driven cavity flow (monophasic flow)
    !Simple shear flow (monophasic flow)
    call fill_f_for_uv_mono_no_body_forces()
elseif (boundary_kind.eq.3) then
    !In this block
    !Pressure-driven flow
    call fill_f_for_uv_pressure_driven()
end if


contains
!==================================================================================================================================
subroutine fill_f_for_uv_mono_no_body_forces()
implicit none

integer :: i,j
real :: d2udx2,d2udy2,d2vdx2,d2vdy2

!u component
do j=1,nju
    do i=1,niu
        d2udx2=(u(i-1,j)-2*u(i,j)+u(i+1,j))/dx**2
        d2udy2=(u(i,j-1)-2*u(i,j)+u(i,j+1))/dy**2
        fu(i,j) = u(i,j)+dt*(-1.5d0*ugradu(i,j,1)+0.5d0*ugradu(i,j,0)+0.5d0*(d2udx2+d2udy2)/Re)
    end do
end do

!v component
do j=1,njv
    do i=1,niv
        d2vdx2=(v(i-1,j)-2*v(i,j)+v(i+1,j))/dx**2
        d2vdy2=(v(i,j-1)-2*v(i,j)+v(i,j+1))/dy**2
        fv(i,j) = v(i,j)+dt*(-1.5d0*ugradv(i,j,1)+0.5d0*ugradv(i,j,0)+0.5d0*(d2vdx2+d2vdy2)/Re)
    end do
end do

return
end subroutine fill_f_for_uv_mono_no_body_forces
!==================================================================================================================================

!==================================================================================================================================
subroutine fill_f_for_uv_pressure_driven()
implicit none

integer :: i,j
real :: d2udx2,d2udy2,d2vdx2,d2vdy2

!The only change in relation to fill_f_for_uv_mono_no_body_forces() is the 12/Re term in u-equation
!u component
do j=1,nju
    do i=1,niu
        d2udx2=(u(i-1,j)-2*u(i,j)+u(i+1,j))/dx**2
        d2udy2=(u(i,j-1)-2*u(i,j)+u(i,j+1))/dy**2
        fu(i,j) = u(i,j)+dt*(-1.5d0*ugradu(i,j,1)+0.5d0*ugradu(i,j,0)+0.5d0*(d2udx2+d2udy2)/Re + 12.0d0/Re)
    end do
end do

!v component
do j=1,njv
    do i=1,niv
        d2vdx2=(v(i-1,j)-2*v(i,j)+v(i+1,j))/dx**2
        d2vdy2=(v(i,j-1)-2*v(i,j)+v(i,j+1))/dy**2
        fv(i,j) = v(i,j)+dt*(-1.5d0*ugradv(i,j,1)+0.5d0*ugradv(i,j,0)+0.5d0*(d2vdx2+d2vdy2)/Re)
    end do
end do

return
end subroutine fill_f_for_uv_pressure_driven
!==================================================================================================================================
end subroutine fill_f_for_uv
