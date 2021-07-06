subroutine flow_update()
implicit none

if     ((boundary_kind.eq.1).or.(boundary_kind.eq.2)) then
    !In this block
    !1 - Lid shear-driven cavity flow
    !2 - Simple shear flow    
    call flow_update_CN()  !CN stands for Cranck-Nicolson
    call get_ucenter()
    call ghosts_update()

elseif(boundary_kind.eq.3) then
    !In this block
    !1 - Pressure-driven flow
    call flow_update_CN_pressure_driven()
    call get_ucenter()
    call ghosts_update()
end if


contains
!==================================================================================================================================
subroutine flow_update_CN()
implicit none
integer :: i,j

!velocity updates
do j=1,nju
    do i=1,niu
        u(i,j) = ustar(i,j)-dt*(chi(i+1,j,1)-chi(i,j,1))/dx
    end do
end do

do j=1,njv
    do i=1,niv
        v(i,j) = vstar(i,j)-dt*(chi(i,j+1,1)-chi(i,j,1))/dy
    end do
end do

!pressure update
do j=1,njp
    do i=1,nip
        press(i,j)=chi(i,j,1)-0.5d0*divustar(i,j)/Re
        divu(i,j) =(u(i,j)-u(i-1,j))/dx+(v(i,j)-v(i,j-1))/dy
    end do
end do

return
end subroutine flow_update_CN
!==================================================================================================================================

!==================================================================================================================================
subroutine flow_update_CN_pressure_driven()
implicit none
integer :: i,j

!Differences in relation to flow_update_CN_pressure_driven() are in u, v and pressure updates.


!velocity updates
do j=1,nju
    do i=1,niu
        u(i,j) = ustar(i,j)-12*dt*(chi(i+1,j,1)-chi(i,j,1))/(Re*dx)
    end do
end do

do j=1,njv
    do i=1,niv
        v(i,j) = vstar(i,j)-12*dt*(chi(i,j+1,1)-chi(i,j,1))/(Re*dy)
    end do
end do

!pressure update
do j=1,njp
    do i=1,nip
        press(i,j)=chi(i,j,1)-divustar(i,j)/24 !Missing Re?
        divu(i,j) =(u(i,j)-u(i-1,j))/dx+(v(i,j)-v(i,j-1))/dy
    end do
end do

return
end subroutine flow_update_CN_pressure_driven
!==================================================================================================================================
end subroutine flow_update
