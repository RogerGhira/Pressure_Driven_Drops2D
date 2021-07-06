subroutine BC()
implicit none

if     (boundary_kind.eq.1) then
    call dirichlet()
elseif ((boundary_kind.eq.2).or.(boundary_kind.eq.3)) then
    call x_periodic()
end if

return
contains
!==================================================================================================================================
subroutine dirichlet()
integer :: i,j

!"A" matrix changings - only at internal points
!"u" component
do i=1,niu
    cu(i,1)=cu(i,1)-bu(i,1)    
    cu(i,nju)=cu(i,nju)-eu(i,nju)
end do

!"v" component
do j=1,njv
    cv(1,j)=cv(1,j)-av(1,j)    
    cv(niv,j)=cv(niv,j)-dv(niv,j)
end do

!chi - always homogeneus Neumann
do i=1,nip
    cchi(i,1)=cchi(i,1)+bchi(i,1)        !dchi/dy(y=0) = 0
    cchi(i,njp)=cchi(i,njp)+echi(i,njp)  !dchi/dy(y=Ly) = 0
end do

do j=1,njp
    cchi(1,j)=cchi(1,j)+achi(1,j)        !dchi/dy(x=0) = 0
    cchi(nip,j)=cchi(nip,j)+dchi(nip,j)  !dchi/dy(x=Lx) = 0
end do

return
end subroutine dirichlet
!==================================================================================================================================
!==================================================================================================================================
subroutine x_periodic()
integer :: i,j

!"A" matrix changings - only at internal points
!"u" component
do i=1,niu
    cu(i,1)=cu(i,1)-bu(i,1)    
    cu(i,nju)=cu(i,nju)-eu(i,nju)
end do

!pressure - always homogeneus Neumann
do i=1,nip
    cchi(i,1)=cchi(i,1)+bchi(i,1)        !dchi/dy(y=0) = 0
    cchi(i,njp)=cchi(i,njp)+echi(i,njp)  !dchi/dy(y=Ly) = 0
end do

return
end subroutine x_periodic
!==================================================================================================================================

end subroutine BC