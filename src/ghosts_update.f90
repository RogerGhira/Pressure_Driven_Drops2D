subroutine ghosts_update()
!ghosts are only at external points
implicit none

if     (boundary_kind.eq.1) then
    !Lid shear-driven cavity flow. Upper surface move with unitary velocity
    call cavity_1st_order()
elseif (boundary_kind.eq.2) then
    !Simple shear flow. Upper and down surfaces move with speed 1 and -1, respectively
    call x_periodic_simple_shear()
elseif (boundary_kind.eq.3) then
    !Pressure driven flow. Upper and down surfaces are at uw=0.0
    call x_periodic_pressure_driven()
end if

return
contains
!==================================================================================================================================
subroutine cavity_2nd_order()
integer :: i,j
real(8) :: uwall,vwall
real(8) :: uwall_prime,vwall_prime

!u ghosts
!tangent components
do i=1,niu
    uwall=0.0d0
    uwall_prime=uwall+dt*(2*(chi(i+1,1,1)-chi(i,1,1))-(chi(i+1,1,0)-chi(i,1,0)))/dx
    ustar(i,0)=2*uwall_prime

    uwall=1.0d0
    uwall_prime=uwall+dt*(2*(chi(i+1,njp,1)-chi(i,njp,1))-(chi(i+1,njp,0)-chi(i,njp,0)))/dx
    ustar(i,nju+1)=2*uwall_prime
end do

!normal components
do j=0,nju+1
    uwall=0.0d0
    ustar(0,j)=uwall

    uwall=0.0d0
    ustar(niu+1,j)=uwall
end do

!v ghosts
!normal components
do i=0,niv+1
    vwall=0.0d0
    v(i,0)=vwall              
    
    vwall=0.0d0
    v(i,njv+1)=vwall   
end do

!tangent components
do j=1,njv    
    vwall=0.0d0
    vwall_prime=vwall+dt*(2*(chi(1,j+1,1)-chi(1,j,1))-(chi(1,j+1,0)-chi(1,j,0)))/dy
    vstar(0,j)=2*vwall
    
    vwall=0.0d0
    vwall_prime=vwall+dt*(2*(chi(nip,j+1,1)-chi(nip,j,1))-(chi(nip,j+1,0)-chi(nip,j,0)))/dy
    vstar(niv+1,j)=2*vwall
end do

!chi ghosts - only lagged edges
do i=0,nip+1
    chi(i,0,:)=0.0d0
    chi(i,njp+1,:)=0.0d0
end do

do j=0,njp+1
    chi(0,j,:)=0.0d0
    chi(nip+1,j,:)=0.0d0
end do

return
end subroutine cavity_2nd_order
!==================================================================================================================================

!==================================================================================================================================
subroutine cavity_1st_order()
implicit none

integer :: i,j

!u ghosts--------------------------------------------------------------------------------------------------------------------------
!lagged boundary surfaces
do i=1,niu
    ustar(i,0) = 0.0d0
    ustar(i,nju+1) =2.0d0
    
    u(i,0)     = 0.0d0-u(i,1)
    u(i,nju+1) = 2.0d0-u(i,nju)
    
    u(i,-1)    = 0.0d0-3*u(i,1)
    u(i,nju+2) = 8.0d0-3*u(i,nju)
    
    u(i,-2)    = 0.0d0-5*u(i,1)
    u(i,nju+3) = 12.0d0-5*u(i,nju)
end do

!exact boudary surfaces (boundary points right on the boundary condition)
do j=0,nju+1 
    ustar(0,j)=0.0d0
    ustar(niu+1,j)=0.0d0
end do
!----------------------------------------------------------------------------------------------------------------------------------

!v ghosts--------------------------------------------------------------------------------------------------------------------------
do j=0,njv+1 !for v component j edges are the lagged ones
    vstar(0,j)=0.0d0
    vstar(niv+1,j)=0.0d0
    
    v(0,j)     = 0.0d0-v(1,j) 
    v(niv+1,j) = 0.0d0-v(niv,j)
    
    v(-1,j)    = 0.0d0-3*v(1,j) 
    v(niv+2,j) = 0.0d0-3*v(niv,j)
    
    v(-2,j)    = 0.0d0-5*v(1,j) 
    v(niv+3,j) = 0.0d0-5*v(niv,j)
end do

!exact boudary surfaces (boundary points right on the boundary condition)
do i=1,niv !for v component i edges are exact (points right on the edge!)
    vstar(i,0)=0.0d0        
    vstar(i,njv+1)=0.0d0
end do
!----------------------------------------------------------------------------------------------------------------------------------


!chi ghosts - only lagged edges
do i=0,nip+1
    chi(i,0,:)=0.0d0
    chi(i,njp+1,:)=0.0d0
end do

do j=0,njp+1
    chi(0,j,:)=0.0d0
    chi(nip+1,j,:)=0.0d0
end do


return
end subroutine cavity_1st_order
!==================================================================================================================================

!==================================================================================================================================
subroutine x_periodic_simple_shear()
implicit none

integer :: i,j
real(8) :: Ut_north,Ut_south,Ut_north_prime,Ut_south_prime

Ut_south= -(ymin+Ly)
Ut_north= +(ymin+Ly)

!u ghosts
do i=1,niu
    
    Ut_south_prime=Ut_south+dt*(2*(chi(i+1,1,1)-chi(i,1,1))-(chi(i+1,1,0)-chi(i,1,0)))/dx
    Ut_north_prime=Ut_north+dt*(2*(chi(i+1,njp,1)-chi(i,njp,1))-(chi(i+1,njp,0)-chi(i,njp,0)))/dx

    ustar(i,0)    = 2*Ut_south_prime      !for u component, i edges are the lagged ones
    ustar(i,nju+1)= 2*Ut_north_prime

    u(i,0)    =2*Ut_south-u(i,1)
    u(i,nju+1)=2*Ut_north-u(i,nju)

    u(i,-1)    = 4*Ut_south-3*u(i,1)
    u(i,nju+2) = 4*Ut_north-3*u(i,nju)
    
    u(i,-2)    = 6*Ut_south-5*u(i,1)
    u(i,nju+3) = 6*Ut_north-5*u(i,nju)
end do

do j=0,nju+1 !periodic BC
    ustar(0,j)=ustar(niu,j)
    ustar(niu+1,j)=ustar(1,j)
    
    u(0,j)=u(niu,j)
    u(niu+1,j)=u(1,j)
    
    u(-1,j)=u(niu-1,j)
    u(niu+2,j)=u(2,j)
    
    u(-2,j)=u(niu-2,j)
    u(niu+3,j)=u(3,j)
end do

!v ghosts
do i=1,niv !for v component i edges are exact (points right on the edge!)
    vstar(i,0)=0.0d0        
    vstar(i,njv+1)=0.0d0

    v(i,0)=0.0d0              
    v(i,njv+1)=0.0d0    

    v(i,-1)=-v(i,1)              
    v(i,njv+2)=-v(i,njv)    
    
    v(i,-2)=-v(i,2)              
    v(i,njv+3)=-v(i,njv-1)    
end do

do j=0,njv+1 !periodic BC
    vstar(0,j)=vstar(niv,j)
    vstar(niv+1,j)=vstar(1,j)
    
    v(0,j)=v(niv,j)  
    v(niv+1,j)=v(1,j)

    v(-1,j)=v(niv-1,j)
    v(niv+2,j)=v(2,j)
    
    v(-2,j)=v(niv-2,j)
    v(niv+3,j)=v(3,j)
end do

!chi ghosts
do i=1,nip
    chi(i,0,:)=0.0d0
    chi(i,njp+1,:)=0.0d0
end do

do j=1,njp !periodic BC
    chi(0,j,:)=chi(nip,j,:)
    chi(nip+1,j,:)=chi(1,j,:)
end do

return
end subroutine x_periodic_simple_shear
!==================================================================================================================================

!==================================================================================================================================
subroutine x_periodic_pressure_driven()
implicit none

integer :: i,j
real(8) :: Uw,Ut_north,Ut_south,Ut_north_prime,Ut_south_prime

Uw = 0.0d0

!u ghosts
do i=1,niu
    Ut_south = -Uw
    Ut_north = +Uw
    
    Ut_south_prime=Ut_south+dt*(2*(chi(i+1,1,1)-chi(i,1,1))-(chi(i+1,1,0)-chi(i,1,0)))/dx
    Ut_north_prime=Ut_north+dt*(2*(chi(i+1,njp,1)-chi(i,njp,1))-(chi(i+1,njp,0)-chi(i,njp,0)))/dx

    ustar(i,0)    = 2*Ut_south_prime      !for u component, i edges are the lagged ones
    ustar(i,nju+1)= 2*Ut_north_prime

    u(i,0)    =2*Ut_south-u(i,1)
    u(i,nju+1)=2*Ut_north-u(i,nju)

    u(i,-1)    = 4*Ut_south-3*u(i,1)
    u(i,nju+2) = 4*Ut_north-3*u(i,nju)
    
    u(i,-2)    = 6*Ut_south-5*u(i,1)
    u(i,nju+3) = 6*Ut_north-5*u(i,nju)
end do

do j=0,nju+1 !periodic BC
    ustar(0,j)=ustar(niu,j)
    ustar(niu+1,j)=ustar(1,j)
    
    u(0,j)=u(niu,j)
    u(niu+1,j)=u(1,j)
    
    u(-1,j)=u(niu-1,j)
    u(niu+2,j)=u(2,j)
    
    u(-2,j)=u(niu-2,j)
    u(niu+3,j)=u(3,j)
end do

!v ghosts
do i=1,niv !for v component i edges are exact (points right on the edge!)
    vstar(i,0)=0.0d0        
    vstar(i,njv+1)=0.0d0

    v(i,0)=0.0d0              
    v(i,njv+1)=0.0d0    

    v(i,-1)=-v(i,1)              
    v(i,njv+2)=-v(i,njv)    
    
    v(i,-2)=-v(i,2)              
    v(i,njv+3)=-v(i,njv-1)    
end do

do j=0,njv+1 !periodic BC
    vstar(0,j)=vstar(niv,j)
    vstar(niv+1,j)=vstar(1,j)
    
    v(0,j)=v(niv,j)  
    v(niv+1,j)=v(1,j)

    v(-1,j)=v(niv-1,j)
    v(niv+2,j)=v(2,j)
    
    v(-2,j)=v(niv-2,j)
    v(niv+3,j)=v(3,j)
end do

!chi ghosts
do i=1,nip
    chi(i,0,:)=0.0d0
    chi(i,njp+1,:)=0.0d0
end do

do j=1,njp !periodic BC
    chi(0,j,:)=chi(nip,j,:)
    chi(nip+1,j,:)=chi(1,j,:)
end do

return
end subroutine x_periodic_pressure_driven
!==================================================================================================================================

end subroutine ghosts_update