subroutine ugradu_update()
implicit none

!call ugradu_update_central()
call ugradu_update_ENO2()

return
contains
!==================================================================================================================================
subroutine ugradu_update_central()
implicit none

integer :: i,j
real :: ubar,vbar,dudx,dudy,dvdx,dvdy

do j=1,nju
    do i=1,niu
        dudx=0.5d0*(u(i+1,j)-u(i-1,j))/dx
        dudy=0.5d0*(u(i,j+1)-u(i,j-1))/dy
        vbar=0.25d0*(v(i,j)+v(i+1,j)+v(i+1,j-1)+v(i,j-1))
        ugradu(i,j,0) = ugradu(i,j,1)
        ugradu(i,j,1) = u(i,j)*dudx+vbar*dudy
    end do
end do

do j=1,njv
    do i=1,niv
        dvdx=0.5d0*(v(i+1,j)-v(i-1,j))/dx
        dvdy=0.5d0*(v(i,j+1)-v(i,j-1))/dy
        ubar=0.25d0*(u(i,j)+u(i,j+1)+u(i-1,j+1)+u(i-1,j))
        ugradv(i,j,0) = ugradv(i,j,1)
        ugradv(i,j,1) = ubar*dvdx+v(i,j)*dvdy
    end do
end do

return
end subroutine ugradu_update_central
!==================================================================================================================================

!==================================================================================================================================
subroutine ugradu_update_ENO2()
implicit none

integer :: i,j,m
real :: ubar,vbar
real :: dudx, dudy, dudx_p, dudx_m, dudy_p, dudy_m
real :: dvdx, dvdy, dvdx_p, dvdx_m, dvdy_p, dvdy_m
real :: Q1prime, Q2prime, c_eno


!building the divided differences for u 
D1x(-1:niu+1, 1:nju)=(u(0:niu+2,1:nju)-u(-1:niu+1,1:nju))/dx   !NOTE: D1x(i,j) means D1x(i+1/2,j+1/2); The same for D1y. We are rounding to the floor
D1y(1:niu, -1:nju+1)=(u(1:niu,0:nju+2)-u(1:niu,-1:nju+1))/dy
D2x(0:niu+1, 1:nju)=0.50d0*(D1x(0:niu+1,1:nju)-D1x(-1:niu,1:nju))/dx
D2y(1:niu, 0:nju+1)=0.50d0*(D1y(1:niu,0:nju+1)-D1y(1:niu,-1:nju))/dy

do j=1,nju
    do i=1,niu
        
        !--------------------------------------------------------------------------------------------------------------------------
        if(u(i,j).ge.0.0d0) then
            m=i-1
        else
            m=i
        end if
        
        Q1prime=D1x(m,j)
        if(abs(D2x(m,j)).le.abs(D2x(m+1,j))) then
            c_eno=D2x(m,j)
        else
            c_eno=D2x(m+1,j)
        end if
        Q2prime=c_eno*(2.0d0*(i-m)-1.0d0)*dx
        dudx=Q1prime+Q2prime
        !--------------------------------------------------------------------------------------------------------------------------
        
        !--------------------------------------------------------------------------------------------------------------------------
        vbar=0.25d0*(v(i,j)+v(i+1,j)+v(i+1,j-1)+v(i,j-1)) ! v-component at u(i,j) position: mean of the four neighbours
        if(vbar.ge.0.0d0) then
            m=j-1
        else
            m=j
        end if

        Q1prime=D1y(i,m)
        if(abs(D2y(i,m)).le.abs(D2y(i,m+1))) then
            c_eno=D2y(i,m)
        else
            c_eno=D2y(i,m+1)
        end if
        Q2prime=c_eno*(2.0d0*(j-m)-1.0d0)*dx
        dudy=Q1prime+Q2prime
        !--------------------------------------------------------------------------------------------------------------------------
        ugradu(i,j,0) = ugradu(i,j,1)  !ugradu update
        ugradu(i,j,1) = u(i,j)*dudx+vbar*dudy
    end do
end do

!building the divided differences for v 
D1x(-1:niv+1, 1:njv)=(v(0:niv+2,1:njv)-v(-1:niv+1,1:njv))/dx    !NOTE: D1x(i,j) means D1x(i+1/2,j+1/2); The same for D1y. We are rounding to the floor
D1y(1:niv, -1:njv+1)=(v(1:niv,0:njv+2)-v(1:niv,-1:njv+1))/dy
D2x(0:niv+1, 1:njv)=0.50d0*(D1x(0:niv+1,1:njv)-D1x(-1:niv,1:njv))/dx
D2y(1:niv, 0:njv+1)=0.50d0*(D1y(1:niv,0:njv+1)-D1y(1:niv,-1:njv))/dy

do j=1,njv
    do i=1,niv
        !--------------------------------------------------------------------------------------------------------------------------
        ubar=0.25d0*(u(i,j)+u(i,j+1)+u(i-1,j+1)+u(i-1,j))  ! u-component at v(i,j) position: mean of the four neighbours
        if(ubar.ge.0.0d0) then
            m=i-1
        else
            m=i
        end if
        Q1prime=D1x(m,j)
        if(abs(D2x(m,j)).le.abs(D2x(m+1,j))) then
            c_eno=D2x(m,j)
        else
            c_eno=D2x(m+1,j)
        end if
        Q2prime=c_eno*(2.0d0*(i-m)-1.0d0)*dx
        dvdx=Q1prime+Q2prime
        !--------------------------------------------------------------------------------------------------------------------------

        !computation of dudy_p and dudy_m------------------------------------------------------------------------------------------
        if(v(i,j).ge.0.0d0) then
            m=j-1
        else
            m=j
        end if
        m=j-1
        Q1prime=D1y(i,m)
        if(abs(D2y(i,m)).le.abs(D2y(i,m+1))) then
            c_eno=D2y(i,m)
        else
            c_eno=D2y(i,m+1)
        end if
        Q2prime=c_eno*(2.0d0*(j-m)-1.0d0)*dy
        dvdy=Q1prime+Q2prime
        !--------------------------------------------------------------------------------------------------------------------------        

        !--------------------------------------------------------------------------------------------------------------------------        
        ugradv(i,j,0) = ugradv(i,j,1)  !ugradv update
        ugradv(i,j,1) = ubar*dvdx + v(i,j)*dvdy
    end do
end do

return
end subroutine ugradu_update_ENO2
!==================================================================================================================================

end subroutine ugradu_update

!**********************************************************************************************************************************
!**********************************************************************************************************************************
subroutine ugradu_update_for_primer()
implicit none

!call ugradu_update_central_for_primer()
call ugradu_update_ENO2_for_primer()

return
contains
!==================================================================================================================================
subroutine ugradu_update_central_for_primer() !basically the same of ugradu_update_central() but with different ugradu update
implicit none

integer :: i,j
real :: ubar,vbar,dudx,dudy,dvdx,dvdy

do j=1,nju
    do i=1,niu
        dudx=0.5d0*(u(i+1,j)-u(i-1,j))/dx
        dudy=0.5d0*(u(i,j+1)-u(i,j-1))/dy
        vbar=0.25d0*(v(i,j)+v(i+1,j)+v(i+1,j-1)+v(i,j-1))
        ugradu(i,j,0) = u(i,j)*dudx+vbar*dudy  !this is the sole diference between this subroutine and the regular one
    end do
end do

do j=1,njv
    do i=1,niv
        dvdx=0.5d0*(v(i+1,j)-v(i-1,j))/dx
        dvdy=0.5d0*(v(i,j+1)-v(i,j-1))/dy
        ubar=0.25d0*(u(i,j)+u(i,j+1)+u(i-1,j+1)+u(i-1,j))
        ugradv(i,j,0) = ubar*dvdx+v(i,j)*dvdy  !this is the sole diference between this subroutine and the regular one
    end do
end do

return
end subroutine ugradu_update_central_for_primer
!==================================================================================================================================

!==================================================================================================================================
subroutine ugradu_update_ENO2_for_primer() !basically the same of ugradu_update_ENO2() but with different ugradu update 
implicit none

integer :: i,j,m
real :: ubar,vbar
real :: dudx, dudy, dudx_p, dudx_m, dudy_p, dudy_m
real :: dvdx, dvdy, dvdx_p, dvdx_m, dvdy_p, dvdy_m
real :: Q1prime, Q2prime, c_eno


!building the divided differences for u 
D1x(-1:niu+1, 1:nju)=(u(0:niu+2,1:nju)-u(-1:niu+1,1:nju))/dx   !NOTE: D1x(i,j) means D1x(i+1/2,j+1/2); The same for D1y. We are rounding to the floor
D1y(1:niu, -1:nju+1)=(u(1:niu,0:nju+2)-u(1:niu,-1:nju+1))/dy
D2x(0:niu+1, 1:nju)=0.50d0*(D1x(0:niu+1,1:nju)-D1x(-1:niu,1:nju))/dx
D2y(1:niu, 0:nju+1)=0.50d0*(D1y(1:niu,0:nju+1)-D1y(1:niu,-1:nju))/dy

do j=1,nju
    do i=1,niu
        
        !--------------------------------------------------------------------------------------------------------------------------
        if(u(i,j).ge.0.0d0) then
            m=i-1
        else
            m=i
        end if
        
        Q1prime=D1x(m,j)
        if(abs(D2x(m,j)).le.abs(D2x(m+1,j))) then
            c_eno=D2x(m,j)
        else
            c_eno=D2x(m+1,j)
        end if
        Q2prime=c_eno*(2.0d0*(i-m)-1.0d0)*dx
        dudx=Q1prime+Q2prime
        !--------------------------------------------------------------------------------------------------------------------------
        
        !--------------------------------------------------------------------------------------------------------------------------
        vbar=0.25d0*(v(i,j)+v(i+1,j)+v(i+1,j-1)+v(i,j-1)) ! v-component at u(i,j) position: mean of the four neighbours
        if(vbar.ge.0.0d0) then
            m=j-1
        else
            m=j
        end if

        Q1prime=D1y(i,m)
        if(abs(D2y(i,m)).le.abs(D2y(i,m+1))) then
            c_eno=D2y(i,m)
        else
            c_eno=D2y(i,m+1)
        end if
        Q2prime=c_eno*(2.0d0*(j-m)-1.0d0)*dx
        dudy=Q1prime+Q2prime
        !--------------------------------------------------------------------------------------------------------------------------
        ugradu(i,j,0) = u(i,j)*dudx+vbar*dudy !this is the sole diference between this subroutine and the regular one
    end do
end do

!building the divided differences for v 
D1x(-1:niv+1, 1:njv)=(v(0:niv+2,1:njv)-v(-1:niv+1,1:njv))/dx    !NOTE: D1x(i,j) means D1x(i+1/2,j+1/2); The same for D1y. We are rounding to the floor
D1y(1:niv, -1:njv+1)=(v(1:niv,0:njv+2)-v(1:niv,-1:njv+1))/dy
D2x(0:niv+1, 1:njv)=0.50d0*(D1x(0:niv+1,1:njv)-D1x(-1:niv,1:njv))/dx
D2y(1:niv, 0:njv+1)=0.50d0*(D1y(1:niv,0:njv+1)-D1y(1:niv,-1:njv))/dy

do j=1,njv
    do i=1,niv
        !--------------------------------------------------------------------------------------------------------------------------
        ubar=0.25d0*(u(i,j)+u(i,j+1)+u(i-1,j+1)+u(i-1,j))  ! u-component at v(i,j) position: mean of the four neighbours
        if(ubar.ge.0.0d0) then
            m=i-1
        else
            m=i
        end if
        Q1prime=D1x(m,j)
        if(abs(D2x(m,j)).le.abs(D2x(m+1,j))) then
            c_eno=D2x(m,j)
        else
            c_eno=D2x(m+1,j)
        end if
        Q2prime=c_eno*(2.0d0*(i-m)-1.0d0)*dx
        dvdx=Q1prime+Q2prime
        !--------------------------------------------------------------------------------------------------------------------------

        !computation of dudy_p and dudy_m------------------------------------------------------------------------------------------
        if(v(i,j).ge.0.0d0) then
            m=j-1
        else
            m=j
        end if
        m=j-1
        Q1prime=D1y(i,m)
        if(abs(D2y(i,m)).le.abs(D2y(i,m+1))) then
            c_eno=D2y(i,m)
        else
            c_eno=D2y(i,m+1)
        end if
        Q2prime=c_eno*(2.0d0*(j-m)-1.0d0)*dy
        dvdy=Q1prime+Q2prime
        !--------------------------------------------------------------------------------------------------------------------------        

        !--------------------------------------------------------------------------------------------------------------------------        
        ugradv(i,j,0) = ubar*dvdx + v(i,j)*dvdy !this is the sole diference between this subroutine and the regular one
    end do
end do

return
end subroutine ugradu_update_ENO2_for_primer
!==================================================================================================================================

end subroutine ugradu_update_for_primer