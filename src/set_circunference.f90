!==================================================================================================================================
subroutine set_circunference(x0,y0,R,nw)
implicit none

real, intent(in) :: x0,y0,R
integer, intent(in) :: nw

integer :: i,j,k
real, allocatable :: theta(:)
real :: dtheta, omega, domega, dlocal, dmin, NRtol, rlocal
real :: p(2),z(2),t(2),v(2),f1,f2,sig

!numerical parameters
domega=1.0d-4
NRtol=1.0d-10

!theta initialization
allocate(theta(1:nw))
dtheta=2*pi/(nw-1)
do k=1,nw
    theta(k)=(k-1)*dtheta
end do

!main loop
do concurrent (i=0:nip+1,j=0:njp+1)
    
    !getting of "p" and translation
    p(1)=0.5d0*(x(i-1,j)+x(i,j))-x0
    p(2)=0.5d0*(y(i,j-1)+y(i,j))-y0
    
    !initial guess for Newton Raphson
    dmin=1.0d10 !initialization of the distance
    do k=1,nw
        z(1)=xpar(R,theta(k))
        z(2)=ypar(R,theta(k))
        dlocal=sqrt((p(1)-z(1))**2+(p(2)-z(2))**2)
        if(dlocal.lt.dmin) then
            dmin=dlocal
            omega=theta(k) !inital guess
        end if
    end do
    
    !Newton Raphson to refine solution
    do
        f1=Fcir(R,omega,p(1),p(2))
        if(abs(f1).lt.NRtol) exit
        f2=Fcir(R,omega+domega,p(1),p(2))
        omega=omega-domega*f1/(f2-f1)
    end do
    
    !deciding for the signal
    z(1)=xpar(R,omega)
    z(2)=ypar(R,omega)
    
    t(1)=xpar(R,omega+domega)-z(1)
    t(2)=ypar(R,omega+domega)-z(2)
    t=t/sqrt(t(1)**2+t(2)**2)
    
    v=(p-z)
    rlocal=sqrt(v(1)**2+v(2)**2)
    if(rlocal.eq.0.0d0) then
        sig=0.0d0
    else
        if(v(1)*t(2)-v(2)*t(1).gt.0.0d0) then 
            sig=1.0d0
        else
            sig=-1.0d0
        end if
    end if
    
    !computing phi
    if(abs(rlocal).lt.abs(phi(i,j,1))) phi(i,j,1) = sig*rlocal
    
end do

return
contains
!----------------------------------------------------------------------------------------------------------------------------------
pure function xpar(R,theta)
implicit none
real, intent(in) :: R,theta
real :: xpar

xpar=R*cos(theta)

end function xpar
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
pure function ypar(R,theta)
implicit none
real, intent(in) :: R,theta
real :: ypar

ypar=R*sin(theta)

end function ypar
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
pure function Fcir(R,theta,p,q)
implicit none
real, intent(in) :: R,theta,p,q
real :: Fcir

Fcir=-(xpar(R,theta)-p)*R*sin(theta)+(ypar(R,theta)-q)*R*cos(theta)

end function Fcir
!----------------------------------------------------------------------------------------------------------------------------------


end subroutine set_circunference
!==================================================================================================================================
