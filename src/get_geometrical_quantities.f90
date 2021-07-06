subroutine get_geometry()
implicit none

call get_geometrical_quantities_tube()  !local geometrical computations are leading to errors in extension procedure

return
contains    
!==================================================================================================================================
subroutine get_geometrical_quantities_tube()
!Computes normal vector by n=grad(psi)/|grad(psi)|, using central differences for psi_x and psi_y
!Computes mean curvature by kappa = div(n), using central differences for nx_x and ny_y
!This is the same procedure done by Croce, Griebel and Schweitzer, International Journal for Numerical Methods in Fluids, 2010
use tools, only: get_machine_eps
implicit none

integer :: i,j,k
real :: grad_phi_norm,phi_x_local,phi_y_local,eps


eps=get_machine_eps()

!interior normal vector - normal vector must be computed overall external tube to allow curvature computation into the interior tube.
do k=1,n_tube
    i=index_i(k)
    j=index_j(k)
    phi_x_local = (phi(i+1,j,1)-phi(i-1,j,1))/(2*dx)
    phi_y_local = (phi(i,j+1,1)-phi(i,j-1,1))/(2*dy)
    grad_phi_norm = sqrt(phi_x_local**2+phi_y_local**2)
    !normalx(i,j) = (psi_x_local)/(grad_psi_norm)
    !normaly(i,j) = (psi_y_local)/(grad_psi_norm)
    normalx(i,j) = (phi_x_local+eps)/(grad_phi_norm+eps)
    normaly(i,j) = (phi_y_local+eps)/(grad_phi_norm+eps)
end do

!interior curvature
do k=1,n_tube
    i=index_i(k)
    j=index_j(k)
    if(mask(i,j).le.2) then
        kappa(i,j) = (normalx(i+1,j)-normalx(i-1,j))/(2*dx)+(normaly(i,j+1)-normaly(i,j-1))/(2*dy)
    else
        kappa(i,j) = 0.0d0
    end if
end do

return
end subroutine get_geometrical_quantities_tube
!==================================================================================================================================

end subroutine get_geometry
