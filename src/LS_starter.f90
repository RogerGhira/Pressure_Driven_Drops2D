!==================================================================================================================================
subroutine LS_starter()
implicit none
integer :: i,j,k


!Memory allocation
allocate(phi(-2:nip+3,-2:njp+3,0:1))  !lots of ghosts to allow high order eno for grad(phi)
allocate(ugradphi(-2:nip+3,-2:njp+3,0:1))
allocate(phi_x(-2:nip+3,-2:njp+3))
allocate(phi_y(-2:nip+3,-2:njp+3))
allocate(sgnphi(-2:nip+3,-2:njp+3))
allocate(normalx(-2:nip+3,-2:njp+3))
allocate(normaly(-2:nip+3,-2:njp+3))
allocate(kappa(-2:nip+3,-2:njp+3))

!tubing structures
allocate(mask(-2:nip+3,-2:njp+3))
allocate(index_i(1:(nip+6)*(njp+6)))
allocate(index_j(1:(nip+6)*(njp+6)))

!variables safety initialization
phi=2.0d0*max(maxval(x)-minval(x),maxval(y)-minval(y)); ugradphi=0.0d0; phi_x=0.0d0; phi_y=0.0d0; sgnphi=0.0d0
normalx=0.0d0; normaly=0.0d0; kappa=0.0d0
n_tube=0; mask=0; index_i=0; index_j=0; alpha=3.0d0*max(dx,dy); beta=6.0d0*max(dx,dy); gamma=9.0d0*max(dx,dy)

!initial level set creation
call set_circunference(2.5d0,0.5d0,3.0d-1,100)

phi(:,:,0)=phi(:,:,1) !all operations are made on phi(:,:,1). Here we are just setting up phi(:,:,0)

!tubing structure and geometric functions initialization!
call make_LS_tubes()

!initialization of "gradient" quantities
do k=1,n_tube
    i=index_i(k)
    j=index_j(k)
    phi_x(i,j) = (phi(i+1,j,1)-phi(i-1,j,1))/(2.0d0*dx)
    phi_y(i,j) = (phi(i,j+1,1)-phi(i,j-1,1))/(2.0d0*dy)
    sgnphi(i,j) = phi(i,j,1)/sqrt(phi(i,j,1)*phi(i,j,1)+(phi_x(i,j)**2+phi_y(i,j)**2)*(0.5d0*(dx+dy))**2)
end do

call get_geometry()

call print_results_to_tecplot_LS('initial_condition.plt')

return
end subroutine LS_starter
!==================================================================================================================================
