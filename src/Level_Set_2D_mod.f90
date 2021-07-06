module Level_Set_2D_mod
use Navier_Stokes_2D_mod !Parent module where the flow came from
implicit none
save  !guarantees that all data values declared in this module will be preserved between references in different proc.

!Level set structure
real(8), allocatable :: phi(:,:,:)      !Level set signaled distance function --> phi(-2:nip+3,-2:njp+3,0:1)
real(8), allocatable :: ugradphi(:,:,:) !u*Grad(phi) --> ugradphi(-2:nip+3,-2:njp+3,0:1)
real(8), allocatable :: phi_x(:,:)      !Grad(phi)*ex. Stored in the cell center (collocated scheme)  --> phi_x(-2:nip+3,-2:njp+3)
real(8), allocatable :: phi_y(:,:)      !Grad(phi)*ey. Stored in the cell center (collocated scheme)  --> phi_y(-2:nip+3,-2:njp+3)
real(8), allocatable :: sgnphi(:,:)     !Signal function. Stored in the cell center (collocated scheme)  --> sgnphi(-2:nip+3,-2:njp+3)
real(8), allocatable :: normalx(:,:)    !Normal vector "x" component. Stored in the cell center (collocated scheme)  --> normalx(-2:nip+3,-2:njp+3)
real(8), allocatable :: normaly(:,:)    !Normal vector "y" component. Stored in the cell center (collocated scheme)  --> normaly(-2:nip+3,-2:njp+3)
real(8), allocatable :: kappa(:,:)      !Mean curvature. Stored in the cell center (collocated scheme)  --> kappa(-2:nip+3,-2:njp+3)

!tubing structure
integer :: n_tube
real(8) :: alpha, beta, gamma       !three tubes
integer, allocatable :: mask(:,:)   !mask for tubing ((narrow band) identification
integer, allocatable :: index_i(:)  !"i" index storage array
integer, allocatable :: index_j(:)  !"j" index storage array

!math constants
real(8), parameter :: pi=3.1415926535897932384626433832795d0

!WENO variables

contains
!==================================================================================================================================
include 'LS_starter.f90'
include 'set_circunference.f90'
include 'make_LS_tubes.f90'
include 'print_results_to_tecplot_LS.f90'
include 'get_geometrical_quantities.f90'
include 'evolve_LS.f90'
!==================================================================================================================================
end module Level_Set_2D_mod
