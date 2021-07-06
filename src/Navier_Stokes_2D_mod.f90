module Navier_Stokes_2D_mod
use Tools  !The only thing we really need from here is the variable "filn_len"
implicit none
save  !guarantees that all data values declared in this module will be preserved between references in different proc.

!Computational domain and grid definitions - MAC: Marker and Cell scheme - Staggered grid to avoid odd-even decoupling 
!                                                                          between the pressure and velocity.
!Rectangular domains only: (Lx,Ly)
!Cartesian grid only: (nx,ny)

!Grid variables
real(8) :: xmin  !"x" starting point for grid construction
real(8) :: ymin  !"y" starting point for grid construction
real(8) :: Lx  !length in "x" direction
real(8) :: Ly  !length in "y" direction
real(8) :: dx  !"x" direction grid spacing - "x" norm
real(8) :: dy  !"y" direction grid spacing - "y" norm

integer :: nx  !number of pressure nodes in "x" direction
integer :: ny  !number of pressure nodes in "y" direction
integer :: niu,nju  !limits of unkwon for "u" component
integer :: niv,njv  !limits of unkwon for "u" component
integer :: nip,njp  !limits of unkwon for pressure

real(8), allocatable :: x(:,:) !"x" coordinates of grid edges: x(0:nx+1,0:ny+1). 
real(8), allocatable :: y(:,:) !"y" coordinates of grid edges: y(0:nx+1,0:ny+1).

!Flow variables - including projection method auxiliars
real(8), allocatable :: u(:,:) !"x" velocity component: u(-2:niu+3,-2:nju+3). NOTE: Wider vectors are necessary to allow straight ENO2 calculations with periodic condictions
real(8), allocatable :: v(:,:) !"y" velocity component: v(-2:niv+3,-2:njv+3). NOTE: Wider vectors are necessary to allow straight ENO2 calculations with periodic condictions
real(8), allocatable :: ugradu(:,:,:) !convective term in x direction ugradu(0:niu+1,0:nuv+1,0:1) --> 0 means "n-1"; 1 means "n";
real(8), allocatable :: ugradv(:,:,:) !convective term in x direction ugradu(0:niu+1,0:nuv+1,0:1) --> 0 means "n-1"; 1 means "n";
real(8), allocatable :: press(:,:) !pressure: p(0:nx,0:ny).
real(8), allocatable :: divustar(:,:)  !divergence of the tentative velocity
real(8), allocatable :: divu(:,:)  !divergence of the velocity vector - must be zero

real(8), allocatable :: uinitial(:,:) !initial velocity field - "x" direction.
real(8), allocatable :: vinitial(:,:) !initial velocity field - "y" direction.
real(8), allocatable :: pressinitial(:,:) !initial pressure field.
real(8), allocatable :: ucenter(:,:)  !u velocity in the cell center (interpolated).
real(8), allocatable :: vcenter(:,:)  !u velocity in the cell center (interpolated).
real(8), allocatable :: ustar(:,:) !"x" tentative velocity component ustar(0:niu+1,0:nju+1).
real(8), allocatable :: vstar(:,:) !"y" tentative velocity component vstar(0:niv+1,0:njv+1).
real(8), allocatable :: chi(:,:,:) !false pressure used in the projection method
real(8), allocatable :: swap(:,:) !auxiliar for variables swapping

!ENO and WENO auxiliary variables
real(8), allocatable :: D1x(:,:) !First divided differences. NOTE: D1x(i,j) means D1x(i+1/2,j+1/2); The same for D1y
real(8), allocatable :: D1y(:,:) !First divided differences.
real(8), allocatable :: D2x(:,:) !First divided differences. 
real(8), allocatable :: D2y(:,:) !First divided differences. 

real(8), allocatable :: au(:,:),bu(:,:),cu(:,:),du(:,:),eu(:,:),fu(:,:)  !linear system components: lots of memory for almost notting!
real(8), allocatable :: av(:,:),bv(:,:),cv(:,:),dv(:,:),ev(:,:),fv(:,:)  !linear system components: lots of memory for almost notting!
real(8), allocatable :: achi(:,:),bchi(:,:),cchi(:,:),dchi(:,:),echi(:,:),fchi(:,:)  !linear system components: lots of memory for almost notting!

!Flow and numerical data
real(8) :: Re      !Reynolds number
real(8) :: dt      !constant time step
real(8) :: time    !physical simulation time
real(8) :: Tmax    !Simulation final time

!Control data
logical :: bck_read_flag=.false.  !flag to read or not the backup
logical :: single_bck_file=.true. !flag to only save the last backup file. If TRUE only the last saved iteration will be kept.
logical :: disp_exec_time=.false. !flag to only save the last backup file. If TRUE only the last saved iteration will be kept.

character(len=filn_len) :: bck_file !backup file

integer :: first_it   !number of the first iteration
integer :: nit        !total number of iterations
integer :: bck_sv_frq !backup saving frequency
integer :: plt_sv_frq !ploting saving frequency

!Boundary condition control data
integer :: boundary_kind !kind of boundary condition. See codes in "Drops2d_2020.dat"
integer :: initial_kind  !flag to use "initial_flow_field" (initial_kind=1) or set u,v, and press to 0 (initial_kind=0)

contains
!==================================================================================================================================
include 'NS_starter.f90'
include 'primer_for_first_iteration.f90'
include 'fill_abc.f90'
include 'ghosts_update.f90' 
include 'BC.f90'
include 'fill_f_for_uv.f90'
include 'fill_f_for_chi.f90'
include 'flow_update.f90'
include 'ugradu_update.f90'
include 'print_results_to_tecplot.f90'
include 'print_backup.f90'
include 'solve_ustar_system.f90'
include 'solve_chi_system.f90'
include 'set_initial_flow_field.f90'
include 'get_ucenter.f90'
!==================================================================================================================================
end module Navier_Stokes_2D_mod
