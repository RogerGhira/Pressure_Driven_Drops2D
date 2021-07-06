module LSys
    
private
public :: CG5,CG5_x_per,PCG5_diag,PCG5_diag_x_per,PCG5_SSOR,PCG5_SSOR_x_per, &
&         BiCGStab5,BiCGStab5_x_per,CG55,CG55_x_per,PCG55_SSOR_x_per,PCG55_SSOR, GSeidel55, &
&         PCG55_SSOR_x_per_v2,PCG55_SSOR_x_per_v3, PCG55_SSOR_v2, PCG55_SSOR_v3, alloc_memory_for_Lsys
public :: CGtol, verb

integer :: verb
real(8), parameter :: pi=3.1415926535897932384626433832795
real(8) :: CGtol

real(8), allocatable :: rglob(:,:)
real(8), allocatable :: pglob(:,:)
real(8), allocatable :: zglob(:,:)
real(8), allocatable :: Auglob(:,:)

contains
!"Elemental" functions
!==================================================================================================================================
function A_dot_u(ni,nj,a,b,c,d,e,u)
!Escalar product [A]{v} for symmetric [A]
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: u
real(8), dimension(0:ni+1,0:nj+1) :: A_dot_u

integer :: i,j

A_dot_u=0.0d0 !initialization

!inner points
do j =1,nj
    do i=1,ni
        A_dot_u(i,j) = a(i,j)*u(i-1,j) &
    &                + b(i,j)*u(i,j-1) &
    &                + c(i,j)*u(i,j)   &
    &                + d(i,j)*u(i+1,j) &
    &                + e(i,j)*u(i,j+1)
    end do
end do

end function A_dot_u
!==================================================================================================================================

!==================================================================================================================================
function A_dot_u_x_per(ni,nj,a,b,c,d,e,u)
!Escalar product [A]{v} for symmetric [A]
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: u
real(8), dimension(0:ni+1,0:nj+1) :: A_dot_u_x_per

integer :: i,j

A_dot_u_x_per=0.0d0 !initialization

!inner points
do concurrent (i=2:ni-1, j=1:nj)
    A_dot_u_x_per(i,j) = a(i,j)*u(i-1,j) &
    &                  + b(i,j)*u(i,j-1) &
    &                  + c(i,j)*u(i,j)   &
    &                  + d(i,j)*u(i+1,j) &
    &                  + e(i,j)*u(i,j+1)
end do

!boundary points
i=1
do concurrent (j=1:nj)
    A_dot_u_x_per(i,j) = a(i,j)*u(ni,j)  &
    &                  + b(i,j)*u(i,j-1) &
    &                  + c(i,j)*u(i,j)   &
    &                  + d(i,j)*u(i+1,j) &
    &                  + e(i,j)*u(i,j+1)
end do

i=ni
do concurrent (j=1:nj)
    A_dot_u_x_per(i,j) = a(i,j)*u(i-1,j) &
    &                  + b(i,j)*u(i,j-1) &
    &                  + c(i,j)*u(i,j)   &
    &                  + d(i,j)*u(1,j)   &
    &                  + e(i,j)*u(i,j+1)
end do

end function A_dot_u_x_per
!==================================================================================================================================

!==================================================================================================================================
function A_dot_u_symm(ni,nj,a,b,c,u)
!Escalar product [A]{v} for symmetric [A]
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: u
real(8), dimension(0:ni+1,0:nj+1) :: A_dot_u_symm

integer :: i,j

A_dot_u_symm=0.0d0 !initialization

!inner points
do concurrent (i=1:ni, j=1:nj)
    A_dot_u_symm(i,j)=a(i,j)*(u(i-1,j)+u(i+1,j))   &
&                    +b(i,j)*(u(i,j-1)+u(i,j+1))   &
&                    +c(i,j)*u(i,j)
end do

end function A_dot_u_symm
!==================================================================================================================================

!==================================================================================================================================
function A_dot_u_symm_x_per(ni,nj,a,b,c,u)
!Escalar product [A]{v} for symmetric [A]
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: u
real(8), dimension(0:ni+1,0:nj+1) :: A_dot_u_symm_x_per

integer :: i,j

A_dot_u_symm_x_per=0.0d0 !initialization

!inner points
do concurrent (i=2:ni-1, j=1:nj)
    A_dot_u_symm_x_per(i,j)=a(i,j)*(u(i-1,j)+u(i+1,j))   &
&                          +b(i,j)*(u(i,j-1)+u(i,j+1))   &
&                          +c(i,j)*u(i,j)
end do

!boundary points
i=1
do concurrent (j=1:nj)
    A_dot_u_symm_x_per(i,j)=a(i,j)*(u(ni,j)+u(i+1,j))    &
&                          +b(i,j)*(u(i,j-1)+u(i,j+1))   &
&                          +c(i,j)*u(i,j)
end do

i=ni
do concurrent (j=1:nj)
    A_dot_u_symm_x_per(i,j)=a(i,j)*(u(i-1,j)+u(1,j))     &
&                          +b(i,j)*(u(i,j-1)+u(i,j+1))   &
&                          +c(i,j)*u(i,j)
end do

end function A_dot_u_symm_x_per
!==================================================================================================================================

!==================================================================================================================================
function u_dot_v(ni,nj,u,v)
!Perform {v} \cdot {u} (inner product)
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: u,v
real(8) :: u_dot_v
integer :: i,j

u_dot_v=0.0d0
do concurrent(i=1:ni,j=1:nj)
        u_dot_v=u_dot_v+u(i,j)*v(i,j)
end do

end function u_dot_v
!==================================================================================================================================

!==================================================================================================================================
pure function get_machine_eps()
!Get upper bound on the relative error due to rounding in floating point arithmetic of the computer ---> machine epsilon
implicit none
real(8) :: get_machine_eps

get_machine_eps=1.0d0
do 
    get_machine_eps=0.5*get_machine_eps
    if((1.0d0+get_machine_eps).eq.1.0d0) exit
end do

end function get_machine_eps
!==================================================================================================================================

!Main subroutines
!==================================================================================================================================
subroutine GSeidel55(ni,nj,a,b,c,d,e,f,u,tol,number_it)
implicit none
integer, intent(in) :: ni,nj
real, dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e    !explicit-shape dummy arrays
real, dimension(0:ni+1,0:nj+1), intent(in) :: f        !explicit-shape dummy arrays
real, dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real, optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: i,j,nit   !number of iterations
real :: alpha,beta,rr,local_tol !local auxiliary variables
real, allocatable :: r(:,:),p(:,:),Au(:,:) !automatic arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(Au(0:ni+1,0:nj+1))
r=0.0d0
Au=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Gauss Seidel working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')
Au=A_dot_u(ni,nj,a,b,c,d,e,u)
r(1:ni,1:nj)=f(1:ni,1:nj)-Au(1:ni,1:nj)
rr=u_dot_v(ni,nj,r,r)

if(abs(rr).lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! The GS residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=0
do
    do j=1,nj
        do i=1,ni
            u(i,j)=(f(i,j)-(a(i,j)*u(i-1,j)+b(i,j)*u(i,j-1)+d(i,j)*u(i+1,j)+e(i,j)*u(i,j+1)))/c(i,j)
        end do
    end do
    
    do j=nj,-1,1
        do i=ni,-1,1
            u(i,j)=(f(i,j)-(a(i,j)*u(i-1,j)+b(i,j)*u(i,j-1)+d(i,j)*u(i+1,j)+e(i,j)*u(i,j+1)))/c(i,j)
        end do
    end do
    
    Au=A_dot_u(ni,nj,a,b,c,d,e,u)
    r(1:ni,1:nj)=f(1:ni,1:nj)-Au(1:ni,1:nj)
    rr=u_dot_v(ni,nj,r,r)
    
    nit=nit+1
    if(abs(rr).lt.local_tol) exit
        
    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)


deallocate(r,Au)

return
end subroutine GSeidel55
!==================================================================================================================================

!==================================================================================================================================
subroutine CG5(ni,nj,a,b,c,f,u,tol,number_it)
!Conjugate Gradient Method for symmetric positve definite LS ([A]{u}={f}). No precontioning is used!
!
!The sparse structure of A must be compatible with a 5 points stencil discretization such that
!                 a(i,j)u(i-1,j) + b(i,j)u(i,j-1) + c(i,j)u(i,j) + a(i,j)u(i+1,j) + b(i,j)u(i,j+1) = f(i,j) (100)
!
!ni, nj: dimensions
!a,b,c: Matrix [A]'s components given in the form of equation (100)
!f: Vector {f}'s components (forcing vector) given in the form of equation (100)
!u: Solution vector. At the entree must be equal to the initial guess
!tol: tolerance. If not present in the argument list, then tol = twice the machine error
!
!Notation: 
!a*b = sum_{i=1,ni} sum_{j=1,nj} a(i,j)*b(i,j)
!A*v = a(i,j)v(i-1,j) + b(i,j)v(i,j-1) + c(i,j)v(i,j) + a(i,j)v(i+1,j) + b(i,j)v(i,j+1)
!
!
!Main algorithm
!1. r_0 := b - A*u_0
!2. p_0 := r_0
!3. from k := 0
!   3.1 alpha_k := (r_k*r_k)/(p_k*A*p_k)
!   3.2 u_{k+1} := u_k + \alpha_k p_k
!   3.3 r_{k+1} := r_k - \alpha_k A*p_k
!   3.4 if r_{k+1} is sufficiently small then exit loop 
!   3.5 beta_k := r_{k+1}*r_{k+1}/(r_k*r_k)
!   3.6 p_{k+1} := r_{k+1} + \beta_k p_k
!end 3
!
!The solution is u_{k+1}
!
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c    !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f    !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u !explicit-shape dummy arrays
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: nit   !number of iterations
real(8) :: alpha,beta,rr,local_tol !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),Au(:,:) !automatic arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(Au(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
Au=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')
Au=A_dot_u_symm(ni,nj,a,b,c,u) !--------------------->{Av} = [A]{u0}: A u_0 computation. Remember, u0 = u at the beginning
r(1:ni,1:nj)=f(1:ni,1:nj)-Au(1:ni,1:nj)  !----------->1. r_0 := b - A u_0
p(1:ni,1:nj)=r(1:ni,1:nj) !-------------------------->2. p_0 := r_0
rr=u_dot_v(ni,nj,r,r) !------------------------------>r_k \cdot z_k: for \alpha_k := (r_k^T r_k)/(p_k^T Ap_k) computation

if(abs(rr).lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! The CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=0
do
    Au=A_dot_u_symm(ni,nj,a,b,c,p) !------------------------->{Av} = [A]{p}: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation
    alpha = rr/u_dot_v(ni,nj,p,Au) !------------------------->3.1 \alpha_k := (r_k^T r_k)/(p_k^T Ap_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj)  !------>3.2 u_{k+1} := u_k + \alpha_k p_k
    r(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Au(1:ni,1:nj) !------>3.3 r_{k+1} := r_k - \alpha_k Ap_k
    
    nit=nit+1
    if(abs(rr).lt.local_tol) exit   !--------------------------->3.4 if r_{k+1} is sufficiently small then exit loop
    
    beta = u_dot_v(ni,nj,r,r)/rr !--------------------------->3.5 \beta_k := \frac{r_{k+1}^T r_{k+1}}{r_k^T r_k}
    p(1:ni,1:nj) = r(1:ni,1:nj) + beta*p(1:ni,1:nj) !-------->3.6 p_{k+1} := r_{k+1} + \beta_k p_k
    rr=u_dot_v(ni,nj,r,r) !---------------------------------->r_k \cdot r_k: for \alpha_k := (r_k^T r_k)/(p_k^T Ap_k) computation for the next iteration
    
    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,Au)

return
end subroutine CG5
!==================================================================================================================================

!==================================================================================================================================
subroutine CG5_x_per(ni,nj,a,b,c,f,u,tol,number_it)
!Conjugate Gradient Method for symmetric positve definite LS ([A]{u}={f}). No precontioning is used!
!
!The sparse structure of A must be compatible with a 5 points stencil discretization such that
!                 a(i,j)u(i-1,j) + b(i,j)u(i,j-1) + c(i,j)u(i,j) + a(i,j)u(i+1,j) + b(i,j)u(i,j+1) = f(i,j) (100)
!
!ni, nj: dimensions
!a,b,c: Matrix [A]'s components given in the form of equation (100)
!f: Vector {f}'s components (forcing vector) given in the form of equation (100)
!u: Solution vector. At the entree must be equal to the initial guess
!tol: tolerance. If not present in the argument list, then tol = twice the machine error
!
!Notation: 
!a*b = sum_{i=1,ni} sum_{j=1,nj} a(i,j)*b(i,j)
!A*v = a(i,j)v(i-1,j) + b(i,j)v(i,j-1) + c(i,j)v(i,j) + a(i,j)v(i+1,j) + b(i,j)v(i,j+1)
!
!
!Main algorithm
!1. r_0 := b - A*u_0
!2. p_0 := r_0
!3. from k := 0
!   3.1 alpha_k := (r_k*r_k)/(p_k*A*p_k)
!   3.2 u_{k+1} := u_k + \alpha_k p_k
!   3.3 r_{k+1} := r_k - \alpha_k A*p_k
!   3.4 if r_{k+1} is sufficiently small then exit loop 
!   3.5 beta_k := \frac{r_{k+1}*r_{k+1}}{r_k*r_k}
!   3.6 p_{k+1} := r_{k+1} + \beta_k p_k
!end 3
!
!The solution is u_{k+1}
!
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c    !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f    !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u !explicit-shape dummy arrays
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: nit   !number of iterations
real(8) :: alpha,beta,rr,local_tol !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),Au(:,:) !automatic arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(Au(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
Au=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if


if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')
Au=A_dot_u_symm_x_per(ni,nj,a,b,c,u) !--------------->{Av} = [A]{u0}: A u_0 computation. Remember, u0 = u at the beginning
r(1:ni,1:nj)=f(1:ni,1:nj)-Au(1:ni,1:nj)  !----------->1. r_0 := b - A u_0
p(1:ni,1:nj)=r(1:ni,1:nj) !-------------------------->2. p_0 := r_0
rr=u_dot_v(ni,nj,r,r) !------------------------------>r_k \cdot z_k: for \alpha_k := (r_k^T r_k)/(p_k^T Ap_k) computation

if(abs(rr).lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! The CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=0
do
    Au=A_dot_u_symm_x_per(ni,nj,a,b,c,p) !------------------->{Av} = [A]{p}: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation
    alpha = rr/u_dot_v(ni,nj,p,Au) !------------------------->3.1 \alpha_k := (r_k^T r_k)/(p_k^T Ap_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj)  !------>3.2 u_{k+1} := u_k + \alpha_k p_k
    r(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Au(1:ni,1:nj) !------>3.3 r_{k+1} := r_k - \alpha_k Ap_k
    
    nit=nit+1
    if(abs(rr).lt.local_tol) exit   !--------------------------->3.4 if r_{k+1} is sufficiently small then exit loop
    
    beta = u_dot_v(ni,nj,r,r)/rr !--------------------------->3.5 \beta_k := \frac{r_{k+1}^T r_{k+1}}{r_k^T r_k}
    p(1:ni,1:nj) = r(1:ni,1:nj) + beta*p(1:ni,1:nj) !-------->3.6 p_{k+1} := r_{k+1} + \beta_k p_k
    rr=u_dot_v(ni,nj,r,r) !---------------------------------->r_k \cdot r_k: for \alpha_k := (r_k^T r_k)/(p_k^T Ap_k) computation for the next iteration
    
    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,Au)

return
end subroutine CG5_x_per
!==================================================================================================================================

!==================================================================================================================================
subroutine PCG5_diag(ni,nj,a,b,c,f,u,tol,number_it)
!Conjugate Gradient Method for symmetric positve definite LS ([A]{u}={f}). Precontioning is carryed out 
!using the diagonal part of matrix [A]
!
!The sparse structure of A must be compatible with a 5 points stencil discretization such that
!                 a(i,j)u(i-1,j) + b(i,j)u(i,j-1) + c(i,j)u(i,j) + a(i,j)u(i+1,j) + b(i,j)u(i,j+1) = f(i,j) (100)
!
!f: Vector {f}'s components (forcing vector) given in the form of equation (100)
!u: Solution vector. At the entree must be equal to the initial guess
!tol: tolerance. If not present in the argument list, then tol = twice the machine error
!
!
!Notation: 
!a*b = sum_{i=1,ni} sum_{j=1,nj} a(i,j)*b(i,j)
!A*v = a(i,j)v(i-1,j) + b(i,j)v(i,j-1) + c(i,j)v(i,j) + a(i,j)v(i+1,j) + b(i,j)v(i,j+1)
!
!Main algorithm
!1. r_0 := b - A*u_0
!2. z_0 := M^{-1}*r_0  !M is the preconditioner matrix. In here, [M] = diag([A])
!3. p_0 := z_0
!
!4. from k := 0
!   4.1 alpha_k := (r_k*z_k)/(p_k*A*p_k)
!   4.2 u_{k+1} := u_k + \alpha_k p_k
!   4.3 r_{k+1} := r_k - \alpha_k A*p_k
!   4.4 if r_{k+1} is sufficiently small then exit loop 
!   4.5 z_{k+1} := M^{-1}*r_{k+1}
!   4.6 beta_k := \frac{z_{k+1}*r_{k+1}}{z_k*r_k}
!   4.7 p_{k+1} := z_{k+1} + \beta_k p_k
!end 4
!
!The solution is u_{k+1}
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f        !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: nit   !number of iterations
real(8) :: alpha,beta,local_tol,rz !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),z(:,:),Au(:,:) !arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(z(0:ni+1,0:nj+1))
allocate(Au(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
z=0.0d0
Au=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

Au= A_dot_u_symm(ni,nj,a,b,c,u) !--------------->{Av} = [A]{u0}: A u_0 computation. Remember, u0 = u at the beginning 
r(1:ni,1:nj)=f(1:ni,1:nj)-Au(1:ni,1:nj)  !------>1. r_0 := b - A u_0
z(1:ni,1:nj)=r(1:ni,1:nj)/c(1:ni,1:nj)  !------->2. z_0 := M^{-1}*r_0
p(1:ni,1:nj)=z(1:ni,1:nj) !--------------------->3. p_0 := z_0

rz=u_dot_v(ni,nj,r,z) !------------------------->r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation
 
if(abs(rz).lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! The CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=0
do
    Au=A_dot_u_symm(ni,nj,a,b,c,p) !------------------------------->{Av} = [A]{p}: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation
    alpha = rz/u_dot_v(ni,nj,p,Au) !------------------------------->4.1 \alpha_k := (r_k^T z_k)/(p_k^T Ap_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj)  !------------>4.2 u_{k+1} := u_k + \alpha_k p_k
    r(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Au(1:ni,1:nj) !------------>4.3 r_{k+1} := r_k - \alpha_k Ap_k

    nit=nit+1
    if(abs(rz).lt.local_tol) exit   !--------------------------------->4.4 if r_{k+1} is sufficiently small then exit loop
    
    z(1:ni,1:nj)=r(1:ni,1:nj)/c(1:ni,1:nj) !----------------------->4.5 z_{k+1} := M^{-1} r_{k+1}. Observe that, here, r is actually r_{k+1}
    beta = u_dot_v(ni,nj,r,z)/rz !--------------------------------->4.6 \beta_k := \frac{z_{k+1}^T r_{k+1}}{z_k^T r_k}
    p(1:ni,1:nj) = z(1:ni,1:nj) + beta*p(1:ni,1:nj) !-------------->4.7 p_{k+1} := z_{k+1} + \beta_k p_k
    rz=u_dot_v(ni,nj,r,z) !---------------------------------------->r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation for the next iteration
    
    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
    end if

end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,z,Au)
return
end subroutine PCG5_diag
!==================================================================================================================================

!==================================================================================================================================
subroutine PCG5_diag_x_per(ni,nj,a,b,c,f,u,tol,number_it)
!Conjugate Gradient Method for symmetric positve definite LS ([A]{u}={f}). Precontioning is carryed out 
!using the diagonal part of matrix [A]
!
!The sparse structure of A must be compatible with a 5 points stencil discretization such that
!                 a(i,j)u(i-1,j) + b(i,j)u(i,j-1) + c(i,j)u(i,j) + a(i,j)u(i+1,j) + b(i,j)u(i,j+1) = f(i,j) (100)
!
!f: Vector {f}'s components (forcing vector) given in the form of equation (100)
!u: Solution vector. At the entree must be equal to the initial guess
!tol: tolerance. If not present in the argument list, then tol = twice the machine error
!
!
!Notation: 
!a*b = sum_{i=1,ni} sum_{j=1,nj} a(i,j)*b(i,j)
!A*v = a(i,j)v(i-1,j) + b(i,j)v(i,j-1) + c(i,j)v(i,j) + a(i,j)v(i+1,j) + b(i,j)v(i,j+1)
!
!Main algorithm
!1. r_0 := b - A*u_0
!2. z_0 := M^{-1}*r_0  !M is the preconditioner matrix. In here, [M] = diag([A])
!3. p_0 := z_0
!
!4. from k := 0
!   4.1 alpha_k := (r_k*z_k)/(p_k*A*p_k)
!   4.2 u_{k+1} := u_k + \alpha_k p_k
!   4.3 r_{k+1} := r_k - \alpha_k A*p_k
!   4.4 if r_{k+1} is sufficiently small then exit loop 
!   4.5 z_{k+1} := M^{-1}*r_{k+1}
!   4.6 beta_k := \frac{z_{k+1}*r_{k+1}}{z_k*r_k}
!   4.7 p_{k+1} := z_{k+1} + \beta_k p_k
!end 4
!
!The solution is u_{k+1}
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f        !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: nit   !number of iterations
real(8) :: alpha,beta,local_tol,rz !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),z(:,:),Au(:,:) !arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(z(0:ni+1,0:nj+1))
allocate(Au(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
z=0.0d0
Au=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

Au= A_dot_u_symm_x_per(ni,nj,a,b,c,u) !--------->{Av} = [A]{u0}: A u_0 computation. Remember, u0 = u at the beginning 
r(1:ni,1:nj)=f(1:ni,1:nj)-Au(1:ni,1:nj)  !------>1. r_0 := b - A u_0
z(1:ni,1:nj)=r(1:ni,1:nj)/c(1:ni,1:nj)  !------->2. z_0 := M^{-1}*r_0
p(1:ni,1:nj)=z(1:ni,1:nj) !--------------------->3. p_0 := z_0

rz=u_dot_v(ni,nj,r,z) !------------------------->r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation
 
if(rz.lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! The CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=0
do
    Au=A_dot_u_symm_x_per(ni,nj,a,b,c,p) !------------------------->{Av} = [A]{p}: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation
    alpha = rz/u_dot_v(ni,nj,p,Au) !------------------------------->4.1 \alpha_k := (r_k^T z_k)/(p_k^T Ap_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj)  !------------>4.2 u_{k+1} := u_k + \alpha_k p_k
    r(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Au(1:ni,1:nj) !------------>4.3 r_{k+1} := r_k - \alpha_k Ap_k

    nit=nit+1
    if(abs(rz).lt.local_tol) exit   !--------------------------------->4.4 if r_{k+1} is sufficiently small then exit loop
    
    z(1:ni,1:nj)=r(1:ni,1:nj)/c(1:ni,1:nj) !----------------------->4.5 z_{k+1} := M^{-1} r_{k+1}. Observe that, here, r is actually r_{k+1}
    beta = u_dot_v(ni,nj,r,z)/rz !--------------------------------->4.6 \beta_k := \frac{z_{k+1}^T r_{k+1}}{z_k^T r_k}
    p(1:ni,1:nj) = z(1:ni,1:nj) + beta*p(1:ni,1:nj) !-------------->4.7 p_{k+1} := z_{k+1} + \beta_k p_k
    rz=u_dot_v(ni,nj,r,z) !---------------------------------------->r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation for the next iteration
    
    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
    end if

end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,z,Au)
return
end subroutine PCG5_diag_x_per
!==================================================================================================================================

!==================================================================================================================================
subroutine PCG5_SSOR(ni,nj,a,b,c,f,u,tol,h,number_it)
!Conjugate Gradient Method for symmetric positve definite LS ([A]{u}={f}). SSOR Preconditioning carried out as in
!Gen, et. al., Appl. Math. Mech. -Engl. Ed. 34(10), 1225-1236 (2013). Exactly as in section 2.2, equation (20).
!
!The sparse structure of A must be compatible with a 5 points stencil discretization such that
!                 a(i,j)u(i-1,j) + b(i,j)u(i,j-1) + c(i,j)u(i,j) + a(i,j)u(i+1,j) + b(i,j)u(i,j+1) = f(i,j) (100)
!
!nx, ny: dimensions
!a,b,c: Matrix [A]'s components given in the form of equation (100)
!f: Vector {f}'s components (forcing vector) given in the form of equation (100)
!u: Solution vector. At the entree must be equal to the initial guess
!tol: tolerance. If not present in the argument list, then tol = twice the machine error
!h: grid spacing for SSOR omega parameter estimation
!
!
!Notation: 
!a*b = sum_{i=1,ni} sum_{j=1,nj} a(i,j)*b(i,j)
!A*v = a(i,j)v(i-1,j) + b(i,j)v(i,j-1) + c(i,j)v(i,j) + a(i,j)v(i+1,j) + b(i,j)v(i,j+1)
!
!
!Main algorithm
!1. r_0 := b - A*u_0
!2. z_0 := M^{-1}*r_0  !M is the preconditioner matrix. 
!3. p_0 := z_0
!
!4. from k := 0
!   4.1 alpha_k := (r_k*z_k)/(p_k*A*p_k)
!   4.2 u_{k+1} := u_k + \alpha_k p_k
!   4.3 r_{k+1} := r_k - \alpha_k A*p_k
!   4.4 if r_{k+1} is sufficiently small then exit loop 
!   4.5 z_{k+1} := M^{-1}*r_{k+1}
!   4.6 beta_k := \frac{z_{k+1}*r_{k+1}}{z_k*r_k}
!   4.7 p_{k+1} := z_{k+1} + \beta_k p_k
!end 4
!
!The solution is u_{k+1}
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), intent(in) :: h
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: nit   !number of iterations
real(8) :: alpha,beta,local_tol,omega,rz_old,rz_new !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),z(:,:),Au(:,:) !automatic arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(z(0:ni+1,0:nj+1))
allocate(Au(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
z=0.0d0
Au=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')

omega=2.0d0/(1.0d0+pi*h)

Au = A_dot_u_symm(ni,nj,a,b,c,u) !Dot put this function directly in the next attribution. It may corrupt the pointers!
r(1:ni,1:nj) = f(1:ni,1:nj)-Au(1:ni,1:nj)!---------------->1. r_0 := f - A u_0
z = SSOR_precon(ni,nj,a,b,c,r,omega) !-------------------->2. z_0 := M^{-1} r_0
p(1:ni,1:nj) = z(1:ni,1:nj) !----------------------------->3. p_0 := z_0

rz_old=u_dot_v(ni,nj,r,z) !---------------------------> r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation

if(abs(rz_old).lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! The CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=1
do
    Au = A_dot_u_symm(ni,nj,a,b,c,p) !Do not put this function directly in the next attribution. It may corrupt the pointers!
    alpha = rz_old/u_dot_v(ni,nj,p,Au) !------------------------------------->4.1 \alpha_k := (r_k.z_k)/(p_k.A.p_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj)         !----------->4.2 u_{k+1} := u_k + \alpha_k p_k
    r(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Au(1:ni,1:nj)        !----------->4.3 r_{k+1} := r_k - \alpha_k Ap_k

    nit=nit+1    
    if(abs(rz_old).lt.local_tol) exit   !----------------------------------->4.4 if r_{k+1} is sufficiently small then exit loop
    
    z=SSOR_precon(ni,nj,a,b,c,r,omega) !-------------->4.5 z_{k+1} := M^{-1} r_{k+1}. Observe that, here, r is actually r_{k+1}
    rz_new=u_dot_v(ni,nj,r,z) !--------->r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation for the next iteration
    beta = rz_new/rz_old !-------------------------------------->4.6 \beta_k := \frac{z_{k+1}^T r_{k+1}}{z_k^T r_k}
    p(1:ni,1:nj) = z(1:ni,1:nj) + beta*p(1:ni,1:nj) !-------------------->4.7 p_{k+1} := z_{k+1} + \beta_k p_k
    rz_old=rz_new
                                                                                                                                !
    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz_old)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz_old)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,z,Au)
return
contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function SSOR_precon(ni,nj,a,b,c,r,omega)
!Porpouse: use SSOR to preconditioning CG. This algorithm might not be easily parallelized.
implicit none
integer, intent(in) :: ni,nj
real(8), intent(in) :: omega !SSOR parameter. Attention: might be necessary to optmized it for non-eliptic problems
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: r
real(8), dimension(0:ni+1,0:nj+1) :: SSOR_precon

integer :: i,j

SSOR_precon=0.0d0

!Lower SSOR
do j=1,nj
    do i=1,ni
        SSOR_precon(i,j) = ((-omega*(a(i,j)*SSOR_precon(i-1,j)+b(i,j)*SSOR_precon(i,j-1)))+omega*(2.0-omega)*r(i,j))/(c(i,j))
    end do
end do

!Upper SSOR
do j=nj,1,-1
    do i=ni,1,-1
        SSOR_precon(i,j) = (-omega*(a(i,j)*SSOR_precon(i+1,j)+b(i,j)*SSOR_precon(i,j+1))+c(i,j)*SSOR_precon(i,j))/c(i,j)
    end do
end do

end function SSOR_precon
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine PCG5_SSOR
!==================================================================================================================================

!==================================================================================================================================
subroutine PCG5_SSOR_x_per(ni,nj,a,b,c,f,u,tol,h,number_it)
!Conjugate Gradient Method for symmetric positve definite LS ([A]{u}={f}). SSOR Preconditioning carried out as in
!Gen, et. al., Appl. Math. Mech. -Engl. Ed. 34(10), 1225-1236 (2013). Exactly as in section 2.2, equation (20).
!
!The sparse structure of A must be compatible with a 5 points stencil discretization such that
!                 a(i,j)u(i-1,j) + b(i,j)u(i,j-1) + c(i,j)u(i,j) + a(i,j)u(i+1,j) + b(i,j)u(i,j+1) = f(i,j) (100)
!
!nx, ny: dimensions
!a,b,c: Matrix [A]'s components given in the form of equation (100)
!f: Vector {f}'s components (forcing vector) given in the form of equation (100)
!u: Solution vector. At the entree must be equal to the initial guess
!tol: tolerance. If not present in the argument list, then tol = twice the machine error
!h: grid spacing for SSOR omega parameter estimation
!
!
!Notation: 
!a*b = sum_{i=1,ni} sum_{j=1,nj} a(i,j)*b(i,j)
!A*v = a(i,j)v(i-1,j) + b(i,j)v(i,j-1) + c(i,j)v(i,j) + a(i,j)v(i+1,j) + b(i,j)v(i,j+1)
!
!
!Main algorithm
!1. r_0 := b - A*u_0
!2. z_0 := M^{-1}*r_0  !M is the preconditioner matrix. In here, [M] = diag([A])
!3. p_0 := z_0
!
!4. from k := 0
!   4.1 alpha_k := (r_k*z_k)/(p_k*A*p_k)
!   4.2 u_{k+1} := u_k + \alpha_k p_k
!   4.3 r_{k+1} := r_k - \alpha_k A*p_k
!   4.4 if r_{k+1} is sufficiently small then exit loop 
!   4.5 z_{k+1} := M^{-1}*r_{k+1}
!   4.6 beta_k := \frac{z_{k+1}*r_{k+1}}{z_k*r_k}
!   4.7 p_{k+1} := z_{k+1} + \beta_k p_k
!end 4
!
!The solution is u_{k+1}
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), intent(in) :: h
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: nit   !number of iterations
real(8) :: alpha,beta,local_tol,omega,rz !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),z(:,:),Au(:,:) !automatic arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(z(0:ni+1,0:nj+1))
allocate(Au(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
z=0.0d0
Au=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')

omega=2.0d0/(1.0d0+pi*h)

Au = A_dot_u_symm_x_per(ni,nj,a,b,c,u) !Dot put this function directly in the next attribution. It may corrupt the pointers!
r(1:ni,1:nj) = f(1:ni,1:nj)-Au(1:ni,1:nj)!---------------->1. r_0 := f - A u_0
z = SSOR_precon(ni,nj,a,b,c,r,omega) !-------------------->2. z_0 := M^{-1} r_0
p(1:ni,1:nj) = z(1:ni,1:nj) !----------------------------->3. p_0 := z_0

rz=u_dot_v(ni,nj,r,z) !-----------------------------------> r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation

if(abs(rz).lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! The CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=1
do
    Au = A_dot_u_symm_x_per(ni,nj,a,b,c,p) !Dot put this function directly in the next attribution. It may corrupt the pointers!
    alpha = rz/u_dot_v(ni,nj,p,Au) !------------------------------------->4.1 \alpha_k := (r_k.z_k)/(p_k.A.p_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj)         !----------->4.2 u_{k+1} := u_k + \alpha_k p_k
    r(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Au(1:ni,1:nj)        !----------->4.3 r_{k+1} := r_k - \alpha_k Ap_k

    nit=nit+1    
    if(abs(rz).lt.local_tol) exit   !--------------------------------------->4.4 if r_{k+1} is sufficiently small then exit loop
    
    z=SSOR_precon(ni,nj,a,b,c,r,omega) !--------------------------------->4.5 z_{k+1} := M^{-1} r_{k+1}. Observe that, here, r is actually r_{k+1}
    beta = u_dot_v(ni,nj,r,z)/rz !--------------------------------------->4.6 \beta_k := \frac{z_{k+1}^T r_{k+1}}{z_k^T r_k}
    p(1:ni,1:nj) = z(1:ni,1:nj) + beta*p(1:ni,1:nj) !-------------------->4.7 p_{k+1} := z_{k+1} + \beta_k p_k
    rz=u_dot_v(ni,nj,r,z) !---------------------------------------------->r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation for the next iteration

    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,z,Au)
return
contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function SSOR_precon(ni,nj,a,b,c,r,omega)
!Porpouse: use SSOR to preconditioning CG. This algorithm might not be easily parallelized.
implicit none
integer, intent(in) :: ni,nj
real(8), intent(in) :: omega !SSOR parameter. Attention: might be necessary to optmized it for non-eliptic problems
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: r
real(8), dimension(0:ni+1,0:nj+1) :: SSOR_precon

integer :: i,j

SSOR_precon=0.0d0

!Lower SSOR
do j=1,nj
    do i=1,ni
        SSOR_precon(i,j) = ((-omega*(a(i,j)*SSOR_precon(i-1,j)+b(i,j)*SSOR_precon(i,j-1)))+omega*(2.0-omega)*r(i,j))/(c(i,j))
    end do
end do

!Upper SSOR
do j=nj,1,-1
    do i=ni,1,-1
        SSOR_precon(i,j) = (-omega*(a(i,j)*SSOR_precon(i+1,j)+b(i,j)*SSOR_precon(i,j+1))+c(i,j)*SSOR_precon(i,j))/c(i,j)
    end do
end do

end function SSOR_precon
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine PCG5_SSOR_x_per
!==================================================================================================================================

!==================================================================================================================================
subroutine BiCGStab5(ni,nj,a,b,c,d,e,f,u,tol,number_it)
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e    !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f        !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it


integer :: nit   !number of iterations
real(8) :: alpha,beta,omega,rr_old,rr_new,local_tol !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),s(:,:),hatr(:,:),Ap(:,:),As(:,:) !automatic arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(s(0:ni+1,0:nj+1))
allocate(hatr(0:ni+1,0:nj+1))
allocate(Ap(0:ni+1,0:nj+1))
allocate(As(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
s=0.0d0
hatr=0.0d0
Ap=0.0d0
As=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        if(present(number_it))number_it=0
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'BiCGStab working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')
Ap=A_dot_u(ni,nj,a,b,c,d,e,u)
r(1:ni,1:nj)=f(1:ni,1:nj)-Ap(1:ni,1:nj)
hatr(1:ni,1:nj)=r(1:ni,1:nj)
p(1:ni,1:nj)=r(1:ni,1:nj)
rr_old=u_dot_v(ni,nj,hatr,r)

if(abs(rr_old).lt.local_tol) then
    if(verb.eq.2) write(*,*)'Caution!!! The BiCGStab residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=0
do
    Ap=A_dot_u(ni,nj,a,b,c,d,e,p)
    alpha = rr_old/u_dot_v(ni,nj,Ap,hatr)
    s(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Ap(1:ni,1:nj)
    As=A_dot_u(ni,nj,a,b,c,d,e,s)
    omega=u_dot_v(ni,nj,As,s)/u_dot_v(ni,nj,As,As)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj) + omega*s(1:ni,1:nj)
    r(1:ni,1:nj) = s(1:ni,1:nj) - omega*As(1:ni,1:nj)
    rr_new=u_dot_v(ni,nj,hatr,r)
    
    nit=nit+1
    if(abs(rr_new).lt.local_tol) exit
    
    beta = (rr_new*alpha)/(rr_old*omega)
    p(1:ni,1:nj) = r(1:ni,1:nj) + beta*(p(1:ni,1:nj) - omega*Ap(1:ni,1:nj))
    rr_old=rr_new
    
    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr_new)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr_new)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,s,hatr,Ap,As)

return
end subroutine BiCGStab5
!==================================================================================================================================

!==================================================================================================================================
subroutine BiCGStab5_x_per(ni,nj,a,b,c,d,e,f,u,tol,number_it)
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e    !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f        !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: nit   !number of iterations
real(8) :: alpha,beta,omega,rr_old,rr_new,local_tol !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),s(:,:),hatr(:,:),Ap(:,:),As(:,:) !automatic arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(s(0:ni+1,0:nj+1))
allocate(hatr(0:ni+1,0:nj+1))
allocate(Ap(0:ni+1,0:nj+1))
allocate(As(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
s=0.0d0
hatr=0.0d0
Ap=0.0d0
As=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        if(present(number_it))number_it=0
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'BiCGStab working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')
Ap=A_dot_u_x_per(ni,nj,a,b,c,d,e,u)
r(1:ni,1:nj)=f(1:ni,1:nj)-Ap(1:ni,1:nj)
hatr(1:ni,1:nj)=r(1:ni,1:nj)
p(1:ni,1:nj)=r(1:ni,1:nj)
rr_old=u_dot_v(ni,nj,hatr,r)

if(abs(rr_old).lt.local_tol) then
    if(verb.eq.2) write(*,*)'Caution!!! The BiCGStab residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=0
do
    Ap=A_dot_u_x_per(ni,nj,a,b,c,d,e,p)
    alpha = rr_old/u_dot_v(ni,nj,Ap,hatr)
    s(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Ap(1:ni,1:nj)
    As=A_dot_u_x_per(ni,nj,a,b,c,d,e,s)
    omega=u_dot_v(ni,nj,As,s)/u_dot_v(ni,nj,As,As)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj) + omega*s(1:ni,1:nj)
    r(1:ni,1:nj) = s(1:ni,1:nj) - omega*As(1:ni,1:nj)
    rr_new=u_dot_v(ni,nj,hatr,r)
    
    nit=nit+1
    if(abs(rr_new).lt.local_tol) exit
    
    beta = (rr_new*alpha)/(rr_old*omega)
    p(1:ni,1:nj) = r(1:ni,1:nj) + beta*(p(1:ni,1:nj) - omega*Ap(1:ni,1:nj))
    rr_old=rr_new
    
    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr_new)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr_new)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,s,hatr,Ap,As)

return
end subroutine BiCGStab5_x_per
!==================================================================================================================================

!==================================================================================================================================
subroutine CG55(ni,nj,a,b,c,d,e,f,u,tol,number_it)
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e    !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f        !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: nit   !number of iterations
real(8) :: alpha,beta,rr,local_tol !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),Au(:,:) !automatic arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(Au(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
Au=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        if(present(number_it))number_it=0
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')
Au=A_dot_u(ni,nj,a,b,c,d,e,u) !--------------------->{Av} = [A]{u0}: A u_0 computation. Remember, u0 = u at the beginning
r(1:ni,1:nj)=f(1:ni,1:nj)-Au(1:ni,1:nj)  !----------->1. r_0 := b - A u_0
p(1:ni,1:nj)=r(1:ni,1:nj) !-------------------------->2. p_0 := r_0
rr=u_dot_v(ni,nj,r,r) !------------------------------>r_k \cdot z_k: for \alpha_k := (r_k^T r_k)/(p_k^T Ap_k) computation

if(abs(rr).lt.local_tol) then
    if(verb.eq.2) write(*,*)'Caution!!! The CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=0
do
    Au=A_dot_u(ni,nj,a,b,c,d,e,p) !-------------------->{Av} = [A]{p}: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation
    alpha = rr/u_dot_v(ni,nj,p,Au) !------------------------->3.1 \alpha_k := (r_k^T r_k)/(p_k^T Ap_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj)  !------>3.2 u_{k+1} := u_k + \alpha_k p_k
    r(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Au(1:ni,1:nj) !------>3.3 r_{k+1} := r_k - \alpha_k Ap_k
    
    nit=nit+1
    if(abs(rr).lt.local_tol) exit   !--------------------------->3.4 if r_{k+1} is sufficiently small then exit loop
    
    beta = u_dot_v(ni,nj,r,r)/rr !--------------------------->3.5 \beta_k := \frac{r_{k+1}^T r_{k+1}}{r_k^T r_k}
    p(1:ni,1:nj) = r(1:ni,1:nj) + beta*p(1:ni,1:nj) !-------->3.6 p_{k+1} := r_{k+1} + \beta_k p_k
    rr=u_dot_v(ni,nj,r,r) !---------------------------------->r_k \cdot r_k: for \alpha_k := (r_k^T r_k)/(p_k^T Ap_k) computation for the next iteration
    
    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,Au)

return
end subroutine CG55
!==================================================================================================================================

!==================================================================================================================================
subroutine CG55_x_per(ni,nj,a,b,c,d,e,f,u,tol,number_it)
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e    !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f        !explicit-shape dummy arrays
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: nit   !number of iterations
real(8) :: alpha,beta,rr,local_tol !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),Au(:,:) !automatic arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(Au(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
Au=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')
Au=A_dot_u_x_per(ni,nj,a,b,c,d,e,u) !--------------------->{Av} = [A]{u0}: A u_0 computation. Remember, u0 = u at the beginning
r(1:ni,1:nj)=f(1:ni,1:nj)-Au(1:ni,1:nj)  !----------->1. r_0 := b - A u_0
p(1:ni,1:nj)=r(1:ni,1:nj) !-------------------------->2. p_0 := r_0
rr=u_dot_v(ni,nj,r,r) !------------------------------>r_k \cdot z_k: for \alpha_k := (r_k^T r_k)/(p_k^T Ap_k) computation

if(abs(rr).lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! The CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=0
do
    Au=A_dot_u_x_per(ni,nj,a,b,c,d,e,p) !-------------------->{Av} = [A]{p}: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation
    alpha = rr/u_dot_v(ni,nj,p,Au) !------------------------->3.1 \alpha_k := (r_k^T r_k)/(p_k^T Ap_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj)  !------>3.2 u_{k+1} := u_k + \alpha_k p_k
    r(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Au(1:ni,1:nj) !------>3.3 r_{k+1} := r_k - \alpha_k Ap_k
    
    nit=nit+1
    if(abs(rr).lt.local_tol) exit   !--------------------------->3.4 if r_{k+1} is sufficiently small then exit loop
    
    beta = u_dot_v(ni,nj,r,r)/rr !--------------------------->3.5 \beta_k := \frac{r_{k+1}^T r_{k+1}}{r_k^T r_k}
    p(1:ni,1:nj) = r(1:ni,1:nj) + beta*p(1:ni,1:nj) !-------->3.6 p_{k+1} := r_{k+1} + \beta_k p_k
    rr=u_dot_v(ni,nj,r,r) !---------------------------------->r_k \cdot r_k: for \alpha_k := (r_k^T r_k)/(p_k^T Ap_k) computation for the next iteration
    
    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rr)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,Au)

return
end subroutine CG55_x_per
!==================================================================================================================================

!==================================================================================================================================
subroutine PCG55_SSOR(ni,nj,a,b,c,d,e,f,u,tol,h,number_it)
!Conjugate Gradient Method for symmetric positve definite LS ([A]{u}={f}). SSOR Preconditioning carried out as in
!Gen, et. al., Appl. Math. Mech. -Engl. Ed. 34(10), 1225-1236 (2013). Exactly as in section 2.2, equation (20).
!
!The sparse structure of A must be compatible with a 5 points stencil discretization such that
!                 a(i,j)u(i-1,j) + b(i,j)u(i,j-1) + c(i,j)u(i,j) + d(i,j)u(i+1,j) + e(i,j)u(i,j+1) = f(i,j) (100)
!
!nx, ny: dimensions
!a,b,c: Matrix [A]'s components given in the form of equation (100)
!f: Vector {f}'s components (forcing vector) given in the form of equation (100)
!u: Solution vector. At the entree must be equal to the initial guess
!tol: tolerance. If not present in the argument list, then tol = twice the machine error
!h: grid spacing for SSOR omega parameter estimation
!
!
!Notation: 
!a*b = sum_{i=1,ni} sum_{j=1,nj} a(i,j)*b(i,j)
!A*v = a(i,j)v(i-1,j) + b(i,j)v(i,j-1) + c(i,j)v(i,j) + a(i,j)v(i+1,j) + b(i,j)v(i,j+1)
!
!
!Main algorithm
!1. r_0 := b - A*u_0
!2. z_0 := M^{-1}*r_0  !M is the preconditioner matrix. In here, [M] = diag([A])
!3. p_0 := z_0
!
!4. from k := 0
!   4.1 alpha_k := (r_k*z_k)/(p_k*A*p_k)
!   4.2 u_{k+1} := u_k + \alpha_k p_k
!   4.3 r_{k+1} := r_k - \alpha_k A*p_k
!   4.4 if r_{k+1} is sufficiently small then exit loop 
!   4.5 z_{k+1} := M^{-1}*r_{k+1}
!   4.6 beta_k := \frac{z_{k+1}*r_{k+1}}{z_k*r_k}
!   4.7 p_{k+1} := z_{k+1} + \beta_k p_k
!end 4
!
!The solution is u_{k+1}
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), intent(in) :: h
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: nit   !number of iterations
real(8) :: alpha,beta,local_tol,omega,rz !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),z(:,:),Au(:,:) !automatic arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(z(0:ni+1,0:nj+1))
allocate(Au(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
z=0.0d0
Au=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')

omega=2.0d0/(1.0d0+pi*h)

Au = A_dot_u(ni,nj,a,b,c,d,e,u) !Dot put this function directly in the next attribution. It may corrupt the pointers!
r(1:ni,1:nj) = f(1:ni,1:nj)-Au(1:ni,1:nj)!---------------->1. r_0 := f - A u_0
z = SSOR_precon(ni,nj,a,b,c,d,e,r,omega) !-------------------->2. z_0 := M^{-1} r_0
p(1:ni,1:nj) = z(1:ni,1:nj) !----------------------------->3. p_0 := z_0

rz=u_dot_v(ni,nj,r,z) !-----------------------------------> r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation

if(abs(rz).lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! The CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=1
do
    Au = A_dot_u(ni,nj,a,b,c,d,e,p) !Dot put this function directly in the next attribution. It may corrupt the pointers!
    alpha = rz/u_dot_v(ni,nj,p,Au) !------------------------------------->4.1 \alpha_k := (r_k.z_k)/(p_k.A.p_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj)         !----------->4.2 u_{k+1} := u_k + \alpha_k p_k
    r(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Au(1:ni,1:nj)        !----------->4.3 r_{k+1} := r_k - \alpha_k Ap_k

    nit=nit+1    
    if(abs(rz).lt.local_tol) exit   !--------------------------------------->4.4 if r_{k+1} is sufficiently small then exit loop
    
    z=SSOR_precon(ni,nj,a,b,c,d,e,r,omega) !--------------------------------->4.5 z_{k+1} := M^{-1} r_{k+1}. Observe that, here, r is actually r_{k+1}
    beta = u_dot_v(ni,nj,r,z)/rz !--------------------------------------->4.6 \beta_k := \frac{z_{k+1}^T r_{k+1}}{z_k^T r_k}
    p(1:ni,1:nj) = z(1:ni,1:nj) + beta*p(1:ni,1:nj) !-------------------->4.7 p_{k+1} := z_{k+1} + \beta_k p_k
    rz=u_dot_v(ni,nj,r,z) !---------------------------------------------->r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation for the next iteration

    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,z,Au)
return
contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function SSOR_precon(ni,nj,a,b,c,d,e,r,omega)
!Porpouse: use SSOR to preconditioning CG. This algorithm might not be easily parallelized.
implicit none
integer, intent(in) :: ni,nj
real(8), intent(in) :: omega !SSOR parameter. Attention: might be necessary to optmized it for non-eliptic problems
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: r
real(8), dimension(0:ni+1,0:nj+1) :: SSOR_precon

integer :: i,j

SSOR_precon=0.0d0

!Lower SSOR
do j=1,nj
    do i=1,ni
        SSOR_precon(i,j) = ((-omega*(a(i,j)*SSOR_precon(i-1,j)+b(i,j)*SSOR_precon(i,j-1)))+omega*(2.0d0-omega)*r(i,j))/(c(i,j))
    end do
end do

!Upper SSOR
do j=nj,1,-1
    do i=ni,1,-1
        SSOR_precon(i,j) = (-omega*(d(i,j)*SSOR_precon(i+1,j)+e(i,j)*SSOR_precon(i,j+1))+c(i,j)*SSOR_precon(i,j))/c(i,j)
    end do
end do

end function SSOR_precon
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine PCG55_SSOR
!==================================================================================================================================

!==================================================================================================================================
subroutine PCG55_SSOR_v2(ni,nj,a,b,c,d,e,f,u,tol,h,number_it)
!Conjugate Gradient Method for symmetric positve definite LS ([A]{u}={f}). SSOR Preconditioning carried out as in
!Gen, et. al., Appl. Math. Mech. -Engl. Ed. 34(10), 1225-1236 (2013). Exactly as in section 2.2, equation (20).
!
!The sparse structure of A must be compatible with a 5 points stencil discretization such that
!                 a(i,j)u(i-1,j) + b(i,j)u(i,j-1) + c(i,j)u(i,j) + d(i,j)u(i+1,j) + e(i,j)u(i,j+1) = f(i,j) (100)
!
!nx, ny: dimensions
!a,b,c: Matrix [A]'s components given in the form of equation (100)
!f: Vector {f}'s components (forcing vector) given in the form of equation (100)
!u: Solution vector. At the entree must be equal to the initial guess
!tol: tolerance. If not present in the argument list, then tol = twice the machine error
!h: grid spacing for SSOR omega parameter estimation
!
!
!Notation: 
!a*b = sum_{i=1,ni} sum_{j=1,nj} a(i,j)*b(i,j)
!A*v = a(i,j)v(i-1,j) + b(i,j)v(i,j-1) + c(i,j)v(i,j) + a(i,j)v(i+1,j) + b(i,j)v(i,j+1)
!
!
!Main algorithm
!1. r_0 := b - A*u_0
!2. z_0 := M^{-1}*r_0  !M is the preconditioner matrix. In here, [M] = diag([A])
!3. p_0 := z_0
!
!4. from k := 0
!   4.1 alpha_k := (r_k*z_k)/(p_k*A*p_k)
!   4.2 u_{k+1} := u_k + \alpha_k p_k
!   4.3 r_{k+1} := r_k - \alpha_k A*p_k
!   4.4 if r_{k+1} is sufficiently small then exit loop 
!   4.5 z_{k+1} := M^{-1}*r_{k+1}
!   4.6 beta_k := \frac{z_{k+1}*r_{k+1}}{z_k*r_k}
!   4.7 p_{k+1} := z_{k+1} + \beta_k p_k
!end 4
!
!The solution is u_{k+1}
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), intent(in) :: h
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: i,j,nit   !number of iterations
real(8) :: alpha,beta,local_tol,omega,rz,rz_old,pAp !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),z(:,:),Au(:,:) !automatic arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(z(0:ni+1,0:nj+1))
allocate(Au(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
z=0.0d0
Au=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')

omega=2.0d0/(1.0d0+pi*h)

!A dot u old style-----------------------------------------------------------------------------------------------------------------
Au=0.0d0
!inner points
do j=1,nj
    do i=1,ni
        Au(i,j) = a(i,j)*u(i-1,j) &
    &           + b(i,j)*u(i,j-1) &
    &           + c(i,j)*u(i,j)   &
    &           + d(i,j)*u(i+1,j) &
    &           + e(i,j)*u(i,j+1)
    end do
end do
!----------------------------------------------------------------------------------------------------------------------------------


r(1:ni,1:nj) = f(1:ni,1:nj)-Au(1:ni,1:nj)!---------------->1. r_0 := f - A u_0
!SSOR pre-conditioning----------->2. z_0 := M^{-1} r_0-----------------------------------------------------------------------------
z=0.0d0

!Lower SSOR
do j=1,nj
    do i=1,ni
        z(i,j) = ((-omega*(a(i,j)*z(i-1,j)+b(i,j)*z(i,j-1)))+omega*(2.0-omega)*r(i,j))/(c(i,j))
    end do
end do

!Upper SSOR
do j=nj,1,-1
    do i=ni,1,-1
        z(i,j) = (-omega*(d(i,j)*z(i+1,j)+e(i,j)*z(i,j+1))+c(i,j)*z(i,j))/c(i,j)
    end do
end do
!----------------------------------------------------------------------------------------------------------------------------------

p(1:ni,1:nj) = z(1:ni,1:nj) !----------------------------->3. p_0 := z_0

!r dot z old style ----------------> r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation---------------------------
rz=0.0d0
do j=1,nj
    do i=1,ni
        rz=rz+r(i,j)*z(i,j)
    end do
end do
!----------------------------------------------------------------------------------------------------------------------------------
rz_old=rz

if(abs(rz).lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! The CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=1
do

    !A dot u old style-----------------------------------------------------------------------------------------------------------------
    Au=0.0d0
    !inner points
    do j=1,nj
        do i=1,ni
            Au(i,j) = a(i,j)*p(i-1,j) &
                &       + b(i,j)*p(i,j-1) &
                &       + c(i,j)*p(i,j)   &
                &       + d(i,j)*p(i+1,j) &
                &       + e(i,j)*p(i,j+1)
            end do
        end do
    !----------------------------------------------------------------------------------------------------------------------------------

    !p dot Ap (Au is just the name of the vector. Here it means Ap!!) old style -----> r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation
    pAp=0.0d0
    do j=1,nj
        do i=1,ni
            pAp=pAp+p(i,j)*Au(i,j)
        end do
    end do
    !----------------------------------------------------------------------------------------------------------------------------------

    alpha = rz/pAp                 !------------------------------------->4.1 \alpha_k := (r_k.z_k)/(p_k.A.p_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj)         !----------->4.2 u_{k+1} := u_k + \alpha_k p_k
    r(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Au(1:ni,1:nj)        !----------->4.3 r_{k+1} := r_k - \alpha_k Ap_k

    nit=nit+1    
    if(abs(rz).lt.local_tol) exit   !--------------------------------------->4.4 if r_{k+1} is sufficiently small then exit loop
    
    !SSOR pre-conditioning----------->4.5 z_{k+1} := M^{-1} r_{k+1}. Observe that, here, r is actually r_{k+1}--------------------------
    z=0.0d0

    !Lower SSOR
    do j=1,nj
        do i=1,ni
            z(i,j) = ((-omega*(a(i,j)*z(i-1,j)+b(i,j)*z(i,j-1)))+omega*(2.0-omega)*r(i,j))/(c(i,j))
        end do
    end do

    !Upper SSOR
    do j=nj,1,-1
        do i=ni,1,-1
            z(i,j) = (-omega*(d(i,j)*z(i+1,j)+e(i,j)*z(i,j+1))+c(i,j)*z(i,j))/c(i,j)
        end do
    end do
    !----------------------------------------------------------------------------------------------------------------------------------
    !r dot z old style ----------------------------------------------------------------------------------------------------------------
    rz=0.0d0
    do j=1,nj
        do i=1,ni
            rz=rz+r(i,j)*z(i,j)
        end do
    end do
    !----------------------------------------------------------------------------------------------------------------------------------
    beta = rz/rz_old !--------------------------------------------------->4.6 \beta_k := \frac{z_{k+1}^T r_{k+1}}{z_k^T r_k}
    p(1:ni,1:nj) = z(1:ni,1:nj) + beta*p(1:ni,1:nj) !-------------------->4.7 p_{k+1} := z_{k+1} + \beta_k p_k
    rz_old=rz !---------------------------------------------->r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation for the next iteration

    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,z,Au)
return
end subroutine PCG55_SSOR_v2
!==================================================================================================================================

!==================================================================================================================================
subroutine PCG55_SSOR_v3(ni,nj,a,b,c,d,e,f,u,tol,h,number_it)
!Conjugate Gradient Method for symmetric positve definite LS ([A]{u}={f}). SSOR Preconditioning carried out as in
!Gen, et. al., Appl. Math. Mech. -Engl. Ed. 34(10), 1225-1236 (2013). Exactly as in section 2.2, equation (20).
!
!The sparse structure of A must be compatible with a 5 points stencil discretization such that
!                 a(i,j)u(i-1,j) + b(i,j)u(i,j-1) + c(i,j)u(i,j) + d(i,j)u(i+1,j) + e(i,j)u(i,j+1) = f(i,j) (100)
!
!nx, ny: dimensions
!a,b,c: Matrix [A]'s components given in the form of equation (100)
!f: Vector {f}'s components (forcing vector) given in the form of equation (100)
!u: Solution vector. At the entree must be equal to the initial guess
!tol: tolerance. If not present in the argument list, then tol = twice the machine error
!h: grid spacing for SSOR omega parameter estimation
!
!
!Notation: 
!a*b = sum_{i=1,ni} sum_{j=1,nj} a(i,j)*b(i,j)
!A*v = a(i,j)v(i-1,j) + b(i,j)v(i,j-1) + c(i,j)v(i,j) + a(i,j)v(i+1,j) + b(i,j)v(i,j+1)
!
!
!Main algorithm
!1. r_0 := b - A*u_0
!2. z_0 := M^{-1}*r_0  !M is the preconditioner matrix. In here, [M] = diag([A])
!3. p_0 := z_0
!
!4. from k := 0
!   4.1 alpha_k := (r_k*z_k)/(p_k*A*p_k)
!   4.2 u_{k+1} := u_k + \alpha_k p_k
!   4.3 r_{k+1} := r_k - \alpha_k A*p_k
!   4.4 if r_{k+1} is sufficiently small then exit loop 
!   4.5 z_{k+1} := M^{-1}*r_{k+1}
!   4.6 beta_k := \frac{z_{k+1}*r_{k+1}}{z_k*r_k}
!   4.7 p_{k+1} := z_{k+1} + \beta_k p_k
!end 4
!
!The solution is u_{k+1}
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), intent(in) :: h
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: i,j,nit   !number of iterations
real(8) :: alpha,beta,local_tol,omega,rz,rz_old,pAp !local auxiliary variables

rglob=0.0d0
pglob=0.0d0
zglob=0.0d0
Auglob=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')

omega=2.0d0/(1.0d0+pi*h)

!A dot u old style-----------------------------------------------------------------------------------------------------------------
Auglob=0.0d0
!inner points
do j=1,nj
    do i=1,ni
        Auglob(i,j) = a(i,j)*u(i-1,j) &
    &           + b(i,j)*u(i,j-1) &
    &           + c(i,j)*u(i,j)   &
    &           + d(i,j)*u(i+1,j) &
    &           + e(i,j)*u(i,j+1)
    end do
end do
!----------------------------------------------------------------------------------------------------------------------------------


rglob(1:ni,1:nj) = f(1:ni,1:nj)-Auglob(1:ni,1:nj)!---------------->1. r_0 := f - A u_0
!SSOR pre-conditioning----------->2. z_0 := M^{-1} r_0-----------------------------------------------------------------------------
zglob=0.0d0

!Lower SSOR
do j=1,nj
    do i=1,ni
        zglob(i,j) = ((-omega*(a(i,j)*zglob(i-1,j)+b(i,j)*zglob(i,j-1)))+omega*(2.0-omega)*rglob(i,j))/(c(i,j))
    end do
end do

!Upper SSOR
do j=nj,1,-1
    do i=ni,1,-1
        zglob(i,j) = (-omega*(d(i,j)*zglob(i+1,j)+e(i,j)*zglob(i,j+1))+c(i,j)*zglob(i,j))/c(i,j)
    end do
end do
!----------------------------------------------------------------------------------------------------------------------------------

pglob(1:ni,1:nj) = zglob(1:ni,1:nj) !----------------------------->3. p_0 := z_0

!r dot z old style ----------------> r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation---------------------------
rz=0.0d0
do j=1,nj
    do i=1,ni
        rz=rz+rglob(i,j)*zglob(i,j)
    end do
end do
!----------------------------------------------------------------------------------------------------------------------------------
rz_old=rz

if(abs(rz).lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! The CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=1
do

    !A dot u old style-----------------------------------------------------------------------------------------------------------------
    Auglob=0.0d0
    !inner points
    do j=1,nj
        do i=1,ni
            Auglob(i,j) = a(i,j)*pglob(i-1,j) &
                &       + b(i,j)*pglob(i,j-1) &
                &       + c(i,j)*pglob(i,j)   &
                &       + d(i,j)*pglob(i+1,j) &
                &       + e(i,j)*pglob(i,j+1)
            end do
        end do
    !----------------------------------------------------------------------------------------------------------------------------------

    !p dot Ap (Au is just the name of the vector. Here it means Ap!!) old style -----> r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation
    pAp=0.0d0
    do j=1,nj
        do i=1,ni
            pAp=pAp+pglob(i,j)*Auglob(i,j)
        end do
    end do
    !----------------------------------------------------------------------------------------------------------------------------------

    alpha = rz/pAp                 !------------------------------------->4.1 \alpha_k := (r_k.z_k)/(p_k.A.p_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*pglob(1:ni,1:nj)         !----------->4.2 u_{k+1} := u_k + \alpha_k p_k
    rglob(1:ni,1:nj) = rglob(1:ni,1:nj) - alpha*Auglob(1:ni,1:nj)        !----------->4.3 r_{k+1} := r_k - \alpha_k Ap_k

    nit=nit+1    
    if(abs(rz).lt.local_tol) exit   !--------------------------------------->4.4 if r_{k+1} is sufficiently small then exit loop
    
    !SSOR pre-conditioning----------->4.5 z_{k+1} := M^{-1} r_{k+1}. Observe that, here, r is actually r_{k+1}--------------------------
    zglob=0.0d0

    !Lower SSOR
    do j=1,nj
        do i=1,ni
            zglob(i,j) = ((-omega*(a(i,j)*zglob(i-1,j)+b(i,j)*zglob(i,j-1)))+omega*(2.0-omega)*rglob(i,j))/(c(i,j))
        end do
    end do

    !Upper SSOR
    do j=nj,1,-1
        do i=ni,1,-1
            zglob(i,j) = (-omega*(d(i,j)*zglob(i+1,j)+e(i,j)*zglob(i,j+1))+c(i,j)*zglob(i,j))/c(i,j)
        end do
    end do
    !----------------------------------------------------------------------------------------------------------------------------------
    !r dot z old style ----------------------------------------------------------------------------------------------------------------
    rz=0.0d0
    do j=1,nj
        do i=1,ni
            rz=rz+rglob(i,j)*zglob(i,j)
        end do
    end do
    !----------------------------------------------------------------------------------------------------------------------------------
    beta = rz/rz_old !--------------------------------------------------->4.6 \beta_k := \frac{z_{k+1}^T r_{k+1}}{z_k^T r_k}
    pglob(1:ni,1:nj) = zglob(1:ni,1:nj) + beta*pglob(1:ni,1:nj) !-------------------->4.7 p_{k+1} := z_{k+1} + \beta_k p_k
    rz_old=rz !---------------------------------------------->r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation for the next iteration

    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

return
end subroutine PCG55_SSOR_v3
!==================================================================================================================================

!==================================================================================================================================
subroutine PCG55_SSOR_x_per(ni,nj,a,b,c,d,e,f,u,tol,h,number_it)
!Conjugate Gradient Method for symmetric positve definite LS ([A]{u}={f}). SSOR Preconditioning carried out as in
!Gen, et. al., Appl. Math. Mech. -Engl. Ed. 34(10), 1225-1236 (2013). Exactly as in section 2.2, equation (20).
!
!The sparse structure of A must be compatible with a 5 points stencil discretization such that
!                 a(i,j)u(i-1,j) + b(i,j)u(i,j-1) + c(i,j)u(i,j) + d(i,j)u(i+1,j) + e(i,j)u(i,j+1) = f(i,j) (100)
!
!nx, ny: dimensions
!a,b,c: Matrix [A]'s components given in the form of equation (100)
!f: Vector {f}'s components (forcing vector) given in the form of equation (100)
!u: Solution vector. At the entree must be equal to the initial guess
!tol: tolerance. If not present in the argument list, then tol = twice the machine error
!h: grid spacing for SSOR omega parameter estimation
!
!
!Notation: 
!a*b = sum_{i=1,ni} sum_{j=1,nj} a(i,j)*b(i,j)
!A*v = a(i,j)v(i-1,j) + b(i,j)v(i,j-1) + c(i,j)v(i,j) + a(i,j)v(i+1,j) + b(i,j)v(i,j+1)
!
!
!Main algorithm
!1. r_0 := b - A*u_0
!2. z_0 := M^{-1}*r_0  !M is the preconditioner matrix. In here, [M] = diag([A])
!3. p_0 := z_0
!
!4. from k := 0
!   4.1 alpha_k := (r_k*z_k)/(p_k*A*p_k)
!   4.2 u_{k+1} := u_k + \alpha_k p_k
!   4.3 r_{k+1} := r_k - \alpha_k A*p_k
!   4.4 if r_{k+1} is sufficiently small then exit loop 
!   4.5 z_{k+1} := M^{-1}*r_{k+1}
!   4.6 beta_k := \frac{z_{k+1}*r_{k+1}}{z_k*r_k}
!   4.7 p_{k+1} := z_{k+1} + \beta_k p_k
!end 4
!
!The solution is u_{k+1}
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), intent(in) :: h
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: nit   !number of iterations
real(8) :: alpha,beta,local_tol,omega,rz !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),z(:,:),Au(:,:) !automatic arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(z(0:ni+1,0:nj+1))
allocate(Au(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
z=0.0d0
Au=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')

omega=2.0d0/(1.0d0+pi*h)

Au = A_dot_u_x_per(ni,nj,a,b,c,d,e,u) !Dot put this function directly in the next attribution. It may corrupt the pointers!
r(1:ni,1:nj) = f(1:ni,1:nj)-Au(1:ni,1:nj)!---------------->1. r_0 := f - A u_0
z = SSOR_precon(ni,nj,a,b,c,d,e,r,omega) !-------------------->2. z_0 := M^{-1} r_0
p(1:ni,1:nj) = z(1:ni,1:nj) !----------------------------->3. p_0 := z_0

rz=u_dot_v(ni,nj,r,z) !-----------------------------------> r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation

if(abs(rz).lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! The CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
   return
end if

nit=1
do
    Au = A_dot_u_x_per(ni,nj,a,b,c,d,e,p) !Dot put this function directly in the next attribution. It may corrupt the pointers!
    alpha = rz/u_dot_v(ni,nj,p,Au) !------------------------------------->4.1 \alpha_k := (r_k.z_k)/(p_k.A.p_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj)         !----------->4.2 u_{k+1} := u_k + \alpha_k p_k
    r(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Au(1:ni,1:nj)        !----------->4.3 r_{k+1} := r_k - \alpha_k Ap_k

    nit=nit+1    
    if(abs(rz).lt.local_tol) exit   !--------------------------------------->4.4 if r_{k+1} is sufficiently small then exit loop
    
    z=SSOR_precon(ni,nj,a,b,c,d,e,r,omega) !--------------------------------->4.5 z_{k+1} := M^{-1} r_{k+1}. Observe that, here, r is actually r_{k+1}
    beta = u_dot_v(ni,nj,r,z)/rz !--------------------------------------->4.6 \beta_k := \frac{z_{k+1}^T r_{k+1}}{z_k^T r_k}
    p(1:ni,1:nj) = z(1:ni,1:nj) + beta*p(1:ni,1:nj) !-------------------->4.7 p_{k+1} := z_{k+1} + \beta_k p_k
    rz=u_dot_v(ni,nj,r,z) !---------------------------------------------->r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation for the next iteration

    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,z,Au)
return
contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function SSOR_precon(ni,nj,a,b,c,d,e,r,omega)
!Porpouse: use SSOR to preconditioning CG. This algorithm might not be easily parallelized.
implicit none
integer, intent(in) :: ni,nj
real(8), intent(in) :: omega !SSOR parameter. Attention: might be necessary to optmized it for non-eliptic problems
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: r
real(8), dimension(0:ni+1,0:nj+1) :: SSOR_precon

integer :: i,j

SSOR_precon=0.0d0

!Lower SSOR
do j=1,nj
    do i=1,ni
        SSOR_precon(i,j) = ((-omega*(a(i,j)*SSOR_precon(i-1,j)+b(i,j)*SSOR_precon(i,j-1)))+omega*(2.0-omega)*r(i,j))/(c(i,j))
    end do
end do

!Upper SSOR
do j=nj,1,-1
    do i=ni,1,-1
        SSOR_precon(i,j) = (-omega*(d(i,j)*SSOR_precon(i+1,j)+e(i,j)*SSOR_precon(i,j+1))+c(i,j)*SSOR_precon(i,j))/c(i,j)
    end do
end do

end function SSOR_precon
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine PCG55_SSOR_x_per
!==================================================================================================================================

!==================================================================================================================================
subroutine PCG55_SSOR_x_per_v2(ni,nj,a,b,c,d,e,f,u,tol,h,number_it)
!Conjugate Gradient Method for symmetric positve definite LS ([A]{u}={f}). SSOR Preconditioning carried out as in
!Gen, et. al., Appl. Math. Mech. -Engl. Ed. 34(10), 1225-1236 (2013). Exactly as in section 2.2, equation (20).
!
!The sparse structure of A must be compatible with a 5 points stencil discretization such that
!                 a(i,j)u(i-1,j) + b(i,j)u(i,j-1) + c(i,j)u(i,j) + d(i,j)u(i+1,j) + e(i,j)u(i,j+1) = f(i,j) (100)
!
!nx, ny: dimensions
!a,b,c: Matrix [A]'s components given in the form of equation (100)
!f: Vector {f}'s components (forcing vector) given in the form of equation (100)
!u: Solution vector. At the entree must be equal to the initial guess
!tol: tolerance. If not present in the argument list, then tol = twice the machine error
!h: grid spacing for SSOR omega parameter estimation
!
!
!Notation: 
!a*b = sum_{i=1,ni} sum_{j=1,nj} a(i,j)*b(i,j)
!A*v = a(i,j)v(i-1,j) + b(i,j)v(i,j-1) + c(i,j)v(i,j) + a(i,j)v(i+1,j) + b(i,j)v(i,j+1)
!
!
!Main algorithm
!1. r_0 := b - A*u_0
!2. z_0 := M^{-1}*r_0  !M is the preconditioner matrix. In here, [M] = diag([A])
!3. p_0 := z_0
!
!4. from k := 0
!   4.1 alpha_k := (r_k*z_k)/(p_k*A*p_k)
!   4.2 u_{k+1} := u_k + \alpha_k p_k
!   4.3 r_{k+1} := r_k - \alpha_k A*p_k
!   4.4 if r_{k+1} is sufficiently small then exit loop 
!   4.5 z_{k+1} := M^{-1}*r_{k+1}
!   4.6 beta_k := \frac{z_{k+1}*r_{k+1}}{z_k*r_k}
!   4.7 p_{k+1} := z_{k+1} + \beta_k p_k
!end 4
!
!The solution is u_{k+1}
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), intent(in) :: h
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: i,j,nit   !number of iterations
real(8) :: alpha,beta,local_tol,omega,rz,rz_old,pAp !local auxiliary variables
real(8), allocatable :: r(:,:),p(:,:),z(:,:),Au(:,:) !automatic arrays for local use only

allocate(r(0:ni+1,0:nj+1))
allocate(p(0:ni+1,0:nj+1))
allocate(z(0:ni+1,0:nj+1))
allocate(Au(0:ni+1,0:nj+1))
r=0.0d0
p=0.0d0
z=0.0d0
Au=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')

omega=2.0d0/(1.0d0+pi*h)

!A dot u old style-----------------------------------------------------------------------------------------------------------------
Au=0.0d0
!inner points
do j=1,nj
    do i=2,ni-1
        Au(i,j) = a(i,j)*u(i-1,j) &
    &           + b(i,j)*u(i,j-1) &
    &           + c(i,j)*u(i,j)   &
    &           + d(i,j)*u(i+1,j) &
    &           + e(i,j)*u(i,j+1)
    end do
end do

!boundary points
do j=1,nj
    Au(1,j)  = a(1,j)*u(ni,j)  &
    &        + b(1,j)*u(1,j-1) &
    &        + c(1,j)*u(1,j)   &
    &        + d(1,j)*u(2,j) &
    &        + e(1,j)*u(1,j+1)
    Au(ni,j) = a(ni,j)*u(ni-1,j) &
    &        + b(ni,j)*u(ni,j-1) &
    &        + c(ni,j)*u(ni,j)   &
    &        + d(ni,j)*u(1,j)    &
    &        + e(ni,j)*u(ni,j+1)
end do
!----------------------------------------------------------------------------------------------------------------------------------


r(1:ni,1:nj) = f(1:ni,1:nj)-Au(1:ni,1:nj)!---------------->1. r_0 := f - A u_0
!SSOR pre-conditioning----------->2. z_0 := M^{-1} r_0-----------------------------------------------------------------------------
z=0.0d0

!Lower SSOR
do j=1,nj
    do i=1,ni
        z(i,j) = ((-omega*(a(i,j)*z(i-1,j)+b(i,j)*z(i,j-1)))+omega*(2.0-omega)*r(i,j))/(c(i,j))
    end do
end do

!Upper SSOR
do j=nj,1,-1
    do i=ni,1,-1
        z(i,j) = (-omega*(d(i,j)*z(i+1,j)+e(i,j)*z(i,j+1))+c(i,j)*z(i,j))/c(i,j)
    end do
end do
!----------------------------------------------------------------------------------------------------------------------------------

p(1:ni,1:nj) = z(1:ni,1:nj) !----------------------------->3. p_0 := z_0

!r dot z old style ----------------> r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation---------------------------
rz=0.0d0
do j=1,nj
    do i=1,ni
        rz=rz+r(i,j)*z(i,j)
    end do
end do
!----------------------------------------------------------------------------------------------------------------------------------
rz_old=rz

if(abs(rz).lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! The CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=1
do

    !A dot u old style-----------------------------------------------------------------------------------------------------------------
    Au=0.0d0
    !inner points
    do j=1,nj
        do i=2,ni-1
            Au(i,j) = a(i,j)*p(i-1,j) &
                &       + b(i,j)*p(i,j-1) &
                &       + c(i,j)*p(i,j)   &
                &       + d(i,j)*p(i+1,j) &
                &       + e(i,j)*p(i,j+1)
            end do
        end do
        
        !boundary points
        do j=1,nj
            Au(1,j)  = a(1,j)*p(ni,j)  &
            &        + b(1,j)*p(1,j-1) &
            &        + c(1,j)*p(1,j)   &
            &        + d(1,j)*p(2,j) &
            &        + e(1,j)*p(1,j+1)
            Au(ni,j) = a(ni,j)*p(ni-1,j) &
            &        + b(ni,j)*p(ni,j-1) &
            &        + c(ni,j)*p(ni,j)   &
            &        + d(ni,j)*p(1,j)    &
            &        + e(ni,j)*p(ni,j+1)
        end do
    !----------------------------------------------------------------------------------------------------------------------------------

    !p dot Ap (Au is just the name of the vector. Here it means Ap!!) old style -----> r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation
    pAp=0.0d0
    do j=1,nj
        do i=1,ni
            pAp=pAp+p(i,j)*Au(i,j)
        end do
    end do
    !----------------------------------------------------------------------------------------------------------------------------------

    alpha = rz/pAp                 !------------------------------------->4.1 \alpha_k := (r_k.z_k)/(p_k.A.p_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*p(1:ni,1:nj)         !----------->4.2 u_{k+1} := u_k + \alpha_k p_k
    r(1:ni,1:nj) = r(1:ni,1:nj) - alpha*Au(1:ni,1:nj)        !----------->4.3 r_{k+1} := r_k - \alpha_k Ap_k

    nit=nit+1    
    if(abs(rz).lt.local_tol) exit   !--------------------------------------->4.4 if r_{k+1} is sufficiently small then exit loop
    
    !SSOR pre-conditioning----------->4.5 z_{k+1} := M^{-1} r_{k+1}. Observe that, here, r is actually r_{k+1}--------------------------
    z=0.0d0

    !Lower SSOR
    do j=1,nj
        do i=1,ni
            z(i,j) = ((-omega*(a(i,j)*z(i-1,j)+b(i,j)*z(i,j-1)))+omega*(2.0-omega)*r(i,j))/(c(i,j))
        end do
    end do

    !Upper SSOR
    do j=nj,1,-1
        do i=ni,1,-1
            z(i,j) = (-omega*(d(i,j)*z(i+1,j)+e(i,j)*z(i,j+1))+c(i,j)*z(i,j))/c(i,j)
        end do
    end do
    !----------------------------------------------------------------------------------------------------------------------------------
    !r dot z old style ----------------------------------------------------------------------------------------------------------------
    rz=0.0d0
    do j=1,nj
        do i=1,ni
            rz=rz+r(i,j)*z(i,j)
        end do
    end do
    !----------------------------------------------------------------------------------------------------------------------------------
    beta = rz/rz_old !--------------------------------------------------->4.6 \beta_k := \frac{z_{k+1}^T r_{k+1}}{z_k^T r_k}
    p(1:ni,1:nj) = z(1:ni,1:nj) + beta*p(1:ni,1:nj) !-------------------->4.7 p_{k+1} := z_{k+1} + \beta_k p_k
    rz_old=rz !---------------------------------------------->r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation for the next iteration

    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

deallocate(r,p,z,Au)
return
end subroutine PCG55_SSOR_x_per_v2
!==================================================================================================================================

!==================================================================================================================================
subroutine PCG55_SSOR_x_per_v3(ni,nj,a,b,c,d,e,f,u,tol,h,number_it)
!Conjugate Gradient Method for symmetric positve definite LS ([A]{u}={f}). SSOR Preconditioning carried out as in
!Gen, et. al., Appl. Math. Mech. -Engl. Ed. 34(10), 1225-1236 (2013). Exactly as in section 2.2, equation (20).
!
!The sparse structure of A must be compatible with a 5 points stencil discretization such that
!                 a(i,j)u(i-1,j) + b(i,j)u(i,j-1) + c(i,j)u(i,j) + d(i,j)u(i+1,j) + e(i,j)u(i,j+1) = f(i,j) (100)
!
!nx, ny: dimensions
!a,b,c: Matrix [A]'s components given in the form of equation (100)
!f: Vector {f}'s components (forcing vector) given in the form of equation (100)
!u: Solution vector. At the entree must be equal to the initial guess
!tol: tolerance. If not present in the argument list, then tol = twice the machine error
!h: grid spacing for SSOR omega parameter estimation
!
!
!Notation: 
!a*b = sum_{i=1,ni} sum_{j=1,nj} a(i,j)*b(i,j)
!A*v = a(i,j)v(i-1,j) + b(i,j)v(i,j-1) + c(i,j)v(i,j) + a(i,j)v(i+1,j) + b(i,j)v(i,j+1)
!
!
!Main algorithm
!1. r_0 := b - A*u_0
!2. z_0 := M^{-1}*r_0  !M is the preconditioner matrix. In here, [M] = diag([A])
!3. p_0 := z_0
!
!4. from k := 0
!   4.1 alpha_k := (r_k*z_k)/(p_k*A*p_k)
!   4.2 u_{k+1} := u_k + \alpha_k p_k
!   4.3 r_{k+1} := r_k - \alpha_k A*p_k
!   4.4 if r_{k+1} is sufficiently small then exit loop 
!   4.5 z_{k+1} := M^{-1}*r_{k+1}
!   4.6 beta_k := \frac{z_{k+1}*r_{k+1}}{z_k*r_k}
!   4.7 p_{k+1} := z_{k+1} + \beta_k p_k
!end 4
!
!The solution is u_{k+1}
implicit none
integer, intent(in) :: ni,nj
real(8), dimension(1:ni,1:nj), intent(in) :: a,b,c,d,e
real(8), dimension(0:ni+1,0:nj+1), intent(in) :: f
real(8), dimension(0:ni+1,0:nj+1), intent(inout) :: u     !explicit-shape dummy arrays
real(8), intent(in) :: h
real(8), optional, intent(in) :: tol
integer, optional, intent(out) :: number_it

integer :: i,j,nit   !number of iterations
real(8) :: alpha,beta,local_tol,omega,rz,rz_old,pAp !local auxiliary variables

rglob=0.0d0
pglob=0.0d0
zglob=0.0d0
Auglob=0.0d0

if(present(tol)) then
    if(tol.eq.0.0d0) then
        local_tol=get_machine_eps()
    else
        local_tol=tol
    end if
else
    local_tol=CGtol
end if

if(verb.eq.2) write(unit=*,fmt=101)
if(verb.eq.2) write(unit=*,fmt='(a45,es10.3e2)')'Conjugate Gradient working. Residual target: ', local_tol

!open(unit=1,file='residual.plt')

omega=2.0d0/(1.0d0+pi*h)

!A dot u old style-----------------------------------------------------------------------------------------------------------------
Auglob=0.0d0
!inner points
do j=1,nj
    do i=2,ni-1
        Auglob(i,j) = a(i,j)*u(i-1,j) &
    &           + b(i,j)*u(i,j-1) &
    &           + c(i,j)*u(i,j)   &
    &           + d(i,j)*u(i+1,j) &
    &           + e(i,j)*u(i,j+1)
    end do
end do

!boundary points
do j=1,nj
    Auglob(1,j)  = a(1,j)*u(ni,j)  &
    &        + b(1,j)*u(1,j-1) &
    &        + c(1,j)*u(1,j)   &
    &        + d(1,j)*u(2,j) &
    &        + e(1,j)*u(1,j+1)
    Auglob(ni,j) = a(ni,j)*u(ni-1,j) &
    &        + b(ni,j)*u(ni,j-1) &
    &        + c(ni,j)*u(ni,j)   &
    &        + d(ni,j)*u(1,j)    &
    &        + e(ni,j)*u(ni,j+1)
end do
!----------------------------------------------------------------------------------------------------------------------------------


rglob(1:ni,1:nj) = f(1:ni,1:nj)-Auglob(1:ni,1:nj)!---------------->1. r_0 := f - A u_0
!SSOR pre-conditioning----------->2. z_0 := M^{-1} r_0-----------------------------------------------------------------------------
zglob=0.0d0

!Lower SSOR
do j=1,nj
    do i=1,ni
        zglob(i,j) = ((-omega*(a(i,j)*zglob(i-1,j)+b(i,j)*zglob(i,j-1)))+omega*(2.0-omega)*rglob(i,j))/(c(i,j))
    end do
end do

!Upper SSOR
do j=nj,1,-1
    do i=ni,1,-1
        zglob(i,j) = (-omega*(d(i,j)*zglob(i+1,j)+e(i,j)*zglob(i,j+1))+c(i,j)*zglob(i,j))/c(i,j)
    end do
end do
!----------------------------------------------------------------------------------------------------------------------------------

pglob(1:ni,1:nj) = zglob(1:ni,1:nj) !----------------------------->3. p_0 := z_0

!r dot z old style ----------------> r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation---------------------------
rz=0.0d0
do j=1,nj
    do i=1,ni
        rz=rz+rglob(i,j)*zglob(i,j)
    end do
end do
!----------------------------------------------------------------------------------------------------------------------------------
rz_old=rz

if(abs(rz).lt.local_tol) then
    if(present(number_it))number_it=0
    if(verb.eq.2) write(*,*)'Caution!!! CG residual is less than tolerance before the first iteration.'
    if(verb.eq.2) write(*,*)'Resuming main program...'
    if(verb.eq.2) write(unit=*,fmt=102)
    return
end if

nit=1
do

    !A dot u old style-----------------------------------------------------------------------------------------------------------------
    Auglob=0.0d0
    !inner points
    do j=1,nj
        do i=2,ni-1
            Auglob(i,j) = a(i,j)*pglob(i-1,j) &
                &       + b(i,j)*pglob(i,j-1) &
                &       + c(i,j)*pglob(i,j)   &
                &       + d(i,j)*pglob(i+1,j) &
                &       + e(i,j)*pglob(i,j+1)
            end do
        end do
        
        !boundary points
        do j=1,nj
            Auglob(1,j)  = a(1,j)*pglob(ni,j)  &
            &            + b(1,j)*pglob(1,j-1) &
            &            + c(1,j)*pglob(1,j)   &
            &            + d(1,j)*pglob(2,j) &
            &            + e(1,j)*pglob(1,j+1)
            Auglob(ni,j) = a(ni,j)*pglob(ni-1,j) &
            &            + b(ni,j)*pglob(ni,j-1) &
            &            + c(ni,j)*pglob(ni,j)   &
            &            + d(ni,j)*pglob(1,j)    &
            &            + e(ni,j)*pglob(ni,j+1)
        end do
    !----------------------------------------------------------------------------------------------------------------------------------

    !p dot Ap (Au is just the name of the vector. Here it means Ap!!) old style -----> r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation
    pAp=0.0d0
    do j=1,nj
        do i=1,ni
            pAp=pAp+pglob(i,j)*Auglob(i,j)
        end do
    end do
    !----------------------------------------------------------------------------------------------------------------------------------

    alpha = rz/pAp                 !------------------------------------->4.1 \alpha_k := (r_k.z_k)/(p_k.A.p_k)
    u(1:ni,1:nj) = u(1:ni,1:nj) + alpha*pglob(1:ni,1:nj)         !----------->4.2 u_{k+1} := u_k + \alpha_k p_k
    rglob(1:ni,1:nj) = rglob(1:ni,1:nj) - alpha*Auglob(1:ni,1:nj)        !----------->4.3 r_{k+1} := r_k - \alpha_k Ap_k

    nit=nit+1    
    if(abs(rz).lt.local_tol) exit   !--------------------------------------->4.4 if r_{k+1} is sufficiently small then exit loop
    
    !SSOR pre-conditioning----------->4.5 z_{k+1} := M^{-1} r_{k+1}. Observe that, here, r is actually r_{k+1}--------------------------
    zglob=0.0d0

    !Lower SSOR
    do j=1,nj
        do i=1,ni
            zglob(i,j) = ((-omega*(a(i,j)*zglob(i-1,j)+b(i,j)*zglob(i,j-1)))+omega*(2.0-omega)*rglob(i,j))/(c(i,j))
        end do
    end do

    !Upper SSOR
    do j=nj,1,-1
        do i=ni,1,-1
            zglob(i,j) = (-omega*(d(i,j)*zglob(i+1,j)+e(i,j)*zglob(i,j+1))+c(i,j)*zglob(i,j))/c(i,j)
        end do
    end do
    !----------------------------------------------------------------------------------------------------------------------------------
    !r dot z old style ----------------------------------------------------------------------------------------------------------------
    rz=0.0d0
    do j=1,nj
        do i=1,ni
            rz=rz+rglob(i,j)*zglob(i,j)
        end do
    end do
    !----------------------------------------------------------------------------------------------------------------------------------
    beta = rz/rz_old !--------------------------------------------------->4.6 \beta_k := \frac{z_{k+1}^T r_{k+1}}{z_k^T r_k}
    pglob(1:ni,1:nj) = zglob(1:ni,1:nj) + beta*pglob(1:ni,1:nj) !-------------------->4.7 p_{k+1} := z_{k+1} + \beta_k p_k
    rz_old=rz !---------------------------------------------->r_k \cdot z_k: for \alpha_k := (r_k^T z_k)/(p_k^T Ap_k) computation for the next iteration

    if(mod(nit,100).eq.0) then 
        if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
    end if
    
end do

if(present(number_it))number_it=nit
if(verb.eq.2) write(unit=*,fmt=100)'Iteration number: ',nit,'Current residual: ',abs(rz)
if(verb.eq.2) write(unit=*,fmt=102)

100 format(1x,a20,i4,3(a20,es10.3e2,1x))
101 format(/,80('-'))
102 format(80('-'),/)

return
end subroutine PCG55_SSOR_x_per_v3
!==================================================================================================================================

!==================================================================================================================================
subroutine alloc_memory_for_Lsys(nimax,njmax)
implicit none
integer, intent(in) :: nimax, njmax

allocate(rglob(0:nimax+1,0:njmax+1))
allocate(pglob(0:nimax+1,0:njmax+1))
allocate(zglob(0:nimax+1,0:njmax+1))
allocate(Auglob(0:nimax+1,0:njmax+1))

end subroutine alloc_memory_for_Lsys
!==================================================================================================================================

end module LSys