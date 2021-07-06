subroutine primer_for_first_iteration()
use LSys
implicit none

integer :: i,j,number_of_chi_it
real(8) :: dudx,dudy,dvdx,dvdy,ubar,vbar

write(unit=*,fmt=*)'Initialization in progress...'
call ghosts_update()
do
    !convective terms for the initial flow field - using central differences here
    call ugradu_update_for_primer()

    !basically a normal time step
    call fill_f_for_uv()
    call solve_ustar_system()
    
    call ghosts_update()
    call fill_f_for_chi()
    call solve_chi_system(number_of_chi_it)
    
!    call ghosts_update()
    call flow_update()
    time=0.0d0+dt
    
    if(number_of_chi_it.eq.0) then !very hard contition to reach.
        call ugradu_update()
        write(unit=*,fmt=*)'Initialization successful!'
        exit
    end if
    write(unit=*,fmt=*)'Number of iterations for pressure convergence:',number_of_chi_it
    
    !setting flow field back to the initial
    if(initial_kind.eq.1) then
        call set_initial_flow_field()
        u=uinitial
        v=vinitial
        press=pressinitial
    else
        u=0.0d0
        v=0.0d0
    end if

    !convective term update considering that adams-bashforth formula is used in fill_f_for_uv. Doing like in the lines beneath
    !one gets ugradu^{1/2}=(ugradu^{1}+ugradu^{0})/2
    ugradu(:,:,1)=ugradu(:,:,1)/3 
    ugradv(:,:,1)=ugradv(:,:,1)/3
end do

return
end subroutine primer_for_first_iteration