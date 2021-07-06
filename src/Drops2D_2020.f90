program Drops2D_2020
!In this version
!1-Crank Nicholson version of the projection version
!2-Pressure free velocity equations algorithm 
!Tentative velocity equation
!u*-1/(2Re)L(u*) = u^{n}+dt*( -[ugrad(u)]^{n+1/2} + 1/(2Re)L(u^{n}) ) !obs: L(u) means Laplacian of u
!u* \cdot \hat{n} = u_wall \cdot \hat{n}
!u* \cdot \hat{t} = u_wall \cdot \hat{t}+dt*grad(chi)^{n+1}  !obs:grad(chi)^{n+1} must be extrapolated
!
!False pressure equation
!L(chi)^{n+1} = div(u*)/dt
!grad(chi^{n+1})\cdot \hat{n} = 0
!
!Projection - update equations
!u^{n+1} = u* - dt grad(chi^{n+1})
!p^{n+1/2} = chi^{n+1} - 1/(2Re) div(u*)
!
!Main reference: Brown, Cortez and Minion, 2001, JCP, 168, pp.464-499
use Level_Set_2D_mod
use Navier_Stokes_2D_mod
use Tools
use LSys
implicit none

integer :: it
character(len=filn_len) :: filn,last_filn='scratch'
integer, dimension(8) :: tini,tfin


call NS_starter()  !memory allocation and variable initialization, including when it is necessary to read it from backup
call fill_abc() !assembling of stencil shaped linear system - filling up of "a", "b", "c", "d", and "e" grid vectors
call BC() !boundary condition application on "a", "b", "c", "d", and "e"

call LS_starter()

!This block can be supressed without significative differences
if(first_it.eq.0) then
    call primer_for_first_iteration() !iterative procedure to perform the first iteration (like in Brown, Cortez and Minion, JCP, 2001, pp. 484)
    time=time+dt
    first_it=first_it+1
end if

call ghosts_update() 

do it=first_it,nit

    !interation number display
    if(verb.ne.0) then
        write(unit=*,fmt=100)'Time step ',it,' of ',nit,' started!'
    end if

       
    call DATE_AND_TIME(values=tini)
    call fill_f_for_uv() !right side of velocity linear systems update
    call solve_ustar_system() !velocity linear system solution
    
    call ghosts_update()
    call fill_f_for_chi() !right side of pressure linear systems update
    call solve_chi_system() !false pressure (chi) linear system solution 

    !velocity and pressure update - including points outside the mesh
    call flow_update() !inside this subroutine ghost_update is called!
    call ugradu_update()
    
    time=time+dt
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!BACKUP SYSTEM AND PLOTTING SUBROUTINES ONLY. NO NUMERICAL ALGORITHMS FROM THIS POINT AHEAD.!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !backup system
    if(mod(it,bck_sv_frq).eq.0) then
        filn = make_file_name('Drops2D_2020','Re',Re,'it',it,'bck')
        if(single_bck_file) then
            call delete_file(last_filn)
            call print_backup(filn)
            last_filn=filn
        else
            call print_backup(filn)
        end if
    end if

    !tecplot printing
    if(mod(it,plt_sv_frq).eq.0.or.it.eq.nit) then
        filn = make_file_name('Drops2D_2020','Re',Re,'it',it,'plt')
        call print_results_to_tecplot_LS(filn)
    end if
    
    !pause system
    call pause_system()
    

    call DATE_AND_TIME(values=tfin)
    if(disp_exec_time) write(unit=*,fmt=1000) 'Total iteration elapsed time :',time_difference(tfin,tini)
1000 format(a60,f8.3)

    !interation number display
    if(verb.ne.0) then
        write(unit=*,fmt=101)'Time step ',it,' completed!'
    end if

end do !main loop end
100 format(/,22('+'),a10,i6,a4,i6,a10,22('+'))
101 format(22('+'),a10,i6,a10,32('+'),/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program Drops2D_2020