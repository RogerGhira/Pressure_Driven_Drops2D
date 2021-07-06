subroutine print_results_to_tecplot_LS(filename)
implicit none
character(len=*), intent(in) :: filename

call print_2D_flow_field(filename)
!call print_2D_flow_field_cell_center(filename)

contains
!==================================================================================================================================
subroutine print_2D_flow_field(filename)
implicit none

character(len=*), intent(in) :: filename
integer :: iostatus  !I/O status: 0 for sucessful
integer :: i,j !loop counters

!Reading data from file++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
open(unit=1, file=filename, status='replace', iostat = iostatus)

fileopen: if(iostatus.eq.0) then
    write(unit=1,fmt=100)
    write(unit=1,fmt=101)time,nx+1,ny+1
    do j=0,ny
        do i=0,nx
                write(unit=1,fmt=102)x(i,j),y(i,j),0.5d0*(u(i,j)+u(i,j+1)),0.5d0*(v(i,j)+v(i+1,j)),&
            &                       0.25d0*(press(i,j)+press(i+1,j)+press(i,j+1)+press(i+1,j+1)),  &                 
            &                       0.25d0*(divu(i,j)+divu(i+1,j)+divu(i,j+1)+divu(i+1,j+1)),      &                
            &                       0.25d0*(phi(i,j,1)+phi(i+1,j,1)+phi(i,j+1,1)+phi(i+1,j+1,1)), &
            &                       0.25d0*(kappa(i,j)+kappa(i+1,j)+kappa(i,j+1)+kappa(i+1,j+1)), &
            &                       0.25d0*(normalx(i,j)+normalx(i+1,j)+normalx(i,j+1)+normalx(i+1,j+1)), &
            &                       0.25d0*(normaly(i,j)+normaly(i+1,j)+normaly(i,j+1)+normaly(i+1,j+1)), &
            &                       0.25d0*(mask(i,j)+mask(i+1,j)+mask(i,j+1)+mask(i+1,j+1))
        end do
    end do
    100 format('variables = "x","y","u","v","p","div","<greek>f</greek>", "<greek>k</greek>", &
        &      "n<sub>x</sub>", "n<sub>y</sub>", "mask"')    
    101 format('zone T="',f6.2,'", i=',i6,', j=',i6,', f=point')    
    102 format(1x,20(e18.10))    
    close(unit=1,status='keep')
else fileopen
    !Impossible to open file. Tell user.
    write(unit=*,fmt=9001) filename,iostatus
    stop
    9001 format('Imposible to open "',a30,'". Status=',i6,/,'INTERRUPTED EXECUTION!')
end if fileopen
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine print_2D_flow_field
!==================================================================================================================================

!==================================================================================================================================
subroutine print_2D_flow_field_cell_center(filename)
implicit none

character(len=*), intent(in) :: filename
integer :: iostatus  !I/O status: 0 for sucessful
integer :: i,j !loop counters

!Reading data from file++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
open(unit=1, file=filename, status='replace', iostat = iostatus)

fileopen: if(iostatus.eq.0) then
    write(unit=1,fmt=100)
    write(unit=1,fmt=101)time,nx,ny
    do j=1,ny
        do i=1,nx
                write(unit=1,fmt=102)0.5d0*(x(i,j)+x(i-1,j)),0.5d0*(y(i,j)+y(i,j-1)),0.5d0*(u(i-1,j)+u(i,j)),0.5d0*(v(i,j-1)+v(i,j)),&
            &                       press(i,j),divu(i,j),phi(i,j,1),kappa(i,j), normalx(i,j), normaly(i,j), 1.0d0*mask(i,j)
        end do
    end do
    100 format('variables = "x","y","u","v","p","div","<greek>f</greek>", "<greek>k</greek>", &
        &      "n<sub>x</sub>", "n<sub>y</sub>", "mask"')    
    101 format('zone T="',f6.2,'", i=',i6,', j=',i6,', f=point')    
    102 format(1x,20(e18.10))    
    close(unit=1,status='keep')
else fileopen
    !Impossible to open file. Tell user.
    write(unit=*,fmt=9001) filename,iostatus
    stop
    9001 format('Imposible to open "',a30,'". Status=',i6,/,'INTERRUPTED EXECUTION!')
end if fileopen
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end subroutine print_2D_flow_field_cell_center
!==================================================================================================================================


end subroutine print_results_to_tecplot_LS