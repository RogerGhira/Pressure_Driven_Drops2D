!Main NSProj_2D printing subroutine for
! 1. Printing field plots;
! 2.
subroutine print_results_to_tecplot(filename)
implicit none
character(len=*), intent(in) :: filename

call print_2D_flow_field(filename)

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
!                write(unit=1,fmt=102)x(i,j),y(i,j),u(i,j),v(i,j),press(i,j),divu(i,j)                
                write(unit=1,fmt=102)x(i,j),y(i,j),0.5d0*(u(i,j)+u(i,j+1)),0.5d0*(v(i,j)+v(i+1,j)),&
            &                       0.25d0*(press(i,j)+press(i+1,j)+press(i,j+1)+press(i+1,j+1)),&                
            &                       0.25d0*(divu(i,j)+divu(i+1,j)+divu(i,j+1)+divu(i+1,j+1))                
        end do
    end do
    100 format('variables = "x","y","u","v","p","div"')    
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
subroutine print_2D_u_field(filename)
implicit none

character(len=*), intent(in) :: filename
integer :: iostatus  !I/O status: 0 for sucessful
integer :: i,j !loop counters
real :: xlocal, ylocal

!Reading data from file++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
open(unit=1, file=filename, status='replace', iostat = iostatus)

fileopen: if(iostatus.eq.0) then
    write(unit=1,fmt=100)
    write(unit=1,fmt=101)time,niu,nju
    do j=1,nju
        do i=1,niu
                xlocal=x(0,j)
                ylocal=0.5d0*(y(0,j-1)+y(0,j))
                write(unit=1,fmt=102)xlocal, ylocal, u(i,j)
        end do
    end do
    100 format('variables = "x","y","u"')    
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

end subroutine print_2D_u_field
!==================================================================================================================================

end subroutine print_results_to_tecplot