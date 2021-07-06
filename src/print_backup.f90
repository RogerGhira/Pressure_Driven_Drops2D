subroutine print_backup(filename)
implicit none
character(len=*), intent(in) :: filename

integer :: iostatus  !I/O status: 0 for sucessful

!Reading data from file++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
open(unit=1, file=filename, status='replace', iostat = iostatus, form='unformatted')

fileopen: if(iostatus.eq.0) then
    write(unit=1)u
    write(unit=1)v
    write(unit=1)press
    close(unit=1)
else fileopen
    !Impossible to open file. Tell user.
    write(unit=*,fmt=9001) filename,iostatus
    stop
    9001 format('Imposible to open "',a60,'". Status=',i6,/,'EXECUTION INTERRUPTED!')
end if fileopen
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


end subroutine print_backup
    