module Tools
implicit none

integer, parameter :: filn_len=60

contains
!==================================================================================================================================
pure function make_file_name(base,real_param_name,real_param,int_param_name,int_param,term)
implicit none

character(len=*), intent(in) :: base  !base
character(len=*), intent(in) :: term   !termination
character(len=*), intent(in) :: real_param_name   !termination
character(len=*), intent(in) :: int_param_name   !termination
real, intent(in) :: real_param
integer, intent(in) :: int_param

character(len=filn_len) :: make_file_name
character(len=filn_len) :: aux
integer :: i,j

write(unit=aux,fmt=001)base,real_param_name,real_param,int_param_name,int_param,term

001 format(a24,'_',a6,'=',es9.3e1,'_',a6,'=',i7,'.',a3)

make_file_name=' '

!This code puts all blank spaces to the left
j=1
do i=1,filn_len
    if(aux(i:i).ne.' ')then
        make_file_name(j:j)=aux(i:i)
        j=j+1
    end if
end do

!bad characters replacement
do i=1,filn_len
    if(make_file_name(i:i).eq.'+') make_file_name(i:i)='p'
    if(make_file_name(i:i).eq.'-') make_file_name(i:i)='-'
end do

return
end function make_file_name
!==================================================================================================================================

!==================================================================================================================================
pure function time_difference(stop_crono,start_crono)
implicit none

integer, dimension(8), intent(in) :: stop_crono,start_crono
real(8) :: time_difference

integer :: delta_day
real(8) :: t_ini,t_fin

delta_day=stop_crono(3)-start_crono(3)

t_ini=real(start_crono(5)*60*60*1000+start_crono(6)*60*1000+start_crono(7)*1000+start_crono(8))
if (delta_day.eq.1) then
   t_fin=real((stop_crono(5)+24)*60*60*1000+stop_crono(6)*60*1000+stop_crono(7)*1000+stop_crono(8))
else
   t_fin=real(stop_crono(5)*60*60*1000+stop_crono(6)*60*1000+stop_crono(7)*1000+stop_crono(8))
end if

time_difference=(t_fin-t_ini)/1.0d3

return
end function time_difference
!==================================================================================================================================

!==================================================================================================================================
subroutine pause_system(probe_file)
implicit none
integer, intent(in), optional :: probe_file !this variable is necessary to properlly close probe file if execution need to be interrupted.
integer :: iostatus
character(len=1) :: kill

open(unit=1, file='pause', status='old', iostat = iostatus)

if(iostatus.eq.0) then
    write(*,*)'Simulation paused. To kill enter "K", to resume, anything else. :) '
    read(*,*)kill
    close(unit=1,status='delete')
    if(kill.eq.'k'.or.kill.eq.'K') then
        if(present(probe_file)) close(unit=probe_file,status='keep')  !closing probe file before interrupt program...
        stop
    end if
    return
end if

end subroutine pause_system
!==================================================================================================================================

!==================================================================================================================================
subroutine delete_file(filename)
implicit none
character(len=*), intent(in) :: filename

open(unit=1,file=filename)
close(unit=1,status='delete')

end subroutine delete_file
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
end module Tools
