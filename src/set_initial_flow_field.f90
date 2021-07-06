subroutine set_initial_flow_field()
implicit none
integer :: i,j

uinitial(:,:)=0.0d0
vinitial(:,:)=0.0d0
pressinitial(:,:)=0.0d0

return
end subroutine set_initial_flow_field