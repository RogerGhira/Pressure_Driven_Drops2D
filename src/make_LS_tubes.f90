subroutine make_LS_tubes()
implicit none

integer :: i,j,p,q
real :: h

h=max(dx,dy)
n_tube=0
mask=4 !out of the last tube, mask is set to 4.
index_i=0
index_j=0

!inner band construction - bisection method like algorithm. If the zero level set crosses any face of an cell, the mask is set to -1.
!It's likely useless.
do j=0,njp+1
    do i=0,nip+1
        if((((phi(i,j,1)*phi(i+1,j,1).lt.0.0d0).or.(phi(i+1,j,1)*phi(i+1,j+1,1).lt.0.0d0).or.(phi(i+1,j+1,1)*phi(i,j+1,1).lt.0.0d0).or.(phi(i,j+1,1)*phi(i,j,1).lt.0.0d0)).or.&
        &   ((phi(i,j,1)*phi(i,j+1,1).lt.0.0d0).or.(phi(i,j+1,1)*phi(i-1,j+1,1).lt.0.0d0).or.(phi(i-1,j+1,1)*phi(i-1,j,1).lt.0.0d0).or.(phi(i-1,j,1)*phi(i,j,1).lt.0.0d0)).or.&
        &   ((phi(i,j,1)*phi(i-1,j,1).lt.0.0d0).or.(phi(i-1,j,1)*phi(i-1,j-1,1).lt.0.0d0).or.(phi(i-1,j-1,1)*phi(i,j-1,1).lt.0.0d0).or.(phi(i,j-1,1)*phi(i,j,1).lt.0.0d0)).or.&
        &   ((phi(i,j,1)*phi(i,j-1,1).lt.0.0d0).or.(phi(i,j-1,1)*phi(i+1,j-1,1).lt.0.0d0).or.(phi(i+1,j-1,1)*phi(i+1,j,1).lt.0.0d0).or.(phi(i+1,j,1)*phi(i,j,1).lt.0.0d0))).and.&
        &   mask(i,j).eq.4)then
            mask(i,j)=-1
            n_tube=n_tube+1
            index_i(n_tube)=i
            index_j(n_tube)=j
        end if
    end do
end do


do j=0,njp+1
    do i=0,nip+1

        !alpha tube
        if(abs(phi(i,j,1)).le.alpha.and.mask(i,j).ne.-1) then
            mask(i,j)=0
            n_tube=n_tube+1
            index_i(n_tube)=i
            index_j(n_tube)=j
        end if
        
        !beta tube
        if(abs(phi(i,j,1)).gt.alpha.and.abs(phi(i,j,1)).le.beta) then
            mask(i,j)=1
            n_tube=n_tube+1
            index_i(n_tube)=i
            index_j(n_tube)=j
        end if
        
        !gamma tube
        if(abs(phi(i,j,1)).gt.beta.and.abs(phi(i,j,1)).le.gamma) then
            mask(i,j)=2
            n_tube=n_tube+1
            index_i(n_tube)=i
            index_j(n_tube)=j
            
            !q and p loops find the border. This is really necessary
            do q=j-1,j+1
                do p=i-1,i+1
                    if(mask(p,q).eq.4.and.abs(phi(p,q,1)).gt.gamma) then
                        mask(p,q)=3
                        n_tube=n_tube+1
                        index_i(n_tube)=p
                        index_j(n_tube)=q
                    end if
                end do
            end do
        end if
    end do
end do

!reseting phi for regions out of the last tube. That's not tubing, in fact, but keeps the phi clean thought the domain.
do j=0,njp+1
    do i=0,nip+1

            if(phi(i,j,1).lt.-gamma-0.5d0*(dx+dy)) then
                phi(i,j,1)=-gamma-0.5d0*(dx+dy)
                phi_x(i,j)=0.0d0
                phi_y(i,j)=0.0d0
            end if
            
            if(phi(i,j,1).gt.gamma+0.5d0*(dx+dy)) then
                phi(i,j,1)=gamma+0.5d0*(dx+dy)
                phi_x(i,j)=0.0d0
                phi_y(i,j)=0.0d0
            end if

    end do
end do

return
end subroutine make_LS_tubes
