!Main NSProj_2D boundary condition subroutine.
subroutine fill_abc()
implicit none

!u and v components at same time
au(:,:)=-dt/(2*Re*dx**2)
bu(:,:)=-dt/(2*Re*dy**2)
cu(:,:)=1.0d0+dt/(Re*dx**2)+dt/(Re*dy**2)
du(:,:)=-dt/(2*Re*dy**2)
eu(:,:)=-dt/(2*Re*dy**2)

av(:,:)=-dt/(2*Re*dx**2)
bv(:,:)=-dt/(2*Re*dy**2)
cv(:,:)=1.0d0+dt/(Re*dx**2)+dt/(Re*dy**2)
dv(:,:)=-dt/(2*Re*dy**2)
ev(:,:)=-dt/(2*Re*dy**2)

!chi component
achi(:,:)=1.0d0/(dx**2)
bchi(:,:)=1.0d0/(dy**2)
cchi(:,:)=-2.0d0/(dx**2)-2.0d0/(dy**2)
dchi(:,:)=1.0d0/(dx**2)
echi(:,:)=1.0d0/(dx**2)

return
end subroutine fill_abc