subroutine writeIEN(mesh ,filename)     

  use class

  implicit none

  type (mesh_class) mesh   
  character*(*) filename
  integer meshf,i,j
      
  meshf = 12     
  write(*,*) "Write graph: ", filename
  open(meshf, file = filename, status = 'unknown')
      
  ! write number of elements and element type
   write(meshf,*) mesh%NElem 
!  if (mesh%NSHLmax.eq.4) then  ! Tetrahedral
!    write(meshf,*) mesh%NElem, 2
!  else if (mesh%NSHL.eq.6) then
!    write(meshf,*) mesh
!  else if (mesh%NSHL.eq.8) then  ! Hexahedral
!    write(meshf,*) mesh%NElem, 3
!  else   
!    write(*,*) "unsupported IEN output"
!    write(*,*) "NSHL: ",mesh%NSHLmax
!    stop
!  endif
   
   
        
  ! write nodal numbering corresponding to the element
  do i = 1, mesh%NElem
    write(meshf,'(9I8)') (mesh%IEN(i,j), j = 1, mesh%NSHL(i))
  end do
      
  close(meshf)
        
end subroutine writeIEN




subroutine writeE2EGraph(mesh ,filename)

 use class

 implicit none

 type (mesh_class) mesh
 character*(*) filename
 integer gfid,i,j,e1,e2,n

 integer :: cnt1(mesh%NNODE),cnt2(mesh%NELEM)
 integer :: flag(mesh%NELEM)

 integer :: INE(mesh%NNODE,99)  ! 99 is arbitrarily big....
 integer :: E2E(mesh%NELEM,4)   ! 4 is based on flag treshold of 4 ==> E2E based on faces
!
!... Invert relevant part of IEN to INE to get node to element map
!

! cnt1 = 0
! do e1 = 1, mesh%NELEM
!   do i = 1, mesh%NSHL
!     n = mesh%IEN(e1,i)
!     cnt1(n) = cnt1(n)+ 1
!     INE(n,cnt1(n)) = e1
!   enddo
! enddo
! write(*,*) "  Max node 2 element connectivity", maxval(cnt1)

!
!... Use INE to generate the element to element (E2E) connectivity graph
!
! cnt2 = 0
! flag = 0
! do e1 = 1, mesh%NELEM

!   do i = 1, mesh%NSHL
!     n = mesh%IEN(e1,i)
!     do j = 1, cnt1(n)
!       e2 = INE(n,j)
!       flag(e2) = flag(e2) + 1
!     enddo
!   enddo

!   do i = 1, mesh%NSHL
!     n = mesh%IEN(e1,i)
!     do j = 1, cnt1(n)
!       e2 = INE(n,j)
!       if ((flag(e2).ge.3).and.(e1.ne.e2)) then
!         cnt2(e1) = cnt2(e1) + 1
!         E2E(e1,cnt2(e1)) = e2
!       endif
!       flag(e2)=0
!     enddo
!   enddo

! enddo

 write(*,*) "  Max element 2 element connectivity", maxval(cnt2)

!
!... write graph
!
 gfid = 12
 write(*,*) "Write graph: ", filename
 open(gfid, file = filename, status = 'unknown')
 write (gfid, *) mesh%NELEM, sum(cnt2)/2
 do e1 = 1, mesh%NELEM
   write(gfid,'(256I8)') (E2E(e1,j),j=1,cnt2(e1))
 enddo
 close(gfid)

end subroutine writeE2EGraph
