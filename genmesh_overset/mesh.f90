!======================================================================
!          
!======================================================================
subroutine deallocMesh(mesh)
 use class 
 implicit none
 type (mesh_class)   mesh
 integer :: i
   
 deallocate(mesh%xg)
 deallocate(mesh%IEN)

 do i = 1, mesh%NBOUND
   deallocate(mesh%bound(i)%FACE_IEN)
   deallocate(mesh%bound(i)%F2E)
   deallocate(mesh%bound(i)%FACE_OR)
   deallocate(mesh%bound(i)%BNODES)
 end do 
 deallocate(mesh%bound)
  
 do i = 1, mesh%NPATCH 
   deallocate(mesh%patch(i)%U_KNOT)
   deallocate(mesh%patch(i)%V_KNOT)
   deallocate(mesh%patch(i)%W_KNOT)
 end do 
 if (mesh%NPATCH.ge.1) then
   deallocate(mesh%patch) 
   deallocate(mesh%wg)
   deallocate(mesh%EPID)
   deallocate(mesh%EIJK)
 endif    
  
end subroutine deallocMesh



!======================================================================
!          
!======================================================================
subroutine readMesh(mesh,filename)
  use class 
  implicit none
  type (mesh_class) :: mesh
  character*(*) :: filename
  integer :: mfid
  integer :: i,j,k
  integer :: ierr

  ! open mesh files
  mfid = 11      
  write(*,*) "Read mesh: ", filename   
  open (mfid, file=filename, status='old', iostat=ierr)
  if (ierr.ne.0) then
    write(*,*) "Can not open file: ", filename
    stop   
  end if
      
!  read(mfid,*) mesh%NSD,  mesh%NSHL, mesh%NSHLB
  read(mfid,*) mesh%NSD, mesh%NSHLmax, mesh%NSHLB
  read(mfid,*) mesh%NNODE, mesh%NELEM, mesh%NBOUND, mesh%NPATCH

  ! read nodes
  allocate(mesh%xg(mesh%NNODE,mesh%NSD), mesh%NID(mesh%NNODE))
  do i = 1, mesh%NNODE
    read(mfid,*) (mesh%xg(i,j), j = 1,mesh%NSD), mesh%NID(i)
  end do  

  ! read elements  
  allocate(mesh%NSHL(mesh%NELEM), mesh%EID(mesh%NELEM), &
           mesh%IEN(mesh%NELEM,mesh%NSHLmax))
  do i = 1, mesh%NELEM
    read(mfid,*) mesh%NSHL(i), (mesh%IEN(i,j), j = 1, mesh%NSHL(i)), &
                 mesh%EID(i)
  end do

  ! read faces
  allocate(mesh%bound(mesh%NBOUND)) 
  do i = 1, mesh%NBOUND
      
    read(mfid,*) mesh%bound(i)%FACE_ID,mesh%bound(i)%NFACE,mesh%bound(i)%NNODE
      
    allocate(mesh%bound(i)%FACE_IEN(mesh%bound(i)%NFACE,mesh%NSHLB))  
    allocate(mesh%bound(i)%F2E(mesh%bound(i)%NFACE))
    allocate(mesh%bound(i)%FACE_OR(mesh%bound(i)%NFACE))     
      
    do j = 1,mesh%bound(i)%NFACE
      read(mfid,*) (mesh%bound(i)%FACE_IEN(j,k), k = 1,mesh%NSHLB)
    end do     
    do j = 1,mesh%bound(i)%NFACE
      read(mfid,*) mesh%bound(i)%F2E(j),mesh%bound(i)%FACE_OR(j)
    end do   
         
    allocate(mesh%bound(i)%BNODES(mesh%bound(i)%NNODE))   
    do j = 1,mesh%bound(i)%NNODE
      read(mfid,*) mesh%bound(i)%BNODES(j)
    end do        
  end do


  ! read nurbs data
  if (mesh%NPATCH.gt.0) then
    allocate(mesh%wg(mesh%NNODE))
    do i = 1, mesh%NNODE
      read(mfid,*) mesh%wg(i)
    end do   
         
    allocate(mesh%EPID(mesh%NELEM), mesh%EIJK(mesh%NELEM,mesh%NSD))
   
    do i = 1, mesh%NELEM
      read(mfid,*) mesh%EPID(i),(mesh%EIJK(i,j), j = 1,mesh%NSD)
    end do  
  end if
   
  ! read patches 
  allocate(mesh%patch(mesh%NPATCH)) 
  do i = 1, mesh%NPATCH

    read(mfid,*) mesh%patch(i)%P,mesh%patch(i)%Q,mesh%patch(i)%R        
    read(mfid,*) mesh%patch(i)%MCP,mesh%patch(i)%NCP,mesh%patch(i)%OCP
        
    allocate(mesh%patch(i)%U_KNOT(mesh%patch(i)%MCP+mesh%patch(i)%P+1))
    allocate(mesh%patch(i)%V_KNOT(mesh%patch(i)%NCP+mesh%patch(i)%Q+1))
    allocate(mesh%patch(i)%W_KNOT(mesh%patch(i)%OCP+mesh%patch(i)%R+1))

    read(mfid,*) (mesh%patch(i)%U_KNOT(j), j=1,mesh%patch(i)%MCP+mesh%patch(i)%P+1)        
    read(mfid,*) (mesh%patch(i)%V_KNOT(j), j=1,mesh%patch(i)%NCP+mesh%patch(i)%Q+1)    
    read(mfid,*) (mesh%patch(i)%W_KNOT(j), j=1,mesh%patch(i)%OCP+mesh%patch(i)%R+1)

  end do      
  
  close(mfid)
    
end subroutine readMesh



!======================================================================
!          
!======================================================================
subroutine writeMesh(mesh, filename)
  use class 
  implicit none
  type(mesh_class) :: mesh
  character*(*) filename
  integer :: mfid
  integer :: i, j, k

  ! open mesh files
  mfid = 11    
  write(*,*) "Write mesh: ", filename      
  open(mfid, file=filename, status='replace')   
  write(mfid,*) mesh%NSD,   mesh%NSHLmax, mesh%NSHLB
  write(mfid,*) mesh%NNODE, mesh%NELEM, mesh%NBOUND, mesh%NPATCH

  ! write nodes
  do i = 1, mesh%NNODE
    write(mfid,*) (mesh%xg(i,j), j = 1, mesh%NSD), mesh%NID(i)
  end do  
    
  ! write elements  
  do i = 1, mesh%NELEM
    write(mfid,'(256I12)') mesh%NSHL(i), &
                           (mesh%IEN(i,j), j = 1, mesh%NSHL(i)), &
                           mesh%EID(i)
  end do
   
  ! write faces
  do i = 1, mesh%NBOUND
    write(mfid,*) mesh%bound(i)%FACE_ID,mesh%bound(i)%NFACE,mesh%bound(i)%NNODE
    do j = 1,mesh%bound(i)%NFACE
      write(mfid,'(256I12)') (mesh%bound(i)%FACE_IEN(j,k), k = 1,mesh%NSHLB)
    end do            
    do j = 1,mesh%bound(i)%NFACE
      write(mfid,*) mesh%bound(i)%F2E(j), mesh%bound(i)%FACE_OR(j)
    end do   
    do j = 1,mesh%bound(i)%NNODE
      write(mfid,*) mesh%bound(i)%BNODES(j)
    end do  
  end do

  ! write nurbs data
  if (mesh%NPATCH.gt.0) then
    do i = 1, mesh%NNODE
      write(mfid,*) mesh%wg(i)
    end do   
         
    do i = 1, mesh%NELEM
      write(mfid,*) mesh%EPID(i),(mesh%EIJK(i,j), j = 1,mesh%NSD)
    end do  
  end if
   
  ! write patches 
  do i = 1, mesh%NPATCH
    write(mfid,*) mesh%patch(i)%P,mesh%patch(i)%Q,mesh%patch(i)%R
    write(mfid,*) mesh%patch(i)%MCP,mesh%patch(i)%NCP,mesh%patch(i)%OCP
        
    write(mfid,'(256E21.12)') (mesh%patch(i)%U_KNOT(j), j=1,mesh%patch(i)%MCP+mesh%patch(i)%P+1)
    write(mfid,'(256E21.12)') (mesh%patch(i)%V_KNOT(j), j=1,mesh%patch(i)%NCP+mesh%patch(i)%Q+1)     
    write(mfid,'(256E21.12)') (mesh%patch(i)%W_KNOT(j), j=1,mesh%patch(i)%OCP+mesh%patch(i)%R+1)
  end do      
      
  close(mfid)

end subroutine writeMesh



!======================================================================
! Output mesh for blade (or boundary)        
!======================================================================
subroutine writeBndMesh(mesh, filename, bb)
  use class 
  implicit none
  integer, intent(in) :: bb
  type (mesh_class) :: mesh
  character(len=*) :: filename
  integer :: mfid
  integer :: i, j, gn, G2L(mesh%NNODE)

  G2L = -1

  ! open mesh files
  mfid = 11    
  write(*,*) "Write boundary mesh: ", filename      
  open(mfid, file=filename, status='replace')   
  write(mfid,*) mesh%NSD, mesh%NSHLB
   
!  ! write faces
!  do bb = 1, mesh%NBOUND
      
!    ! only output the blade
!    if (mesh%bound(bb)%FACE_ID == 7)  then

      write(mfid,*) mesh%bound(bb)%NNODE, mesh%bound(bb)%NFace, &
                    bb, mesh%bound(bb)%FACE_ID

      do i = 1, mesh%bound(bb)%NNODE
        gn = mesh%bound(bb)%BNODES(i)
        G2L(gn) = i
        write(mfid,*) (mesh%xg(gn,j), j = 1, mesh%NSD)
      end do

      do i = 1, mesh%bound(bb)%NFace
        write(mfid,'(4I8)') (G2L(mesh%bound(bb)%FACE_IEN(i,j)), j = 1, mesh%NSHLB)
      end do  

!    end if
!  end do   

  close(mfid)
   
end subroutine writeBndMesh
