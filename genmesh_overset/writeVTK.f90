subroutine writeMeshVTK(mesh,fname)

  use class
  
  implicit none
  
  character*(*) fname  
  type (mesh_class) mesh   
  
  integer ifil,i,j
  
  integer  idxmap (8)  
  
  idxmap(1) = 1 
  idxmap(2) = 2 
  idxmap(3) = 4
  idxmap(4) = 3
  idxmap(5) = 5
  idxmap(6) = 6
  idxmap(7) = 8
  idxmap(8) = 7   
       
  !==============================
  ! output solution in vtk
  !==============================   
  ifil = 99
  write(*,*) "Write mesh: ", fname  
  open(ifil, file=fname, status='replace', form='formatted')
      
  write(ifil,'(a)') '# vtk DataFile Version 3.0'
  write(ifil,'(a)') 'vtk output'
  write(ifil,'(a)') 'ASCII'
  write(ifil,'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(ifil,'(a,1x,I8,1x,a)') 'POINTS ',mesh%NNode, 'double'

  do i = 1, mesh%NNode
    write(ifil,'(3E17.8)') (real(mesh%xg(i,j),4), j = 1, mesh%NSD)    
  end do

  if (mesh%NSHLmax.le.6) then   ! Linear Tetrahedral and/or Prism
    write(ifil,'(a,1x,I8,1x,I8)') 'CELLS ',mesh%NElem, mesh%NElem + sum(mesh%NSHL)
    do i = 1, mesh%NElem
      write(ifil,'(5I8)') mesh%NSHL(i), (mesh%IEN(i,j)-1, j = 1, mesh%NSHL(i))
    end do

    write(ifil,'(a,1x,I8)') 'CELL_TYPES',mesh%NElem
    do i = 1, mesh%NElem
      if (mesh%NSHL(i).eq.4) then
        write(ifil,'(I8)') 10
      else  
        write(ifil,'(I8)') 13
      endif
    enddo       

  else if (mesh%NSHLmax.eq.8) then  ! Linear Hexahedral
    write(ifil,'(a,1x,I8,1x,I8)') 'CELLS ',mesh%NElem, mesh%NElem*(mesh%NSHLmax+1)
    do i = 1, mesh%NElem
      write(ifil,'(9I8)') mesh%NSHLmax, (mesh%IEN(i,idxmap(j))-1, j = 1, mesh%NSHLmax)
    end do

    write(ifil,'(a,1x,I8)') 'CELL_TYPES',mesh%NElem
    do i = 1, mesh%NElem
      write(ifil,'(I8)') 12
    enddo   
  else if (mesh%NSHLmax.eq.20) then  ! Quadratic Hexahedral
    write(ifil,'(a,1x,I8,1x,I8)') 'CELLS ',mesh%NElem, mesh%NElem*(mesh%NSHLmax+1)
    do i = 1, mesh%NElem
      write(ifil,'(21I8)') mesh%NSHLmax, (mesh%IEN(i,j)-1, j = 1, mesh%NSHLmax)
    end do

    write(ifil,'(a,1x,I8)') 'CELL_TYPES',mesh%NElem
    do i = 1, mesh%NElem
      write(ifil,'(I8)') 25
    enddo     
  endif 
    
  close(ifil)

end subroutine

!------------------------------------------------------------------------
!                 
!------------------------------------------------------------------------
subroutine writeVTK(mesh,fname,dg,ug,pg,phig)

   use class
   
   implicit none
   
   character*(*) fname  
   type (mesh_class) mesh
   
   real(8) dg(mesh%NNODE,mesh%NSD)
   real(8) ug(mesh%NNODE,mesh%NSD)
   real(8) pg(mesh%NNODE)
   real(8) phig(mesh%NNODE)
   
   integer ifil,i,j
   
   integer  idxmap (8)  
   
   idxmap(1) = 1 
   idxmap(2) = 2 
   idxmap(3) = 4
   idxmap(4) = 3
   idxmap(5) = 5
   idxmap(6) = 6
   idxmap(7) = 8
   idxmap(8) = 7   
       
  !==============================
  ! output solution in vtk
  !==============================   
  ifil = 99  
  write(*,*) "Write boundary: ", fname  
  open(ifil, file=fname, status='replace', form='formatted')
      
  write(ifil,'(a)') '# vtk DataFile Version 3.0'
  write(ifil,'(a)') 'vtk output'
  write(ifil,'(a)') 'ASCII'
  write(ifil,'(a)') 'DATASET UNSTRUCTURED_GRID'

  write(ifil,'(a,1x,I8,1x,a)') 'POINTS ',mesh%NNode, 'double'

  do i = 1, mesh%NNode
    write(ifil,'(3E17.8)') (real(mesh%xg(i,j) + dg(i,j),4), j = 1, mesh%NSD)    
  end do

  if (mesh%NSHLmax.le.6) then   ! Linear Tetrahedral and/or Prism
    write(ifil,'(a,1x,I8,1x,I8)') 'CELLS ',mesh%NElem, mesh%NElem + sum(mesh%NSHL)
    do i = 1, mesh%NElem
      write(ifil,'(5I8)') mesh%NSHL(i), (mesh%IEN(i,j)-1, j = 1, mesh%NSHL(i))
    end do

    write(ifil,'(a,1x,I8)') 'CELL_TYPES',mesh%NElem
    do i = 1, mesh%NElem
      if (mesh%NSHL(i).eq.4) then
        write(ifil,'(I8)') 10
      else  
        write(ifil,'(I8)') 13
      endif
    enddo       
  else if (mesh%NSHLmax.eq.8) then  ! Linear Hexahedral
    write(ifil,'(a,1x,I8,1x,I8)') 'CELLS ',mesh%NElem, mesh%NElem*(mesh%NSHLmax+1)
    do i = 1, mesh%NElem
      write(ifil,'(9I8)') mesh%NSHLmax, (mesh%IEN(i,idxmap(j))-1, j = 1, mesh%NSHLmax)
    end do

    write(ifil,'(a,1x,I8)') 'CELL_TYPES',mesh%NElem
    do i = 1, mesh%NElem
      write(ifil,'(I8)') 12
    enddo   
  else if (mesh%NSHLmax.eq.20) then  ! Quadratic Hexahedral
    write(ifil,'(a,1x,I8,1x,I8)') 'CELLS ',mesh%NElem, mesh%NElem*(mesh%NSHLmax+1)
    do i = 1, mesh%NElem
      write(ifil,'(21I8)') mesh%NSHLmax, (mesh%IEN(i,j)-1, j = 1, mesh%NSHLmax)
    end do

    write(ifil,'(a,1x,I8)') 'CELL_TYPES',mesh%NElem
    do i = 1, mesh%NElem
      write(ifil,'(I8)') 25
    enddo     
  endif 

  write(ifil,'(a,1x,I8)') 'POINT_DATA', mesh%NNode
     
  write(ifil,'(a)') 'VECTORS u double'	  
  do i = 1, mesh%NNode
    write(ifil,'(3E17.8)') real(ug(i,1),4),real(ug(i,2),4),real(ug(i,3),4)
  enddo

  write(ifil,'(a)') 'SCALARS p double'
  write(ifil,'(a)') 'LOOKUP_TABLE default'
  do i = 1, mesh%NNode
    write(ifil,'(E17.8)') real(pg(i),4)
  enddo 

  write(ifil,'(a)') 'SCALARS phi double'
  write(ifil,'(a)') 'LOOKUP_TABLE default'
  do i = 1, mesh%NNode
    write(ifil,'(E17.8)') real(phig(i),4)
  enddo 
    
  close(ifil)

end subroutine
!------------------------------------------------------------------------
!                 
!------------------------------------------------------------------------
subroutine writeBndMeshVTK(mesh,fname,bn)

  use class
      
  implicit none
  
  character*(*) fname  
  type (mesh_class) mesh  
  integer bn

  integer gn,meshf,i,j
  integer G2L(mesh%NNODE)    

  write(*,*) "Writting boundary ", bn," : " ,fname 
  meshf = 99
  open(meshf, file=fname, status='unknown', form='formatted')
      
  write(meshf,'(a)') '# vtk DataFile Version 3.0'
  write(meshf,'(a)') 'vtk output'
  write(meshf,'(a)') 'ASCII'
  write(meshf,'(a)') 'DATASET UNSTRUCTURED_GRID'
  
  write(meshf,'(a,1x,I8,1x,a)') 'POINTS ',mesh%bound(bn)%NNODE, 'double'

  do i = 1, mesh%bound(bn)%NNODE
    gn      = mesh%bound(bn)%BNODES(i)
    G2L(gn) = i     
    write(meshf,'(3E17.8)') (real(mesh%xg(gn,j)), j = 1, mesh%NSD)    
  end do

  write(meshf,'(a,1x,I8,1x,I8)') 'CELLS ',mesh%bound(bn)%NFace, mesh%bound(bn)%NFace*(mesh%NSHLB+1)
  do i = 1,mesh%bound(bn)%NFace
    write(meshf,'(5I8)') mesh%NSHLB, (G2L(mesh%bound(bn)%FACE_IEN(i,j))-1, j = 1, mesh%NSHLB)
  end do

  if (mesh%NSHLB.eq.3) then  ! Linear triangle
    write(meshf,'(a,1x,I8)') 'CELL_TYPES',mesh%bound(bn)%NFace
    do i = 1, mesh%bound(bn)%NFace
      write(meshf,'(I8)') 5
    enddo
  else if (mesh%NSHLB.eq.4) then  ! Linear quad
    write(meshf,'(a,1x,I8)') 'CELL_TYPES',mesh%bound(bn)%NFace
    do i = 1, mesh%bound(bn)%NFace
      write(meshf,'(I8)') 9
    enddo
  endif  
    
  close(meshf)

end subroutine writeBndMeshVTK

!------------------------------------------------------------------------
!                 
!------------------------------------------------------------------------
subroutine writeBndVTK(mesh,fname,bn,dg)

  use class
      
  implicit none
  
  character*(*) fname  
  type (mesh_class) mesh  
  integer bn
  real(8) dg(mesh%NNODE, mesh%NSD)
  
  integer gn,meshf,i,j
  integer G2L(mesh%NNODE)    

  write(*,*) "Writting boundary ", bn," : " ,fname 
  meshf = 99
  open(meshf, file=fname, status='unknown', form='formatted')
      
  write(meshf,'(a)') '# vtk DataFile Version 3.0'
  write(meshf,'(a)') 'vtk output'
  write(meshf,'(a)') 'ASCII'
  write(meshf,'(a)') 'DATASET UNSTRUCTURED_GRID'
  
  write(meshf,'(a,1x,I8,1x,a)') 'POINTS ',mesh%bound(bn)%NNODE, 'double'

  do i = 1, mesh%bound(bn)%NNODE
    gn      = mesh%bound(bn)%BNODES(i)
    G2L(gn) = i     
    write(meshf,'(3E17.8)') (real(mesh%xg(gn,j) + dg(gn,j)), j = 1, mesh%NSD)    
  end do

  write(meshf,'(a,1x,I8,1x,I8)') 'CELLS ',mesh%bound(bn)%NFace, mesh%bound(bn)%NFace*4
  do i = 1,mesh%bound(bn)%NFace
    write(meshf,'(5I8)') 3, (G2L(mesh%bound(bn)%FACE_IEN(i,j))-1, j = 1, mesh%NSHLB)
  end do

  write(meshf,'(a,1x,I8)') 'CELL_TYPES',mesh%bound(bn)%NFace
  do i = 1, mesh%bound(bn)%NFace
    write(meshf,'(I8)') 5
  enddo



  close(meshf)

end subroutine writeBndVTK


