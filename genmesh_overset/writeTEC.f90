!======================================================================
! subroutine to output as Tecplot format                 
!======================================================================
subroutine writeTEC(mesh, fname, ug, pg, Tg, phig, rphig, rflag)

  use class
   
  implicit none
   
  character*(*) :: fname  
  type (mesh_class) :: mesh
   
  real(8), intent(in) :: ug(mesh%NNODE,mesh%NSD), &
                         pg(mesh%NNODE), &
						 Tg(mesh%NNODE), &
                         phig(mesh%NNODE), &
                         rphig(mesh%NNODE)
  integer, intent(in) :: rflag

  real(8) :: lxg(mesh%NSD), ldg(mesh%NSD), lug(mesh%NSD), lugm(mesh%NSD), umag, cs
  real(8) :: lacgm(mesh%NSD)
  integer :: ifil, i, j   
  integer :: idxmap(8)  
   
  idxmap(1) = 1 
  idxmap(2) = 2 
  idxmap(3) = 4
  idxmap(4) = 3
  idxmap(5) = 5
  idxmap(6) = 6
  idxmap(7) = 8
  idxmap(8) = 7   
       
  !==============================
  ! output solution in tec
  !==============================   
  ifil = 99  
  write(*,*) "Write solution: ", fname
  ! if (rflag == 1) write(*,*) "Rotate back to ref.; theta =", theta 
  open(ifil, file=fname, status='replace', form='formatted')
      
  write(ifil,*) 'VARIABLES =  "X" "Y" "Z" "U" "V" "W" &
      "P" "PHI" "RPHI" "T" '
  write(ifil,*) 'ZONE N=', mesh%NNode, ', E=', mesh%NElem

  ! Linear Tetrahedral
  if (mesh%NSHLmax == 4) then
    write(ifil,*) 'DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON'   
  ! Linear Tetrahedral and/or Prism
  else if (mesh%NSHLmax == 6) then
    write(ifil,*) 'DATAPACKING=POINT, ZONETYPE=FEBRICK'    
  ! Linear Hexahedral
  else if (mesh%NSHLmax == 8) then
    write(ifil,*) 'DATAPACKING=POINT, ZONETYPE=FEBRICK'
  else
    write(*,*) 'ERROR: undefined mesh%NSHL =', mesh%NSHLmax
    stop
  end if

  do i = 1, mesh%NNode

    if (rflag == 1) then
      ! call rot_vec(mesh%NSD, theta, mesh%xg(i,:), lxg)
      ! call rot_vec(mesh%NSD, theta,      dg(i,:), ldg)
      ! call rot_vec(mesh%NSD, theta,      ug(i,:), lug)
      ! call rot_vec(mesh%NSD, theta,     ugm(i,:), lugm)

    else
      lxg  = mesh%xg(i,:)
      lug  = ug(i,:)
    end if
    umag = sqrt(sum(lug*lug))
    cs = sqrt(1.4d0*286.62d0*Tg(i))
    
    write(ifil,100) lxg, &      ! xg
                    lug, &            ! ug
                    pg(i), &          ! p
                    phig(i), &        ! phi
                    rphig(i), &        ! rphi
                    Tg(i)             ! T
  end do


  if (mesh%NSHLmax == 4) then

    do i = 1, mesh%NElem
      write(ifil,'(4I8)') (mesh%IEN(i,j), j = 1, mesh%NSHLmax)
    end do

  ! Linear Tetrahedral and/or Prism
  else if (mesh%NSHLmax == 6) then  

    do i = 1, mesh%NElem
      if (mesh%NSHL(i) == 4) then
        write(ifil,'(8I8)') mesh%IEN(i,1), mesh%IEN(i,2), mesh%IEN(i,4), mesh%IEN(i,4), &
                            mesh%IEN(i,3), mesh%IEN(i,3), mesh%IEN(i,3), mesh%IEN(i,3)
      else if (mesh%NSHL(i) == 6) then
        write(ifil,'(8I8)') mesh%IEN(i,1), mesh%IEN(i,2), mesh%IEN(i,3), mesh%IEN(i,3), &
                            mesh%IEN(i,4), mesh%IEN(i,5), mesh%IEN(i,6), mesh%IEN(i,6)
      else
        write(*,*) "ERROR: Undefined mesh%NSHL(i) in writeTEC.f90"
        stop
      end if
    end do  

  else if (mesh%NSHLmax == 8) then

    do i = 1, mesh%NElem
      write(ifil,'(8I8)') (mesh%IEN(i,idxmap(j)), j = 1, mesh%NSHLmax)
    end do

  end if

  close(ifil)

  100 format (10E13.5)

end subroutine writeTEC



!======================================================================
!                 
!======================================================================
subroutine writeBndTEC(mesh, fname, bn, ug, pg, Tg, rflag)

  use class
  implicit none
  
  character*(*) :: fname  
  type (mesh_class) :: mesh  
  integer :: bn
  real(8), intent(in) :: ug(mesh%NNODE,mesh%NSD), &
                         pg(mesh%NNODE), &
						 Tg(mesh%NNODE)
  integer, intent(in) :: rflag
  integer :: gn, meshf, i, j
  integer :: G2L(mesh%NNODE)   
  real(8) :: lxg(mesh%NSD), ldg(mesh%NSD), lug(mesh%NSD), lugm(mesh%NSD) 

  write(*,*) "Writting boundary ", bn," : " ,fname 
  meshf = 99
  open(meshf, file=fname, status='replace', form='formatted')

  write(meshf,*) 'VARIABLES =  "X" "Y" "Z" "U" "V" "W" "P" "Umag" '
  write(meshf,*) 'ZONE N=', mesh%bound(bn)%NNODE, ', E=', mesh%bound(bn)%NFace
  write(meshf,*) 'DATAPACKING=POINT, ZONETYPE=FETRIANGLE' 
  
  do i = 1, mesh%bound(bn)%NNODE

    gn      = mesh%bound(bn)%BNODES(i)
    G2L(gn) = i   

    if (rflag == 1) then
      ! call rot_vec(mesh%NSD, theta, mesh%xg(gn,:), lxg)
      ! call rot_vec(mesh%NSD, theta,      dg(gn,:), ldg)
      ! call rot_vec(mesh%NSD, theta,      ug(gn,:), lug)
      ! call rot_vec(mesh%NSD, theta,     ugm(gn,:), lugm)

      ! lug = lug - lugm

    else
      lxg  = mesh%xg(gn,:)
      lug  =      ug(gn,:)
    end if

    write(meshf,'(9E13.5)') lxg, &              ! xg
                            lug, &                    ! ug
                            pg(gn), &                 ! p
                            sqrt(sum(lug**2d0))        ! |ug|
  end do

  do i = 1, mesh%bound(bn)%NFace
    write(meshf,'(3I8)') (G2L(mesh%bound(bn)%FACE_IEN(i,j)), j = 1, mesh%NSHLB)
  end do

  close(meshf)

end subroutine writeBndTEC



!======================================================================
!                 
!======================================================================
subroutine rot_vec(nsd, theta, sol, rot)
  implicit none
  integer, intent(in)  :: nsd
  real(8), intent(in)  :: theta, sol(nsd)
  real(8), intent(out) :: rot(nsd)

  real(8) :: tmpx, tmpy

  tmpx = sol(1)
  tmpy = sol(2)
  rot(1) = cos(-theta)*tmpx - sin(-theta)*tmpy
  rot(2) = sin(-theta)*tmpx + cos(-theta)*tmpy
  rot(3) = sol(3)
end subroutine rot_vec




!======================================================================
! Output boundary meshes for Tecplot        
!======================================================================
subroutine writeBndMeshTEC(mesh, filename, bb)
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
   
  write(mfid,*) 'VARIABLES =  "X" "Y" "Z"'
  write(mfid,*) 'ZONE N=', mesh%bound(bb)%NNODE, ', E=', mesh%bound(bb)%NFace
  write(mfid,*) 'DATAPACKING=POINT, ZONETYPE=FETRIANGLE' 

  do i = 1, mesh%bound(bb)%NNODE
    gn = mesh%bound(bb)%BNODES(i)
    G2L(gn) = i
    write(mfid,*) (mesh%xg(gn,j), j = 1, mesh%NSD)
  end do

  do i = 1, mesh%bound(bb)%NFace
    write(mfid,'(4I8)') (G2L(mesh%bound(bb)%FACE_IEN(i,j)), j = 1, mesh%NSHLB)
  end do  

  close(mfid)
   
end subroutine writeBndMeshTEC
