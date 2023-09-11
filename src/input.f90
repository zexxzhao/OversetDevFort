!======================================================================
! 
!======================================================================
subroutine input(id)
  implicit none
  integer, intent(in) :: id
  call input_mesh(id)
  call input_bmesh()
end subroutine input


!======================================================================
!
!======================================================================
subroutine input_mesh(id)
  use aAdjKeep
  use commonvars
  use mpi
  implicit none

  integer, intent(in) :: id

  character(len=30) :: fname, iname
  character(len=10) :: cname

  integer :: mfid, i, j, k, itmp1, itmp2, itmp3, nshl, counter

  counter = 1
  ! open mesh files
  mfid = 11
  fname = 'part'//trim(cname(id))//'.dat'
  open (mfid, file=fname, status='old')

  read (mfid, *) NSD, NSHL, NSHLBmax
  read (mfid, *) NNODE, NELEM, NBOUND, NPATCH

!  if (myid + 1 == 11) then
  !       write(*,*) counter, NSD, NSHL, NSHLBmax
  !      counter = counter + 1
  !     write(*,*) counter, NNODE, NELEM, NBOUND, NPATCH
  !    counter = counter + 1
!endif

  ! read nodes
  allocate (xg(NNODE, NSD), NodeID(NNODE))
  do i = 1, NNODE
    if (fem_flag == 1) then
      read (mfid, *) (xg(i, j), j=1, NSD), NodeID(i)
      ! if (myid + 1 == 1) then
      !   write (*, *) counter, xg(i, :)
      !   counter = counter + 1
      ! end if
    else
      read (mfid, *) (xg(i, j), j=1, NSD)
    end if
  end do
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

  ! read elements
  allocate (ELMNSHL(NELEM), IEN(NELEM, NSHL))
  allocate (ELM_ID(NELEM))
  do i = 1, NELEM
    if (fem_flag == 1) then
      read (mfid, *) ELMNSHL(i), (IEN(i, j), j=1, ELMNSHL(i)), ELM_ID(i)
    else
      read (mfid, *) ELMNSHL(i), (IEN(i, j), j=1, ELMNSHL(i))
    end if
    ! if ((id == 120) .and. (i == NELEM)) then
    !   write (*, *), ELMNSHL(i), (IEN(i, j), j=1, ELMNSHL(i))
    ! end if
  end do
  maxNSHL = maxval(ELMNSHL)
  ! read faces
  allocate (bound(NBOUND))
  do i = 1, NBOUND

    read (mfid, *) bound(i)%FACE_ID, bound(i)%NFACE, bound(i)%NNODE

    allocate (bound(i)%FACE_IEN(bound(i)%NFACE, NSHLBmax))
    allocate (bound(i)%F2E(bound(i)%NFACE))
    allocate (bound(i)%FACE_OR(bound(i)%NFACE))
    allocate (bound(i)%NSHLB(bound(i)%NFACE))
    do j = 1, bound(i)%NFACE
      bound(i)%NSHLB(j) = NSHLBmax
      if (fem_flag == 1) then
        read (mfid, *) (bound(i)%FACE_IEN(j, k), k=1, NSHLBmax)
      else
        read (mfid, *) bound(i)%NSHLB(j), (bound(i)%FACE_IEN(j, k), k=1, NSHLBmax)
      end if
    end do
    do j = 1, bound(i)%NFACE
      read (mfid, *) bound(i)%F2E(j), bound(i)%FACE_OR(j)
    end do

    allocate (bound(i)%BNODES(bound(i)%NNODE))
    do j = 1, bound(i)%NNODE
      read (mfid, *) bound(i)%BNODES(j)
    end do
  end do
  ! read nurbs data
  allocate (wg(NNODE))

  if (NPATCH > 0) then
    iga = .true.
    do i = 1, NNODE
      read (mfid, *) wg(i)
    end do

    allocate (EPID(NELEM), EIJK(NELEM, NSD))

    do i = 1, NELEM
      read (mfid, *) EPID(i), (EIJK(i, j), j=1, NSD)
    end do
  else
    iga = .false.
    wg = 1.0d0
  end if

  ! read patches
  ! allocate (patch(NPATCH))
  ! do i = 1, NPATCH

  !   read (mfid, *) patch(i)%P, patch(i)%Q, patch(i)%R
  !   read (mfid, *) patch(i)%MCP, patch(i)%NCP, patch(i)%OCP

  !   allocate (patch(i)%U_KNOT(patch(i)%MCP + patch(i)%P + 1))
  !   allocate (patch(i)%V_KNOT(patch(i)%NCP + patch(i)%Q + 1))
  !   allocate (patch(i)%W_KNOT(patch(i)%OCP + patch(i)%R + 1))

  !   read (mfid, *) (patch(i)%U_KNOT(j), j=1, patch(i)%MCP + patch(i)%P + 1)
  !   read (mfid, *) (patch(i)%V_KNOT(j), j=1, patch(i)%NCP + patch(i)%Q + 1)
  !   read (mfid, *) (patch(i)%W_KNOT(j), j=1, patch(i)%OCP + patch(i)%R + 1)
  ! end do

  close (mfid)

  ! Init Flags
  allocate (IPER(NNODE))
  do i = 1, NNODE
    IPER(i) = i
  end do

  allocate (IBC(NNODE, 2*NSD + 2))
  IBC = 0
  ! allocate (IS_SOLID_NODE(NNODE))
  ! IS_SOLID_NODE_ASSIGNED = .false.

  allocate (EL_TYP(NELEM))
  EL_TYP = 0

  allocate (P_Flag(NNODE))
  P_Flag = 1

  allocate (D_Flag(NNODE))
  D_Flag = 0

  allocate (ELMNGAUSS(NELEM))
  ELMNGAUSS = ELMNSHL
  do i = 1, NBOUND
    allocate (bound(i)%NGAUSSB(bound(i)%NFACE))
    bound(i)%NGAUSSB = bound(i)%NSHLB
!     if(ismaster) write(*,*) i, bound(i)%NNODE
  end do

  !----------------------------------------------------------
  ! read the mapping between partitioned volume boundary
  ! and unpartitioned shell mesh (bmesh.*.dat)
  !----------------------------------------------------------
  write (iname, '(I8)') id
  fname = 'l2b.'//trim(adjustl(iname))//'.dat'
  open (mfid, file=fname, status='old')

  do i = 1, NBOUND

    read(mfid, *) itmp1, itmp2, itmp3
    if (itmp1 /= bound(i)%FACE_ID) then
      write(*,*) "ERROR in l2b: Boundary ID doesn't match"
      stop
    end if
    if (itmp2 /= bound(i)%NNODE) then
      write(*,*) "ERROR in l2b: Boundary node number doesn't match"
      stop
    end if
    if (itmp3 /= bound(i)%NFACE) then
      write(*,*) "ERROR in l2b: Boundary face number doesn't match"
      stop
    end if

    !if (bound(i)%NNODE > 0) then
      allocate (bound(i)%L2SNODE(bound(i)%NNODE))
      do j = 1, bound(i)%NNODE
        read (mfid, *) itmp1, bound(i)%L2SNODE(j) ! this is local node ID to the global bound node ID
      end do
    !end if

    !if (bound(i)%NFACE > 0) then
      allocate (bound(i)%L2SELEM(bound(i)%NFACE))

      do j = 1, bound(i)%NFACE
        ! read partitioned boundary element number and the corresponding
        ! shell mesh element number
        read(mfid,*) itmp1, bound(i)%L2SELEM(j)
        if (itmp1 /= bound(i)%F2E(j)) then
          write(*,*) "ERROR in l2b: Boundary element doesn't match"
          stop
        end if
      end do
    !end if

  end do
  close (mfid)
  call get_mesh_statistics()

end subroutine input_mesh
!======================================================================
!
!======================================================================
subroutine input_bmesh()
  use aAdjKeep
  use commonvars
  use mpi
  implicit none
  ! read the coordinates of bmesh only
  character(len=32) :: fname, cname
  integer :: i, j, k, ifac
  integer :: ib, mfid, tmp
  integer,allocatable :: owned_local(:)
  logical, parameter :: debug = .true.
  real(8) :: xb(3), xs(3)
  do ib = 1, NRES
    write(cname, '(I8)') bmesh_id_set(ib)
    fname = 'bmesh.'//trim(adjustl(cname))//'.dat'
    mfid = 11
    open (mfid, file=fname, status='old')
    read (mfid, *) tmp, tmp
    read (mfid, *) bmesh(ib)%NNODE, tmp, tmp, bmesh(ib)%FACE_ID
    allocate (bmesh(ib)%xg(bmesh(ib)%NNODE, NSD))
    allocate (bmesh(ib)%owned(bmesh(ib)%NNODE))
    do i = 1, bmesh(ib)%NNODE
      read (mfid, *) (bmesh(ib)%XG(i, j), j=1, NSD)
    end do
    if(ismaster) write(*,*) "Reading bmesh file: ", fname, bmesh(ib)%NNODE, bmesh(ib)%FACE_ID
    ! find the boundary ID
    do i = 1, NBOUND
      if(bmesh(ib)%FACE_ID == bound(i)%FACE_ID) then
        ifac = i
        exit
      end if
    enddo
    ! check if the node ID is correct
    if(debug) then
      do i = 1, bound(ifac)%NNODE
        xb(:) = xg(bound(ifac)%BNODES(i), :)
        k = bound(ifac)%L2SNODE(i)
        xs(:) = bmesh(ib)%xg(k, :)
        if(norm2(xb - xs) > 1d-6) then
          write(*,*) "ERROR: bmesh/bound node doesn't match", myid
          write(*,*) "bmesh node", k, bmesh(ib)%xg(k, :)
          write(*,*) "bound node", i, xg(bound(ifac)%BNODES(i), :)
          call MPI_Abort(MPI_COMM_WORLD, 201, mpi_err)
        endif
      enddo
    endif
    allocate(owned_local(bmesh(ib)%NNODE))
    owned_local(:) = 0
    ! assert
    if(debug .and. maxval(bound(ifac)%L2SNODE(:)) > bmesh(ib)%NNODE) then
      write(*,*) "ERROR: bmesh node ID is larger than boundary node ID", myid
      !write(*,*) "bmesh node ID", maxval(bound(ifac)%L2SNODE(:))
      !write(*,*) "boundary node ID", bmesh(ib)%NNODE
      ! call MPI_Abort(MPI_COMM_WORLD, 201, mpi_err)
    end if
    do i = 1, bound(ifac)%NNODE
      owned_local(bound(ifac)%L2SNODE(i)) = myid
    end do
    ! get the maximum of owned_local
    call MPI_AllREDUCE(owned_local(1), bmesh(ib)%owned(1), bmesh(ib)%NNODE, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpi_err)
    do i = 1, bmesh(ib)%NNODE
      if (bmesh(ib)%owned(i) /= myid) then
        bmesh(ib)%owned(i) = 0
      else
        bmesh(ib)%owned(i) = 1
      end if
    end do
    if(debug) then
      write(*,*) "owned", myid, bmesh(ib)%NNODE
    endif
    deallocate(owned_local)
    close (mfid)
  enddo

end subroutine input_bmesh



!======================================================================
!
!======================================================================
subroutine get_mesh_statistics()
  use commonvars
  use aAdjKeep
  use mpi
  implicit none

  integer :: i
  integer :: nodes(2)

  nodes(:) = 0

  do i = 1, NNODE
    if(NODEID(i) == ELM_ID_SET(1)) then
      nodes(1) = nodes(1) + 1
    else if(NODEID(i) == ELM_ID_SET(2)) then
      nodes(2) = nodes(2) + 1
    end if
  end do

  write(*,*) "Node distribution of", myid, nodes(:), nodes(1) + nodes(2)
end subroutine get_mesh_statistics
