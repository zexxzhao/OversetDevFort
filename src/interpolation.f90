!==========================================================
! Brutal force interpolation without any optimization
!==========================================================

subroutine interpolation_tet(mesh_id)
  use aAdjkeep
  use commonvars
  use commonpars
  use mpi
  implicit none

  integer, intent(IN) :: mesh_id

  integer :: iel, nshl, i, j, k

  real(8) :: xl(4, NSD), dl(4, NSD)
  real(8) :: Amat(3, 3), Ainv(3, 3), bvec(3), xvec(4), tmp
  integer :: found(bmesh(mesh_id)%NNODE)
  real(8) :: ug_local(bmesh(mesh_id)%NNODE, NSD)
  found(:) = 0
  bmesh(mesh_id)%ug(:, :) = 0d0
  ug_local(:, :) = 0d0
  ! loop over elements for collision detection
  do iel = 1, NELEM
    if(ELM_ID(iel) /= SRC_ID_SET(mesh_id)) cycle
    nshl = ELMNSHL(iel)
    if(nshl > 4) cycle ! only for tetrahedron
    
    ! get the coordinates of the element
    do j = 1, NSHL
      xl(j, :) = xg(IEN(iel, j), :)
      dl(j, :) = dg(IEN(iel, j), :)
    enddo
    ! use xl(1, :) as the origin and store the relative coordinates in Amat
    do k = 1, NSD
      do j = 2, nshl
        Amat(j-1, k) = xl(j, k) - xl(1, k) + dl(j, k) - dl(1, k) 
      enddo
    enddo
    Amat = transpose(Amat)
    ! calculate the inverse of Amat
    call get_inverse_3x3(Amat, Ainv, tmp)

    bmesh(mesh_id)%ug(:, :) = 0
    do j = 1, bmesh(mesh_id)%NNODE
      if(found(j) > 0) cycle
      ! bvec(:) = xg(j, :) - xl(1, :) + dg(j, :) - dl(1, :)
      bvec(:) = bmesh(mesh_id)%xg(j, :) - xl(1, :) &
        + bmesh(mesh_id)%dg(j, :)
      xvec(1:3) = matmul(Ainv, bvec)
      xvec(4) = 1d0 - sum(xvec(1:3))
      if(all((/(xvec(i) >= 0, i = 1, 4)/))) then
        ug_local(j, :) = ug_local(j, :) + xvec(1) * ug(IEN(iel, 1), :)
        ug_local(j, :) = ug_local(j, :) + xvec(2) * ug(IEN(iel, 2), :)
        ug_local(j, :) = ug_local(j, :) + xvec(3) * ug(IEN(iel, 3), :)
        ug_local(j, :) = ug_local(j, :) + xvec(4) * ug(IEN(iel, 4), :)
        found(j) = found(j) + 1
      endif
    enddo
  enddo
  
  call MPI_ALLREDUCE(ug_local(1, 1), bmesh(mesh_id)%ug(1, 1), &
    bmesh(mesh_id)%NNODE * NSD, MPI_DOUBLE_PRECISION, &
    MPI_SUM, MPI_COMM_WORLD)
  
end subroutine interpolation_tet
!==========================================================
!
!==========================================================
subroutine update_bmesh_displacement(mesh_id)
  use aAdjkeep
  use commonvars
  use commonpars
  use mpi
  implicit none

  integer, intent(IN) :: mesh_id
  real(8) :: dg_local(bmesh(mesh_id)%NNODE, NSD)
  integer :: i, j
  dg_local(:, :) = 0
  ! insert dg into dg_local
  do i = 1, bmesh(mesh_id)%NNODE
    if(bmesh(mesh_id)%owned(i) > 0) then
      j = bound(ELM_ID_SET(mesh_id))%BNODES(i)
      dg_local(i, :) = dg(j, :) 
    endif
  enddo
  call MPI_Allreduce(&
    dg_local(1,1), bmesh(mesh_id)%dg(1,1), &
    bmesh(mesh_id)%NNODE * NSD, MPI_DOUBLE_PRECISION, &
    MPI_SUM, MPI_COMM_WORLD)
  
end subroutine update_bmesh_displacement
