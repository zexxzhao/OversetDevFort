!======================================================================
!
!======================================================================
subroutine solveOverset_stag(istep)

  use aAdjKeep
  use mpi
  use commonvars
  use commonpars
  implicit none

  integer, intent(in) :: istep

  integer :: inewt, i, NL_max
  real(8) :: residual0(NRES*NSubdomain), residual(NRES*NSubdomain)
  integer :: converged(NRES*NSubdomain)
  integer :: subdomain_flag

  real(8) :: sol(NNODE, NRES + NSD - 1)
  real(8) :: t1, t2
  real(8) :: utol(NRES*NSubdomain)
  
  utol = (/NS_NL_UTOL, NS_NL_PTOL, NS_NL_UTOL, NS_NL_PTOL/)
  inewt = 0
  ! momres0 = -1.0d0
  ! conres0 = -1.0d0
  ! convres0 = -1.0d0
  ! lsres0 = -1.0d0

  residual(:) = 0d0
  residual0(:) = 0d0

  subdomain_flag = 1 + 2
  IBC(:, :) = 0
  call setBCs_overset(subdomain_flag)
  call assembleNavSto(ASSEMBLE_TENSOR_VEC, subdomain_flag)

  call calculate_residual_overset( &
    residual0, rhsgu, rhsgp, nodeid, elm_id_set, &
    subdomain_flag, nnode, nsd, nres, nsubdomain)

  residual(:) = residual0(:)
  call check_convergence(converged, residual, residual0, utol)
  call print_residual(residual, residual0, utol, inewt)

  if(istep <= 1) then 
    NL_max = 1
  else 
    NL_max = NS_NL_itermax
  endif
  do inewt = 1, NL_max
    !---------------------------
    ! Solve Mesh 1
    !---------------------------
    subdomain_flag = 1
    IBC(:, :) = 0
    call setBCs_overset(subdomain_flag)
    ! rotate Mesh 1 since it contains the rotating blades
    call rotate_mesh(subdomain_flag)
    call assembleNavSto(ASSEMBLE_TENSOR_VEC + ASSEMBLE_TENSOR_MAT, subdomain_flag)

    if (myid .eq. 0) then
      call CPU_TIME(t1)
    end if
    sol(:, 1:4) = 0d0
    call SparseGMRES_up(col, row, IBC, IPER, D_FLAG, P_FLAG, &
                        rhsGu, rhsGp, sol(:, 1:3), sol(:, 4), &
                        lhsK11, lhsG, lhsD1, lhsM, icnt, &
                        NS_GMRES_tol, NS_GMRES_itermax, &
                        NS_GMRES_itermin, &
                        NNODE, maxNSHL, NSD)
    if (myid .eq. 0) then
      call CPU_TIME(t2)
      write (*, *) "Total time solve GMRES Mesh#1:", t2 - t1, "seconds"
    end if
    acg = acg + sol(:, 1:3)
    ug = ug + gami*Delt*sol(:, 1:3)
    pg = pg + alfi*gami*Delt*sol(:, 4)
    !---------------------------
    ! Solve Mesh 2
    !---------------------------
    subdomain_flag = 2
    IBC(:, :) = 0
    call setBCs_overset(subdomain_flag)
    call assembleNavSto(ASSEMBLE_TENSOR_VEC + ASSEMBLE_TENSOR_MAT, subdomain_flag)

    if (myid .eq. 0) then
      call CPU_TIME(t1)
    end if
    sol(:, 1:4) = 0d0
    call SparseGMRES_up(col, row, IBC, IPER, D_FLAG, P_FLAG, &
                        rhsGu, rhsGp, sol(:, 1:3), sol(:, 4), &
                        lhsK11, lhsG, lhsD1, lhsM, icnt, &
                        NS_GMRES_tol, NS_GMRES_itermax, &
                        NS_GMRES_itermin, &
                        NNODE, maxNSHL, NSD)
    if (myid .eq. 0) then
      call CPU_TIME(t2)
      write (*, *) "Total time solve GMRES Mesh#2:", t2 - t1, "seconds"
    end if
    acg = acg + sol(:, 1:3)
    ug = ug + gami*Delt*sol(:, 1:3)
    pg = pg + alfi*gami*Delt*sol(:, 4)

    !---------------------------
    ! Check Convergence
    !---------------------------
    subdomain_flag = 1 + 2
    IBC(:, :) = 0
    call setBCs_overset(subdomain_flag)
    call assembleNavSto(ASSEMBLE_TENSOR_VEC, subdomain_flag)


    call calculate_residual_overset( &
      residual, RHSGu, RHSGp, NODEID, ELM_ID_SET, &
      subdomain_flag, NNODE, NSD, NRES, NSubdomain)
    call check_convergence(converged, residual, residual0, utol)
    call print_residual(residual, residual0, utol, inewt)

    if(ismaster ) write (*, *) "Convergence:", converged
    if (size(converged) == sum(converged)) exit

    if (numnodes > 1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
  end do

end subroutine solveOverset_stag

!============================================================================
!
!============================================================================
subroutine calculate_residual_overset(&
  r, RHSGu, RHSGp, NODE_ID, NODE_ID_SET, &
  subdomain_flag, NNODE, NSD, NRES, NSubdomain)
  implicit none
  include "mpif.h"

  real(8), intent(in) :: RHSGu(NNODE, NSD), RHSGp(NNODE)
  integer, intent(in) :: NODE_ID(NNODE), NODE_ID_SET(NRES)
  integer, intent(in) :: subdomain_flag, NNODE, NSD, NRES, NSubdomain
  real(8), intent(out) :: r(NRES * NSubdomain)

  real(8) :: r_local(NRES * NSubdomain)
  integer :: j, mpi_err

  r_local(:) = 0d0
  ! sum of RHSGu(j, :) ** 2 where NODE_ID_SET(1) == NODE_ID(j)
  if(iand(subdomain_flag, 1) > 0) then
    do j = 1, NNODE
      if (NODE_ID_SET(1) == NODE_ID(j)) then
        r_local(1) = r_local(1) + sum(RHSGu(j, :)**2)
        r_local(2) = r_local(2) + RHSGp(j)**2
      end if
    end do
  end if

  ! sum of RHSGu(j, :) ** 2 where NODE_ID_SET(2) == NODE_ID(j)
  if(iand(subdomain_flag, 2) > 0) then
    do j = 1, NNODE
      if (NODE_ID_SET(2) == NODE_ID(j)) then
        r_local(3) = r_local(3) + sum(RHSGu(j, :)**2)
        r_local(4) = r_local(4) + RHSGp(j)**2
      end if
    end do
  end if  


  call MPI_ALLREDUCE(r_local(1), r(1), NRES*NSubdomain, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)

  do j = 1, NRES*NSubdomain
    r(j) = sqrt(r(j))
  end do

end subroutine calculate_residual_overset

!============================================================================
!
!============================================================================
subroutine check_convergence(converged, residual, residual0, tol)
  use commonvars
  use commonpars
  implicit none

  integer, intent(inout) :: converged(NRES*NSubdomain)

  real(8), intent(in) :: residual(NRES*NSubdomain), residual0(NRES*NSubdomain), tol(NRES*NSubdomain)
  integer :: j

  converged(:) = 0
  do j = 1, NRES*NSubdomain
    if (residual(j) < tol(j) * max(1e-15, residual0(j))) then
      converged(j) = 1
    end if
  end do
end subroutine check_convergence

!============================================================================
!
!============================================================================
subroutine print_residual(residual, residual0, tol, inewt)
  use commonpars
  use mpi
  implicit none

  real(8), intent(in) :: residual(NRES*NSubdomain), residual0(NRES*NSubdomain), tol(NRES*NSubdomain)
  integer, intent(in) :: inewt

  if(ismaster) then
    write(*, 9999) inewt, "Mesh #1, MOM", residual(1), residual(1)/max(residual0(1), 1d-15), tol(1)
    write(*, 9999) inewt, "Mesh #1, CON", residual(2), residual(2)/max(residual0(2), 1d-15), tol(2)
    write(*, 9999) inewt, "Mesh #2, MOM", residual(3), residual(3)/max(residual0(3), 1d-15), tol(3)
    write(*, 9999) inewt, "Mesh #2, CON", residual(4), residual(4)/max(residual0(4), 1d-15), tol(4)
  endif

9999 format(i2, ') ', a12, '  abs=', e10.4, '  rel=', e10.4, '  rtol=', e10.4)
end subroutine print_residual

!============================================================================
!
!============================================================================
subroutine print_sol()
  use aAdjKeep
  use commonvars
  use commonpars
  use mpi

  implicit none

  real(8) :: umagmax_local, umagmin_local, umagmax, umagmin
  real(8) :: pmax_local, pmin_local, pmax, pmin
  real(8) :: umag(NNODE)
  ! integer :: mpi_err

  umag = sqrt(ug(:, 1)**2 + ug(:, 2)**2 + ug(:, 3)**2)

  umagmax_local = maxval(umag)
  umagmin_local = minval(umag)
  pmax_local = maxval(pg)
  pmin_local = minval(pg)

  call MPI_ALLREDUCE(umagmax_local, umagmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpi_err)
  call MPI_ALLREDUCE(umagmin_local, umagmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpi_err)
  call MPI_ALLREDUCE(pmax_local, pmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpi_err)
  call MPI_ALLREDUCE(pmin_local, pmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpi_err)
  if (ismaster) then
    write(*,*) "|U|min, |U|max [m/s]: ", umagmin, umagmax
    write(*,*) "Pmin, Pmax [Pa]: ", pmin, pmax
  end if

end subroutine print_sol

!============================================================================
! Calculate acgm, ugm, dg using specified `omega_x`
!============================================================================

subroutine rotate_mesh(meshid)
  use commonvars
  use commonpars
  use mpi
  use aAdjKeep
  implicit none
  
  integer, intent(in) :: meshid
  integer :: i

  real(8) :: xold(NSD), xnew(NSD)

  theta = omega_x * Delt
  do i = 1,NNODE
    if(NODEID(i) == ELM_ID_SET(meshid)) then
      ! rotate the mesh along x-axis by theta
      xold(:) = xg(i, :) + dgold(i, :)
      xnew(1) = xold(1)
      xnew(2) = xold(2) * cos(theta) - xold(3) * sin(theta)
      xnew(3) = xold(2) * sin(theta) + xold(3) * cos(theta)
      dg(i, :) = xnew - xg(i, :)
      ! calculate the acgm and ugm based on generalized alpha and newmark-beta method
      ! where acgm is the n+1 acceleration, ugm is the n+1 velocity, dg is the n+1 displacement
      ! acgmold is the n acceleration, ugmold is the n velocity, dgold is the n displacement

      acgm(i, :) = (dg(i,:) - dgold(i, :) - delt * ugmold(i, :)) / (delt**2 * beti) &
        - (1-2*beti) / (2*beti) * acgmold(i, :)
      ugm(i, :) = ugmold(i, :) + delt * (1-gami) * acgmold(i, :) + delt * gami * acgm(i, :)
    else
      dg(i, :) = 0d0
      acgm(i, :) = 0d0
      ugm(i, :) = 0d0
    endif
  enddo

  
end subroutine rotate_mesh
