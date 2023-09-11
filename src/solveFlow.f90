!======================================================================
!
!======================================================================
subroutine solmultiphasethermofluid_stag(istep)

  use aAdjKeep
  use mpi
  use commonvars
  use commonpars
  implicit none

  integer, intent(in) :: istep

  integer :: inewt, i, NL_max
  ! real(8) :: momres0, conres0, convres0, meshres0, lsres0
  real(8) :: residual0(4), residual(4)
  integer :: converged(4)

  real(8) :: sol(NNODE, 6)
  real(8) :: t1, t2
  real(8) :: utol(4)
  
  utol = (/NS_NL_UTOL, NS_NL_PTOL, LSC_NL_TOL, TEM_NL_TOL/)
  inewt = 0
  ! momres0 = -1.0d0
  ! conres0 = -1.0d0
  ! convres0 = -1.0d0
  ! lsres0 = -1.0d0

  residual(:) = 0d0
  residual0(:) = 0d0
  converged(:) = 0

  IBC(:, :) = 0
  call setBCs_NSVOF()
  call setBCs_Tem()
  !call assembleNavStoVOFTem(ASSEMBLE_TENSOR_VEC, &
  !                          ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)
  call assembleQuenching(ASSEMBLE_TENSOR_VEC, &
                         ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)

  ! write(*,*) "myid=", myid, "GET RHS0"
  call calculate_residual(residual0, RHSGu, RHSGp, RHSGls, RHSGtem, NNODE, NSD)
  residual(:) = residual0(:)
  call check_convergence(converged, residual, residual0, utol, &
                         ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)
  call print_residual(residual, residual0, utol, &
                      ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM, &
                      inewt)

  if(istep <= 1) then 
    NL_max = 1
  else 
    NL_max = NS_NL_itermax
  endif
  do inewt = 1, NL_max
    IBC(:, :) = 0
    call setBCs_NSVOF()
    call setBCs_Tem()
 
    !---------------------------
    ! Solve NavStoVOF
    !---------------------------
    ! IBC(:, :) = 0
    ! call setBCs_NSVOF()
    ! IBC(:, 5) = 1
    ! IBC(:, 6:8) = 1
    ! call assembleNavStoVOFTem(ASSEMBLE_TENSOR_MAT + ASSEMBLE_TENSOR_VEC, &
    !                           ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF)
 
    call assembleQuenching(ASSEMBLE_TENSOR_MAT + ASSEMBLE_TENSOR_VEC, &
                           ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF)
    call calculate_residual(residual, RHSGu, RHSGp, RHSGls, RHSGtem, NNODE, NSD)
    !call check_convergence(converged, residual, residual0, utol, &
    !                      ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF)
    !call print_residual(residual, residual0, utol, &
    !                    ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF, &
    !                    inewt)

    if (myid .eq. 0) then
      call CPU_TIME(t1)
    end if
    sol(:, 1:5) = 0d0
    !if((residual(1)< 1d-16 ).and. (residual(2)< 1d-16 ).and. (residual(3)< 1d-16 )) then
    call SparseGMRES_up(col, row, IBC, IPER, D_FLAG, P_FLAG, &
                        rhsGu, rhsGp, sol(:, 1:3), sol(:, 4), &
                        lhsK11, lhsG, lhsD1, lhsM, icnt, &
                        NS_GMRES_tol, NS_GMRES_itermax, &
                        NS_GMRES_itermin, &
                        NNODE, maxNSHL, NSD, &
                        sol(:, 5), lhsLS, &
                        lhsLSu, lhsUls, lhsPls, rhsGls)
    !endif

    if (myid .eq. 0) then
      call CPU_TIME(t2)
      write (*, *) "Total time solve NSVOF GMRES:", t2 - t1, "seconds"
    end if
    acg = acg + sol(:, 1:3)
    ug = ug + gami*Delt*sol(:, 1:3)
    pg = pg + alfi*gami*Delt*sol(:, 4)

    rphig = rphig + sol(:, 5)
    phig = phig + gami*Delt*sol(:, 5)

    !-----------------------------
    ! Solve Temperature
    !-----------------------------
    if (istep > 0) then
      ! IBC(:, :) = 0
      ! call setBCs_Tem()
      ! IBC(:, 1:5) = 1
      ! IBC(:, 7:8) = 1
      ! call assembleNavStoVOFTem(ASSEMBLE_TENSOR_MAT + ASSEMBLE_TENSOR_VEC, &
      !                           ASSEMBLE_FIELD_TEM)
      call assembleQuenching(ASSEMBLE_TENSOR_MAT + ASSEMBLE_TENSOR_VEC, &
                             ASSEMBLE_FIELD_TEM)
      ! call calculate_residual(residual, RHSGu, RHSGp, RHSGls, RHSGtem, NNODE, NSD)
      ! call check_convergence(converged, residual, residual0, utol, &
      !                        ASSEMBLE_FIELD_TEM)
      ! call print_residual(residual, residual0, utol, &
      !                     ASSEMBLE_FIELD_TEM, &
      !                     inewt)

      !write(*,*) "Residual T = ", residual0(4), residual(4)
      if (myid .eq. 0) then
        call CPU_TIME(t1)
      end if

      sol(:, 6) = 0d0
      call SparseGMRES_tem(LHStem, TEM_GMRES_tol, col, row, &
                           rhsgtem, sol(:, 6), TEM_GMRES_itermax, TEM_GMRES_itermin, &
                           NNODE, maxNSHL, icnt, NSD)
      if (myid .eq. 0) then
        call CPU_TIME(t2)
        write (*, *) "Total time solve TEM GMRES:", t2 - t1, "seconds"
      end if
      !write(*,*) "solT =", sum(sol(:, 6)**2) 
      rTg = rTg + sol(:, 6)
      Tg = Tg + gami*Delt*sol(:, 6)

    end if
    ! IBC(:, :) = 0
    ! call setBCs_NSVOF()
    ! call setBCs_Tem()
    ! call assembleNavStoVOFTem(ASSEMBLE_TENSOR_VEC, &
    !                        ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)
    call assembleQuenching(ASSEMBLE_TENSOR_VEC, &
                           ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)

    call calculate_residual(residual, RHSGu, RHSGp, RHSGls, RHSGtem, NNODE, NSD)
    call check_convergence(converged, residual, residual0, utol, &
                           ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)
    call print_residual(residual, residual0, utol, &
                        ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM, &
                        inewt)

    if(ismaster ) write (*, *) "Convergence:", converged
    if (size(converged) == sum(converged)) exit

    if (numnodes > 1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
  end do

end subroutine solmultiphasethermofluid_stag

!============================================================================
!
!============================================================================
subroutine calculate_residual(r, RHSGu, RHSGp, RHSGls, RHSGtem, NNODE, NSD)
  implicit none
  include "mpif.h"

  real(8), intent(in) :: RHSGu(NNODE, NSD), RHSGp(NNODE), RHSGls(NNODE), RHSGtem(NNODE)
  integer, intent(in) :: NNODE, NSD
  real(8), intent(out) :: r(4)

  real(8) :: r_local(4)
  integer :: mpi_err

  r_local(1) = sum(RHSGu(:, :) * RHSGu(:, :))
  r_local(2) = sum(RHSGp(:) * RHSGp(:))
  r_local(3) = sum(RHSGls(:) * RHSGls(:))
  r_local(4) = sum(RHSGtem(:) * RHSGtem(:))

  call MPI_ALLREDUCE(r_local(1), r(1), 4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)

end subroutine calculate_residual

!============================================================================
!
!============================================================================
subroutine check_convergence(converged, residual, residual0, tol, assemble_field_flag)
  use commonvars
  use commonpars
  implicit none

  integer, intent(out) :: converged(4)

  real(8), intent(in) :: residual(4), residual0(4), tol(4)
  integer, intent(in) :: assemble_field_flag

  converged(:) = 0
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_NS) > 0) then
    if ((residual(1) < tol(1) * residual0(1)).or. ( residual(1)<NS_GMRES_atol))then
      converged(1) = 1
    end if
    if ((residual(2) < tol(2) * residual0(2)).or. ( residual(2)<NS_GMRES_atol)) then
      converged(2) = 1
    end if
  end if
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_VOF) > 0)then
    if ((residual(3) < tol(3)*residual0(3)) .or. ( residual(3)< LSC_GMRES_atol))then
      converged(3) = 1
    end if
  end if
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0) then
    if ((residual(4) < tol(4) * residual0(4) ).or.  ( residual(4)<TEM_GMRES_atol))then
      converged(4) = 1
    end if
  end if

end subroutine check_convergence

!============================================================================
!
!============================================================================
subroutine print_residual(residual, residual0, tol, assemble_field_flag, inewt)
  use commonpars
  use mpi
  implicit none

  real(8), intent(in) :: residual(4), residual0(4), tol(4)
  integer, intent(in) :: assemble_field_flag, inewt

  if (iand(assemble_field_flag, ASSEMBLE_FIELD_NS) > 0 .and. ismaster) then
    write(*,9999) inewt, "MOM", residual(1), residual(1)/max(residual0(1), 1d-15), tol(1)
    write(*,9999) inewt, "CON", residual(2), residual(2)/max(residual0(2), 1d-15), tol(2)
  end if
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_VOF) > 0 .and. ismaster) then
    write(*,9999) inewt, "VOF", residual(3), residual(3)/max(residual0(3), 1d-15), tol(3)
  end if
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0 .and. ismaster) then
    write(*,9999) inewt, "TEM", residual(4), residual(4)/max(residual0(4), 1d-15), tol(4)
  end if

9999 format(i2, ') ', a3, '  abs = ', e10.4, '  rel = ', e10.4, '  rtol = ', e10.4)
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

  real(8) :: phimin, phimax, phimin_local, phimax_local
  real(8) :: temmin, temmax, temmin_local, temmax_local

  ! integer :: mpi_err

  phimin_local = minval(phig)
  phimax_local = maxval(phig)

  temmin_local = minval(Tg)
  temmax_local = maxval(Tg)

  call MPI_ALLREDUCE(phimin_local, phimin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpi_err)
  call MPI_ALLREDUCE(phimax_local, phimax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpi_err)

  call MPI_ALLREDUCE(temmin_local, temmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpi_err)
  call MPI_ALLREDUCE(temmax_local, temmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpi_err)

  if (ismaster) then
    write(*,*) "phimin, phimax [ ] = ", phimin, phimax
    write(*,*) "temmin, temmax [C] = ", temmin-273d0, temmax-273d0
  end if

end subroutine print_sol
