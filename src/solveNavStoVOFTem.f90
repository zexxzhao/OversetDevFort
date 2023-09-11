!======================================================================
!
!======================================================================

! subroutine calculate_residual(residual)
!   use mpi
!   use aAdjKeep
!   implicit none
!   integer, parameter :: NRES = 4
! 
!   real(8), intent(out) :: residual(4)
! 
!   real(8) :: rl(NRES)
!   integer :: i
! 
!   rl(1) = sum(RHSGu(:, :)*RHSGu(:, :))
!   rl(2) = sum(RHSGp(:)*RHSGp(:))
!   rl(3) = sum(RHSGls(:)*RHSGls(:))
!   rl(4) = sum(RHSGtem(:)*RHSGtem(:))
!   if (numnodes > 1) then
!     call MPI_ALLREDUCE(rl, residual, NRES, MPI_DOUBLE_PRECISION, &
!                        MPI_SUM, MPI_COMM_WORLD, mpi_err)
!   else
!     residual(:) = rl(:)
!   end if
! 
!   do i = 1, NRES
!     residual(i) = sqrt(residual(i))
!   end do
! 
! end subroutine calculate_residual

!======================================================================
!
!======================================================================

! subroutine check_convergence(r0, r, verbose, inewt, assemble_field_flag, converged)
!   use mpi
!   use aAdjKeep
!   use commonvars
!   use commonpars
!   implicit none
! 
!   integer, parameter :: NRES = 4
! 
!   real(8), intent(in) :: r0(NRES), r(NRES)
!   integer, intent(in) :: inewt
!   integer, intent(in) :: verbose
!   integer, intent(in) :: assemble_field_flag
! 
!   integer, intent(out) :: converged(NRES)
!   character(len=80) :: fomt
! 
!   converged(:) = 0
! 
!   if (r0(1)*NS_NL_UTOL >= r(1)) converged(1) = 1
!   if (r0(2)*NS_NL_PTOL >= r(2)) converged(2) = 1
!   if (r0(3)*LSC_NL_TOL >= r(3)) converged(3) = 1
!   if (r0(4)*LSC_NL_TOL >= r(4)) converged(4) = 1
!   if (ismaster .and. verbose > 0) then
!   fomt = "(I3,a,x,ES13.6,x,F12.6)"
!   if (iand(assemble_field_flag, ASSEMBLE_FIELD_NS) > 0) then
!     write (*, fomt) inewt, ") Total Mom. Res. Norm = ", &
!       r(1), 1.0d2*r(1)/max(r0(1), 1d-15)
!     write (*, fomt) inewt, ") Continuity Res. Norm = ", &
!       r(2), 1.0d2*r(2)/max(r0(2), 1d-15)
!   end if
! 
!   if (iand(assemble_field_flag, ASSEMBLE_FIELD_LS) > 0) then
!     write (*, fomt) inewt, ") LS Res. Norm = ", &
!       r(3), 1.0d2*r(3)/max(r0(3), 1d-15)
!   end if
!   if (iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0) then
!     write (*, fomt) inewt, ") TEM Res. Norm = ", &
!       r(4), 1.0d2*r(4)/max(r0(4), 1d-15)
!   end if
!   write (*, *)
!   end if
! 
! end subroutine check_convergence

!======================================================================
!
!======================================================================
! subroutine print_residual(r, r0, utol, assemble_field_flag, inewt)
! 
!   use commonpars
!   use mpi
! 
!   implicit none
! 
!   real(8), intent(in) :: r(4), r0(4), utol(4)
!   integer, intent(in) :: assemble_field_flag
!   integer, intent(in) :: inewt
! 
!   character(len=80) :: fomt
!   fomt = "(I3,a,x,ES13.6,x,F12.6,a,F12.6,a)"
!   if(ismaster) then
!     if(iand(assemble_field_flag, ASSEMBLE_FIELD_NS) > 0) then
!       write(*,fomt) inewt, ") Total Mom. Res. Norm = ", &
!         r(1), 1.0d2*r(1)/max(r0(1), 1d-15), "%", 1.0d2*utol(1), "%"
!       write(*,fomt) inewt, ") Continuity Res. Norm = ", &
!         r(2), 1.0d2*r(2)/max(r0(2), 1d-15), "%", 1.0d2*utol(2), "%"
!     endif
!     if(iand(assemble_field_flag, ASSEMBLE_FIELD_LS) > 0) then
!       write(*,fomt) inewt, ") LS Res. Norm = ", &
!         r(3), 1.0d2*r(3)/max(r0(3), 1d-15), "%", 1.0d2*utol(3), "%"
!     endif
!     if(iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0) then
!       write(*,fomt) inewt, ") TEM Res. Norm = ", &
!         r(4), 1.0d2*r(4)/max(r0(4), 1d-15), "%", 1.0d2*utol(4), "%"
!     endif
!   endif
! end subroutine print_residual

!======================================================================
!
!======================================================================
subroutine assembleNavStoVOFTem(assemble_tensor_flag, assemble_field_flag)
  use aAdjKeep
  use mpi
  use commonvars
  use commonpars

  implicit none

  integer, intent(in) :: assemble_tensor_flag ! assembe Jacobian mat or vec
  integer, intent(in) :: assemble_field_flag ! assemble NS + LS/VOF + Tem

  real(8) :: dgAlpha(NNODE, NSD), &
             ugAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
             acgAlpha(NNODE, NSD), acgmAlpha(NNODE, NSD), &
             pgAlpha(NNODE), phigAlpha(NNODE), rphigAlpha(NNODE), &
             rTgAlpha(NNODE), TgAlpha(NNODE)

  real(8) :: t1, t2

  !---------------------------------
  ! Alpha stage
  !---------------------------------
  acgAlpha = acgold + almi*(acg - acgold)
  acgmAlpha = acgmold + almi*(acgm - acgmold)
  ugAlpha = ugold + alfi*(ug - ugold)
  ugmAlpha = ugmold + alfi*(ugm - ugmold)
  dgAlpha = dgold + alfi*(dg - dgold)
  pgAlpha = pg

  phigAlpha = phigold + alfi*(phig - phigold)
  rphigAlpha = rphigold + almi*(rphig - rphigold)

  TgAlpha = Tgold + alfi*(Tg - Tgold)
  rTgAlpha = rTgold + almi*(rTg - rTgold)

  !---------------------------------
  ! zero out RHS and LHS
  !---------------------------------

  if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
    RHSGu = 0.0d0
    RHSGp = 0.0d0
    RHSGls = 0.0d0
    RHSGtem = 0.0d0
  end if

  if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
    LHSK11 = 0.0d0
    LHSG = 0.0d0
    LHSD1 = 0.0d0
    LHSM = 0.0d0

    LHSLS = 0.0d0
    LHSULS = 0.0d0
    LHSLSU = 0.0d0
    LHSPLS = 0.0d0

    LHSTem = 0.0d0
  endif

  if (myid .eq. 0) then
    call CPU_TIME(t1)
  endif
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF) > 0) then
    ! write(*,*) myid, "ug-alpha", sum(ugAlpha(:, :) ** 2), assemble_tensor_flag 
    call IntElmAss_NSVOF(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                         acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
                         TgAlpha, rTgAlpha, &
                         assemble_tensor_flag)
    ! write(*,*) myid, "RHSgu1", sum(RHSGu(:, :) ** 2), &
    !       sum(LHSK11**2), sum(lhsG**2), sum(lhsD1**2), sum(lhsM**2)
    call FaceAssembly_NS_weak(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                              acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
                              assemble_tensor_flag)
    ! write(*,*) myid, "RHSgu2", sum(RHSGu(:, :) ** 2), &
    !      sum(LHSK11**2), sum(lhsG**2), sum(lhsD1**2), sum(lhsM**2)
  end if
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0) then
    call IntElmAss_Tem(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                       acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
                       TgAlpha, rTgAlpha, &
                       assemble_tensor_flag)

  end if

  if (myid .eq. 0) then
    call CPU_TIME(t2)
    write (*, *) "Total time assemble:", t2 - t1, "seconds"
  end if

  if (numnodes > 1 .and. iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
    if(iand(assemble_field_flag, ASSEMBLE_FIELD_NS) > 0) then
      call commu(RHSGp, 1, 'in ')
      call commu(RHSGu, NSD, 'in ')
    endif
    if(iand(assemble_field_flag, ASSEMBLE_FIELD_VOF) > 0) then
      call commu(RHSGls, 1, 'in ')
    endif
    if(iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0) then
      call commu(RHSGTem, 1, 'in ')
    endif
  end if

end subroutine assembleNavStoVOFTem
