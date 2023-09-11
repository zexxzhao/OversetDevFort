!======================================================================
!
!======================================================================
subroutine solveNavSto(inewt, MomRes0, ConRes0, NSConverged, lsres0, LSConverged)

  use aAdjKeep
  use mpi
  use commonvars

  implicit none

  integer, intent(in)    :: inewt
  real(8), intent(inout) :: MomRes0, ConRes0, LSRes0
  logical, intent(out)   :: NSConverged, LSConverged

  integer :: i, j, b, p

  real(8) :: dgAlpha(NNODE, NSD), &
             ugAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
             acgAlpha(NNODE, NSD), acgmAlpha(NNODE, NSD), &
             pgAlpha(NNODE), phigAlpha(NNODE), rphigAlpha(NNODE)

  real(8) :: MomRes, ConRes, LSRes, MomResl, ConResl, LSResl, gfor(NSD), gmom(NSD)
  real(8) :: tmp1(NNODE, NSD), vtmp(NSD), TRACT(NNODE, NSD)
  real(8) :: t1, t2
  character(len=80) :: fomt

  call setBCs_CFD()
  acgAlpha = acgold + almi*(acg - acgold)
  acgmAlpha = acgmold + almi*(acgm - acgmold)
  ugAlpha = ugold + alfi*(ug - ugold)
  ugmAlpha = ugmold + alfi*(ugm - ugmold)
  dgAlpha = dgold + alfi*(dg - dgold)
  pgAlpha = pg

  phigAlpha = phigold + alfi*(phig - phigold)
  rphigAlpha = rphigold + almi*(rphig - rphigold)

  ! Assemble Navier-Stokes LHS/RHS
  RHSGu = 0.0d0
  RHSGp = 0.0d0
  RHSGls = 0.0d0

  LHSK11 = 0.0d0
  LHSG = 0.0d0
  LHSD1 = 0.0d0
  LHSM = 0.0d0

  LHSLS = 0.0d0
  LHSULS = 0.0d0
  LHSLSU = 0.0d0
  LHSPLS = 0.0d0
  lhsgq = 0.0d0
  rhsgq = 0.0d0
!!!  call FaceAssembly_NS_conv(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!!!                            acgmAlpha, pgAlpha, phigAlpha)

!!!  call FaceAssembly_NS(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!!!                       acgmAlpha, pgAlpha, phigAlpha)
  if (myid .eq. 0) then
    call CPU_TIME(t1)
  end if
  call IntElmAss_3D(inewt, dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                    acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, 0)

  call FaceAssembly_NS_weak(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                            acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, 0)
  if (myid .eq. 0) then
    call CPU_TIME(t2)
    write (*, *) "Total time assemble:", t2 - t1, "seconds"
  end if
  ! Compute residuals
  if (numnodes > 1) then
    call commu(RHSGp, 1, 'in ')
    call commu(RHSGu, NSD, 'in ')
    call commu(RHSGls, 1, 'in ')
    call commu(rhsgq, 1, 'in ')
    call commu(rhsgq, 1, 'out ')
    call commu(lhsgq, 1, 'in ')
    call commu(lhsgq, 1, 'out ')
  end if

  MomResl = sum(RHSGu(:, 1)**2) + sum(RHSGu(:, 2)**2) + sum(RHSGu(:, 3)**2)
  ConResl = sum(RHSGp**2)
  LSResl = sum(RHSGls**2)

  MomRes = MomResl
  ConRes = ConResl
  LSRes = LSResl

  if (numnodes > 1) then
    call MPI_ALLREDUCE(MomResl, MomRes, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
    call MPI_ALLREDUCE(ConResl, ConRes, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
    call MPI_ALLREDUCE(LSResl, LSRes, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
  end if

  MomRes = sqrt(MomRes)
  ConRes = sqrt(ConRes)
  LSRes = sqrt(LSRes)

  ! Print residual
  if (inewt == 1) then
    MomRes0 = MomRes
    ConRes0 = ConRes
    LSRes0 = LSRes
  end if

  if (ismaster) then
    fomt = "(I3,a,x,ES13.6,x,F12.6)"
    write (*, fomt) inewt, ") Total Mom. Res. Norm = ", &
      MomRes, 1.0d2*MomRes/MomRes0

    write (*, fomt) inewt, ") Continuity Res. Norm = ", &
      ConRes, 1.0d2*ConRes/ConRes0

    write (*, fomt) inewt, ") LS Res. Norm = ", &
      LSRes, 1.0d2*LSRes/LSRes0
    write (*, *)
  end if

  if (isnan(MomRes + ConRes + LSRes)) then
    write (*, *) myid, inewt, '!=== Norm is NaN ===!'
    stop
  end if

  ! Check convergence and solve linear system
  if (((MomRes/MomRes0) < NS_NL_Utol) .and. &
      ((ConRes/ConRes0) < NS_NL_Ptol) .and. &
      ((LSRes/LSRes0) < LSC_NL_tol)) then
    NSConverged = .true.
    if (ismaster) then
      write (*, *)
      write (*, *) "    Converged: skip linear solve"
      write (*, *)
    end if

  else
    NSConverged = .false.
    acgAlpha = 0.0d0; pgAlpha = 0.0d0; rphigAlpha = 0.0d0
!    lhsLSu = 0.0d0!; lhsUls = 0.0d0; lhsPls = 0.0d0
    if (myid .eq. 0) then
      call CPU_TIME(t1)
    end if
    call SparseGMRES_up(col, row, IBC, IPER, D_FLAG, P_FLAG, &
                        rhsGu, rhsGp, acgAlpha, pgAlpha, &
                        lhsK11, lhsG, lhsD1, lhsM, icnt, &
                        NS_GMRES_tol, NS_GMRES_itermax, &
                        NS_GMRES_itermin, &
                        NNODE, maxNSHL, NSD, &
                        rphigAlpha, lhsLS, &
                        lhsLSu, lhsUls, lhsPls, rhsGls)

    if (myid .eq. 0) then
      call CPU_TIME(t2)
      write (*, *) "Total time solve:", t2 - t1, "seconds"
    end if

    acg = acg + acgAlpha
    ug = ug + gami*Delt*acgAlpha
    pg = pg + alfi*gami*Delt*pgAlpha

    rphig = rphig + rphigAlpha
    phig = phig + gami*Delt*rphigAlpha
  end if

  !--------------------------------------------------------------------
  ! compute torque (weak BC/Surface Int.)
  !--------------------------------------------------------------------
  acgAlpha = acgold + almi*(acg - acgold)
  acgmAlpha = acgmold + almi*(acgm - acgmold)
  ugAlpha = ugold + alfi*(ug - ugold)
  ugmAlpha = ugmold + alfi*(ugm - ugmold)
  dgAlpha = dgold + alfi*(dg - dgold)
  pgAlpha = pg
  phigAlpha = phigold + alfi*(phig - phigold)
  rphigAlpha = rphigold + alfi*(rphig - rphigold)

!!$  ! compute RHS and the torque for individual blade
!!$  call compute_RHS_PhVI(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!!$                        acgmAlpha, pgAlpha, phigAlpha)

  ! compute RHS and the torque (RHSGu is used afterwards for FSI)
!  call compute_RHS(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!                   acgmAlpha, pgAlpha, phigAlpha, 4,     &  ! 1-nrel; 4-5mw
!                   torque3, torque4, force_trac, force_cons)

!!$!!$  ! output the force to compare with shell
!!$!!$  RHSGm = 0.0d0
!!$!!$  do b = 1, NBOUND
!!$!!$    if (bound(b)%Face_ID == 4) then
!!$!!$      do i = 1, bound(b)%NNODE
!!$!!$        j = bound(b)%BNODES(i)
!!$!!$        RHSGm(j,:) = RHSGu(j,:)
!!$!!$      end do
!!$!!$    end if
!!$!!$  end do
!!$!!$
!!$!!$  MomResl = sum(RHSGm(:,1)**2) + sum(RHSGm(:,2)**2) + sum(RHSGm(:,3)**2)
!!$!!$  MomRes  = MomResl
!!$!!$
!!$!!$  if (numnodes > 1) then
!!$!!$    call MPI_ALLREDUCE(MomResl, MomRes, 1, MPI_DOUBLE_PRECISION, &
!!$!!$                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
!!$!!$  end if
!!$!!$
!!$!!$  MomRes = sqrt(MomRes)
!!$!!$
!!$!!$  if (ismaster) then
!!$!!$    write(*,*) "*) Total Force Norm = ", MomRes
!!$!!$  end if
!!$
!!$  !-- compute conservative force and torque (Strong BC/Volume Int.) ---
!!$  if (BCugType(4,1) == 3 .or. BCugType(5,1) == 3 .or. &
!!$      BCugType(6,1) == 3 .or. BCugType(7,1) == 3) then
!!$
!!$    RHSGu = 0.0d0
!!$    call IntElmAss_3D(inewt, dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!!$                      acgmAlpha, pgAlpha, phigAlpha, 1)
!!$
!!$    if (numnodes > 1) call commu(RHSGu, NSD, 'in ')
!!$
!!$    ! only leave vaules on the boundary
!!$    RHSGm = 0.0d0
!!$    do b = 1, NBOUND
!!$      if (bound(b)%Face_ID == 21 .or. bound(b)%Face_ID == 22 .or. &
!!$          bound(b)%Face_ID == 23 .or. bound(b)%Face_ID == 24) then
!!$        do i = 1, bound(b)%NNODE
!!$          j = bound(b)%BNODES(i)
!!$          RHSGm(j,:) = RHSGu(j,:)
!!$        end do
!!$      end if
!!$    end do
!!$
!!$    call compute_moment3(NNODE, NSD, RHSGm, xg, dgAlpha, torque2)
!!$  end if
!!$  !-- end -----------------------------------------------------------

  if (numnodes > 1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

end subroutine solveNavSto

