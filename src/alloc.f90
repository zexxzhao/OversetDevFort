!======================================================================
!
!======================================================================
subroutine allocMatVec()

  use aAdjKeep
  use commonvars
  use defs_shell

  implicit none

  integer :: NGAUSSMAX, i

  ! Allocate solution vectors
  allocate (dg(NNODE, NSD), dgold(NNODE, NSD))

  allocate (ug(NNODE, NSD), ugold(NNODE, NSD), &
            ugm(NNODE, NSD), ugmold(NNODE, NSD))
  do i = 1, NRES
    allocate (bmesh(i)%ug(bmesh(i)%NNODE, NSD))
    allocate (bmesh(i)%dg(bmesh(i)%NNODE, NSD))
    ! allocate (bmesh(i)%owned(bmesh(i)%NNODE))
  enddo

  allocate (lhsgq(NNODE), rhsgq(NNODE))
  allocate (acg(NNODE, NSD), acgold(NNODE, NSD), &
            acgm(NNODE, NSD), acgmold(NNODE, NSD))

  allocate (phig(NNODE), phigold(NNODE), &
            rphig(NNODE), rphigold(NNODE))

  allocate (pg(NNODE), pgold(NNODE))

  allocate (rTg(NNODE), rTgold(NNODE))
  allocate (Tg(NNODE), Tgold(NNODE))

  ! allocate (FPKS(NELEM, NSD, NSD))

  ! allocate averaged solution vectors
  allocate (uavg(NNODE, NSD), pavg(NNODE))

  ! Allocate Residuals
  allocate (RHSGu(NNODE, NSD), RHSGm(NNODE, NSD), &
            RHSGp(NNODE), RHSGls(NNODE), RHSGTem(NNODE))

  ! Allocate Matrices
  allocate (lhsK11(NSD*NSD, icnt))
  allocate (lhsK12(NSD*NSD, icnt), lhsK22(NSD*NSD, icnt))

  allocate (lhsG(NSD, icnt), lhsD1(NSD, icnt), &
            lhsD2(NSD, icnt))

  allocate (LHSM(icnt), LHSLS(icnt), LHSmass(icnt))

  allocate (LHSlsu(NSD, icnt), LHSUls(NSD, icnt), LHSPLS(icnt))

  allocate (LHSTem(icnt))

  NGAUSSMAX = maxval(ELMNGAUSS)

  !--------------------------------------------
  ! allocation for shell
  !--------------------------------------------
  ! if (solshell) then
  !   ! allocate the shell RHS and LHS
! !  allocate(SH%RHSG(SH%TSP%NNODE,NSD), &
! !           SH%RHSG_EXT(SH%TSP%NNODE,NSD), &
! !           SH%RHSG_GRA(SH%TSP%NNODE,NSD))

  !   allocate (SH%RHSG(SH%NRB%NNODE, NSD), &
  !             SH%RHSG_EXT(SH%NRB%NNODE, NSD), &
  !             SH%RHSG_GRA(SH%NRB%NNODE, NSD))

  !   allocate (SH%LHSK(NSD*NSD, SH%icnt))
  !   SH%RHSG = 0.0d0; SH%RHSG_EXT = 0.0d0; SH%RHSG_GRA = 0.0d0
  !   SH%LHSK = 0.0d0

  !   ! Allocate solution variables for SHELL
  !   allocate (SH%FEM%dsh(SH%FEM%NNODE, NSD), SH%FEM%dshOld(SH%FEM%NNODE, NSD), &
  !             SH%FEM%ush(SH%FEM%NNODE, NSD), SH%FEM%ushOld(SH%FEM%NNODE, NSD), &
  !             SH%FEM%ash(SH%FEM%NNODE, NSD), SH%FEM%ashOld(SH%FEM%NNODE, NSD))
  !   SH%FEM%dsh = 0.0d0; SH%FEM%dshOld = 0.0d0
  !   SH%FEM%ush = 0.0d0; SH%FEM%ushOld = 0.0d0
  !   SH%FEM%ash = 0.0d0; SH%FEM%ashOld = 0.0d0

! !  allocate(SH%TSP%dsh(SH%TSP%NNODE,NSD), SH%TSP%dshOld(SH%TSP%NNODE,NSD), &
! !           SH%TSP%ush(SH%TSP%NNODE,NSD), SH%TSP%ushOld(SH%TSP%NNODE,NSD), &
! !           SH%TSP%ash(SH%TSP%NNODE,NSD), SH%TSP%ashOld(SH%TSP%NNODE,NSD))
! !  SH%TSP%dsh = 0.0d0; SH%TSP%dshOld = 0.0d0
! !  SH%TSP%ush = 0.0d0; SH%TSP%ushOld = 0.0d0
! !  SH%TSP%ash = 0.0d0; SH%TSP%ashOld = 0.0d0

  !   allocate (SH%NRB%dsh(SH%NRB%NNODE, NSD), SH%NRB%dshOld(SH%NRB%NNODE, NSD), &
  !             SH%NRB%ush(SH%NRB%NNODE, NSD), SH%NRB%ushOld(SH%NRB%NNODE, NSD), &
  !             SH%NRB%ash(SH%NRB%NNODE, NSD), SH%NRB%ashOld(SH%NRB%NNODE, NSD))
  !   SH%NRB%dsh = 0.0d0; SH%NRB%dshOld = 0.0d0
  !   SH%NRB%ush = 0.0d0; SH%NRB%ushOld = 0.0d0
  !   SH%NRB%ash = 0.0d0; SH%NRB%ashOld = 0.0d0
  ! end if

  ! !--------------------------------------------
  ! ! Allocate solution variables for Non-matching boundaries
  ! !--------------------------------------------
  ! if (nonmatch) then
  !   do i = 1, 2
  !     allocate (NM%FEM(i)%dsh(NM%FEM(i)%NNODE, NSD), NM%FEM(i)%dshOld(NM%FEM(i)%NNODE, NSD), &
  !               NM%FEM(i)%ush(NM%FEM(i)%NNODE, NSD), NM%FEM(i)%ushOld(NM%FEM(i)%NNODE, NSD), &
  !               NM%FEM(i)%ash(NM%FEM(i)%NNODE, NSD), NM%FEM(i)%ashOld(NM%FEM(i)%NNODE, NSD))
  !     NM%FEM(i)%dsh = 0.0d0; NM%FEM(i)%dshOld = 0.0d0
  !     NM%FEM(i)%ush = 0.0d0; NM%FEM(i)%ushOld = 0.0d0
  !     NM%FEM(i)%ash = 0.0d0; NM%FEM(i)%ashOld = 0.0d0
  !   end do
  ! end if
end subroutine allocMatVec

!======================================================================
!
!======================================================================
subroutine deallocMatVec

  use aAdjKeep
  use mpi
  use commonvars

  implicit none

  integer :: i

  ! deallocate mesh and flags
  deallocate (xg, wg, IEN)
  deallocate (ELM_ID)
  deallocate (IPER, IBC)
  ! deallocate (IS_SOLID_NODE)

  deallocate (EL_TYP, P_Flag, D_Flag)
  do i = 1, NRES
    deallocate (bmesh(i)%xg)
    deallocate (bmesh(i)%ug)
    deallocate (bmesh(i)%dg)
    deallocate (bmesh(i)%owned)
  enddo

  !if (NPATCH .gt. 0) deallocate (EPID, EIJK)
  !do i = 1, NPATCH
  !  deallocate (patch(i)%U_KNOT)
  !  deallocate (patch(i)%V_KNOT)
  !  deallocate (patch(i)%W_KNOT)
  !end do
  !deallocate (patch)

  do i = 1, NBOUND
    deallocate (bound(i)%FACE_IEN)
    deallocate (bound(i)%F2E)
    deallocate (bound(i)%FACE_OR)
    deallocate (bound(i)%BNODES)
  end do
  deallocate (bound)

  ! deallocate sparseStruct
  deallocate (col, row)

  ! dellocate solution vectors
  deallocate (dg, dgold)
  deallocate (ug, ugold, ugm, ugmold)
  deallocate (acg, acgold, acgm, acgmold)
  deallocate (phig, phigold, rphig, rphigold)
  deallocate (pg, pgold)
  deallocate (rTg, rTgold)
  deallocate (Tg, Tgold)
  deallocate (uavg, pavg)
  ! deallocate (FPKS)

  ! deallocate Residuals
  deallocate (RHSGu, RHSGm, RHSGp, RHSGls, RHSGtem)

  ! deallocate Matrices
  deallocate (lhsK11, lhsK12, lhsK22)
  deallocate (lhsG, lhsD1, lhsD2)
  deallocate (LHSM, LHSLS, LHSmass)
  deallocate (LHSlsu, LHSUls, lhsPls)
  deallocate (LHStem)

  ! Conditional deallocations
  if (numnodes .gt. 1) deallocate (ilworku)

end subroutine deallocMatVec
