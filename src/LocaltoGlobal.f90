!======================================================================
!
!======================================================================
subroutine LocaltoGlobal_3D(nshl, iel, Rhsu, Rhsp)

  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl

  integer, intent(in) :: iel
  real(8), intent(in) :: Rhsu(NSD, NSHL), Rhsp(NSHL)
  integer :: bb

  do bb = 1, NSHL
    RHSGu(IEN(iel, bb), :) = RHSGu(IEN(iel, bb), :) + Rhsu(:, bb)
    RHSGp(IEN(iel, bb)) = RHSGp(IEN(iel, bb)) + Rhsp(bb)
  end do

end subroutine LocaltoGlobal_3D

!======================================================================
!
!======================================================================
subroutine LocaltoGlobalNSVOF_3D(RHSGu, RHSGp, RHSGls, &
                                 NNODE, NSD, &
                                 NSHL, maxNSHL, &
                                 ix, &
                                 Rhsu, Rhsp, Rhsls)

  !use aAdjKeep
  !use commonvars
  implicit none

  integer, intent(in) :: NNODE, NSD, NSHL, maxNSHL
  real(8), intent(in) :: Rhsu(NSD, NSHL), Rhsp(NSHL), Rhsls(NSHL)

  integer, intent(in) :: ix(maxNSHL)

  real(8), intent(inout) :: RHSGu(NNODE, NSD), RHSGp(NNODE), RHSGls(NNODE)
  integer :: bb, gix

  do bb = 1, NSHL
    gix = ix(bb)
    RHSGu(gix, :) = RHSGu(gix, :) + Rhsu(:, bb)
    RHSGp(gix) = RHSGp(gix) + Rhsp(bb)
    RHSGls(gix) = RHSGls(gix) + Rhsls(bb)
  end do

end subroutine LocaltoGlobalNSVOF_3D

!======================================================================
!
!======================================================================

subroutine LocaltoGlobal_3D_mesh(nshl, iel, Rhsm)

  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl

  integer, intent(in) :: iel
  real(8), intent(in) :: Rhsm(NSD, NSHL)
  integer :: bb

  do bb = 1, NSHL
    RHSGm(IEN(iel, bb), :) = RHSGm(IEN(iel, bb), :) + Rhsm(:, bb)
  end do

end subroutine LocaltoGlobal_3D_mesh

!======================================================================
!
!======================================================================
subroutine LocaltoGlobal_ls(nshl, iel, Rhs)

  use aAdjKeep
  use commonvars

  implicit none

  integer, intent(in) :: nshl

  integer, intent(in) :: iel
  real(8), intent(in) :: Rhs(NSHL)
  integer :: bb

  do bb = 1, NSHL
    RHSGls(IEN(iel, bb)) = RHSGls(IEN(iel, bb)) + Rhs(bb)
  end do

end subroutine LocaltoGlobal_ls

!======================================================================
!
!======================================================================
subroutine LocaltoGlobal_NS_conv(nshl, iel, Rhsu, Rhsp, Rhsls)

  use aAdjKeep
  use commonvars

  implicit none

  integer, intent(in) :: nshl

  integer, intent(in) :: iel
  real(8), intent(in) :: Rhsu(NSD, NSHL), Rhsp(NSHL), Rhsls(NSHL)
  integer :: bb

  do bb = 1, NSHL
    RHSGu(IEN(iel, bb), :) = RHSGu(IEN(iel, bb), :) + Rhsu(:, bb)
    RHSGp(IEN(iel, bb)) = RHSGp(IEN(iel, bb)) + Rhsp(bb)
    RHSGls(IEN(iel, bb)) = RHSGls(IEN(iel, bb)) + Rhsls(bb)
  end do

end subroutine LocaltoGlobal_NS_conv

!======================================================================
!
!======================================================================
subroutine LocaltoGlobal(nshl, iel, Rhsloc, RHSglob)

  use aAdjKeep
  use commonvars

  implicit none

  integer, intent(in) :: nshl

  integer :: iel
  integer :: aa, bb
  real(8) :: Rhsloc(NSHL), RHSglob(NNODE)

  do bb = 1, NSHL
    RHSglob(IEN(iel, bb)) = RHSglob(IEN(iel, bb)) + Rhsloc(bb)
  end do

end subroutine LocaltoGlobal

