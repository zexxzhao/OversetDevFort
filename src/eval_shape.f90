!======================================================================
! subroutine to evaluate shape or basis functions for volume
!======================================================================
subroutine eval_shape(nshl, iel, gp, xl, dl, wl, shlu, shgradgu, shhessg, &
                      dxidx, Gij, Ginv, hess_flag)
  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl

  integer, intent(in)  :: iel, hess_flag
  real(8), intent(in)  :: gp(NSD), xl(NSHL, NSD), dl(NSHL, NSD), wl(NSHL)
  real(8), intent(out) :: shlu(NSHL), shgradgu(NSHL, NSD), &
                          shhessg(NSHL, NSD, NSD), dxidx(NSD, NSD), &
                          Gij(NSD, NSD), Ginv(NSD, NSD)
  integer :: pn, ni, nj, nk

  if (iga) then
    ! pn = EPID(iel)
    ! ni = EIJK(iel, 1)
    ! nj = EIJK(iel, 2)
    ! nk = EIJK(iel, 3)

    ! call eval_shape_nurbs(nshl, gp, ni, nj, nk, patch(pn), xl + dl, wl, &
    !                       shlu, shgradgu, shhessg, dxidx, Gij, Ginv, &
    !                       hess_flag)
  else

    if (nshl == 4 .or. nshl == 6) then
      call eval_shape_fem(nshl, gp, xl, dl, shlu, shgradgu, dxidx, Gij, &
                          Ginv)

    else
      write (*, *) "ERROR: Undefined NSHL"
      stop
    end if
  end if

end subroutine eval_shape

!======================================================================
! subroutine to evaluate shape or basis functions for surface
!======================================================================
subroutine eval_faceshape(nshl, iel, gp, mgp, faceor, xl, dl, wl, shlu, &
                          shgradgu, dxidx, Gij, Ginv, nor)
  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl

  integer, intent(in)  :: iel
  real(8), intent(in)  :: gp(2), mgp(3), xl(NSHL, NSD), dl(NSHL, NSD), wl(NSHL)
  real(8), intent(out) :: shlu(NSHL), shgradgu(NSHL, NSD), &
                          dxidx(NSD, NSD), nor(NSD), Gij(NSD, NSD), &
                          Ginv(NSD, NSD)
  integer :: pn, ni, nj, nk, faceor

  if (iga) then
    ! pn = EPID(iel)
    ! ni = EIJK(iel, 1)
    ! nj = EIJK(iel, 2)
    ! nk = EIJK(iel, 3)
    ! call eval_faceshape_nurbs(nshl, gp, faceor, ni, nj, nk, patch(pn), &
    !                           xl + dl, wl, shlu, shgradgu, dxidx, &
    !                           Gij, Ginv, nor)

  else
    if (nshl == 4) then
      call eval_faceshape_tet(nshl, gp, mgp, faceor, xl, dl, shlu, shgradgu, &
                              dxidx, nor)

    else if (nshl == 6) then
      call eval_faceshape_pri(nshl, gp, mgp, faceor, xl, dl, shlu, shgradgu, &
                              dxidx, nor)
    end if
  end if

end subroutine eval_faceshape

!======================================================================
! subroutine to invert an 3x3 matrix
!======================================================================
subroutine get_inverse_3x3(Amat, Ainv, DetJ)
  implicit none
  real(8), intent(in)  :: Amat(3, 3)
  real(8), intent(out) :: Ainv(3, 3), DetJ
  real(8) :: tmp

  Ainv = 0.0d0

  Ainv(1, 1) = Amat(2, 2)*Amat(3, 3) - Amat(3, 2)*Amat(2, 3)
  Ainv(1, 2) = Amat(3, 2)*Amat(1, 3) - Amat(1, 2)*Amat(3, 3)
  Ainv(1, 3) = Amat(1, 2)*Amat(2, 3) - Amat(1, 3)*Amat(2, 2)

  tmp = 1.0d0/(Ainv(1, 1)*Amat(1, 1) + Ainv(1, 2)*Amat(2, 1) + &
               Ainv(1, 3)*Amat(3, 1))

  Ainv(1, 1) = Ainv(1, 1)*tmp
  Ainv(1, 2) = Ainv(1, 2)*tmp
  Ainv(1, 3) = Ainv(1, 3)*tmp

  Ainv(2, 1) = (Amat(2, 3)*Amat(3, 1) - Amat(2, 1)*Amat(3, 3))*tmp
  Ainv(2, 2) = (Amat(1, 1)*Amat(3, 3) - Amat(3, 1)*Amat(1, 3))*tmp
  Ainv(2, 3) = (Amat(2, 1)*Amat(1, 3) - Amat(1, 1)*Amat(2, 3))*tmp
  Ainv(3, 1) = (Amat(2, 1)*Amat(3, 2) - Amat(2, 2)*Amat(3, 1))*tmp
  Ainv(3, 2) = (Amat(3, 1)*Amat(1, 2) - Amat(1, 1)*Amat(3, 2))*tmp
  Ainv(3, 3) = (Amat(1, 1)*Amat(2, 2) - Amat(1, 2)*Amat(2, 1))*tmp

  DetJ = 1.0d0/tmp
end subroutine get_inverse_3x3
