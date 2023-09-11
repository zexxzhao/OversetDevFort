!======================================================================
!
!======================================================================
subroutine SparseMatLoc_3D(list, n, goal, locat)
  implicit none
  integer, intent(in) :: n, list(n), goal
  integer, intent(out):: locat
  integer :: rowvl, rowvh, rowv

  ! Initialize
  rowvl = 1
  rowvh = n + 1

  ! do a binary search
  do
    if (rowvh - rowvl > 1) then
      rowv = (rowvh + rowvl)/2
      if (list(rowv) > goal) then
        rowvh = rowv
      else
        rowvl = rowv
      end if
    else
      exit
    end if
  end do
  locat = rowvl
end subroutine SparseMatLoc_3D

!======================================================================
!
!======================================================================
subroutine FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe)
  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl

  integer, intent(in) :: iel
  real(8), intent(in) :: xKebe11(NSD*NSD, NSHL, NSHL), &
                         xGebe(NSD, NSHL, NSHL), &
                         xDebe1(NSD, NSHL, NSHL), &
                         xMebe(NSHL, NSHL)
  integer :: a, b, c, d, ee, n, k, locn, i

  do a = 1, NSHL
    i = IEN(iel, a)
    c = col(i)
    n = col(i + 1) - c
    do b = 1, NSHL

      call SparseMatLoc_3D(row(c:c + n - 1), n, IEN(iel, b), locn)

      k = locn + c - 1

      LHSK11(1:NSD*NSD, k) = LHSK11(1:NSD*NSD, k) + xKebe11(1:NSD*NSD, a, b)

      LHSG(1:NSD, k) = LHSG(1:NSD, k) + xGebe(1:NSD, a, b)
      LHSD1(1:NSD, k) = LHSD1(1:NSD, k) + xDebe1(1:NSD, a, b)

      LHSM(k) = LHSM(k) + xMebe(a, b)

    end do
  end do

end subroutine FillSparseMat_3D

!======================================================================
!
!======================================================================
subroutine FillSparseMesh_3D(nshl, iel, xKebe22)

  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl

  integer, intent(in) :: iel
  real(8), intent(in) :: xKebe22(NSD*NSD, NSHL, NSHL)

  integer :: a, b, c, d, ee, n, k, locn, i

  do a = 1, NSHL
    i = IEN(iel, a)
    c = col(i)
    n = col(i + 1) - c
    do b = 1, NSHL

      call SparseMatLoc_3D(row(c:c + n - 1), n, IEN(iel, b), locn)

      k = locn + c - 1

      LHSK22(1:NSD*NSD, k) = LHSK22(1:NSD*NSD, k) + xKebe22(1:NSD*NSD, a, b)

    end do

  end do

end subroutine FillSparseMesh_3D

!======================================================================
!
!======================================================================
subroutine GetSparseMat_3D(nshl, iel, xKebe, xGebe, xDebet, xMebe)

  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl

  integer :: iel
  integer :: a, b, c, d, ee, n, k, locn, i

  real(8) :: xKebe(NSD*NSD, NSHL, NSHL), xGebe(NSD, NSHL, NSHL), &
             xDebet(NSD, NSHL, NSHL), xMebe(NSHL, NSHL)

  do a = 1, NSHL
    i = IEN(iel, a)
    c = col(i)
    n = col(i + 1) - c
    do b = 1, NSHL

      call SparseMatLoc_3D(row(c:c + n - 1), n, IEN(iel, b), locn)

      k = locn + c - 1

      xKebe(:, a, b) = LHSK11(:, k)
      xGebe(:, a, b) = LHSG(:, k)
      xDebet(:, a, b) = LHSD1(:, k)

      xMebe(a, b) = LHSM(k)

    end do
  end do

end subroutine GetSparseMat_3D

!======================================================================
!
!======================================================================
subroutine FillSparseMat_ls(nshl, iel, xMebe)

  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl

  integer, intent(in) :: iel
  real(8), intent(in) :: xMebe(NSHL, NSHL)

  integer :: a, b, c, d, ee, n, k, locn, i

  do a = 1, NSHL
    i = IEN(iel, a)
    c = col(i)
    n = col(i + 1) - c
    do b = 1, NSHL

      call SparseMatLoc_3D(row(c:c + n - 1), n, IEN(iel, b), locn)

      k = locn + c - 1
      LHSls(k) = LHSls(k) + xMebe(a, b)
    end do
  end do

end subroutine FillSparseMat_ls

!======================================================================
!
!======================================================================
subroutine FillSparseMat_tem(nshl, iel, xTebe)

  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl

  integer, intent(in) :: iel
  real(8), intent(in) :: xTebe(NSHL, NSHL)

  integer :: a, b, c, d, ee, n, k, locn, i

  do a = 1, NSHL
    i = IEN(iel, a)
    c = col(i)
    n = col(i + 1) - c
    do b = 1, NSHL

      call SparseMatLoc_3D(row(c:c + n - 1), n, IEN(iel, b), locn)

      k = locn + c - 1
      LHStem(k) = LHStem(k) + xTebe(a, b)
    end do
  end do

end subroutine FillSparseMat_tem
!======================================================================
!
!======================================================================
subroutine FillSparseMat_fm(nshl, iel, xMebe)

  use aAdjKeep
  use commonvars

  implicit none

  integer, intent(in) :: nshl

  integer, intent(in) :: iel
  real(8), intent(in) :: xMebe(NSHL, NSHL)

  integer :: a, b, c, d, ee, n, k, locn, i

  do a = 1, NSHL
    i = IEN(iel, a)
    c = col(i)
    n = col(i + 1) - c
    do b = 1, NSHL

      call SparseMatLoc_3D(row(c:c + n - 1), n, IEN(iel, b), locn)

      k = locn + c - 1
      LHSMass(k) = LHSMass(k) + xMebe(a, b)
    end do

  end do
end subroutine FillSparseMat_fm

!======================================================================
!
!======================================================================
subroutine FillSparseMat_NS_conv(nshl, iel, xKebe11, xGebe, xDebe1, &
                                 xDebe2, xMebe, xMebels, xLSUebe, &
                                 xPLSebe, xULSebe)
  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl

  integer, intent(in) :: iel
  real(8), intent(in) :: xKebe11(NSD*NSD, NSHL, NSHL), &
                         xGebe(NSD, NSHL, NSHL), &
                         xDebe1(NSD, NSHL, NSHL), &
                         xDebe2(NSD, NSHL, NSHL), &
                         xMebe(NSHL, NSHL), xMebels(NSHL, NSHL), &
                         xLSUebe(NSD, NSHL, NSHL), xPLSebe(NSHL, NSHL), &
                         xULSebe(NSD, NSHL, NSHL)

  integer :: a, b, c, d, ee, n, k, locn, i

  do a = 1, NSHL
    i = IEN(iel, a)
    c = col(i)
    n = col(i + 1) - c
    do b = 1, NSHL

      call SparseMatLoc_3D(row(c:c + n - 1), n, IEN(iel, b), locn)

      k = locn + c - 1

      LHSK11(1:NSD*NSD, k) = LHSK11(1:NSD*NSD, k) + xKebe11(1:NSD*NSD, a, b)

      LHSG(1:NSD, k) = LHSG(1:NSD, k) + xGebe(1:NSD, a, b)

      LHSD1(1:NSD, k) = LHSD1(1:NSD, k) + xDebe1(1:NSD, a, b)

      LHSD2(1:NSD, k) = LHSD2(1:NSD, k) + xDebe2(1:NSD, a, b)

      LHSM(k) = LHSM(k) + xMebe(a, b)

      LHSls(k) = LHSls(k) + xMebels(a, b)

      LHSlsu(1:NSD, k) = LHSlsu(1:NSD, k) + xLSUebe(1:NSD, a, b)

      LHSPls(k) = LHSPls(k) + xPLSebe(a, b)

      LHSUls(1:NSD, k) = LHSUls(1:NSD, k) + xULSebe(1:NSD, a, b)

    end do
  end do

end subroutine FillSparseMat_NS_conv
