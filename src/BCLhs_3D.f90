!======================================================================
!
!======================================================================
subroutine BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, &
                    Rhsu, Rhsp)
  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl
  real(8), intent(inout) :: xKebe11(NSD*NSD, NSHL, NSHL), &
                            xGebe(NSD, NSHL, NSHL), &
                            xDebe1(NSD, NSHL, NSHL), &
                            xMebe(NSHL, NSHL), &
                            Rhsu(NSD, NSHL), Rhsp(NSHL)
  integer :: iel, aa, bb, cc

  ! loop through local nodes
  do aa = 1, NSHL
    ! global node number of local node aa
    cc = IEN(iel, aa)

    ! Physics constraint in x direction
    if (IBC(cc, 1) == 1) then

      ! update K-matrix to account for constraint
      ! project away row (X-dir)
      do bb = 1, NSHL
        xKebe11(1, aa, bb) = 0.0d0 ! Phys-Phys
        xKebe11(2, aa, bb) = 0.0d0
        xKebe11(3, aa, bb) = 0.0d0

        xGebe(1, aa, bb) = 0.0d0
      end do

      !project away column (X-dir)
      do bb = 1, NSHL
        xKebe11(1, bb, aa) = 0.0d0 ! Phys-Phys
        xKebe11(4, bb, aa) = 0.0d0
        xKebe11(7, bb, aa) = 0.0d0

        xDebe1(1, bb, aa) = 0.0d0
      end do

      xKebe11(1, aa, aa) = 1.0d0

      Rhsu(1, aa) = 0.0d0   ! rhs
    end if

    ! Physics constraint in y direction
    if (IBC(cc, 2) == 1) then

      ! update K-matrix to account for constraint
      ! project away row (X-dir)
      do bb = 1, NSHL
        xKebe11(4, aa, bb) = 0.0d0 ! Phys-Phys
        xKebe11(5, aa, bb) = 0.0d0
        xKebe11(6, aa, bb) = 0.0d0

        xGebe(2, aa, bb) = 0.0d0
      end do

      ! project away column (X-dir)
      do bb = 1, NSHL
        xKebe11(2, bb, aa) = 0.0d0 ! Phys-Phys
        xKebe11(5, bb, aa) = 0.0d0
        xKebe11(8, bb, aa) = 0.0d0

        xDebe1(2, bb, aa) = 0.0d0
      end do

      xKebe11(5, aa, aa) = 1.0d0

      Rhsu(2, aa) = 0.0d0   ! rhs
    end if

    ! Physics constraint in z direction
    if (IBC(cc, 3) == 1) then

      ! update K-matrix to account for constraint
      ! project away row (X-dir)
      do bb = 1, NSHL
        xKebe11(7, aa, bb) = 0.0d0 ! Phys-Phys
        xKebe11(8, aa, bb) = 0.0d0
        xKebe11(9, aa, bb) = 0.0d0

        xGebe(3, aa, bb) = 0.0d0
      end do

      !project away column (X-dir)
      do bb = 1, NSHL
        xKebe11(3, bb, aa) = 0.0d0 ! Phys-Phys
        xKebe11(6, bb, aa) = 0.0d0
        xKebe11(9, bb, aa) = 0.0d0

        xDebe1(3, bb, aa) = 0.0d0
      end do

      xKebe11(9, aa, aa) = 1.0d0

      Rhsu(3, aa) = 0.0d0   ! rhs
    end if

    ! Pressure Constraint
    if (IBC(cc, 4) == 1) then

      ! update K-matrix to account for constraint
      ! project away row (X-dir)
      do bb = 1, NSHL
        xMebe(aa, bb) = 0.0d0 ! Phys-Phys
        xDebe1(:, aa, bb) = 0.0d0
        xGebe(:, bb, aa) = 0.0d0
      end do

      xMebe(aa, aa) = 1.0d0

      Rhsp(aa) = 0.0d0     ! rhs
    end if

  end do

end subroutine BCLhs_3D

!======================================================================
!
!======================================================================
subroutine BCMesh_3D(nshl, iel, xKebe22, Rhsm)

  use aAdjKeep
  use commonvars
  implicit none

  integer :: iel, aa, bb, cc
  integer, intent(in)    :: nshl
  real(8) :: xKebe22(NSD*NSD, NSHL, NSHL)
  real(8) :: Rhsm(NSD, NSHL)

  ! loop through local nodes
  do aa = 1, NSHL
    ! global node number of local node aa
    cc = IEN(iel, aa)

    ! Mesh constraint in x direction
    if (IBC(cc, 4) == 1) then
      ! project away row (X-dir)
      do bb = 1, NSHL
        xKebe22(1, aa, bb) = 0.0d0 ! Phys-Phys
        xKebe22(2, aa, bb) = 0.0d0
        xKebe22(3, aa, bb) = 0.0d0
      end do

      ! project away column (X-dir)
      do bb = 1, NSHL
        xKebe22(1, bb, aa) = 0.0d0 ! Phys-Phys
        xKebe22(4, bb, aa) = 0.0d0
        xKebe22(7, bb, aa) = 0.0d0
      end do

      xKebe22(1, aa, aa) = 1.0d0

      Rhsm(1, aa) = 0.0d0   ! rhs
    end if

    ! Mesh constraint in y direction
    if (IBC(cc, 5) == 1) then

      ! project away row (X-dir)
      do bb = 1, NSHL
        xKebe22(4, aa, bb) = 0.0d0 ! Phys-Phys
        xKebe22(5, aa, bb) = 0.0d0
        xKebe22(6, aa, bb) = 0.0d0
      end do

      ! project away column (X-dir)
      do bb = 1, NSHL
        xKebe22(2, bb, aa) = 0.0d0 ! Phys-Phys
        xKebe22(5, bb, aa) = 0.0d0
        xKebe22(8, bb, aa) = 0.0d0
      end do

      xKebe22(5, aa, aa) = 1.0d0

      Rhsm(2, aa) = 0.0d0   ! rhs

    end if

    ! Mesh constraint in z direction
    if (IBC(cc, 6) == 1) then
      ! project away row (X-dir)
      do bb = 1, NSHL
        xKebe22(7, aa, bb) = 0.0d0 ! Phys-Phys
        xKebe22(8, aa, bb) = 0.0d0
        xKebe22(9, aa, bb) = 0.0d0
      end do

      ! project away column (X-dir)
      do bb = 1, NSHL
        xKebe22(3, bb, aa) = 0.0d0 ! Phys-Phys
        xKebe22(6, bb, aa) = 0.0d0
        xKebe22(9, bb, aa) = 0.0d0
      end do

      xKebe22(9, aa, aa) = 1.0d0

      Rhsm(3, aa) = 0.0d0   ! rhs

    end if

  end do

end subroutine BCMesh_3D

!======================================================================
!
!======================================================================
subroutine BCLhs_ls(nshl, iel, xMebe, Rhs)

  use aAdjKeep
  use commonvars
  implicit none

  integer :: iel, aa, bb, cc
  integer, intent(in) :: nshl
  real(8) :: xMebe(NSHL, NSHL), Rhs(NSHL)

  ! loop through local nodes
  do aa = 1, NSHL

    cc = IEN(iel, aa)  ! global node number of local node aa

    if (IBC(cc, 8) == 1) then

      xMebe(aa, :) = 0.0d0
      xMebe(:, aa) = 0.0d0

      xMebe(aa, aa) = 1.0d0
      Rhs(aa) = 0.0d0    ! rhs

    end if
  end do

end subroutine BCLhs_ls

!======================================================================
!
!======================================================================
subroutine BCLhs_tem(nshl, iel, xTebe, Rhs)

  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl, iel
  real(8), intent(inout) :: xTebe(NSHL, NSHL), Rhs(NSHL)

  integer :: aa, bb, cc

  ! loop through local nodes
  do aa = 1, NSHL

    cc = IEN(iel, aa)  ! global node number of local node aa

    if (IBC(cc, 6) == 1) then

      xTebe(aa, :) = 0.0d0
      xTebe(:, aa) = 0.0d0

      xTebe(aa, aa) = 1.0d0
      Rhs(aa) = 0.0d0    ! rhs

    end if
  end do

end subroutine BCLhs_tem
!======================================================================
!
!======================================================================
subroutine BCLhs_NS_conv(nshl, iel, xKebe11, xGebe, xDebe1, &
                         xDebe2, xMebe, xMebels, xLSUebe, xPLSebe, &
                         xULSebe, Rhsu, Rhsp, Rhsls)
  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in)    :: nshl
  integer :: iel, aa, bb, cc

  real(8) :: xKebe11(NSD*NSD, NSHL, NSHL), &
             xGebe(NSD, NSHL, NSHL), &
             xDebe1(NSD, NSHL, NSHL), &
             xDebe2(NSD, NSHL, NSHL), &
             xMebe(NSHL, NSHL), &
             xMebels(NSHL, NSHL), &
             xLSUebe(NSD, NSHL, NSHL), &
             xPLSebe(NSHL, NSHL), &
             xULSebe(NSD, NSHL, NSHL)
  real(8) :: Rhsu(NSD, NSHL), Rhsp(NSHL), Rhsls(NSHL)

  ! loop through local nodes
  do aa = 1, NSHL
    ! global node number of local node aa
    cc = IEN(iel, aa)

    ! Physics constraint in x direction
    if (IBC(cc, 1) == 1) then

      ! update K-matrix to account for constraint
      ! project away row (X-dir)
      do bb = 1, NSHL
        xKebe11(1, aa, bb) = 0.0d0 ! Phys-Phys
        xKebe11(2, aa, bb) = 0.0d0
        xKebe11(3, aa, bb) = 0.0d0

        xGebe(1, aa, bb) = 0.0d0
        xDebe2(1, aa, bb) = 0.0d0
        xLSUebe(1, aa, bb) = 0.0d0
      end do

      !project away column (X-dir)
      do bb = 1, NSHL
        xKebe11(1, bb, aa) = 0.0d0 ! Phys-Phys
        xKebe11(4, bb, aa) = 0.0d0
        xKebe11(7, bb, aa) = 0.0d0

        xDebe1(1, bb, aa) = 0.0d0
        xLSUebe(1, bb, aa) = 0.0d0
      end do

      xKebe11(1, aa, aa) = 1.0d0

      Rhsu(1, aa) = 0.0d0   ! rhs
    end if

    ! Physics constraint in y direction
    if (IBC(cc, 2) == 1) then

      ! update K-matrix to account for constraint
      ! project away row (X-dir)
      do bb = 1, NSHL
        xKebe11(4, aa, bb) = 0.0d0 ! Phys-Phys
        xKebe11(5, aa, bb) = 0.0d0
        xKebe11(6, aa, bb) = 0.0d0

        xGebe(2, aa, bb) = 0.0d0
        xDebe2(2, aa, bb) = 0.0d0
        xLSUebe(2, aa, bb) = 0.0d0
      end do

      ! project away column (X-dir)
      do bb = 1, NSHL
        xKebe11(2, bb, aa) = 0.0d0 ! Phys-Phys
        xKebe11(5, bb, aa) = 0.0d0
        xKebe11(8, bb, aa) = 0.0d0

        xDebe1(2, bb, aa) = 0.0d0
        xLSUebe(2, bb, aa) = 0.0d0
      end do

      xKebe11(5, aa, aa) = 1.0d0

      Rhsu(2, aa) = 0.0d0   ! rhs
    end if

    ! Physics constraint in z direction
    if (IBC(cc, 3) == 1) then

      ! update K-matrix to account for constraint
      ! project away row (X-dir)
      do bb = 1, NSHL
        xKebe11(7, aa, bb) = 0.0d0 ! Phys-Phys
        xKebe11(8, aa, bb) = 0.0d0
        xKebe11(9, aa, bb) = 0.0d0

        xGebe(3, aa, bb) = 0.0d0
        xDebe2(3, aa, bb) = 0.0d0
        xLSUebe(3, aa, bb) = 0.0d0
      end do

      !project away column (X-dir)
      do bb = 1, NSHL
        xKebe11(3, bb, aa) = 0.0d0 ! Phys-Phys
        xKebe11(6, bb, aa) = 0.0d0
        xKebe11(9, bb, aa) = 0.0d0

        xDebe1(3, bb, aa) = 0.0d0
        xLSUebe(3, bb, aa) = 0.0d0
      end do

      xKebe11(9, aa, aa) = 1.0d0

      Rhsu(3, aa) = 0.0d0   ! rhs
    end if

    ! Pressure Constraint
    if (IBC(cc, 7) == 1) then

      ! update K-matrix to account for constraint
      ! project away row (X-dir)
      do bb = 1, NSHL
        xMebe(aa, bb) = 0.0d0 ! Phys-Phys
        xDebe1(:, aa, bb) = 0.0d0
        xGebe(:, bb, aa) = 0.0d0

        xPLSebe(aa, :) = 0.0d0
      end do

      xMebe(aa, aa) = 1.0d0

      Rhsp(aa) = 0.0d0     ! rhs
    end if

    ! Levelset Constraint
    if (IBC(cc, 8) .eq. 1) then

      xMebels(aa, :) = 0.0d0
      xMebels(:, aa) = 0.0d0
      xLSUebe(:, aa, :) = 0.0d0
      xPLSebe(:, aa) = 0.0d0
      xULSebe(:, :, aa) = 0.0d0

      xMebels(aa, aa) = 1.0d0
      Rhsls(aa) = 0.0d0    ! rhs

    end if

  end do

end subroutine BCLhs_NS_conv
