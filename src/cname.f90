!======================================================================
!
!======================================================================
function cname(i)

  implicit none

  character(len=10) :: cname
  integer, intent(in) :: i

  integer :: il(0:8), ic0, ii, k
  logical :: beg
  character(len=10) :: cc

  ic0 = ICHAR("0")
  cc = " "
  ii = i

  il(0) = mod(ii, 10)
  do k = 1, 8
    ii = (ii - il(k - 1))/10
    il(k) = mod(ii, 10)
  end do

  beg = .false.

  do k = 8, 1, -1
    if (il(k) .ne. 0 .or. beg) then
      beg = .true.
      cc = TRIM(cc)//CHAR(ic0 + il(k))
    end if
  end do

  cc = TRIM(cc)//CHAR(ic0 + il(0))
  cname = "."//cc

end function cname

!======================================================================
!
!======================================================================
function cname2(i)

  implicit none

  character(len=10) :: cname2
  integer, intent(in) :: i

  integer :: ic0, ii, k, il(0:8)
  logical :: beg
  character(len=10) :: cc

  ic0 = ICHAR("0")
  cc = " "
  ii = i

  il(0) = mod(ii, 10)
  do k = 1, 8
    ii = (ii - il(k - 1))/10
    il(k) = mod(ii, 10)
  end do

  beg = .false.

  do k = 8, 1, -1
    if (il(k) .ne. 0 .or. beg) then
      beg = .true.
      cc = TRIM(cc)//CHAR(ic0 + il(k))
    end if
  end do

  cname2 = TRIM(cc)//CHAR(ic0 + il(0))

end function cname2
