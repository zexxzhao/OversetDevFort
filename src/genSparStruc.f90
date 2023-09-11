!======================================================================
!
!======================================================================
subroutine genSparStruc

  use aAdjKeep
  use commonvars

  implicit none

  integer :: tmpr(NNODE)
  integer :: adjcnt(NNODE), mloc(1)
  integer :: i, j, k, imin, ibig, ncol
  integer :: row_fill_list(NNODE, 27*maxNSHL)

  ! compute sparse matrix data structures
  row_fill_list = 0
  adjcnt = 0
  call Asadj(row_fill_list, adjcnt)

  ! allocate
  icnt = sum(adjcnt)

  allocate (col(NNODE + 1), row(icnt))

  ! build the colm array
  col(1) = 1
  do i = 1, NNODE
    col(i + 1) = col(i) + adjcnt(i)
  end do

  ! sort the rowp into increasing order
  ibig = 10*NNODE
  icnt = 0
  do i = 1, NNODE
    ncol = adjcnt(i)
    tmpr(1:ncol) = row_fill_list(i, 1:ncol)
    do j = 1, ncol
      icnt = icnt + 1
      imin = minval(tmpr(1:ncol))
      mloc = minloc(tmpr(1:ncol))
      row(icnt) = imin
      tmpr(mloc(1)) = ibig
    end do
  end do

end subroutine genSparStruc

!======================================================================
!
!======================================================================
subroutine Asadj(row_fill_list, adjcnt)

  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(out) :: row_fill_list(NNODE, 27*maxNSHL), adjcnt(NNODE)

  integer :: k, i, j, ni, nj, nk, l, ibroke, knd, nshl, &
             jnd, jlngth

  integer, allocatable :: ndlist(:)

  do i = 1, NELEM

    nshl = ELMNSHL(i)
    allocate (ndlist(nshl))

    do j = 1, nshl
      ndlist(j) = IEN(i, j)  ! gen list of global "nodes" for this element
    end do

    do j = 1, nshl
      jnd = ndlist(j)   ! jnd is the global "node" we are working on
      jlngth = adjcnt(jnd)  ! current length of j's list
      do k = 1, nshl
        knd = ndlist(k)
        ibroke = 0

        do l = 1, jlngth   ! row_fill_list is, for each node, the
          ! list of nodes that I have already
          ! detected interaction with
          if (row_fill_list(jnd, l) == knd) then
            ibroke = 1
            exit
          end if
        end do

        ! to get here k was not in  j's list so add it
        if (ibroke == 0) then
          jlngth = jlngth + 1 ! lengthen list
          if (jlngth .gt. 27*maxNSHL) write (*, *) "Stop the Music!!!"

          row_fill_list(jnd, jlngth) = knd ! add unique entry to list
        end if
      end do       ! finished checking all the k's for this j
      adjcnt(jnd) = jlngth  ! update the counter
    end do      ! done with j's

    deallocate (ndlist)

  end do         ! done with elements in this block

end subroutine Asadj
