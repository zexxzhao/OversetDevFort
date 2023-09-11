!======================================================================
!
!======================================================================
subroutine ctypes

  use mpi
  use commonvars

  implicit none

  integer :: sizeofdouble, is, itask, workf, kdof
  integer :: isbegin(maxseg), lenseg(maxseg), ioffset(maxseg)
  character(len=30) :: fname
  character(len=10) :: cname

  !-- READ LOCAL WORK FOR U AND P --------------------

  workf = 17
  fname = trim('lwork'//cname(myid + 1))//'.dat'
  open (workf, file=fname, status='old')

  read (workf, *) numtask, nlworku

  allocate (ilworku(nlworku))
  ilworku = 0
  ilworku(1) = numtask

  itkbeg = 1

  do itask = 1, numtask
    read (workf, *) ilworku(itkbeg + 1) ! tag
    read (workf, *) ilworku(itkbeg + 2) ! iacc - role in communication
    read (workf, *) ilworku(itkbeg + 3) ! partner node
    read (workf, *) ilworku(itkbeg + 4) ! number of segments
    do is = 1, ilworku(itkbeg + 4)
      read (workf, *) ilworku(itkbeg + 3 + 2*is) ! isgbegin
      read (workf, *) ilworku(itkbeg + 4 + 2*is) ! lenseg
    end do
    itkbeg = itkbeg + 4 + 2*ilworku(itkbeg + 4)
  end do

  close (workf)

  !-- DONE READING LOCAL WORK FOR D, U and P --------------

  call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, sizeofdouble, mpi_err)
  lstrideu = NNODE*sizeofdouble

  ! maxfront is a common variable being set in this routine
  maxfrontu = 0
  numtask = ilworku(1)
  itkbeg = 1

  if (numtask > maxtask) then
    write (*, *) "ERROR: numtask > maxtask", numtask, maxtask
    call MPI_abort(MPI_COMM_WORLD, 911, mpi_err)
  end if

  do itask = 1, numtask

    ! iacc = 0 ==> this task is a send (slave)
    !  = 1 ==> this task is a recieve (master)
    iacc = ilworku(itkbeg + 2)

    ! numseg : number of data segments to be communicated
    numseg = ilworku(itkbeg + 4)

    if (numseg > maxseg) then
      write (*, *) "ERROR: numseg > maxseg at task ", itask
      call MPI_abort(MPI_COMM_WORLD, 911, mpi_err)
    end if

    ! adjust the number of the other processor, since processors
    ! are numbered here starting from 0, not 1.
    ilworku(itkbeg + 3) = ilworku(itkbeg + 3) - 1

    ! lfront = total number of nodes involved in this task
    isbegin = 0  ! Initialize
    lenseg = 0
    lfrontu = 0
    do is = 1, numseg

      ! isbegin(is): starting node number for each segment
      isbegin(is) = ilworku(itkbeg + 3 + 2*is)

      ! lenseg(is): length of each segment (number of nodes)
      lenseg(is) = ilworku(itkbeg + 4 + 2*is)

      ! increment the total node counter
      lfrontu = lfrontu + lenseg(is)

    end do

    ! maxfront: number of nodes which will be communicated, including
    !   all segments. Note that after the loop over tasks
    !   is complete, maxfront will contain the maximum number
    !   of nodes for any of the tasks.
    maxfrontu = max(maxfrontu, lfrontu)

    ! ioffset: array offset from the first node in the first segment
    ioffset = 0     ! Initialize
    do is = 1, numseg
      ioffset(is) = isbegin(is) - isbegin(1)
    end do

    !-- now set up the MPI data types which will be used in commu.f.
    ! These data types represent the indexed sets that will be sent
    ! and recieved.
    !
    ! the following call to MPI_TYPE_INDEXED will create a new data
    ! type which will represent the blocks of data we wish to transfer
    ! for this task. A handle to the new type is returned
    ! (sevsegtype(itask,1)). This data type describes the blocks of
    ! data to be transferred in terms of segments.
    ! Input to this routine:
    !     numseg: number of segments in this task
    !     lenseg: length of each segment (number of nodes)
    !     ioffset: where to begin each block with respect to the
    !      first segment
    !     MPI_DOUBLE_PRECISION: type to set for each of the blocks --
    call MPI_TYPE_INDEXED(numseg, lenseg, ioffset, &
                          MPI_DOUBLE_PRECISION, sevsegtypeu(itask, 1), &
                          mpi_err)

    ! now create a new data type for each of the types of arrays we
    ! may wish to communicate with. For example ndof will be used when
    ! communicating the residual vector. Each one of these is derived
    ! from the first data type defined above, sevsegtype(itask,1).
    call MPI_TYPE_HVECTOR(NSD, 1, lstrideu, sevsegtypeu(itask, 1), &
                          sevsegtypeu(itask, 2), mpi_err)

    call MPI_TYPE_HVECTOR(NSD*NSD, 1, lstrideu, sevsegtypeu(itask, 1), &
                          sevsegtypeu(itask, 3), mpi_err)

    call MPI_TYPE_HVECTOR((NSD + 1)*(NSD + 1), 1, lstrideu, &
                          sevsegtypeu(itask, 1), &
                          sevsegtypeu(itask, 4), mpi_err)

    call MPI_TYPE_INDEXED(numseg, lenseg, ioffset, &
                          MPI_INTEGER, sevsegtypeu(itask, 5), mpi_err)

    ! now this must be done to make MPI recognize each of the data
    ! types that were just defined
    do kdof = 1, 5
      call MPI_TYPE_COMMIT(sevsegtypeu(itask, kdof), mpi_err)
    end do

    ! set the counter to the index in ilwork where the next task
    ! begins
    itkbeg = itkbeg + 4 + 2*numseg

  end do    !end loop over tasks

end subroutine ctypes
