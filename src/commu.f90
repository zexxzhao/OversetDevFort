!======================================================================
! This subroutine is responsible for interprocessor communication of
! the residual and solution vectors.
!
! input:
!     global(nshg,n): global vector to be communicated. Note that
!          this vector is local to the processor, (i.e.
!          not distributed across processors)
!     n:       second dimension of the array to be communicated
!     code:    = 'in' for communicating with the residual
!          = 'out' for cummunicating the solution
!---------------------------------------------------------------------
! The array ilwork describes the details of the communications.
! Each communication step (call of this routine) consists of a
! sequence of "tasks", where a task is defined as a communication
! between two processors where data is exchanged. This would imply
! that for a given processor, there will be as many tasks as there
! are processors with which it must communicate. Details of the
! ilwork array appear below.
!======================================================================

  subroutine commu(global, n, code)

    use mpi
    use commonvars

    implicit none

    character(len=3), intent(in) :: code

    integer :: stat(MPI_STATUS_SIZE, 2*numnodes), req(2*numnodes)
!!!  integer stat(MPI_STATUS_SIZE), req(2*numnodes)
    integer :: idl, jdl, itask, n, m, is, itemp, idof, kdof, j
    integer :: lenseg

    real(8) :: global(NNODE, n), rtemp(maxfrontu*n, numnodes)

    if (n == 1) then      ! like a scalar
      kdof = 1
    else if (n == NSD) then   ! like the normal vectors
      kdof = 2
    else if (n == (NSD*NSD)) then ! bdiag
      kdof = 3
    else if (n == ((NSD + 1)*(NSD + 1))) then ! bdiag
      kdof = 4
    end if

    rtemp = 0.0d0      ! Initialize

!... Note that when adding another kdof to the above set, we must
!... also make changes in ctypes.f and auxmpi.h

!---------------------------------------------------------------------
!  ilwork(1): number of tasks
!
!  The following information is contained in ilwork for each task:
!     itag: tag of the communication
!     iacc: == 0 if task is a send
!       == 1 if task is a recieve
!     iother: rank of processor with which this communication occurs
!     numseg: number of data "segments" to be sent or recieved. A
!     segment is defined as a continuous section of the global
!     vector to be communicated, (i.e. a group of nodes (or,
!     rather, "shape function coefficients") which occur
!     sequentially in the array global(nshg,n)).
!     isbeg:  location of the first segment in the array owned by the
!     current processor.
!
! The two types of communication are 'in', where the residual is being
! communicated, and 'out', where the solution is being communicated.
! Note that when the type is 'out', senders recieve and recievers send.
!
! The following comment pertains to a communication of type 'in':
!
!     If the task is a send, then all of the numseg segments are
!     sent with a single call to MPI_SEND. Where these segments live in
!     the array is built into the array sevsegtype, which is a common
!     array constructed in the subroutine "ctypes.f". In other words,
!     sevsegtype is a data type that describes the indices of the blocks
!     to be sent, in terms of there beginning index, and the length of
!     each segment. Using this, we can make a single send to take care of
!     all the segments for this task.
!
!     If the task is a recieve, then once the vector is recieved, the
!     recieved segments must be added to the correct locations in the
!     current array. These locations are described in ilwork as the
!     beginning position, then the length of the segment.
!---------------------------------------------------------------------

    numtask = ilworku(1)

    itkbeg = 1
    m = 0
    idl = 0

    ! loop through all the communication tasks
    do itask = 1, numtask
      m = m + 1
      itag = ilworku(itkbeg + 1)
      iacc = ilworku(itkbeg + 2)
      iother = ilworku(itkbeg + 3)

      if (iother < 0) iother = MPI_PROC_NULL

      numseg = ilworku(itkbeg + 4)

      ! when there's no communication,
      ! there should be no isgbeg..
      if (iother < 0) then
        isgbeg = 1
      else
        isgbeg = ilworku(itkbeg + 5)
      end if

      !--------------------------------------------------
      ! if iacc == 0, then this task is a send. (slave)
      !--------------------------------------------------
      if (iacc == 0) then

        ! residual communication
        if (code == 'in ') then

          call MPI_ISEND(global(isgbeg, 1), 1, &
                         sevsegtypeu(itask, kdof), iother, &
                         itag, MPI_COMM_WORLD, req(m), mpi_err)

!$$$       do is = 1, numseg
!$$$      isgbeg = ilworku (itkbeg + 3 + 2*is)
!$$$      lenseg = ilworku (itkbeg + 4 + 2*is)
!$$$      isgend = isgbeg + lenseg - 1
!$$$      global(isgbeg:isgend,:) = 0d+0
!$$$       enddo

!$$$       call MPI_SEND(global(isgbeg, 1),1,
!$$$     &      sevsegtypeu(itask,kdof),
!$$$     &      iother,itag,MPI_COMM_WORLD,
!$$$     &      mpi_err)

        end if

        ! solution communication
        if (code == 'out') then

          call MPI_IRECV(global(isgbeg, 1), 1, &
                         sevsegtypeu(itask, kdof), iother, &
                         itag, MPI_COMM_WORLD, req(m), mpi_err)

!$$$       call MPI_RECV(global(isgbeg,1), 1,
!$$$     &      sevsegtypeu(itask,kdof),
!$$$     &      iother, itag,MPI_COMM_WORLD,status,
!$$$     &      mpi_err)

        end if

        !------------------------------------------------------
        ! if iacc == 1, then this task is a recieve. (master)
        !------------------------------------------------------
      else

        ! masters receive the residuals
        if (code == 'in ') then

          ! determine the number of total number of nodes involved
          ! in this communication (lfront), including all segments
          lfrontu = 0
          do is = 1, numseg
            lenseg = ilworku(itkbeg + 4 + 2*is)
            lfrontu = lfrontu + lenseg
          end do

          ! recieve all segments for this task in a single step
          idl = idl + 1  ! stands for i Do Later, the number to fix later

          call MPI_IRECV(rtemp(1, idl), lfrontu*n, &
                         MPI_DOUBLE_PRECISION, iother, &
                         itag, MPI_COMM_WORLD, req(m), mpi_err)

!$$$       write(*,*) myid, "Receiving...", idl
!$$$       call MPI_RECV(rtemp(1,idl),lfrontu*n,
!$$$     &      MPI_DOUBLE_PRECISION,
!$$$     &      iother,itag,MPI_COMM_WORLD,status,
!$$$     &      mpi_err)
!$$$       write(*,*) myid, "Received", idl

        end if

        ! masters send out the solutions
        if (code == 'out') then

          call MPI_ISEND(global(isgbeg, 1), 1, &
                         sevsegtypeu(itask, kdof), iother, &
                         itag, MPI_COMM_WORLD, req(m), mpi_err)

!$$$       call MPI_SEND(global(isgbeg, 1),1,
!$$$     &      sevsegtypeu(itask,kdof),
!$$$     &      iother,itag,MPI_COMM_WORLD,
!$$$     &      mpi_err)

        end if
      end if

      itkbeg = itkbeg + 4 + 2*numseg

      call MPI_WAIT(req(m), stat(:, m), mpi_err)

    end do          ! end tasks loop

!!!  call MPI_WAITALL(m, req, stat, mpi_err)

    !-----------------------------------------------------------------
    ! Stuff added below is a delayed assembly of that which was
    ! communicated above but due to the switch to non-blocking
    ! receives could not be assembled until after the waitall.
    ! Only necessary for commu "in"
    !-----------------------------------------------------------------
    if (code == 'in ') then

      itkbeg = 1
      jdl = 0

      ! time to do all the segments that needed to be
      ! assembled into the global vector
      do j = 1, numtask

        iacc = ilworku(itkbeg + 2)
        numseg = ilworku(itkbeg + 4)

        ! if it's a master
        if (iacc == 1) then

          ! keep track of order of rtemp's
          jdl = jdl + 1

          ! add the recieved data to the global array on the current processor.
          ! Note that this involves splitting up the chunk of recieved data
          ! into its correct segment locations for the current processor.
          itemp = 1
          do idof = 1, n
            do is = 1, numseg
              isgbeg = ilworku(itkbeg + 3 + 2*is)
              lenseg = ilworku(itkbeg + 4 + 2*is)
              isgend = isgbeg + lenseg - 1

              global(isgbeg:isgend, idof) = global(isgbeg:isgend, idof) &
                                            + rtemp(itemp:itemp + lenseg - 1, jdl)

              itemp = itemp + lenseg
            end do
          end do

        else  ! end of receive (iacc=1), now do iacc=0 (zero sent)
!!!      else if (iacc == 0) then

          do is = 1, numseg
            isgbeg = ilworku(itkbeg + 3 + 2*is)
            lenseg = ilworku(itkbeg + 4 + 2*is)
            isgend = isgbeg + lenseg - 1
            global(isgbeg:isgend, :) = 0.0d0
          end do

        end if

        itkbeg = itkbeg + 4 + 2*numseg

      end do
    end if    ! end of commu "in"

  end subroutine commu

!======================================================================
! This subroutine is responsible for interprocessor communication.
! it sets the slave dofs to zero vectors.
!======================================================================
  subroutine zeroslaves(global, n)
    use mpi
    use commonvars
    implicit none

    integer, intent(in)    :: n
    real(8), intent(inout) :: global(NNODE, n)

    integer :: idl, jdl, itask, m, is, itemp, idof, kdof, j
    integer :: lenseg

    numtask = ilworku(1)
    itkbeg = 1
    jdl = 0
    do j = 1, numtask

      iacc = ilworku(itkbeg + 2)
      numseg = ilworku(itkbeg + 4)

      if (iacc == 0) then

        do is = 1, numseg
          isgbeg = ilworku(itkbeg + 3 + 2*is)
          lenseg = ilworku(itkbeg + 4 + 2*is)
          isgend = isgbeg + lenseg - 1
          global(isgbeg:isgend, 1:n) = 0.0d0
        end do

      end if

      itkbeg = itkbeg + 4 + 2*numseg
    end do
  end subroutine zeroslaves

!======================================================================
! This subroutine is responsible for interprocessor communication of
! the residual and solution vectors.
!
! input:
!     global(nshg,n): global vector to be communicated. Note that
!          this vector is local to the processor, (i.e.
!          not distributed across processors)
!     n:       second dimension of the array to be communicated
!     code:    = 'in' for communicating with the residual
!          = 'out' for cummunicating the solution
!
!---------------------------------------------------------------------
!
! The array ilwork describes the details of the communications.
! Each communication step (call of this routine) consists of a
! sequence of "tasks", where a task is defined as a communication
! between two processors where data is exchanged. This would imply
! that for a given processor, there will be as many tasks as there
! are processors with which it must communicate. Details of the
! ilwork array appear below.
!======================================================================
  subroutine icommu(global)

    use mpi
    use commonvars

    implicit none

    integer :: stat(MPI_STATUS_SIZE, 2*numnodes), req(2*numnodes)
    integer :: idl, jdl, itask, n, m, is, itemp, idof, kdof, j
    integer :: lenseg
    integer, dimension(NNODE) :: global, rtemp

    n = 1
    kdof = 5

    numtask = ilworku(1)

    itkbeg = 1
    m = 0
    idl = 0

    do itask = 1, numtask
      m = m + 1
      itag = ilworku(itkbeg + 1)
      iacc = ilworku(itkbeg + 2)
      iother = ilworku(itkbeg + 3)
      numseg = ilworku(itkbeg + 4)

      ! when there's no communication,
      ! there should be no isgbeg..
      if (iother < 0) then
        isgbeg = 1
        iother = MPI_PROC_NULL
      else
        isgbeg = ilworku(itkbeg + 5)
      end if
      rtemp = 0
      if (iacc == 0) then

        call MPI_ISEND(global(isgbeg), 1, &
                       sevsegtypeu(itask, kdof), &
                       iother, itag, MPI_COMM_WORLD, req(m), &
                       mpi_err)

      else

        call MPI_IRECV(rtemp(isgbeg), 1, &
                       sevsegtypeu(itask, kdof), &
                       iother, itag, MPI_COMM_WORLD, req(m), &
                       mpi_err)

      end if

      call MPI_WAIT(req(m), status, mpi_err)

      global = global + rtemp

      itkbeg = itkbeg + 4 + 2*numseg
    end do

    itkbeg = 1
    m = 0
    idl = 0

    do itask = 1, numtask
      m = m + 1
      itag = ilworku(itkbeg + 1)
      iacc = ilworku(itkbeg + 2)
      iother = ilworku(itkbeg + 3)
      numseg = ilworku(itkbeg + 4)

      ! when there's no communication,
      ! there should be no isgbeg..
      if (iother < 0) then
        isgbeg = 1
        iother = MPI_PROC_NULL
      else
        isgbeg = ilworku(itkbeg + 5)
      end if

      if (iacc == 0) then

        call MPI_IRECV(global(isgbeg), 1, &
                       sevsegtypeu(itask, kdof), &
                       iother, itag, MPI_COMM_WORLD, req(m), &
                       mpi_err)
      else
        call MPI_ISEND(global(isgbeg), 1, &
                       sevsegtypeu(itask, kdof), &
                       iother, itag, MPI_COMM_WORLD, req(m), &
                       mpi_err)
      end if
      call MPI_WAIT(req(m), status, mpi_err)

      itkbeg = itkbeg + 4 + 2*numseg

    end do          ! end tasks loop

  end subroutine icommu
