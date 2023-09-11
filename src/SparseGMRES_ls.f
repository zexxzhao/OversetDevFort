      subroutine SparseGMRES_ls_diag(
     &        lhsM,
     &        Utol, col, row,
     &        rhsG,
     &        sol,
     &        Kspaceu, Kspaceu_mn, 
     &        NNODZu, NSHLu, icntu, NSD
     &        )

      use mpi     
      implicit none

      integer NNODZu, NSHLu, icntu, NSD

      real(8) unorm_ref, rru, rrul

      integer n, i, j, k, iK, iKs, jK, lK, Kspaceu,
     &     Kspaceu_mn, 
     &     is, lenseg

      integer kk
      integer col(NNODZu+1), row(icntu)

      real(8)  lhsM (icntu)

      real(8) HBrg(Kspaceu+1,Kspaceu)
      real(8) Rcos(Kspaceu), Rsin(Kspaceu)

      real(8) uBrg1(NNODZu,Kspaceu+1)
      real(8) rhstmp1(NNODZu)
      real(8) rhsG(NNODZu),     sol(NNODZu)
      real(8) eBrg(Kspaceu+1),
     &     yBrg(Kspaceu),
     &     rr, unorm, epsnrm, beta,rrglob,
     &     ercheck, tmp, tmp1, Utol
      real(8)  lhsMdiag(NNODZu),temp1(NNODZu)
  
                                ! Algorithm from Shakib/Johan's theses
      rhstmp1 = rhsG

! Precondition Resudual

      lhsMdiag   = 0.0d0      
      
      do i = 1, NNODZu         
         do j = col(i), col(i+1)-1
            n = row(j)
            if (n == i) then               
                lhsMdiag(i) = LHSM(j)
            endif           
         enddo         
      enddo
 
      if (numnodes .gt. 1) then
            call commu(lhsMdiag,1,'in ')
            call commu(lhsMdiag,1,'out')      
      endif
      
      do i = 1, NNODZu
        if (lhsMdiag(i).ne.0d0) then
          lhsMdiag(i) = 1.0d0/lhsMdiag(i)
        endif
      enddo
  
      rr = sum(rhstmp1*rhstmp1)
      rrglob = rr
       if (numnodes .gt. 1) then
        call MPI_ALLREDUCE (rrglob,rr, 1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, 
     &        MPI_COMM_WORLD,mpi_err)
       endif

      if (rrglob == 0d0) then
        sol =0d0
        return
      endif

!------------- zero out common values -------------- 
      Rcos  = 0d0
      Rsin  = 0d0
      HBrg  = 0d0
      uBrg1 = 0d0

!---------------   Element-by-element Inverse Preconditioner   ---------------
        !  matrix-vector product
      rhstmp1 = lhsMdiag*rhstmp1
  
      uBrg1(:,1)   = rhstmp1
                                ! calculate norm of rhs
      rr = sum(rhstmp1*rhstmp1)
      rrglob = rr
      if (numnodes .gt. 1) then
          call MPI_ALLREDUCE (rrglob, rr, 1,
     &        MPI_DOUBLE_PRECISION,MPI_SUM, 
     &        MPI_COMM_WORLD,mpi_err)
      endif
      unorm   = sqrt(rr)

      epsnrm  = Utol * unorm    ! set up tolerance
                                ! set up RHS of the Hessenberg's problem
      eBrg    = 0d0
      eBrg(1) = unorm

      unorm = 1d0/unorm         ! normalize the first Krylov vector
      uBrg1(:,1) = uBrg1(:,1)   * unorm

      unorm_ref = unorm*1d2
    
                                ! loop through GMRES iterations
      do iK = 1, Kspaceu
     
         iKs = iK

         ! matrix-vector product
         temp1 = uBrg1(:,iKs)
c--------------------------------------------------------------------c             
         ! Periodicity (Slave = Master) - GL
         if (numnodes .gt. 1) call commu(temp1,1,'out')
                  
!         do i = 1, NNODZu
!           if ((IBC(i,4)==3).or.(IBC(i,5)==3).or.(IBC(i,6)==3)) then
!             temp1(i,:) = temp1(IPER(i),:) ! Slave = Master 
!           end if
!         end do
        
         ! Product        
         call SparseProd_LS(lhsM,col, row, temp1, uBrg1(:,iKs+1),
     &                      NNODZu, NSHLu, icntu, NSD)
     
         !Communicate to Masters, Zero out Slaves - LG         
!         do i = 1, NNODZu
!           if ((IBC(i,4)==3).or.(IBC(i,5)==3).or.(IBC(i,6)==3)) then
!               
!             uBrg1(IPER(i),:,iKs+1) = uBrg1(IPER(i),:,iKs+1) + 
!     &                               uBrg1(i,:,iKs+1) ! Master = Master+Slave
!             uBrg1(i,:,iKs+1) = 0d+0 ! Slave = zero
!               
!           end if          
!         end do
                  
         if (numnodes .gt. 1) call commu(uBrg1(:,iKs+1),1,'in ')         
c-------------------------------------------------------------------------c
                               
         uBrg1(:,iKs+1) = lhsMdiag*uBrg1(:,iKs+1) 

         do jK = 1, iKs+1       ! orthogonalize and get the norm

            if (jK .eq. 1) then

               rr = sum(uBrg1(:,iKs+1)*uBrg1(:,1)) !{u_{i+1}*u_1} vector
               rrglob=rr
               if (numnodes .gt. 1) 	            
     &              call MPI_ALLREDUCE (rrglob, rr, 1,
     &                 MPI_DOUBLE_PRECISION,MPI_SUM, 
     &                 MPI_COMM_WORLD,mpi_err)

               beta = rr

            else
                                !     project off jK-1 vector
               uBrg1(:  ,iKs+1) =
     &              uBrg1(:  ,iKs+1) - beta * uBrg1(:  ,jK-1)

               rr = sum(uBrg1(:,iKs+1)*uBrg1(:,jK)) !{u_{i+1}*u_1} vector
               rrglob = rr
             
               if (numnodes .gt. 1)
     &              call MPI_ALLREDUCE (rrglob, rr, 1,
     &                 MPI_DOUBLE_PRECISION,MPI_SUM, 
     &                 MPI_COMM_WORLD,mpi_err)

               beta = rr

            endif

            HBrg(jK,iKs) = beta ! put this in the Hessenberg Matrix

         enddo

         unorm           = sqrt(beta)
         HBrg(iKs+1,iKs) = unorm ! this fills the 1 sub diagonal band

!     normalize the Krylov vector
         unorm = 1d0/unorm      
                                ! normalize the next Krylov vector
         uBrg1(:  ,iKs+1) = uBrg1(:,iKs+1)   * unorm

c.... construct and reduce the Hessenberg Matrix
c     since there is only one subdiagonal we can use a Givens rotation
c     to rotate off each subdiagonal AS IT IS FORMED. We do this because it
c     allows us to check progress of solution and quit when satisfied.  Note
c     that all future K vects will put a subdiagonal in the next colmn so
c     there is no penalty to work ahead as  the rotation for the next vector
c     will be unaffected by this rotation.

c     H Y = E ========>   R_i H Y = R_i E

         do jK = 1, iKs-1
            tmp =  Rcos(jK) * HBrg(jK,iKs) + Rsin(jK) * HBrg(jK+1,iKs)
            HBrg(jK+1,iKs) = 
     &           - Rsin(jK) * HBrg(jK,  iKs)
     &           + Rcos(jK) * HBrg(jK+1,iKs)
            HBrg(jK,  iKs) =  tmp
         enddo

         tmp  = sqrt(HBrg(iKs,iKs)**2 + HBrg(iKs+1,iKs)**2)
         tmp1 = 1d0/tmp
         Rcos(iKs) = HBrg(iKs,  iKs) * tmp1
         Rsin(iKs) = HBrg(iKs+1,iKs) * tmp1
         HBrg(iKs,  iKs) = tmp
         HBrg(iKs+1,iKs) = 0d0
                                ! rotate eBrg    R_i E
         tmp         = + Rcos(iKs) * eBrg(iKs) + Rsin(iKs) * eBrg(iKs+1)
         eBrg(iKs+1) = - Rsin(iKs) * eBrg(iKs) + Rcos(iKs) * eBrg(iKs+1)
         eBrg(iKs)   = tmp
                                ! check for convergence
         ercheck=abs(eBrg(iKs+1))

         if ((ismaster) .and.(mod(iKs,400).eq.0)) then
               write(*,*) "Iteration =  ", iKs
               write(*,'(a,x,E12.4,x,E12.4,x,F12.6)')
     &	         "Residual, Goal, reduction (%)=  ",
     &              ercheck, epsnrm, ercheck*unorm_ref
         endif 

         if (ercheck .le. epsnrm .and. iK >= Kspaceu_mn) exit

      enddo                ! end of GMRES iteration loop

!     ------------------------->   Solution   <------------------------
!     if converged or end of Krylov space
                                ! solve for yBrg
      do jK = iKs, 1, -1
         yBrg(jK)     = eBrg(jK) / HBrg(jK,jK)
         eBrg(1:jK-1) = eBrg(1:jK-1) - yBrg(jK) * HBrg(1:jK-1,jK)
      enddo
                                ! update Dy
      do jK = 1, iKs
         sol = sol + yBrg(jK) * uBrg1(:  ,jK)
      enddo

      if (numnodes .gt. 1) then      
          call commu(sol,1,'out')    	    
      endif

      if ((ismaster) .and.(iKs.eq.Kspaceu))  
     &     write(*,9000) iKs, ercheck, ercheck*unorm_ref

      return

 2000 format(i10, f17.12)
 9000 format(
     &     10x,' number of GMRES iterations:', 2x,i10,/,
     &     10x,'       residual final value:', 2x,e10.4,/,
     &     10x,'       reduction percentage:', 2x,f10.7,/ )
      
      end
