
      subroutine getInterpolationPoints(gp,NGAUSS,NSD)
     
      implicit none
      
      integer :: NGAUSS,NSD
      real(8) :: gp(NGAUSS,NSD)
      
      if (NGAUSS .eq.8 ) then
        gp(1,:) = (/-1d0,-1d0,-1d0 /)
        gp(2,:) = (/ 1d0,-1d0,-1d0 /)
        gp(3,:) = (/-1d0, 1d0,-1d0 /)
        gp(4,:) = (/ 1d0, 1d0,-1d0 /) 
      
        gp(5,:) = (/-1d0,-1d0, 1d0 /)
        gp(6,:) = (/ 1d0,-1d0, 1d0 /)
        gp(7,:) = (/-1d0, 1d0, 1d0 /)
        gp(8,:) = (/ 1d0, 1d0, 1d0 /)       
      else if (NGAUSS .eq.20 ) then

        gp(1,:) = (/-1d0,-1d0,-1d0 /)
        gp(2,:) = (/ 1d0,-1d0,-1d0 /)
        gp(3,:) = (/ 1d0, 1d0,-1d0 /)
        gp(4,:) = (/-1d0, 1d0,-1d0 /)
      
        gp(5,:) = (/-1d0,-1d0, 1d0 /)
        gp(6,:) = (/ 1d0,-1d0, 1d0 /)
        gp(7,:) = (/ 1d0, 1d0, 1d0 /)
        gp(8,:) = (/-1d0, 1d0, 1d0 /)
   
        
        gp(9, :) = (/ 0d0,-1d0,-1d0 /)
        gp(10,:) = (/ 1d0, 0d0,-1d0 /)
        gp(11,:) = (/ 0d0, 1d0,-1d0 /)
        gp(12,:) = (/-1d0, 0d0,-1d0 /)  
            
        gp(13,:) = (/ 0d0,-1d0, 1d0 /)
        gp(14,:) = (/ 1d0, 0d0, 1d0 /)
        gp(15,:) = (/ 0d0, 1d0, 1d0 /)
        gp(16,:) = (/-1d0, 0d0, 1d0 /)          
        
        gp(17,:) = (/-1d0,-1d0, 0d0 /)
        gp(18,:) = (/ 1d0,-1d0, 0d0 /)        
        gp(19,:) = (/ 1d0, 1d0, 0d0 /) 
        gp(20,:) = (/-1d0, 1d0, 0d0 /)

                   
      else   
        write(*,*) 'Wrong NGAUSS: ', NGAUSS  
        stop
      endif           
           
      end subroutine
!-------------------------------------------------               
!           
!-------------------------------------------------                     
      subroutine eval_shape_nurbs(gp,ni,nj,nk,patch,xl,wl,shlu, NSD, NSHL)
!,shgradgu,shhessg,dxidx,hess_flag)

      use class
      
      implicit none
      
!...  Input Variables 
      integer NSD, NSHL       
      real*8 gp(NSD)
      integer ni, nj, nk
      type (NURBSpatch) :: patch
      real*8 xl(NSHL,NSD),wl(NSHL)
      
      real*8 shlu(NSHL),shgradgu(NSHL,NSD),shhessg(NSHL,NSD,NSD)
      real*8 dxidx(NSD,NSD)
      !integer hess_flag

!...  Local Variables           
      integer nders
      real*8 u_hat, v_hat, w_hat, du, dv, dw,da
      integer Pu,Qu,Ru
      
      real*8 shgradlu(NSHL,NSD), tempshl(NSHL),tempshgradl(NSHL,NSD),shhessl(NSHL,6), tempshhessl(NSHL,6)

      real*8 dxdxi(NSD,NSD),dxdxixj(NSD,6), locLHS(6,6)

      real*8 Nnu(1,patch%P+1), Mmu(1,patch%Q+1)
      real*8  Ou(1,patch%R+1)    ! Mu = Mmu, since mu is viscosity
      
      real*8 u, v, w, denom_sum, derv_sum_U, derv_sum_V
      real*8   derv_sum_W, derv_sum_UU
      real*8   derv_sum_UV, derv_sum_UW,derv_sum_VV, derv_sum_VW,derv_sum_WW

      integer i, j, k, icount
      real*8 tmp
      
      !hess_flag = 0
      
!            icount=0
!	    do k = 0, Ru	
!               do j = 0, Qu	
!                  do i = 0, Pu
!		   icount = icount + 1
!	    	     B_NETu(ni-i,nj-j,nk-k,1:NSD) = xl(icount,:)
!	           enddo
!	        enddo
!	    enddo
                  
! ------------------------------------------------------------------

      u_hat = gp(1)
      v_hat = gp(2)
      w_hat = gp(3)
      
      Pu = patch%P
      Qu = patch%Q
      Ru = patch%R
  
!     Get u and v coordinates of integration point
      
      u = ((patch%U_KNOT(ni+1) - patch%U_KNOT(ni))*u_hat + &
           patch%U_KNOT(ni+1) + patch%U_KNOT(ni))/2d+0
      v = ((patch%V_KNOT(nj+1) - patch%V_KNOT(nj))*v_hat + &
           patch%V_KNOT(nj+1) + patch%V_KNOT(nj))/2d+0
      w = ((patch%W_KNOT(nk+1) - patch%W_KNOT(nk))*w_hat + &
           patch%W_KNOT(nk+1) + patch%W_KNOT(nk))/2d+0

!     Get knot span sizes
      
      du = patch%U_KNOT(ni+1)-patch%U_KNOT(ni)
      dv = patch%V_KNOT(nj+1)-patch%V_KNOT(nj)
      dw = patch%W_KNOT(nk+1)-patch%W_KNOT(nk)
      da = du*dv*dw

      
!     Evaluate 1D shape functions and derivatives each direction
      
      nders = 0

      
      call dersbasisfuns(ni,Pu,patch%MCP,u,nders,patch%U_KNOT,Nnu) ! calculate in u dir.
      call dersbasisfuns(nj,Qu,patch%NCP,v,nders,patch%V_KNOT,Mmu) ! calculate in v dir.
      call dersbasisfuns(nk,Ru,patch%OCP,w,nders,patch%W_KNOT,Ou ) ! calculate in v dir.

      
!     Form basis functions and derivatives dR/du and dR/dv
      
      icount = 0
      denom_sum = 0d+0
      derv_sum_U = 0d+0
      derv_sum_V = 0d+0
      derv_sum_W = 0d+0
      shgradlu = 0d+0
      shhessl = 0d+0
      shhessg = 0d+0
      tempshl = 0d+0
      tempshgradl = 0d+0
      tempshhessl = 0d+0
            
      do k = 0, Ru
        do j = 0, Qu
          do i = 0, Pu
            icount = icount+1
            
!...        basis functions
            shlu(icount) = Nnu(1,Pu+1-i)*Mmu(1,Qu+1-j)*Ou(1,Ru+1-k)*wl(icount)  !B_NET(ni-i,nj-j,nk-k,NSD+1)
            denom_sum = denom_sum + shlu(icount)
           
          enddo
        enddo
      enddo
           
      do i = 1,NSHL
        shlu(i) = shlu(i)/denom_sum
      enddo
           
      return
      end


!--------------------------      
      subroutine dersbasisfuns(i,pl,ml,u,nders,u_knotl,ders)
      
      IMPLICIT NONE

!c     --------------VARIABLE DECLARATIONS--------------------------------
!c...  knot span, degree of curve, number of control points, counters
      integer i, pl, ml, j, r, k, j1, j2, s1, s2, rk, pk,nders
!c     parameter value, vector of knots, derivative matrix
      real*8 u, u_knotl(pl+ml+1), ders(nders+1,pl+1), ndu(pl+1,pl+1),d

      real*8 left(pl+1), right(pl+1), saved, temp, a(2,pl+1)

!c     -------------------------------------------------------------------

      
      ndu(1,1) = 1d0
      do j = 1,pl
         left(j+1) = u - u_knotl(i+1-j)
         right(j+1) = u_knotl(i+j) - u
         saved = 0d0
         do r = 0,j-1
            ndu(j+1,r+1) = right(r+2) + left(j-r+1)
            temp = ndu(r+1,j)/ndu(j+1,r+1)
            ndu(r+1,j+1) = saved + right(r+2)*temp
            saved = left(j-r+1)*temp
         enddo
         ndu(j+1,j+1) = saved
      enddo
      
                                ! load basis functions
      do j = 0,pl
         ders(1,j+1) = ndu(j+1,pl+1)
      enddo
                                ! compute derivatives
      do r = 0,pl ! loop over function index
         s1 = 0
         s2 = 1                 ! alternate rows in array a
         a(1,1) = 1d0
                                ! loop to compute kth derivative
         do k = 1,nders
            d = 0d+0
            rk = r-k
            pk = pl-k
            if (r >= k) then
               a(s2+1,1) = a(s1+1,1)/ndu(pk+2,rk+1)
               d = a(s2+1,1)*ndu(rk+1,pk+1)
            endif
            if (rk >= -1) then
               j1 = 1
            else 
               j1 = -rk
            endif
            if ((r-1) <= pk) then
               j2 = k-1
            else 
               j2 = pl-r
            endif
            do j = j1,j2
               a(s2+1,j+1) = (a(s1+1,j+1) - a(s1+1,j))/ndu(pk+2,rk+j+1)
               d = d + a(s2+1,j+1)*ndu(rk+j+1,pk+1)
            enddo
            if (r <= pk) then
               a(s2+1,k+1) = -a(s1+1,k)/ndu(pk+2,r+1)
               d = d + a(s2+1,k+1)*ndu(r+1,pk+1)
            endif
            ders(k+1,r+1) = d
            j = s1
            s1 = s2
            s2 = j              ! switch rows
         enddo
      enddo
      
!c     Multiply through by the correct factors
      r = pl
      do k = 1,nders
         do j = 0,pl
            ders(k+1,j+1) = ders(k+1,j+1)*r
         enddo
         r = r*(pl-k)
      enddo

      return
      end
