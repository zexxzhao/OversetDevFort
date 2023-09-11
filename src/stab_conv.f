      subroutine e3STAB_TAU(ui, Gij, tauM)
      
      use aAdjKeep      
      use commonvars
      implicit none 
      
      real(8), intent(in) :: ui(NSD), Gij(NSD,NSD)
      real(8) :: tauM
      real(8) :: taua, taut, dtfact
      integer :: i, j

c...  build remaining diagonal entries of tauM  
    
      dtfact = 4d0

      taua = 0d0
      do j = 1, NSD
	do i = 1, NSD
	  taua  = taua  + ui(i)*Gij(i,j)*ui(j)
	enddo
      enddo    

      taut = dtfact*Dtgl*Dtgl      
      tauM = 1d0/sqrt(taut+taua) 

      return
      end
!------------------------------------------------------------------     
      
      subroutine e3STAB_TAU_nt(ui, Gij, tauM)
      
      use aAdjKeep     
      use commonvars
      implicit none 
      
      real(8), intent(in) :: ui(NSD), Gij(NSD,NSD)
      real(8) :: tauM
      real(8) :: taua
      integer :: i, j

c...  build remaining diagonal entries of tauM  
    
      taua = 0d0
      do j = 1, NSD
	do i = 1, NSD
	  taua  = taua  + ui(i)*Gij(i,j)*ui(j)
	enddo
      enddo    
    
      tauM = 1d0/sqrt(taua+1d-15) 

      return
      end

!------------------------------------------------------------------
      subroutine e3DC_CAU(dfidxi,ui,Gij,res,tauM,kdc)

      use aAdjKeep
      use commonvars
      implicit none
      
      real(8), intent(in) :: dfidxi(NSD),ui(NSD), Gij(NSD,NSD),
     &                       res

      real(8) :: tauM, kdc(NSD,NSD)
      real(8) :: temp
      integer :: i, j

      temp = 0d0
      do j = 1, NSD
         do i = 1, NSD
            temp  = temp  + dfidxi(i)*Gij(i,j)*dfidxi(j)
         enddo
      enddo
      
      kdc = 0d0
      do i = 1, NSD      
        kdc(i,i) = abs(res)/sqrt(temp+1d-15)
      enddo      
      
      return
      end
      
      subroutine e3DC_CAUs(dfidxi,ui,Ginv,res,tauM,kdc)

      use aAdjKeep
      use commonvars
      implicit none
      
      real(8), intent(in) :: dfidxi(NSD), ui(NSD), Ginv(NSD,NSD),
     &                       res
      real(8) :: temp, tauM, kdc
      integer :: i, j

      temp = 0d0
      do j = 1, NSD
         do i = 1, NSD
            temp  = temp  + dfidxi(i)*Ginv(i,j)*dfidxi(j)
         enddo
      enddo
                
      kdc = abs(res)*sqrt(temp)/(sum(dfidxi*dfidxi)  +1d-10)    
      
      return
      end
      
