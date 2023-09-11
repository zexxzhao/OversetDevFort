!------------------------------------------------------------------------
!                                                                        
!        Main routine to call all the subroutines                        
!                                                                        
!------------------------------------------------------------------------
subroutine interpolate(omesh, odg, oug, ougm, opg, ophig, interp,         &
                       imesh, idg, iug, iugm, ipg, iphig)

  use class 
  implicit none
  
  type (mesh_class) :: omesh
  real(8), intent(in) :: odg(omesh%NNODE,omesh%NSD),  &
                         oug(omesh%NNODE,omesh%NSD),  &
                         ougm(omesh%NNODE,omesh%NSD), &
                         opg(omesh%NNODE), &
                         ophig(omesh%NNODE) 
   
  integer :: interp 
   
  type (mesh_class) :: imesh
  real(8) :: idg(imesh%NNODE,imesh%NSD), iug(imesh%NNODE,imesh%NSD), &
             iugm(imesh%NNODE,imesh%NSD), ipg(imesh%NNODE), &
             iphig(imesh%NNODE)  

  integer :: pn,ni,nj,nk,idx,i,j,k,ip,e1,e2,ln1,ln2,gn,gn1,e,n,NGAUSS 
  real(8) :: gp(imesh%NSHLmax,imesh%NSD),shlu(omesh%NSHLmax)
  real(8) :: xl(omesh%NSHLmax,omesh%NSD),wl(omesh%NSHLmax)  
  real(8) :: dl(omesh%NSHLmax,omesh%NSD),ul(omesh%NSHLmax,omesh%NSD), &
             uml(omesh%NSHLmax,omesh%NSD),pl(omesh%NSHLmax), phil(omesh%NSHLmax)
  real(8) :: xi(omesh%NSD),dx(omesh%NSD),di(omesh%NSD),ui(omesh%NSD), &
             umi(omesh%NSD), pi, phi
  logical :: flag(imesh%NNODE)
         
! 
!... Set interpolation
!       
  if (interp .eq.0) then
    NGAUSS      = 8   
  else if (interp .eq.1) then     
    NGAUSS      = 20   
  else
    write(*,*) 'Wrong interpolation_mode: ', interp     
    write(*,*) ' choose:  0 / 1'   
    stop
  endif 
   
  call getInterpolationPoints(gp,NGAUSS,omesh%NSD)        

!
!... Loop over element and Interpolate values 
!

  flag = .true.
  do e1 = 1, omesh%NELEM

    pn = omesh%EPID(e1)
    ni = omesh%EIJK(e1,1)
    nj = omesh%EIJK(e1,2)
    nk = omesh%EIJK(e1,3)  
   
    do j = 1, omesh%NSHLmax 
      idx = omesh%IEN(e1,j)
      xl(j,:)  = omesh%xg(idx,:)
      wl(j)    = omesh%wg(idx)
      dl(j,:)  = odg(idx,:)
      ul(j,:)  = oug(idx,:)
      uml(j,:) = ougm(idx,:)
      pl(j)    = opg(idx)
      phil(j)  = ophig(idx)

    enddo

! Interpolate values   
    do ln1 = 1, NGAUSS
      gn1 = imesh%IEN(e1,ln1)
      if (flag(gn1)) then
        call eval_shape_nurbs(gp(ln1,:),ni,nj,nk,omesh%patch(pn),xl,wl,shlu, omesh%NSD, omesh%NSHL)

        do j = 1, omesh%NSD
          xi(j)  = sum(xl(:,j)*shlu)
          di(j)  = sum(dl(:,j)*shlu)
          ui(j)  = sum(ul(:,j)*shlu) 
          umi(j) = sum(uml(:,j)*shlu) 
        enddo

        pi  = sum(pl*shlu)
        phi = sum(phil*shlu) 
                               
        imesh%xg(gn1,:) = xi
        idg  (gn1,:) = di
        iug  (gn1,:) = ui
        iugm (gn1,:) = umi
        ipg  (gn1)   = pi 
        iphig(gn1)   = phi
        
        flag(gn1)= .false.
      endif 
    enddo 
     
  enddo     
  
end subroutine  
