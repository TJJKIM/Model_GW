!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for checking convergency of green function !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ConvCheckG (Dimen,Greenf_k_iWn,dummy_k_iWn,DiffElement,NumDiffElement,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,ConvTH)
Implicit none
!argument
integer,intent (in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim
complex(8) ,intent(in) ::  Greenf_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim), dummy_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)
integer,intent (inout) :: NumDiffElement
real(8),intent (inout) :: DiffElement, ConvTH

!local variable
integer kx,ky,kz,n,i,j

DiffElement=0.0;
NumDiffElement=0
if (Dimen==3) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do kz=1,K3GridNum
                    do n=1,iWnGridNum
                         do i=1,MatDim
                              do j=1,MatDim
                                   if ( abs( real(Greenf_k_iWn(kx,ky,kz,n,i,j)) - real(dummy_k_iWn(kx,ky,kz,n,i,j)) ) >ConvTH) then
                                           NumDiffElement=NumDiffElement+1
                                   endif
 
                                   if ( abs( aimag(Greenf_k_iWn(kx,ky,kz,n,i,j)) - aimag(dummy_k_iWn(kx,ky,kz,n,i,j)) ) >ConvTH ) then
                                           NumDiffElement=NumDiffElement+1
                                   endif
 
                                   DiffElement=DiffElement+abs( real(Greenf_k_iWn(kx,ky,kz,n,i,j)) - real(dummy_k_iWn(kx,ky,kz,n,i,j)) ) + abs( aimag(Greenf_k_iWn(kx,ky,kz,n,i,j)) - aimag(dummy_k_iWn(kx,ky,kz,n,i,j)) )
 
 
                              enddo
                         enddo
                    enddo
               enddo
          enddo
     enddo
endif

if (Dimen==2) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
              do n=1,iWnGridNum
                   do i=1,MatDim
                        do j=1,MatDim
                             if ( abs( real(Greenf_k_iWn(kx,ky,1,n,i,j)) - real(dummy_k_iWn(kx,ky,1,n,i,j)) ) >ConvTH) then
                                     NumDiffElement=NumDiffElement+1
                             endif

                             if ( abs( aimag(Greenf_k_iWn(kx,ky,1,n,i,j)) - aimag(dummy_k_iWn(kx,ky,1,n,i,j)) ) >ConvTH ) then
                                     NumDiffElement=NumDiffElement+1
                             endif

                             DiffElement=DiffElement+abs( real(Greenf_k_iWn(kx,ky,1,n,i,j)) - real(dummy_k_iWn(kx,ky,1,n,i,j)) ) + abs( aimag(Greenf_k_iWn(kx,ky,1,n,i,j)) - aimag(dummy_k_iWn(kx,ky,1,n,i,j)) )


                        enddo
                   enddo
              enddo
          enddo
     enddo
endif


if (Dimen==1) then
     do kx=1,K1GridNum
         do n=1,iWnGridNum
              do i=1,MatDim
                   do j=1,MatDim
                        if ( abs( real(Greenf_k_iWn(kx,1,1,n,i,j)) - real(dummy_k_iWn(kx,1,1,n,i,j)) ) >ConvTH) then
                                NumDiffElement=NumDiffElement+1
                        endif

                        if ( abs( aimag(Greenf_k_iWn(kx,1,1,n,i,j)) - aimag(dummy_k_iWn(kx,1,1,n,i,j)) ) >ConvTH ) then
                                NumDiffElement=NumDiffElement+1
                        endif

                        DiffElement=DiffElement+abs( real(Greenf_k_iWn(kx,1,1,n,i,j)) - real(dummy_k_iWn(kx,1,1,n,i,j)) ) + abs( aimag(Greenf_k_iWn(kx,1,1,n,i,j)) - aimag(dummy_k_iWn(kx,1,1,n,i,j)) )


                   enddo
              enddo
         enddo
     enddo
endif




return
end subroutine

