!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for fourier transformation from k_iWn to k_tau !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine k_iWn2k_tau(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Greenf_k_iWn,dummy_k_tau)
Implicit none
! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim
real(8), intent(in) :: K1list(K1GridNum,3), K2list(K2GridNum,3), K3list(K3GridNum,3),Taulist(TauGridNum), R1list(R1GridNum,3), R2list(R2GridNum,3),R3list(R3GridNum,3),beta
complex(8), intent(in) :: iWnlist(iWnGridNum), Greenf_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)
complex(8), intent(inout) :: dummy_k_tau(K1GridNum,K2GridNum,K3GridNum,TauGridNum,MatDim,MatDim)

! local variable
integer i,j,kx,ky,kz,l,E,n
complex(8) Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim),dummy3(MatDim,MatDim)
complex(8) gtilda2(MatDim,MatDim), gtilda3(MatDim,MatDim), MatGcut(MatDim,MatDim), asymtoticG(MatDim,MatDim)


if (Dimen==3) then

     do l=1,TauGridNum !time index sum
          gtilda2(:,:)=(0.0d0,0.0d0)
          gtilda3(:,:)=(0.0d0,0.0d0) !initalize gtilda for calculating asymtotic part of matsubar green function
          do kx=1,(K1GridNum) ! kx-index sum
               do ky=1,(K2GridNum) ! ky-index sum
                    do kz=1,(K3GridNum) ! kz-index sum
                         dummy(:,:)=(0.0d0,0.0d0) !initialize
                         call Identity(MatDim,Iden)
                         gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( Greenf_k_iWn(kx,ky,kz,iWnGridNum,:,:)+Greenf_k_iWn(kx,ky,kz,1,:,:)  )/(2.0d0)
                         gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( Greenf_k_iWn(kx,ky,kz,iWnGridNum,:,:)-Greenf_k_iWn(kx,ky,kz,1,:,:) )/(2.0d0) - Iden/iWnlist(iWnGridNum) )
     
           !    do n=1,(iWnGridNum) ! Matsubara frequency iWn sumation 
            !        dummy= dummy + ((1.0d0)/(beta))* ( Greenf_k_iWn(k,n,:,:) - Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) - gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(-iWnlist(n)*(taulist(l)))
     
     
                         do n=(int(iWnGridNum/2)+1),(iWnGridNum) ! Matsubara frequency iWn sumation
                              dummy= dummy + ((1.0d0)/(beta))* ( Greenf_k_iWn(kx,ky,kz,n,:,:) - Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) - gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(-iWnlist(n)*(taulist(l)))
                              dummy= dummy + ((1.0d0)/(beta))* ( Greenf_k_iWn(kx,ky,kz,iWnGridNum+1-n,:,:) + Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) + gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(iWnlist(n)*(taulist(l)))
                         end do
     
          !     end do
     
                         asymtoticG = -  ( Iden/(2.0d0) ) + ( gtilda2*beta*( (taulist(l))/((2.0d0)*(beta))-(0.25d0)) ) - ( gtilda3*beta*beta/(4.0d0)*( ( (taulist(l))*(taulist(l)) )/(beta*beta) - (taulist(l))/beta ) )
                         dummy=dummy+asymtoticG
                         dummy_k_tau(kx,ky,kz,l,:,:)=dummy
                   end do
               end do
          end do
     end do

endif



if (Dimen==2) then

     do l=1,TauGridNum !time index sum
          gtilda2(:,:)=(0.0d0,0.0d0)
          gtilda3(:,:)=(0.0d0,0.0d0) !initalize gtilda for calculating asymtotic part of matsubar green function
          do kx=1,(K1GridNum) ! kx-index sum
               do ky=1,(K2GridNum) ! ky-index sum
                   dummy(:,:)=(0.0d0,0.0d0) !initialize
                   call Identity(MatDim,Iden)
                   gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( Greenf_k_iWn(kx,ky,1,iWnGridNum,:,:)+Greenf_k_iWn(kx,ky,1,1,:,:)  )/(2.0d0)
                   gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( Greenf_k_iWn(kx,ky,1,iWnGridNum,:,:)-Greenf_k_iWn(kx,ky,1,1,:,:) )/(2.0d0) - Iden/iWnlist(iWnGridNum) )
     
           !   1,(iWnGridNum) ! Matsubara frequency iWn sumation 
            !  dummy= dummy + ((1.0d0)/(beta))* ( Greenf_k_iWn(k,n,:,:) - Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) - gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(-iWnlist(n)*(taulist(l)))

     
                   do n=(int(iWnGridNum/2)+1),(iWnGridNum) ! Matsubara frequency iWn sumation
                        dummy= dummy + ((1.0d0)/(beta))* ( Greenf_k_iWn(kx,ky,1,n,:,:) - Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) - gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(-iWnlist(n)*(taulist(l)))
                        dummy= dummy + ((1.0d0)/(beta))* ( Greenf_k_iWn(kx,ky,1,iWnGridNum+1-n,:,:) + Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) + gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(iWnlist(n)*(taulist(l)))
                   end do
     
          !    o
               
                   asymtoticG = -  ( Iden/(2.0d0) ) + ( gtilda2*beta*( (taulist(l))/((2.0d0)*(beta))-(0.25d0)) ) - ( gtilda3*beta*beta/(4.0d0)*( ( (taulist(l))*(taulist(l)) )/(beta*beta) - (taulist(l))/beta ) )
                   dummy=dummy+asymtoticG
                   dummy_k_tau(kx,ky,1,l,:,:)=dummy
               end do
          end do
     end do

endif




if (Dimen==1) then

     do l=1,TauGridNum !time index sum
          gtilda2(:,:)=(0.0d0,0.0d0)
          gtilda3(:,:)=(0.0d0,0.0d0) !initalize gtilda for calculating asymtotic part of matsubar green function
          do kx=1,(K1GridNum) ! kx-index sum
              dummy(:,:)=(0.0d0,0.0d0) !initialize
              call Identity(MatDim,Iden)
              gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( Greenf_k_iWn(kx,1,1,iWnGridNum,:,:)+Greenf_k_iWn(kx,1,1,1,:,:)  )/(2.0d0)
              gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( Greenf_k_iWn(kx,1,1,iWnGridNum,:,:)-Greenf_k_iWn(kx,1,1,1,:,:) )/(2.0d0) - Iden/iWnlist(iWnGridNum) )
     

     
              do n=(int(iWnGridNum/2)+1),(iWnGridNum) ! Matsubara frequency iWn sumation
                   dummy= dummy + ((1.0d0)/(beta))* ( Greenf_k_iWn(kx,1,1,n,:,:) - Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) - gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(-iWnlist(n)*(taulist(l)))
                   dummy= dummy + ((1.0d0)/(beta))* ( Greenf_k_iWn(kx,1,1,iWnGridNum+1-n,:,:) + Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) + gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(iWnlist(n)*(taulist(l)))
              end do
     
          
              asymtoticG = -  ( Iden/(2.0d0) ) + ( gtilda2*beta*( (taulist(l))/((2.0d0)*(beta))-(0.25d0)) ) - ( gtilda3*beta*beta/(4.0d0)*( ( (taulist(l))*(taulist(l)) )/(beta*beta) - (taulist(l))/beta ) )
              dummy=dummy+asymtoticG
              dummy_k_tau(kx,1,1,l,:,:)=dummy
          end do
     end do

endif



return
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for fourier transformation from k_tau to k_iWn !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine k_tau2k_iWn(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_k_iWn)

Implicit none
! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum, MatDim
real(8), intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GRidNum,3), Taulist(TauGridNum), R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3), beta
complex(8), intent(in) :: iWnlist(iWnGridNum), dummy_k_tau(K1GridNum,K2GridNum,K3GridNum,TauGridNum,MatDim,MatDim)
complex(8), intent(inout) :: dummy_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)

! local variable
integer i,j,kx,ky,kz,l,E,n
real (8) fR(TauGridNum,MatDim,MatDim), fIm(TauGridNum,MatDim,MatDim), gR(TauGridNum,MatDim,MatDim),gIM(TauGridNum,MatDim,MatDim)
complex(8) Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG(MatDim,MatDim),dummy(MatDim,MatDim)

if (Dimen==3) then
     do n=1,iWnGridNum !iWn loop
          do kx=1,(K1GridNum) ! kx loop
               do ky=1,(K2GridNum) !ky loop
                    do kz=1,(K3GridNum) !kz loop
                    !dummy(:,:)=(0.0,0.0)
                         do l=1,(TauGridNum) ! time integration for FT
                              fR(l,1,1)= real( dummy_k_tau(kx,ky,kz,l,1,1) )*cos( aimag(iWnlist(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,kz,l,1,1) )*sin( aimag(iWnlist(n)) * Taulist(l) )
                              fIm(l,1,1)= real( dummy_k_tau(kx,ky,kz,l,1,1) )*sin( aimag(iWnlist(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,kz,l,1,1) )*cos( aimag(iWnlist(n)) * Taulist(l))
 
                              fR(l,1,2)= real( dummy_k_tau(kx,ky,kz,l,1,2) )*cos( aimag(iWnlist(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,kz,l,1,2) )*sin( aimag(iWnlist(n)) * Taulist(l) )
                              fIm(l,1,2)= real( dummy_k_tau(kx,ky,kz,l,1,2) )*sin( aimag(iWnlist(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,kz,l,1,2) )*cos( aimag(iWnlist(n)) * Taulist(l))
 
 
                              fR(l,2,1)= real( dummy_k_tau(kx,ky,kz,l,2,1) )*cos( aimag(iWnlist(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,kz,l,2,1) )*sin( aimag(iWnlist(n)) * Taulist(l) )
                              fIm(l,2,1)= real( dummy_k_tau(kx,ky,kz,l,2,1) )*sin( aimag(iWnlist(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,kz,l,2,1) )*cos( aimag(iWnlist(n)) * Taulist(l))
 
 
                              fR(l,2,2)= real( dummy_k_tau(kx,ky,kz,l,2,2) )*cos( aimag(iWnlist(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,kz,l,2,2) )*sin( aimag(iWnlist(n)) * Taulist(l) )
                              fIm(l,2,2)= real( dummy_k_tau(kx,ky,kz,l,2,2) )*sin( aimag(iWnlist(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,kz,l,2,2) )*cos( aimag(iWnlist(n)) * Taulist(l))
 
 
                         end do

                         call fderiv(-1,TauGridNum,Taulist,fR(:,1,1),gR(:,1,1))
                         call fderiv(-1,TauGridNum,Taulist,fIm(:,1,1),gIm(:,1,1))
                         call fderiv(-1,TauGridNum,Taulist,fR(:,1,2),gR(:,1,2))
                         call fderiv(-1,TauGridNum,Taulist,fIm(:,1,2),gIm(:,1,2))
                         call fderiv(-1,TauGridNum,Taulist,fR(:,2,1),gR(:,2,1))
                         call fderiv(-1,TauGridNum,Taulist,fIm(:,2,1),gIm(:,2,1))
                         call fderiv(-1,TauGridNum,Taulist,fR(:,2,2),gR(:,2,2))
                         call fderiv(-1,TauGridNum,Taulist,fIm(:,2,2),gIm(:,2,2))
 
                         dummy_k_iWn(kx,ky,kz,n,:,:)=gR(TauGridNum,:,:)+gIm(TauGridNum,:,:)*(0.0d0,1.0d0)
                    end do
               end do
          end do
     end do
endif


if (Dimen==2) then
     do n=1,iWnGridNum !iWn loop
          do kx=1,(K1GridNum) ! kx loop
               do ky=1,(K2GridNum) !ky loop
                         do l=1,(TauGridNum) ! time integration for FT
                              fR(l,1,1)= real( dummy_k_tau(kx,ky,1,l,1,1) )*cos( aimag(iWnlist(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,1,l,1,1) )*sin( aimag(iWnlist(n)) * Taulist(l) )
                              fIm(l,1,1)= real( dummy_k_tau(kx,ky,1,l,1,1) )*sin( aimag(iWnlist(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,1,l,1,1) )*cos( aimag(iWnlist(n)) * Taulist(l))

                              fR(l,1,2)= real( dummy_k_tau(kx,ky,1,l,1,2) )*cos( aimag(iWnlist(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,1,l,1,2) )*sin( aimag(iWnlist(n)) * Taulist(l) )
                              fIm(l,1,2)= real( dummy_k_tau(kx,ky,1,l,1,2) )*sin( aimag(iWnlist(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,1,l,1,2) )*cos( aimag(iWnlist(n)) * Taulist(l))


                              fR(l,2,1)= real( dummy_k_tau(kx,ky,1,l,2,1) )*cos( aimag(iWnlist(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,1,l,2,1) )*sin( aimag(iWnlist(n)) * Taulist(l) )
                              fIm(l,2,1)= real( dummy_k_tau(kx,ky,1,l,2,1) )*sin( aimag(iWnlist(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,1,l,2,1) )*cos( aimag(iWnlist(n)) * Taulist(l))


                              fR(l,2,2)= real( dummy_k_tau(kx,ky,1,l,2,2) )*cos( aimag(iWnlist(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,1,l,2,2) )*sin( aimag(iWnlist(n)) * Taulist(l) )
                              fIm(l,2,2)= real( dummy_k_tau(kx,ky,1,l,2,2) )*sin( aimag(iWnlist(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,1,l,2,2) )*cos( aimag(iWnlist(n)) * Taulist(l))


                         end do

                         call fderiv(-1,TauGridNum,Taulist,fR(:,1,1),gR(:,1,1))
                         call fderiv(-1,TauGridNum,Taulist,fIm(:,1,1),gIm(:,1,1))
                         call fderiv(-1,TauGridNum,Taulist,fR(:,1,2),gR(:,1,2))
                         call fderiv(-1,TauGridNum,Taulist,fIm(:,1,2),gIm(:,1,2))
                         call fderiv(-1,TauGridNum,Taulist,fR(:,2,1),gR(:,2,1))
                         call fderiv(-1,TauGridNum,Taulist,fIm(:,2,1),gIm(:,2,1))
                         call fderiv(-1,TauGridNum,Taulist,fR(:,2,2),gR(:,2,2))
                         call fderiv(-1,TauGridNum,Taulist,fIm(:,2,2),gIm(:,2,2))

                         dummy_k_iWn(kx,ky,1,n,:,:)=gR(TauGridNum,:,:)+gIm(TauGridNum,:,:)*(0.0d0,1.0d0)
               end do
          end do
     end do
endif



if (Dimen==1) then
     do n=1,iWnGridNum !iWn loop
          do kx=1,(K1GridNum) ! kx loop
              do l=1,(TauGridNum) ! time integration for FT
                   fR(l,1,1)= real( dummy_k_tau(kx,1,1,l,1,1) )*cos( aimag(iWnlist(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,1,1,l,1,1) )*sin( aimag(iWnlist(n)) * Taulist(l) )
                   fIm(l,1,1)= real( dummy_k_tau(kx,1,1,l,1,1) )*sin( aimag(iWnlist(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,1,1,l,1,1) )*cos( aimag(iWnlist(n)) * Taulist(l))

                   fR(l,1,2)= real( dummy_k_tau(kx,1,1,l,1,2) )*cos( aimag(iWnlist(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,1,1,l,1,2) )*sin( aimag(iWnlist(n)) * Taulist(l) )
                   fIm(l,1,2)= real( dummy_k_tau(kx,1,1,l,1,2) )*sin( aimag(iWnlist(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,1,1,l,1,2) )*cos( aimag(iWnlist(n)) * Taulist(l))


                   fR(l,2,1)= real( dummy_k_tau(kx,1,1,l,2,1) )*cos( aimag(iWnlist(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,1,1,l,2,1) )*sin( aimag(iWnlist(n)) * Taulist(l) )
                   fIm(l,2,1)= real( dummy_k_tau(kx,1,1,l,2,1) )*sin( aimag(iWnlist(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,1,1,l,2,1) )*cos( aimag(iWnlist(n)) * Taulist(l))


                   fR(l,2,2)= real( dummy_k_tau(kx,1,1,l,2,2) )*cos( aimag(iWnlist(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,1,1,l,2,2) )*sin( aimag(iWnlist(n)) * Taulist(l) )
                   fIm(l,2,2)= real( dummy_k_tau(kx,1,1,l,2,2) )*sin( aimag(iWnlist(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,1,1,l,2,2) )*cos( aimag(iWnlist(n)) * Taulist(l))


              end do

              call fderiv(-1,TauGridNum,Taulist,fR(:,1,1),gR(:,1,1))
              call fderiv(-1,TauGridNum,Taulist,fIm(:,1,1),gIm(:,1,1))
              call fderiv(-1,TauGridNum,Taulist,fR(:,1,2),gR(:,1,2))
              call fderiv(-1,TauGridNum,Taulist,fIm(:,1,2),gIm(:,1,2))
              call fderiv(-1,TauGridNum,Taulist,fR(:,2,1),gR(:,2,1))
              call fderiv(-1,TauGridNum,Taulist,fIm(:,2,1),gIm(:,2,1))
              call fderiv(-1,TauGridNum,Taulist,fR(:,2,2),gR(:,2,2))
              call fderiv(-1,TauGridNum,Taulist,fIm(:,2,2),gIm(:,2,2))

              dummy_k_iWn(kx,1,1,n,:,:)=gR(TauGridNum,:,:)+gIm(TauGridNum,:,:)*(0.0d0,1.0d0)
          end do
     end do
endif




return
end subroutine





























!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for fourier transformation from k_iWn to k_tau !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine k_iWn2k_tau_Even(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist2,Taulist,R1list,R2list,R3list,dummy_k_iWn2,dummy_k_tau)
Implicit none
! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim
real(8), intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3), Taulist(TauGridNum), R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3),beta
complex(8), intent(in) :: iWnlist2(iWnGridNum+1), dummy_k_iWn2(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim)
complex(8), intent(inout) :: dummy_k_tau(K1GridNum,K2GridNum,K3GridNum,TauGridNum,MatDim,MatDim)

! local variable
integer i,j,kx,ky,kz,l,E,n,Z
complex(8) Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim),dummy3(MatDim,MatDim)
complex(8) w2tilda(MatDim,MatDim), MatGcut(MatDim,MatDim), asymtoticW(MatDim,MatDim)


Z=(iWnGridNum/2)+1

if (Dimen==3) then
     do l=1,(TauGridNum) !time index sum
     w2tilda(:,:)=(0.0d0,0.0d0) !initalize gtilda for calculating asymtotic part of matsubar green function
          do kx=1,(K1GridNum) ! kx-index sum
               do ky=1,(K2GridNum) ! ky-index sum
                    do kz=1,(K3GridNum) ! kz-index sum
                         dummy(:,:)=(0.0d0,0.0d0) !initialize
                         call Identity(MatDim,Iden)
                         w2tilda(:,:)=iWnlist2(iWnGridNum+1)*iWnlist2(iWnGridNum+1)/(2.0d0)*( dummy_k_iWn2(kx,ky,kz,iWnGridNum+1,:,:)+dummy_k_iWn2(kx,ky,kz,1,:,:) )
                         dummy=dummy+((1.0d0)/beta)*dummy_k_iWn2(kx,ky,kz,Z,:,:)
                         do n=Z+1,(iWnGridNum+1) ! Matsubara frequency iWn sumation
 
                              dummy=dummy+((1.0d0)/beta)*(dummy_k_iWn2(kx,ky,kz,n,:,:)-w2tilda/(iWnlist2(n)*iWnlist2(n)))*exp(-iWnlist2(n)*Taulist(l))
                              dummy=dummy+((1.0d0)/beta)*(dummy_k_iWn2(kx,ky,kz,iWnGridNum+2-n,:,:)-w2tilda/(iWnlist2(iWnGridNum+2-n)*iWnlist2(iWnGridNum+2-n)))*exp(-iWnlist2(iWnGridNum+2-n)*Taulist(l))
                         end do
 
                         asymtoticW =  -w2tilda*(beta/(2.0d0))*( (Taulist(l)/beta)*(Taulist(l)/beta) - (Taulist(l)/beta) + (1.0d0)/(6.0d0) )
                         dummy=dummy+asymtoticW
                         dummy_k_tau(kx,ky,kz,l,:,:)=dummy
                    end do
               end do
          end do
     end do
endif


if (Dimen==2) then
     do l=1,(TauGridNum) !time index sum
     w2tilda(:,:)=(0.0d0,0.0d0) !initalize gtilda for calculating asymtotic part of matsubar green function
          do kx=1,(K1GridNum) ! kx-index sum
               do ky=1,(K2GridNum) ! ky-index sum
                   dummy(:,:)=(0.0d0,0.0d0) !initialize
                   call Identity(MatDim,Iden)
                   w2tilda(:,:)=iWnlist2(iWnGridNum+1)*iWnlist2(iWnGridNum+1)/(2.0d0)*( dummy_k_iWn2(kx,ky,1,iWnGridNum+1,:,:)+dummy_k_iWn2(kx,ky,1,1,:,:) )
                   dummy=dummy+((1.0d0)/beta)*dummy_k_iWn2(kx,ky,1,Z,:,:)
                   do n=Z+1,(iWnGridNum+1) ! Matsubara frequency iWn sumation

                        dummy=dummy+((1.0d0)/beta)*(dummy_k_iWn2(kx,ky,1,n,:,:)-w2tilda/(iWnlist2(n)*iWnlist2(n)))*exp(-iWnlist2(n)*Taulist(l))
                        dummy=dummy+((1.0d0)/beta)*(dummy_k_iWn2(kx,ky,1,iWnGridNum+2-n,:,:)-w2tilda/(iWnlist2(iWnGridNum+2-n)*iWnlist2(iWnGridNum+2-n)))*exp(-iWnlist2(iWnGridNum+2-n)*Taulist(l))
                   end do

                   asymtoticW =  -w2tilda*(beta/(2.0d0))*( (Taulist(l)/beta)*(Taulist(l)/beta) - (Taulist(l)/beta) + (1.0d0)/(6.0d0) )
                   dummy=dummy+asymtoticW
                   dummy_k_tau(kx,ky,1,l,:,:)=dummy
               end do
          end do
     end do
endif


if (Dimen==1) then
     do l=1,(TauGridNum) !time index sum
     w2tilda(:,:)=(0.0d0,0.0d0) !initalize gtilda for calculating asymtotic part of matsubar green function
          do kx=1,(K1GridNum) ! kx-index sum
              dummy(:,:)=(0.0d0,0.0d0) !initialize
              call Identity(MatDim,Iden)
              w2tilda(:,:)=iWnlist2(iWnGridNum+1)*iWnlist2(iWnGridNum+1)/(2.0d0)*( dummy_k_iWn2(kx,1,1,iWnGridNum+1,:,:)+dummy_k_iWn2(kx,1,1,1,:,:) )
              dummy=dummy+((1.0d0)/beta)*dummy_k_iWn2(kx,1,1,Z,:,:)
              do n=Z+1,(iWnGridNum+1) ! Matsubara frequency iWn sumation

                   dummy=dummy+((1.0d0)/beta)*(dummy_k_iWn2(kx,1,1,n,:,:)-w2tilda/(iWnlist2(n)*iWnlist2(n)))*exp(-iWnlist2(n)*Taulist(l))
                   dummy=dummy+((1.0d0)/beta)*(dummy_k_iWn2(kx,1,1,iWnGridNum+2-n,:,:)-w2tilda/(iWnlist2(iWnGridNum+2-n)*iWnlist2(iWnGridNum+2-n)))*exp(-iWnlist2(iWnGridNum+2-n)*Taulist(l))
              end do

              asymtoticW =  -w2tilda*(beta/(2.0d0))*( (Taulist(l)/beta)*(Taulist(l)/beta) - (Taulist(l)/beta) + (1.0d0)/(6.0d0) )
              dummy=dummy+asymtoticW
              dummy_k_tau(kx,1,1,l,:,:)=dummy
          end do
     end do
endif



return
end subroutine
                 











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for fourier transformation from k_tau to k_iWn !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine k_tau2k_iWn_Even(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist2,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_k_iWn2)

Implicit none
! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum, R3GridNum,MatDim
real(8), intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3), Taulist(TauGridNum), R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3), beta
complex(8), intent(in) :: iWnlist2(iWnGridNum+1), dummy_k_tau(K1GridNum,K2GridNum,K3GridNum,TauGridNum,MatDim,MatDim)
complex(8), intent(inout) :: dummy_k_iWn2(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim)

! local variable
integer i,j,kx,ky,kz,l,E,n
real (8) fR(TauGridNum,MatDim,MatDim), fIm(TauGridNum,MatDim,MatDim), gR(TauGridNum,MatDim,MatDim),gIM(TauGridNum,MatDim,MatDim)
complex(8) Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim),dummy3(MatDim,MatDim)




if (Dimen==3) then
     do n=1,(iWnGridNum+1) !iWn loop
          do kx=1,(K1GridNum) ! kx loop
               do ky=1,(K2GridNum) ! ky loop
                    do kz=1,(K3GridNum) ! kz loop
                    !dummy(:,:)=(0.0,0.0) !initialize
                    !dummy(:,:)= dummy(:,:) + (  dummy_k_tau(k,1,:,:)*exp(iWnlist2(n)*Taulist(1))  )*(Taulist(1))
                         do l=1,TauGridNum ! time integration for FT
 
                              fR(l,1,1)= real( dummy_k_tau(kx,ky,kz,l,1,1) )*cos( aimag(iWnlist2(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,kz,l,1,1) )*sin( aimag(iWnlist2(n)) * Taulist(l) )
                              fIm(l,1,1)= real( dummy_k_tau(kx,ky,kz,l,1,1) )*sin( aimag(iWnlist2(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,kz,l,1,1) )*cos( aimag(iWnlist2(n)) * Taulist(l))
 
                              fR(l,1,2)= real( dummy_k_tau(kx,ky,kz,l,1,2) )*cos( aimag(iWnlist2(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,kz,l,1,2) )*sin( aimag(iWnlist2(n)) * Taulist(l) )
                              fIm(l,1,2)= real( dummy_k_tau(kx,ky,kz,l,1,2) )*sin( aimag(iWnlist2(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,kz,l,1,2) )*cos( aimag(iWnlist2(n)) * Taulist(l))
 
 
                              fR(l,2,1)= real( dummy_k_tau(kx,ky,kz,l,2,1) )*cos( aimag(iWnlist2(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,kz,l,2,1) )*sin( aimag(iWnlist2(n)) * Taulist(l) )
                              fIm(l,2,1)= real( dummy_k_tau(kx,ky,kz,l,2,1) )*sin( aimag(iWnlist2(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,kz,l,2,1) )*cos( aimag(iWnlist2(n)) * Taulist(l))
 
 
                              fR(l,2,2)= real( dummy_k_tau(kx,ky,kz,l,2,2) )*cos( aimag(iWnlist2(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,kz,l,2,2) )*sin( aimag(iWnlist2(n)) * Taulist(l) )
                              fIm(l,2,2)= real( dummy_k_tau(kx,ky,kz,l,2,2) )*sin( aimag(iWnlist2(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,kz,l,2,2) )*cos( aimag(iWnlist2(n)) * Taulist(l))
 
 
 
                         end do
 
                         call fderiv(-1,TauGridNum,Taulist,fR(:,1,1),gR(:,1,1))
                         call fderiv(-1,TauGridNum,Taulist,fIm(:,1,1),gIm(:,1,1))
                         call fderiv(-1,TauGridNum,Taulist,fR(:,1,2),gR(:,1,2))
                         call fderiv(-1,TauGridNum,Taulist,fIm(:,1,2),gIm(:,1,2))
                         call fderiv(-1,TauGridNum,Taulist,fR(:,2,1),gR(:,2,1))
                         call fderiv(-1,TauGridNum,Taulist,fIm(:,2,1),gIm(:,2,1))
                         call fderiv(-1,TauGridNum,Taulist,fR(:,2,2),gR(:,2,2))
                         call fderiv(-1,TauGridNum,Taulist,fIm(:,2,2),gIm(:,2,2))
 
                         dummy_k_iWn2(kx,ky,kz,n,:,:)=gR(TauGridNum,:,:)+gIm(TauGridNum,:,:)*(0.0d0,1.0d0)
 
 
 
 
                    end do
               end do
          end do
     end do
endif

if (Dimen==2) then
     do n=1,(iWnGridNum+1) !iWn loop
          do kx=1,(K1GridNum) ! kx loop
               do ky=1,(K2GridNum) ! ky loop
                     do l=1,TauGridNum ! time integration for FT
                  
                          fR(l,1,1)= real( dummy_k_tau(kx,ky,1,l,1,1) )*cos( aimag(iWnlist2(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,1,l,1,1) )*sin( aimag(iWnlist2(n)) * Taulist(l) )
                          fIm(l,1,1)= real( dummy_k_tau(kx,ky,1,l,1,1) )*sin( aimag(iWnlist2(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,1,l,1,1) )*cos( aimag(iWnlist2(n)) * Taulist(l))
                  
                          fR(l,1,2)= real( dummy_k_tau(kx,ky,1,l,1,2) )*cos( aimag(iWnlist2(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,1,l,1,2) )*sin( aimag(iWnlist2(n)) * Taulist(l) )
                          fIm(l,1,2)= real( dummy_k_tau(kx,ky,1,l,1,2) )*sin( aimag(iWnlist2(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,1,l,1,2) )*cos( aimag(iWnlist2(n)) * Taulist(l))
                  
                  
                          fR(l,2,1)= real( dummy_k_tau(kx,ky,1,l,2,1) )*cos( aimag(iWnlist2(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,1,l,2,1) )*sin( aimag(iWnlist2(n)) * Taulist(l) )
                          fIm(l,2,1)= real( dummy_k_tau(kx,ky,1,l,2,1) )*sin( aimag(iWnlist2(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,1,l,2,1) )*cos( aimag(iWnlist2(n)) * Taulist(l))
                  
                  
                          fR(l,2,2)= real( dummy_k_tau(kx,ky,1,l,2,2) )*cos( aimag(iWnlist2(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,ky,1,l,2,2) )*sin( aimag(iWnlist2(n)) * Taulist(l) )
                          fIm(l,2,2)= real( dummy_k_tau(kx,ky,1,l,2,2) )*sin( aimag(iWnlist2(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,ky,1,l,2,2) )*cos( aimag(iWnlist2(n)) * Taulist(l))
                  
                  
                  
                     end do
                  
                     call fderiv(-1,TauGridNum,Taulist,fR(:,1,1),gR(:,1,1))
                     call fderiv(-1,TauGridNum,Taulist,fIm(:,1,1),gIm(:,1,1))
                     call fderiv(-1,TauGridNum,Taulist,fR(:,1,2),gR(:,1,2))
                     call fderiv(-1,TauGridNum,Taulist,fIm(:,1,2),gIm(:,1,2))
                     call fderiv(-1,TauGridNum,Taulist,fR(:,2,1),gR(:,2,1))
                     call fderiv(-1,TauGridNum,Taulist,fIm(:,2,1),gIm(:,2,1))
                     call fderiv(-1,TauGridNum,Taulist,fR(:,2,2),gR(:,2,2))
                     call fderiv(-1,TauGridNum,Taulist,fIm(:,2,2),gIm(:,2,2))
                  
                     dummy_k_iWn2(kx,ky,1,n,:,:)=gR(TauGridNum,:,:)+gIm(TauGridNum,:,:)*(0.0d0,1.0d0)




               end do
          end do
     end do
endif

if (Dimen==1) then
     do n=1,(iWnGridNum+1) !iWn loop
          do kx=1,(K1GridNum) ! kx loop
               do l=1,TauGridNum ! time integration for FT

                    fR(l,1,1)= real( dummy_k_tau(kx,1,1,l,1,1) )*cos( aimag(iWnlist2(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,1,1,l,1,1) )*sin( aimag(iWnlist2(n)) * Taulist(l) )
                    fIm(l,1,1)= real( dummy_k_tau(kx,1,1,l,1,1) )*sin( aimag(iWnlist2(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,1,1,l,1,1) )*cos( aimag(iWnlist2(n)) * Taulist(l))

                    fR(l,1,2)= real( dummy_k_tau(kx,1,1,l,1,2) )*cos( aimag(iWnlist2(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,1,1,l,1,2) )*sin( aimag(iWnlist2(n)) * Taulist(l) )
                    fIm(l,1,2)= real( dummy_k_tau(kx,1,1,l,1,2) )*sin( aimag(iWnlist2(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,1,1,l,1,2) )*cos( aimag(iWnlist2(n)) * Taulist(l))


                    fR(l,2,1)= real( dummy_k_tau(kx,1,1,l,2,1) )*cos( aimag(iWnlist2(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,1,1,l,2,1) )*sin( aimag(iWnlist2(n)) * Taulist(l) )
                    fIm(l,2,1)= real( dummy_k_tau(kx,1,1,l,2,1) )*sin( aimag(iWnlist2(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,1,1,l,2,1) )*cos( aimag(iWnlist2(n)) * Taulist(l))


                    fR(l,2,2)= real( dummy_k_tau(kx,1,1,l,2,2) )*cos( aimag(iWnlist2(n)) * Taulist(l) ) - aimag( dummy_k_tau(kx,1,1,l,2,2) )*sin( aimag(iWnlist2(n)) * Taulist(l) )
                    fIm(l,2,2)= real( dummy_k_tau(kx,1,1,l,2,2) )*sin( aimag(iWnlist2(n)) * Taulist(l) ) + aimag( dummy_k_tau(kx,1,1,l,2,2) )*cos( aimag(iWnlist2(n)) * Taulist(l))



               end do

               call fderiv(-1,TauGridNum,Taulist,fR(:,1,1),gR(:,1,1))
               call fderiv(-1,TauGridNum,Taulist,fIm(:,1,1),gIm(:,1,1))
               call fderiv(-1,TauGridNum,Taulist,fR(:,1,2),gR(:,1,2))
               call fderiv(-1,TauGridNum,Taulist,fIm(:,1,2),gIm(:,1,2))
               call fderiv(-1,TauGridNum,Taulist,fR(:,2,1),gR(:,2,1))
               call fderiv(-1,TauGridNum,Taulist,fIm(:,2,1),gIm(:,2,1))
               call fderiv(-1,TauGridNum,Taulist,fR(:,2,2),gR(:,2,2))
               call fderiv(-1,TauGridNum,Taulist,fIm(:,2,2),gIm(:,2,2))

               dummy_k_iWn2(kx,1,1,n,:,:)=gR(TauGridNum,:,:)+gIm(TauGridNum,:,:)*(0.0d0,1.0d0)




          end do
     end do
endif



return
end subroutine














!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for fourier transformation from k_tau to R_tau !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine k_tau2R_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_R_tau)

Implicit none
! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim
real(8), intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3), Taulist(TauGridNum), R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3),beta
complex(8), intent(in) ::  dummy_k_tau(K1GridNum,K2GridNum,K3GridNum,TauGridNum,MatDim,MatDim)
complex(8), intent(inout) :: dummy_R_tau(R1GridNum,R2GridNum,R3GridNum,TauGridNum,MatDim,MatDim)

! local variable
integer l



if (Dimen==3) then
      do l=1,TauGridNum ! time loop
           call fftw3d(dummy_k_tau(:,:,:,l,1,1),dummy_R_tau(:,:,:,l,1,1),K1GridNum,K2GridNum,K3GridNum,-1)
           call fftw3d(dummy_k_tau(:,:,:,l,1,2),dummy_R_tau(:,:,:,l,1,2),K1GridNum,K2GridNum,K3GridNum,-1)
           call fftw3d(dummy_k_tau(:,:,:,l,2,1),dummy_R_tau(:,:,:,l,2,1),K1GridNum,K2GridNum,K3GridNum,-1)
           call fftw3d(dummy_k_tau(:,:,:,l,2,2),dummy_R_tau(:,:,:,l,2,2),K1GridNum,K2GridNum,K3GridNum,-1)
     end do
endif


if (Dimen==2) then
      do l=1,TauGridNum ! time loop
           call fftw2d(dummy_k_tau(:,:,1,l,1,1),dummy_R_tau(:,:,1,l,1,1),K1GridNum,K2GridNum,-1)
           call fftw2d(dummy_k_tau(:,:,1,l,1,2),dummy_R_tau(:,:,1,l,1,2),K1GridNum,K2GridNum,-1)
           call fftw2d(dummy_k_tau(:,:,1,l,2,1),dummy_R_tau(:,:,1,l,2,1),K1GridNum,K2GridNum,-1)
           call fftw2d(dummy_k_tau(:,:,1,l,2,2),dummy_R_tau(:,:,1,l,2,2),K1GridNum,K2GridNum,-1)

     end do
endif

if (Dimen==1) then
      do l=1,TauGridNum ! time loop
           call fftw(dummy_k_tau(:,1,1,l,1,1),dummy_R_tau(:,1,1,l,1,1),K1GridNum,-1)
           call fftw(dummy_k_tau(:,1,1,l,1,2),dummy_R_tau(:,1,1,l,1,2),K1GridNum,-1)
           call fftw(dummy_k_tau(:,1,1,l,2,1),dummy_R_tau(:,1,1,l,2,1),K1GridNum,-1)
           call fftw(dummy_k_tau(:,1,1,l,2,2),dummy_R_tau(:,1,1,l,2,2),K1GridNum,-1)

     end do
endif




return
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for fourier transformation from R_tau to k_tau !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine R_tau2k_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_R_tau,dummy_k_tau)

Implicit none
! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim
real(8), intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3), Taulist(TauGridNum), R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3),beta
complex(8), intent(in) :: iWnlist(iWnGridNum), dummy_R_tau(R1GridNum,R2GridNum,R3GridNum,TauGridNum,MatDim,MatDim)
complex(8), intent(inout) :: dummy_k_tau(K1GridNum,K2GridNum,K3GridNum,TauGridNum,MatDim,MatDim)

! local variable
integer l

if (Dimen==3) then
     do l=1,TauGridNum ! time loop
           call fftw3d(dummy_R_tau(:,:,:,l,1,1),dummy_k_tau(:,:,:,l,1,1),R1GridNum,R2GridNum,R3GridNum,1)
           call fftw3d(dummy_R_tau(:,:,:,l,1,2),dummy_k_tau(:,:,:,l,1,2),R1GridNum,R2GridNum,R3GridNum,1)
           call fftw3d(dummy_R_tau(:,:,:,l,2,1),dummy_k_tau(:,:,:,l,2,1),R1GridNum,R2GridNum,R3GridNum,1)
           call fftw3d(dummy_R_tau(:,:,:,l,2,2),dummy_k_tau(:,:,:,l,2,2),R1GridNum,R2GridNum,R3GRidNum,1)
     end do
endif


if (Dimen==2) then
     do l=1,TauGridNum ! time loop
           call fftw2d(dummy_R_tau(:,:,1,l,1,1),dummy_k_tau(:,:,1,l,1,1),R1GridNum,R2GridNum,1) 
           call fftw2d(dummy_R_tau(:,:,1,l,1,2),dummy_k_tau(:,:,1,l,1,2),R1GridNum,R2GridNum,1) 
           call fftw2d(dummy_R_tau(:,:,1,l,2,1),dummy_k_tau(:,:,1,l,2,1),R1GridNum,R2GridNum,1)
           call fftw2d(dummy_R_tau(:,:,1,l,2,2),dummy_k_tau(:,:,1,l,2,2),R1GridNum,R2GridNum,1)
     end do
endif


if (Dimen==1) then
     do l=1,TauGridNum ! time loop
           call fftw(dummy_R_tau(:,1,1,l,1,1),dummy_k_tau(:,1,1,l,1,1),R1GridNum,1) 
           call fftw(dummy_R_tau(:,1,1,l,1,2),dummy_k_tau(:,1,1,l,1,2),R1GridNum,1) 
           call fftw(dummy_R_tau(:,1,1,l,2,1),dummy_k_tau(:,1,1,l,2,1),R1GridNum,1)
           call fftw(dummy_R_tau(:,1,1,l,2,2),dummy_k_tau(:,1,1,l,2,2),R1GridNum,1)
     end do
endif




return
end subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for fast fourier transform !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine fftw3d(inarray,outarray,xGridNum,yGridNum,zGridNum,direction)
!argument
integer,intent(in) :: xGridNum,yGridNum,zGridNum,direction
complex(8), intent(in) :: inarray(xGridNum,yGridNum,zGridNum)
complex(8), intent(inout) :: outarray(xGridNum,yGridNum,zGridNum)

!local variabe
integer*8 plan

if (direction .ne. 1) then   ! k->i
     call dfftw_plan_dft_3d(plan,xGridNum,yGridNum,zGridNum,inarray,outarray,-1,FFTW_ESTIMATE)
     call dfftw_execute_dft(plan,inarray, outarray)
     outarray=outarray/dble(xGridNum)/dble(yGridNum)/dble(zGridNum)
else    ! i->k
     call dfftw_plan_dft_3d(plan,xGridNum,yGridNum,zGridNum,inarray,outarray,1,FFTW_ESTIMATE)
     call dfftw_execute_dft(plan,inarray, outarray)
endif

call dfftw_destroy_plan(plan)
!outarray=cmplx(outarray)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for fast fourier transform !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine fftw2d(inarray,outarray,xGridNum,yGridNum,direction)
!argument
integer,intent(in) :: xGridNum,yGridNum,direction
complex(8), intent(in) :: inarray(xGridNum,yGridNum)
complex(8), intent(inout) :: outarray(xGridNum,yGridNum)

!local variabe
integer*8 plan

if (direction .ne. 1) then   ! k->i
     call dfftw_plan_dft_2d(plan,xGridNum,yGridNum,inarray,outarray,-1,FFTW_ESTIMATE)
     call dfftw_execute_dft(plan,inarray, outarray)
     outarray=outarray/dble(xGridNum)/dble(yGridNum)
else    ! i->k
     call dfftw_plan_dft_2d(plan,xGridNum,yGridNum,inarray,outarray,1,FFTW_ESTIMATE)
     call dfftw_execute_dft(plan,inarray, outarray)
endif

call dfftw_destroy_plan(plan)
!outarray=cmplx(outarray)

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for fast fourier transform !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine fftw(inarray,outarray,GridNum,direction)
!argument
integer,intent(in) :: GridNum,direction
complex(8), intent(in) :: inarray(GridNum)
complex(8), intent(inout) :: outarray(GridNum)

!local variabe
integer*8 plan

if (direction .ne. 1) then   ! k->i
     call dfftw_plan_dft_1d(plan,GridNum,inarray,outarray,-1,FFTW_ESTIMATE)
     call dfftw_execute_dft(plan,inarray, outarray)
     outarray=outarray/dble(GridNum)
else    ! i->k
     call dfftw_plan_dft_1d(plan,GridNum,inarray,outarray,1,FFTW_ESTIMATE)
     call dfftw_execute_dft(plan,inarray, outarray)
endif

call dfftw_destroy_plan(plan)
!outarray=cmplx(outarray)

end subroutine





! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: fderiv
! !INTERFACE:
subroutine fderiv(m,n,x,f,g)
! !INPUT/OUTPUT PARAMETERS:
!   m : order of derivative (in,integer)
!   n : number of points (in,integer)
!   x : abscissa array (in,real(n))
!   f : function array (in,real(n))
!   g : (anti-)derivative of f (out,real(n))
! !DESCRIPTION:
!   Given function $f$ defined on a set of points $x_i$ then if $m\ge 0$ this
!   routine computes the $m$th derivative of $f$ at each point. If $m<0$ the
!   anti-derivative of $f$ given by
!   $$ g(x_i)=\int_{x_1}^{x_i} f(x)\,dx $$
!   is calculated. If $m=-1$ then an accurate integral is computed by fitting
!   the function to a clamped cubic spline. When $m=-3$ the fast but low
!   accuracy trapezoidal integration method is used. Simpson's integration,
!   which is slower but more accurate than the trapezoidal method, is used if
!   $m=-2$.
!
! !REVISION HISTORY:
!   Created May 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: m,n
real(8), intent(in) :: x(n),f(n)
real(8), intent(out) :: g(n)
! local variables
integer i
real(8) x0,x1,x2,dx
! automatic arrays
real(8) cf(3,n)
if (n.le.0) then
  write(*,*)
  write(*,'("Error(fderiv): invalid number of points : ",I8)') n
  write(*,*)
  stop
end if
select case(m)
case(-3)
! low accuracy trapezoidal integration
  g(1)=0.d0
  do i=1,n-1
    g(i+1)=g(i)+0.5d0*(x(i+1)-x(i))*(f(i+1)+f(i))
  end do
  return
case(-2)
! medium accuracy Simpson integration
  g(1)=0.d0
  do i=1,n-2
    x0=x(i)
    x1=x(i+1)
    x2=x(i+2)
    g(i+1)=g(i)+(x0-x1)*(f(i+2)*(x0-x1)**2+f(i+1)*(x2-x0)*(x0+2.d0*x1-3.d0*x2) &
     +f(i)*(x2-x1)*(2.d0*x0+x1-3.d0*x2))/(6.d0*(x0-x2)*(x1-x2))
  end do
  x0=x(n)
  x1=x(n-1)
  x2=x(n-2)
  g(n)=g(n-1)+(x1-x0)*(f(n-2)*(x1-x0)**2+f(n)*(x1-x2)*(3.d0*x2-x1-2.d0*x0) &
   +f(n-1)*(x0-x2)*(3.d0*x2-2.d0*x1-x0))/(6.d0*(x2-x1)*(x2-x0))
  return
case(0)
  g(:)=f(:)
  return
case(4:)
  g(:)=0.d0
  return
end select
! high accuracy integration/differentiation from spline interpolation
call spline(n,x,f,cf)
select case(m)
case(:-1)
  g(1)=0.d0
  do i=1,n-1
    dx=x(i+1)-x(i)
    g(i+1)=g(i)+(((0.25d0*cf(3,i)*dx+0.3333333333333333333d0*cf(2,i))*dx &
     +0.5d0*cf(1,i))*dx+f(i))*dx
  end do
case(1)
  g(:)=cf(1,:)
case(2)
  g(:)=2.d0*cf(2,:)
case(3)
  g(:)=6.d0*cf(3,:)
end select
return
end subroutine
!EOC

! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: spline
! !INTERFACE:
subroutine spline(n,x,f,cf)
! !INPUT/OUTPUT PARAMETERS:
!   n  : number of points (in,integer)
!   x  : abscissa array (in,real(n))
!   f  : input data array (in,real(n))
!   cf : cubic spline coefficients (out,real(3,n))
! !DESCRIPTION:
!   Calculates the coefficients of a cubic spline fitted to input data. In other
!   words, given a set of data points $f_i$ defined at $x_i$, where
!   $i=1\ldots n$, the coefficients $c_j^i$ are determined such that
!   $$ y_i(x)=f_i+c_1^i(x-x_i)+c_2^i(x-x_i)^2+c_3^i(x-x_i)^3, $$
!   is the interpolating function for $x\in[x_i,x_{i+1})$. The coefficients are
!   determined piecewise by fitting a cubic polynomial to adjacent points.
!
! !REVISION HISTORY:
!   Created November 2011 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: x(n),f(n)
real(8), intent(out) :: cf(3,n)
! local variables
integer i
real(8) x0,x1,x2,x3,y0,y1,y2,y3
real(8) c1,c2,c3,t0,t1,t2,t3,t4,t5,t6
if (n.le.0) then
  write(*,*)
  write(*,'("Error(spline): n <= 0 : ",I8)') n
  write(*,*)
  stop
end if
if (n.eq.1) then
  cf(:,1)=0.d0
  return
end if
if (n.eq.2) then
  cf(1,1)=(f(2)-f(1))/(x(2)-x(1))
  cf(2:3,1)=0.d0
  cf(1,2)=cf(1,1)
  cf(2:3,2)=0.d0
  return
end if
if (n.eq.3) then
  x0=x(1)
  x1=x(2)-x0
  x2=x(3)-x0
  y0=f(1)
  y1=f(2)-y0
  y2=f(3)-y0
  t0=1.d0/(x1*x2*(x2-x1))
  t1=x1*y2
  t2=x2*y1
  c1=t0*(x2*t2-x1*t1)
  c2=t0*(t1-t2)
  cf(1,1)=c1
  cf(2,1)=c2
  cf(3,1)=0.d0
  t3=2.d0*c2
  cf(1,2)=c1+t3*x1
  cf(2,2)=c2
  cf(3,2)=0.d0
  cf(1,3)=c1+t3*x2
  cf(2,3)=c2
  cf(3,3)=0.d0
  return
end if
y0=f(1)
y1=f(2)-y0
y2=f(3)-y0
y3=f(4)-y0
x0=x(1)
x1=x(2)-x0
x2=x(3)-x0
x3=x(4)-x0
t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
t1=x1*x2*y3
t2=x2*x3*y1
t3=x3*x1*y2
t4=x1**2
t5=x2**2
t6=x3**2
y1=t3*t6-t1*t5
y3=t2*t5-t3*t4
y2=t1*t4-t2*t6
c1=t0*(x1*y1+x2*y2+x3*y3)
c2=-t0*(y1+y2+y3)
c3=t0*(t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1))
cf(1,1)=c1
cf(2,1)=c2
cf(3,1)=c3
cf(1,2)=c1+2.d0*c2*x1+3.d0*c3*t4
cf(2,2)=c2+3.d0*c3*x1
cf(3,2)=c3
if (n.eq.4) then
  cf(1,3)=c1+2.d0*c2*x2+3.d0*c3*t5
  cf(2,3)=c2+3.d0*c3*x2
  cf(3,3)=c3
  cf(1,4)=c1+2.d0*c2*x3+3.d0*c3*t6
  cf(2,4)=c2+3.d0*c3*x3
  cf(3,4)=c3
  return
end if
do i=3,n-2
  y0=f(i)
  y1=f(i-1)-y0
  y2=f(i+1)-y0
  y3=f(i+2)-y0
  x0=x(i)
  x1=x(i-1)-x0
  x2=x(i+1)-x0
  x3=x(i+2)-x0
  t1=x1*x2*y3
  t2=x2*x3*y1
  t3=x3*x1*y2
  t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
  c3=t0*(t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1))
  t4=x1**2
  t5=x2**2
  t6=x3**2
  y1=t3*t6-t1*t5
  y2=t1*t4-t2*t6
  y3=t2*t5-t3*t4
  cf(1,i)=t0*(x1*y1+x2*y2+x3*y3)
  cf(2,i)=-t0*(y1+y2+y3)
  cf(3,i)=c3
end do
c1=cf(1,n-2)
c2=cf(2,n-2)
c3=cf(3,n-2)
cf(1,n-1)=c1+2.d0*c2*x2+3.d0*c3*t5
cf(2,n-1)=c2+3.d0*c3*x2
cf(3,n-1)=c3
cf(1,n)=c1+2.d0*c2*x3+3.d0*c3*t6
cf(2,n)=c2+3.d0*c3*x3
cf(3,n)=c3
return
end subroutine
!EOC

