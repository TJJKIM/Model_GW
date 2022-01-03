!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating local green function !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CallocG(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,dummy_k_iWn,locGreenf)
Implicit none

! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,Dimen
complex(8), intent(in) :: dummy_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)
complex(8), intent(inout) :: locGreenf(iWnGridNum,MatDim,MatDim)

!local variable
integer kx,ky,kz,n
complex(8) dummy(MatDim,MatDim)

print *," "
print *,"Calculating local Greenfunction Gloc(iWn)..."

if (Dimen==3) then
     do n=1,iWnGridNum
          dummy(:,:)=(0.0,0.0)
          do kx=1,(K1GridNum)
               do ky=1,(K2GridNum)
                    do kz=1,(K3GridNum)
                         dummy=dummy+dummy_k_iWn(kx,ky,kz,n,:,:)
                    enddo
               enddo
          enddo
          locGreenf(n,:,:)=dummy/(K1GridNum)/(K2GridNum)/(K3GridNum)
     enddo
endif

if (Dimen==2) then
     do n=1,iWnGridNum
          dummy(:,:)=(0.0,0.0)
          do kx=1,(K1GridNum)
               do ky=1,(K2GridNum)
                    dummy=dummy+dummy_k_iWn(kx,ky,1,n,:,:)
               enddo
          enddo
          locGreenf(n,:,:)=dummy/(K1GridNum)/(K2GridNum)
     enddo
endif


if (Dimen==1) then
     do n=1,iWnGridNum
          dummy(:,:)=(0.0,0.0)
          do kx=1,(K1GridNum)
               dummy=dummy+dummy_k_iWn(kx,1,1,n,:,:)
          enddo
          locGreenf(n,:,:)=dummy/(K1GridNum)
     enddo
endif


return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating local slef-energy  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalocSelfE(Dimen,MatDim,GWSelfE_k_iWn,locSelfE,K1GridNum,K2GridNum,K3GridNum,iWnGridNum)
Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum,iWnGridNum, Dimen,MatDim
complex(8), intent(in) :: GWSelfE_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)
complex(8), intent(inout) :: locSelfE(iWnGridNum,MatDim,MatDim)
! local variable
integer kx,ky,kz,n,Z

locSelfE(:,:,:)=(0.0d0,0.0d0) ! initialize

if(Dimen==3) then
     do kx=1,(K1GridNum) !k space loop
          do ky=1,(K2GridNum) !k space loop
               do kz=1,(K3GridNum) !k space loop
                    do n=1,(iWnGridNum)
                          locSelfE(n,:,:)= locSelfE(n,:,:)+ ( (1.0d0/(K1GridNum*K2GridNum*K3GridNum))*( GWSelfE_k_iWn(kx,ky,kz,n,:,:)  )  )
                    enddo
               enddo
          enddo
     enddo
endif


if(Dimen==2) then
     do kx=1,(K1GridNum) !k space loop
          do ky=1,(K2GridNum) !k space loop
              do n=1,(iWnGridNum)
                    locSelfE(n,:,:)= locSelfE(n,:,:)+ ( (1.0d0/(K1GridNum*K2GridNum))*( GWSelfE_k_iWn(kx,ky,1,n,:,:)  )  )
              enddo
          enddo
     enddo
endif


if(Dimen==1) then
     do kx=1,(K1GridNum) !k space loop
         do n=1,(iWnGridNum)
               locSelfE(n,:,:)= locSelfE(n,:,:)+ ( (1.0d0/(K1GridNum))*( GWSelfE_k_iWn(kx,1,1,n,:,:)  )  )
         enddo
     enddo
endif



return
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating local polarization !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalocPol(Dimen,MatDim,GWPol_k_iWn,locPol,K1GridNum,K2GridNum,K3GridNum,iWnGridNum)
Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum,iWnGridNum,Dimen,MatDim
complex(8), intent(in) :: GWPol_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim)
complex(8), intent(inout) :: locPol(iWnGridNum+1,MatDim,MatDim)
! local variable
integer kx,ky,kz,n,Z

locPol(:,:,:)=(0.0d0,0.0d0) ! initialize


if (Dimen==3) then
     do kx=1,(K1GridNum) !k space loop
          do ky=1,(K2GridNum) !k space loop
               do kz=1,(K3GridNum) !k space loop
                    do n=1,(iWnGridNum+1)
                         locPol(n,:,:)= locPol(n,:,:)+ ( (1.0d0/(K1GridNum*K2GridNum*K3GridNum))*(  GWPol_k_iWn(kx,ky,kz,n,:,:) ) )
                    enddo
               enddo
          enddo
     enddo
endif


if (Dimen==2) then
     do kx=1,(K1GridNum) !k space loop
          do ky=1,(K2GridNum) !k space loop
              do n=1,(iWnGridNum+1)
                   locPol(n,:,:)= locPol(n,:,:)+ ( (1.0d0/(K1GridNum*K2GridNum))*(  GWPol_k_iWn(kx,ky,1,n,:,:) ) )
              enddo
          enddo
     enddo
endif


if (Dimen==1) then
     do kx=1,(K1GridNum) !k space loop
         do n=1,(iWnGridNum+1)
              locPol(n,:,:)= locPol(n,:,:)+ ( (1.0d0/(K1GridNum))*(  GWPol_k_iWn(kx,1,1,n,:,:) ) )
         enddo
     enddo
endif


return
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating local green function !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CallocW(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,dummy_k_iWn,locW)
Implicit none

! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim
complex(8), intent(in) :: dummy_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim)
complex(8), intent(inout) :: locW(iWnGridNum+1,MatDim,MatDim)

!local variable
integer kx,ky,kz,n
complex(8) dummy(MatDim,MatDim)

print *," "
print *,"Calculating local Greenfunction Gloc(iWn)..."

if (Dimen==3) then
     do n=1,(iWnGridNum+1)
          dummy(:,:)=(0.0,0.0)
          do kx=1,(K1GridNum)
               do ky=1,(K2GridNum)
                    do kz=1,(K3GridNum)
                         dummy=dummy+dummy_k_iWn(kx,ky,kz,n,:,:)
                    enddo
               enddo
          enddo
          locW(n,:,:)=dummy/(K1GridNum)/(K2GridNum)/(K3GridNum)
     enddo
endif 


if (Dimen==2) then
     do n=1,(iWnGridNum+1)
          dummy(:,:)=(0.0,0.0)
          do kx=1,(K1GridNum)
               do ky=1,(K2GridNum)
                    dummy=dummy+dummy_k_iWn(kx,ky,1,n,:,:)
               enddo
          enddo
          locW(n,:,:)=dummy/(K1GridNum)/(K2GridNum)
     enddo
endif


if (Dimen==1) then
     do n=1,(iWnGridNum+1)
          dummy(:,:)=(0.0,0.0)
          do kx=1,(K1GridNum)
               dummy=dummy+dummy_k_iWn(kx,ky,kz,n,:,:)
          enddo
          locW(n,:,:)=dummy/(K1GridNum)
     enddo
endif





return
end subroutine

