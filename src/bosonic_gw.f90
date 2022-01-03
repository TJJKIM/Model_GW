!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating polarization !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine PolCal(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_R_tau,dummy_R_tauInv,Polarization_R_tau)


Implicit none
! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim
real(8), intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3), Taulist(TauGridNum), R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3),beta
complex(8), intent(in) :: iWnlist(iWnGridNum), dummy_R_tau(R1GridNum,R2GridNum,R3GridNum,TauGridNum,MatDim,MatDim), dummy_R_tauInv(R1GridNum,R2GridNum,R3GridNum,TauGridNum,MatDim,MatDim)
complex(8), intent(inout) :: Polarization_R_tau(R1GridNum,R2GridNum,R3GridNum,TauGridNum,MatDim,MatDim)

! local variable
integer ix,iy,iz,j,k,l,E,n,o1,o2
complex(8) Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim),dummy3(MatDim,MatDim)

Polarization_R_tau=(0.0d0,0.0d0)

if (Dimen==3) then
     do ix=1,(R1GridNum) ! real Rx space site loop
          do iy=1,(R2GridNum) ! real Ry space site loop
               do iz=1,(R3GridNum) ! real Rz space site loop
                    do l=1,(TauGridNum) ! time space loop
                         do o1=1,MatDim
                              do o2=1,MatDim
                                   Polarization_R_tau(ix,iy,iz,l,o1,o2)=-(1.0d0)*dummy_R_tau(ix,iy,iz,l,o1,o2)*dummy_R_tauInv(ix,iy,iz,TauGridNum+1-l,o2,o1)
                                   !Polarization_R_tau(ix,iy,iz,l,1,2)=-(1.0d0)*dummy_R_tau(ix,iy,iz,l,1,2)*dummy_R_tauInv(ix,iy,iz,TauGridNum+1-l,2,1)
                                   !Polarization_R_tau(ix,iy,iz,l,2,1)=-(1.0d0)*dummy_R_tau(ix,iy,iz,l,2,1)*dummy_R_tauInv(ix,iy,iz,TauGridNum+1-l,1,2)
                                   !Polarization_R_tau(ix,iy,iz,l,2,2)=-(1.0d0)*dummy_R_tau(ix,iy,iz,l,2,2)*dummy_R_tauInv(ix,iy,iz,TauGridNum+1-l,2,2)
                              enddo
                         enddo
                    enddo
               enddo
          enddo
     enddo
endif



if (Dimen==2) then
     do ix=1,(R1GridNum) ! real Rx space site loop
          do iy=1,(R2GridNum) ! real Ry space site loop
              do l=1,(TauGridNum) ! time space loop
                   do o1=1,(MatDim)
                        do o2=1,(MatDim)
                              Polarization_R_tau(ix,iy,1,l,o1,o2)=-(1.0d0)*dummy_R_tau(ix,iy,1,l,o1,o2)*dummy_R_tauInv(ix,iy,1,TauGridNum+1-l,o2,o1)
                              !Polarization_R_tau(ix,iy,1,l,1,2)=-(1.0d0)*dummy_R_tau(ix,iy,1,l,1,2)*dummy_R_tauInv(ix,iy,1,TauGridNum+1-l,2,1)
                              !Polarization_R_tau(ix,iy,1,l,2,1)=-(1.0d0)*dummy_R_tau(ix,iy,1,l,2,1)*dummy_R_tauInv(ix,iy,1,TauGridNum+1-l,1,2)
                              !Polarization_R_tau(ix,iy,1,l,2,2)=-(1.0d0)*dummy_R_tau(ix,iy,1,l,2,2)*dummy_R_tauInv(ix,iy,1,TauGridNum+1-l,2,2)
                        enddo 
                   enddo
              enddo
          enddo
     enddo
endif



if (Dimen==1) then
     do ix=1,(R1GridNum) ! real Rx space site loop
         do l=1,(TauGridNum) ! time space loop
              do o1=1,(MatDim)
                   do o2=1,(MatDim)
                        Polarization_R_tau(ix,1,1,l,o1,o2)=-(1.0d0)*dummy_R_tau(ix,1,1,l,o1,o2)*dummy_R_tauInv(ix,1,1,TauGridNum+1-l,o2,o1)
                        !Polarization_R_tau(ix,1,1,l,1,2)=-(1.0d0)*dummy_R_tau(ix,1,1,l,1,2)*dummy_R_tauInv(ix,1,1,TauGridNum+1-l,2,1)
                        !Polarization_R_tau(ix,1,1,l,2,1)=-(1.0d0)*dummy_R_tau(ix,1,1,l,2,1)*dummy_R_tauInv(ix,1,1,TauGridNum+1-l,1,2)
                        !Polarization_R_tau(ix,1,1,l,2,2)=-(1.0d0)*dummy_R_tau(ix,1,1,l,2,2)*dummy_R_tauInv(ix,1,1,TauGridNum+1-l,2,2)
                   enddo
              enddo
         enddo
     enddo
endif



return
end subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating screend coulomb interactin W !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine WcorrCal(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,U,Vq,MatDim,beta,K1list,K2list,K3list,iWnlist2,Taulist,R1list,R2list,R3list,Polarization_k_iWn,ScreenedW_k_iWn)


Implicit none
! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim
real(8), intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3) ,K3list(K3GridNum,3),Taulist(TauGridNum), R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3),beta,U(MatDim,MatDim)
complex(8), intent(in) :: iWnlist2(iWnGridNum+1), Polarization_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim), Vq(K1GridNum,K2GridNum,K3GridNum,MatDim,MatDim)
complex(8), intent(inout) :: ScreenedW_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim)

! local variable
integer i,j,kx,ky,kz,l,E,n,o1,o2
complex(8) PolSum,Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim),dummy3(MatDim,MatDim)

ScreenedW_k_iWn = (0.0d0,0.0d0)

!call Identity(MatDim,Iden)
if (Dimen==3) then
     do kx =1,(K1GridNum) ! kx space loop
          do ky =1,(K2GridNum) ! ky space loop
               do kz =1,(K3GridNum) ! kz space loop
                    do n=1,(iWnGridNum+1) !  iWn loop
                         PolSum=sum(Polarization_k_iWn(kx,ky,kz,n,:,:))
                         !print *, PolSum
                         do o1=1,(MatDim)
                              do o2=1,(MatDim)
      !                            dummy(:,:)=(0.0d0,0.0d0) !initialize
      !                            dummy2(:,:)=(0.0d0,0.0d0) !initialize
      !                            dummy(1,1)=1.0d0-U/(2.0d0)*( Polarization_k_iWn(k,n,1,1)+Polarization_k_iWn(k,n,2,1) )
      !                            dummy(1,2)=-U/(2.0d0)*( Polarization_k_iWn(k,n,1,2)+Polarization_k_iWn(k,n,2,2) )
      !                            dummy(2,1)=-U/(2.0d0)*( Polarization_k_iWn(k,n,1,1)+Polarization_k_iWn(k,n,2,1) )
      !                            dummy(2,2)=1.0d0-U/(2.0d0)*( Polarization_k_iWn(k,n,1,2)+Polarization_k_iWn(k,n,2,2) )
      !                            call MatInv(MatDim,dummy,dummy2)
      !                            Vq=2*V*( cos( Kxlist(kx)*ax )+cos( Kylist(ky)*ay )+cos( Kzlist(kz)*az ) )
                                  
                                   ScreenedW_k_iWn(kx,ky,kz,n,o1,o2)=(U(o1,o2)+Vq(kx,ky,kz,o1,o2))*(U(o1,o2)+Vq(kx,ky,kz,o1,o2))*PolSum / ( 1.0d0 - (U(o1,o2)+Vq(kx,ky,kz,o1,o2))*PolSum )
                   
 
                              enddo 
                         enddo 
                    enddo
               enddo
          enddo
     enddo
endif

if (Dimen==2) then
     do kx =1,(K1GridNum) ! kx space loop
          do ky =1,(K2GridNum) ! ky space loop
              do n=1,(iWnGridNum+1) !  iWn loop
                   PolSum=sum(Polarization_k_iWn(kx,ky,1,n,:,:))
                   !if (n==(iWnGridNum+1)/2) then
                   !      print *, PolSum 
                   !endif
                   do o1=1,(MatDim)
                        do o2=1,(MatDim)
      !                      dummy(:,:)=(0.0d0,0.0d0) !initialize
      !                      dummy2(:,:)=(0.0d0,0.0d0) !initialize
      !                      dummy(1,1)=1.0d0-U/(2.0d0)*( Polarization_k_iWn(k,n,1,1)+Polarization_k_iWn(k,n,2,1) )
      !                      dummy(1,2)=-U/(2.0d0)*( Polarization_k_iWn(k,n,1,2)+Polarization_k_iWn(k,n,2,2) )
      !                      dummy(2,1)=-U/(2.0d0)*( Polarization_k_iWn(k,n,1,1)+Polarization_k_iWn(k,n,2,1) )
      !                      dummy(2,2)=1.0d0-U/(2.0d0)*( Polarization_k_iWn(k,n,1,2)+Polarization_k_iWn(k,n,2,2) )
      !                      call MatInv(MatDim,dummy,dummy2)
      !                      Vq=2*V*( cos( Kxlist(kx)*ax )+cos( Kylist(ky)*ay )+cos( Kzlist(kz)*az ) )
          
                             ScreenedW_k_iWn(kx,ky,1,n,o1,o2)=(U(o1,o2)+Vq(kx,ky,1,o1,o2))*(U(o1,o2)+Vq(kx,ky,1,o1,o2))*PolSum / ( 1.0d0 - (U(o1,o2)+Vq(kx,ky,1,o1,o2))*PolSum )


                        enddo
                   enddo
              enddo
          enddo
     enddo
endif


if (Dimen==1) then
     do kx =1,(K1GridNum) ! kx space loop
         do n=1,(iWnGridNum+1) !  iWn loop
              PolSum=sum(Polarization_k_iWn(kx,1,1,n,:,:))
              do o1=1,(MatDim)
                   do o2=1,(MatDim)
      !                 dummy(:,:)=(0.0d0,0.0d0) !initialize
      !                 dummy2(:,:)=(0.0d0,0.0d0) !initialize
      !                 dummy(1,1)=1.0d0-U/(2.0d0)*( Polarization_k_iWn(k,n,1,1)+Polarization_k_iWn(k,n,2,1) )
      !                 dummy(1,2)=-U/(2.0d0)*( Polarization_k_iWn(k,n,1,2)+Polarization_k_iWn(k,n,2,2) )
      !                 dummy(2,1)=-U/(2.0d0)*( Polarization_k_iWn(k,n,1,1)+Polarization_k_iWn(k,n,2,1) )
      !                 dummy(2,2)=1.0d0-U/(2.0d0)*( Polarization_k_iWn(k,n,1,2)+Polarization_k_iWn(k,n,2,2) )
      !                 call MatInv(MatDim,dummy,dummy2)
      !                 Vq=2*V*( cos( Kxlist(kx)*ax )+cos( Kylist(ky)*ay )+cos( Kzlist(kz)*az ) )
        
                        ScreenedW_k_iWn(kx,1,1,n,o1,o2)=(U(o1,o2)+Vq(kx,1,1,o1,o2))*(U(o1,o2)+Vq(kx,1,1,o1,o2))*PolSum / ( 1.0d0 - (U(o1,o2)+Vq(kx,1,1,o1,o2))*PolSum )


                   enddo
              enddo
         enddo
     enddo
endif






return
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for Screend W sum check loop !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Corrsum(Dimen,ScreenedW_k_iWn,beta,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,U)
Implicit none
! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum, K3GridNum,iWnGridNum,MatDim
real(8), intent(in) :: beta,U(MatDim,MatDim)
complex(8), intent(in) :: ScreenedW_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim)

! local variable
integer kx,ky,kz,n,Z,i
complex(8) Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim),dummy3(MatDim,MatDim)

dummy(:,:)=(0.0d0,0.0d0)

if (Dimen ==3) then
     do kx =1,(K1GridNum) ! kx space loop
          do ky =1,(K2GridNum) ! ky space loop
               do kz =1,(K3GridNum) ! kz space loop
                    Z=(iWnGridNum/2+1)
                    dummy=dummy+( (1.0d0)/dble(K1GridNum)/dble(K2GridNum)/dble(K3GridNum) )*ScreenedW_k_iWn(kx,ky,kz,Z,:,:)
               enddo
          enddo
     enddo
endif

if (Dimen ==2) then
     do kx =1,(K1GridNum) ! kx space loop
          do ky =1,(K2GridNum) ! ky space loop
               Z=(iWnGridNum/2+1)
               dummy=dummy+( (1.0d0)/dble(K1GridNum)/dble(K2GridNum) )*ScreenedW_k_iWn(kx,ky,1,Z,:,:)
          enddo
     enddo
endif

if (Dimen ==1) then
     do kx =1,(K1GridNum) ! kx space loop
         Z=(iWnGridNum/2+1)
         !print *, ScreenedW_k_iWn(kx,1,1,Z,:,:)
         dummy=dummy+( (1.0d0)/dble(K1GridNum) )*ScreenedW_k_iWn(kx,1,1,Z,:,:)
     enddo
endif




print *, "Bare interaction : "
print *, "[Onsite]"
do i=1,MatDim
     print 15, U(i,1:MatDim)
     15 format(30(f8.4))
enddo
print *, " "
print *, "Screened Interaction :"
print *, "[Wc]"
do i=1,MatDim
     print 34, dummy(i,1:MatDim)
     34 format(30(f8.4))
enddo


return
end subroutine

