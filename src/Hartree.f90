!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating hartree self energy in (k,iWn) space !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Hartree_self(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,U,V0,MatDim,beta,delta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Greenf_k_iWn,Hartree_k_iWn)


Implicit none
! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim
real(8), intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3), Taulist(TauGridNum), R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3),beta,U(MatDim,MatDim),delta
complex(8), intent(in) :: iWnlist(iWnGridNum),V0(MatDim,MatDim), Greenf_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)
complex(8), intent(inout) :: Hartree_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)

! local variable
integer i,j,kx,ky,kz,l,E,n,m,o1,o2
complex(8) Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim),dummy3(MatDim,MatDim)
complex(8) gtilda2(MatDim,MatDim), gtilda3(MatDim,MatDim), MatGcut(MatDim,MatDim), asymtoticG(MatDim,MatDim)


call Identity(MatDim,Iden)
dummy(:,:)=(0.0d0,0.0d0) !initialize
Hartree_k_iWn(:,:,:,:,:,:)=(0.0d0,0.0d0) !initialize

if (Dimen==3) then
     do kx =1,(K1GridNum) ! kx space loop
          do ky=1,(K2GridNum) ! ky space loop
               do kz=1,(K3GridNum) ! ky space loop
                    gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( Greenf_k_iWn(kx,ky,kz,iWnGridNum,:,:)+Greenf_k_iWn(kx,ky,kz,1,:,:) )/(2.0)
                    gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( Greenf_k_iWn(kx,ky,kz,iWnGridNum,:,:)-Greenf_k_iWn(kx,ky,kz,1,:,:) )/(2.0) - Iden/iWnlist(iWnGridNum) )


                    do n=(int(iWnGridNum/2)+1),(iWnGridNum) !  iWn loop

                         dummy= dummy - ((1.0d0)/(beta))* ( Greenf_k_iWn(kx,ky,kz,n,:,:) - Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) - gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(-iWnlist(n)*(beta-delta))
                         dummy= dummy - ((1.0d0)/(beta))* ( Greenf_k_iWn(kx,ky,kz,iWnGridNum+1-n,:,:) + Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) + gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(iWnlist(n)*(beta-delta))


                    enddo
                    asymtoticG =  ( Iden/(2.0d0) ) - ( gtilda2*beta*( (beta-delta)/((2.0d0)*(beta))-(0.25d0)) ) + ( gtilda3*beta*beta/(4.0d0)*( (beta-delta)*(beta-delta)/(beta*beta) - (beta-delta)/beta) )
                    dummy=dummy+asymtoticG
               enddo
          enddo
     enddo
      
     !dummy=dummy/(K1GridNum)/(K2GridNum)/(K3GridNum) !test
     !print *, "Kpts: ",K1list(kx),K2list(ky),K3list(kz)
     !print *, dummy
     do o1=1,(MatDim)
          do o2=1,(MatDim)
               !if( U(o1,o2)/=0.0 ) then ! to find same position orbital. ( because onsite interaction almoste finite )
                    Hartree_k_iWn(:,:,:,:,o1,o1)=Hartree_k_iWn(:,:,:,:,o1,o1)+(U(o1,o2)+V0(o1,o2))*( dummy(o2,o2) )/(K1GridNum)/(K2GridNum)/(K3GridNum)
                    !Hartree_k_iWn(:,:,:,:,(MatDim/2)+o1,(MatDim/2)+o2)=Hartree_k_iWn(:,:,:,:,o1,o2)+(U(o1,o2)+V0(o1,o2))*( dummy(o2,o2) + dummy( (MatDim/2)+o2, (MatDim/2)+o2 ) )/(K1GridNum)/(K2GridNum)/(K3GridNum)
               !endif
          enddo
     enddo

endif




if (Dimen==2) then
     do kx =1,(K1GridNum) ! kx space loop
          do ky=1,(K2GridNum) ! ky space loop
               gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( Greenf_k_iWn(kx,ky,1,iWnGridNum,:,:)+Greenf_k_iWn(kx,ky,1,1,:,:) )/(2.0)
               gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( Greenf_k_iWn(kx,ky,1,iWnGridNum,:,:)-Greenf_k_iWn(kx,ky,1,1,:,:) )/(2.0) - Iden/iWnlist(iWnGridNum) )


               do n=(int(iWnGridNum/2)+1),(iWnGridNum) !  iWn loop

                    dummy= dummy - ((1.0d0)/(beta))* ( Greenf_k_iWn(kx,ky,1,n,:,:) - Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) - gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(-iWnlist(n)*(beta-delta))
                    dummy= dummy - ((1.0d0)/(beta))* ( Greenf_k_iWn(kx,ky,1,iWnGridNum+1-n,:,:) + Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) + gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(iWnlist(n)*(beta-delta))


               enddo
               asymtoticG =  ( Iden/(2.0d0) ) - ( gtilda2*beta*( (beta-delta)/((2.0d0)*(beta))-(0.25d0)) ) + ( gtilda3*beta*beta/(4.0d0)*( (beta-delta)*(beta-delta)/(beta*beta) - (beta-delta)/beta) )
               dummy=dummy+asymtoticG
          enddo
     enddo

     !dummy=dummy/(K1GridNum)/(K2GridNum) !test
     !print *, "Kpts: ",K1list(kx),K2list(ky),K3list(1)
     !print *, dummy
     !do o=1,(MatDim/2)
     !     Hartree_k_iWn(:,:,:,:,o,o)=(U+V0)*(dummy(1,1)+dummy(2,2))/(KxGridNum)/(KyGridNum)/(KzGridNum)
     !     Hartree_k_iWn(:,:,:,:,2,2)=(U+V0)*(dummy(1,1)+dummy(2,2))/(KxGridNum)/(KyGridNum)/(KzGridNum)
     !enddo

     do o1=1,(MatDim)
          do o2=1,(MatDim)
               !if( U(o1,o2)/=0.0 ) then ! to find same position orbital. ( because onsite interaction almoste finite )
                    Hartree_k_iWn(:,:,1,:,o1,o1)=Hartree_k_iWn(:,:,1,:,o1,o1)+(U(o1,o2)+V0(o1,o2))*( dummy(o2,o2) )/(K1GridNum)/(K2GridNum)
                    !Hartree_k_iWn(:,:,1,:,(MatDim/2)+o1,(MatDim/2)+o2)=Hartree_k_iWn(:,:,1,:,o1,o2)+(U(o1,o2)+V0(o1,o2))*( dummy(o2,o2) + dummy( (MatDim/2)+o2, (MatDim/2)+o2 ) )/(K1GridNum)/(K2GridNum)
               !endif
          enddo
     enddo




endif





if (Dimen==1) then
     do kx =1,(K1GridNum) ! kx space loop
          gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( Greenf_k_iWn(kx,1,1,iWnGridNum,:,:)+Greenf_k_iWn(kx,1,1,1,:,:) )/(2.0)
          gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( Greenf_k_iWn(kx,1,1,iWnGridNum,:,:)-Greenf_k_iWn(kx,1,1,1,:,:) )/(2.0) - Iden/iWnlist(iWnGridNum) )


          do n=(int(iWnGridNum/2)+1),(iWnGridNum) !  iWn loop

               dummy= dummy - ((1.0d0)/(beta))* ( Greenf_k_iWn(kx,1,1,n,:,:) - Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) - gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(-iWnlist(n)*(beta-delta))
               dummy= dummy - ((1.0d0)/(beta))* ( Greenf_k_iWn(kx,1,1,iWnGridNum+1-n,:,:) + Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) + gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(iWnlist(n)*(beta-delta))


          enddo
          asymtoticG =  ( Iden/(2.0d0) ) - ( gtilda2*beta*( (beta-delta)/((2.0d0)*(beta))-(0.25d0)) ) + ( gtilda3*beta*beta/(4.0d0)*( (beta-delta)*(beta-delta)/(beta*beta) - (beta-delta)/beta) )
          dummy=dummy+asymtoticG
     enddo
     
     !dummy=dummy/(K1GridNum) !test
     !kx=kx-1
     !print *,"kx :", kx
     !print *, "Kpts: ",K1list(kx),K2list(1),K3list(1)
     !print *, "Occu(1,1:4):", dummy(1,:)
     !print *, "Occu(2:1:4):", dummy(2,:)
     !print *, "Occu(3:1:4):", dummy(3,:)
     !print *, "Occu(4:1:4):", dummy(4,:)
     !do o=1,(MatDim/2)
     !     Hartree_k_iWn(:,:,:,:,o,o)=(U+V0)*(dummy(1,1)+dummy(2,2))/(KxGridNum)/(KyGridNum)/(KzGridNum)
     !     Hartree_k_iWn(:,:,:,:,2,2)=(U+V0)*(dummy(1,1)+dummy(2,2))/(KxGridNum)/(KyGridNum)/(KzGridNum)
     !enddo

     do o1=1,(MatDim)
          do o2=1,(MatDim)
               !if( U(o1,o2)/=0.0 ) then ! to find same position orbital. ( because onsite interaction almoste finite )
                    Hartree_k_iWn(:,1,1,:,o1,o1)=Hartree_k_iWn(:,1,1,:,o1,o1)+(U(o1,o2)+V0(o1,o2))*( dummy(o2,o2)  )/(K1GridNum)
                    !Hartree_k_iWn(:,1,1,:,(MatDim/2)+o1,(MatDim/2)+o2)=Hartree_k_iWn(:,1,1,:,(MatDim/2)+o1,(MatDim/2)+o2)+(U((MatDim/2)+o1,(MatDim/2)+o2)+V0((MatDim/2)+o1,(MatDim/2)+o2))*( dummy(o2,o2) + dummy( (MatDim/2)+o2, (MatDim/2)+o2 ) )/(K1GridNum)
               !endif
          enddo
     enddo

endif



!if (KzGridNum==1 .and. KyGridNum/=1) then        ! for 2d case
!     Hartree_k_iWn(:,:,:,:,1,1)=(U+4*V)*(dummy(1,1)+dummy(2,2))/(KxGridNum)/(KyGridNum)/(KzGridNum)
!     Hartree_k_iWn(:,:,:,:,2,2)=(U+4*V)*(dummy(1,1)+dummy(2,2))/(KxGridNum)/(KyGridNum)/(KzGridNum)
!endif


!if (KzGridNum==1 .and. KyGridNum==1) then        ! for 1d case
!     Hartree_k_iWn(:,:,:,:,1,1)=(U+2*V)*(dummy(1,1)+dummy(2,2))/(KxGridNum)/(KyGridNum)/(KzGridNum)
!     Hartree_k_iWn(:,:,:,:,2,2)=(U+2*V)*(dummy(1,1)+dummy(2,2))/(KxGridNum)/(KyGridNum)/(KzGridNum)
!endif



return
end subroutine



