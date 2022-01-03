
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating Fock self energy in (k,iWn) space !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Fock_self(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,U,Vq,MatDim,beta,delta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Greenf_k_iWn,Fock_k_iWn)

Implicit none
! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim
real(8), intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3), Taulist(TauGridNum), R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3),beta,U(MatDim,MatDim),delta
complex(8), intent(in) :: iWnlist(iWnGridNum), Greenf_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim),Vq(K1GridNum,K2GridNum,K3GridNum,MatDim,MatDim)
complex(8), intent(inout) :: Fock_k_iWn(K1GridNum,K2GridNum,K3GRidNum,iWnGridNum,MatDim,MatDim)

! local variable
integer i,j,kx,ky,kz,qx,qy,qz,l,E,n,n2,m,o1,o2
complex(8) Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim),dummy3(MatDim,MatDim)
complex(8) gtilda2(MatDim,MatDim), gtilda3(MatDim,MatDim), MatGcut(MatDim,MatDim), asymtoticG(MatDim,MatDim)


call Identity(MatDim,Iden)
dummy(:,:)=(0.0,0.0) !initialize

Fock_k_iWn(:,:,:,:,:,:)=(0.0,0.0) ! Initialize Fock_self Energy

if (Dimen==3) then
     do qx =1,(K1GridNum) ! kx space loop
          do qy =1,(K2GridNum) ! ky space loop
               do qz =1,(K3GridNum) ! kz space loop

                    dummy2(:,:)=(0.0,0.0) !initialize for every kpts
                    gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( Greenf_k_iWn(qx,qy,qz,iWnGridNum,:,:)+Greenf_k_iWn(qx,qy,qz,1,:,:) )/(2.0d0)
                    gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( Greenf_k_iWn(qx,qy,qz,iWnGridNum,:,:)-Greenf_k_iWn(qx,qy,qz,1,:,: ) )/(2.0d0) - Iden/iWnlist(iWnGridNum) )


                    do n=(int(iWnGridNum/2)+1),(iWnGridNum) !  iWn loop
 
                         dummy2= dummy2 - ((1.0d0)/(beta))* ( Greenf_k_iWn(qx,qy,qz,n,:,:) - Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) - gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(-iWnlist(n)*(beta-delta))
                         dummy2= dummy2 - ((1.0d0)/(beta))* ( Greenf_k_iWn(qx,qy,qz,iWnGridNum+1-n,:,:) + Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) + gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(iWnlist(n)*(beta-delta))


                    enddo
                    asymtoticG =  ( Iden/(2.0d0) ) - ( gtilda2*beta*( (beta-delta)/((2.0d0)*(beta))-(0.25d0)) ) + ( gtilda3*beta*beta/(4.0d0)*( (beta-delta)*(beta-delta)/(beta*beta) - (beta-delta)/beta) )
                    dummy2=dummy2+asymtoticG
                   
                    
                    do kx=1,(K1GridNum)
                         do ky=1,(K2GridNum)
                              do kz=1,(K3GridNum)
                                   do n2=1,(iWnGridNum)
                                       do o1=1,(MatDim)
                                            do o2=1,(MatDim)
                                                 if (kx>=qx .and. ky>=qy .and. kz>=qz) then
                                                      Fock_k_iWn(kx,ky,kz,n2,o1,o2)= Fock_k_iWn(kx,ky,kz,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(kx-qx+1,ky-qy+1,kz-qz+1,o1,o2))/(K1GridNum*K2GridNum*K3GridNum)
                                                 endif
                                                 if (kx<qx .and. ky>=qy .and. kz>=qz) then
                                                      Fock_k_iWn(kx,ky,kz,n2,o1,o2)= Fock_k_iWn(kx,ky,kz,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(k1GridNum+kx-qx+1,ky-qy+1,kz-qz+1,o1,o2))/(K1GridNum*K2GridNum*K3GridNum)
                                                 endif
                                                 if (kx>=qx .and. ky<qy .and. kz>=qz) then
                                                      Fock_k_iWn(kx,ky,kz,n2,o1,o2)= Fock_k_iWn(kx,ky,kz,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(kx-qx+1,K2GridNum+ky-qy+1,kz-qz+1,o1,o2))/(K1GridNum*K2GridNum*K3GridNum)
                                                 endif
                                                 if (kx>=qx .and. ky>=qy .and. kz<qz) then
                                                      Fock_k_iWn(kx,ky,kz,n2,o1,o2)= Fock_k_iWn(kx,ky,kz,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(kx-qx+1,ky-qy+1,K3GridNum+kz-qz+1,o1,o2))/(K1GridNum*K2GridNum*K3GridNum)
                                                 endif
                                                 if (kx>=qx .and. ky<qy .and. kz<qz) then
                                                      Fock_k_iWn(kx,ky,kz,n2,o1,o2)= Fock_k_iWn(kx,ky,kz,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(kx-qx+1,K2GridNum+ky-qy+1,K3GridNum+kz-qz+1,o1,o2))/(K1GridNum*K2GridNum*K3GridNum)
                                                 endif
                                                 if (kx<qx .and. ky>=qy .and. kz<qz) then
                                                      Fock_k_iWn(kx,ky,kz,n2,o1,o2)= Fock_k_iWn(kx,ky,kz,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(K1GridNum+kx-qx+1,ky-qy+1,K3GridNum+kz-qz+1,o1,o2))/(K1GridNum*K2GridNum*K3GridNum)
                                                 endif
                                                 if (kx<qx .and. ky<qy .and. kz>=qz) then
                                                      Fock_k_iWn(kx,ky,kz,n2,o1,o2)= Fock_k_iWn(kx,ky,kz,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(K1GridNum+kx-qx+1,K2GridNum+ky-qy+1,kz-qz+1,o1,o2))/(K1GridNum*K2GridNum*K3GridNum)
                                                 endif
                                                 if (kx<qx .and. ky<qy .and. kz<qz) then
                                                      Fock_k_iWn(kx,ky,kz,n2,o1,o2)= Fock_k_iWn(kx,ky,kz,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(K1GridNum+kx-qx+1,K2GridNum+ky-qy+1,K3GridNum+kz-qz+1,o1,o2))/(K1GridNum*K2GridNum*K3GridNum)
                                                 endif


                                             enddo
                                        enddo
                                   enddo
                              enddo
                         enddo
                    enddo


               enddo
          enddo
     enddo
endif





if (Dimen==2) then
     do qx =1,(K1GridNum) ! kx space loop
          do qy =1,(K2GridNum) ! ky space loop

              dummy2(:,:)=(0.0,0.0) !initialize for every kpts
              gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( Greenf_k_iWn(qx,qy,1,iWnGridNum,:,:)+Greenf_k_iWn(qx,qy,1,1,:,:) )/(2.0d0)
              gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( Greenf_k_iWn(qx,qy,1,iWnGridNum,:,:)-Greenf_k_iWn(qx,qy,1,1,:,: ) )/(2.0d0) - Iden/iWnlist(iWnGridNum) )


              do n=(int(iWnGridNum/2)+1),(iWnGridNum) !  iWn loop

                   dummy2= dummy2 - ((1.0d0)/(beta))* ( Greenf_k_iWn(qx,qy,1,n,:,:) - Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) - gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(-iWnlist(n)*(beta-delta))
                   dummy2= dummy2 - ((1.0d0)/(beta))* ( Greenf_k_iWn(qx,qy,1,iWnGridNum+1-n,:,:) + Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) + gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(iWnlist(n)*(beta-delta))


              enddo
              asymtoticG =  ( Iden/(2.0d0) ) - ( gtilda2*beta*( (beta-delta)/((2.0d0)*(beta))-(0.25d0)) ) + ( gtilda3*beta*beta/(4.0d0)*( (beta-delta)*(beta-delta)/(beta*beta) - (beta-delta)/beta) )
              dummy2=dummy2+asymtoticG
  
              !print *,"qx,qy:", qx,qy      
              !print *, "nom",dummy2

              do kx=1,(K1GridNum)
                   do ky=1,(K2GridNum)
                       do n2=1,(iWnGridNum)
                           do o1=1,(MatDim)
                                do o2=1,(MatDim)
                                     if (kx>=qx .and. ky>=qy) then
                                          Fock_k_iWn(kx,ky,1,n2,o1,o2)= Fock_k_iWn(kx,ky,1,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(kx-qx+1,ky-qy+1,1,o1,o2))/(K1GridNum*K2GridNum)
                                     endif
                                     if (kx<qx .and. ky>=qy ) then
                                          Fock_k_iWn(kx,ky,1,n2,o1,o2)= Fock_k_iWn(kx,ky,1,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(k1GridNum+kx-qx+1,ky-qy+1,1,o1,o2))/(K1GridNum*K2GridNum)
                                     endif
                                     if (kx>=qx .and. ky<qy ) then
                                          Fock_k_iWn(kx,ky,1,n2,o1,o2)= Fock_k_iWn(kx,ky,1,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(kx-qx+1,K2GridNum+ky-qy+1,1,o1,o2))/(K1GridNum*K2GridNum)
                                     endif
                                     if (kx<qx .and. ky<qy ) then
                                          Fock_k_iWn(kx,ky,1,n2,o1,o2)= Fock_k_iWn(kx,ky,1,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(K1GridNum+kx-qx+1,K2GridNum+ky-qy+1,1,o1,o2))/(K1GridNum*K2GridNum)
                                     endif


                                 enddo
                            enddo
                       enddo
                   enddo
              enddo


          enddo
     enddo
endif










if (Dimen==1) then
     do qx =1,(K1GridNum) ! kx space loop

           dummy2(:,:)=(0.0,0.0) !initialize for every kpts
           gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( Greenf_k_iWn(qx,1,1,iWnGridNum,:,:)+Greenf_k_iWn(qx,1,1,1,:,:) )/(2.0d0)
           gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( Greenf_k_iWn(qx,1,1,iWnGridNum,:,:)-Greenf_k_iWn(qx,1,1,1,:,: ) )/(2.0d0) - Iden/iWnlist(iWnGridNum) )
         
         
           do n=(int(iWnGridNum/2)+1),(iWnGridNum) !  iWn loop
         
                dummy2= dummy2 - ((1.0d0)/(beta))* ( Greenf_k_iWn(qx,1,1,n,:,:) - Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) - gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(-iWnlist(n)*(beta-delta))
                dummy2= dummy2 - ((1.0d0)/(beta))* ( Greenf_k_iWn(qx,1,1,iWnGridNum+1-n,:,:) + Iden/( iWnlist(n) ) - gtilda2/( iWnlist(n)*iWnlist(n) ) + gtilda3/( iWnlist(n)*iWnlist(n)*iWnlist(n) ) ) *exp(iWnlist(n)*(beta-delta))
         
         
           enddo
           asymtoticG =  ( Iden/(2.0d0) ) - ( gtilda2*beta*( (beta-delta)/((2.0d0)*(beta))-(0.25d0)) ) + ( gtilda3*beta*beta/(4.0d0)*( (beta-delta)*(beta-delta)/(beta*beta) - (beta-delta)/beta) )
           dummy2=dummy2+asymtoticG
         
         
         
           do kx=1,(K1GridNum)
              do n2=1,(iWnGridNum)
                  do o1=1,(MatDim)
                       do o2=1,(MatDim)
                            if (kx>=qx ) then
                                 Fock_k_iWn(kx,1,1,n2,o1,o2)= Fock_k_iWn(kx,1,1,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(kx-qx+1,1,1,o1,o2))/(K1GridNum)
                            endif
                            if (kx<qx  ) then
                                 Fock_k_iWn(kx,1,1,n2,o1,o2)= Fock_k_iWn(kx,1,1,n2,o1,o2)+dummy2(o1,o2)*(-U(o1,o2)-Vq(k1GridNum+kx-qx+1,1,1,o1,o2))/(K1GridNum)
                            endif
         
         
                        enddo
                   enddo
              enddo
           enddo


     enddo
endif
      



return
end subroutine


