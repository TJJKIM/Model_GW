
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating local slef-energy  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalGWDMFTSelfE(Dimen,GWSelfE_k_iWn,locSelfE,ImpuritySelfE,K1GridNum,K2GridNum,K3GridNum,MatDim,iWnGridNum,GWCTQMCSelfE_k_iWn)
Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,Dimen
complex(8), intent(in) :: GWSelfE_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim),locSelfE(iWnGridNum,MatDim,MatDim),ImpuritySelfE(iWnGridNum,MatDim)
complex(8), intent(inout) :: GWCTQMCSelfE_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)
! local variable
integer kx,ky,kz,n,m,Z
complex(8) GW_m_loc


if (Dimen==3) then
     do kx=1,(K1GridNum) !k space loop
          do ky=1,(K2GridNum) !k space loop
               do kz=1,(K3GridNum) !k space loop
                    do n=1,(iWnGridNum)


                         GWCTQMCSelfE_k_iWn(kx,ky,kz,n,:,:)= ( GWSelfE_k_iWn(kx,ky,kz,n,:,:) - locSelfE(n,:,:) )

                         do m=1,(MatDim)
                              GW_m_loc=GWCTQMCSelfE_k_iWn(kx,ky,kz,n,m,m)
                              GWCTQMCSelfE_k_iWn(kx,ky,kz,n,m,m)= GW_m_loc + ImpuritySelfE(n,m)
                         enddo


                    enddo
               enddo
          enddo
     enddo

endif


if (Dimen==2) then
     do kx=1,(K1GridNum) !k space loop
          do ky=1,(K2GridNum) !k space loop
               do n=1,(iWnGridNum)
            
            
                    GWCTQMCSelfE_k_iWn(kx,ky,1,n,:,:)= ( GWSelfE_k_iWn(kx,ky,1,n,:,:) - locSelfE(n,:,:) )
            
                    do m=1,(MatDim)
                         GW_m_loc=GWCTQMCSelfE_k_iWn(kx,ky,1,n,m,m)
                         GWCTQMCSelfE_k_iWn(kx,ky,1,n,m,m)= GW_m_loc + ImpuritySelfE(n,m)
                    enddo
            
            
               enddo
          enddo
     enddo

endif

if (Dimen==1) then
     do kx=1,(K1GridNum) !k space loop
           do n=1,(iWnGridNum)


                GWCTQMCSelfE_k_iWn(kx,1,1,n,:,:)= ( GWSelfE_k_iWn(kx,1,1,n,:,:) - locSelfE(n,:,:) )

                do m=1,(MatDim)
                     GW_m_loc=GWCTQMCSelfE_k_iWn(kx,1,1,n,m,m)
                     GWCTQMCSelfE_k_iWn(kx,1,1,n,m,m)= GW_m_loc + ImpuritySelfE(n,m)
                enddo


           enddo
     enddo

endif






return
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating local slef-energy  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalGWDMFTPol(Dimen,GWPol_k_iWn,locPol,ImpurityPol,K1GridNum,K2GridNum,K3GridNum,MatDim,iWnGridNum,GWCTQMCPol_k_iWn)
Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,Dimen
complex(8), intent(in) :: GWPol_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim),locPol(iWnGridNum,MatDim,MatDim),ImpurityPol(iWnGridNum,MatDim)
complex(8), intent(inout) :: GWCTQMCPol_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)
! local variable
integer kx,ky,kz,n,m,Z
complex(8) GW_m_loc


if (Dimen==3) then
     do kx=1,(K1GridNum) !k space loop
          do ky=1,(K2GridNum) !k space loop
               do kz=1,(K3GridNum) !k space loop
                    do n=1,(iWnGridNum)


                         GWCTQMCPol_k_iWn(kx,ky,kz,n,:,:)= ( GWPol_k_iWn(kx,ky,kz,n,:,:) - locPol(n,:,:) )

                         do m=1,(MatDim)
                              GW_m_loc=GWCTQMCPol_k_iWn(kx,ky,kz,n,m,m)
                              GWCTQMCPol_k_iWn(kx,ky,kz,n,m,m)= GW_m_loc + ImpurityPol(n,m)
                         enddo


                    enddo
               enddo
          enddo
     enddo

endif



if (Dimen==2) then
     do kx=1,(K1GridNum) !k space loop
          do ky=1,(K2GridNum) !k space loop
              do n=1,(iWnGridNum)


                   GWCTQMCPol_k_iWn(kx,ky,1,n,:,:)= ( GWPol_k_iWn(kx,ky,1,n,:,:) - locPol(n,:,:) )

                   do m=1,(MatDim)
                        GW_m_loc=GWCTQMCPol_k_iWn(kx,ky,1,n,m,m)
                        GWCTQMCPol_k_iWn(kx,ky,1,n,m,m)= GW_m_loc + ImpurityPol(n,m)
                   enddo


              enddo
          enddo
     enddo

endif




if (Dimen==3) then
     do kx=1,(K1GridNum) !k space loop
          do n=1,(iWnGridNum)


               GWCTQMCPol_k_iWn(kx,1,1,n,:,:)= ( GWPol_k_iWn(kx,1,1,n,:,:) - locPol(n,:,:) )

               do m=1,(MatDim)
                    GW_m_loc=GWCTQMCPol_k_iWn(kx,1,1,n,m,m)
                    GWCTQMCPol_k_iWn(kx,1,1,n,m,m)= GW_m_loc + ImpurityPol(n,m)
               enddo


          enddo
     enddo

endif






return
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating screend coulomb interactin Wc in GWDMFT !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GWDMFTWCal(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,U,Vq,MatDim,beta,K1list,K2list,K3list,iWnlist2,Taulist,R1list,R2list,R3list,Polarization_k_iWn,ScreenedW_k_iWn)


Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim, Dimen
real(8), intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3), Taulist(TauGridNum), R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3),beta,U(MatDim,MatDim)
complex(8), intent(in) :: iWnlist2(iWnGridNum+1),Vq(K1GridNum,K2GridNum,K3GridNum,MatDim,MatDim), Polarization_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim)
complex(8), intent(inout) :: ScreenedW_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim)

! local variable
integer i,j,kx,ky,kz,l,E,n
complex(8) Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim),dummy3(MatDim,MatDim)

if (Dimen==3) then
     call Identity(MatDim,Iden)
     do kx =1,(K1GridNum) ! k space loop
          do ky =1,(K2GridNum) ! k space loop
               do kz =1,(K3GridNum) ! k space loop
                    do n=1,(iWnGridNum+1) !  iWn loop
 
                         dummy=Iden-matmul( (U(:,:)+Vq(kx,ky,kz,:,:)), Polarization_k_iWn(kx,ky,kz,n,:,:) )
                         call MatInv(MatDim,dummy,dummy2)
                         dummy3=matmul( (U(:,:)+Vq(kx,ky,kz,:,:)), dummy2 )
                         ScreenedW_k_iWn(kx,ky,kz,n,:,:)=dummy3
                    enddo
               enddo
          enddo
     enddo

endif



if (Dimen==2) then
     call Identity(MatDim,Iden)
     do kx =1,(K1GridNum) ! k space loop
          do ky =1,(K2GridNum) ! k space loop
               do n=1,(iWnGridNum+1) !  iWn loop

                    dummy=Iden-matmul( (U(:,:)+Vq(kx,ky,kz,:,:)), Polarization_k_iWn(kx,ky,kz,n,:,:) )
                    call MatInv(MatDim,dummy,dummy2)
                    dummy3=matmul( (U(:,:)+Vq(kx,ky,kz,:,:)), dummy2 )
                    ScreenedW_k_iWn(kx,ky,kz,n,:,:)=dummy3
               enddo
          enddo
     enddo

endif



if (Dimen==1) then
     call Identity(MatDim,Iden)
     do kx =1,(K1GridNum) ! k space loop
          do n=1,(iWnGridNum+1) !  iWn loop

               dummy=Iden-matmul( (U(:,:)+Vq(kx,ky,kz,:,:)), Polarization_k_iWn(kx,ky,kz,n,:,:) )
               call MatInv(MatDim,dummy,dummy2)
               dummy3=matmul( (U(:,:)+Vq(kx,ky,kz,:,:)), dummy2 )
               ScreenedW_k_iWn(kx,ky,kz,n,:,:)=dummy3
          enddo
     enddo

endif


return
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for reading GWDMFT loop Num !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine WriteGWDMFTloopNum(GWDMFTloopNum)
Implicit none
!arguments
integer, intent(in) :: GWDMFTloopNum

open (unit=2, file ="GWDMFT.loopNum")
write(2,*) ,int(GWDMFTloopNum)
close(2)

return
end subroutine


