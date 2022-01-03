
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for searching Chemical potential from bi-sect root method with given G and self energy as fixing Ntot through G=inv( inv(G) -selfenergy + ChemP) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine GWCTQMC_ChemPgivenN(Dimen,t1,U0,U1,ChemP,iWnlist,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,iWnGridNum,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,Noccup,ChemOccup,MatDim,delta,beta,Ntot,ChemTH,NonIntGreenf_k_iWn,GWCTQMCSelfE_k_iWn,dummy_k_iWn)
Implicit none
!argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum, MatDim
real(8) ,intent(in) ::  delta, beta,Noccup(MatDim),Ntot,ChemTH,t1(MatDim,MatDim),U0(MatDim,MatDim),U1(MatDim,MatDim)
real(8) ,intent(in) ::  K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3),Taulist(TauGridNum), R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3)
complex(8), intent(in) :: iWnlist(iWnGridNum),NonIntGreenf_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim),GWCTQMCSelfE_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)
real(8) ,intent(inout) :: ChemP
complex(8), intent(inout) :: ChemOccup(MatDim)
complex(8), intent(inout) ::  dummy_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)

!variable
integer i,j
real(8) lowoccupsum, upoccupsum, midoccupsum
real(8) lowene, upene, midene
complex(8) testOccup(MatDim)
print *, " "
print *, " "
print *, "-----------------------------------------------------------------"
print *, "*** Determining Chemical Potential from Calculated G, SelfE  ***"
print *," "

!ChemP=-2*abs(maxval(t1))-2*abs(maxval(U0))-2*abs(maxval(U1))
ChemP=-2*Dimen*abs(maxval(t1))-2*Dimen*abs(maxval(U0))-2*Dimen*abs(maxval(U1))

call  NewG_GWCTQMC(Dimen,ChemP,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,NonIntGreenf_k_iWn,GWCTQMCSelfE_k_iWn,dummy_k_iWn)
call G_k_iWn_Occup2(Dimen,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,Noccup,ChemOccup,MatDim,delta,ChemP,beta,dummy_k_iWn)



print 11,"[",1,"]   ", " Guess of Chemical Potential : ",ChemP
11 format(1x,A,I3,1x,A,1x,A,f)
print 12," Occupancy : ", ChemOccup
12 format(1x,A,30(f10.6))
print 13," total : ", real(sum(ChemOccup))
13 format(1x,A,f10.6)



lowoccupsum=real(sum(ChemOccup))
lowene=ChemP


if ( abs(lowoccupsum-Ntot)<ChemTH ) then
     print *," "
     print *," ~ Chempcial potential founded -- ChemP :", ChemP
     print 1," Given Occupancy (for Hamiltonian) :", sum(Noccup)
     1 format(1x,A,(f10.6))
     print 2," Calculated Occupancy :", real(ChemOccup)
     2 format(1x,A,30(f10.6))

     return
endif



!ChemP=2*abs(maxval(t1))+2*abs(maxval(U0))+2*abs(maxval(U1))
ChemP=2*Dimen*abs(maxval(t1))+2*Dimen*abs(maxval(U0))+2*Dimen*abs(maxval(U1))

call  NewG_GWCTQMC(Dimen,ChemP,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,NonIntGreenf_k_iWn,GWCTQMCSelfE_k_iWn,dummy_k_iWn)
call G_k_iWn_Occup2(Dimen,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,Noccup,ChemOccup,MatDim,delta,ChemP,beta,dummy_k_iWn)



print 14,"[",2,"]   ", " Guess of Chemical Potential : ",ChemP
14 format(1x,A,I3,1x,A,1x,A,f)
print 15," Occupancy : ", ChemOccup
15 format(1x,A,30(f10.6))
print 16," total : ", real(sum(ChemOccup))
16 format(1x,A,f10.6)

upoccupsum=real(sum(ChemOccup))
upene=ChemP

if ( abs(upoccupsum-Ntot)<ChemTH ) then
     print *," "
     print *," ~ Chempcial potential founded -- ChemP :", ChemP
     print 3," Given Occupancy (for Hamiltonian) :", sum(Noccup)
     3 format(1x,A,(f10.6))
     print 4," Calculated Occupancy :", real(ChemOccup)
     4 format(1x,A,30(f10.6))
     return
endif




ChemP=(upene+lowene)/(2.0)

call  NewG_GWCTQMC(Dimen,ChemP,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,NonIntGreenf_k_iWn,GWCTQMCSelfE_k_iWn,dummy_k_iWn)
call G_k_iWn_Occup2(Dimen,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,Noccup,ChemOccup,MatDim,delta,ChemP,beta,dummy_k_iWn)




print 17,"[",3,"]   ", " Guess of Chemical Potential : ",ChemP
17 format(1x,A,I3,1x,A,1x,A,f)
print 18," Occupancy : ", ChemOccup
18 format(1x,A,30(f10.6))
print 19," total : ", real(sum(ChemOccup))
19 format(1x,A,f10.6)



midoccupsum=real(sum(ChemOccup))
midene=ChemP

if ( abs(midoccupsum-Ntot)<ChemTH ) then
     print *," "
     print *," ~ Chempcial potential founded -- ChemP :", ChemP
     print 5," Given Occupancy (for Hamiltonian) :", sum(Noccup)
     5 format(1x,A,(f10.6))
     print 6," Calculated Occupancy :", real(ChemOccup)
     6 format(1x,A,30(f10.6))
     return
endif


j=3;
do i =1,100

     if( ( lowoccupsum<Ntot ) .and. ( Ntot<midoccupsum  ) ) then

          ChemP=(lowene+midene)/(2.0)



          call  NewG_GWCTQMC(Dimen,ChemP,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,NonIntGreenf_k_iWn,GWCTQMCSelfE_k_iWn,dummy_k_iWn)
          call G_k_iWn_Occup2(Dimen,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,Noccup,ChemOccup,MatDim,delta,ChemP,beta,dummy_k_iWn)

          j=j+1

          print 20,"[",j,"]   ", " Guess of Chemical Potential : ",ChemP
          20 format(1x,A,I3,1x,A,1x,A,f)
          print 21," Occupancy : ", ChemOccup
          21 format(1x,A,30(f10.6))
          print 22," total : ", real(sum(ChemOccup))
          22 format(1x,A,f10.6)


          upoccupsum=midoccupsum
          upene=midene
          midoccupsum=real(sum(ChemOccup))
          midene=ChemP

          if ( abs(midoccupsum-Ntot)<ChemTH ) then
               print *," "
               print *," ~ Chempcial potential founded -- ChemP :", ChemP
               print 7," Given Occupancy (for Hamiltonian) :", sum(Noccup)
               7 format(1x,A,(f10.6))
               print 8," Calculated Occupancy :", real(ChemOccup)
               8 format(1x,A,30(f10.6))
               return
          endif


     endif

     if( ( midoccupsum<Ntot ) .and. (Ntot<upoccupsum ) ) then



          ChemP=(midene+upene)/(2.0d0)


          call  NewG_GWCTQMC(Dimen,ChemP,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,NonIntGreenf_k_iWn,GWCTQMCSelfE_k_iWn,dummy_k_iWn)
          call G_k_iWn_Occup2(Dimen,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,Noccup,ChemOccup,MatDim,delta,ChemP,beta,dummy_k_iWn)

          j=j+1

          print 23,"[",j,"]   ", " Guess of Chemical Potential : ",ChemP
          23 format(1x,A,I3,1x,A,1x,A,f)
          print 24," Occupancy : ", ChemOccup
          24 format(1x,A,30(f10.6))
          print 25," total : ", real(sum(ChemOccup))
          25 format(1x,A,f10.6)


          lowoccupsum=midoccupsum
          lowene=midene
          midoccupsum=real(sum(ChemOccup))
          midene=ChemP
          if ( abs(midoccupsum-Ntot)<ChemTH ) then
               print *," "
               print *," ~ Chempcial potential founded -- ChemP :", ChemP
               print 9," Given Occupancy (for Hamiltonian) :", sum(Noccup)
               9 format(1x,A,(f10.6))
               print 10," Calculated Occupancy :", real(ChemOccup)
               10 format(1x,A,30(f10.6))
               return
          endif


     endif
end do

return
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for constructing new green function !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine NewG_GWCTQMC(Dimen,ChemP,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,NonIntGreenf_k_iWn,GWCTQMCSelfE_k_iWn,dummy_k_iWn)

Implicit none

! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim
real(8), intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3), Taulist(TauGridNum), R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3),beta,ChemP
complex(8), intent(in) :: iWnlist(iWnGridNum), NonIntGreenf_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim), GWCTQMCSelfE_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)
complex(8), intent(inout) :: dummy_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)

! local variable
integer i,j,kx,ky,kz,l,E,n
complex(8) Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim),dummy3(MatDim,MatDim)


if (Dimen==3) then
     do kx =1,(K1GridNum)
          do ky =1,(K2GridNum)
               do kz =1,(K3GridNum)
                    do n=1,iWnGridNum
                         call Identity(MatDim,Iden)
                         dummy(:,:)=(0.0d0,0.0d0)
                         dummy2(:,:)=(0.0d0,0.0d0)
                         dummy3(:,:)=(0.0d0,0.0d0)
                         !call MatInv(MatDim,NonIntGreenf_k_iWn(k,n,:,:),dummy)
                         call MatInv(MatDim,NonIntGreenf_k_iWn(kx,ky,kz,n,:,:),dummy)
                         dummy2(:,:)=dummy(:,:)-GWCTQMCSelfE_k_iWn(kx,ky,kz,n,:,:)+ChemP*Iden
                         !call MatInv(MatDim,dummy2,dummy3)      
                         call MatInv(MatDim,dummy2,dummy3)
                         dummy_k_iWn(kx,ky,kz,n,:,:)=dummy3
                    enddo
               enddo
          enddo
     enddo
endif



if (Dimen==2) then
     do kx =1,(K1GridNum)
          do ky =1,(K2GridNum)
               do n=1,iWnGridNum
                    call Identity(MatDim,Iden)
                    dummy(:,:)=(0.0d0,0.0d0)
                    dummy2(:,:)=(0.0d0,0.0d0)
                    dummy3(:,:)=(0.0d0,0.0d0)
                    !call MatInv(MatDim,NonIntGreenf_k_iWn(k,n,:,:),dummy)
                    call MatInv(MatDim,NonIntGreenf_k_iWn(kx,ky,1,n,:,:),dummy)
                    dummy2(:,:)=dummy(:,:)--GWCTQMCSelfE_k_iWn(kx,ky,1,n,:,:)+ChemP*Iden
                    !call MatInv(MatDim,dummy2,dummy3)      
                    call MatInv(MatDim,dummy2,dummy3)
                    dummy_k_iWn(kx,ky,1,n,:,:)=dummy3
               enddo
          enddo
     enddo
endif


if (Dimen==1) then
     do kx =1,(K1GridNum)
          do n=1,iWnGridNum
               call Identity(MatDim,Iden)
               dummy(:,:)=(0.0d0,0.0d0)
               dummy2(:,:)=(0.0d0,0.0d0)
               dummy3(:,:)=(0.0d0,0.0d0)
               !call MatInv(MatDim,NonIntGreenf_k_iWn(k,n,:,:),dummy)
               call MatInv(MatDim,NonIntGreenf_k_iWn(kx,1,1,n,:,:),dummy)
               dummy2(:,:)=dummy(:,:)-GWCTQMCSelfE_k_iWn(kx,1,1,n,:,:)+ChemP*Iden
               !call MatInv(MatDim,dummy2,dummy3)      
               call MatInv(MatDim,dummy2,dummy3)
               dummy_k_iWn(kx,1,1,n,:,:)=dummy3
          enddo
     enddo
endif







return

end subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating occupation for a given grenn function set of G(k,iWn) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine G_k_iWn_Occup2(Dimen,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,Noccup,ChemOccup,MatDim,delta,ChemP,beta,Greenf_k_iWn)
Implicit none
! arguments
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum, MatDim
real(8) ,intent(in) :: K1list(K1GridNum,3), K2list(K2GridNum,3),K3list(K3GridNum,3),delta, beta
real(8) ,intent(in) :: ChemP,Noccup(MatDim)
complex(8), intent(in) :: iWnlist(iWnGridNum)
complex(8), intent(in) ::  Greenf_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)
complex(8) ,intent(inout) :: ChemOccup(MatDim)

! local variable
integer i,jx,jy,jz,k,l,E,o
complex(8) Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim),dummy3(MatDim,MatDim)
complex(8) gtilda2(MatDim,MatDim), gtilda3(MatDim,MatDim), MatGcut(MatDim,MatDim), asymtoticG(MatDim,MatDim)



call Identity(MatDim,Iden)
gtilda2(:,:)=(0.0d0,0.0d0)
gtilda3(:,:)=(0.0d0,0.0d0) !initalize gtilda for calculating asymtotic part of matsubar green function
dummy(:,:)=(0.0d0,0.0d0) !initialize



if(Dimen==3) then
     do jx=1,(K1GridNum) ! H components constructure in kx-space, that is H_k
          do jy=1,(K2GridNum) ! H components constructure in ky-space, that is H_k
               do jz=1,(K3GridNum) ! H components constructure in kz-space, that is H_k
                    gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( Greenf_k_iWn(jx,jy,jz,iWnGridNum,:,:)+Greenf_k_iWn(jx,jy,jz,1,:,:) )/(2.0d0)
                    gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( Greenf_k_iWn(jx,jy,jz,iWnGridNum,:,:)-Greenf_k_iWn(jx,jy,jz,1,:,:) )/(2.0d0) - Iden/iWnlist(iWnGridNum) )
 
                    do k=(int(iWnGridNum/2)+1),(iWnGridNum) ! Matsubara frequency iWn sumation
                         dummy= dummy - ((1.0d0)/(beta))* ( Greenf_k_iWn(jx,jy,jz,k,:,:) - Iden/( iWnlist(k) ) - gtilda2/( iWnlist(k)*iWnlist(k) ) - gtilda3/( iWnlist(k)*iWnlist(k)*iWnlist(k) ) ) *exp(-iWnlist(k)*(beta-delta))
                         dummy= dummy - ((1.0d0)/(beta))* ( Greenf_k_iWn(jx,jy,jz,iWnGridNum+1-k,:,:) + Iden/( iWnlist(k) ) - gtilda2/( iWnlist(k)*iWnlist(k) ) + gtilda3/( iWnlist(k)*iWnlist(k)*iWnlist(k) ) ) *exp(iWnlist(k)*(beta-delta))
                    end do
 
                    asymtoticG =  ( Iden/(2.0d0) ) - ( gtilda2*beta*( (beta-delta)/((2.0d0)*(beta))-(0.25d0)) ) + ( gtilda3*beta*beta/(4.0d0)*( (beta-delta)*(beta-delta)/(beta*beta) - (beta-delta)/beta) )
                    dummy=dummy+asymtoticG
 
               end do
          end do
     end do

     do o=1,MatDim
          ChemOccup(o)=dummy(o,o)/dble(K1GridNum)/dble(K2GridNum)/dble(K3GridNum)
     enddo

end if




if(Dimen==2) then
     do jx=1,(K1GridNum) ! H components constructure in kx-space, that is H_k
          do jy=1,(K2GridNum) ! H components constructure in ky-space, that is H_k
               gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( Greenf_k_iWn(jx,jy,1,iWnGridNum,:,:)+Greenf_k_iWn(jx,jy,1,1,:,:) )/(2.0d0)
               gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( Greenf_k_iWn(jx,jy,1,iWnGridNum,:,:)-Greenf_k_iWn(jx,jy,1,1,:,:) )/(2.0d0) - Iden/iWnlist(iWnGridNum) )
          
               do k=(int(iWnGridNum/2)+1),(iWnGridNum) ! Matsubara frequency iWn sumation
                    dummy= dummy - ((1.0d0)/(beta))* ( Greenf_k_iWn(jx,jy,1,k,:,:) - Iden/( iWnlist(k) ) - gtilda2/( iWnlist(k)*iWnlist(k) ) - gtilda3/( iWnlist(k)*iWnlist(k)*iWnlist(k) ) ) *exp(-iWnlist(k)*(beta-delta))
                    dummy= dummy - ((1.0d0)/(beta))* ( Greenf_k_iWn(jx,jy,1,iWnGridNum+1-k,:,:) + Iden/( iWnlist(k) ) - gtilda2/( iWnlist(k)*iWnlist(k) ) + gtilda3/( iWnlist(k)*iWnlist(k)*iWnlist(k) ) ) *exp(iWnlist(k)*(beta-delta))
               end do
          
               asymtoticG =  ( Iden/(2.0d0) ) - ( gtilda2*beta*( (beta-delta)/((2.0d0)*(beta))-(0.25d0)) ) + ( gtilda3*beta*beta/(4.0d0)*( (beta-delta)*(beta-delta)/(beta*beta) - (beta-delta)/beta) )
               dummy=dummy+asymtoticG
     
          end do
     end do

     do o=1,MatDim
          ChemOccup(o)=dummy(o,o)/dble(K1GridNum)/dble(K2GridNum)
     enddo


end if



if(Dimen==1) then
     do jx=1,(K1GridNum) ! H components constructure in kx-space, that is H_k
          gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( Greenf_k_iWn(jx,1,1,iWnGridNum,:,:)+Greenf_k_iWn(jx,1,1,1,:,:) )/(2.0d0)
          gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( Greenf_k_iWn(jx,1,1,iWnGridNum,:,:)-Greenf_k_iWn(jx,1,1,1,:,:) )/(2.0d0) - Iden/iWnlist(iWnGridNum) )
     
          do k=(int(iWnGridNum/2)+1),(iWnGridNum) ! Matsubara frequency iWn sumation
               dummy= dummy - ((1.0d0)/(beta))* ( Greenf_k_iWn(jx,1,1,k,:,:) - Iden/( iWnlist(k) ) - gtilda2/( iWnlist(k)*iWnlist(k) ) - gtilda3/( iWnlist(k)*iWnlist(k)*iWnlist(k) ) ) *exp(-iWnlist(k)*(beta-delta))
               dummy= dummy - ((1.0d0)/(beta))* ( Greenf_k_iWn(jx,1,1,iWnGridNum+1-k,:,:) + Iden/( iWnlist(k) ) - gtilda2/( iWnlist(k)*iWnlist(k) ) + gtilda3/( iWnlist(k)*iWnlist(k)*iWnlist(k) ) ) *exp(iWnlist(k)*(beta-delta))
          end do
     
          asymtoticG =  ( Iden/(2.0d0) ) - ( gtilda2*beta*( (beta-delta)/((2.0d0)*(beta))-(0.25d0)) ) + ( gtilda3*beta*beta/(4.0d0)*( (beta-delta)*(beta-delta)/(beta*beta) - (beta-delta)/beta) )
          dummy=dummy+asymtoticG
     
     end do

     do o=1,MatDim
          ChemOccup(o)=dummy(o,o)/dble(K1GridNum)
     enddo

end if





return
end subroutine

































