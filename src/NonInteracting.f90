!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for searching Chemical potential & converged occupancy Noccup !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine NonIntChemP(MinV,MaxV,Dimen,OrbN,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,tq,t0,U,Noccup,ChemOccup,MatDim,delta,ChemP,beta,Ntot,ChemTH,lwork,DoCTQMC,NonIntHeigVal,NonIntGreenf_k_iWn,NonIntGreenf0ChemP_k_iWn)
Implicit none
!argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum,iWnGridNum, MatDim, lwork, DoCTQMC , Dimen,OrbN
real(8) ,intent(in) :: MinV,MaxV,K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3),  delta, beta,Ntot,ChemTH ,U(MatDim,MatDim) ,t0(MatDim,MatDim)
complex(8), intent(in) :: iWnlist(iWnGridNum), tq(K1GridNum,K2GridNum,K3GridNum,MatDim,MatDim)
real(8) ,intent(inout) :: ChemP
real (8), intent(inout) :: Noccup(MatDim)
complex(8), intent(inout) :: ChemOccup(MatDim)
complex(8), intent(inout) ::  NonIntHeigVal(K1GridNum,K2GridNum,K3GridNum,MatDim,1),NonIntGreenf_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim),NonIntGreenf0ChemP_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)

!variable
integer i,j
real(8) lowoccupsum, upoccupsum, midoccupsum
real(8) lowene, upene, midene



do i=1,100
     call NonIntChemPgivenN(MinV,MaxV,Dimen,OrbN,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,tq,t0,U,Noccup,ChemOccup,MatDim,delta,ChemP,beta,Ntot,ChemTH,lwork,DoCTQMC,NonIntHeigVal,NonIntGreenf_k_iWn,NonIntGreenf0ChemP_k_iWn)


     if ( abs( sum(Noccup(:)-real(ChemOccup(:)) ) )<ChemTH*MatDim ) then
         Noccup(:)=real(ChemOccup(:))
         print *," "
         print *, "------------------------------------------------------------------------"
         print *, "Input occupancy for constract hamiltonian and output Occupancy Converged"
         print 1, " ChemP : ",ChemP
         1 format(1x,A,f)
         print 2, "Noccup : ",Noccup
         2 format(1x,A,30(f10.6))
         print *, "------------------------------------------------------------------------"
         return
     endif

     Noccup(:)=real(ChemOccup(:))
enddo





return
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for searching Chemical potential from bi-sect root method with hamiltonian constructred by given occupancy (Noccup) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine NonIntChemPgivenN(MinV,MaxV,Dimen,OrbN,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,tq,t0,U,Noccup,ChemOccup,MatDim,delta,ChemP,beta,Ntot,ChemTH,lwork,DoCTQMC,NonIntHeigVal,NonIntGreenf_k_iWn,NonIntGreenf0ChemP_k_iWn)
Implicit none
!argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum,iWnGridNum, MatDim, lwork,DoCTQMC,Dimen,OrbN
real(8) ,intent(in) :: MinV,MaxV,K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3), delta, beta,Noccup(MatDim),Ntot,ChemTH, t0(MatDim,MatDim),U(MatDim,MatDim)
complex(8), intent(in) :: iWnlist(iWnGridNum), tq(K1GridNum,K2GridNum,K3GridNum,MatDim,MatDim)
real(8) ,intent(inout) :: ChemP
complex(8), intent(inout) :: ChemOccup(MatDim)
complex(8), intent(inout) ::  NonIntHeigVal(K1GridNum,K2GridNum,K3GridNum,MatDim,1), NonIntGreenf_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim),NonIntGreenf0ChemP_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)

!variable
integer i,j
real(8) lowoccupsum, upoccupsum, midoccupsum
real(8) lowene, upene, midene
complex(8) testOccup(2,1)
print *, " "
print *, " "
print *, "------------------------------------------------------------------------------"
print *, "*** Determining Chemical Potential for Non-interacting Hamiltonian with n:",int(sum(noccup))," ***"
print *," "

!print *,"**********HERE**********" 
ChemP=MinV



call NonIntOccupancy2(Dimen,OrbN,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,tq,t0,U,Noccup,ChemOccup,MatDim,delta,ChemP,beta,lwork,DoCTQMC,NonIntHeigVal,NonIntGreenf_k_iWn,NonIntGreenf0ChemP_k_iWn) ! Noccup is initial occupation(total occupation fixed with sum of initial occupation), ChemOccup is Occupation for a given ChemP within Hamiltonian constaructed by Noccup

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


ChemP=MaxV

call NonIntOccupancy2(Dimen,OrbN,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,tq,t0,U,Noccup,ChemOccup,MatDim,delta,ChemP,beta,lwork,DoCTQMC,NonIntHeigVal,NonIntGreenf_k_iWn,NonIntGreenf0ChemP_k_iWn) ! Noccup is initial occupation(total occupation fixed with sum of initial occupation), ChemOccup is Occupation for a given ChemP within Hamiltonian constaructed by Noccup

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




ChemP=(upene+lowene)/(2.0d0)

call NonIntOccupancy2(Dimen,OrbN,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,tq,t0,U,Noccup,ChemOccup,MatDim,delta,ChemP,beta,lwork,DoCTQMC,NonIntHeigVal,NonIntGreenf_k_iWn,NonIntGreenf0ChemP_k_iWn) ! Noccup is initial occupation(total occupation fixed with sum of initial occupation), ChemOccup is Occupation for a given ChemP within Hamiltonian constaructed by Noccup

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

          ChemP=(lowene+midene)/(2.0d0)

          call NonIntOccupancy2(Dimen,OrbN,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,tq,t0,U,Noccup,ChemOccup,MatDim,delta,ChemP,beta,lwork,DoCTQMC,NonIntHeigVal,NonIntGreenf_k_iWn,NonIntGreenf0ChemP_k_iWn) ! Noccup is initial occupation(total occupation fixed with sum of initial occupation), ChemOccup is Occupation for a given ChemP within Hamiltonian constaructed by Noccup
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

          call NonIntOccupancy2(Dimen,OrbN,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,tq,t0,U,Noccup,ChemOccup,MatDim,delta,ChemP,beta,lwork,DoCTQMC,NonIntHeigVal,NonIntGreenf_k_iWn,NonIntGreenf0ChemP_k_iWn) ! Noccup is initial occupation(total occupation fixed with sum of initial occupation), ChemOccup is Occupation for a given ChemP within Hamiltonian constaructed by Noccup
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




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for making a model hamiltonian !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine NonIntOccupancy2(Dimen,OrbN,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,tq,t0,U,Noccup,ChemOccup,MatDim,delta,ChemP,beta,lwork,DoCTQMC,NonIntHeigVal,NonIntGreenf_k_iWn,NonIntGreenf0ChemP_k_iWn)
Implicit none
! arguments
integer, intent(in) :: K1GridNum, K2GridNum,K3GridNum, iWnGridNum, MatDim, DoCTQMC, lwork, Dimen, OrbN
real(8) ,intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3),  delta, beta
real(8) ,intent(in) :: ChemP,Noccup(MatDim), U(MatDim,MatDim), t0(MatDim,MatDim)
complex(8), intent(in) :: iWnlist(iWnGridNum), tq(K1GridNum,K2GridNum,K3GridNum,MatDim,MatDim)
complex(8) ,intent(inout) :: ChemOccup(MatDim)
complex(8), intent(inout) :: NonIntHeigVal(K1GridNum,K2GridNum,K3GridNum,MatDim,1), NonIntGreenf_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim),NonIntGreenf0ChemP_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)

! local variable
integer i,jx,jy,jz,k,l,E,p,n
complex(8) Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG_noChem1(MatDim,MatDim),dummyG_noChem2(MatDim,MatDim),dummyG1(MatDim,MatDim),dummyG2(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim),dummy3(MatDim,MatDim)
complex(8) gtilda2(MatDim,MatDim), gtilda3(MatDim,MatDim), MatGcut1(MatDim,MatDim), MatGcut2(MatDim,MatDim) ,asymtoticG(MatDim,MatDim)



          gtilda2(:,:)=(0.0d0,0.0d0)
          gtilda3(:,:)=(0.0d0,0.0d0) !initalize gtilda for calculating asymtotic part of matsubar green function
          dummy(:,:)=(0.0d0,0.0d0) !initialize


          if (Dimen==3) then
               do jx=1,(K1GridNum) ! H components constructure in kx-space, that is H_k
                    do jy=1,(K2GridNum) ! H components constructure in ky-space, that is H_k
                         do jz=1,(K3GridNum) ! H components constructure in kz-space, that is H_k
                         call Identity(MatDim,Iden)
                         call SingleTightHam(MatDim,OrbN,tq(jx,jy,jz,:,:),t0,U,dummyH)
 
                              if ( DoCTQMC == 1 .or. DoCTQMC==2) then
                                     call Eigen_valvec(dummyH,NonIntHeigVal(jx,jy,jz,:,1),dummy2,MatDim,lwork)
                              endif


                              call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(1),MatGcut1) ! Matsubaragreen function at cutoff matsubar frequency ,w0, and klist(j)  
                              call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(iWnGridNum),MatGcut2) ! Matsubaragreen function at cutoff matsubar frequency ,w0, and klist(j)  
                             !call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(1),MatGcut1) ! Matsubaragreen function at cutoff matsubar frequency ,w0, and klist(j)  
                             !call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(iWnGridNum),MatGcut2) ! Matsubaragreen function at cutoff matsubar frequency ,w0, and klist(j)  
                              gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( MatGcut2+MatGcut1 )/(2.0d0)
                              gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( MatGCut2-MatGCut1 )/(2.0d0) - Iden/iWnlist(iWnGridNum) )

                              do k=(int(iWnGridNum/2)+1),(iWnGridNum) ! Matsubara frequency iWn sumation
                                   dummyG1(:,:)=(0.0d0,0.0d0)
                                   dummyG2(:,:)=(0.0d0,0.0d0)
                                   call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(k),dummyG1)
                                   call MatsubaraG(MatDim,dummyH,0.0d0,iWnlist(k),dummyG_noChem1)
                                   call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(iWnGridNum+1-k),dummyG2)
                                   call MatsubaraG(MatDim,dummyH,0.0d0,iWnlist(iWnGridNum+1-k),dummyG_noChem2)
                                  !call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(k),dummyG1)
                                  !call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(iWnGridNum+1-k),dummyG2)
                                   NonIntGreenf_k_iWn(jx,jy,jz,k,:,:)=dummyG1

                                   !TEST
                                   !print *, ""
                                   !print 1,dummyG1(1,1),dummyG1(1,2),dummyG1(1,3), dummyG1(1,4), dummyG1(1,5), dummyG1(1,6)
                                   !1 format (6(f,f))
                                   !print 2,dummyG1(2,1),dummyG1(2,2),dummyG1(2,3), dummyG1(2,4), dummyG1(2,5), dummyG1(2,6)
                                   !2 format (6(f,f))
                                   !print 3,dummyG1(3,1),dummyG1(3,2),dummyG1(3,3), dummyG1(3,4), dummyG1(3,5), dummyG1(3,6)
                                   !3 format (6(f,f))
                                   !print 4,dummyG1(4,1),dummyG1(4,2),dummyG1(4,3), dummyG1(4,4), dummyG1(4,5), dummyG1(4,6)
                                   !4 format (6(f,f))
                                   !print 5,dummyG1(5,1),dummyG1(5,2),dummyG1(5,3), dummyG1(5,4), dummyG1(5,5), dummyG1(5,6)
                                   !5 format (6(f,f))
                                   !print 6,dummyG1(6,1),dummyG1(6,2),dummyG1(6,3), dummyG1(6,4), dummyG1(6,5), dummyG1(6,6)
                                   !6 format (6(f,f))
                                   !print *, ""


                                   NonIntGreenf0ChemP_k_iWn(jx,jy,jz,k,:,:)=dummyG_noChem1
                                   NonIntGreenf_k_iWn(jx,jy,jz,iWnGridNum+1-k,:,:)=dummyG2
                                   NonIntGreenf0ChemP_k_iWn(jx,jy,jz,iWnGridNum+1-k,:,:)=dummyG_noChem2
                                   !!!print *,"(",j,",",k,")",Greenf_k_iWn(j,k,:,:)
                                   dummy= dummy - ((1.0d0)/(beta))* ( dummyG1 - Iden/( iWnlist(k) ) - gtilda2/( iWnlist(k)*iWnlist(k) ) - gtilda3/( iWnlist(k)*iWnlist(k)*iWnlist(k) ) ) *exp(-iWnlist(k)*(beta-delta))
                                   dummy= dummy - ((1.0d0)/(beta))* ( dummyG2 + Iden/( iWnlist(k) ) - gtilda2/( iWnlist(k)*iWnlist(k) ) + gtilda3/( iWnlist(k)*iWnlist(k)*iWnlist(k) ) ) *exp(iWnlist(k)*(beta-delta))
                              end do

                              asymtoticG =  ( Iden/(2.0d0) ) - ( gtilda2*beta*( (beta-delta)/((2.0d0)*(beta))-(0.25d0)) ) + ( gtilda3*beta*beta/(4.0d0)*( (beta-delta)*(beta-delta)/(beta*beta) - (beta-delta)/beta) )
                              dummy=dummy+asymtoticG

                          end do
                     end do
               end do
               do n=1,MatDim
                    ChemOccup(n)=dummy(n,n)/dble(K1GridNum)/dble(K2GridNum)/dble(K3GridNum)
               enddo
          endif



          if(Dimen==2) then

               do jx=1,(K1GridNum) ! H components constructure in kx-space, that is H_k
                    do jy=1,(K2GridNum) ! H components constructure in ky-space, that is H_k
                         call Identity(MatDim,Iden)
                         call SingleTightHam(MatDim,OrbN,tq(jx,jy,1,:,:),t0,U,dummyH)

                         if ( DoCTQMC == 1 .or. DoCTQMC==2) then
                              call Eigen_valvec(dummyH,NonIntHeigVal(jx,jy,1,:,1),dummy2,MatDim,lwork)
                         endif


                         call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(1),MatGcut1) ! Matsubaragreen function at cutoff matsubar frequency ,w0, and klist(j)  
                         call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(iWnGridNum),MatGcut2) ! Matsubaragreen function at cutoff matsubar frequency ,w0, and klist(j)  
                    !call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(1),MatGcut1) ! Matsubaragreen function at cutoff matsubar frequency ,w0, and klist(j)  
                    !call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(iWnGridNum),MatGcut2) ! Matsubaragreen function at cutoff matsubar frequency ,w0, and klist(j)  
                         gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( MatGcut2+MatGcut1 )/(2.0d0)
                         gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( MatGCut2-MatGCut1 )/(2.0d0) - Iden/iWnlist(iWnGridNum) )

                         do k=(int(iWnGridNum/2)+1),(iWnGridNum) ! Matsubara frequency iWn sumation
                              dummyG1(:,:)=(0.0d0,0.0d0)
                              dummyG2(:,:)=(0.0d0,0.0d0)
                              call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(k),dummyG1)
                              call MatsubaraG(MatDim,dummyH,0.0d0,iWnlist(k),dummyG_noChem1)
                              call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(iWnGridNum+1-k),dummyG2)
                              call MatsubaraG(MatDim,dummyH,0.0d0,iWnlist(iWnGridNum+1-k),dummyG_noChem2)
                              !call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(k),dummyG1)
                              !call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(iWnGridNum+1-k),dummyG2)
                              NonIntGreenf_k_iWn(jx,jy,1,k,:,:)=dummyG1
                              NonIntGreenf0ChemP_k_iWn(jx,jy,1,k,:,:)=dummyG_noChem1
                              NonIntGreenf_k_iWn(jx,jy,1,iWnGridNum+1-k,:,:)=dummyG2
                              NonIntGreenf0ChemP_k_iWn(jx,jy,1,iWnGridNum+1-k,:,:)=dummyG_noChem2
                              !!!print *,"(",j,",",k,")",Greenf_k_iWn(j,k,:,:)
                              dummy= dummy - ((1.0d0)/(beta))* ( dummyG1 - Iden/( iWnlist(k) ) - gtilda2/( iWnlist(k)*iWnlist(k) ) - gtilda3/( iWnlist(k)*iWnlist(k)*iWnlist(k) ) ) *exp(-iWnlist(k)*(beta-delta))
                              dummy= dummy - ((1.0d0)/(beta))* ( dummyG2 + Iden/( iWnlist(k) ) - gtilda2/( iWnlist(k)*iWnlist(k) ) + gtilda3/( iWnlist(k)*iWnlist(k)*iWnlist(k) ) ) *exp(iWnlist(k)*(beta-delta))
                         end do

                         asymtoticG =  ( Iden/(2.0d0) ) - ( gtilda2*beta*( (beta-delta)/((2.0d0)*(beta))-(0.25d0)) ) + ( gtilda3*beta*beta/(4.0d0)*( (beta-delta)*(beta-delta)/(beta*beta) - (beta-delta)/beta) )
                         dummy=dummy+asymtoticG

                    end do
               end do
               do n=1,MatDim
                    ChemOccup(n)=dummy(n,n)/dble(K1GridNum)/dble(K2GridNum)
               enddo

          endif
         
          if (Dimen==1) then

               do jx=1,(K1GridNum) ! H components constructure in k-space, that is H_k

                    call Identity(MatDim,Iden)
                    call SingleTightHam(MatDim,OrbN,tq(jx,1,1,:,:),t0,U,dummyH)

                    if (DoCTQMC == 1 .or. DoCTQMC==2) then
                         call Eigen_valvec(dummyH,NonIntHeigVal(jx,1,1,:,1),dummy2,MatDim,lwork)
                    endif


                    call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(1),MatGcut1) ! Matsubaragreen function at cutoff matsubar frequency ,w0, and klist(j)  
                    call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(iWnGridNum),MatGcut2) ! Matsubaragreen function at cutoff matsubar frequency ,w0, and klist(j)  
                    !call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(1),MatGcut1) ! Matsubaragreen function at cutoff matsubar frequency ,w0, and klist(j)  
                    !call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(iWnGridNum),MatGcut2) ! Matsubaragreen function at cutoff matsubar frequency ,w0, and klist(j)  
                    gtilda2(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( MatGcut2+MatGcut1 )/(2.0d0)
                    gtilda3(:,:)=iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*iWnlist(iWnGridNum)*( ( MatGCut2-MatGCut1 )/(2.0d0) - Iden/iWnlist(iWnGridNum) )

                    do k=(int(iWnGridNum/2)+1),(iWnGridNum) ! Matsubara frequency iWn sumation
                         dummyG1(:,:)=(0.0d0,0.0d0)
                         dummyG2(:,:)=(0.0d0,0.0d0)
                         call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(k),dummyG1)
                         call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(iWnGridNum+1-k),dummyG2)
                         call MatsubaraG(MatDim,dummyH,0.0,iWnlist(k),dummyG_noChem1)
                         call MatsubaraG(MatDim,dummyH,0.0,iWnlist(iWnGridNum+1-k),dummyG_noChem2)
                         !call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(k),dummyG1)
                         !call MatsubaraG(MatDim,dummyH,ChemP,iWnlist(iWnGridNum+1-k),dummyG2)
                         NonIntGreenf_k_iWn(jx,1,1,k,:,:)=dummyG1
                         NonIntGreenf0ChemP_k_iWn(jx,1,1,k,:,:)=dummyG_noChem1
                         NonIntGreenf_k_iWn(jx,1,1,iWnGridNum+1-k,:,:)=dummyG2
                         NonIntGreenf0ChemP_k_iWn(jx,1,1,iWnGridNum+1-k,:,:)=dummyG_noChem2
                         !!!print *,"(",j,",",k,")",Greenf_k_iWn(j,k,:,:)
                         dummy= dummy - ((1.0d0)/(beta))* ( dummyG1 - Iden/( iWnlist(k) ) - gtilda2/( iWnlist(k)*iWnlist(k) ) - gtilda3/( iWnlist(k)*iWnlist(k)*iWnlist(k) ) ) *exp(-iWnlist(k)*(beta-delta))
                         dummy= dummy - ((1.0d0)/(beta))* ( dummyG2 + Iden/( iWnlist(k) ) - gtilda2/( iWnlist(k)*iWnlist(k) ) + gtilda3/( iWnlist(k)*iWnlist(k)*iWnlist(k) ) ) *exp(iWnlist(k)*(beta-delta))
                    end do

                    asymtoticG =  ( Iden/(2.0d0) ) - ( gtilda2*beta*( (beta-delta)/((2.0d0)*(beta))-(0.25d0)) ) + ( gtilda3*beta*beta/(4.0d0)*( (beta-delta)*(beta-delta)/(beta*beta) - (beta-delta)/beta) )
                    dummy=dummy+asymtoticG

               end do
               do n=1,MatDim
                    ChemOccup(n)=dummy(n,n)/dble(K1GridNum)
               enddo

          endif



return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for construction of identity matrix for a given dimension !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine Identity(MatDim,IdenM)

Implicit none
!arguments
integer, intent(in) :: MatDim
complex(8), intent(inout) :: IdenM(MatDim,MatDim)

!local variable
integer i
IdenM(:,:)=(0.0d0,0.0d0)
do i=1,MatDim
   IdenM(i,i)=(1.0d0,0.0d0)
end do

return
end subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for making Hartree Hamiltonian !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SingleTightHam(MatDim,OrbN,tk,t0,U,HartreeH)
Implicit none

!arguments
integer, intent(in) :: MatDim,OrbN
real(8), intent(in) :: t0(MatDim,MatDim),U(MatDim,MatDim)
complex(8), intent(in) :: tk(MatDim,MatDim)
complex(8), intent(inout) :: HartreeH(MatDim,MatDim)


!local variable
integer i,j,k
complex(8) Iden(MatDim,MatDim)
 
HartreeH=tk
do i=1,MatDim
     HartreeH(i,i)=-HartreeH(i,i)+t0(i,i)-U(i,i)/(2.0d0)

enddo
!HartreeH(1:OrbN,1:OrbN)=HartreeH(1:OrbN,1:OrbN)+t0(1:OrbN,1:OrbN)-U(1:OrbN,1:OrbN)/(2.0d0)
!HartreeH(OrbN+1:2*OrbN,OrbN+1:2*OrbN)=HartreeH(OrbN+1:2*OrbN,OrbN+1:2*OrbN)+t0(1:OrbN,1:OrbN)-U(1:OrbN,1:OrbN)/(2.0d0)

!print *,""
!print 1, HartreeH(1,1), HartreeH(1,2), HartreeH(1,3), HartreeH(1,4)
!1 format(4(f8.4,f8.4))
!print 2, HartreeH(2,1), HartreeH(2,2), HartreeH(2,3), HartreeH(2,4)
!2 format(4(f8.4,f8.4))
!print 3, HartreeH(3,1), HartreeH(3,2), HartreeH(3,3), HartreeH(3,4)
!3 format(4(f8.4,f8.4))
!print 4, HartreeH(4,1), HartreeH(4,2), HartreeH(4,3), HartreeH(4,4)
!4 format(4(f8.4,f8.4))
!print *,""


return
end subroutine






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for making matsubara green function !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MatsubaraG(MatDim,Ham,ChemP,Matsubf,MatG)
Implicit none

!arguments
integer, intent(in) :: MatDim
real(8), intent(in) :: ChemP
complex(8), intent(in) :: Matsubf, Ham(MatDim,MatDim)
complex(8), intent(inout) :: MatG(MatDim,MatDim)

!local variable 
complex(8) IdenM(MatDim,MatDim),dumG(MatDim,MatDim)

call Identity(MatDim,IdenM)
call MatInv(MatDim,((Matsubf+ChemP)*IdenM-Ham),dumG)
MatG=dumG


return
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for obtaining inverse matrix of given complex matrix !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MatInv(MatDim,A,InvA)   ! for complex matrix inverse 
Implicit none
! arguments
integer, intent(in) :: MatDim
complex(8), intent(in) :: A(MatDim,MatDim)
complex(8), intent(inout) :: InvA(MatDim,MatDim)

!local variable 
integer ipiv(MatDim), info
complex(8)  work(MatDim*MatDim)  ! work array for LAPACK


  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  InvA = A

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call zgetrf(MatDim, MatDim, InvA, MatDim, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call zgetri(MatDim, InvA, MatDim, ipiv, work, MatDim, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if

return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for obtaining eigenvalue and -vector of given complex matrix hamiltonian !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Eigen_valvec(H,eigval,eigvec,MatDim,lwork)
Implicit none
!arguments
integer, intent(in) :: MatDim
integer, intent(in) :: lwork
complex(8), intent(in) :: H(MatDim,MatDim)
complex(8), intent(inout) :: eigval(MatDim)
complex(8), intent(inout) :: eigvec(MatDim,MatDim)

!local variable
double precision rwork(2*MatDim)
complex(8) work(lwork)
integer info
complex(8) dum(MatDim,MatDim)

!lwork=3*MatDim-1

call zgeev('N','V',MatDim,H,MatDim,eigval,dum,MatDim,eigvec,MatDim,work,lwork,rwork,info)

return
end subroutine








!!!!!!!!!!!!!!!1 ETC



subroutine FindMaxValInMat(tq,D,val)
Implicit none
!arguments
integer, intent(in) :: D
real(8), intent(in) :: tq(D,D)
real(8), intent(inout) :: val

!local variable
integer i,j
real(8) maxv

maxv=0.0d0
!lwork=3*MatDim-1
do i=1,D
     do j=1,D
          if(tq(i,j)>maxv) then
               maxv=tq(i,j)
          endif
     enddo
enddo
val=maxv
return
end subroutine


subroutine FindMinValInMat(tq,D,val)
Implicit none
!arguments
integer, intent(in) :: D
real(8), intent(in) :: tq(D,D)
real(8), intent(inout) :: val

!local variable
integer i,j
real(8) minv

minv=0.0d0
!lwork=3*MatDim-1
do i=1,D
     do j=1,D
          if(tq(i,j)<minv) then
               minv=tq(i,j)
          endif
     enddo
enddo
val=minv
return
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for initializing occupation number !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initializeOccupancy(Ntot,MatDim,InitOccup,Noccup)
Implicit none

!arguments
integer, intent(in) :: InitOccup, MatDim
real(8), intent(in) :: Ntot
real(8), intent(inout) :: Noccup(MatDim)

!loacal variable
integer i,j

!para
if ( InitOccup == 1 ) then
     do i=1,MatDim

           Noccup(i)=Ntot/(MatDim)

     end do
endif

!ferro
if ( InitOccup == 2 ) then
     do i=1,MatDim

           if (i < (MatDim/2+1) ) then
               Noccup(i)=Ntot/(MatDim)*(1.9d0)
           else
               Noccup(i)=Ntot/(MatDim)*(0.1d0)
           endif

     end do
endif

!antiferro
!if ( InitOccup == 3 ) then
!       do i=1,2
!         do j=1,MatDim
!           if( mod(i+j,2) == 0) then
!               Noccup(i,j)=Ntot/(MatDim*2)*(1.9)
!           else 
!               Noccup(i,j)=Ntot/(MatDim*2)*(0.1)
!           endif
!
!         end do
!     end do
!
!endif


return
end subroutine


