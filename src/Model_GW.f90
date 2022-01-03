program HFpGW


use modmain
Implicit none

integer iii,jjj,G,Z,TE,NumDiffElement,C, GWDMFTloopNum,i,j,k,m,n 
real(8) ChemP,ChemPOld,DiffElement,diffChemP,NonIntChemPval, maxv,minv,dum
complex(8) testOccup(2,1)
complex(8) test(2,2), test2(2,2),test3(2,2)

call readinput
allocate(Noccup(MatDim))
allocate(ChemOccup(MatDim))
allocate(localE(MatDim))
allocate(NeighborInfo(OrbN,3,MaxOrbitN,4)) ! NeighborInfo(i,j,k,l) : j th neighbor of i th orbital, l=1,2,3 for cartesian coord, l=4 for orbital number
allocate(NeighborDist(OrbN,3))
allocate(NeighborN(OrbN,3))
allocate(tq(K1GridNum,K2GridNum,K3GridNum,MatDim,MatDim))
allocate(Vq(K1GridNum,K2GridNum,K3GridNum,MatDim,MatDim))
allocate(R1list(R1GridNum,3))
allocate(R2list(R2GridNum,3))
allocate(R3list(R3GridNum,3))
allocate(K1list(K1GridNum,3))
allocate(K2list(K1GridNum,3))
allocate(K3list(K1GridNum,3))
allocate(Taulist(TauGridNum))
allocate(Taulist2(TauGridNum+1))
allocate(iWnlist(iWnGridNum))
allocate(iWnlist2(iWnGridNum+1))
!original allocateval

allocate(Greenf_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim))
allocate(NonIntGreenf_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim))
allocate(NonIntGreenf0ChemP_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim))
allocate(Polarization_R_tau(R1GridNum,R2GridNum,R3GridNum,TauGridNum,MatDim,MatDim))
allocate(SelfEcorr_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim))
allocate(Hartree_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim))
allocate(Fock_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim))
allocate(ScreenedWcorr_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim))
allocate(W_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim))
allocate(dummy_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim))
allocate(dummy_k_iWn2(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim))
allocate(dummy_k_tau(K1GridNum,K2GridNum,K3GridNum,TauGridNum,MatDim,MatDim))
allocate(dummy_k_tauII(K1GridNum,K2GridNum,K3GridNum,TauGridNum,MatDim,MatDim))
allocate(dummy_R_iWn(R1GridNum,R2GridNum,R3GridNum,iWnGridNum,MatDim,MatDim))
allocate(dummy_R_tau(R1GridNum,R2GridNum,R3GridNum,TauGridNum,MatDim,MatDim))
allocate(dummy_R_tauInv(R1GridNum,R2GridNum,R3GridNum,TauGridNum,MatDim,MatDim))
allocate(dummy_R_tau2(R1GridNum,R2GridNum,R3GridNum,TauGridNum+1,MatDim,MatDim))
allocate(dummy_R_tauII(R1GridNum,R2GridNum,R3GridNum,TauGridNum,MatDim,MatDim))
allocate(dummy_R_tauIII(R1GridNum,R2GridNum,R3GridNum,TauGridNum,MatDim,MatDim))
allocate(locGreenf(iWnGridNum,MatDim,MatDim))
allocate(locW(iWnGridNum+1,MatDim,MatDim))
allocate(locPol(iWnGridNum+1,MatDim,MatDim))
allocate(locSelfE(iWnGridNum,MatDim,MatDim))
allocate(NonIntHeigVal(K1GridNum,K2GridNum,K3GridNum,MatDim,1))
allocate(ImpuritySelfE(iWnGridNum,MatDim))
allocate(ImpurityPol(iWnGridNum+1,MatDim))
allocate(GWSelfE_k_iWn_Old(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim))
allocate(GWPol_k_iWn_Old(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim))
allocate(GWCTQMCSelfE_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim))
allocate(GWCTQMCPol_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim))
allocate(Hartree_loc(iWnGridNum,MatDim,MatDim))
allocate(Fock_loc(iWnGridNum,MatDim,MatDim))







print *,""
print *, "--- Diemnsion ---"
print *,Dimen

print *,""
print *, "--- Orbital # ---"
print *, OrbN

print *,""
print *, "--- MatDim ---"
print *, MatDim

print *,""
print *, "--- MaxOrbit ---"
print *, MaxOrbitN


print *,""
print *, "--- Unit Vector ---"
     print 3,"a1 = ", a1(:)
     3 format(1x,A,f8.4,f8.4,f8.4)
     print 4,"a2 = ", a2(:)
     4 format(1x,A,f8.4,f8.4,f8.4)
     print 5,"a3 = ", a3(:)
     5 format(1x,A,f8.4,f8.4,f8.4)

call ReciprocalV(a1,a2,a3,b1,b2,b3,PI)
print *,""
print *, "--- Reciprocal Vector ---"
     print 6,"b1 = ", b1(:)
     6 format(1x,A,f8.4,f8.4,f8.4)
     print 7,"b2 = ", b2(:)
     7 format(1x,A,f8.4,f8.4,f8.4)
     print 8,"b3 = ", b3(:)
     8 format(1x,A,f8.4,f8.4,f8.4)



print *,""
print *, "--- Orbital position (Frac)---"
do i=1,OrbN
     print 9,"[",int(i),"]:", OrbitPosi(i,:)
     9 format(1x,A,I3,1x,A,f8.4,f8.4,f8.4)
enddo

print *,""
print *, "--- Orbital position (Cart)---"
do i=1,OrbN
     print 10,"[",int(i),"]:", OrbitPosi_car(i,:)
     10 format(1x,A,I3,1x,A,f8.4,f8.4,f8.4)
enddo


call FindSameAtomOrb(OrbitPosi_car,OrbN,distingPos,distingAtom)
print *,"--- Distinguished position ---"
do i=1,OrbN
     print 60,"[",int(i),"]:", distingPos(i,:)
     60 format(1x,A,I3,1x,A,f8.4,f8.4,f8.4)
enddo

print *,"--- Distinguished Atom ---"
do i=1,OrbN
     print 70,"[",int(i),"]:", distingAtom(i,:)
     70 format(1x,A,I3,1x,A,30(I3))
enddo


print *,""
print *, "--- Spin Symm Breaking (chemical shift val) ---"
print 80, SpinSplit(:) 
      80 format(30(f8.4))

print *,""
print *, "--- Grid ---"
print *, "Tau: ", TauGridNum
print *, "iWn: ", iWnGridNum
print *, "K1: ", K1GridNum
print *, "K2: ", K2GridNum
print *, "K3: ", K3GridNum
print *, "R1: ", R1GridNum
print *, "R2: ", R2GridNum
print *, "R3: ", R3GridNum

print *,""
print *, "--- MaxloopN ---"
print *, SCFloop

print *,""
print *, "--- OnsiteE ---"
do i=1,MatDim
     print 11, t0(i,1:MatDim)
     11 format(30(f8.4))
enddo


print *,""
print *, "--- Hopping ---"
print *, "[NN]"
do i=1,MatDim
     print 12, t1(i,1:MatDim)
     12 format(30(f8.4))
enddo
print *, "[2NN]"
do i=1,MatDim
     print 13, t2(i,1:MatDim)
     13 format(30(f8.4))
enddo
print *, "[3NN]"
do i=1,MatDim
     print 14, t3(i,1:MatDim)
     14 format(30(f8.4))
enddo

print *, ""
print *, "--- Coulomb Int ---"
print *, "[Onsite]"
do i=1,MatDim
     print 15, U(i,1:MatDim)
     15 format(30(f8.4))
enddo
print *, "[NN]"
do i=1,MatDim
     print 16, U1(i,1:MatDim)
     16 format(30(f8.4))
enddo
print *, "[2NN]"
do i=1,MatDim
     print 17, U2(i,1:MatDim)
     17 format(30(f8.4))
enddo
print *, "[3NN]"
do i=1,MAtDim
     print 18, U3(i,1:MatDim)
     18 format(30(f8.4))
enddo


print *,""
print *, "--- Temperature ---"
print 19, "T (in K) : ",Tem
19 format(1x,A,f8.4)
print 20, "beta : ",beta
20 format(1x,A,f8.4)

print *,""
print *, "--- Occupation ---"
print 21, Ntot
      21 format((f8.4))

print *,""
print *, "--- Criterion ---"
print *, "Mu threshold : ", ChemTH
print *, "Greenf diff threshold : ", ConvTH


print *,""
print *, "--- Output file ---"
print 25, "G(k,iWn) out : ", allgreenfSave
25 format(1x,A,I3)
print 26, "LocSelfE out : ", locselfESave
26 format(1x,A,I3)
print 27, "Every step log out : ", LogOutOn
27 format(1x,A,I3)


print *,""
print *, "--- Mixing parameters (New portion) ---"
print 28, "GW Polarization : ", GWMixing_Pol
28 format(1x,A,f8.4)
print 22, "GW Self Energy : ", GWMixing_SelfE
22 format(1x,A,f8.4)
print 23, "Imp Polarization : ", DMFTMixing_Pol
23 format(1x,A,f8.4)
print 24, "Imp Self Energy : ", DMFTMixing_SelfE
24 format(1x,A,f8.4)


print *,""
print *, "--- ++ method ---"
print *, " * CTQMC =1 : DMFT, CTQMC=2 : EDMFT, GWCTQMC =1 : GW+EDMFT       * "
print *, " * HartreeFock =1 : Self-consitent HF, all=0 : Self-consitent GW * "
print 90, "Hartree-Fock : ", DoHF
90 format(1x,A,I3)
print 29, "CTQMC : ", DoCTQMC
29 format(1x,A,I3)
print 31, "GWCTQMC : ", DoGW_CTQMC
31 format(1x,A,I3)

print *,""
print *, "--- DMFT restart with Mu shift ---"
print 32, " Mu shifted at 1st step of DMFT :", MuShift
32 format(1x,A,f8.4)


call LatticeConst(Dimen,a1,a2,a3, OrbN, MaxOrbitN, OrbitPosi_car, NeighborInfo, NeighborDist,NeighborN)

print *,""
print *,""
print *, "--- Neighbor Info ---"
do i=1,OrbN 
    print 30,"[Orbital :",i,"] :", NeighborN(i,1),", ",NeighborN(i,2),", ", NeighborN(i,3)
    30 format(1x,A,I3,1x,A,I3,1x,A,I3,1x,A,I3)
    print *,""
    do j=1,3
        print *,""
        print 1,j,"-th neighbor : Dist =",NeighborDist(i,j)
        1 format(I3,1x,A,f8.4)
        do n=1,NeighborN(i,j)
             print 2,"[",n,"]- ","position: ", NeighborInfo(i,j,n,1:3),",","atom :", int(NeighborInfo(i,j,n,4))
             2 format(1x,A,I3,1x,A,1x,A,f8.4,f8.4,f8.4,1x,A1,3x,A10,I3)
        enddo
    enddo
print *, ""
print *,""
enddo
print *, ""


call initializeOccupancy(Ntot,MatDim,InitOccup,Noccup)
Call Kgrid(b1,K1GridNum,K1list)
Call Kgrid(b2,K2GridNum,K2list)
Call Kgrid(b3,K3GridNum,K3list)



Call Rgrid(a1,R1GridNum,R1list)
Call Rgrid(a2,R2GridNum,R2list)
Call Rgrid(a3,R3GridNum,R3list)




call gridfromZ4(beta,TauGridNum,Taulist)
call gridfromZEven(beta,TauGridNum,Taulist2)
call MatsubfGrid(PI,iWnGridNum,iWnlist,beta)
call MatsubfGridEven(PI,iWnGridNum,iWnlist2,beta)



call parameterR2K(Dimen,OrbN,MatDim,MaxOrbitN,K1GridNum,K2GridNum,K3GridNum,K1list,K2list,K3list,R1GridNum,R2GridNum,R3GridNum,R1list,R2list,R3list,t0,t1,t2,t3,U,U1,U2,U3,OrbitPosi_car,NeighborN,NeighborInfo,tq,Vq)
call Save_tqVq(Dimen,K1GridNum,K2GridNum,K3GridNum,K1list,K2list,K3list,MatDim,tq,Vq)





call FindMaxValInMat(t0,OrbN,dum)
maxv=maxv+dum
call FindMaxValInMat(t1,OrbN,dum)
maxv=maxv+6*dum
call FindMaxValInMat(t2,OrbN,dum)
maxv=maxv+6*dum
call FindMaxValInMat(t3,OrbN,dum)
maxv=maxv+4*dum
call FindMaxValInMat(U,OrbN,dum)
maxv=maxv+dum
call FindMaxValInMat(U1,OrbN,dum)
maxv=maxv+3*dum
call FindMaxValInMat(U2,OrbN,dum)
maxv=maxv+6*dum
call FindMaxValInMat(U3,OrbN,dum)
maxv=maxv+4*dum


call FindMinValInMat(t0,OrbN,dum)
minv=minv+dum
call FindMaxValInMat(t1,OrbN,dum)
minv=minv-6*dum
call FindMaxValInMat(t2,OrbN,dum)
minv=minv-6*dum
call FindMaxValInMat(t3,OrbN,dum)
minv=minv-4*dum
call FindMaxValInMat(U,OrbN,dum)
minv=minv-dum
call FindMaxValInMat(U1,OrbN,dum)
minv=minv-3*dum
call FindMaxValInMat(U2,OrbN,dum)
minv=minv-6*dum
call FindMaxValInMat(U3,OrbN,dum)
minv=minv-4*dum


call FindMaxValInMat(t0,OrbN,dum)

!!! If Chemical potential of Noninteracting part can't be founded change these value which determin the chemical potential scan range
maxv=maxv-(MatDim*dum)
minv=minv+(MatDim*dum)

maxv=10.0
minv=-10.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! Main part !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!-----------------------------------------!
!---------- Noniteracting part -----------!
!-----------------------------------------!


print *," "
print *," "
print *," *********************************** "
print *," ******* Noninteracting part ******* "
print *," *********************************** "
print *," "
print *," "

call NonIntChemP(minv,maxv,Dimen,OrbN,K1list,K2list,K3list,iWnlist,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,tq,t0,U,Noccup,ChemOccup,MatDim,delta,ChemP,beta,Ntot,ChemTH,lwork,DoCTQMC,NonIntHeigVal,NonIntGreenf_k_iWn,NonIntGreenf0ChemP_k_iWn)

NonIntChemPval=ChemP
Greenf_k_iWn=NonIntGreenf_k_iWn


if (DoCTQMC==1 .or. DoCTQMC==2 ) then
 call SaveNonInteigV(Dimen,K1GridNum,K2GridNum,K3GridNum,K1list,K2list,K3list,MatDim,ChemP,localE,NonIntHeigVal)
endif


if (SCFloop==0 .or. DoCTQMC ==1 .or. DoCTQMC==2 ) then
 call CallocG(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,Greenf_k_iWn,locGreenf)
 call SavelocG(OrbN,DoCTQMC,localE,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,ChemP,NonIntChemPval,Noccup,beta,MatDim,iWnlist,locGreenf)
 stop
endif



!---------------------------------------------------!
!---------- Self-consitent GW+EDMFT loop -----------!
!---------------------------------------------------!


if (DoGW_CTQMC==1) then

     call ReadGWDMFTloopNum(GWDMFTloopNum)
     print *," "
     print *," "
     print *," ****************************"
     print *," ******* GW+DMFT PART *******"
     print *," ****************************"
     print *," "
     print *," "
     print *," "


     if (GWDMFTloopNum==0) then
          print *," *---------------------------------------------------------*"
          print *," *-- Saving Initilize SelfE = 0.0 and Polarization = 0.0 --*"
          print *," *---------------------------------------------------------*"

          call PrintSelfE_k_iWn_Init(Dimen,OrbN,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist,MatDim)
          call PrintPol_k_iWn_Init(Dimen,OrbN,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist2,MatDim)

          print *," Reading GW SelfE_k_iWn, Pol_k_iWn ... "
          call ReadSelfE_k_iWn(Dimen,OrbN,MatDim,iWnGridNum,K1GridNum,K2GridNum,K3GridNum,GWSelfE_k_iWn_Old)
          call ReadPol_k_iWn(Dimen,OrbN,MatDim,iWnGridNum,K1GridNum,K2GridNum,K3GridNum,GWPol_k_iWn_Old)
 
          !print *, GWPol_k_iWn_Old(1,1,1,1,1,1)
          !print *, GWPol_k_iWn_Old(1,1,1,(iWnGridNum/2)-1,1,1)
          !print *, GWPol_k_iWn_Old(1,1,1,(iWnGridNum/2),1,1)
          !print *, GWPol_k_iWn_Old(1,1,1,(iWnGridNum/2)+1,1,1)
          !print *, GWPol_k_iWn_Old(1,1,1,(iWnGridNum/2)+1,1,2)
          !print *, GWPol_k_iWn_Old(1,1,1,(iWnGridNum/2)+2,1,1)

          print *," Calculating local selfE, local Pol ... "
          call CalocSelfE(Dimen,MatDim,GWSelfE_k_iWn_Old,locSelfE,K1GridNum,K2GridNum,K3GridNum,iWnGridNum)
          call CalocPol(Dimen,MatDim,GWPol_k_iWn_Old,locPol,K1GridNum,K2GridNum,K3GridNum,iWnGridNum)

          print *," Reading Impurity self-energy, polarization ... "
          call ReadImpuritySelfE(MatDim,iWnGridNum,ImpuritySelfE)
          call ReadImpurityPol(MatDim,iWnGridNum,ImpurityPol)


          print *," Calculating GW+DMFT SelfE_k_iWn, Pol_k_iWn ... "
          call CalGWDMFTSelfE(Dimen,GWSelfE_k_iWn_Old,locSelfE,ImpuritySelfE,K1GridNum,K2GridNum,K3GridNum,MatDim,iWnGridNum,GWCTQMCSelfE_k_iWn)
          call CalGWDMFTPol(Dimen,GWPol_k_iWn_Old,locPol,ImpurityPol,K1GridNum,K2GridNum,K3GridNum,MatDim,iWnGridNum,GWCTQMCPol_k_iWn)


          print *," Print GW+DMFT SelfE_k_iWn, Pol_k_iWn ... "
          call  PrintGWEDMFTSelfE(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist,MatDim,GWCTQMCSelfE_k_iWn)
          call  PrintGWEDMFTPol(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist2,MatDim,GWCTQMCPol_k_iWn)


          print *,"*----- Cal ChemP and GreenF with GW+DMFT selfE for hyb...-----"
          call GWCTQMC_ChemPgivenN(Dimen,t1,U,U1,ChemP,iWnlist,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,iWnGridNum,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,Noccup,ChemOccup,MatDim,delta,beta,Ntot,ChemTH,NonIntGreenf0ChemP_k_iWn,GWCTQMCSelfE_k_iWn,Greenf_k_iWn) ! In this step, ChemP are determined, and dummy_k_iWn will become Next Greenfunction calculated with determined ChemP -- Here Greenf_k_iWn changed as a New Greenfunction

          print *," Calculating W_k_iWn from Pol_k_iWn ... "
          call GWDMFTWCal(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,U,Vq,MatDim,beta,K1list,K2list,K3list,iWnlist2,Taulist,R1list,R2list,R3list,GWCTQMCPol_k_iWn,W_k_iWn) ! Here W_k_iWn is the screened W part



          !!!!!!!!!!!!!!!!!!! Wiess  field out for DMFT calculation!!!!!!!!!!!!!!!!!!!!!!!!!

          print *," Calculating local Greenfunction from G_k_iWn ..."
          call CallocG(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,Greenf_k_iWn,locGreenf)
          print *," Calculating local W from W_k_iWn ... "
          call CallocW(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,W_k_iWn,locW)
          print *," Calculating local SelfE again with New GWCTQMCSelfE_k_iWn ..."
          call CalocSelfE(Dimen,MatDim,GWCTQMCSelfE_k_iWn,locSelfE,K1GridNum,K2GridNum,K3GridNum,iWnGridNum)
          print *," Calculating local Polarization again with New GWCTQMCPol_k_iWn ..."
          call CalocPol(Dimen,MatDim,GWCTQMCPol_k_iWn,locPol,K1GridNum,K2GridNum,K3GridNum,iWnGridNum)




          print *," Printing hybridization, hybridization.dat...."
          call Printhybridization(locGreenf,locSelfE,MatDim,iWnGridNum,iWnlist,ChemP,beta,Noccup,U)
          print *," Printing Nonlocal Retarded Interaction, NonlocRetarded.dat ...."
          call PrintRetardedInt(locW,locPol,iWnGridNum,MatDim,U,iWnlist2,GWDMFTloopNum)


          !call PrintlocW_iWn_wPolMat(dummy_k_iWn2,beta,KxGridNum,KyGridNum,iWnGridNum,Kxlist,Kylist,iWnlist2,ax,ay,MatDim,U,V,1)
          call Printdyn_iWn(locW,locPol,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,iWnlist2,U,MatDim,1)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






          call Hartree_self(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,U,Vq(1,1,1,:,:),MatDim,beta,delta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Greenf_k_iWn,Hartree_k_iWn)

          print *,"Fourier transform of G : (k,iWn)->(k,tau) "
          call k_iWn2k_tau(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Greenf_k_iWn,dummy_k_tau)


          print *,"Fourier transform of G : (k,tau)->(R,tau) "
          call k_tau2R_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_R_tau)
          call R_tau2k_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_R_tauInv)
          if (Dimen==3) then
               dummy_R_tauInv=dummy_R_tauInv/dble(K1GridNum)/dble(K2GridNum)/dble(K3GridNum) ! because there is no normalization in inverse trainsfomation
          endif
          if (Dimen==2) then
               dummy_R_tauInv=dummy_R_tauInv/dble(K1GridNum)/dble(K2GridNum) ! because there is no normalization in inverse trainsfomation
          endif
          if (Dimen==1) then
               dummy_R_tauInv=dummy_R_tauInv/dble(K1GridNum) ! because there is no normalization in inverse trainsfomation
          endif

          print *,"Polarization, (Chi) calculation in (R,tau)  "
          call PolCal(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_R_tau,dummy_R_tauInv,Polarization_R_tau)


          print *,"Fourier transform of Chi : (R,tau)->(k,tau) "
          call R_tau2k_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Polarization_R_tau,dummy_k_tau)

          print *,"Fourier transform of Chi : (k,tau)->(k,iWn) "
          call k_tau2k_iWn_Even(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist2,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_k_iWn2)



          print *,"Fourier transform of W : (k,iWn)->(k,tau) "
          call k_iWn2k_tau_Even(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist2,Taulist,R1list,R2list,R3list,W_k_iWn,dummy_k_tau)

          print *,"Fourier transform of W : (k,tau)->(R,tau) "
          call k_tau2R_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_R_tauII)




          print *,"Calculating GW-self energy in (R,tau) space "
          call GWselfenergy(Dimen,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,dummy_R_tau,dummy_R_tauII,dummy_R_tauIII)

          print *,"Fourier transform of Self energy : (R,tau)->(k,tau) "
          call R_tau2k_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_R_tauIII,dummy_k_tau)


          print *,"Fourier transform of Self energy : (k,tau)->(k,iWn) "
          call k_tau2k_iWn(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_k_tau,SelfEcorr_k_iWn)

          ChemPOld = ChemP

          call PrintSelfE_k_iWn(Dimen,OrbN,Hartree_k_iWn,SelfEcorr_k_iWn,GWSelfE_k_iWn_Old,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3GridNum,iWnlist,MatDim,GWMixing_SelfE)
          call PrintPol_k_iWn(Dimen,OrbN,dummy_k_iWn2,GWPol_k_iWn_Old,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist2,MatDim,GWMixing_Pol)


          print *,"---------------------------------------------------------"
          print *,"-[GW] SelfE_k_iWn.dat, [GW] Pol_k_iWn.dat files are made-"
          print *,"---------------------------------------------------------"
          GWDMFTloopNum=GWDMFTloopNum+1
          call WriteGWDMFTloopNum(GWDMFTloopNum)
          stop



     endif



     print *," Reading GW SelfE_k_iWn, Pol_k_iWn ... "
     call ReadSelfE_k_iWn(Dimen,OrbN,MatDim,iWnGridNum,K1GridNum,K2GridNum,K3GridNum,GWSelfE_k_iWn_Old)
     call ReadPol_k_iWn(Dimen,OrbN,MatDim,iWnGridNum,K1GridNum,K2GridNum,K3GridNum,GWPol_k_iWn_Old)


     print *," Calculating local selfE, local Pol ... "
     call CalocSelfE(Dimen,MatDim,GWSelfE_k_iWn_Old,locSelfE,K1GridNum,K2GridNum,K3GridNum,iWnGridNum)
     call CalocPol(Dimen,MatDim,GWPol_k_iWn_Old,locPol,K1GridNum,K2GridNum,K3GridNum,iWnGridNum)


     print *," Reading Impurity self-energy, polarization ... "
     call ReadImpuritySelfE(MatDim,iWnGridNum,ImpuritySelfE)
     call ReadImpurityPol(MatDim,iWnGridNum,ImpurityPol)


     print *," Calculating GW+DMFT SelfE_k_iWn, Pol_k_iWn ... "
     call CalGWDMFTSelfE(Dimen,GWSelfE_k_iWn_Old,locSelfE,ImpuritySelfE,K1GridNum,K2GridNum,K3GridNum,MatDim,iWnGridNum,GWCTQMCSelfE_k_iWn)
     call CalGWDMFTPol(Dimen,GWPol_k_iWn_Old,locPol,ImpurityPol,K1GridNum,K2GridNum,K3GridNum,MatDim,iWnGridNum,GWCTQMCPol_k_iWn)


     print *," Print GW+DMFT SelfE_k_iWn, Pol_k_iWn ... "
     call  PrintGWEDMFTSelfE(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist,MatDim,GWCTQMCSelfE_k_iWn)
     call  PrintGWEDMFTPol(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist2,MatDim,GWCTQMCPol_k_iWn)



     print *,"*----- Cal ChemP and GreenF with GW+DMFT selfE for hyb...-----"
     call GWCTQMC_ChemPgivenN(Dimen,t1,U,U1,ChemP,iWnlist,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,iWnGridNum,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,Noccup,ChemOccup,MatDim,delta,beta,Ntot,ChemTH,NonIntGreenf0ChemP_k_iWn,GWCTQMCSelfE_k_iWn,Greenf_k_iWn) ! In this step, ChemP are determined, and dummy_k_iWn will become Next Greenfunction calculated with determined ChemP -- Here Greenf_k_iWn changed as a New Greenfunction

     print *," Calculating W_k_iWn from Pol_k_iWn ... "
     call GWDMFTWCal(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,U,Vq,MatDim,beta,K1list,K2list,K3list,iWnlist2,Taulist,R1list,R2list,R3list,GWCTQMCPol_k_iWn,W_k_iWn) ! Here W_k_iWn is the screened W part


     !!!!!!!!!!!!!!!!!!! Wiess  field out for DMFT calculation!!!!!!!!!!!!!!!!!!!!!!!!!
     print *," Calculating local Greenfunction from G_k_iWn ..."
     call CallocG(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,Greenf_k_iWn,locGreenf)
     print *," Calculating local W from W_k_iWn ... "
     call CallocW(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,W_k_iWn,locW)
     print *," Calculating local SelfE again with New GWCTQMCSelfE_k_iWn ..."
     call CalocSelfE(Dimen,MatDim,GWCTQMCSelfE_k_iWn,locSelfE,K1GridNum,K2GridNum,K3GridNum,iWnGridNum)
     print *," Calculating local Polarization again with New GWCTQMCPol_k_iWn ..."
     call CalocPol(Dimen,MatDim,GWCTQMCPol_k_iWn,locPol,K1GridNum,K2GridNum,K3GridNum,iWnGridNum)




     print *," Printing hybridization, hybridization.dat...."
     call Printhybridization(locGreenf,locSelfE,MatDim,iWnGridNum,iWnlist,ChemP,beta,Noccup,U)
     print *," Printing Nonlocal Retarded Interaction, NonlocRetarded.dat ...."
     call PrintRetardedInt(locW,locPol,iWnGridNum,MatDim,U,iWnlist2,GWDMFTloopNum)


     !call PrintlocW_iWn_wPolMat(dummy_k_iWn2,beta,KxGridNum,KyGridNum,iWnGridNum,Kxlist,Kylist,iWnlist2,ax,ay,MatDim,U,V,1)
     call Printdyn_iWn(locW,locPol,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,iWnlist2,U,MatDim,GWDMFTloopNum+1)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     call Hartree_self(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,U,Vq(1,1,1,:,:),MatDim,beta,delta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Greenf_k_iWn,Hartree_k_iWn)

     print *,"Fourier transform of G : (k,iWn)->(k,tau) "
     call k_iWn2k_tau(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Greenf_k_iWn,dummy_k_tau)

     print *,"Fourier transform of G : (k,tau)->(R,tau) "
     call k_tau2R_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_R_tau)
     call R_tau2k_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_R_tauInv)
     if (Dimen==3) then
          dummy_R_tauInv=dummy_R_tauInv/dble(K1GridNum)/dble(K2GridNum)/dble(K3GridNum) ! because there is no normalization in inverse trainsfomation
     endif
     if (Dimen==2) then
          dummy_R_tauInv=dummy_R_tauInv/dble(K1GridNum)/dble(K2GridNum) ! because there is no normalization in inverse trainsfomation
     endif
     if (Dimen==1) then
          dummy_R_tauInv=dummy_R_tauInv/dble(K1GridNum) ! because there is no normalization in inverse trainsfomation
     endif

     print *,"Polarization, (Chi) calculation in (R,tau)  "
     call PolCal(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_R_tau,dummy_R_tauInv,Polarization_R_tau)


     print *,"Fourier transform of Chi : (R,tau)->(k,tau) "
     call R_tau2k_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Polarization_R_tau,dummy_k_tau)

     print *,"Fourier transform of Chi : (k,tau)->(k,iWn) "
     call k_tau2k_iWn_Even(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist2,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_k_iWn2)



     print *,"Fourier transform of W : (k,iWn)->(k,tau) "
     call k_iWn2k_tau_Even(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist2,Taulist,R1list,R2list,R3list,W_k_iWn,dummy_k_tau)

     print *,"Fourier transform of W : (k,tau)->(R,tau) "
     call k_tau2R_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_R_tauII)




     print *,"Calculating GW-self energy in (R,tau) space "
     call GWselfenergy(Dimen,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,dummy_R_tau,dummy_R_tauII,dummy_R_tauIII)

     print *,"Fourier transform of Self energy : (R,tau)->(k,tau) "
     call R_tau2k_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_R_tauIII,dummy_k_tau)


     print *,"Fourier transform of Self energy : (k,tau)->(k,iWn) "
     call k_tau2k_iWn(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_k_tau,SelfEcorr_k_iWn)

     ChemPOld = ChemP

     call PrintSelfE_k_iWn(Dimen,OrbN,Hartree_k_iWn,SelfEcorr_k_iWn,GWSelfE_k_iWn_Old,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3GridNum,iWnlist,MatDim,GWMixing_SelfE)
     call PrintPol_k_iWn(Dimen,OrbN,dummy_k_iWn2,GWPol_k_iWn_Old,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist2,MatDim,GWMixing_Pol)


     print *,"---------------------------------------------------------"
     print *,"-[GW] SelfE_k_iWn.dat, [GW] Pol_k_iWn.dat files are made-"
     print *,"---------------------------------------------------------"
     GWDMFTloopNum=GWDMFTloopNum+1
     call WriteGWDMFTloopNum(GWDMFTloopNum)
     stop




endif








!---------------------------------------------!
!---------- Self-consitent GW loop -----------!
!---------------------------------------------!


if(DoGW_CTQMC==0 .and. DoCTQMC==0 .and. DoHF==0) then
     print *," "
     print *," "
     print *," ***************************************"
     print *," ******* Self-consistent GW PART *******"
     print *," ***************************************"
     print *," "
     print *," "
     print *," "

     do G=1,SCFloop           ! GW-self consistent loop, Max iteration :  SCFloop




          print *," "
          print *," ----- GW SCF # : ", G ," ----- "
          print *," "
          print *,"Calculating Hartree Self-energy "
          call Hartree_self(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,U,Vq(1,1,1,:,:),MatDim,beta,delta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Greenf_k_iWn,Hartree_k_iWn)
 
          print *,Hartree_k_iWn(1,1,1,iWnGridNum/2+1,1,1)


          print *,"Calculating Fock Self-energy "
          call Fock_self(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,U,Vq,MatDim,beta,delta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Greenf_k_iWn,Fock_k_iWn)


          print *,Fock_k_iWn(1,1,1,iWnGridNum/2+1,1,1)


          print *,"Fourier transform of G : (k,iWn)->(k,tau) "
          call k_iWn2k_tau(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Greenf_k_iWn,dummy_k_tau)



          print *,"Fourier transform of G : (k,tau)->(R,tau) "
          call k_tau2R_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_R_tau)
          call R_tau2k_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_R_tauInv)
          if (Dimen==3) then
               dummy_R_tauInv=dummy_R_tauInv/dble(K1GridNum)/dble(K2GridNum)/dble(K3GridNum) ! because there is no normalization in inverse trainsfomation
          endif
          if (Dimen==2) then
               dummy_R_tauInv=dummy_R_tauInv/dble(K1GridNum)/dble(K2GridNum) ! because there is no normalization in inverse trainsfomation
          endif
          if (Dimen==1) then
               dummy_R_tauInv=dummy_R_tauInv/dble(K1GridNum) ! because there is no normalization in inverse trainsfomation
          endif



          print *,"Polarization, (Chi) calculation in (R,tau)  "
          call PolCal(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_R_tau,dummy_R_tauInv,Polarization_R_tau)


          print *,"Fourier transform of Chi : (R,tau)->(k,tau) "
          call R_tau2k_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Polarization_R_tau,dummy_k_tau)

          print *,"Fourier transform of Chi : (k,tau)->(k,iWn) "
          call k_tau2k_iWn_Even(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist2,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_k_iWn2)
 
          print *,"Calculating Correlation part, (U(chi)U/1-(chi)U) in (k,iWn) "
          call WcorrCal(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,U,Vq,MatDim,beta,K1list,K2list,K3list,iWnlist2,Taulist,R1list,R2list,R3list,dummy_k_iWn2,ScreenedWcorr_k_iWn)



          print *,"---------Screening of coulomb interaction-----------"
          call Corrsum(Dimen,ScreenedWcorr_k_iWn,beta,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,U)
          print *,"----------------------------------------------------"

          print *,"Fourier transform of Wc : (k,iWn)->(k,tau) "
          call k_iWn2k_tau_Even(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist2,Taulist,R1list,R2list,R3list,ScreenedWcorr_k_iWn,dummy_k_tau)

          print *,"Fourier transform of Wc : (k,tau)->(R,tau) "
          call k_tau2R_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,dummy_k_tau,dummy_R_tauII)

          print *,"Calculating GW-self energy in (R,tau) space "
          call GWselfenergy(Dimen,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,dummy_R_tau,dummy_R_tauII,dummy_R_tauIII)


          print *,"Fourier transform of Self energy : (R,tau)->(k,tau) "
          call R_tau2k_taufftw(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_R_tauIII,dummy_k_tau)




          print *,"Fourier transform of Self energy : (k,tau)->(k,iWn) "
          call k_tau2k_iWn(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,dummy_k_tau,SelfEcorr_k_iWn)



          ChemPOld = ChemP
          call GWChemPgivenN(Dimen,t1,U,U1,ChemP,tq,Vq,iWnlist,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,iWnGridNum,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,Noccup,ChemOccup,MatDim,delta,beta,Ntot,ChemTH,NonIntGreenf_k_iWn,Hartree_k_iWn,Fock_k_iWn,SelfEcorr_k_iWn,dummy_k_iWn)

          call ConvCheckG (Dimen,Greenf_k_iWn,dummy_k_iWn,DiffElement,NumDiffElement,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,ConvTH)




          diffChemp=abs(ChemP-ChemPOld)
          if(  ( diffChemp < (ChemTH) .and. (   (NumDiffElement < 1)  )  )   .or. (G==SCFloop)                   )then
               print *," "
               print *,"  **************************************  "
               print *,"  ***    Green function Converged    ***  "
               print *,"  **************************************  "
               print *," "
               print *," "
               print *,"... Saving local Green function G(k,iWn) ..."
               print *," "
               print *," "
               call CallocG(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,dummy_k_iWn,locGreenf)
               call SavelocG(OrbN,DoCTQMC,localE,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,ChemP,NonIntChemPval,Noccup,beta,MatDim,iWnlist,locGreenf)

               !call SavelocselfE(locselfESave,Hartree_k_iWn,Fock_k_iWn,SelfEcorr_k_iWn,KxGridNum,KyGridNum,KzGridNum,iWnGridNum,iWnlist,MatDim)
               !endif
               print *, "!!! Files are saved as LocalG_*.dat and inout_inform !!!"
               print *," "

               if (G==SCFloop .and. SCFloop.ne.1 ) then
                    print *, " "
                    print *," **************************************** "
                    print *," ***   Convergency dosen't achieved   *** "
                    print *," **************************************** "
                    print *, " "

               endif



               exit
          endif

          ChemPOld = ChemP
          Noccup=real(ChemOccup)
          !Greenf_k_iWn=(0.5d0)*dummy_k_iWn+(0.5d0)*(Greenf_k_iWn) !Mixing
          Greenf_k_iWn=dummy_k_iWn
          print *, " "
          print *, " diffChemP : ",diffChemp, ", TH : ",(ChemTH)
          print *, " DiffElement : ",DiffElement
          print *, " NumDiffElement : ",NumDiffElement,", TH : < 1  (",ConvTH,")"
          print *, "------------------------------------------------"
          print *, " "

          if (G==SCFloop .and. SCFloop.ne.1 ) then
               print *, " "
               print *," **************************************** "
               print *," ***   Convergency dosen't achieved   *** "
               print *," **************************************** "
               print *, " "

          endif



     enddo
endif















!---------------------------------------------!
!---------- Self-consitent HF loop -----------!
!---------------------------------------------!

if(DoGW_CTQMC==0 .and. DoCTQMC==0 .and. DoHF==1) then
     print *," "
     print *," "
     print *," ***************************************"
     print *," ******* Self-consistent HF PART *******"
     print *," ***************************************"
     print *," "
     print *," "
     print *," "


     do G=1,SCFloop           ! GW-self consistent loop, Max iteration :  SCFloop

          print *," "
          print *," ----- HF SCF # : ", G ," ----- "
          print *," "
          print *,"Calculating Hartree Self-energy "
          call Hartree_self(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,U,Vq(1,1,1,:,:),MatDim,beta,delta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Greenf_k_iWn,Hartree_k_iWn)
 
          print *,"Calculating Fock Self-energy "
          call Fock_self(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,U,Vq,MatDim,beta,delta,K1list,K2list,K3list,iWnlist,Taulist,R1list,R2list,R3list,Greenf_k_iWn,Fock_k_iWn)
          print *, "-------------------------------"
          print *, "Hartree E :", maxval(real(Hartree_k_iWn(1,1,1,1,:,:)))
          print *, "Fock E    :", minval(real(Fock_k_iWn(1,1,1,1,:,:)))
          print *, "-------------------------------"

          call HFChemPgivenN(G,SpinSplit,Dimen,t1,U1,ChemP,tq,U,Vq,iWnlist,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,iWnGridNum,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,Noccup,ChemOccup,MatDim,delta,beta,Ntot,ChemTH,NonIntGreenf_k_iWn,Hartree_k_iWn,Fock_k_iWn,dummy_k_iWn)
          call ConvCheckG(Greenf_k_iWn,dummy_k_iWn,DiffElement,NumDiffElement,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,ConvTH)
          diffChemp=abs(ChemP-ChemPOld)
          if(  ( diffChemp < (ChemTH) .and. (   (NumDiffElement < 1)  )  )   .or. (G==SCFloop)                   )then
               print *," "
               print *,"  **************************************  "
               print *,"  ***    Green function Converged    ***  "
               print *,"  **************************************  "
               print *," "
               print *," "
               print *,"... Saving local Green function G(k,iWn) ..."
               print *," "
               print *," "
               call CallocG(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,dummy_k_iWn,locGreenf)
               call SavelocG(OrbN,DoCTQMC,localE,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,ChemP,NonIntChemPval,Noccup,beta,MatDim,iWnlist,locGreenf)
               
               !call SavelocselfE(locselfESave,Hartree_k_iWn,Fock_k_iWn,SelfEcorr_k_iWn,KxGridNum,KyGridNum,KzGridNum,iWnGridNum,iWnlist,MatDim)
               !endif
               print *, "!!! Files are saved as LocalG_*.dat and inout_inform !!!"
               print *," "
 
               if (G==SCFloop .and. SCFloop.ne.1 ) then
                    print *, " "
                    print *," **************************************** "
                    print *," ***   Convergency dosen't achieved   *** "
                    print *," **************************************** "
                    print *, " "
 
               endif
 
 
 
               exit
          endif

          ChemPOld = ChemP
          Noccup=real(ChemOccup)
          !Greenf_k_iWn=(0.5d0)*dummy_k_iWn+(0.5d0)*(Greenf_k_iWn) !Mixing
          Greenf_k_iWn=dummy_k_iWn
          print *, " "
          print *, " diffChemP : ",diffChemp, ", TH : ",(ChemTH)
          print *, " DiffElement : ",DiffElement
          print *, " NumDiffElement : ",NumDiffElement,", TH : < 1  (",ConvTH,")"
          print *, "------------------------------------------------"
          print *, " "


          if (G==SCFloop .and. SCFloop.ne.1 ) then
               print *, " "
               print *," **************************************** "
               print *," ***   Convergency dosen't achieved   *** "
               print *," **************************************** "
               print *, " "
 
          endif


     enddo
endif














deallocate(NeighborInfo) ! NeighborInfo(i,j,k,l) : j th neighbor of i th orbital, l=1,2,3 for cartesian coord, l=4 for orbital number
deallocate(SpinSplit) 
deallocate(NeighborDist)
deallocate(NeighborN)
deallocate(tq)
deallocate(Vq)
deallocate(R1list)
deallocate(R2list)
deallocate(R3list)
deallocate(K1list)
deallocate(K2list)
deallocate(K3list)
deallocate(Taulist)
deallocate(Taulist2)
deallocate(iWnlist)
deallocate(iWnlist2)
deallocate(t0)
deallocate(t1)
deallocate(t2)
deallocate(t3)
deallocate(U)
deallocate(U1)
deallocate(U2)
deallocate(U3)



!original allocatable
deallocate(Greenf_k_iWn)
deallocate(NonIntGreenf_k_iWn)
deallocate(NonIntGreenf0ChemP_k_iWn)
deallocate(Polarization_R_tau)
deallocate(SelfEcorr_k_iWn)
deallocate(Hartree_k_iWn)
deallocate(Fock_k_iWn)
deallocate(ScreenedWcorr_k_iWn)
deallocate(W_k_iWn)
deallocate(dummy_k_iWn)
deallocate(dummy_k_iWn2)
deallocate(dummy_k_tau)
deallocate(dummy_k_tauII)
deallocate(dummy_R_iWn)
deallocate(dummy_R_tau)
deallocate(dummy_R_tauInv)
deallocate(dummy_R_tau2)
deallocate(dummy_R_tauII)
deallocate(dummy_R_tauIII)
deallocate(locGreenf)
deallocate(locW)
deallocate(locPol)
deallocate(locSelfE)
deallocate(NonIntHeigVal)
deallocate(ImpuritySelfE)
deallocate(ImpurityPol)
deallocate(GWSelfE_k_iWn_Old)
deallocate(GWPol_k_iWn_Old)
deallocate(GWCTQMCSelfE_k_iWn)
deallocate(GWCTQMCPol_k_iWn)
deallocate(Hartree_loc)
deallocate(Fock_loc)


end program HFpGW



