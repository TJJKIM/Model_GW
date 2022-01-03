subroutine readinput


use modmain



implicit none
! local variables
logical nosym,highq
integer is,ia,ias,iostat
integer i,j,k,l,n,p
integer itind,itind2
real(8) sc,sc1,sc2,sc3
real(8) solscf,zn,b
real(8) axang(4),rot(3,3)
real(8) v1(3),v2(3)
real(8) Lfactor
character(256) block,symb,str





!!!!!!!!!!!!!!!!!!
! default values !
!!!!!!!!!!!!!!!!!!
Dimen = 2
MatDim = 2 
TauGridNum = 2000
iWnGridCut = 20.0d0
K1GridNum = 128
K2GridNum = 128
K3GridNum = 128
R1GridNum = 128
R2GridNum = 128
R3GridNum = 128
SCFloop  = 30
a1(1) = 1.0d0
a1(2) = 0.0d0
a1(3) = 0.0d0
a2(1) = 0.0d0
a2(2) = 1.0d0
a2(3) = 0.0d0
a3(1) = 0.0d0
a3(2) = 0.0d0
a3(3) = 1.0d0
Tem = 100.0d0
Ntot = 1.0d0
delta = 0.0d0
ChemTH = 0.0000000000001d0 
ConvTH = 0.0000001d0
Erange = 15.0d0
DoCTQMC = 0
DoGW_CTQMC = 0
DoHF = 0
locselfESave = 0
allgreenfSave = 0
iWnGridCheck = 0
TauGridCheck1 = 0
TauGridCheck2 = 0
LogOuton = 0
GWMixing_Pol = 1.0
GWMixing_SelfE = 1.0
DMFTMixing_Pol = 1.0
DMFTMixing_SelfE = 1.0
MuShift = 0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Read from 1dhubbard_GW.in  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(50, file='Model_GW.in', action='READ', status='OLD', form='FORMATTED',iostat=iostat)
if (iostat.ne.0) then
  write(*,*)
  write(*,'("Error(readinput): error opening 3dhubbard_GW.in")')
  write(*,*)
  stop
end if
10 continue
read(50,*,end=30) block
! check for a comment
if ((scan(trim(block),'!').eq.1).or.(scan(trim(block),'#').eq.1)) goto 10
select case(trim(block))
case('Dimension')
  read(50,*,err=20) Dimen
case('Grid')
  read(50,*,err=20) TauGridNum
  read(50,*,err=20) iWnGridCut
  read(50,*,err=20) K1GridNum
  read(50,*,err=20) K2GridNum
  read(50,*,err=20) K3GridNum
  read(50,*,err=20) R1GridNum
  read(50,*,err=20) R2GridNum
  read(50,*,err=20) R3GridNum
case('OrbitalN')
  read(50,*,err=20) OrbN
MatDim=OrbN*2
allocate(t0(MatDim,MatDim)) ! Onsite energy (orbital energy)
allocate(t1(MatDim,MatDim)) ! NN hopping parametera
allocate(t2(MatDim,MatDim)) ! 2NN hopping parameter
allocate(t3(MatDim,MatDim)) ! 3NN hopping parameter
allocate(U(MatDim,MatDim)) ! Onsite U
allocate(U1(MatDim,MatDim)) ! NN intersite U
allocate(U2(MatDim,MatDim)) ! 2NN intersite U
allocate(U3(MatDim,MatDim)) ! 3NN intersite U
allocate(OrbitPosi(OrbN,3)) !fractional
allocate(OrbitPosi_car(OrbN,3)) !cartesian
allocate(SpinSplit(OrbN))

allocate(distingAtom(OrbN,OrbN))
allocate(distingPos(OrbN,3))
case('SpinSymBreaking')
  read(50,*,err=20) SpinSplit(:)

case('OrbitalPosition')
  do itind=1,OrbN
       read(50,*,err=20) OrbitPosi(itind,:)
  enddo
case('MaxLoopNum')
  read(50,*,err=20) SCFloop
case('UnitcellVector')
  read(50,*,err=20) Lfactor
  read(50,*,err=20) a1(:)
  read(50,*,err=20) a2(:)
  read(50,*,err=20) a3(:)
case('OnsiteE')
  do itind=1,MatDim
            read(50,*,err=20) t0(itind,1:MatDim)
  enddo
case('NNHopping')
  do itind=1,MatDim
            read(50,*,err=20) t1(itind,1:MatDim)
  enddo
case('2NNHopping')
  do itind=1,MatDim
            read(50,*,err=20) t2(itind,1:MatDim)
  enddo
case('3NNHopping')
  do itind=1,MatDim
            read(50,*,err=20) t3(itind,1:MatDim)
  enddo
case('OnsiteU')
  do itind=1,MatDim
            read(50,*,err=20) U(itind,1:MatDim)
  enddo
case('NNU')
  do itind=1,MatDim
            read(50,*,err=20) U1(itind,1:MatDim)
  enddo
case('2NNU')
  do itind=1,MatDim
            read(50,*,err=20) U2(itind,1:MatDim)
  enddo
case('3NNU')
  do itind=1,MatDim
            read(50,*,err=20) U3(itind,1:MatDim)
  enddo
case('Temperature')
  read(50,*,err=20) Tem
case('Occupation')
  read(50,*,err=20) Ntot
case('delta')
  read(50,*,err=20) delta
case('ThresholdChemP')
  read(50,*,err=20) ChemTH
case('ThresholdGreenf')
  read(50,*,err=20) ConvTH
case('EnergyRange')
  read(50,*,err=20) Erange
case('GridNumCheck')
  read(50,*,err=20) iWnGridCheck
  read(50,*,err=20) TauGridCheck1
  read(50,*,err=20) TauGridCheck2
case('CTQMC')
  read(50,*,err=20) DoCTQMC
case('HartreeFock')
  read(50,*,err=20) DoHF
case('GW_CTQMC')
  read(50,*,err=20) DoGW_CTQMC
case('GreenfOutAllk')
  read(50,*,err=20) allgreenfSave
case('LocalSelfEOut')
  read(50,*,err=20) locselfESave
case('LogOutOn')
  read(50,*,err=20) LogOutOn
case('GWMixing')
  read(50,*,err=20) GWMixing_Pol
  read(50,*,err=20) GWMixing_SelfE
case('DMFTMixing')
  read(50,*,err=20) DMFTMixing_Pol
  read(50,*,err=20) DMFTMixing_SelfE
case('MuShiftInDMFT_1stloop')
  read(50,*,err=20) MuShift
case('')
  goto 10
case default
  write(*,*)
  write(*,'("Error(readinput): invalid block name : ",A)') trim(block)
  write(*,*)
  stop


end select
goto 10
20 continue
write(*,*)
write(*,'("Error(readinput): error reading from elk.in")')
write(*,'("Problem occurred in ''",A,"'' block")') trim(block)
write(*,'("Check input convention in manual")')
write(*,*)
stop
30 continue
close(50)
beta=(1.0/(boltzmanC*Tem))
dumGNum=((iWnGridCut*beta/PI)+1)/2 ! # of grid for matsubara frequency
iWnGridNum=dumGNum*2 ! # of grid for matsubara frequency
lwork=5*MatDim-1

do itind=1,OrbN
      OrbitPosi_car(itind,:)=a1(:)*OrbitPosi(itind,1) + a2(:)*OrbitPosi(itind,2) + a3(:)*OrbitPosi(itind,3)
enddo



return
end subroutine
