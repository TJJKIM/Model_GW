!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for saving eigen values of NonInt Hamiltonian !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SaveNonInteigV(Dimen,K1GridNum,K2GridNum,K3GridNum,K1list,K2list,K3list,MatDim,ChemP,localE,NonIntHeigVal)

Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum,MatDim,Dimen
real(8), intent(in) :: K1list(K1GridNum,3), K2list(K2GridNum,3), K3list(K3GridNum,3), ChemP
complex(8), intent(in) ::  NonIntHeigVal(K1GridNum,K2GridNum,K3GridNum,MatDim,1)
complex(8), intent(inout) :: localE(MatDim)

! local variable
integer kx,ky,kz,n

localE(:)=(0.0d0,0.0d0) ! initialize
open (unit=2, file ="NonIntHeigV.dat")
!open (unit=3, file ="NonIntHeigV.gnuplot")
!write(2,*), "plot "NonIntHeigV" "

if (Dimen==3) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do kz=1,K3GridNum
                    localE(:)=localE(:)+(NonIntHeigVal(kx,ky,kz,:,1)/(K1GridNum*K2GridNum*K3GridNum*1.0d0))
               enddo
          enddo
     enddo

     write(2,3) , real(localE(:))
     3 format (30(f))

     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do kz=1,K3GridNum
                    write(2,4) , real(K1list(kx,1)+K2list(ky,1)+K3list(kz,1)), real(K1list(kx,2)+K2list(ky,2)+K3list(kz,2)), real(K1list(kx,3)+K2list(ky,3)+K3list(kz,3)), real(NonIntHeigVal(kx,ky,kz,:,1))
                    4 format (30(f))
               enddo
          enddo
     enddo


endif

if (Dimen==2) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               localE(:)=localE(:)+(NonIntHeigVal(kx,ky,1,:,1)/(K1GridNum*K2GridNum*1.0d0))
          enddo
     enddo

     write(2,5) , real(localE(:))
     5 format (30(f))


     do kx=1,K1GridNum
          do ky=1,K2GridNum
               write(2,6) , real(K1list(kx,1)+K2list(ky,1)), real(K1list(kx,2)+K2list(ky,2)), real(K1list(kx,3)+K2list(ky,3)), real(NonIntHeigVal(kx,ky,1,:,1))
               6 format (30(f))
          enddo
     enddo

endif

if (Dimen==1) then
     do kx=1,K1GridNum
          localE(:)=localE(:)+(NonIntHeigVal(kx,1,1,:,1)/(K1GridNum*1.0d0))
     enddo


     write(2,7) , real(localE(:))
     7 format (30(f))


     do kx=1,K1GridNum
          write(2,8) , real(K1list(kx,1)),real(K1list(kx,2)),real(K1list(kx,3)), real(NonIntHeigVal(kx,1,1,:,1))
          8 format (30(f))
     enddo

endif



print *, " "
print *, "************ local energy of orbitals ************"
print 10, localE(:)
10 format (30(f10.6,f10.6))
print *, "************ local energy of orbitals ************"
print *, " "
close(2)
!close(3)
return
end subroutine










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for saving local greenfunction as files !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine SavelocG(OrbN,DoCTQMC,localE,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,ChemP,NonIntChemPval,Noccup,beta,MatDim,iWnlist,locGreenf)
Implicit none

! argument
integer, intent(in) :: OrbN,MatDim,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,DoCTQMC
real(8), intent(in) :: ChemP,beta,Noccup(MatDim),NonIntChemPval
complex(8), intent(in) :: iWnlist(iWnGridNum),locGreenf(iWnGridNum,MatDim,MatDim),localE(MatDim)

!local variable
integer i,kx,ky,kz,n,Z,j,k,i1,i2
complex(8) dummy(MatDim,MatDim),dummy2(MatDim,MatDim),localEMat(MatDim,MatDim),hybrid(MatDim,MatDim), Iden(MatDim,MatDim), instG(OrbN), insthyb(OrbN), testG(MatDim,MatDim), testhyb(MatDim,MatDim)

Z=(iWnGridNum/2)

open (unit=2, file ="LocalG.dat")
open (unit=3, file ="LocalG_diag.dat")
if (DoCTQMC==1 .or. DoCTQMC==2) then
     open (unit=8, file ="hybridization.dat")
          write(8,12,advance="no"), "iWn"
          12 format(1x,A6)
          write(8,13,advance="no"), ""
          13 format(1x,A3)

          do i1=1,MatDim !  
               do i2=1,MatDim !  
                    write(8,14,advance="no"), "hyb_",i1,"_",i2,"_Re"
                    14 format(1x,A14,I2,1x,A,I2,1x,A)
                         write(8,15,advance="no"), "hyb_",i1,"_",i2,"_Im"
                         15 format(1x,A14,I2,1x,A,I2,1x,A)
               enddo
          enddo
          write(8,*) 


          write(8,16), beta, Noccup, NonIntChemPval
          16 format (30(f))
endif

call Identity(MatDim,Iden)

if (DoCTQMC==1 .or. DoCTQMC==2) then
     do i=1,MatDim
          localEMat(i,i)=localE(i)
     enddo
endif



write(2,17,advance="no"), "iWn"
17 format(1x,A6)

write(2,18,advance="no"), ""
18 format(1x,A3)


write(3,21,advance="no"), "iWn"
21 format(1x,A6)

write(3,22,advance="no"), ""
22 format(1x,A3)


!Not TEST
!do i1=1,OrbN ! In formal hubbard model, Only diagonal term survive in local Greenfunction, also spin up,down part exactly same 
!     write(2,4,advance="no"), "GreenF_",i1,"_",i1,"_Re"
!     4 format(1x,A14,I2,1x,A,I2,1x,A)
!          write(2,15,advance="no"), "GreenF_",i1,"_",i1,"_Im"
!          15 format(1x,A14,I2,1x,A,I2,1x,A)
!enddo
!Not TEST

!TEST
do i1=1,MatDim ! In formal hubbard model, Only diagonal term survive in local Greenfunction, also spin up,down part exactly same 
     write(3,23,advance="no"), "GreenF_",i1,"_",i1,"_Re"
     23 format(1x,A14,I2,1x,A,I2,1x,A)
          write(3,24,advance="no"), "GreenF_",i1,"_",i1,"_Im"
          24 format(1x,A14,I2,1x,A,I2,1x,A)

     do i2=1,MatDim ! In formal hubbard model, Only diagonal term survive in local Greenfunction, also spin up,down part exactly same 
          write(2,25,advance="no"), "GreenF_",i1,"_",i2,"_Re"
          25 format(1x,A14,I2,1x,A,I2,1x,A)
               write(2,26,advance="no"), "GreenF_",i1,"_",i2,"_Im"
               26 format(1x,A14,I2,1x,A,I2,1x,A)
     enddo
enddo
!TEST


write(2,*)
write(3,*)

do n=1,iWnGridNum ! iWn Sum
     
     !Not TEST    
     !do j=1,OrbN
     !     instG(j)=locGreenf(n,j,j)
     !enddo
     !if (n>Z) then ! -iWn correspond to conj(iWn)
     !     write(2,2) ,aimag(iWnlist(n)), instG
     !     2 format (f10.6,30(f,f))
     !endif
     !Not TEST

     !TEST
     testG=locGreenf(n,:,:)
     if(n>Z) then
          write(2,27,advance="no"), aimag(iWnlist(n))
          27 format (f)
          write(3,28,advance="no"), aimag(iWnlist(n))
          28 format (f)

          do i1=1,MatDim
               write(3,31, advance="no"), testG(i1,i1)
               31 format (30(f,f))
               do i2=1,MatDim
                    write(2,32, advance="no"), testG(i1,i2)
                    32 format (30(f,f))
               enddo
          enddo
          write(2,*)
          write(3,*)
     endif
     !TEST
 


     if ( (DoCTQMC==1 .or. DoCTQMC==2) .and. n>Z ) then ! -iWn = conj(f(iWn))
          call MatInv(MatDim,locGreenf(n,:,:),dummy2)
          hybrid(:,:)=(iWnlist(n)*Iden)+(ChemP*Iden)-dummy2-localEMat ! It is easy to know this is diagonal, so we consider only diagonal term and spin up = down 
          write(8,33,advance="no"), aimag(iWnlist(n))
          33 format (f)
          do i1=1,MatDim
               do i2=1,MatDim
               write(8,34, advance="no"), hybrid(i1,i2)
               34 format (30(f,f))
               enddo
          enddo
          write(8,*)
     endif
enddo

close(3)
close(2)
if (DoCTQMC==1 .or. DoCTQMC==2) then
     close(8)
endif


return
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for printing SelfE_k_iWn, Pol_k_iWn !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PrintSelfE_k_iWn_Init(Dimen,OrbN,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist,MatDim)

Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum, iWnGridNum, MatDim, Dimen, OrbN
real(8), intent(in) :: K1list(K1GridNum,3), K2list(K2GridNum,3), K3list(K3GridNum,3)
complex(8), intent(in) :: iWnlist(iWnGridNum)

!local variable
integer kx,ky,kz, n, Z,i1,i2

open (unit=2, file ="SelfE_k_iWn.dat")

write(2,35, advance="no"), "kx", "ky", "kz", "iWn",""
35 format(1x,A6,1x,A9,1x,A9,1x,A15,1x,A11)
do i1=1,MatDim
     do i2=1,MatDim
          write(2,36, advance="no"), "Self_",i1,"_",i2,"_Re"
          36 format(1x,A14,I2,1x,A,I2,1x,A) 
          write(2,37, advance="no"), "Self_",i1,"_",i2,"_Im"
          37 format(1x,A14,I2,1x,A,I2,1x,A)
     enddo
enddo

write(2,*),

Z=(iWnGridNum/2)

if(Dimen==3) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do kz=1,K3GridNum
                    do n=1,iWnGridNum
                         if (n>Z) then
                              write(2,38,advance="no") , real(K1list(kx,1)+K2list(ky,1)+K3list(kz,1)), real(K1list(kx,2)+K2list(ky,2)+K3list(kz,2)), real(K1list(kx,3)+K2list(ky,3)+K3list(kz,3))
                              38 format (f10.6,f10.6,f10.6)
                              write(2,41,advance="no"), aimag(iWnlist(n))
                              41 format (f)
                              do i1=1,MatDim
                                   do i2=1,MatDim
                                        write(2,42,advance="no"), 0.0d0, 0.0d0
                                        42 format (f,f)
                                   enddo
                              enddo
                              write(2,*)
                         endif
                    enddo
               enddo
          enddo
     enddo
endif


if(Dimen==2) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do n=1,iWnGridNum
                    if(n>z) then
                         !write(2,*)
                         write(2,43,advance="no") , real(K1list(kx,1)+K2list(ky,1)), real(K1list(kx,2)+K2list(ky,2)), real(K1list(kx,3)+K2list(ky,3))
                         43 format (f10.6,f10.6,f10.6)
                         write(2,44,advance="no"), aimag(iWnlist(n))
                         44 format (f)
                         do i1=1,MatDim
                              do i2=1,MatDim
                                   write(2,45,advance="no"), (0.0d0, 0.0d0)
                                   45 format (f,f)
                              enddo
                         enddo
                         write(2,*)
                    endif
               enddo
          enddo
     enddo
endif


if(Dimen==1) then
     do kx=1,K1GridNum
          do n=1,iWnGridNum
               if(n>z) then
                    write(2,46,advance="no") , real(K1list(kx,1)), real(K1list(kx,2)), real(K1list(kx,3))
                    46 format (f10.6,f10.6,f10.6)
                    write(2,47,advance="no"), aimag(iWnlist(n))
                    47 format (f)
                    do i1=1,MatDim
                         do i2=1,MatDim
                              write(2,48,advance="no"), 0.0d0, 0.0d0
                              48 format (f,f)
                         enddo
                    enddo
                    write(2,*)
               endif
          enddo
     enddo
endif


close(2)
return

end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for printing SelfE_k_iWn, Pol_k_iWn !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PrintPol_k_iWn_Init(Dimen,OrbN,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist2,MatDim)

Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum, iWnGridNum, MatDim,Dimen, OrbN
real(8), intent(in) :: K1list(K1GridNum,3), K2list(K2GridNum,3), K3list(K3GridNum,3)
complex(8), intent(in) :: iWnlist2(iWnGridNum+1)

!local variable
integer kx,ky,kz, m, Z,i1,i2

open (unit=2, file ="Pol_k_iWn.dat")
write(2,49, advance="no"), "kx", "ky", "kz", "iWn",""
49 format(1x,A6,1x,A9,1x,A9,1x,A15,1x,A11)
do i1=1,MatDim
     do i2=1,MatDim
          write(2,51, advance="no"), "Pol__",i1,"_",i2,"_Re"
          51 format(1x,A14,I2,1x,A,I2,1x,A)
          write(2,52, advance="no"), "Pol__",i1,"_",i2,"_Im"
          52 format(1x,A14,I2,1x,A,I2,1x,A)
     enddo
enddo

write(2,*),


Z=(iWnGridNum/2)

if(Dimen==3) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do kz=1,K3GridNum
                    do m=1,iWnGridNum+1
                         if (m>Z) then
                              write(2,53,advance="no") , real(K1list(kx,1)+K2list(ky,1)+K3list(kz,1)), real(K1list(kx,2)+K2list(ky,2)+K3list(kz,2)), real(K1list(kx,3)+K2list(ky,3)+K3list(kz,3))
                              53 format (f10.6,f10.6,f10.6)
                              write(2,54,advance="no"), aimag(iWnlist2(m))
                              54 format (f)
                              do i1=1,MatDim
                                   do i2=1,MatDim
                                        write(2,55,advance="no"), 0.0d0, 0.0d0
                                        55 format (f,f)
                                   enddo
                              enddo
                              write(2,*)
                         endif
                    enddo
               enddo
          enddo
     enddo
endif


if(Dimen==2) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do m=1,iWnGridNum+1
                    if(m>Z) then
                         write(2,56,advance="no") , real(K1list(kx,1)+K2list(ky,1)), real(K1list(kx,2)+K2list(ky,2)), real(K1list(kx,3)+K2list(ky,3))
                         56 format (f10.6,f10.6,f10.6)
                         write(2,57,advance="no"), aimag(iWnlist2(m))
                         57 format (f)
                         do i1=1,MatDim
                              do i2=1,MatDim
                                   write(2,58,advance="no"), 0.0d0, 0.0d0
                                   58 format (f,f)
                              enddo
                         enddo
                         write(2,*)
                    endif
               enddo
          enddo
     enddo
endif


if(Dimen==1) then
     do kx=1,K1GridNum
          do m=1,iWnGridNum+1
               if(m>Z) then
                    write(2,59,advance="no") , real(K1list(kx,1)), real(K1list(kx,2)), real(K1list(kx,3))
                    59 format (f10.6,f10.6,f10.6)
                    write(2,60,advance="no"), aimag(iWnlist2(m))
                    60 format (f)
                    do i1=1,MatDim
                         do i2=1,MatDim
                              write(2,61,advance="no"), 0.0d0, 0.0d0
                              61 format (f,f)
                         enddo
                    enddo
                    write(2,*)
               endif
          enddo
     enddo
endif



close(2)
return
end subroutine



subroutine Save_tqVq(Dimen,K1GridNum,K2GridNum,K3GridNum,K1list,K2list,K3list,MatDim,tq,Vq)

Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum,MatDim,Dimen
real(8), intent(in) :: K1list(K1GridNum,3), K2list(K2GridNum,3), K3list(K3GridNum,3)
complex(8), intent(in) ::  tq(K1GridNum,K2GridNum,K3GridNum,MatDim,MatDim), Vq(K1GridNum,K2GridNum,K3GridNum,MatDim,MatDim)

! local variable
integer kx,ky,kz,o1,o2

open (unit=100, file ="tq.dat")
open (unit=101, file ="Vq.dat")


do kx=1,K1GridNum
     do ky=1,K2GridNum
          do kz=1,K3GridNum
               write(100,62), "Kpt :", K1list(kx,:)+K2list(ky,:)+K3list(kz,:)
               write(101,63), "Kpt :", K1list(kx,:)+K2list(ky,:)+K3list(kz,:)
               62 format(1x,A,f8.4,f8.4,f8.4)
               63 format(1x,A,f8.4,f8.4,f8.4)
               write(100, *), "tq :"              
               write(101, *), "Vq :"              
               do o1=1,MatDim
                         write(100,64) (tq(kx,ky,kz,o1,o2), ",", o2=1,MatDim)
                         64 format(50((f8.4,f8.4),1x,A))
                         write(101,65) (Vq(kx,ky,kz,o1,o2), ",", o2=1,MatDim)
                         65 format(50((f8.4,f8.4),1x,A))
               enddo
               write(100,*) 
               write(101,*) 
          enddo
     enddo
enddo



close(100)
close(101)
return
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for printing SelfE_k_iWn, Pol_k_iWn !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PrintGWEDMFTSelfE(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist,MatDim,GWCTQMCSelfE_k_iWn)

Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum, iWnGridNum, MatDim, Dimen
real(8), intent(in) :: K1list(K1GridNum,3), K2list(K2GridNum,3), K3list(K3GridNum,3)
complex(8), intent(in) :: iWnlist(iWnGridNum), GWCTQMCSelfE_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)

!local variable
integer i1,i2,kx,ky,kz, n, Z
complex(8) SelfE

open (unit=2, file ="GWEDMFTSelfE.dat")


write(2,66, advance="no"), "kx", "ky", "kz", "iWn",""
66 format(1x,A6,1x,A9,1x,A9,1x,A15,1x,A11)
do i1=1,MatDim
     do i2=1,MatDim
          write(2,67, advance="no"), "Self_",i1,"_",i2,"_Re"
          67 format(1x,A14,I2,1x,A,I2,1x,A)
          write(2,68, advance="no"), "Self_",i1,"_",i2,"_Im"
          68 format(1x,A14,I2,1x,A,I2,1x,A)
     enddo
enddo

write(2,*),

Z=(iWnGridNum/2)



!Z=(iWnGridNum/2)+1
!print *, "HI"




if(Dimen==3) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do kz=1,K3GridNum
                    do n=1,iWnGridNum
                         if (n>Z) then
                              write(2,69,advance="no") , real(K1list(kx,1)+K2list(ky,1)+K3list(kz,1)), real(K1list(kx,2)+K2list(ky,2)+K3list(kz,2)), real(K1list(kx,3)+K2list(ky,3)+K3list(kz,3))
                              69 format (f10.6,f10.6,f10.6)
                              write(2,71,advance="no"), aimag(iWnlist(n))
                              71 format (f)
                              do i1=1,MatDim
                                   do i2=1,MatDim
                                        write(2,72,advance="no"), real(GWCTQMCSelfE_k_iWn(kx,ky,kz,n,i1,i2)), aimag(GWCTQMCSelfE_k_iWn(kx,ky,kz,n,i1,i2))
                                        72 format (f,f)
                                   enddo
                              enddo
                              write(2,*)
                         endif
                    enddo
               enddo
          enddo
     enddo
endif



if(Dimen==2) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do n=1,iWnGridNum
                    if (n>Z) then
                         write(2,73,advance="no") , real(K1list(kx,1)+K2list(ky,1)), real(K1list(kx,2)+K2list(ky,2)), real(K1list(kx,3)+K2list(ky,3))
                         73 format (f10.6,f10.6,f10.6)
                         write(2,74,advance="no"), aimag(iWnlist(n))
                         74 format (f)
                         do i1=1,MatDim
                              do i2=1,MatDim
                                   write(2,75,advance="no"), real(GWCTQMCSelfE_k_iWn(kx,ky,1,n,i1,i2)), aimag(GWCTQMCSelfE_k_iWn(kx,ky,1,n,i1,i2))
                                   75 format (f,f)
                              enddo
                         enddo
                         write(2,*)
                    endif
               enddo
          enddo
     enddo
endif



if(Dimen==1) then
     do kx=1,K1GridNum
          do n=1,iWnGridNum
               if (n>Z) then
                    write(2,76,advance="no") , real(K1list(kx,1)), real(K1list(kx,2)), real(K1list(kx,3))
                    76 format (f10.6,f10.6,f10.6)
                    write(2,77,advance="no"), aimag(iWnlist(n))
                    77 format (f)
                    do i1=1,MatDim
                         do i2=1,MatDim
                              write(2,78,advance="no"), real(GWCTQMCSelfE_k_iWn(kx,1,1,n,i1,i2)), aimag(GWCTQMCSelfE_k_iWn(kx,1,1,n,i1,i2))
                              78 format (f,f)
                         enddo
                    enddo
                    write(2,*)
               endif
          enddo
     enddo
endif



close(2)
return
end subroutine




subroutine PrintGWEDMFTPol(Dimen,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist,MatDim,GWCTQMCPol_k_iWn)

Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum, iWnGridNum, MatDim, Dimen
real(8), intent(in) :: K1list(K1GridNum,3), K2list(K2GridNum,3), K3list(K3GridNum,3)
complex(8), intent(in) :: iWnlist(iWnGridNum), GWCTQMCPol_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)

!local variable
integer i1,i2,kx,ky,kz, n, Z
complex(8) Pol

open (unit=2, file ="GWEDMFTPol.dat")


write(2,79, advance="no"), "kx", "ky", "kz", "iWn",""
79 format(1x,A6,1x,A9,1x,A9,1x,A15,1x,A11)
do i1=1,MatDim
     do i2=1,MatDim
          write(2,81, advance="no"), "Pol__",i1,"_",i2,"_Re"
          81 format(1x,A14,I2,1x,A,I2,1x,A)
          write(2,82, advance="no"), "Pol__",i1,"_",i2,"_Im"
          82 format(1x,A14,I2,1x,A,I2,1x,A)
     enddo
enddo

write(2,*),

Z=(iWnGridNum/2)



!Z=(iWnGridNum/2)+1
!print *, "HI"

if(Dimen==3) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do kz=1,K3GridNum
                    do n=1,iWnGridNum
                         if (n>Z) then
                              write(2,83,advance="no") , real(K1list(kx,1)+K2list(ky,1)+K3list(kz,1)), real(K1list(kx,2)+K2list(ky,2)+K3list(kz,2)), real(K1list(kx,3)+K2list(ky,3)+K3list(kz,3))
                              83 format (f10.6,f10.6,f10.6)
                              write(2,84,advance="no"), aimag(iWnlist(n))
                              84 format (f)
                              do i1=1,MatDim
                                   do i2=1,MatDim
                                        write(2,85,advance="no"), real(GWCTQMCPol_k_iWn(kx,ky,kz,n,i1,i2)), aimag(GWCTQMCPol_k_iWn(kx,ky,kz,n,i1,i2))
                                        85 format (f,f)
                                   enddo
                              enddo
                              write(2,*)
                         endif
                    enddo
               enddo
          enddo
     enddo
endif



if(Dimen==2) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do n=1,iWnGridNum
                    if (n>Z) then
                         write(2,86,advance="no") , real(K1list(kx,1)+K2list(ky,1)), real(K1list(kx,2)+K2list(ky,2)), real(K1list(kx,3)+K2list(ky,3))
                         86 format (f10.6,f10.6,f10.6)
                         write(2,87,advance="no"), aimag(iWnlist(n))
                         87 format (f)
                         do i1=1,MatDim
                              do i2=1,MatDim
                                   write(2,88,advance="no"), real(GWCTQMCPol_k_iWn(kx,ky,1,n,i1,i2)), aimag(GWCTQMCPol_k_iWn(kx,ky,1,n,i1,i2))
                                   88 format (f,f)
                              enddo
                         enddo
                         write(2,*)
                    endif
               enddo
          enddo
     enddo
endif


if(Dimen==1) then
     do kx=1,K1GridNum
          do n=1,iWnGridNum
               if (n>Z) then
                    write(2,89,advance="no") , real(K1list(kx,1)), real(K1list(kx,2)), real(K1list(kx,3))
                    89 format (f10.6,f10.6,f10.6)
                    write(2,91,advance="no"), aimag(iWnlist(n))
                    91 format (f)
                    do i1=1,MatDim
                         do i2=1,MatDim
                              write(2,92,advance="no"), real(GWCTQMCPol_k_iWn(kx,1,1,n,i1,i2)), aimag(GWCTQMCPol_k_iWn(kx,1,1,n,i1,i2))
                              92 format (f,f)
                         enddo
                    enddo
                    write(2,*)
               endif
          enddo
     enddo
endif



close(2)
return
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating hybridization and priting at "hybridization.dat" !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine  Printhybridization(locGreenf,locSelfE,MatDim,iWnGridNum,iWnlist,ChemP,beta,Noccup,U)
Implicit none
! argument
integer, intent(in) :: iWnGridNum,MatDim
real(8), intent(in) :: ChemP, beta, Noccup(MatDim),U(MatDim,MatDim)
complex(8), intent(in) :: iWnlist(iWnGridNum), locGreenf(iWnGridNum,MatDim,MatDim), locSelfE(iWnGridNum,MatDim,MatDim)


!local variable
integer i,k,n,i1,i2
complex(8) dummy(MatDim,MatDim),dummy2(MatDim,MatDim),localEMat(MatDim,MatDim),hybrid(MatDim,MatDim), Iden(MatDim,MatDim)






do i=1,MatDim
     localEMat(i,i)=-U(i,i)/(2.0d0)
enddo

call Identity(MatDim,Iden)
open (unit=8, file ="hybridization.dat")

write(8,93), beta, Noccup, ChemP
93 format (30(f))


write(8,94, advance="no"), "iWn",""
94 format(1x,A15,1x,A7)
do i1=1,MatDim
     do i2=1,MatDim
          write(8,95, advance="no"), "hyb__",i1,"_",i2,"_Re"
          95 format(1x,A14,I2,1x,A,I2,1x,A)
          write(8,96, advance="no"), "hyb__",i1,"_",i2,"_Im"
          96 format(1x,A14,I2,1x,A,I2,1x,A)
     enddo
enddo

write(8,*),

do n=1,iWnGridNum ! iWn Sum
     if ( aimag(iWnlist(n))>0 ) then
  
          call MatInv(MatDim,locGreenf(n,:,:),dummy2)
          hybrid(:,:)=(iWnlist(n)*Iden)+(ChemP*Iden)-dummy2-(matmul(locSelfE(n,:,:),Iden))-localEMat
  
          write(8,97,advance="no") ,aimag(iWnlist(n))
          97 format (f)
  
          do i1=1,MatDim
               do i2=1,MatDim
                    write(8,98,advance="no"), real(hybrid(i1,i2)), aimag(hybrid(i1,i2))
                    98 format (f,f)
               enddo
          enddo
          write(8,*)
  
  
     end if
enddo






close(8)

return
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating nonlocal Retarded Interaction at "NonlocRtareded.dat" !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine PrintRetardedInt(locW,locPol,iWnGridNum,MatDim,U,iWnlist2,GWDMFTloopNum)
Implicit none
! argument
integer, intent(in) :: iWnGridNum,MatDim,GWDMFTloopNum
real(8), intent(in) :: U(MatDim,MatDim)
complex(8), intent(in) :: iWnlist2(iWnGridNum+1), locW(iWnGridNum+1,MatDim,MatDim), locPol(iWnGridNum+1,MatDim,MatDim)


!local variable
integer n,i1,i2
complex(8) Iden(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim), dummy3(MatDim,MatDim),  NonlocRetarded(MatDim,MatDim), RetardedInt(MatDim,MatDim)


call Identity(MatDim,Iden)


open (unit=8, file ="NonlocRetarded.dat")
open (unit=9, file ="RetrdInteraction.dat")


if (GWDMFTloopNum>0 .or. GWDMFTloopNum==0) then
     do n=1,(iWnGridNum+1) ! iWn Sum
           if ( (aimag(iWnlist2(n))>0.0) .or. (aimag(iWnlist2(n))==0.0) ) then
               
                call MatInv(MatDim,locW(n,:,:),dummy)
                call MatInv(MatDim,dummy+locPol(n,:,:),dummy2)

                NonlocRetarded= dummy2 - U
                RetardedInt= dummy2


                write(8,101,advance="no") ,aimag(iWnlist2(n))
                101 format (f)
 
                do i1=1,MatDim
                     do i2=1,MatDim
                          write(8,102,advance="no"), real(NonLocRetarded(i1,i2)), aimag(NonlocRetarded(i1,i2))
                          102 format (f,f)
                     enddo
                enddo
                write(8,*)


                write(9,103,advance="no") ,aimag(iWnlist2(n))
                103 format (f)

                do i1=1,MatDim
                     do i2=1,MatDim
                          write(9,104,advance="no"), real(RetardedInt(i1,i2)), aimag(RetardedInt(i1,i2))
                          104 format (f,f)
                     enddo
                enddo
                write(9,*)




           endif
     enddo
endif

close(8)
close(9)
return
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for printing SelfE_k_iWn, Pol_k_iWn !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Printdyn_iWn(locW,locPol,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,iWnlist2,U,MatDim,logNum)

Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum, iWnGridNum, MatDim,logNum
real(8), intent(in) :: U(MatDim,MatDim)
complex(8), intent(in) ::locW(iWnGridNum+1,MatDim,MatDim),locPol(iWnGridNum+1,MatDim,MatDim), iWnlist2(iWnGridNum+1)

!local variable
integer m, Z,i1, i2
complex(8)  Iden(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim), dummy3(MatDim,MatDim), dyn(MatDim,MatDim) 
real (8) Vq
character*20 :: filename
WRITE(filename,'(a,i4.4,a)') "GWdyn_",logNum,".dat"

open (unit=2, file =filename)
!Z=(iWnGridNum/2)+1

if (logNum>1 .or. logNum==1) then
     do m=1,(iWnGridNum+1) ! iWn Sum
           if ( (aimag(iWnlist2(m))>0.0) .or. (aimag(iWnlist2(m))==0.0) ) then

                call MatInv(MatDim,locW(m,:,:),dummy)
                call MatInv(MatDim,dummy+locPol(m,:,:),dummy2)

                dyn= dummy2 - U


                write(2,105,advance="no") ,aimag(iWnlist2(m))
                105 format (f)

                do i1=1,MatDim
                     do i2=1,MatDim
                          write(2,106,advance="no"), real(dyn(i1,i2)), aimag(dyn(i1,i2))
                          106 format (f,f)
                     enddo
                enddo


           endif
     enddo
endif




close(2)
return

end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for printing SelfE_k_iWn, Pol_k_iWn !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PrintSelfE_k_iWn(Dimen,OrbN,Hartree_k_iWn,SelfEcorr_k_iWn,GWSelfE_k_iWn_Old,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist,MatDim,GWMixing_SelfE)

Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum, iWnGridNum, MatDim, Dimen, OrbN
real(8), intent(in) :: K1list(K1GridNum,3), K2list(K2GridNum,3), K3list(K3GridNum,3),GWMixing_SelfE
complex(8), intent(in) :: Hartree_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim),SelfEcorr_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim),GWSelfE_k_iWn_Old(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim), iWnlist(iWnGridNum)

!local variable
integer kx,ky,kz, n, Z,i1,i2
complex(8) SelfE_New(MatDim,MatDim), SelfE_Old(MatDim,MatDim), SelfE(MatDim,MatDim)


open (unit=2, file ="SelfE_k_iWn.dat")

write(2,107, advance="no"), "kx", "ky", "kz", "iWn",""
107 format(1x,A6,1x,A9,1x,A9,1x,A15,1x,A11)
do i1=1,MatDim
     do i2=1,MatDim
          write(2,108, advance="no"), "Self_",i1,"_",i2,"_Re"
          108 format(1x,A14,I2,1x,A,I2,1x,A)
          write(2,109, advance="no"), "Self_",i1,"_",i2,"_Im"
          109 format(1x,A14,I2,1x,A,I2,1x,A)
     enddo
enddo

write(2,*),

Z=(iWnGridNum/2)

if(Dimen==3) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do kz=1,K3GridNum
                    do n=1,iWnGridNum
                         if (n>Z) then
                              write(2,111,advance="no") , real(K1list(kx,1)+K2list(ky,1)+K3list(kz,1)), real(K1list(kx,2)+K2list(ky,2)+K3list(kz,2)), real(K1list(kx,3)+K2list(ky,3)+K3list(kz,3))
                              111 format (f10.6,f10.6,f10.6)
                              write(2,112,advance="no"), aimag(iWnlist(n))
                              112 format (f)

                              SelfE_Old=GWSelfE_k_iWn_Old(kx,ky,kz,n,:,:)
                              SelfE_New=Hartree_k_iWn(kx,ky,kz,n,:,:)+SelfEcorr_k_iWn(kx,ky,kz,n,:,:)
                              SelfE= SelfE_Old*(1.0d0-GWMixing_SelfE) + SelfE_New*(GWMixing_SelfE)


                              do i1=1,MatDim
                                   do i2=1,MatDim
                                        write(2,113,advance="no"), real(SelfE(i1,i2)), aimag(SelfE(i1,i2))
                                        113 format (f,f)
                                   enddo
                              enddo
                              write(2,*)
                         endif
                    enddo
               enddo
          enddo
     enddo
endif



if(Dimen==2) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
              do n=1,iWnGridNum
                   if (n>Z) then
                        write(2,114,advance="no") , real(K1list(kx,1)+K2list(ky,1)), real(K1list(kx,2)+K2list(ky,2)), real(K1list(kx,3)+K2list(ky,3))
                        114 format (f10.6,f10.6,f10.6)
                        write(2,115,advance="no"), aimag(iWnlist(n))
                        115 format (f)

                        SelfE_Old=GWSelfE_k_iWn_Old(kx,ky,1,n,:,:)
                        SelfE_New=Hartree_k_iWn(kx,ky,1,n,:,:)+SelfEcorr_k_iWn(kx,ky,1,n,:,:)
                        SelfE= SelfE_Old*(1.0d0-GWMixing_SelfE) + SelfE_New*(GWMixing_SelfE)


                        do i1=1,MatDim
                             do i2=1,MatDim
                                  write(2,116,advance="no"), real(SelfE(i1,i2)), aimag(SelfE(i1,i2))
                                  116 format (f,f)
                             enddo
                        enddo
                        write(2,*)
                   endif
              enddo
          enddo
     enddo
endif



if(Dimen==1) then
     do kx=1,K1GridNum
         do n=1,iWnGridNum
              if (n>Z) then
                   write(2,117,advance="no") , real(K1list(kx,1)), real(K1list(kx,2)), real(K1list(kx,3))
                   117 format (f10.6,f10.6,f10.6)
                   write(2,118,advance="no"), aimag(iWnlist(n))
                   118 format (f)

                   SelfE_Old=GWSelfE_k_iWn_Old(kx,1,1,n,:,:)
                   SelfE_New=Hartree_k_iWn(kx,1,1,n,:,:)+SelfEcorr_k_iWn(kx,1,1,n,:,:)
                   SelfE= SelfE_Old*(1.0d0-GWMixing_SelfE) + SelfE_New*(GWMixing_SelfE)


                   do i1=1,MatDim
                        do i2=1,MatDim
                             write(2,121,advance="no"), real(SelfE(i1,i2)), aimag(SelfE(i1,i2))
                             121 format (f,f)
                        enddo
                   enddo
                   write(2,*)
              endif
         enddo
     enddo
endif



close(2)
return

end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for printing SelfE_k_iWn, Pol_k_iWn !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PrintPol_k_iWn(Dimen,OrbN,Polarization_k_iWn,GWPol_k_iWn_Old,K1GridNum,K2GridNum,K3GridNum,iWnGridNum,K1list,K2list,K3list,iWnlist2,MatDim,GWMixing_Pol)

Implicit none
! argument
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum, iWnGridNum, MatDim, Dimen, OrbN
real(8), intent(in) :: K1list(K1GridNum,3), K2list(K2GridNum,3), K3list(K3GridNum,3),GWMixing_Pol
complex(8), intent(in) :: Polarization_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim),GWPol_k_iWn_Old(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim), iWnlist2(iWnGridNum+1)

!local variable
integer kx,ky,kz, n, Z,i1,i2
complex(8) Pol_New(MatDim,MatDim), Pol_Old(MatDim,MatDim), Pol(MatDim,MatDim)


open (unit=2, file ="Pol_k_iWn.dat")

write(2,122, advance="no"), "kx", "ky", "kz", "iWn",""
122 format(1x,A6,1x,A9,1x,A9,1x,A15,1x,A11)
do i1=1,MatDim
     do i2=1,MatDim
          write(2,123, advance="no"), "Pol__",i1,"_",i2,"_Re"
          123 format(1x,A14,I2,1x,A,I2,1x,A)
          write(2,124, advance="no"), "Pol__",i1,"_",i2,"_Im"
          124 format(1x,A14,I2,1x,A,I2,1x,A)
     enddo
enddo

write(2,*),

Z=(iWnGridNum/2)

if(Dimen==3) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do kz=1,K3GridNum
                    do n=1,iWnGridNum
                         if (n>Z) then
                              write(2,125,advance="no") , real(K1list(kx,1)+K2list(ky,1)+K3list(kz,1)), real(K1list(kx,2)+K2list(ky,2)+K3list(kz,2)), real(K1list(kx,3)+K2list(ky,3)+K3list(kz,3))
                              125 format (f10.6,f10.6,f10.6)
                              write(2,126,advance="no"), aimag(iWnlist2(n))
                              126 format (f)

                              Pol_Old=GWPol_k_iWn_Old(kx,ky,kz,n,:,:)
                              Pol_New= sum(Polarization_k_iWn(kx,ky,kz,n,:,:))
                              Pol= Pol_Old*(1.0d0-GWMixing_Pol) + Pol_New*(GWMixing_Pol)

                              do i1=1,MatDim
                                   do i2=1,MatDim
                                        write(2,127,advance="no"), real(Pol(i1,i2)), aimag(Pol(i1,i2))
                                        127 format (f,f)
                                   enddo
                              enddo
                              write(2,*)
                         endif
                    enddo
               enddo
          enddo
     enddo
endif



if(Dimen==2) then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do n=1,iWnGridNum
                    if (n>Z) then
                         write(2,128,advance="no") , real(K1list(kx,1)+K2list(ky,1)), real(K1list(kx,2)+K2list(ky,2)), real(K1list(kx,3)+K2list(ky,3))
                         128 format (f10.6,f10.6,f10.6)
                         write(2,131,advance="no"), aimag(iWnlist2(n))
                         131 format (f)
            
                         Pol_Old=GWPol_k_iWn_Old(kx,ky,1,n,:,:)
                         Pol_New= sum(Polarization_k_iWn(kx,ky,1,n,:,:))
                         Pol= Pol_Old*(1.0d0-GWMixing_Pol) + Pol_New*(GWMixing_Pol)
            
                         do i1=1,MatDim
                              do i2=1,MatDim
                                   write(2,132,advance="no"), real(Pol(i1,i2)), aimag(Pol(i1,i2))
                                   132 format (f,f)
                              enddo
                         enddo
                         write(2,*)
                    endif
               enddo
          enddo
     enddo
endif



if(Dimen==1) then
     do kx=1,K1GridNum
          do n=1,iWnGridNum
               if (n>Z) then
                    write(2,133,advance="no") , real(K1list(kx,1)), real(K1list(kx,2)), real(K1list(kx,3))
                    133 format (f10.6,f10.6,f10.6)
                    write(2,134,advance="no"), aimag(iWnlist2(n))
                    134 format (f)
       
                    Pol_Old=GWPol_k_iWn_Old(kx,1,1,n,:,:)
                    Pol_New= sum(Polarization_k_iWn(kx,1,1,n,:,:))
                    Pol= Pol_Old*(1.0d0-GWMixing_Pol) + Pol_New*(GWMixing_Pol)
       
                    do i1=1,MatDim
                         do i2=1,MatDim
                              write(2,135,advance="no"), real(Pol(i1,i2)), aimag(Pol(i1,i2))
                              135 format (f,f)
                         enddo
                    enddo
                    write(2,*)
               endif
          enddo
     enddo
endif







close(2)
return

end subroutine



