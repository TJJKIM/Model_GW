!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for reading GWDMFT loop Num !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine ReadGWDMFTloopNum(GWDMFTloopNum)
Implicit none
!arguments
integer, intent(inout) :: GWDMFTloopNum


open (2, file ='GWDMFT.loopNum', status = 'old')
read(2,*) GWDMFTloopNum

close(2)

return
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for reading Impurity self energy from "SelfE_k_iWn.dat" !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine  ReadSelfE_k_iWn(Dimen,OrbN,MatDim,iWnGridNum,K1GridNum,K2GridNum,K3GridNum,GWSelfE_k_iWn)
Implicit none
! arguments
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum, iWnGridNum, Dimen, MatDim, OrbN
complex(8) ,intent(inout) :: GWSelfE_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum,MatDim,MatDim)

! local variable
integer Z,kx,ky,kz,n,o,o2,i
real kxval, kyval,kzval, iWnval 
real(8) SelfE(MatDim*MatDim*2)

Z=(iWnGridNum/2)

open (2, file ='SelfE_k_iWn.dat', status = 'old')
read(2,*)

if (Dimen==3)then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do kz=1,K3GridNum
                    do n=Z+1,(iWnGridNum)
                         read(2,*) kxval,kyval,kzval,iWnval,SelfE
                         !3 format(f,f,f,f,30(f))
                         !read(2,4) SelfE
                         !4 format (30(f,f))
                         i=1
                         do o=1,MatDim
                              do o2=1,MatDim
                                   GWSelfE_k_iWn(kx,ky,kz,n,o,o2)=SelfE(i)*(1.0d0,0.0d0)+SelfE(i+1)*(0.0d0,1.0d0)
                                   GWSelfE_k_iWn(kx,ky,kz,iWnGridNum+1-n,o,o2)=SelfE(i)*(1.0d0,0.0d0)-SelfE(i+1)*(0.0d0,1.0d0)
                                   i=i+2
                              enddo
                         enddo 
                                   !print *,kxval
                                   !print 3, SelfE(1)*1.0d0,SelfE(2)
                                   !3 format (f,f)
                                   !print *, SelfE(3),",",SelfE(4)
                                   !print *, SelfE(5),",",SelfE(6)
                          !         print *, GwSelfE_k_iWn(kx,ky,kz,n,1,1)
                          !         print *, GwSelfE_k_iWn(kx,ky,kz,n,1,2)
                          !         print *, GwSelfE_k_iWn(kx,ky,kz,n,1,3)
                          !         print *, GwSelfE_k_iWn(kx,ky,kz,n,3,1)
                              
                         !stop
                         !GWSelfE_k_iWn(kx,ky,kz,n)=(1.0d0,0.0d0)*SR+(0.0d0,1.0d0)*SI
                    enddo
               enddo
          enddo
     enddo
endif

if (Dimen==2)then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do n=Z+1,(iWnGridNum)
                    read(2,*) kxval,kyval,kzval,iWnval,SelfE
                    i=1
                    do o=1,MatDim
                         do o2=1,MatDim
                              GWSelfE_k_iWn(kx,ky,1,n,o,o2)=SelfE(i)*(1.0d0,0.0d0)+SelfE(i+1)*(0.0d0,1.0d0)
                              GWSelfE_k_iWn(kx,ky,1,iWnGridNum+1-n,o,o2)=SelfE(i)*(1.0d0,0.0d0)-SelfE(i+1)*(0.0d0,1.0d0)
                              i=i+2
                         enddo
                    enddo
               enddo
          enddo
     enddo
endif

if (Dimen==1)then
     do kx=1,K1GridNum
         do n=Z+1,(iWnGridNum)
              read(2,*) kxval,kyval,kzval,iWnval,SelfE
              i=1
              do o=1,MatDim
                   do o2=1,MatDim
                        GWSelfE_k_iWn(kx,1,1,n,o,o2)=SelfE(i)*(1.0d0,0.0d0)+SelfE(i+1)*(0.0d0,1.0d0)
                        GWSelfE_k_iWn(kx,1,1,iWnGridNum+1-n,o,o2)=SelfE(i)*(1.0d0,0.0d0)-SelfE(i+1)*(0.0d0,1.0d0)
                        i=i+2
                   enddo
              enddo
         enddo
     enddo
endif

close(2)
return
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for reading Impurity self energy from "SelfE_k_iWn.dat" !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine  ReadPol_k_iWn(Dimen,OrbN,MatDim,iWnGridNum,K1GridNum,K2GridNum,K3GridNum,GWPol_k_iWn)
Implicit none
! arguments
integer, intent(in) :: K1GridNum,K2GridNum,K3GridNum, iWnGridNum, Dimen, MatDim, OrbN
complex(8) ,intent(inout) :: GWPol_k_iWn(K1GridNum,K2GridNum,K3GridNum,iWnGridNum+1,MatDim,MatDim)

! local variable
integer Z,kx,ky,kz,n,o,o2,i
real kxval, kyval,kzval, iWnval
real(8) Pol(MatDim*MatDim*2)

Z=(iWnGridNum/2)

open (2, file ='Pol_k_iWn.dat', status = 'old')
read(2,*)

if (Dimen==3)then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
               do kz=1,K3GridNum
                    read(2,*) kxval,kyval,kzval,iWnval,Pol
                    i=1
                    do o=1,MatDim
                         do o2=1,MatDim
                              GWPol_k_iWn(kx,ky,kz,Z+1,o,o2)=Pol(i)*(1.0d0,0.0d0)+Pol(i+1)*(0.0d0,1.0d0)
                              i=i+2
                         enddo
                    enddo

                    do n=Z+2,(iWnGridNum+1)
                         read(2,*) kxval,kyval,kzval,iWnval,Pol
                         !3 format(f,f,f,f,30(f))
                         !read(2,4) SelfE
                         !4 format (30(f,f))
                         i=1
                         do o=1,MatDim
                              do o2=1,MatDim
                                   GWPol_k_iWn(kx,ky,kz,n,o,o2)=Pol(i)*(1.0d0,0.0d0)+Pol(i+1)*(0.0d0,1.0d0)
                                   GWPol_k_iWn(kx,ky,kz,iWnGridNum+2-n,o,o2)=Pol(i)*(1.0d0,0.0d0)-Pol(i+1)*(0.0d0,1.0d0)
                                   i=i+2
                              enddo
                         enddo
                                   !print *,kxval
                                   !print 3, SelfE(1)*1.0d0,SelfE(2)
                                   !3 format (f,f)
                                   !print *, SelfE(3),",",SelfE(4)
                                   !print *, SelfE(5),",",SelfE(6)
                          !         print *, GwSelfE_k_iWn(kx,ky,kz,n,1,1)
                          !         print *, GwSelfE_k_iWn(kx,ky,kz,n,1,2)
                          !         print *, GwSelfE_k_iWn(kx,ky,kz,n,1,3)
                          !         print *, GwSelfE_k_iWn(kx,ky,kz,n,3,1)

                         !stop
                         !GWSelfE_k_iWn(kx,ky,kz,n)=(1.0d0,0.0d0)*SR+(0.0d0,1.0d0)*SI
                    enddo
               enddo
          enddo
     enddo
endif


if (Dimen==2)then
     do kx=1,K1GridNum
          do ky=1,K2GridNum
              read(2,*) kxval,kyval,kzval,iWnval,Pol
              i=1
              do o=1,MatDim
                   do o2=1,MatDim
                        GWPol_k_iWn(kx,ky,1,Z+1,o,o2)=Pol(i)*(1.0d0,0.0d0)+Pol(i+1)*(0.0d0,1.0d0)
                        i=i+2
                   enddo
              enddo

              do n=Z+2,(iWnGridNum+1)
                   read(2,*) kxval,kyval,kzval,iWnval,Pol
                   i=1
                   do o=1,MatDim
                        do o2=1,MatDim
                             GWPol_k_iWn(kx,ky,1,n,o,o2)=Pol(i)*(1.0d0,0.0d0)+Pol(i+1)*(0.0d0,1.0d0)
                             GWPol_k_iWn(kx,ky,1,iWnGridNum+2-n,o,o2)=Pol(i)*(1.0d0,0.0d0)-Pol(i+1)*(0.0d0,1.0d0)
                             i=i+2
                        enddo
                   enddo
              enddo
          enddo
     enddo
endif


if (Dimen==1)then
     do kx=1,K1GridNum
         read(2,*) kxval,kyval,kzval,iWnval,Pol
         i=1
         do o=1,MatDim
              do o2=1,MatDim
                   GWPol_k_iWn(kx,1,1,Z+1,o,o2)=Pol(i)*(1.0d0,0.0d0)+Pol(i+1)*(0.0d0,1.0d0)
                   i=i+2
              enddo
         enddo

         do n=Z+2,(iWnGridNum+1)
              read(2,*) kxval,kyval,kzval,iWnval,Pol
              i=1
              do o=1,MatDim
                   do o2=1,MatDim
                        GWPol_k_iWn(kx,1,1,n,o,o2)=Pol(i)*(1.0d0,0.0d0)+Pol(i+1)*(0.0d0,1.0d0)
                        GWPol_k_iWn(kx,1,1,iWnGridNum+2-n,o,o2)=Pol(i)*(1.0d0,0.0d0)-Pol(i+1)*(0.0d0,1.0d0)
                        i=i+2
                   enddo
              enddo
         enddo
     enddo
endif





close(2)
return
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for reading Impurity self energy from "ImpuritySelfE.dat" !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine  ReadImpuritySelfE(MatDim,iWnGridNum,ImpuritySelfE)
Implicit none
! arguments
integer, intent(in) :: iWnGridNum, MatDim
complex(8) ,intent(inout) :: ImpuritySelfE(iWnGridNum,MatDim)

! local variable
integer n,i1
real(8) iWn,SelfE(MatDim*2)

open (2, file ='ImpuritySelfE.dat', status = 'old')

do n=1,(iWnGridNum/2)
     read(2,*) iWn, SelfE
     do i1=1,MatDim
          ImpuritySelfE(iWnGridNum/2-n+1,i1)=(1.0d0,0.0d0)*SelfE(2*(i1-1)+1)-(0.0d0,1.0d0)*SelfE(2*(i1-1)+2)
          ImpuritySelfE(iWnGridNum/2+n,i1)=(1.0d0,0.0d0)*SelfE(2*(i1-1)+1)+(0.0d0,1.0d0)*SelfE(2*(i1-1)+2)

     enddo
enddo



close(2)
return
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for reading Impurity self energy from "ImpurityPol.dat" !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine ReadImpurityPol(MatDim,iWnGridNum,ImpurityPol)

! arguments
integer, intent(in) :: iWnGridNum, MatDim
complex(8) ,intent(inout) :: ImpurityPol(iWnGridNum+1,MatDim)

! local variable
integer n,i1
real iWn, Pol(MatDim*2)

open (2, file ='ImpurityPol.dat', status = 'old')


read(2,*) iWn, Pol
do i1=1,MatDim
     ImpurityPol(iWnGridNum/2+1,i1)=(1.0d0,0.0d0)*Pol(2*(i1-1)+1)+(0.0d0,1.0d0)*Pol(2*(i1-1)+2)
enddo
do n=1,(iWnGridNum/2-2)
     read(2,*) iWn, Pol
     do i1=1,MatDim
          ImpurityPol(iWnGridNum/2-n+1,i1)=(1.0d0,0.0d0)*Pol(2*(i1-1)+1)-(0.0d0,1.0d0)*Pol(2*(i1-1)+2)
          ImpurityPol(iWnGridNum/2+n+1,i1)=(1.0d0,0.0d0)*Pol(2*(i1-1)+1)+(0.0d0,1.0d0)*Pol(2*(i1-1)+2)
     enddo
enddo



close(2)



return
end subroutine

