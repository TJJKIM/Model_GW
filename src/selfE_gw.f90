
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for calculating GW-self energy in (R,tau) space !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine GWselfenergy(Dimen,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim,beta,K1list,K2list,K3list,Taulist,R1list,R2list,R3list,dummy_R_tau,dummy_R_tauII,dummy_R_tauIII)
Implicit none


! argument
integer, intent(in) :: Dimen,K1GridNum,K2GridNum,K3GridNum,TauGridNum,R1GridNum,R2GridNum,R3GridNum,MatDim
real(8), intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3), K3list(K3GridNum,3),Taulist(TauGridNum), R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3),beta
complex(8), intent(in) :: dummy_R_tau(R1GridNum,R2GridNum,R3GridNum,TauGridNum,MatDim,MatDim), dummy_R_tauII(R1GridNum,R2GridNum,R3GridNum,TauGridNum,MatDim,MatDim)
complex(8), intent(inout) :: dummy_R_tauIII(R1GridNum,R2GridNum,R3GridNum,TauGridNum,MatDim,MatDim)

! local variable
integer ix,iy,iz,j,k,l,E,n
complex(8) Iden(MatDim,MatDim),dummyH(MatDim,MatDim),dummyG(MatDim,MatDim), dummy(MatDim,MatDim),dummy2(MatDim,MatDim),dummy3(MatDim,MatDim)


if (Dimen==3) then
     do ix = 1, (R1GridNum)
          do iy = 1, (R2GridNum)
               do iz = 1, (R3GridNum)
                    do l = 1, (TauGridNum)
                         dummy_R_tauIII(ix,iy,iz,l,:,:)=dummy_R_tau(ix,iy,iz,l,:,:)*(-1.0d0)*dummy_R_tauII(ix,iy,iz,l,:,:)
                         !dummy_R_tauIII(ix,iy,iz,l,1,2)=dummy_R_tau(ix,iy,iz,l,1,2)*(-1.0d0)*dummy_R_tauII(ix,iy,iz,l,1,2)
                         !dummy_R_tauIII(ix,iy,iz,l,2,1)=dummy_R_tau(ix,iy,iz,l,2,1)*(-1.0d0)*dummy_R_tauII(ix,iy,iz,l,2,1)
                         !dummy_R_tauIII(ix,iy,iz,l,2,2)=dummy_R_tau(ix,iy,iz,l,2,2)*(-1.0d0)*dummy_R_tauII(ix,iy,iz,l,2,2)
                    enddo
               enddo
          enddo
     enddo
endif

if (Dimen==2) then
     do ix = 1, (R1GridNum)
          do iy = 1, (R2GridNum)
               do l = 1, (TauGridNum)
                   dummy_R_tauIII(ix,iy,1,l,:,:)=dummy_R_tau(ix,iy,1,l,:,:)*(-1.0d0)*dummy_R_tauII(ix,iy,1,l,:,:)
                   !dummy_R_tauIII(ix,iy,iz,l,1,2)=dummy_R_tau(ix,iy,iz,l,1,2)*(-1.0d0)*dummy_R_tauII(ix,iy,iz,l,1,2)
                   !dummy_R_tauIII(ix,iy,iz,l,2,1)=dummy_R_tau(ix,iy,iz,l,2,1)*(-1.0d0)*dummy_R_tauII(ix,iy,iz,l,2,1)
                   !dummy_R_tauIII(ix,iy,iz,l,2,2)=dummy_R_tau(ix,iy,iz,l,2,2)*(-1.0d0)*dummy_R_tauII(ix,iy,iz,l,2,2)
               enddo
          enddo
     enddo
endif


if (Dimen==1) then
     do ix = 1, (R1GridNum)
          do l = 1, (TauGridNum)
              dummy_R_tauIII(ix,1,1,l,:,:)=dummy_R_tau(ix,1,1,l,:,:)*(-1.0d0)*dummy_R_tauII(ix,1,1,l,:,:)
              !dummy_R_tauIII(ix,iy,iz,l,1,2)=dummy_R_tau(ix,iy,iz,l,1,2)*(-1.0d0)*dummy_R_tauII(ix,iy,iz,l,1,2)
              !dummy_R_tauIII(ix,iy,iz,l,2,1)=dummy_R_tau(ix,iy,iz,l,2,1)*(-1.0d0)*dummy_R_tauII(ix,iy,iz,l,2,1)
              !dummy_R_tauIII(ix,iy,iz,l,2,2)=dummy_R_tau(ix,iy,iz,l,2,2)*(-1.0d0)*dummy_R_tauII(ix,iy,iz,l,2,2)
          enddo
     enddo
endif



return
end subroutine


