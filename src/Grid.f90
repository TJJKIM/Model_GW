subroutine Kgrid(b,KGridNum,Klist)
Implicit none
!argument
real(8), intent(in) :: b(3)
integer, intent(in) :: KGridNum
real(8), intent(inout) :: Klist(KGridNum,3)
! local variable
integer k

do k=1,(KGridNum)
     Klist(k,:)=b*(dble(k-1)/dble(KGridNum))
enddo
return
end subroutine



subroutine Rgrid(a,RGridNum,Rlist)
Implicit none
!argument
real(8), intent(in) :: a(3)
integer, intent(in) :: RGridNum
real(8), intent(inout) :: Rlist(RGridNum,3)
! local variable
integer i

do i=1,(RGridNum)
     Rlist(i,:)=a*dble(i-1)
enddo
return
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for making grid from 0 to given number !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gridfromZ4(beta,NumGrid,Taulist)
Implicit none

!arguments
integer, intent(in) :: NumGrid
real(8), intent(in) :: beta
real(8), intent(inout) :: Taulist(NumGrid)

! local variable 
integer i
real (8) r

r= ((8.d0)*beta)**(1.d0/4.d0)/dble(NumGrid-1)

do i =1,(NumGrid)
     if ( i < (NumGrid/2) .or. i == (NumGrid/2) )   then
          Taulist(i)= ( r*dble(i-1) )**(4.d0)
     endif
     if  ( i > (NumGrid/2)  ) then
          Taulist(i)= -( r*dble(i-NumGrid) )**(4.d0) + beta
     endif
end do



return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for making grid from 0 to given number !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gridfromZEven(beta,NumGrid,Taulist)
Implicit none

!arguments
integer, intent(in) :: NumGrid
real(8), intent(in) :: beta
real(8), intent(inout) :: Taulist(NumGrid)

! local variable 
integer i

do i =1,(NumGrid+1)
     Taulist(i)= (dble(i-1)*(beta))/dble(NumGrid)
end do

return
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for obtaining matsubara frequency grid !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MatsubfGrid(PI,NumGrid,iWnlist,beta)
Implicit none
!arguments
integer,intent(in) :: NumGrid
real(8), intent(in) ::  beta, PI
complex(8),intent(inout) :: iWnlist(NumGrid)

!local variables
integer i
complex(8) dummy

!print *,(0.0, (2*(-(NumGrid/2)+20)+1)*PI/beta)
do i = 1,(NumGrid)
     dummy=((2.0d0)*(-(NumGrid/2)+i-1)+1)*PI/beta*cmplx(0.0d0,1.0d0)
     iWnlist(i)= dummy
end do

return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine for obtaining matsubara frequency grid !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MatsubfGridEven(PI,NumGrid,iWnlist2,beta)
Implicit none
!arguments
integer,intent(in) :: NumGrid
real(8), intent(in) ::  beta, PI
complex(8),intent(inout) :: iWnlist2(NumGrid+1)

!local variables
integer i
complex(8) dummy

!print *,(0.0, (2*(-(NumGrid/2)+20)+1)*PI/beta)
do i = 1,(NumGrid+1)
     dummy=((2.0d0)*(-((NumGrid)/2)+i-1))*PI/beta*(0.0d0,1.0d0)
     iWnlist2(i)= dummy
end do

return
end subroutine

