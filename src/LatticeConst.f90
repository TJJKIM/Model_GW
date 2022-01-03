



subroutine LatticeConst(Dimen,a1,a2,a3, OrbN, MaxOrbitN, OrbitPosi_car, NeighborInfo, NeighborDist, NeighborN)
!use modmain
Implicit none
!argument
real(8), intent(in) :: a1(3),a2(3),a3(3)
integer, intent(in) :: OrbN, MaxOrbitN,Dimen
real(8), intent(in) :: OrbitPosi_car(OrbN,3)
real(8), intent(inout) :: NeighborInfo(OrbN,3,MaxOrbitN,4), NeighborDist(OrbN,3)
integer, intent(inout) :: NeighborN(OrbN,3)

!local variable
integer i1,i2,i3,lao,o,n,m,che(3)
real(8) Dist, targetorb_pos(3) ,diffcrit, dummy(3)

diffcrit=0.0000000001
NeighborDist(:,:)=100.0
NeighborInfo=0.0d0

if(Dimen==3) then
!loop for finding the distance of neighbors
     do o=1,OrbN  
          do i1=-3,3
               do i2=-3,3
                    do i3=-3,3
                         do lao=1,OrbN  
                              targetorb_pos= (a1(:)*i1 + a2(:)*i2 + a3(:)*i3) + OrbitPosi_car(lao,:)
                              call VecDist(OrbitPosi_car(o,:),targetorb_pos,Dist)
                              dummy=NeighborDist(o,:)

                              if  ( (NeighborDist(o,3)-Dist) > diffcrit .and. Dist >diffcrit) then
                                   if ( abs(Dist-dummy(2))>diffcrit .and. abs(Dist-dummy(1))>diffcrit ) then
                                         dummy(3)=Dist
                                   endif
                              end if
                      
                              if ( (NeighborDist(o,2)-Dist) > diffcrit .and. Dist >diffcrit) then
                                   if ( abs(Dist-dummy(1)) >diffcrit) then
                                        dummy(3)=NeighborDist(o,2)
                                        dummy(2)=Dist
                                   endif
                              end if

                              if ( (NeighborDist(o,1)-Dist) > diffcrit .and. Dist >diffcrit) then
                                   dummy(3)=NeighborDist(o,2)
                                   dummy(2)=NeighborDist(o,1)
                                   dummy(1)=Dist
                              end if

                              NeighborDist(o,:)=dummy
             
                         enddo
                    enddo
                enddo
          enddo
     enddo     

!loop for finding the neighbors whose distance is corresponding to the distance we founded in upper loop
     do o=1,OrbN
          che(:)=1
          do i1=-3,3
               do i2=-3,3
                    do i3=-3,3
                         do lao=1,OrbN
                              targetorb_pos= (a1(:)*i1 + a2(:)*i2 + a3(:)*i3) + OrbitPosi_car(lao,:)
                              call VecDist(OrbitPosi_car(o,:),targetorb_pos,Dist)
                              do n=1,3 
                                   if (abs(Dist-NeighborDist(o,n))<diffcrit) then
                                        NeighborInfo(o,n,che(n),1:3)=targetorb_pos
                                        NeighborInfo(o,n,che(n),4)=lao
                                        che(n)=che(n)+1
                                   endif
                              enddo 
                         enddo
                    enddo
                enddo
          enddo
     enddo

! loop for findind out the number of neighbor
     do o=1,OrbN
          do n=1,3
               che(n)=0
               do m=1,MaxOrbitN
                    if( int(NeighborInfo(o,n,m,4))/=0) then
                         che(n)=che(n)+1
                         NeighborN(o,n)=che(n)
                    endif
               enddo
          enddo
     enddo

endif


if(Dimen==2) then
!loop for finding the distance of neighbors
     do o=1,OrbN  
          do i1=-3,3
               do i2=-3,3
                    do lao=1,OrbN  
                         targetorb_pos= (a1(:)*i1 + a2(:)*i2) + OrbitPosi_car(lao,:)
                         call VecDist(OrbitPosi_car(o,:),targetorb_pos,Dist)
                         dummy=NeighborDist(o,:)

                         if  ( (NeighborDist(o,3)-Dist) > diffcrit .and. Dist >diffcrit) then
                              if ( abs(Dist-dummy(2))>diffcrit .and. abs(Dist-dummy(1))>diffcrit ) then
                                    dummy(3)=Dist
                              endif
                         end if
                      
                         if ( (NeighborDist(o,2)-Dist) > diffcrit .and. Dist >diffcrit) then
                              if ( abs(Dist-dummy(1)) >diffcrit) then
                                   dummy(3)=NeighborDist(o,2)
                                   dummy(2)=Dist
                              endif
                         end if

                         if ( (NeighborDist(o,1)-Dist) > diffcrit .and. Dist >diffcrit) then
                              dummy(3)=NeighborDist(o,2)
                              dummy(2)=NeighborDist(o,1)
                              dummy(1)=Dist
                         end if

                         NeighborDist(o,:)=dummy
                    enddo
               enddo
          enddo
     enddo

!loop for finding the neighbors whose distance is corresponding to the distance we founded in upper loop
     do o=1,OrbN
          che(:)=1
          do i1=-3,3
               do i2=-3,3
                    do lao=1,OrbN
                         targetorb_pos= (a1(:)*i1 + a2(:)*i2 ) + OrbitPosi_car(lao,:)
                         call VecDist(OrbitPosi_car(o,:),targetorb_pos,Dist)
                         do n=1,3 
                              if (abs(Dist-NeighborDist(o,n))<diffcrit) then
                                   NeighborInfo(o,n,che(n),1:3)=targetorb_pos
                                   NeighborInfo(o,n,che(n),4)=lao
                                   che(n)=che(n)+1
                              endif
                         enddo 
                    enddo
               enddo
           enddo
     enddo

! loop for findind out the number of neighbor
     do o=1,OrbN
          do n=1,3
               che(n)=0
               do m=1,MaxOrbitN
                    if( int(NeighborInfo(o,n,m,4))/=0) then
                         che(n)=che(n)+1
                         NeighborN(o,n)=che(n)
                    endif
               enddo
          enddo
     enddo

endif



if(Dimen==1) then
!loop for finding the distance of neighbors
     do o=1,OrbN  
          do i1=-3,3
               do lao=1,OrbN  
                    targetorb_pos= (a1(:)*i1) + OrbitPosi_car(lao,:)
                    call VecDist(OrbitPosi_car(o,:),targetorb_pos,Dist)
                    dummy=NeighborDist(o,:)

                    if  ( (NeighborDist(o,3)-Dist) > diffcrit .and. Dist >diffcrit) then
                         if ( abs(Dist-dummy(2))>diffcrit .and. abs(Dist-dummy(1))>diffcrit ) then
                               dummy(3)=Dist
                         endif
                    end if
                      
                    if ( (NeighborDist(o,2)-Dist) > diffcrit .and. Dist >diffcrit) then
                         if ( abs(Dist-dummy(1)) >diffcrit) then
                              dummy(3)=NeighborDist(o,2)
                              dummy(2)=Dist
                         endif
                    end if

                    if ( (NeighborDist(o,1)-Dist) > diffcrit .and. Dist >diffcrit) then
                         dummy(3)=NeighborDist(o,2)
                         dummy(2)=NeighborDist(o,1)
                         dummy(1)=Dist
                    end if

                    NeighborDist(o,:)=dummy
             
               enddo
          enddo
      enddo

!loop for finding the neighbors whose distance is corresponding to the distance we founded in upper loop
     do o=1,OrbN
          che(:)=1
          do i1=-3,3
               do lao=1,OrbN
                    targetorb_pos= (a1(:)*i1) + OrbitPosi_car(lao,:)
                    call VecDist(OrbitPosi_car(o,:),targetorb_pos,Dist)
                    do n=1,3 
                         if (abs(Dist-NeighborDist(o,n))<diffcrit) then
                              NeighborInfo(o,n,che(n),1:3)=targetorb_pos
                              NeighborInfo(o,n,che(n),4)=lao
                              che(n)=che(n)+1
                         endif
                    enddo 
               enddo
          enddo
      enddo

! loop for findind out the number of neighbor
     do o=1,OrbN
          do n=1,3
               che(n)=0
               do m=1,MaxOrbitN
                    if( int(NeighborInfo(o,n,m,4))/=0) then
                         che(n)=che(n)+1
                         NeighborN(o,n)=che(n)
                    endif
               enddo
          enddo
     enddo

endif





!print *, NeighborDist(1,:)
!print *, NeighborDist(2,:)

!print *, "*----------------------------------------------*"
!print *, "1-1"
!do n=1,MaxOrbitN
!     print *, NeighborInfo(1,1,n,:)
!enddo
!print *, "1-2"
!do n=1,MaxOrbitN
!print *, NeighborInfo(1,2,n,:)
!enddo
!print *, "1-3"
!do n=1,MaxOrbitN
!print *, NeighborInfo(1,3,n,:)
!enddo

!print *, "2-1"
!do n=1,MaxOrbitN
!print *, NeighborInfo(2,1,n,:)
!enddo
!print *, "2-2"
!do n=1,MaxOrbitN
!print *, NeighborInfo(2,2,n,:)
!enddo
!print *, "2-3"
!do n=1,MaxOrbitN
!print *, NeighborInfo(2,3,n,:)
!enddo

!print *, NEighborInfo(1,3,:,4)
!print *, NEighborInfo(2,3,:,4)
!print *, NeighborN(1,:)
!print *, NeighborN(2,:)

return
end subroutine


subroutine ReciprocalV(a1,a2,a3,b1,b2,b3,PI)
Implicit none
!argument
real(8), intent (in) :: PI
real(8),intent (in) :: a1(3),a2(3),a3(3)
real(8), intent (inout) :: b1(3), b2(3), b3(3)

!local variable
real(8) cross1(3),cross2(3),cross3(3)
real(8) fact(3)

cross1(1) = a2(2)*a3(3)-a2(3)*a3(2)
cross1(2) = -(a2(1)*a3(3)-a2(3)*a3(1))
cross1(3) = a2(1)*a3(2)-a2(2)*a3(1)

cross2(1) = a3(2)*a1(3)-a3(3)*a1(2)
cross2(2) = -(a3(1)*a1(3)-a3(3)*a1(1))
cross2(3) = a3(1)*a1(2)-a3(2)*a1(1)

cross3(1) = a1(2)*a2(3)-a1(3)*a2(2)
cross3(2) = -(a1(1)*a2(3)-a1(3)*a2(1))
cross3(3) = a1(1)*a2(2)-a1(2)*a2(1)

fact(1)=2*PI/(a1(1)*cross1(1) + a1(2)*cross1(2) + a1(3)*cross1(3))
fact(2)=2*PI/(a2(1)*cross2(1) + a2(2)*cross2(2) + a2(3)*cross2(3))
fact(3)=2*PI/(a3(1)*cross3(1) + a3(2)*cross3(2) + a3(3)*cross3(3))

b1(:)=fact(1)*cross1(:)
b2(:)=fact(2)*cross2(:)
b3(:)=fact(3)*cross3(:)
return
end subroutine



subroutine VecDist(v1,v2,val)
Implicit none
!argument
real(8), intent(in) :: v1(3), v2(3)
real(8), intent(inout) :: val

!local variable
real(8) diffvec(3)

diffvec=v1-v2
val=sqrt( diffvec(1)*diffvec(1) + diffvec(2)*diffvec(2) + diffvec(3)*diffvec(3) )

return

end subroutine 







subroutine FindSameAtomOrb(OrbitPosi_car,OrbN,distingPos,distingAtom)

Implicit none
!argument
integer,intent(in) :: OrbN
real(8), intent(in) :: OrbitPosi_car(OrbN,3)
integer, intent(inout) :: distingAtom(OrbN,OrbN)
real(8), intent(inout) :: distingPos(OrbN,3)

!local variable
integer i1,i2,i,j(OrbN),m 
real crit

i=1
j(1)=1
crit=0.0000001

distingPos(1,:)=OrbitPosi_car(1,:)
distingAtom(1,1)=1
j(1)=j(1)+1

do i1=2,OrbN
          m=1
          do i2=1,i
               if ( sum(abs(OrbitPosi_car(i1,:)-distingPos(i2,:))) <crit ) then
                    distingAtom(i2,j(i2))=i1
                    j(i2)=j(i2)+1
                    m=0
               endif
          enddo
          if (m==1) then 
               i=i+1
               distingPos(i,:)=OrbitPosi_car(i1,:)
               distingAtom(i,1)=i1
               j(i)=2
          endif
 
enddo

return
endsubroutine
