subroutine parameterR2K(Dimen,OrbN,MatDim,MaxOrbitN,K1GridNum,K2GridNum,K3GridNum,K1list,K2list,K3list,R1GridNum,R2GridNum,R3GridNum,R1list,R2list,R3list,t0,t1,t2,t3,U,U1,U2,U3,OrbitPosi_car,NeighborN,NeighborInfo,tq,Vq)
Implicit none
!arguments
integer, intent(in) :: Dimen,OrbN, MatDim, MaxOrbitN,NeighborN(OrbN,3), K1GridNum, K2GridNum, K3GridNum, R1GridNum, R2GridNum, R3GridNum 
real(8), intent(in) :: K1list(K1GridNum,3),K2list(K2GridNum,3),K3list(K3GridNum,3),R1list(R1GridNum,3),R2list(R2GridNum,3),R3list(R3GridNum,3)
real(8), intent(in) :: t0(MatDim,MatDim), t1(MatDim,MatDim), t2(MatDim,MatDim), t3(MatDim,MatDim), U(MatDim,MatDim), U1(MatDim,MatDim), U2(MatDim,MatDim), U3(MatDim,MatDim)
real(8), intent(in) :: NeighborInfo(OrbN,3,MaxOrbitN,4), OrbitPosi_car(OrbN,3)
complex(8), intent(inout) :: tq(K1GridNum,K2GridNum,K3GridNum,MatDim,MatDim),Vq(K1GridNum,K2GridNum,K3GridNum,MatDim,MatDim)


!local variable
integer k1,k2,k3,o,o2,n,m

tq(:,:,:,:,:)=(0.0d0,0.0d0)
Vq(:,:,:,:,:)=(0.0d0,0.0d0)

if(Dimen==3) then
     do k1=1,K1GridNum
          do k2=1,K2GridNum
               do k3=1,K3GridNum
                    do o=1,OrbN
                         do n=1,3 
                              do m=1,NeighborN(o,n) 
                                   do o2=1,OrbN

                                        if( int(NeighborInfo(o,n,m,4))==o2 ) then
                                             if (n==1) then
                                                  tq(k1,k2,k3,o,o2)=tq(k1,k2,k3,o,o2)+ ( t1(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  tq(k1,k2,k3,OrbN+o,o2)=tq(k1,k2,k3,OrbN+o,o2)+ ( t1(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  tq(k1,k2,k3,o,OrbN+o2)=tq(k1,k2,k3,o,OrbN+o2)+ ( t1(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  tq(k1,k2,k3,OrbN+o,OrbN+o2)=tq(k1,k2,k3,OrbN+o,OrbN+o2)+ ( t1(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )



                                                  Vq(k1,k2,k3,o,o2)=Vq(k1,k2,k3,o,o2)+ ( U1(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  Vq(k1,k2,k3,o,OrbN+o2)=Vq(k1,k2,k3,o,OrbN+o2)+ ( U1(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  Vq(k1,k2,k3,OrbN+o,o2)=Vq(k1,k2,k3,OrbN+o,o2)+ ( U1(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  Vq(k1,k2,k3,OrbN+o,OrbN+o2)=Vq(k1,k2,k3,OrbN+o,OrbN+o2)+ ( U1(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )


                                             endif

                                             if (n==2) then
                                                  tq(k1,k2,k3,o,o2)=tq(k1,k2,k3,o,o2)+ ( t2(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  tq(k1,k2,k3,OrbN+o,o2)=tq(k1,k2,k3,OrbN+o,o2)+ ( t2(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  tq(k1,k2,k3,o,OrbN+o2)=tq(k1,k2,k3,o,OrbN+o2)+ ( t2(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  tq(k1,k2,k3,OrbN+o,OrbN+o2)=tq(k1,k2,k3,OrbN+o,OrbN+o2)+ ( t2(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )           


                                                  Vq(k1,k2,k3,o,o2)=Vq(k1,k2,k3,o,o2)+ ( U2(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  Vq(k1,k2,k3,o,OrbN+o2)=Vq(k1,k2,k3,o,OrbN+o2)+ ( U2(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  Vq(k1,k2,k3,OrbN+o,o2)=Vq(k1,k2,k3,OrbN+o,o2)+ ( U2(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )

                                                  Vq(k1,k2,k3,OrbN+o,OrbN+o2)=Vq(k1,k2,k3,OrbN+o,OrbN+o2)+ ( U2(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )


       
                                             endif

                                             if (n==3) then
                                                  tq(k1,k2,k3,o,o2)=tq(k1,k2,k3,o,o2)+ ( t3(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  tq(k1,k2,k3,OrbN+o,o2)=tq(k1,k2,k3,OrbN+o,o2)+ ( t3(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  tq(k1,k2,k3,o,OrbN+o2)=tq(k1,k2,k3,o,OrbN+o2)+ ( t3(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  tq(k1,k2,k3,OrbN+o,OrbN+o2)=tq(k1,k2,k3,OrbN+o,OrbN+o2)+ ( t3(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )       


                                                  Vq(k1,k2,k3,o,o2)=Vq(k1,k2,k3,o,o2)+ ( U3(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )

                                                  Vq(k1,k2,k3,o,OrbN+o2)=Vq(k1,k2,k3,o,OrbN+o2)+ ( U3(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )
                                                  Vq(k1,k2,k3,OrbN+o,o2)=Vq(k1,k2,k3,OrbN+o,o2)+ ( U3(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )

                                                  Vq(k1,k2,k3,OrbN+o,OrbN+o2)=Vq(k1,k2,k3,OrbN+o,OrbN+o2)+ ( U3(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)+K3list(k3,:)))  ) )

           
                                             endif
                             
                                        endif
                                   enddo
                              enddo 
                         enddo
                    enddo
               enddo
          enddo
     enddo
endif



if(Dimen==2) then
     do k1=1,K1GridNum
          do k2=1,K2GridNum
               do o=1,OrbN
                    do n=1,3 
                         do m=1,NeighborN(o,n) 
                              do o2=1,OrbN
                                   if( int(NeighborInfo(o,n,m,4))==o2 ) then
                                        if (n==1) then
                                             tq(k1,k2,1,o,o2)=tq(k1,k2,1,o,o2)+ ( t1(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             tq(k1,k2,1,OrbN+o,o2)=tq(k1,k2,1,OrbN+o,o2)+ ( t1(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             tq(k1,k2,1,o,OrbN+o2)=tq(k1,k2,1,o,OrbN+o2)+ ( t1(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             tq(k1,k2,1,OrbN+o,OrbN+o2)=tq(k1,k2,1,OrbN+o,OrbN+o2)+ ( t1(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )



                                             Vq(k1,k2,1,o,o2)=Vq(k1,k2,1,o,o2)+ ( U1(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             Vq(k1,k2,1,o,OrbN+o2)=Vq(k1,k2,1,o,OrbN+o2)+ ( U1(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             Vq(k1,k2,1,OrbN+o,o2)=Vq(k1,k2,1,OrbN+o,o2)+ ( U1(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             Vq(k1,k2,1,OrbN+o,OrbN+o2)=Vq(k1,k2,1,OrbN+o,OrbN+o2)+ ( U1(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )


                                        endif

                                        if (n==2) then
                                             tq(k1,k2,1,o,o2)=tq(k1,k2,1,o,o2)+ ( t2(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             tq(k1,k2,1,OrbN+o,o2)=tq(k1,k2,1,OrbN+o,o2)+ ( t2(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             tq(k1,k2,1,o,OrbN+o2)=tq(k1,k2,1,o,OrbN+o2)+ ( t2(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             tq(k1,k2,1,OrbN+o,OrbN+o2)=tq(k1,k2,1,OrbN+o,OrbN+o2)+ ( t2(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )           


                                             Vq(k1,k2,1,o,o2)=Vq(k1,k2,1,o,o2)+ ( U2(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             Vq(k1,k2,1,o,OrbN+o2)=Vq(k1,k2,1,o,OrbN+o2)+ ( U2(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             Vq(k1,k2,1,OrbN+o,o2)=Vq(k1,k2,1,OrbN+o,o2)+ ( U2(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             Vq(k1,k2,1,OrbN+o,OrbN+o2)=Vq(k1,k2,1,OrbN+o,OrbN+o2)+ ( U2(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )


       
                                        endif

                                        if (n==3) then
                                             tq(k1,k2,1,o,o2)=tq(k1,k2,1,o,o2)+ ( t3(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             tq(k1,k2,1,OrbN+o,o2)=tq(k1,k2,1,OrbN+o,o2)+ ( t3(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             tq(k1,k2,1,o,OrbN+o2)=tq(k1,k2,1,o,OrbN+o2)+ ( t3(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             tq(k1,k2,1,OrbN+o,OrbN+o2)=tq(k1,k2,1,OrbN+o,OrbN+o2)+ ( t3(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )       


                                             Vq(k1,k2,1,o,o2)=Vq(k1,k2,1,o,o2)+ ( U3(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )

                                             Vq(k1,k2,1,o,OrbN+o2)=Vq(k1,k2,1,o,OrbN+o2)+ ( U3(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )
                                             Vq(k1,k2,1,OrbN+o,o2)=Vq(k1,k2,1,OrbN+o,o2)+ ( U3(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )

                                             Vq(k1,k2,1,OrbN+o,OrbN+o2)=Vq(k1,k2,1,OrbN+o,OrbN+o2)+ ( U3(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)+K2list(k2,:)))  ) )

           
                                        endif
                             
                                   endif
                              enddo
                         enddo 
                    enddo
               enddo
          enddo
     enddo
endif



if(Dimen==1) then
     do k1=1,K1GridNum
          do o=1,OrbN
               do n=1,3 
                    do m=1,NeighborN(o,n) 
                         do o2=1,OrbN

                              if( int(NeighborInfo(o,n,m,4))==o2 ) then
                                   if (n==1) then
                                        tq(k1,1,1,o,o2)=tq(k1,1,1,o,o2)+ ( t1(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        tq(k1,1,1,OrbN+o,o2)=tq(k1,1,1,OrbN+o,o2)+ ( t1(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        tq(k1,1,1,o,OrbN+o2)=tq(k1,1,1,o,OrbN+o2)+ ( t1(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        tq(k1,1,1,OrbN+o,OrbN+o2)=tq(k1,1,1,OrbN+o,OrbN+o2)+ ( t1(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )



                                        Vq(k1,1,1,o,o2)=Vq(k1,1,1,o,o2)+ ( U1(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        Vq(k1,1,1,o,OrbN+o2)=Vq(k1,1,1,o,OrbN+o2)+ ( U1(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        Vq(k1,1,1,OrbN+o,o2)=Vq(k1,1,1,OrbN+o,o2)+ ( U1(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        Vq(k1,1,1,OrbN+o,OrbN+o2)=Vq(k1,1,1,OrbN+o,OrbN+o2)+ ( U1(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )

                                   endif

                                   if (n==2) then
                                        tq(k1,1,1,o,o2)=tq(k1,1,1,o,o2)+ ( t2(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        tq(k1,1,1,OrbN+o,o2)=tq(k1,1,1,OrbN+o,o2)+ ( t2(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        tq(k1,1,1,o,OrbN+o2)=tq(k1,1,1,o,OrbN+o2)+ ( t2(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        tq(k1,1,1,OrbN+o,OrbN+o2)=tq(k1,1,1,OrbN+o,OrbN+o2)+ ( t2(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )           

                                        Vq(k1,1,1,o,o2)=Vq(k1,1,1,o,o2)+ ( U2(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        Vq(k1,1,1,o,OrbN+o2)=Vq(k1,1,1,o,OrbN+o2)+ ( U2(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        Vq(k1,1,1,OrbN+o,o2)=Vq(k1,1,1,OrbN+o,o2)+ ( U2(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )

                                        Vq(k1,1,1,OrbN+o,OrbN+o2)=Vq(k1,1,1,OrbN+o,OrbN+o2)+ ( U2(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )


       
                                   endif
                                   if (n==3) then
                                        tq(k1,1,1,o,o2)=tq(k1,1,1,o,o2)+ ( t3(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        tq(k1,1,1,OrbN+o,o2)=tq(k1,1,1,OrbN+o,o2)+ ( t3(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        tq(k1,1,1,o,OrbN+o2)=tq(k1,1,1,o,OrbN+o2)+ ( t3(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        tq(k1,1,1,OrbN+o,OrbN+o2)=tq(k1,1,1,OrbN+o,OrbN+o2)+ ( t3(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )       


                                        Vq(k1,1,1,o,o2)=Vq(k1,1,1,o,o2)+ ( U3(o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )

                                        Vq(k1,1,1,o,OrbN+o2)=Vq(k1,1,1,o,OrbN+o2)+ ( U3(o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        Vq(k1,1,1,OrbN+o,o2)=Vq(k1,1,1,OrbN+o,o2)+ ( U3(OrbN+o,o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )
                                        Vq(k1,1,1,OrbN+o,OrbN+o2)=Vq(k1,1,1,OrbN+o,OrbN+o2)+ ( U3(OrbN+o,OrbN+o2) * exp( (0.0d0,1.0d0)*dot_product( OrbitPosi_car(o,:)-NeighborInfo(o,n,m,1:3) ,(K1list(k1,:)))  ) )

           
                                   endif
                            
                              endif
                         enddo
                    enddo 
               enddo
          enddo
     enddo
endif


return 

end subroutine
