Subroutine AN2MOM(A,N,MOMTEN)
implicit none
!	Starting with the A and N axis, calculates the elements
!	of the moment tensor with unit scalar moment.
!	Coordinate system:  X = North, Y = East, Z = Down
!	Convention used is that of Dziewonski & Woodhouse
!	(JGR 88, 3247-3271, 1983) and Aki & Richards (p 118)
!	If an element is < 0.000001 (ABS), set to zero
!	Moment tensor components:  M(I,j) = A(I)*N(J)+A(J)*N(I)

integer j
real A(3), N(3), MOMTEN(6)
!Moment tensor components:  M(I,j) = A(I)*N(J)+A(J)*N(I)
MOMTEN(1) = 2.0*A(3)*N(3)     !  MRR = M(3,3)
MOMTEN(2) = 2.0*A(1)*N(1)     !  MTT = M(1,1)
MOMTEN(3) = 2.0*A(2)*N(2)     !  MPP = M(2,2)
MOMTEN(4) = A(1)*N(3)+A(3)*N(1)     !  MRT = M(1,3)
MOMTEN(5) = -A(2)*N(3)-A(3)*N(2)!  MRP = -M(2,3)
MOMTEN(6) = -A(2)*N(1)-A(1)*N(2)!  MTP = -M(2,1)   
do j=1,6
   if (ABS(MOMTEN(j)) .LT. 0.000001) MOMTEN(j) = 0.0
end do
return
end subroutine AN2MOM	
