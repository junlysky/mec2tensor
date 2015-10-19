Subroutine DSRIN (ANGS,ANBTP,ANGS2,PTTP,MOMTEN,eigen)
  implicit none

!	Calculates other representations of fault planes with
!		dip, strike and rake (A&R convinention) input.  All
!		angles are in radians.
!       Added moment tensor output (D&W convention)
!       normalized to unit scalar moment
!	    When ported to the PC, on one compiler there
!		was a roundoff problem for limiting cases.  Included
!		a fix.

real N(3), MOMTEN(6)
real PTTP(4),ANGS(3),ANGS2(3),ANBTP(6),P(3),T(3),A(3),B(3)

real sr2,rake,str,dip
real eigen(3),aa(3,3), U(3,3),w(3),isotrop
real, parameter :: pi = 3.1415926536
integer j
DATA SR2/0.707107/

RAKE= ANGS(3)
STR = ANGS(2)
DIP = ANGS(1)
A(1) = COS(RAKE)*COS(STR) + SIN(RAKE)*COS(DIP)*SIN(STR)
A(2) = COS(RAKE)*SIN(STR) - SIN(RAKE)*COS(DIP)*COS(STR)
A(3) = -SIN(RAKE)*SIN(DIP)
N(1) = -SIN(STR)*SIN(DIP)
N(2) = COS(STR)*SIN(DIP)
N(3) = -COS(DIP)
	do j=1,3
          if (abs(A(j)) .le. 0.0001) A(j) = 0.0
          if ((ABS(A(j))-1.0) .gt. 0.0) A(j)=A(j)/abs(A(j))
          if (abs(N(j)) .le. 0.0001) N(j) = 0.0
          if ((ABS(N(j))-1.0) .gt. 0.0) N(j)=N(j)/abs(N(j))
        end do
	call V2TRPL(A,ANBTP(1),PI)
	call V2TRPL(N,ANBTP(3),PI)
	do 100 J=1,3
	T(J) = SR2*(A(J) + N(J))
100	P(J) = SR2*(A(J) - N(J))
	B(1) = P(2)*T(3) - P(3)*T(2)
	B(2) = P(3)*T(1) - P(1)*T(3)
	B(3) = P(1)*T(2) - P(2)*T(1)
	CALL V2TRPL(P,PTTP(1),PI)
	CALL V2TRPL(T,PTTP(3),PI)
	CALL V2TRPL(B,ANBTP(5),PI)
	CALL AN2DSR(N,A,ANGS2,PI)
    CALL AN2MOM(A,N,MOMTEN) 

    AA(3,3) = MOMTEN(1)
    AA(1,1) = MOMTEN(2)
    AA(2,2) = MOMTEN(3)
    AA(1,3) = MOMTEN(4)
    AA(3,1) = MOMTEN(4)
    AA(2,3) = -MOMTEN(5)
    AA(3,2) = -MOMTEN(5) 
    AA(1,2) = -MOMTEN(6)
    AA(2,1) = -MOMTEN(6)
    call eig(aa,u,w)

    isotrop = 0.0      ! M_all = M_isotrop + M_deviatoric
    do j=1,3
       isotrop = isotrop + w(j)
    end do
  
    do j=1,3
       w(j) = w(j) - isotrop/3.0
       eigen(j) = w(j)
    end do

	return
end subroutine DSRIN
