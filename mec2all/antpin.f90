Subroutine ANTPIN (ANBTP,ANGS,ANGS2,PTTP,MOMTEN,eigen)
implicit none

!	Calculates other representations of fault planes with
!		trend and plunge of A and N as input.  All
!		angles are in radians.
!	    Added moment tensor output.

real N(3), MOMTEN(6)
real PTTP(4),ANGS(3),ANGS2(3),ANBTP(6),P(3),T(3),A(3),B(3)
real pi,sr2

real aa(3,3),U(3,3),w(3),isotrop,eigen(3)

integer j
pi = 3.141592653

	DATA SR2/0.707107/
	call TRPL2V(ANBTP(1),A)
	call TRPL2V(ANBTP(3),N)
 do 100 J=1,3
	  T(J) = SR2*(A(J) + N(J))
	  P(J) = SR2*(A(J) - N(J))
100	CONTINUE
      B(1) = P(2)*T(3) - P(3)*T(2)
      B(2) = P(3)*T(1) - P(1)*T(3)
      B(3) = P(1)*T(2) - P(2)*T(1)
      call V2TRPL(P,PTTP(1),PI)
      call V2TRPL(T,PTTP(3),PI)
      call V2TRPL(B,ANBTP(5),PI)
      call AN2DSR(A,N,ANGS,PI)
      call AN2DSR(N,A,ANGS2,PI)
      call AN2MOM(A,N,MOMTEN)


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
end subroutine ANTPIN
