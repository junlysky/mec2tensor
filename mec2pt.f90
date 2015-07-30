! junlysky.info@gmail.com
program main
implicit none

Real,Parameter :: PI = 3.1415926
Real lon,lat,ANGS(3)
Real ANBTP(6),ANGS2(3),PTTP(4),MOMTEN(6)

Character(10) line(5)
Integer total,i,j,stat


open(521,file='mecinfo.dat')
open(520,file='mec.d',status='old',iostat=stat)
if(stat/=0) then
   write(*,*) "open mec.dat fail ."
   stop
   endif
write(521,"(18A8)")"lon","lat","str1","dip1","rake1","str2","dip2","rake2","mrr","mtt","mpp","mrt","mrp","mtp",&
&"p-trend","p-plu","t-trend","t-plu"
Call Getline(520,total) 
read(520,*) line(1:5)
write(*,*)line    
do i=1,total-1
   read(520,*)lon,lat,ANGS(1:3)
   write(*,81)lon,lat,ANGS(1:3)
   Call mec2pt(ANGS,ANBTP,ANGS2,PTTP,MOMTEN)
   write(*,81)lon,lat,ANGS,ANGS2(2),ANGS2(1),ANGS2(3),MOMTEN,PTTP
   write(521,81)lon,lat,ANGS,ANGS2(2),ANGS2(1),ANGS2(3),MOMTEN,PTTP 
81 format(18F8.2)
end do

close(520)
close(521)
end program main

! subroutine 1
subroutine mec2pt(ANGS,ANBTP,ANGS2,PTTP,MOMTEN)
implicit none

Real,Parameter :: PI = 3.1415926
Real N(3), MOMTEN(6)
Real PTTP(4),ANGS(3),ANGS2(3),ANBTP(6),P(3),T(3),A(3),B(3)
Real STR,DIP,RAKE,SR2,RDEG
Integer I,J
DATA SR2/0.707107/
RDEG = 45.0/ATAN(1.0)

DO 300 J=1,3
300  ANGS(J) = ANGS(J)/RDEG

STR = ANGS(1)
DIP = ANGS(2)
RAKE= ANGS(3)
A(1) = COS(RAKE)*COS(STR) + SIN(RAKE)*COS(DIP)*SIN(STR)
A(2) = COS(RAKE)*SIN(STR) - SIN(RAKE)*COS(DIP)*COS(STR)
A(3) = -SIN(RAKE)*SIN(DIP)
N(1) = -SIN(STR)*SIN(DIP)
N(2) = COS(STR)*SIN(DIP)
N(3) = -COS(DIP)
do j=1,3
   if (abs(A(j)) .le. 0.0001) A(j) = 0.0
   IF ((ABS(A(j))-1.0) .gt. 0.0) A(j)=A(j)/abs(A(j))
   if (abs(N(j)) .le. 0.0001) N(j) = 0.0
   IF ((ABS(N(j))-1.0) .gt. 0.0) N(j)=N(j)/abs(N(j))
   end do

CALL V2TRPL(A,ANBTP(1),PI)
CALL V2TRPL(N,ANBTP(3),PI)
DO 101 J=1,3
    T(J) = SR2*(A(J) + N(J))
101 P(J) = SR2*(A(J) - N(J))
    
B(1) = P(2)*T(3) - P(3)*T(2)
B(2) = P(3)*T(1) - P(1)*T(3)
B(3) = P(1)*T(2) - P(2)*T(1)
      
CALL V2TRPL(P,PTTP(1),PI)
CALL V2TRPL(T,PTTP(3),PI)
CALL V2TRPL(B,ANBTP(5),PI)
CALL AN2DSR(N,A,ANGS2,PI)
CALL AN2MOM(A,N,MOMTEN)

DO 400 I=1,3
   ANGS(I) = ANGS(I)*RDEG
   ANGS2(I) = ANGS2(I)*RDEG
   PTTP(I) = PTTP(I)*RDEG
   ANBTP(I) = ANBTP(I)*RDEG
400   CONTINUE
   ANBTP(4) = ANBTP(4)*RDEG
   ANBTP(5) = ANBTP(5)*RDEG
   ANBTP(6) = ANBTP(6)*RDEG
   PTTP(4) = PTTP(4)*RDEG
   Return
end subroutine mec2pt

! subroutine 2
subroutine V2TRPL(XYZ,TRPL,PI)
implicit none
      REAL XYZ(3),TRPL(2)
      REAL PI,C,S
      INTEGER J
      do j=1,3
        if (abs(xyz(j)) .le. 0.0001) xyz(j) = 0.0
        IF (ABS(ABS(XYZ(j))-1.0).LT.0.0001) xyz(j)=xyz(j)/abs(xyz(j))
      end do
      IF (ABS(XYZ(3)) .eq. 1.0) THEN 

!     plunge is 90 degrees

        if (xyz(3) .lt. 0.0) then
          trpl(1) = PI
        else
          TRPL(1) = 0.0
        end if
        TRPL(2) = 0.5*PI
        RETURN
      END IF
      IF (ABS(XYZ(1)) .LT. 0.0001) THEN
        IF (XYZ(2) .GT. 0.0) THEN
          TRPL(1) = PI/2.
        ELSE IF (XYZ(2) .LT. 0.0) THEN
          TRPL(1) = 3.0*PI/2.0
        ELSE
          TRPL(1) = 0.0
        END IF
      ELSE
        TRPL(1) = ATAN2(XYZ(2),XYZ(1))
      END IF
      C = COS(TRPL(1))
      S = SIN(TRPL(1))
      IF (ABS(C) .GE. 0.1) TRPL(2) = ATAN2(XYZ(3),XYZ(1)/C)
      IF (ABS(C) .LT. 0.1) TRPL(2) = ATAN2(XYZ(3),XYZ(2)/S)
      IF (TRPL(2) .LT. 0.0) THEN
        TRPL(2) = -TRPL(2)
        TRPL(1) = TRPL(1) - PI
        END IF
      IF (TRPL(1) .LT. 0.0) TRPL(1) = TRPL(1) + 2.0*PI
      RETURN
END subroutine V2TRPL

! subroutine 3
subroutine AN2DSR(A,N,ANGS,PI)
implicit none
      REAL N(3),A(3),ANGS(3)
      REAL PI,A1,a2,acosarg 
      if (N(3) .eq. -1.0) then
        angs(2) = atan2(a(2),a(1))
        angs(1) = 0.0
      else
        ANGS(2) = ATAN2(-N(1),N(2))
        if (N(3) .eq. 0.0) then
          angs(1) = 0.5*PI
        else IF (ABS(SIN(ANGS(2))) .ge. 0.1) then
          ANGS(1) = ATAN2(-N(1)/SIN(ANGS(2)),-N(3))
        else
          ANGS(1) = ATAN2(N(2)/COS(ANGS(2)),-N(3))
        end if
      end if
      A1 = A(1)*COS(ANGS(2)) + A(2)*SIN(ANGS(2))
      if (abs(a1) .lt. 0.0001) a1 = 0.0
      if (a(3) .ne. 0.0) then
        if (angs(1) .ne. 0.0) then
          ANGS(3) = ATAN2(-A(3)/SIN(ANGS(1)),A1)
        else
          ANGS(3) = atan2(-1000000.0*A(3),A1)
        end if
      else
        a2 = a(1)*sin(angs(2)) - a(2)*cos(angs(2))
        if (abs(a2) .lt. 0.0001) a2 = 0.0
        if (abs(sin(2*angs(2))) .ge. 0.0001) then
          angs(3) = atan2(a2/sin(2*angs(2)),a1)
        else if (abs(sin(angs(2))) .ge. 0.0001) then
          acosarg = amin1(1.0,amax1(-1.0,a(2)/sin(angs(2))))
          angs(3) = acos(acosarg)
        else
          acosarg = amin1(1.0,amax1(-1.0,a1))
          angs(3) = acos(a1)
        end if
      end if
      IF (ANGS(1) .lt. 0.0) then
        ANGS(1) = ANGS(1) + PI
        ANGS(3) = PI - ANGS(3)
        IF (ANGS(3) .GT. PI) ANGS(3) = ANGS(3) - 2*PI
      end if
      IF(ANGS(1) .gt. 0.5*PI) then
        ANGS(1)=PI-ANGS(1)
        ANGS(2)=ANGS(2)+PI
        ANGS(3)=-ANGS(3)
        IF (ANGS(2) .GE. 2*PI) ANGS(2) = ANGS(2) - 2*PI
      end if
      IF (ANGS(2) .LT. 0.0) ANGS(2) = ANGS(2) + 2.0*PI
      RETURN
END subroutine AN2DSR

! subroutine 4
subroutine AN2MOM(A,N,MOMTEN)
implicit none
      INTEGER J
      REAL*4 A(3), N(3), MOMTEN(6)
!           Moment tensor components:  M(I,j) = A(I)*N(J)+A(J)*N(I)
      MOMTEN(1) = 2.0*A(3)*N(3)     !  MRR = M(3,3)
      MOMTEN(2) = 2.0*A(1)*N(1)     !  MTT = M(1,1)
      MOMTEN(3) = 2.0*A(2)*N(2)     !  MPP = M(2,2)
      MOMTEN(4) = A(1)*N(3)+A(3)*N(1)     !  MRT = M(1,3)
      MOMTEN(5) = -A(2)*N(3)-A(3)*N(2)!  MRP = -M(2,3)
      MOMTEN(6) = -A(2)*N(1)-A(1)*N(2)!  MTP = -M(2,1)   
      DO 105 J=1,6
        IF (ABS(MOMTEN(J)) .LT. 0.000001) MOMTEN(J) = 0.0
105   CONTINUE
      RETURN
END subroutine AN2MOM

! subroutine 5
subroutine Getline(iFileUnit,GetFileN)
implicit none

Integer ,Intent(IN) :: iFileUnit 
Integer :: stat,GetFileN
Character(Len=1) :: cDummy
GetFileN = 0
Rewind( iFileUnit )
    Do
      Read( iFileUnit , * , iostat = stat ) cDummy
      if ( stat /= 0 ) Exit
      GetFileN = GetFileN + 1
    End Do
    Rewind( iFileUnit )
    Return 
End subroutine Getline