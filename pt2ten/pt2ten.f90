!junlysky.info@gail.com
! usage: compute the tensor elements using p-azi,p-plunge,t-azi and t-plunge
!input file format:lon lat p-azi p-plu t-azi t-plu 

program main
implicit none

  real,parameter :: PI = 3.1415926
  real :: lon,lat,PTTP(4)
  
  real :: ANGS(3),ANGS2(3),MOMTEN(6)
  real :: N(3),ANBTP(6),P(3),T(3),A(3),B(3)
  character(100) :: line(6)
  integer i,j,total,stat,k
real SR2,RDEG
DATA SR2/0.707107/
RDEG  = 45.0/ATAN(1.0)

  open(521,file='pt2tensor.out')
  open(520,file='ptin.dat',status='old',iostat=stat)
  if(stat/=0) then
     write(*,*) "open ptin.dat fail."
     stop
  endif
  call getline(520,total)
  write(521,"(22A8)")"lon","lat","str1","dip1","rake1","str2","dip2","rake2",&
       "mrr","mtt","mpp","mrt","mrp","mtp",&
       "p-trend","p-plu","b-trend","b-plu","t-trend","t-plu"

  do i=1,total
!     read(520,*) lon,lat,PTTP(3),PTTP(4),PTTP(1),PTTP(2)  !pttp(1) p-trend, pttp(2) p-plu
     read(520,*)lon,lat,PTTP
     write(*,88)lon,lat,PTTP
     do j = 1,4
        PTTP(j)=PTTP(j)/RDEG
     end do
     call pt2mec(PTTP,ANGS,ANGS2,ANBTP,MOMTEN,PI)
     do k = 1,3
        ANGS(k) = ANGS(k)*RDEG
        ANGS2(k) = ANGS2(k)*RDEG
        PTTP(k) = PTTP(k)*RDEG
        ANBTP(k) = ANBTP(k)*RDEG
     end do
     ANBTP(4) = ANBTP(4)*RDEG
     ANBTP(5) = ANBTP(5)*RDEG
     ANBTP(6) = ANBTP(6)*RDEG
     PTTP(4) = PTTP(4)*RDEG
!     write(521,81)lon,lat,PTTP(1),PTTP(2),PTTP(3),PTTP(4),MOMTEN
   write(521,81)lon,lat,ANGS(2),ANGS(1),ANGS(3),ANGS2(2),ANGS2(1),ANGS2(3),&
        MOMTEN,PTTP(1),PTTP(2),ANBTP(5),ANBTP(6),PTTP(3),PTTP(4)

!     write(521,81)lon,lat,ANGS(2),ANGS(1),ANGS(3),ANGS2(2),ANGS2(1),ANGS2(3),MOMTEN,ANBTP(5),ANBTP(6),PTTP
81 format(23F8.2)
88 format(6F8.2)
  end do
  close(520)
  close(521)
end program main

subroutine pt2mec(PTTP,ANGS,ANGS2,ANBTP,MOMTEN,PI)
implicit none
! 1. all angles are in radians
! 2. the input data are p-azi p-plunge t-azi t-plunge
! 3. PTTP(1) is p-trend(also p-azi), PTTP(2) is p-plunge,
!    similar PTTP(3),PTTP(4)
! 4. ANBTP(5) is b-azi, ANBTP(6) is b-plu
real N(3),MOMTEN(6),pi
real PTTP(4),ANGS(3),ANGS2(3),ANBTP(6),P(3),T(3),A(3),B(3)
real SR2,RDEG
integer i,j
DATA SR2/0.707107/
RDEG  = 45.0/ATAN(1.0)

call TRPL2V(PTTP(1),P)     ! P(1:3) ->  x y z
call TRPL2V(PTTP(3),T)     ! T(1:3) ->  x y z
do j=1,3
   A(j) = SR2*(P(j) + T(j)) 
   N(j) = SR2*(T(j) - P(j))
end do
B(1) = P(2)*T(3) - P(3)*T(2)
B(2) = P(3)*T(1) - P(1)*T(3)
B(3) = P(1)*T(2) - P(2)*T(1)
CALL V2TRPL(A,ANBTP(1))
CALL V2TRPL(N,ANBTP(3))
CALL V2TRPL(B,ANBTP(5))
CALL AN2DSR(A,N,ANGS,PI)
CALL AN2DSR(N,A,ANGS2,PI)
CALL AN2MOM(A,N,MOMTEN)
return
end subroutine pt2mec

subroutine TRPL2V(TRPL,XYZ)
implicit none
! Transforms to XYZ components of a unit vector from
! the trend and plunge for the vector.
! trend is azimuth ( clockwise from north looking down)
! plunge is the downward dip measured from the horizontal
! all angles in radians
! X is north，Y is east,Z is down
! XYZ(1) is x(north),XYZ(2) is y(east),XYZ(3) is z(down)
! trpl(1) is azimuth angle
! trpl(2) is plunge angle

  real XYZ(3),TRPL(2)
  integer j
  XYZ(1) = COS(TRPL(1))*COS(TRPL(2))
  XYZ(2) = SIN(TRPL(1))*COS(TRPL(2))
  XYZ(3) = SIN(TRPL(2))
  do j=1,3
     if (abs(xyz(j)) .lt. 0.0001) xyz(j) = 0.0
     if (abs(abs(xyz(j))-1.0).lt.0.0001) xyz(j)=xyz(j)/abs(xyz(j))
  end do
  return
end subroutine TRPL2V

subroutine V2TRPL(XYZ,TRPL)
implicit none
!If the component of Z is negative (up), the plunge,TRPL(2),
! is replaced by its negative and the trend, TRPL(1),
! Is changed by PI.
! The trend is returned between 0 and 2*PI, the plunge between 0 and PI/2.
! If xyz(3) = -1.0, make the trend PI,made consistency in the roundoff -- all are now 0.0001
real XYZ(3),TRPL(2)
real c,s
real,parameter :: PI = 3.141592653
integer j

do j=1,3
   if (abs(xyz(j)) .le. 0.0001) xyz(j) = 0.0
   IF (ABS(ABS(XYZ(j))-1.0).LT.0.0001) xyz(j)=xyz(j)/abs(xyz(j))
end do
IF (ABS(XYZ(3)) .eq. 1.0) THEN 
!	plunge is 90 degrees
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
end subroutine V2TRPL


subroutine AN2DSR(A,N,ANGS,PI)
  
! AN2DSR(1) is dip,AN2DSR（2）is strike,AN2DSR(3) is rake
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
end subroutine AN2DSR

subroutine AN2MOM(A,N,MOMTEN)
implicit none
integer j
real A(3), N(3), MOMTEN(6)
!Moment tensor components:  M(I,j) = A(I)*N(J)+A(J)*N(I)
MOMTEN(1) = 2.0*A(3)*N(3)     !  MRR = M(3,3)
MOMTEN(2) = 2.0*A(1)*N(1)     !  MTT = M(1,1)
MOMTEN(3) = 2.0*A(2)*N(2)     !  MPP = M(2,2)
MOMTEN(4) = A(1)*N(3)+A(3)*N(1)     !  MRT = M(1,3)
MOMTEN(5) = -A(2)*N(3)-A(3)*N(2)!  MRP = -M(2,3)
MOMTEN(6) = -A(2)*N(1)-A(1)*N(2)!  MTP = -M(2,1)   
do J=1,6
   IF (ABS(MOMTEN(J)) .LT. 0.000001) MOMTEN(J) = 0.0
end do
return
end subroutine AN2MOM

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
 