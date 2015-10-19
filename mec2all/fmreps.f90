Subroutine FMREPS(ANBTP,ANGS,PTTP,ANGS2,AN,PT,DSR,mt,MOMTEN,eigen)
  implicit none

!	Input: A and N (trned and plunge), or P and T or dip, strike
!	  and rake (A&R convention) or moment tensor (A&R)
!	Output: the other representations plus the auxiliary plane.
!	PTTP:  4 parameters, trend and plunge for P and T
!	  P trend and plunge, T trend and plunge
!	ANBTP  6 parameters, t and p for A, N and B respectively.
!	ANGS   3 parameters,  dip, strike and rake for first plane
!	ANGS2  3 parameters, dip strike and rake for auxiliary plane
!	AN, PT, DSR are LOGICAL variables which are true if
!	  A and N, P and T or dip-strike-rake are input
!	MOMTEN  6 parameters:  the moment tensor (unit scalar
!	  magnitude) D & W convention
!	Angles come in and go out in degrees.
!	If LUNIT1 and/or LUNIT2 are positive, the representations
!	  are written on those logical unit numbers.

logical AN,PT,DSR,MT
real*4 MOMTEN(6)
integer LUNIT1, LUNIT2
real PTTP(4),ANGS(3),ANGS2(3),ANBTP(6)
real*4 MRR, MTT, MPP, MRT, MRP, MTP
integer j,i
real rdeg,pi
real eigen(3)

RDEG = 45.0/ATAN(1.0)
PI = 4.0*ATAN(1.0)


if (MT) then
   mrr=MOMTEN(1)
   mtt=MOMTEN(2)
   mpp=MOMTEN(3)
   mrt=MOMTEN(4)
   mrp=MOMTEN(5)
   mtp=MOMTEN(6)
   call mt_in(MRR,MTT,MPP,MRT,MRP,MTP,PTTP,eigen)
   call PTTPIN(PTTP,ANGS,ANGS2,ANBTP,MOMTEN)

else IF (PT) THEN
   do 100 J=1,4
100	  PTTP(J) = PTTP(J)/RDEG
      call PTTPIN(PTTP,ANGS,ANGS2,ANBTP,MOMTEN,eigen)

elseif (DSR) then
do 300 J=1,3
300	  ANGS(J) = ANGS(J)/RDEG
   call DSRIN(ANGS,ANBTP,ANGS2,PTTP,MOMTEN,eigen)

end if
do 400 I=1,3
   ANGS(I) = ANGS(I)*RDEG
   ANGS2(I) = ANGS2(I)*RDEG
   PTTP(I) = PTTP(I)*RDEG
   ANBTP(I) = ANBTP(I)*RDEG
400	CONTINUE
   ANBTP(4) = ANBTP(4)*RDEG
   ANBTP(5) = ANBTP(5)*RDEG
   ANBTP(6) = ANBTP(6)*RDEG
   PTTP(4) = PTTP(4)*RDEG

   RETURN

END Subroutine FMREPS
