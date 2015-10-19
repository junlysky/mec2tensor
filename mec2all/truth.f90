Logical function TRUTH(MSG)
!implicit none
!PURPOSE:
!		ROUTINE ACCEPTS A MESSAGE (QUESTION) REQUIRING A Y/N RESPONSE
!		AND THE RETURN "TRUTH" IS SET:
!			TRUTH = .TRUE.	  IF RESPONSE IS Y(y)
!			TRUTH = .FALSE.                 N(n)
!
! ROUTINES CALLED:
!			IYESNO
!				WHICH CALLS
!						NSTRNG
!						PRINTX
!
!
! USE:
!	I=TRUTH('ANSWER Y OR N')
!	IF (I) ......
!
! OR
!	IF (TRUTH('REPLY Y OR N')) ....
!
!
! AUTHOR:			ALAN LINDE ... AUGUST 1980
!
!  ENTRY
!			ILOGIC (Alan's original name)
!-
Character*(*) MSG
Logical	ILOGIC

!integer IANS

Entry   ILOGIC(MSG)
TRUTH=.FALSE.
call IYESNO(MSG,IANS)
if (IANS.EQ.1) TRUTH=.TRUE.
ILOGIC = TRUTH
return
end function TRUTH
