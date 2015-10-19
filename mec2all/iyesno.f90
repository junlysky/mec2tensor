SUBROUTINE IYESNO(MSG,IANS)

!     PURPOSE:
!	     THIS LITTLE SUBROUTINE ASKS A QUESTION AND RETURNS A
!	     RESPONSE TO THAT QUESTION. THE ANSWER TO THE QUESTION
!	     MUST BE EITHER 'Y' FOR YES, 'N' FOR NO, OR NOTHING
!	     (i.e. simply hitting carrage return) FOR THE DEFAULT
!	     REPONSE TO THE QUESTION.
!
!     ON INPUT:
!	    MSG = BYTE STRING CONTAINING THE QUESTION
!
!     ON OUTPUT:
!	    IANS = THE LOGICAL REPONSE TO THE QUESTION (1 or 0)
!     EXTRA FEATURES:
!	    DEFAULT SITUATION IS:
!	    IF LAST 3 CHARACTERS IN 'MSG' ARE
!	  	     [Y]  OR  [N]
!	    THEN 'IANS' = 1   OR   0
!
!	    IF LAST 3 CHARACTERS ARE NOT ONE OF ABOVE PAIRS
!	    THEN 'IANS' = 0
!	    (i.e. default for no supplied default is N)
!	30 JULY 1989:  IF ENTERED CHARACTER IS A BLANK OR A TAB, 
!	    TREATS AS A NULL ENTRY.
!       27 July 1993: Did input read through cstring so can have 
!         comment lines

CHARACTER*1 DELIM/'$'/,CHARIN,BLANK/' '/
CHARACTER*3 TEST,UCY,LCY
character*80 string_in
CHARACTER*(*) MSG



DATA UCY/'[Y]'/,LCY/'[y]'/
KK = LEN(MSG)
IF (MSG(KK:KK) .EQ. DELIM) KK = KK - 1
TEST = MSG(KK-2:KK)
CALL PRINTX(MSG)
call cstring(string_in,nchar)
IF ((NCHAR.GT.0) .AND. (string_in(1:1).EQ.BLANK)) NCHAR = 0
IF (NCHAR .EQ. 0) THEN
   IF ((TEST .EQ. UCY) .OR. (TEST .EQ. LCY)) THEN
      IANS = 1
   ELSE 
      IANS = 0
   END IF
ELSE
   charin = string_in(1:1)
   IF (CHARIN .EQ. UCY(2:2) .OR. CHARIN .EQ. LCY(2:2)) THEN
      IANS = 1
   ELSE
      IANS = 0
   END IF
END IF
RETURN
END SUBROUTINE IYESNO
