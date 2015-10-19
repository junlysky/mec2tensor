Subroutine Printx(line)
  implicit none
!  OUTPUTS A MESSAGE TO THE TERMINAL
!  PRINTX STARTS WITH A LINE FEED BUT DOES NOT END WITH A CARRIAGE RETURN
!  THE PRINT HEAD REMAINS AT THE END OF THE MESSAGE
!
!  IF THE MESSAGE LENGTH IS LESS THAN 40,
!	DOTS ARE INSERTED UP TO COL. 39
!	AND A COLON IS PUT IN COL. 40.

  Character*(*) line
  Character(60) BUF
  Character(2) COLON
  Character(1)  DOT,DELIM

  integer j,kk
  integer,external:: lenc
  
  DATA DELIM/'$'/,DOT/'.'/,COLON/': '/
  KK = lenc(line)	!  length minus right-hand blanks
  if (line(kk:kk) .EQ. DELIM) kk= kk - 1
  if (kk .GT. 58) kk = 59
  BUF(1:kk) = line(1:kk)
  if (kk .LT. 49) then
     do J=kk+1,49
	    BUF(J:J) = DOT
  end do
  kk = 49
end if
BUF(kk:kk+1) = COLON
  kk = kk + 1
 write(*,'(1x,A,$)') BUF(1:kk)
 return
end Subroutine Printx



integer function lenc(string)
implicit none
!	Returns length of character variable STRING excluding right-hand
!	  most blanks or nulls

character*(*) string

integer j,length

length = len(string)	! total length
if (length .eq. 0) then
   lenc = 0
   return
end if

if(ichar(string(length:length)).eq.0)string(length:length) = ' '
do j=length,1,-1
   lenc = j
   if (string(j:j).ne.' ' .and. ichar(string(j:j)).ne.0) return
end do
lenc = 0
return
end function lenc
