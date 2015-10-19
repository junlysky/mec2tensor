subroutine cstring(string,nstring)

!	Input a character string with a read(*,'(A)') string
!	If first two characters are /* it will read the next entry
!	Tab is a delimiter.
!	Returns string and nstring, number of characters to tab.
!	string stars with first non-blank character.
!       25 May 2001.  Took out parameter statement for tab.

  logical more
  CHARACTER*1 TAB
  CHARACTER*(*) string
!
  tab = char(9)
  more = .true.
  do while (more)
     read(*,'(A)') string
     nstring = lenc(string)
     more = (nstring.ge.2 .and. string(1:2).eq.'/*')
  end do
  IF (nstring .GT. 0) THEN
     NTAB = INDEX(string(1:nstring),TAB)
     IF (NTAB .GT. 0) nstring = NTAB - 1
  end if
  return
end subroutine cstring
