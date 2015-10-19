character*(*) function getstring(prompt)

!  outputs 'prompt' using PRINTX
!  and accepts input character string
!				Alan Linde ... Aug 1986
!       27 July 1993: Did input read through cstring so can have 
!         comment lines
!	12 February 95:  Kill leading blanks
!
  character*(*) prompt
  character*80 temp
  getstring = ' '

  call printx(prompt)
  kk=lenc(prompt)
  if (prompt(kk:kk).eq.']') then
     ll=0
     do i=kk-1,1,-1
	    if (prompt(i:i).eq.'['.and.ll.eq.0) ll=i+1
  end do
  if (ll.ne.0) getstring=prompt(ll:kk-1)
end if
!  get the response
call cstring(temp,nout)
!  Kill leading blanks
do while (nout.gt.1 .and. temp(1:1).eq.' ')
   nout = nout - 1
   temp(1:nout) = temp(2:nout+1)
   temp(nout+1:nout+1) = ' '
end do
if (nout .gt. 0) getstring=temp(1:nout)
return
end function getstring
