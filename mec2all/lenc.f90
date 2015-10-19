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
