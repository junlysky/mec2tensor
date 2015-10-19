Subroutine TIMDAT(NOUT,PROGNM)

! time and date 
!      INPUT IS PROGRAM NAME
!         A LINE TO NOUT WITH DATE AND TIME AND PROGNM
!	unix version:  25 June 1991  jas/vtso

  character*24 fdate
  CHARACTER*(*) PROGNM
  if (nout .le. 0) return
  write(nout,*) fdate(),' for program ',prognm(1:lenc(prognm))
  RETURN
END SUBROUTINE TIMDAT
