program mec2all
  implicit none
!	  Based on the Aki & Richards convention, gives all 
!	  representations of a focal mechanism for input
!	  Dip, Strike, Rake    or
!	  A and N trend and plunge    or
!	  P and T trend and plunge.
!
!	  The code is Fortran 90
!	  Subroutine specific to this program are
!	  FMREPS,PTTPIN,ANTPIN,DSRIN,AN2MOM,V2TRPL,TRPL2V,AN2DSR, MT_IN

Real,Parameter :: PI = 3.1415926
Logical AN,PT,DSR,TRUTH,first,MT
Real ANGS(3),eigen(3)
Real ANBTP(6),ANGS2(3),PTTP(4),MOMTEN(6)
Character(80) getstring,commnt
integer j
real*4 MRR, MTT, MPP, MRT, MRP, MTP

integer,external::lenc
character(30) usrname
integer sta
character(10) eyear,emon,eday,ehour,emin,esec
character(30) evtime
real lon,lat,depth,mag          !mag is Mw
real Mo                        ! Mo =(10**(1.5 * Mw)) * 10**(9.1)
integer EE                   ! EE = log(Mo)



!open(unit=2,file='new.dat',status='unknown')
open(unit=2,file='new.dat')

write(*,*) "Input the usr name:"
read(*,*) usrname
write(*,*) "input the origin time(year month day hour min sec):"
read(*,*) eyear,emon,eday,ehour,emin,esec
write(*,*)"input the epcenter location(lon  lat): " 
read(*,*) lon,lat
write(*,*)"input the depth and magnitude of the earthquake (depth  magnitude) :"
read(*,*) depth,mag
write(*,*) "input the number of stations"
read(*,*) sta

!501 format("Timing and location information")
500 format(A4,"/",A2,"/",a2,"  ",a2,":",a2,":",a2)
501 format("Epicenter:",F6.2," ",F6.2)
502 format("CENC MOMENT TENSOR SOLUTION")
503 format(A15)
504 format("Depth=",F4.1,2X,"No. of sta=",I3)
505 format("Mw=",F3.1)

write(2,500) eyear,emon,eday,ehour,emin,esec
write(2,501) lat,lon
write(2,*) " "
write(2,502) 
write(2,503) usrname
write(2,*) " "
write(2,504) depth,sta
write(2,505) mag     
write(2,*) " "

Mo =(10**(1.5 * mag)) * 10**(16.1)
ee = 1.5*mag + 16.1
ee = int(ee)

  PT = .FALSE.
  DSR = .FALSE.
  MT = .FALSE.
 
  if (truth('strike, dip and rake?..[Y]')) then
     DSR = .TRUE.
     call PRINTX('Enter strike, dip and Rake (degrees)')
     read(*,*) ANGS(2),ANGS(1),ANGS(3)

  elseif (truth('P and T axes trend and plunge?..[Y]')) THEN
     PT = .TRUE.
     call PRINTX ('Enter trend and plunge for P and T (t,p,t,p)')
     read(*,*) (PTTP(J),J=1,4)

  elseif  (truth('Moment-tensor input?...[Y]')) then
     MT = .true.
     call PRINTX ('Enter tensor (mrr,mtt,mpp,mrt,mrp,mtp)')
     read(*,*) MRR, MTT, MPP, MRT, MRP, MTP
        MOMTEN(1)=mrr
        MOMTEN(2)=mtt
        MOMTEN(3)=mpp
        MOMTEN(4)=mrt
        MOMTEN(5)=mrp
        MOMTEN(6)=mtp
      
  end if

506 FORMAT("Moment Tensor;  Scale 10**",I2,"dynes-cm")
507 format("MRR = ",F5.2,"  ","MTT = ",F5.2)
508 format("MPP = ",F5.2,"  ","MRT = ",F5.2)
509 format("MRP = ",F5.2,"  ","MTP = ",F5.2)
!510 format(*,*) " "
511 format("Principal axes:")
512 format(6X,"azi      plg    eigen")
513 FORMAT(X,"T ",3F8.2)
514 format(X,"N ",3F8.2)
515 format(X,"P ",3F8.2)
!516 format(*,*) " "
517 format("Best Double Couple:Mo=",E10.3,"dynes-cm")  
518 format("NP1:strike=",F8.2,2X,"dip=",F8.2,2X,"rake=",F8.2)
519 format("NP2:strike=",F8.2,2X,"dip=",F8.2,2X,"rake=",F8.2)

write(*,*) " "
write(*,*) " "
write(*,500) eyear,emon,eday,ehour,emin,esec
write(*,501) lat,lon
write(*,*) " "
write(*,502) 
write(*,503) usrname
write(*,*) " "
write(*,504) depth,sta
write(*,505) mag
write(*,*) " "
  
!  call FMREPS(ANBTP,ANGS,PTTP,ANGS2,AN,PT,DSR,mt,MOMTEN,eigen,2,6)
  call FMREPS(ANBTP,ANGS,PTTP,ANGS2,AN,PT,DSR,mt,MOMTEN,eigen)

  write(2,506) ee
  write(2,507) MOMTEN(1),MOMTEN(2)
  write(2,508) MOMTEN(3),MOMTEN(4)
  write(2,509) MOMTEN(5),MOMTEN(6)
  write(2,*) " "
  write(2,511)
  write(2,512)
  write(2,513) (PTTP(j),j=3,4),eigen(3)
  write(2,514) (ANBTP(j),j=5,6),eigen(2)
  write(2,515) (PTTP(j),j=1,2),eigen(1)
  write(2,*) " "
  write(2,517)  Mo
  write(2,518)  ANGS(2),ANGS(1),ANGS(3)
  WRITE(2,519)  ANGS2(2),ANGS2(1),ANGS2(3)

  write(*,506) ee
  write(*,507) MOMTEN(1),MOMTEN(2)
  write(*,508) MOMTEN(3),MOMTEN(4)
  write(*,509) MOMTEN(5),MOMTEN(6)
  write(*,*) " "
  write(*,511)
  write(*,512)
  write(*,513) (PTTP(j),j=3,4),eigen(3)
  write(*,514) (ANBTP(j),j=5,6),eigen(2)
  write(*,515) (PTTP(j),j=1,2),eigen(1)
  write(*,*) " "
  write(*,517) Mo
  write(*,518)  ANGS(2),ANGS(1),ANGS(3)
  WRITE(*,519)  ANGS2(2),ANGS2(1),ANGS2(3)

  first = .true.

  call bball(momten,pttp(1),pttp(2),pttp(3),pttp(4),2,first)     ! plot beach ball in outfile
  call bball(momten,pttp(1),pttp(2),pttp(3),pttp(4),6,first)     ! plot beach ball in screem


  first = .true.


end program mec2all

!!!!!!!!!!!!!!!!!!!!!!!!!!subroutine strat !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bball(g,pazim,pplng,tazim,tplng,unit,first)

! ...... generate printer plot rendition of lower hemisphere 
!        equal area projection
!	g has the six elements of the moment tensor, the rest are the
!	  plunge and trends of the P and T axes in degrees. unit is the output
!	  unit.
!	From Stuart Sipkin and Bob Uhrhammer 1993
!	1 October 2001: Replaced his sind, etc, with sin as not all compilers
!		know about degree versions of sin, cos, etc.
!-
  dimension g(6)
  integer unit
  character*1 ach(39,72),aplus,aminus,apaxis,ataxis,ablank
  logical first
!
  data aplus,aminus,apaxis,ataxis,ablank /'#','-','P','T',' '/
  data radius /1.41/
!
! ...... construct lower hemisphere fps 
!
      save
      rd = 45.0/atan(1.0)
      r0=radius
      x0=r0+0.250
      y0=r0+0.500
      ix0=12.*x0
      iy0=6.5*y0
      do 3 i=1,2*ix0
      do 2 j=1,2*iy0
      dx=real(i-ix0)/12.
      dy=-real(j-iy0)/6.5
      dd=dx*dx+dy*dy
      if(dd.gt.0.) then
        del=sqrt(dd)
      else
        del=0.
      endif
      if((dx.eq.0.).and.(dy.eq.0.)) then
        theta=0.
      else
        theta=rd*atan2(dx,dy)
      endif
      if(del.gt.r0) then
        ach(j,i)=ablank
        go to 1
      endif
      if(del.ge.r0) then
        aoi=90.0
      else
        aoi=90.*del/r0
      endif
      if(polar(g,aoi,theta,first).gt.0.) then
        ach(j,i)=aplus
      else
        ach(j,i)=aminus
      endif
1     continue
2     continue
3     continue
!
! ...... add P & T axis
!
      ixp=nint(r0*12.*(90.-pplng)*sin(pazim/rd)/90.+real(ix0))
      iyp=nint(-r0*6.5*(90.-pplng)*cos(pazim/rd)/90.+real(iy0))
      do 5 i=ixp-1,ixp+1
         do 4 j=iyp-1,iyp+1
      ach(j,i)=ablank
4     continue
5     continue
      ach(iyp,ixp)=apaxis
      ixt=nint(r0*12.*(90.-tplng)*sin(tazim/rd)/90.+real(ix0))
      iyt=nint(-r0*6.5*(90.-tplng)*cos(tazim/rd)/90.+real(iy0))
      do 7 i=ixt-1,ixt+1
      do 6 j=iyt-1,iyt+1
      ach(j,i)=ablank
6     continue
7     continue
      ach(iyt,ixt)=ataxis
!
! ...... add fps plot
!
      do 8 i=1,2*iy0-2
      write(unit,'(72a1)') (ach(i,j),j=1,2*ix0)
8     continue
!
      return
      end

      real*4 function polar(g,aoi,theta,first)
!
! ...... compute first motion podsretc_CMT.lstlarity as a function of aoi & theta
!        for a moment tensor for a double-couple solution.
!	Conventions differ slightly from Sipkin.  My moment tensor is derived
!	  from the outer product of two vectors and is hence normalized.  The
!	  order is also different from his, apparently.  I also did not know
!	  cosd and sind existed.
!
      dimension g(6)
      real mxx,mxy,mxz,myy,myz,mzz
      logical first
!
      save
      rd = 45.0/atan(1.0)
      if(first) then
        mxx= g(2)
        mxy=-g(6)
        mxz= g(4)
        myy= g(3)
        myz=-g(5)
        mzz= g(1)
        first = .false.
      endif
	x = cos(theta/rd)*sin(aoi/rd)
	y = sin(theta/rd)*sin(aoi/rd)
	z = cos(aoi/rd)
!
      polar = x*mxx*x + 2*x*mxy*y + 2*x*mxz*z + 2*y*myz*z +y*myy*y+z*mzz*z
!
      return
      end


subroutine ang(anbtp)
!
!	Given A and B trend and plunge, find and print angle A makes
!	with the plane formed by the vertical and the trend of B.
!	See subroutines srchfm.f and fltsol.f for explanation
!-
	dimension anbtp(6), btrpl(2), xtrpl(2), atrpl(2), a(3),x(3)
	rddeg = 45.0/atan(1.0)
	atrpl(1) = anbtp(1)
	atrpl(2) = anbtp(2)
	btrpl(1) = anbtp(5)
	btrpl(2) = anbtp(6)
	xtrpl(1) = btrpl(1) + 180.0
	xtrpl(2) = 90.0 - btrpl(2)
	do j=1,2
	  atrpl(j) = atrpl(j)/rddeg
	  xtrpl(j) = xtrpl(j)/rddeg
	enddo
	call trpl2v(atrpl,a)
	call trpl2v(xtrpl,x)
	cosangle = 0.0
	do j=1,3
	  cosangle = cosangle + a(j)*x(j)
	enddo
	angle = rddeg *acos(cosangle)
	write(*,1) angle,90-angle
	write(2,1) angle,90-angle
1       FORMAT('"A" angle with plane defined by vertical & B trend:',&
        & F5.1,/,'"N" angle with plane defined by vertical & B trend:',F5.1)
	return
end subroutine ang
