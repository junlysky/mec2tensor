!junlsky.info@gmail.com

! the subroutines are: 
! 1. eig(a,u,w)
! 2. tred2(nm,n,a,d,e,z)
! 3. imtql2(nm,n,d,e,z,ierr)
! 4. v2trpl(XYZ,PTTP(3),PI)  
! 5. tensor2pt(PTTP,MOMTEN,PI)
! 6. pt2mec(PTTP,ANGS,ANGS2,ANBTP,MOMTEN,PI)

Program main
implicit none

real,parameter :: PI=3.141592653
real :: lon,lat,MOMTEN(6)
character(50) line(8)
real ANGS(3),ANGS2(3),ANBTP(6),PTTP(4)
integer i,j,total,stat
real RDEG

RDEG  = 45.0/ATAN(1.0)

open(521,file='ten2mec.out')
open(520,file='tensor.dat',status='old',iostat=stat)
if(stat/=0) then
   write(*,*) "open tensor.dat fail ."
   stop
endif
Call Getline(520,total)
write(521,"(22A8)")"lon","lat","str1","dip1","rake1","str2","dip2","rake2","mrr","mtt","mpp","mrt","mrp","mtp",&
&"p-trend","p-plu","b-trend","b-plu","t-trend","t-plu"

do i=1,total
   read(520,*)lon,lat,MOMTEN(1:6)
   write(*,81)lon,lat,MOMTEN
   call tensor2pt(MOMTEN,PTTP)
   call pt2mec(PTTP,ANGS,ANGS2,ANBTP,PI)

   do j=1,3
      ANGS(j) = ANGS(j)*RDEG
      ANGS2(j) = ANGS2(j)*RDEG
      PTTP(j) = PTTP(j)*RDEG
      ANBTP(j) = ANBTP(j)*RDEG
   end do
   ANBTP(4) = ANBTP(4)*RDEG
   ANBTP(5) = ANBTP(5)*RDEG
   ANBTP(6) = ANBTP(6)*RDEG
   PTTP(4) = PTTP(4)*RDEG
   
   write(521,81)lon,lat,ANGS2(2),ANGS2(1),ANGS2(3),ANGS(2),ANGS(1),ANGS(3),MOMTEN,PTTP(1),PTTP(2),ANBTP(5),ANBTP(6),PTTP(3),PTTP(4)
81 format(22F8.2)
end do

close(520)
close(521)
end program main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tensor2pt(MOMTEN,PTTP)
implicit none
!	eigenvalues/vectors using EISPACK routines from www.netlib.no
!	Much of code adapted from Jost/Herrmann mteig.f and mtdec.f
!	Uses original EISPACK routines for TRED2 and IMTQL2, not NR
!	Also includes subroutine eig, which calls TRED2 and IMTQL2.
!	15 March 2002
!-
real,parameter :: PI=3.141592653
real :: A(3,3), U(3,3), W(3), PTTP(4), XYZ(3),MOMTEN(6)
real ::  MRR, MTT, MPP, MRT, MRP, MTP
integer :: i,j

MRR=MOMTEN(1)
MTT=MOMTEN(2)
MPP=MOMTEN(3)
MRT=MOMTEN(4)
MRP=MOMTEN(5)
MTP=MOMTEN(6)

!Convention is X north, Y east, Z down
A(3,3) = MRR
A(1,1) = MTT
A(2,2) = MPP
A(1,3) = MRT
A(3,1) = A(1,3)
A(2,3) = -MRP
A(3,2) = A(2,3)
A(1,2) = -MTP
A(2,1) = A(1,2)
call eig(a,u,w)

!Get trend and plunge for P
do j=1,3
   xyz(j) = u(j,1)
end do
call V2TRPL(XYZ,PTTP(1),PI)
!Get trend and plunge for T
do j=1,3
   xyz(j) = u(j,3)
end do
call V2TRPL(XYZ,PTTP(3),PI)

end subroutine tensor2pt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pt2mec(PTTP,ANGS,ANGS2,ANBTP,PI)
!	Calculates other representations of fault planes with
!		trend and plunge of P and T as input.  All
!		angles are in radians.
real N(3),MOMTEN(6)
real PTTP(4),ANGS(3),ANGS2(3),ANBTP(6),P(3),T(3),A(3),B(3)
DATA SR2/0.707107/
call TRPL2V(PTTP(1),P)
call TRPL2V(PTTP(3),T)
do j=1,3
   A(j) = SR2*(P(j) + T(j))
   N(j) = SR2*(T(j) - P(j))
end do
B(1) = P(2)*T(3) - P(3)*T(2)
B(2) = P(3)*T(1) - P(1)*T(3)
B(3) = P(1)*T(2) - P(2)*T(1)
call V2TRPL(A,ANBTP(1),PI)
call V2TRPL(N,ANBTP(3),PI)
call V2TRPL(B,ANBTP(5),PI)
call AN2DSR(A,N,ANGS,PI)
call AN2DSR(N,A,ANGS2,PI)
return 

end subroutine pt2mec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eig (a,u,w)
real ::  A(3,3), U(3,3), W(3), work(3)
integer i,j,n,np
np=3
do i=1,3
   do j=1,3
      u(i,j) = a(i,j)
   end do
end do
n = 3
call tred2(np,n,a,w,work,u)
call imtql2(np,n,w,work,u,ierr)
!This system has P, B, T as a right-hand coordinate system

do j=1,3
   u(j,1) = -u(j,1)
end do
return
end subroutine eig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine v2trpl(XYZ,TRPL,PI)
!
!    the trend and plunge for the vector.
!  Trend is the azimuth (clockwise from north looking down)
!  Plunge is the downward dip measured from the horizontal.
!  All angles in radians
!  X is north, Y is east, Z is down
!  If the component of Z is negative (up), the plunge,TRPL(2),
!    is replaced by its negative and the trend, TRPL(1),
!    Is changed by PI.
!  The trend is returned between 0 and 2*PI, the plunge
!    between 0 and PI/2.
!  12 January 2000: If xyz(3) = -1.0, make the trend PI.  Made
!    consistency in the roundoff -- all are now 0.0001
!-
real :: XYZ(3),TRPL(2)
do j=1,3
   if (abs(xyz(j)) .le. 0.0001) xyz(j) = 0.0
   if (abs(abs(xyz(j))-1.0).LT.0.0001) xyz(j)=xyz(j)/abs(xyz(j))
end do
if (abs(XYZ(3)) .eq. 1.0) then
!
!  plunge is 90 degrees
!
   if (xyz(3) .lt. 0.0) then
      trpl(1) = PI
   else
      trpl(1) = 0.0
   end if
   trpl(2) = 0.5*PI
   return
end if
if (abs(XYZ(1)) .LT. 0.0001) then
   if (XYZ(2) .GT. 0.0) then
      trpl(1) = PI/2.
   elseif  (XYZ(2) .LT. 0.0) then
      trpl(1) = 3.0*PI/2.0
   else
      trpl(1) = 0.0
   end if
else
   trpl(1) = ATAN2(XYZ(2),XYZ(1))
end if
C = COS(trpl(1))
S = SIN(trpl(1))
if (ABS(C) .GE. 0.1) trpl(2) = ATAN2(XYZ(3),XYZ(1)/C)
if (ABS(C) .LT. 0.1) trpl(2) = ATAN2(XYZ(3),XYZ(2)/S)
if (trpl(2) .LT. 0.0) THEN
   trpl(2) = -trpl(2)
   trpl(1) = trpl(1) - PI
end if
if (trpl(1) .LT. 0.0) trpl(1) = trpl(1) + 2.0*PI
return

end subroutine v2trpl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
SUBROUTINE TRPL2V(TRPL,XYZ)
!	Transforms to XYZ components of a unit vector from
!		the trend and plunge for the vector.
!	Trend is the azimuth (clockwise from north looking down)
!	Plunge is the downward dip measured from the horizontal.
!	All angles in radians
!	X is north, Y is east, Z is down
!-
real XYZ(3),TRPL(2)
XYZ(1) = COS(TRPL(1))*COS(TRPL(2))
XYZ(2) = SIN(TRPL(1))*COS(TRPL(2))
XYZ(3) = SIN(TRPL(2))
do j=1,3
   if (abs(xyz(j)) .lt. 0.0001) xyz(j) = 0.0
   if (abs(abs(xyz(j))-1.0).lt.0.0001) xyz(j)=xyz(j)/abs(xyz(j))
end do
RETURN
END SUBROUTINE TRPL2V



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!subroutine
subroutine tred2(nm,n,a,d,e,z)

integer i,j,k,l,n,nm
real a(nm,n),d(n),e(n),z(nm,n)
real f,g,h,hh,scale
!
!     this subroutine is a translation of the algol procedure tred2,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix to a
!     symmetric tridiagonal matrix using and accumulating
!     orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        z contains the orthogonal transformation matrix
!          produced in the reduction.
!
!        a and z may coincide.  if distinct, a is unaltered.
!
!     Questions and comments should be directed to Alan K. Cline,
!     Pleasant Valley Software, 8603 Altus Cove, Austin, TX 78759.
!     Electronic mail to cline@cs.utexas.edu.
!
!     this version dated january 1989. (for the IBM 3090vf)
!
!     ------------------------------------------------------------------
do i =1,n
   do j=1,n
      z(j,i) = a(j,i)
   end do
end do
do i = 1,n
   d(i) = a(n,i)
end do

! do 300 i = n, 2, -1
do i = n,2,-1
   l = i - 1
   h = 0.0e0
   scale = 0.0e0
   if (l .lt. 2) go to 130
   !     .......... scale row (algol tol then not needed) ..........
   do k = 1, l
      scale = scale + abs(d(k))
   end do
   if (scale .ne. 0.0e0) go to 140
130 e(i) = d(l)
   !
   !"    ( ignore recrdeps
   do j=1,l
      d(j) = z(l,j)
      z(i,j) = 0.0e0
      z(j,i) = 0.0e0
   end do
   !
   go to 290
   !
140 do     k = 1, l
      d(k) = d(k) / scale
      h = h + d(k) * d(k)
   end do          !
   f = d(l)
   g = -sign(sqrt(h),f)
   e(i) = scale * g
   h = h - f * g
   d(l) = f - g
   !     .......... form a*u ..........
   do     j = 1, l
      e(j) = 0.0e0
   end do
   !
   do     j = 1, l
      f = d(j)
      z(j,i) = f
      g = e(j) + z(j,j) * f
      !
      do     k = j+1, l
         g = g + z(k,j) * d(k)
         e(k) = e(k) + z(k,j) * f
      end do
      !
      e(j) = g
   end do
   !     .......... form p ..........
   f = 0.0e0
   !
   do     j = 1, l
      e(j) = e(j) / h
      f = f + e(j) * d(j)
   end do
   hh = -f / (h + h)
   !     .......... form q ..........
   do     j = 1, l
      e(j) = e(j) + hh * d(j)
   end do
   !     .......... form reduced a ..........
   do     j = 1, l
      f = -d(j)
      g = -e(j)
      !
      do     k = j, l
         z(k,j) = z(k,j) + f * e(k) + g * d(k)
      end do
      d(j) = z(l,j)
      z(i,j) = 0.0e0
   end do
290 d(i) = h
end do

!     .......... accumulation of transformation matrices ..........
do     i = 2, n
   l = i - 1
   z(n,l) = z(l,l)
   z(l,l) = 1.0e0
   h = d(i)
   if (h .eq. 0.0e0) go to 380
   !
   do     k = 1, l
      d(k) = z(k,i) / h
   end do
   !"    ( ignore recrdeps
   !"    ( prefer vector
   do     j = 1, l
      g = 0.0e0
      do     k = 1, l
         g = g + z(k,i) * z(k,j)
      end do
      g = -g
      do     k = 1, l
         z(k,j) = z(k,j) + g * d(k)
      end do
   end do
380 do     k = 1, l
      z(k,i) = 0.0e0
   end do
end do

do     i = 1, n
   d(i) = z(n,i)
   z(n,i) = 0.0e0
end do
!
z(n,n) = 1.0e0
e(1) = 0.0e0
return
end subroutine  tred2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!subroutine
subroutine imtql2(nm,n,d,e,z,ierr)

integer i,j,k,l,m,n,nm,ierr
real d(n),e(n),z(nm,n)
real b,c,f,g,p,r,s,tst1,tst2
!
!     this subroutine is a translation of the algol procedure imtql2,
!     num. math. 12, 377-383(1968) by martin and wilkinson,
!     as modified in num. math. 15, 450(1970) by dubrulle.
!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the implicit ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     Questions and comments should be directed to Alan K. Cline,
!     Pleasant Valley Software, 8603 Altus Cove, Austin, TX 78759.
!     Electronic mail to cline@cs.utexas.edu.
!
!     this version dated january 1989. (for the IBM 3090vf)
!
!     ------------------------------------------------------------------
!

ierr = 0
if (n .eq. 1) go to 1001

do     i = 2, n
   e(i-1) = e(i)
end do
e(n) = 0.0e0
do  240   l = 1, n
   j = 0
!     .......... look for small sub-diagonal element ..........
105 do     m = l, n-1
      tst1 = abs(d(m)) + abs(d(m+1))
      tst2 = tst1 + abs(e(m))
      if (tst2 .eq. tst1) go to 120
   end do
120 p = d(l)
   if (m .eq. l) go to 240
   if (j .eq. 30) go to 1000
   j = j + 1
   !     .......... form shift ..........
   g = (d(l+1) - p) / (2.0e0 * e(l))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! *      r = pythag(g,1.0d0)
   !ccccccccccccccccccccccccccccccccccccccccccccccccccc
   if (abs(g).le.1.0e0) then
      r = sqrt(1.0e0 + g*g)
   else
      r = g * sqrt(1.0e0 + (1.0e0/g)**2)
   endif
   !ccccccccccccccccccccccccccccccccccccccccccccccccccc
   g = d(m) - p + e(l) / (g + sign(r,g))
   s = 1.0e0
   c = 1.0e0
   p = 0.0e0
   !.......... for i=m-1 step -1 until l do -- ..........
   do  i = m-1, l, -1
      f = s * e(i)
      b = c * e(i)
      !ccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! *         r = pythag(f,g)
      !ccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (abs(f).ge.abs(g)) then
         r = abs(f) * sqrt(1.0e0 + (g/f)**2)
      else if (g .ne. 0.0e0) then
         r = abs(g) * sqrt((f/g)**2 + 1.0e0)
      else
         r = abs(f)
      endif
      !ccccccccccccccccccccccccccccccccccccccccccccccccccc
      e(i+1) = r
      if (r .eq. 0.0e0) then
         !     .......... recover from underflow ..........
         d(i+1) = d(i+1) - p
         e(m) = 0.0e0
         go to 105
      endif
      s = f / r
      c = g / r
      g = d(i+1) - p
      r = (d(i) - g) * s + 2.0e0 * c * b
      p = s * r
      d(i+1) = g + p
      g = c * r - b
      !   .......... form vector ..........
      do     k = 1, n
         f = z(k,i+1)
         z(k,i+1) = s * z(k,i) + c * f
         z(k,i) = c * z(k,i) - s * f
      end do
   end do

   d(l) = d(l) - p
   e(l) = g
   e(m) = 0.0e0
   go to 105
240 continue 
   !     .......... order eigenvalues and eigenvectors ..........
   do     i = 1, n-1
      k = i
      p = d(i)

      do  260   j = i+1, n
         if (d(j) .ge. p) go to 260
         k = j
         p = d(j)
260      continue 
         d(k) = d(i)
         d(i) = p

         do     j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
         end do
      end do

      go to 1001
      !     .......... set error -- no convergence to an
      !                eigenvalue after 30 iterations ..........
1000  ierr = l
1001  return

end subroutine imtql2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine AN2DSR(A,N,ANGS,PI)
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
END subroutine AN2DSR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
