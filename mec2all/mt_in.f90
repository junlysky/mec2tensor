subroutine mt_in(MRR,MTT,MPP,MRT,MRP,MTP,pttp,eigen)
implicit none
!	eigenvalues/vectors using EISPACK routines from www.netlib.no
!	Much of code adapted from Jost/Herrmann mteig.f and mtdec.f
!	Uses original EISPACK routines for TRED2 and IMTQL2, not NR
!	Also includes subroutine eig, which calls TRED2 and IMTQL2.
!	15 March 2002
!
real, parameter :: pi = 3.1415926536 
real A(3,3), U(3,3), W(3), PTTP(4), XYZ(3)
real*4 MRR, MTT, MPP, MRT, MRP, MTP, isotrop

real eigen(3)
real devmom,eps,eps1,eps2
integer j

!write(*,*) 'Input MRR MTT MPP MRT MRP MTP (free format)'
!read(*,*) MRR, MTT, MPP, MRT, MRP, MTP

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

isotrop = 0.0      ! M_all = M_isotrop + M_deviatoric
do j=1,3
   isotrop = isotrop + w(j)
end do
devmom = 0.0
do j=1,3
   w(j) = w(j) - isotrop/3.0
   devmom = devmom + w(j)*w(j)
   eigen(j) = w(j)
end do
devmom = sqrt(0.5*devmom)
!if (devmom .lt. 0.001*isotrop) devmom = 0.0
!write(*,*) ' '

if (devmom .eq. 0.0) then
   write(2,*) 'Exiting because purely isotropic source'
   write(*,*) 'Exiting because purely isotropic source'
   stop
end if

eps = abs(w(2)/amax1(-w(1),w(3)))
eps1 = eps*200.0
eps2 = 100.0 - eps1

if (eps .ge. 0.25) then
   write(2,*) ' Exiting because less than 50% double couple'
   write(*,*) ' Exiting because less than 50% double couple'
   stop
end if

do j=1,3
   xyz(j) = u(j,1)
end do
call V2TRPL(XYZ,PTTP(1),PI)

!	Get trend and plunge for T

do j=1,3
   xyz(j) = u(j,3)
end do
call V2TRPL(XYZ,PTTP(3),PI)

return
end subroutine mt_in



!!!!!!!!!!!!!!!!!!!!!!subroutine start!!!!!!!!!!!!!!!!!!!
subroutine eig (a,u,w)
  real  A(3,3), U(3,3), W(3), work(3)
  np = 3
  do i=1,np
     do j=1,3
        u(i,j) = a(i,j)
     end do
  end do
  n = 3
  call tred2(np,n,a,w,work,u)
  call imtql2(np,n,w,work,u,ierr)
!
!	This system has P, B, T as a right-hand coordinate system.
!	I prefer P, T, B
!
  do j=1,3
     u(j,1) = -u(j,1)
  end do
  return
end subroutine eig

!---------------
subroutine tred2(nm,n,a,d,e,z)
!
  integer i,j,k,l,n,nm
  real a(nm,n),d(n),e(n),z(nm,n)
  real f,g,h,hh,scale
!
!     this subroutine is a translation of the algol procedure tred2,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!c
!c     this subroutine reduces a real symmetric matrix to a
!c     symmetric tridiagonal matrix using and accumulating
!c     orthogonal similarity transformations.
!c
!c     on input
!c
!c        nm must be set to the row dimension of two-dimensional
!c          array parameters as declared in the calling program
!c          dimension statement.
!c
!c        n is the order of the matrix.
!c
!c        a contains the real symmetric input matrix.  only the
!c          lower triangle of the matrix need be supplied.
!c
!c     on output
!c
!c        d contains the diagonal elements of the tridiagonal matrix.
!c
!c        e contains the subdiagonal elements of the tridiagonal
!c          matrix in its last n-1 positions.  e(1) is set to zero.
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
!
!?      call xuflow(0)
  do 100 i = 1, n
     do 100 j = i, n
100     z(j,i) = a(j,i)
!
  do 110 i = 1, n
110  d(i) = a(n,i)
!
  do 300 i = n, 2, -1
     l = i - 1
     h = 0.0e0
     scale = 0.0e0
     if (l .lt. 2) go to 130
!     .......... scale row (algol tol then not needed) ..........
  do 120 k = 1, l
120  scale = scale + abs(d(k))
!
     if (scale .ne. 0.0e0) go to 140
130  e(i) = d(l)
!
!"    ( ignore recrdeps
     do 135 j = 1, l
        d(j) = z(l,j)
        z(i,j) = 0.0e0
        z(j,i) = 0.0e0
135     continue

  go to 290

140 do 150 k = 1, l
     d(k) = d(k) / scale
     h = h + d(k) * d(k)
150  continue

     f = d(l)
     g = -sign(sqrt(h),f)
     e(i) = scale * g
     h = h - f * g
     d(l) = f - g
!     .......... form a*u ..........
  do 170 j = 1, l
170  e(j) = 0.0e0

  do 240 j = 1, l
     f = d(j)
     z(j,i) = f
     g = e(j) + z(j,j) * f

  do 200 k = j+1, l
     g = g + z(k,j) * d(k)
     e(k) = e(k) + z(k,j) * f
200  continue

  e(j) = g
240 continue
!     .......... form p ..........
  f = 0.0e0

  do 245 j = 1, l
     e(j) = e(j) / h
     f = f + e(j) * d(j)
245  continue

   hh = -f / (h + h)
!     .......... form q ..........
   do 250 j = 1, l
250   e(j) = e(j) + hh * d(j)
!     .......... form reduced a ..........
   do 280 j = 1, l
      f = -d(j)
      g = -e(j)

   do 260 k = j, l
260   z(k,j) = z(k,j) + f * e(k) + g * d(k)
!
      d(j) = z(l,j)
      z(i,j) = 0.0e0
280    continue

290   d(i) = h
300   continue
!     .......... accumulation of transformation matrices ..........
   do 500 i = 2, n
      l = i - 1
      z(n,l) = z(l,l)
      z(l,l) = 1.0e0
      h = d(i)
      if (h .eq. 0.0e0) go to 380

      do 330 k = 1, l
330      d(k) = z(k,i) / h
!"    ( ignore recrdeps
!"    ( prefer vector
      do 360 j = 1, l
      g = 0.0e0

      do 340 k = 1, l
340      g = g + z(k,i) * z(k,j)

  g = -g

  do 350 k = 1, l
350  z(k,j) = z(k,j) + g * d(k)
360  continue

380  do 400 k = 1, l
400     z(k,i) = 0.0e0

500     continue

!"    ( prefer vector
   do 520 i = 1, n
      d(i) = z(n,i)
      z(n,i) = 0.0e0
520   continue

      z(n,n) = 1.0e0
      e(1) = 0.0e0
      return
end subroutine tred2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine imtql2(nm,n,d,e,z,ierr)
implicit none
integer i,j,k,l,m,n,nm,ierr
real d(n),e(n),z(nm,n)
real b,c,f,g,p,r,s,tst1,tst2

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

!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     Questions and comments should be directed to Alan K. Cline,
!     Pleasant Valley Software, 8603 Altus Cove, Austin, TX 78759.
!     Electronic mail to cline@cs.utexas.edu.

!     this version dated january 1989. (for the IBM 3090vf)

!     ------------------------------------------------------------------

!?      call xuflow(0)
      ierr = 0
      if (n .eq. 1) go to 1001

      do 100 i = 2, n
100      e(i-1) = e(i)

      e(n) = 0.0e0

      do 240 l = 1, n
         j = 0
!     .......... look for small sub-diagonal element ..........
105      do 110 m = l, n-1
         tst1 = abs(d(m)) + abs(d(m+1))
         tst2 = tst1 + abs(e(m))
         if (tst2 .eq. tst1) go to 120
110      continue

120      p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
!     .......... form shift ..........
         g = (d(l+1) - p) / (2.0e0 * e(l))
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c *      r = pythag(g,1.0d0)
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (abs(g).le.1.0e0) then
            r = sqrt(1.0e0 + g*g)
         else
            r = g * sqrt(1.0e0 + (1.0e0/g)**2)
         endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
         g = d(m) - p + e(l) / (g + sign(r,g))
         s = 1.0e0
         c = 1.0e0
         p = 0.0e0
!c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 i = m-1, l, -1
            f = s * e(i)
            b = c * e(i)
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
!c *         r = pythag(f,g)
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
            if (abs(f).ge.abs(g)) then
               r = abs(f) * sqrt(1.0e0 + (g/f)**2)
            else if (g .ne. 0.0e0) then
               r = abs(g) * sqrt((f/g)**2 + 1.0e0)
            else
               r = abs(f)
            endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccc
            e(i+1) = r
            if (r .eq. 0.0e0) then
!c     .......... recover from underflow ..........
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
!c     .......... form vector ..........
            do 180 k = 1, n
               f = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * f
               z(k,i) = c * z(k,i) - s * f
180            continue
!c
200            continue
!c
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0e0
         go to 105
240      continue
!c     .......... order eigenvalues and eigenvectors ..........
      do 300 i = 1, n-1
         k = i
         p = d(i)
!c
         do 260 j = i+1, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
260         continue
!c
         d(k) = d(i)
         d(i) = p
!c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
280         continue
!c
300         continue
!c
      go to 1001
!c     .......... set error -- no convergence to an
!c                eigenvalue after 30 iterations ..........
1000  ierr = l
1001  return
end subroutine imtql2
