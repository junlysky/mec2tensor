Subroutine TRPL2V(TRPL,XYZ)
implicit none
! Transforms to XYZ components of a unit vector from
! the trend and plunge for the vector.
! trend is azimuth ( clockwise from north looking down)
! plunge is the downward dip measured from the horizontal
! all angles in radians
! X is north，Y is east,Z is down
! XYZ(1) is x(north),XYZ(2) is y(east),XYZ(3) is z(down)
! trpl(1) is azimuth angle
! trpl(2) is plunge angle

  real XYZ(3),TRPL(2)
  integer j
  XYZ(1) = COS(TRPL(1))*COS(TRPL(2))
  XYZ(2) = SIN(TRPL(1))*COS(TRPL(2))
  XYZ(3) = SIN(TRPL(2))
  do j=1,3
     if (abs(xyz(j)) .lt. 0.0001) xyz(j) = 0.0
     if (abs(abs(xyz(j))-1.0).lt.0.0001) xyz(j)=xyz(j)/abs(xyz(j))
  end do
  return
end subroutine TRPL2V