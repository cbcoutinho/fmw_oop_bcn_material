module blas_mod
! Implementation or interfaces to BLAS functions
  use types_mod
  implicit none
  private 
  public :: idamax

contains

  function idamax ( n, dx, incx )

    !*****************************************************************************80
    !
    !! IDAMAX indexes the array element of maximum absolute value.
    !
    !  Discussion:
    !
    !    This routine uses double precision real arithmetic.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    08 April 1999
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Charles Lawson, Richard Hanson, 
    !    David Kincaid, Fred Krogh.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
    !    LINPACK User's Guide,
    !    SIAM, 1979,
    !    ISBN13: 978-0-898711-72-1,
    !    LC: QA214.L56.
    !
    !    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
    !    Algorithm 539, 
    !    Basic Linear Algebra Subprograms for Fortran Usage,
    !    ACM Transactions on Mathematical Software,
    !    Volume 5, Number 3, September 1979, pages 308-323.
    !
    !  Parameters:
    !
    !    Input, integer(ip) :: N, the number of entries in the vector.
    !
    !    Input, real(RP) DX(*), the vector to be examined.
    !
    !    Input, integer(ip) :: INCX, the increment between successive 
    !    entries of SX.
    !
    !    Output, integer(ip) :: IDAMAX, the index of the element of SX of 
    !    maximum absolute value.
    !
    implicit none

    integer(IP), intent(in) :: n
    real(RP)   , intent(in) :: dx(*)
    integer(IP), intent(in) :: incx
    integer(IP)             :: idamax
    
    real(RP) :: dmax
    integer(IP) :: i
    integer(IP) :: ix

    idamax = 0

    if ( n < 1 .or. incx <= 0 ) then
       return
    end if

    idamax = 1

    if ( n == 1 ) then
       return
    end if

    if ( incx == 1 ) then

       dmax = abs ( dx(1) )

       do i = 2, n
          if ( dmax < abs ( dx(i) ) ) then
             idamax = i
             dmax = abs ( dx(i) )
          end if
       end do

    else

       ix = 1
       dmax = abs ( dx(1) )
       ix = ix + incx

       do i = 2, n
          if ( dmax < abs ( dx(ix) ) ) then
             idamax = i
             dmax = abs ( dx(ix) )
          end if
          ix = ix + incx
       end do

    end if

    return
  end function idamax

end module blas_mod
