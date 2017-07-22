#include "mcheck.i90"
module band_matrix_mod
  use types_mod
  use blas_mod
  use matrix_mod
  implicit none
  private

  type, extends(matrix_t) :: band_matrix_t
     private
     real(rp), allocatable :: a(:,:)
     integer(ip) :: ml  = -1
     integer(ip) :: mu  = -1
     integer(ip) :: lda = -1
   contains
     procedure :: create    => band_matrix_create
     procedure :: assembly  => band_matrix_assembly
     procedure :: apply     => band_matrix_apply
     procedure :: factorize => band_matrix_factorize
     procedure :: backsolve => band_matrix_backsolve
     procedure :: free      => band_matrix_free
  end type band_matrix_t
  public :: band_matrix_t

contains
  
  subroutine band_matrix_create(this,n,ml,mu,nz)
    class(band_matrix_t) , intent(inout) :: this
    integer(ip)          , intent(in)    :: n
    integer(ip), optional, intent(in)    :: ml,mu,nz
    call this%free()
    call this%set_size(n)
    mcheck(present(ml),"ml dummy argument required by band_matrix_create")
    mcheck(present(mu),"mu dummy argument required by band_matrix_create")
    this%ml = ml
    this%mu = mu
    this%lda = 2*ml+mu+1
    allocate ( this%a(1:this%lda,1:n) )
    this%a = 0.0_rp
  end subroutine band_matrix_create

  subroutine band_matrix_assembly(this,i,j,a) 
    implicit none
    class(band_matrix_t), intent(inout) :: this
    integer(ip)         , intent(in)    :: i,j
    real(rp)            , intent(in)    :: a
    this%a(i-j+this%ml+this%mu+1,j) = this%a(i-j+this%ml+this%mu+1,j) + a
  end subroutine band_matrix_assembly

  subroutine band_matrix_apply(this,x,y) 
    implicit none
    class(band_matrix_t), intent(in)    :: this
    real(rp)            , intent(in)    :: x(:)
    real(rp)            , intent(inout) :: y(:)
    call mv_gb ( this%get_size(), this%get_size(), this%ml, this%mu, this%a, x, y )
  end subroutine band_matrix_apply

  subroutine band_matrix_factorize(this, factors, pivots) 
    implicit none
    class(band_matrix_t)        , intent(in)    :: this
    class(matrix_t), allocatable, intent(inout) :: factors
    integer(ip)                 , intent(inout) :: pivots(:)
    integer(ip) :: info, lda
    
    if (allocated(factors)) then
      call factors%free()
      deallocate(factors)
    end if
    
    allocate(band_matrix_t :: factors)
    select type(factors)
    type is(band_matrix_t)
       call factors%create(this%get_size(),this%ml,this%mu)
       factors%a = this%a
       call dgbfa ( factors%a, this%lda, this%get_size(), this%ml, this%mu, pivots, info )
       mcheck(info==0, 'Error in band factorization')
    class default
    end select
  end subroutine band_matrix_factorize

  subroutine band_matrix_backsolve(this,pivots,rhs,x) 
    class(band_matrix_t), intent(in)    :: this
    integer(ip)         , intent(in)    :: pivots(:)
    real(rp)            , intent(in)    :: rhs(:)
    real(rp)            , intent(inout) :: x(:)
    x = rhs
    call dgbsl ( this%a, this%lda, this%get_size(), this%ml, this%mu, pivots, x, 0 )
  end subroutine band_matrix_backsolve

  subroutine band_matrix_free(this)
    class(band_matrix_t), intent(inout) :: this
    call this%set_size(-1)
    this%ml = -1
    this%mu = -1
    this%lda = -1
    if (allocated(this%a)) deallocate(this%a)
  end subroutine band_matrix_free
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Legacy code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mv_gb ( m, n, ml, mu, a, x, b )

    !*****************************************************************************80
    !
    !! MV_GB multiplies a banded matrix by an R8VEC.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    04 June 2014
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer(ip) :: M, the number of rows of the matrix.
    !    M must be positive.
    !
    !    Input, integer(ip) :: N, the number of columns of the matrix.
    !    N must be positive.
    !
    !    Input, integer(ip) :: ML, MU, the lower and upper bandwidths.
    !
    !    Input, real(RP) A(2*ML+MU+1,N), the matrix.
    !
    !    Input, real(RP) X(N), the vector to be multiplied by A.
    !
    !    Output, real(RP) B(M), the product A * x.
    !
    implicit none

    integer(IP), intent(in)  :: m
    integer(IP), intent(in)  :: n
    integer(IP), intent(in)  :: ml
    integer(IP), intent(in)  :: mu
    real(RP)   , intent(in)  :: a(2*ml+mu+1,n)
    real(RP)   , intent(in)  :: x(n)
    real(RP)   , intent(out) :: b(m)
    
    integer(IP) :: i
    integer(IP) :: j
    integer(IP) :: jhi
    integer(IP) :: jlo
    
    b(1:m) = 0.0D+00

    do i = 1, n
       jlo = max ( 1, i - ml )
       jhi = min ( n, i + mu )
       do j = jlo, jhi
          b(i) = b(i) + a(i-j+ml+mu+1,j) * x(j)
       end do
    end do

    return
  end subroutine mv_gb

  subroutine dgbfa ( abd, lda, n, ml, mu, ipvt, info )

    !*****************************************************************************80
    !
    !! DGBFA factors a real band matrix by elimination.
    !
    !  Discussion:
    !
    !    DGBFA is usually called by DGBCO, but it can be called
    !    directly with a saving in time if RCOND is not needed.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 May 2005
    !
    !  Author:
    !
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
    !  Parameters:
    !
    !    Input/output, real(RP) ABD(LDA,N).  On input, the matrix in band
    !    storage.  The columns of the matrix are stored in the columns of ABD
    !    and the diagonals of the matrix are stored in rows ML+1 through
    !    2*ML+MU+1 of ABD.  On output, an upper triangular matrix in band storage
    !    and the multipliers which were used to obtain it.  The factorization
    !    can be written A = L*U where L is a product of permutation and unit lower
    !    triangular matrices and U is upper triangular.
    !
    !    Input, integer(ip) :: LDA, the leading dimension of the array ABD.
    !    2*ML + MU + 1 <= LDA is required.
    !
    !    Input, integer(ip) :: N, the order of the matrix.
    !
    !    Input, integer(ip) :: ML, MU, the number of diagonals below and above
    !    the main diagonal.  0 <= ML < N, 0 <= MU < N.
    !
    !    Output, integer(ip) :: IPVT(N), the pivot indices.
    !
    !    Output, integer(ip) :: INFO, error flag.
    !    0, normal value.
    !    K, if U(K,K) == 0.0D+00.  This is not an error condition for this
    !      subroutine, but it does indicate that DGBSL will divide by zero if
    !      called.  Use RCOND in DGBCO for a reliable indication of singularity.
    !
    implicit none

    real(RP)   , intent(inout) :: abd(lda,n)
    integer(IP), intent(in)    :: lda
    integer(IP), intent(in)    :: n
    integer(IP), intent(out)   :: ipvt(n)
    integer(IP), intent(out)   :: info
    
    integer(IP) :: i
    integer(IP) :: i0    
    integer(IP) :: j
    integer(IP) :: j0
    integer(IP) :: j1
    integer(IP) :: ju
    integer(IP) :: jz
    integer(IP) :: k
    integer(IP) :: l
    integer(IP) :: lm
    integer(IP) :: m
    integer(IP) :: ml
    integer(IP) :: mm
    integer(IP) :: mu
    real(RP) t

    m = ml + mu + 1
    info = 0
    !
    !  Zero initial fill-in columns.
    !
    j0 = mu + 2
    j1 = min ( n, m ) - 1

    do jz = j0, j1
       i0 = m + 1 - jz
       do i = i0, ml
          abd(i,jz) = 0.0D+00
       end do
    end do

    jz = j1
    ju = 0
    !
    !  Gaussian elimination with partial pivoting.
    !
    do k = 1, n-1
       !
       !  Zero out the next fill-in column.
       !
       jz = jz + 1
       if ( jz <= n ) then
          abd(1:ml,jz) = 0.0D+00
       end if
       !
       !  Find L = pivot index.
       !
       lm = min ( ml, n-k )
       l = idamax ( lm+1, abd(m,k), 1 ) + m - 1
       ipvt(k) = l + k - m
       !
       !  Zero pivot implies this column already triangularized.
       !
       if ( abd(l,k) == 0.0D+00 ) then

          info = k
          !
          !  Interchange if necessary.
          !
       else

          if ( l /= m ) then
             t = abd(l,k)
             abd(l,k) = abd(m,k)
             abd(m,k) = t
          end if
          !
          !  Compute multipliers.
          !
          t = -1.0D+00 / abd(m,k)
          call dscal ( lm, t, abd(m+1,k), 1 )
          !
          !  Row elimination with column indexing.
          !
          ju = min ( max ( ju, mu + ipvt(k) ), n )
          mm = m

          do j = k+1, ju
             l = l - 1
             mm = mm - 1
             t = abd(l,j)
             if ( l /= mm ) then
                abd(l,j) = abd(mm,j)
                abd(mm,j) = t
             end if
             call daxpy ( lm, t, abd(m+1,k), 1, abd(mm+1,j), 1 )
          end do

       end if

    end do

    ipvt(n) = n

    if ( abd(m,n) == 0.0D+00 ) then
       info = n
    end if

    return
  end subroutine dgbfa

    subroutine dgbsl ( abd, lda, n, ml, mu, ipvt, b, job )

    !*****************************************************************************80
    !
    !! DGBSL solves a real banded system factored by DGBCO or DGBFA.
    !
    !  Discussion:
    !
    !    DGBSL can solve either A * X = B  or  A' * X = B.
    !
    !    A division by zero will occur if the input factor contains a
    !    zero on the diagonal.  Technically this indicates singularity
    !    but it is often caused by improper arguments or improper
    !    setting of LDA.  It will not occur if the subroutines are
    !    called correctly and if DGBCO has set 0.0 < RCOND
    !    or DGBFA has set INFO == 0.
    !
    !    To compute inverse(A) * C  where C is a matrix with P columns:
    !
    !      call dgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )
    !
    !      if ( rcond is too small ) then
    !        exit
    !      end if
    !
    !      do j = 1, p
    !        call dgbsl ( abd, lda, n, ml, mu, ipvt, c(1,j), 0 )
    !      end do
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    17 May 2005
    !
    !  Author:
    !
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
    !  Parameters:
    !
    !    Input, real(RP) ABD(LDA,N), the output from DGBCO or DGBFA.
    !
    !    Input, integer(ip) :: LDA, the leading dimension of the array ABD.
    !
    !    Input, integer(ip) :: N, the order of the matrix.
    !
    !    Input, integer(ip) :: ML, MU, the number of diagonals below and above
    !    the main diagonal.  0 <= ML < N, 0 <= MU < N.
    !
    !    Input, integer(ip) :: IPVT(N), the pivot vector from DGBCO or DGBFA.
    !
    !    Input/output, real(RP) B(N).  On input, the right hand side.
    !    On output, the solution.
    !
    !    Input, integer(ip) :: JOB, job choice.
    !    0, solve A*X=B.
    !    nonzero, solve A'*X=B.
    !
    implicit none

    real(RP)   , intent(in)     :: abd(lda,n)
    integer(IP), intent(in)    :: lda
    integer(IP), intent(in)    :: n
    integer(IP), intent(in)    :: ml
    integer(IP), intent(in)    :: mu
    integer(IP), intent(in)    :: ipvt(n)
    real(RP)   , intent(inout) :: b(n)
    integer(IP), intent(in)    :: job
    
    integer(IP) :: k
    integer(IP) :: l
    integer(IP) :: la
    integer(IP) :: lb
    integer(IP) :: lm
    integer(IP) :: m
    
    real(RP) :: t

    m = mu + ml + 1
    !
    !  JOB = 0, Solve A * x = b.
    !
    !  First solve L * y = b.
    !
    if ( job == 0 ) then

       if ( 0 < ml ) then

          do k = 1, n-1
             lm = min ( ml, n-k )
             l = ipvt(k)
             t = b(l)
             if ( l /= k ) then
                b(l) = b(k)
                b(k) = t
             end if
             call daxpy ( lm, t, abd(m+1,k), 1, b(k+1), 1 )
          end do

       end if
       !
       !  Now solve U * x = y.
       !
       do k = n, 1, -1
          b(k) = b(k) / abd(m,k)
          lm = min ( k, m ) - 1
          la = m - lm
          lb = k - lm
          t = -b(k)
          call daxpy ( lm, t, abd(la,k), 1, b(lb), 1 )
       end do
       !
       !  JOB nonzero, solve A' * x = b.
       !
       !  First solve U' * y = b.
       !
    else

       do k = 1, n
          lm = min ( k, m ) - 1
          la = m - lm
          lb = k - lm
          t = ddot ( lm, abd(la,k), 1, b(lb), 1 )
          b(k) = ( b(k) - t ) / abd(m,k)
       end do
       !
       !  Now solve L' * x = y.
       !
       if ( 0 < ml ) then

          do k = n - 1, 1, -1
             lm = min ( ml, n - k )
             b(k) = b(k) + ddot ( lm, abd(m+1,k), 1, b(k+1), 1 )
             l = ipvt(k)
             if ( l /= k ) then
                t = b(l)
                b(l) = b(k)
                b(k) = t
             end if
          end do

       end if

    end if

    return
  end subroutine dgbsl
  subroutine daxpy ( n, da, dx, incx, dy, incy )
    !*****************************************************************************80
    !
    !! DAXPY computes constant times a vector plus a vector.
    !
    !  Discussion:
    !
    !    This routine uses double precision real arithmetic.
    !
    !    This routine uses unrolled loops for increments equal to one.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    16 May 2005
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
    !    Input, integer(ip) :: N, the number of elements in DX and DY.
    !
    !    Input, real(RP) DA, the multiplier of DX.
    !
    !    Input, real(RP) DX(*), the first vector.
    !
    !    Input, integer(ip) :: INCX, the increment between successive 
    !    entries of DX.
    !
    !    Input/output, real(RP) DY(*), the second vector.
    !    On output, DY(*) has been replaced by DY(*) + DA * DX(*).
    !
    !    Input, integer(ip) :: INCY, the increment between successive 
    !    entries of DY.
    !
    implicit none
    integer(IP), intent(in)    :: n
    real(RP)   , intent(in)    :: da
    real(RP)   , intent(in)    :: dx(*)
    integer(IP), intent(in)    :: incx
    real(RP)   , intent(inout) :: dy(*)
    integer(IP), intent(in)    :: incy
    
    integer(IP) :: i
    integer(IP) :: ix
    integer(IP) :: iy
    integer(IP) :: m

    if ( n <= 0 ) then
       return
    end if

    if ( da == 0.0D+00 ) then
       return
    end if
    !
    !  Code for unequal increments or equal increments
    !  not equal to 1.
    !
    if ( incx /= 1 .or. incy /= 1 ) then

       if ( 0 <= incx ) then
          ix = 1
       else
          ix = ( - n + 1 ) * incx + 1
       end if

       if ( 0 <= incy ) then
          iy = 1
       else
          iy = ( - n + 1 ) * incy + 1
       end if

       do i = 1, n
          dy(iy) = dy(iy) + da * dx(ix)
          ix = ix + incx
          iy = iy + incy
       end do
       !
       !  Code for both increments equal to 1.
       !
    else

       m = mod ( n, 4 )

       dy(1:m) = dy(1:m) + da * dx(1:m)

       do i = m + 1, n, 4
          dy(i  ) = dy(i  ) + da * dx(i  )
          dy(i+1) = dy(i+1) + da * dx(i+1)
          dy(i+2) = dy(i+2) + da * dx(i+2)
          dy(i+3) = dy(i+3) + da * dx(i+3)
       end do

    end if

    return
  end subroutine daxpy
  function ddot ( n, dx, incx, dy, incy )
    !*****************************************************************************80
    !
    !! DDOT forms the dot product of two vectors.
    !
    !  Discussion:
    !
    !    This routine uses double precision real arithmetic.
    !
    !    This routine uses unrolled loops for increments equal to one.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    16 May 2005
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
    !    Input, integer(ip) :: N, the number of entries in the vectors.
    !
    !    Input, real(RP) DX(*), the first vector.
    !
    !    Input, integer(ip) :: INCX, the increment between successive 
    !    entries in DX.
    !
    !    Input, real(RP) DY(*), the second vector.
    !
    !    Input, integer(ip) :: INCY, the increment between successive 
    !    entries in DY.
    !
    !    Output, real(RP) DDOT, the sum of the product of the 
    !    corresponding entries of DX and DY.
    !
    implicit none

    integer(IP), intent(in) :: n
    real(RP)   , intent(in) :: dx(*)
    integer(IP), intent(in) :: incx
    real(RP)   , intent(in) :: dy(*)
    integer(IP), intent(in) :: incy
    real(RP)                :: ddot
    
    real(RP) :: dtemp
    integer(IP) :: i
    integer(IP) :: ix
    integer(IP) :: iy
    integer(IP) :: m
    
    ddot = 0.0D+00
    dtemp = 0.0D+00

    if ( n <= 0 ) then
       return
    end if
    !
    !  Code for unequal increments or equal increments
    !  not equal to 1.
    !
    if ( incx /= 1 .or. incy /= 1 ) then

       if ( 0 <= incx ) then
          ix = 1
       else
          ix = ( - n + 1 ) * incx + 1
       end if

       if ( 0 <= incy ) then
          iy = 1
       else
          iy = ( - n + 1 ) * incy + 1
       end if

       do i = 1, n
          dtemp = dtemp + dx(ix) * dy(iy)
          ix = ix + incx
          iy = iy + incy
       end do
       !
       !  Code for both increments equal to 1.
       !
    else

       m = mod ( n, 5 )

       do i = 1, m
          dtemp = dtemp + dx(i) * dy(i)
       end do

       do i = m + 1, n, 5

          dtemp = dtemp + dx(i  ) * dy(i  ) &
               + dx(i+1) * dy(i+1) &
               + dx(i+2) * dy(i+2) &
               + dx(i+3) * dy(i+3) &
               + dx(i+4) * dy(i+4)
       end do

    end if

    ddot = dtemp

    return
  end function ddot
  subroutine dscal ( n, sa, x, incx )
    !*****************************************************************************80
    !
    !! DSCAL scales a vector by a constant.
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
    !    Input, real(RP) SA, the multiplier.
    !
    !    Input/output, real(RP) X(*), the vector to be scaled.
    !
    !    Input, integer(ip) :: INCX, the increment between successive 
    !    entries of X.
    !
    implicit none
    integer(IP), intent(in)    :: n
    real(RP)   , intent(in)    :: sa
    real(RP)   , intent(inout) :: x(*)
    integer(IP), intent(in)    :: incx
    
    integer(IP) :: i
    integer(IP) :: ix
    integer(IP) :: m
    
    if ( n <= 0 ) then

    else if ( incx == 1 ) then

       m = mod ( n, 5 )

       x(1:m) = sa * x(1:m)

       do i = m+1, n, 5
          x(i)   = sa * x(i)
          x(i+1) = sa * x(i+1)
          x(i+2) = sa * x(i+2)
          x(i+3) = sa * x(i+3)
          x(i+4) = sa * x(i+4)
       end do

    else

       if ( 0 <= incx ) then
          ix = 1
       else
          ix = ( - n + 1 ) * incx + 1
       end if

       do i = 1, n
          x(ix) = sa * x(ix)
          ix = ix + incx
       end do

    end if

    return
  end subroutine dscal

end module band_matrix_mod
