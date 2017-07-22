module direct_solver_mod
  use types_mod
  implicit none
  private
  public :: dgbfa, dgbsl, dgefa, dgesl, dpbufa, dpbusl
contains

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
  subroutine dgefa ( a, lda, n, ipvt, info )

    !*****************************************************************************80
    !
    !! DGEFA factors a real general matrix.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 March 2001
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
    !    Input/output, real(RP) A(LDA,N).
    !    On intput, the matrix to be factored.
    !    On output, an upper triangular matrix and the multipliers used to obtain
    !    it.  The factorization can be written A=L*U, where L is a product of
    !    permutation and unit lower triangular matrices, and U is upper triangular.
    !
    !    Input, integer(ip) :: LDA, the leading dimension of A.
    !
    !    Input, integer(ip) :: N, the order of the matrix A.
    !
    !    Output, integer(ip) :: IPVT(N), the pivot indices.
    !
    !    Output, integer(ip) :: INFO, singularity indicator.
    !    0, normal value.
    !    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
    !    but it does indicate that DGESL or DGEDI will divide by zero if called.
    !    Use RCOND in DGECO for a reliable indication of singularity.
    !
    implicit none
    real(RP)    , intent(inout) :: a(lda,n)
    integer(IP) , intent(in)    :: lda
    integer(IP) , intent(in)    :: n
    integer(IP) , intent(out)   :: ipvt(n)
    integer(IP) , intent(out)   :: info
    
    integer(IP) :: j
    integer(IP) :: k
    integer(IP) :: l
    real(RP)    :: t
    !
    !  Gaussian elimination with partial pivoting.
    !
    info = 0

    do k = 1, n - 1
       !
       !  Find L = pivot index.
       !
       l = idamax ( n-k+1, a(k,k), 1 ) + k - 1
       ipvt(k) = l
       !
       !  Zero pivot implies this column already triangularized.
       !
       if ( a(l,k) == 0.0D+00 ) then
          info = k
          cycle
       end if
       !
       !  Interchange if necessary.
       !
       if ( l /= k ) then
          t = a(l,k)
          a(l,k) = a(k,k)
          a(k,k) = t
       end if
       !
       !  Compute multipliers.
       !
       t = -1.0D+00 / a(k,k)
       a(k+1:n,k) = a(k+1:n,k) * t
       !
       !  Row elimination with column indexing.
       !
       do j = k+1, n
          t = a(l,j)
          if ( l /= k ) then
             a(l,j) = a(k,j)
             a(k,j) = t
          end if
          a(k+1:n,j) = a(k+1:n,j) + t * a(k+1:n,k)
       end do

    end do

    ipvt(n) = n

    if ( a(n,n) == 0.0D+00 ) then
       info = n
    end if

    return
  end subroutine dgefa
  subroutine dgesl ( a, lda, n, ipvt, b, job )
    !*****************************************************************************80
    !
    !! DGESL solves a real general linear system A * X = B.
    !
    !  Discussion:
    !
    !    DGESL can solve either of the systems A * X = B or A' * X = B.
    !
    !    The system matrix must have been factored by DGECO or DGEFA.
    !
    !    A division by zero will occur if the input factor contains a
    !    zero on the diagonal.  Technically this indicates singularity
    !    but it is often caused by improper arguments or improper
    !    setting of LDA.  It will not occur if the subroutines are
    !    called correctly and if DGECO has set 0.0 < RCOND
    !    or DGEFA has set INFO == 0.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 March 2001
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
    !    Input, real(RP) A(LDA,N), the output from DGECO or DGEFA.
    !
    !    Input, integer(ip) :: LDA, the leading dimension of A.
    !
    !    Input, integer(ip) :: N, the order of the matrix A.
    !
    !    Input, integer(ip) :: IPVT(N), the pivot vector from DGECO or DGEFA.
    !
    !    Input/output, real(RP) B(N).
    !    On input, the right hand side vector.
    !    On output, the solution vector.
    !
    !    Input, integer(ip) :: JOB.
    !    0, solve A * X = B;
    !    nonzero, solve A' * X = B.
    !
    implicit none
    real(RP)   , intent(in)    :: a(lda,n)
    integer(IP), intent(in)    :: lda
    integer(IP), intent(in)    :: n
    integer(IP), intent(in)    :: ipvt(n)
    real(RP)   , intent(inout) :: b(n)
    integer(IP), intent(in)    :: job
    
    integer(IP) :: k
    integer(IP) :: l
    real(RP)    :: t
    !
    !  Solve A * X = B.
    !
    if ( job == 0 ) then

       do k = 1, n-1

          l = ipvt(k)
          t = b(l)

          if ( l /= k ) then
             b(l) = b(k)
             b(k) = t
          end if

          b(k+1:n) = b(k+1:n) + t * a(k+1:n,k)

       end do

       do k = n, 1, -1
          b(k) = b(k) / a(k,k)
          t = -b(k)
          b(1:k-1) = b(1:k-1) + t * a(1:k-1,k)
       end do

    else
       !
       !  Solve A' * X = B.
       !
       do k = 1, n
          t = dot_product ( a(1:k-1,k), b(1:k-1) )
          b(k) = ( b(k) - t ) / a(k,k)
       end do

       do k = n-1, 1, -1

          b(k) = b(k) + dot_product ( a(k+1:n,k), b(k+1:n) )
          l = ipvt(k)

          if ( l /= k ) then
             t = b(l)
             b(l) = b(k)
             b(k) = t
          end if

       end do

    end if

    return
  end subroutine dgesl
  
  subroutine  dpbufa ( n, mu, a, info )
!*****************************************************************************80
!
!! PBU_FA factors an PBU matrix.
!
!  Discussion:
!
!    The PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!    The matrix A must be a positive definite symmetric band matrix.
!
!    Once factored, linear systems A*x=b involving the matrix can be solved
!    by calling PBU_SL.  No pivoting is performed.  Pivoting is not necessary
!    for positive definite symmetric matrices.  If the matrix is not positive
!    definite, the algorithm may behave correctly, but it is also possible
!    that an illegal divide by zero will occur.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 1998
!
!  Author:
!
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
!    FORTRAN90 version by John Burkardt
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
!    Input, integer(IP) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer(IP) MU, the number of superdiagonals of the matrix.
!    MU must be at least 0, and no more than N-1.
!
!    Input/output, real(RP) A(MU+1,N), the N by N matrix, stored 
!    in LINPACK positive definite symmetric band matrix storage.
!    On output, A contains information describing a factored form
!    of the matrix, that can be used to solve linear systems
!    A*x=b, using PBU_SL.
!
!    Output, integer(IP) INFO, singularity flag.
!    0, the matrix is nonsingular.
!    nonzero, the matrix is singular.
!
  implicit none

  integer(IP), intent(in)    :: n
  integer(IP), intent(in)    :: mu
  real(RP)   , intent(inout) :: a(mu+1,n)
  integer(IP) :: info
  
  integer(IP) :: ik
  
  integer(IP) :: j
  integer(IP) :: jk
  integer(IP) :: k
  integer(IP) :: mm
  real(RP) :: s

  info = 0

  do j = 1, n

    ik = mu + 1
    jk = max ( j - mu, 1 )
    mm = max ( mu + 2 - j, 1 )

    s = 0.0D+00

    do k = mm, mu

      a(k,j) = ( a(k,j) - sum ( a(ik:ik+k-mm-1,jk) * a(mm:k-1,j) ) ) &
        / a(mu+1,jk)

      s = s + a(k,j) * a(k,j)

      ik = ik - 1
      jk = jk + 1

    end do

    s = a(mu+1,j) - s

    if ( s <= 0.0D+00 ) then
      info = j
      return
    end if

    a(mu+1,j) = sqrt ( s )

  end do

  return
end subroutine dpbufa
  
  subroutine dpbusl ( n, mu, a_lu, b )
!*****************************************************************************80
!
!! PBU_SL solves an PBU system factored by PBU_FA.
!
!  Discussion:
!
!    The PBU storage format is for a symmetric positive definite band matrix.
!
!    To save storage, only the diagonal and upper triangle of A is stored,
!    in a compact diagonal format that preserves columns.
!
!    The diagonal is stored in row MU+1 of the array.
!    The first superdiagonal in row MU, columns 2 through N.
!    The second superdiagonal in row MU-1, columns 3 through N.
!    The MU-th superdiagonal in row 1, columns MU+1 through N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 October 1998
!
!  Author:
!
!    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
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
!    Input, integer(IP) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer(IP) MU, the number of superdiagonals of the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real(RP) A_LU(MU+1,N), the LU factors from PBU_FA.
!
!    Input/output, real(RP) B(N).
!    On input, B contains the right hand side of the linear system
!    to be solved.
!    On output, B contains X, the solution vector.
!
  implicit none

  integer(IP), intent(in)    :: n
  integer(IP), intent(in)    :: mu
  real(RP)   , intent(in)    :: a_lu(mu+1,n)
  real(RP)   , intent(inout) :: b(n)
  
  integer(IP) :: i
  integer(IP) :: ilo
  integer(IP) :: k
!
!  Solve L * Y = B.
!
  do k = 1, n
    ilo = max ( 1, k - mu )
    b(k) = ( b(k) - sum ( b(ilo:k-1) * a_lu(mu+1+ilo-k:mu,k) ) ) &
      / a_lu(mu+1,k)
  end do
!
!  Solve U * X = Y.
!
  do k = n, 1, -1

    b(k) = b(k) / a_lu(mu+1,k)

    ilo = max ( 1, k - mu )
    do i = ilo, k - 1
      b(i) = b(i) - b(k) * a_lu(mu+1+i-k,k)
    end do

  end do

  return
end subroutine dpbusl
  
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

  
  
end module direct_solver_mod
