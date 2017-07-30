#include "mcheck.i90"
module dense_matrix_mod
  use types_mod,  only: IP, RP
  use blas_mod,   only: idamax
  use matrix_mod, only: matrix_t
  implicit none
  private

  public :: dense_matrix_t
  type, extends(matrix_t) :: dense_matrix_t
     private
     real(rp), allocatable :: a(:,:)
   contains
     procedure :: create    => dense_matrix_create
     procedure :: assembly  => dense_matrix_assembly
     procedure :: apply     => dense_matrix_apply
     procedure :: factorize => dense_matrix_factorize
     procedure :: backsolve => dense_matrix_backsolve
     procedure :: free      => dense_matrix_free
  end type dense_matrix_t

contains

  subroutine dense_matrix_create(this,n,ml,mu,nz)
    class(dense_matrix_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: n
    integer(ip), optional, intent(in)    :: ml,mu,nz
    call this%free()
    call this%set_size(n)
    allocate(this%a(n,n))
    this%a = 0.0_rp
  end subroutine dense_matrix_create

  subroutine dense_matrix_assembly(this,i,j,a)
    implicit none
    class(dense_matrix_t), intent(inout) :: this
    integer(ip)    , intent(in)    :: i,j
    real(rp)       , intent(in)    :: a
    this%a(i,j) = this%a(i,j) + a
  end subroutine dense_matrix_assembly

  subroutine dense_matrix_apply(this,x,y)
    implicit none
    class(dense_matrix_t), intent(in) :: this
    real(rp)       , intent(in)       :: x(:)
    real(rp)       , intent(inout)    :: y(:)
    call mv_ge ( this%get_size(), this%get_size(), this%a, x, y )
  end subroutine dense_matrix_apply

  subroutine dense_matrix_factorize(this, factors, pivots)
    implicit none
    class(dense_matrix_t)       , intent(in)    :: this
    class(matrix_t), allocatable, intent(inout) :: factors
    integer(ip)                 , intent(inout) :: pivots(:)
    integer(ip) :: info

    if (allocated(factors)) then
      call factors%free()
      deallocate(factors)
    end if

    allocate(dense_matrix_t :: factors)
    select type(factors)
    type is(dense_matrix_t)
       call factors%create(this%get_size())
       factors%a = this%a
       call dgefa ( factors%a, this%get_size(), this%get_size(), pivots, info )
       mcheck(info==0,'Error in dense factorization')
    class default
    end select
  end subroutine dense_matrix_factorize

  subroutine dense_matrix_backsolve(this,pivots,rhs,x)
    class(dense_matrix_t), intent(in)    :: this
    integer(ip)          , intent(in)    :: pivots(:)
    real(rp)             , intent(in)    :: rhs(:)
    real(rp)             , intent(inout) :: x(:)
    x = rhs
    call dgesl ( this%a, this%get_size(), this%get_size(), pivots, x, 0 )
  end subroutine dense_matrix_backsolve

  subroutine dense_matrix_free(this)
    class(dense_matrix_t), intent(inout) :: this
    call this%set_size(-1)
    if (allocated(this%a)) deallocate(this%a)
  end subroutine dense_matrix_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Legacy code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mv_ge ( m, n, a, x, b )

    !*****************************************************************************80
    !
    !! MV_GE multiplies an R8GE matrix by an R8VEC.
    !
    !  Discussion:
    !
    !    The R8GE storage format is used for a general M by N matrix.  A storage
    !    space is made for each entry.  The two dimensional logical
    !    array can be thought of as a vector of M*N entries, starting with
    !    the M entries in the column 1, then the M entries in column 2
    !    and so on.  Considered as a vector, the entry A(I,J) is then stored
    !    in vector location I+(J-1)*M.
    !
    !    R8GE storage is used by LINPACK and LAPACK.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 January 1999
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
    !    Input, real(RP) A(M,N), the R8GE matrix.
    !
    !    Input, real(RP) X(N), the vector to be multiplied by A.
    !
    !    Output, real(RP) B(M), the product A * x.
    !
    implicit none
    integer(IP), intent(in) :: m
    integer(IP), intent(in) :: n
    real(RP)   , intent(in)  :: a(m,n)
    real(RP)   , intent(in)  :: x(n)
    real(RP)   , intent(out) :: b(m)
    b(1:m) = matmul ( a(1:m,1:n), x(1:n) )
    return
  end subroutine mv_ge

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
       if ( a(l,k) == 0.0_RP ) then
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
       t = -1.0_RP / a(k,k)
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

    if ( a(n,n) == 0.0_RP ) then
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

end module dense_matrix_mod
