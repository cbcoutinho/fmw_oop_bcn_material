#include "mcheck.i90"
module symmetric_band_matrix_mod
  use types_mod
  use blas_mod
  use matrix_mod
  implicit none
  private

  type, extends(matrix_t) :: symmetric_band_matrix_t
     private
     real(rp), allocatable :: a(:,:)
     integer(ip) :: mu  = -1
     integer(ip) :: lda = -1
   contains
     procedure :: create    => symmetric_band_matrix_create
     procedure :: assembly  => symmetric_band_matrix_assembly
     procedure :: apply     => symmetric_band_matrix_apply
     procedure :: factorize => symmetric_band_matrix_factorize
     procedure :: backsolve => symmetric_band_matrix_backsolve
     procedure :: free      => symmetric_band_matrix_free
  end type symmetric_band_matrix_t
  public :: symmetric_band_matrix_t

contains
  
  subroutine symmetric_band_matrix_create(this,n,ml,mu,nz)
    class(symmetric_band_matrix_t) , intent(inout) :: this
    integer(ip)          , intent(in)    :: n
    integer(ip), optional, intent(in)    :: ml,mu,nz
    call this%free()
    call this%set_size(n)
    mcheck(present(mu),"mu dummy argument required by symmetric_band_matrix_create")
    this%mu = mu
    this%lda = this%mu+1
    allocate ( this%a(1:this%lda,1:n) )
    this%a = 0.0_rp
  end subroutine symmetric_band_matrix_create

  subroutine symmetric_band_matrix_assembly(this,i,j,a) 
    implicit none
    class(symmetric_band_matrix_t), intent(inout) :: this
    integer(ip)         , intent(in)    :: i,j
    real(rp)            , intent(in)    :: a
    if (j>=i) this%a(i-j+this%mu+1,j) = this%a(i-j+this%mu+1,j) + a
  end subroutine symmetric_band_matrix_assembly

  subroutine symmetric_band_matrix_apply(this,x,y) 
    implicit none
    class(symmetric_band_matrix_t), intent(in)    :: this
    real(rp)            , intent(in)    :: x(:)
    real(rp)            , intent(inout) :: y(:)
    call mv_pbu ( this%get_size(), this%mu, this%a, x, y )
  end subroutine symmetric_band_matrix_apply

  subroutine symmetric_band_matrix_factorize(this, factors, pivots) 
    implicit none
    class(symmetric_band_matrix_t)        , intent(in)    :: this
    class(matrix_t), allocatable, intent(inout) :: factors
    integer(ip)                 , intent(inout) :: pivots(:)
    integer(ip) :: info, lda
    
    if (allocated(factors)) then
      call factors%free()
      deallocate(factors)
    end if
    
    allocate(symmetric_band_matrix_t :: factors)
    select type(factors)
    type is(symmetric_band_matrix_t)
       call factors%create(this%get_size(),mu=this%mu)
       factors%a = this%a
       call dpbufa ( this%get_size(), this%mu, factors%a, info )
       mcheck(info==0, 'Error in band factorization')
    class default
    end select
  end subroutine symmetric_band_matrix_factorize

  subroutine symmetric_band_matrix_backsolve(this,pivots,rhs,x) 
    class(symmetric_band_matrix_t), intent(in)    :: this
    integer(ip)         , intent(in)    :: pivots(:)
    real(rp)            , intent(in)    :: rhs(:)
    real(rp)            , intent(inout) :: x(:)
    x = rhs
    call dpbusl ( this%get_size(), this%mu, this%a, x )
  end subroutine symmetric_band_matrix_backsolve

  subroutine symmetric_band_matrix_free(this)
    class(symmetric_band_matrix_t), intent(inout) :: this
    call this%set_size(-1)
    this%mu = -1
    this%lda = -1
    if (allocated(this%a)) deallocate(this%a)
  end subroutine symmetric_band_matrix_free
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Legacy code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mv_pbu ( n, mu, a, x, b )

!*****************************************************************************80
!
!! PBU_MV multiplies an PBU matrix by an R8VEC.
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
!    15 October 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer(IP) N, the order of the matrix.
!    N must be positive.
!
!    Input, integer(IP) MU, the number of superdiagonals in the matrix.
!    MU must be at least 0 and no more than N-1.
!
!    Input, real(RP) A(MU+1,N), the PBU matrix.
!
!    Input, real(RP) X(N), the vector to be multiplied by A.
!
!    Output, real(RP) B(N), the result vector A * x.
!
  implicit none

  integer(IP), intent(in)    :: n
  integer(IP), intent(in)    :: mu
  real(RP)   , intent(in)    :: a(mu+1,n)
  real(RP)   , intent(in)    :: x(n)
  real(RP)   , intent(inout) :: b(n)
  
  integer(IP) :: i
  integer(IP) :: ieqn
  integer(IP) :: j
  
!
!  Multiply X by the diagonal of the matrix.
!
  b(1:n) = a(mu+1,1:n) * x(1:n)
!
!  Multiply X by the superdiagonals of the matrix.
!
  do i = mu, 1, -1
    do j = mu + 2 - i, n
      ieqn = i + j - mu - 1
      b(ieqn) = b(ieqn) + a(i,j) * x(j)
      b(j) = b(j) + a(i,j) * x(ieqn)
    end do
  end do
  return
end subroutine mv_pbu

subroutine dpbufa ( n, mu, a, info )
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


end module symmetric_band_matrix_mod
