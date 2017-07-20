module sparse_matrix_names
  use types
  use blas_names
  use matrix_names
  implicit none
  private

  type, extends(matrix_t) :: sparse_matrix_t
     private
     integer(ip) :: nz, k
     integer(ip), allocatable :: row(:), col(:)
     real(rp)   , allocatable :: a(:)
   contains
     procedure :: create    => sparse_matrix_create
     procedure :: assembly  => sparse_matrix_assembly
     procedure :: apply     => sparse_matrix_apply
     procedure :: factorize => sparse_matrix_factorize
     procedure :: backsolve => sparse_matrix_backsolve
  end type sparse_matrix_t
  public :: sparse_matrix_t

contains
 
  subroutine sparse_matrix_create(this,n,ml,mu,nz)
    class(sparse_matrix_t), intent(inout) :: this
    integer(ip)          , intent(in)    :: n
    integer(ip), optional, intent(in)    :: ml,mu,nz
    call this%set_size(n)
    this%nz = nz
    allocate ( this%row(1:this%nz) )
    allocate ( this%col(1:this%nz) )
    allocate ( this%a(1:this%nz) )
    this%a = 0.0_rp
    this%k = 0
  end subroutine sparse_matrix_create

  subroutine sparse_matrix_assembly(this,i,j,a) 
    implicit none
    class(sparse_matrix_t), intent(inout) :: this
    integer(ip)           , intent(in)    :: i,j
    real(rp)              , intent(in)    :: a
    this%k = this%k + 1
    this%row(this%k) = i
    this%col(this%k) = j
    this%a(this%k)   = a
  end subroutine sparse_matrix_assembly

  subroutine sparse_matrix_apply(this,x,y) 
    implicit none
    class(sparse_matrix_t), intent(in) :: this
    real(rp)       , intent(in)       :: x(:)
    real(rp)       , intent(inout)    :: y(:)
    call mv_st ( this%get_size(), this%get_size(), this%nz, this%row, this%col, this%a, x, y)
  end subroutine sparse_matrix_apply

  subroutine sparse_matrix_factorize(this, factors, pivots) 
    implicit none
    class(sparse_matrix_t)       , intent(in)    :: this
    class(matrix_t), allocatable, intent(inout) :: factors
    integer(ip)                 , intent(inout) :: pivots(:)
    write(*,*) 'Sparse matrix factorization not implemented'
  end subroutine sparse_matrix_factorize

  subroutine sparse_matrix_backsolve(this,pivots,rhs,x) 
    class(sparse_matrix_t), intent(in)    :: this
    integer(ip)    , intent(in)    :: pivots(:)
    real(rp)       , intent(in)    :: rhs(:)
    real(rp)       , intent(inout) :: x(:)
    write(*,*) 'Sparse matrix backsolve not implemented'
  end subroutine sparse_matrix_backsolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Legacy code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mv_st ( m, n, nz_num, row, col, a, x, b )

    !*****************************************************************************80
    !
    !! MV_ST multiplies a sparse triple matrix times a vector.
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
    !    Input, integer(ip) :: M, N, the number of rows and columns.
    !
    !    Input, integer(ip) :: NZ_NUM, the number of nonzero values.
    !
    !    Input, integer(ip) :: ROW(NZ_NUM), COL(NZ_NUM), the row and 
    !    column indices.
    !
    !    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero values in the matrix.
    !
    !    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
    !
    !    Output, real ( kind = 8 ) B(M), the product A*X.
    !
    implicit none

    integer(ip) :: m
    integer(ip) :: n
    integer(ip) :: nz_num

    real ( kind = 8 ) a(nz_num)
    real ( kind = 8 ) b(m)
    integer(ip) :: col(nz_num)
    integer(ip) :: k
    integer(ip) :: row(nz_num)
    real ( kind = 8 ) x(n)

    b(1:m) = 0.0D+00
    do k = 1, nz_num
       b(row(k)) = b(row(k)) + a(k) * x(col(k))
    end do

    return
  end subroutine mv_st


end module sparse_matrix_names
