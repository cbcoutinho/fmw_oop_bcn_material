module wathen_problem_mod
  use types_mod,          only: IP, RP
  use random_numbers_mod, only: r8_uniform_01, r8vec_uniform_01
  use matrix_mod,         only: matrix_t
  implicit none
  private

  real(RP), parameter :: EM(8,8) =  reshape ( [ &
  6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0, &
  -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, &
  2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0, &
  -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, &
  3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0, &
  -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, &
  2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0, &
  -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0 ], &
  [8, 8] )


  public :: wathen_problem_t
  type wathen_problem_t
    private
    integer(IP) :: nx,ny
  contains
    procedure :: setup    => wathen_setup
    procedure :: assembly => wathen_fill
  end type wathen_problem_t

  public :: mv_pbu, dpbufa, dpbusl

contains

#include "../symmetric_band_matrix_subroutines/symmetric_band_matrix_subroutines.f90"

  subroutine wathen_setup( this, nx, ny, n, a )
    class(wathen_problem_t), intent(inout) :: this
    integer(IP)    , intent(in)    :: nx
    integer(IP)    , intent(in)    :: ny
    integer(IP)    , intent(out)   :: n
    class(matrix_t), intent(inout) :: a
    integer(IP) :: ml,md,mu,nz
    this%nx = nx
    this%ny = ny
    call wathen_order ( nx, ny, n )
    call wathen_bandwidth ( nx, ny, ml, md, mu )
    call wathen_st_size ( nx, ny, nz )
    call a%create(n,ml,mu,nz)
  end subroutine wathen_setup

  subroutine wathen_fill ( this, seed, a )

    !*****************************************************************************80
    !
    ! WATHEN_FILL returns the Wathen matrix as an abstract matrix matrix_t,
    ! see the discutssion in the legacy WATHEN_GE
    !
    implicit none
    class(wathen_problem_t), intent(inout) :: this
    integer(IP)            , intent(inout) :: seed
    class(matrix_t)        , intent(inout) :: a

    ! Esta mierda hay que hacerla bien!!!
    real ( kind = 8 ), dimension ( 8, 8 ), save :: em =  reshape ( [ &
    6.0, -6.0,  2.0, -8.0,  3.0, -8.0,  2.0, -6.0, &
    -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, &
    2.0, -6.0,  6.0, -6.0,  2.0, -8.0,  3.0, -8.0, &
    -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, -8.0, 16.0, &
    3.0, -8.0,  2.0, -6.0,  6.0, -6.0,  2.0, -8.0, &
    -8.0, 16.0, -8.0, 20.0, -6.0, 32.0, -6.0, 20.0, &
    2.0, -8.0,  3.0, -8.0,  2.0, -6.0,  6.0, -6.0, &
    -6.0, 20.0, -8.0, 16.0, -8.0, 20.0, -6.0, 32.0 ], &
    [ 8, 8 ] )
    integer(IP) :: i
    integer(IP) :: j
    integer(IP) :: kcol
    integer(IP) :: krow
    integer(IP) :: node(8)
    !real ( kind = 8 ) r8_uniform_01
    real ( kind = 8 ) rho

    do j = 1, this%ny
      do i = 1, this%nx

        node(1) = 3 * j * this%nx + 2 * j + 2 * i + 1
        node(2) = node(1) - 1
        node(3) = node(1) - 2
        node(4) = ( 3 * j - 1 ) * this%nx + 2 * j + i - 1
        node(5) = ( 3 * j - 3 ) * this%nx + 2 * j + 2 * i - 3
        node(6) = node(5) + 1
        node(7) = node(5) + 2
        node(8) = node(4) + 1

        rho = 100.0D+00 * r8_uniform_01 ( seed )

        do krow = 1, 8
          do kcol = 1, 8
            call a%assembly(node(krow),node(kcol),rho * em(krow,kcol))
          end do
        end do

      end do
    end do

    return
  end subroutine wathen_fill

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Legacy code
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wathen_bandwidth ( nx, ny, l, d, u )

    !*****************************************************************************80
    !
    !! WATHEN_BANDWIDTH returns the bandwidth of the WATHEN matrix.
    !
    !  Discussion:
    !
    !    The bandwidth measures the minimal number of contiguous diagonals,
    !    including the central diagonal, which contain all the nonzero elements
    !    of a matrix.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    03 June 2014
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Nicholas Higham,
    !    Algorithm 694: A Collection of Test Matrices in MATLAB,
    !    ACM Transactions on Mathematical Software,
    !    Volume 17, Number 3, September 1991, pages 289-305.
    !
    !    Andrew Wathen,
    !    Realistic eigenvalue bounds for the Galerkin mass matrix,
    !    IMA Journal of Numerical Analysis,
    !    Volume 7, 1987, pages 449-457.
    !
    !  Parameters:
    !
    !    Input, integer(ip) :: NX, NY, values which determine the size of A.
    !
    !    Output, integer(ip) :: L, D, U, the lower, diagonal, and upper
    !    bandwidths of the matrix,
    !
    implicit none
    integer(IP), intent(in)  :: nx
    integer(IP), intent(in)  :: ny
    integer(IP), intent(out) :: l
    integer(IP), intent(out) :: d
    integer(IP), intent(out) :: u
    l = 3 * nx + 4
    d = 1
    u = 3 * nx + 4
    return
  end subroutine wathen_bandwidth

  subroutine wathen_gb ( nx, ny, n, seed, a )

    !*****************************************************************************80
    !
    !! WATHEN_GB returns the Wathen matrix, using general banded (GB) storage.
    !
    !  Discussion:
    !
    !    The Wathen matrix is a finite element matrix which is sparse.
    !
    !    The entries of the matrix depend in part on a physical quantity
    !    related to density.  That density is here assigned random values between
    !    0 and 100.
    !
    !    The matrix order N is determined by the input quantities NX and NY,
    !    which would usually be the number of elements in the X and Y directions.
    !    The value of N is
    !
    !      N = 3*NX*NY + 2*NX + 2*NY + 1,
    !
    !    The matrix is the consistent mass matrix for a regular NX by NY grid
    !    of 8 node serendipity elements.
    !
    !    The local element numbering is
    !
    !      3--2--1
    !      |     |
    !      4     8
    !      |     |
    !      5--6--7
    !
    !    Here is an illustration for NX = 3, NY = 2:
    !
    !     23-24-25-26-27-28-29
    !      |     |     |     |
    !     19    20    21    22
    !      |     |     |     |
    !     12-13-14-15-16-17-18
    !      |     |     |     |
    !      8     9    10    11
    !      |     |     |     |
    !      1--2--3--4--5--6--7
    !
    !    For this example, the total number of nodes is, as expected,
    !
    !      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
    !
    !    The matrix is symmetric positive definite for any positive values of the
    !    density RHO(X,Y).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 July 2014
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Nicholas Higham,
    !    Algorithm 694: A Collection of Test Matrices in MATLAB,
    !    ACM Transactions on Mathematical Software,
    !    Volume 17, Number 3, September 1991, pages 289-305.
    !
    !    Andrew Wathen,
    !    Realistic eigenvalue bounds for the Galerkin mass matrix,
    !    IMA Journal of Numerical Analysis,
    !    Volume 7, Number 4, October 1987, pages 449-457.
    !
    !  Parameters:
    !
    !    Input, integer(ip) :: NX, NY, values which determine the size
    !    of the matrix.
    !
    !    Input, integer(ip) :: N, the number of rows and columns.
    !
    !    Input/output, integer(ip) :: SEED, the random number seed.
    !
    !    Output, real(RP) A(9*NX+13,N), the matrix.
    !
    implicit none
    integer(IP), intent(in)    :: nx
    integer(IP), intent(in)    :: ny
    integer(IP), intent(in)    :: n
    integer(IP), intent(inout) :: seed
    real(RP)   , intent(out)   :: a(9*nx+13,n)
    integer(IP) :: i
    integer(IP) :: ii
    integer(IP) :: j
    integer(IP) :: jj
    integer(IP) :: kcol
    integer(IP) :: krow
    integer(IP) :: ml
    integer(IP) :: mu
    integer(IP) :: node(8)
    real(RP) :: rho


    ml = 3 * nx + 4
    mu = 3 * nx + 4

    a(1:9*nx+13,1:n) = 0.0D+00

    do j = 1, ny
      do i = 1, nx

        node(1) = 3 * j * nx + 2 * j + 2 * i + 1
        node(2) = node(1) - 1
        node(3) = node(1) - 2
        node(4) = ( 3 * j - 1 ) * nx + 2 * j + i - 1
        node(5) = ( 3 * j - 3 ) * nx + 2 * j + 2 * i - 3
        node(6) = node(5) + 1
        node(7) = node(5) + 2
        node(8) = node(4) + 1

        rho = 100.0D+00 * r8_uniform_01 ( seed )

        do krow = 1, 8
          do kcol = 1, 8
            ii = node(krow);
            jj = node(kcol);
            a(ii-jj+ml+mu+1,jj) = a(ii-jj+ml+mu+1,jj) &
            + rho * EM(krow,kcol)
          end do
        end do

      end do
    end do

    return
  end subroutine wathen_gb

  subroutine wathen_ge ( nx, ny, n, seed, a )

    !*****************************************************************************80
    !
    !! WATHEN_GE returns the Wathen matrix as a general storage (GE) matrix.
    !
    !  Discussion:
    !
    !    The Wathen matrix is a finite element matrix which is sparse.
    !
    !    The entries of the matrix depend in part on a physical quantity
    !    related to density.  That density is here assigned random values between
    !    0 and 100.
    !
    !    The matrix order N is determined by the input quantities NX and NY,
    !    which would usually be the number of elements in the X and Y directions.
    !    The value of N is
    !
    !      N = 3*NX*NY + 2*NX + 2*NY + 1,
    !
    !    The matrix is the consistent mass matrix for a regular NX by NY grid
    !    of 8 node serendipity elements.
    !
    !    The local element numbering is
    !
    !      3--2--1
    !      |     |
    !      4     8
    !      |     |
    !      5--6--7
    !
    !    Here is an illustration for NX = 3, NY = 2:
    !
    !     23-24-25-26-27-28-29
    !      |     |     |     |
    !     19    20    21    22
    !      |     |     |     |
    !     12-13-14-15-16-17-18
    !      |     |     |     |
    !      8     9    10    11
    !      |     |     |     |
    !      1--2--3--4--5--6--7
    !
    !    For this example, the total number of nodes is, as expected,
    !
    !      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
    !
    !    The matrix is symmetric positive definite for any positive values of the
    !    density RHO(X,Y).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 July 2014
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Nicholas Higham,
    !    Algorithm 694: A Collection of Test Matrices in MATLAB,
    !    ACM Transactions on Mathematical Software,
    !    Volume 17, Number 3, September 1991, pages 289-305.
    !
    !    Andrew Wathen,
    !    Realistic eigenvalue bounds for the Galerkin mass matrix,
    !    IMA Journal of Numerical Analysis,
    !    Volume 7, Number 4, October 1987, pages 449-457.
    !
    !  Parameters:
    !
    !    Input, integer(ip) :: NX, NY, values which determine the size
    !    of the matrix.
    !
    !    Input, integer(ip) :: N, the number of rows and columns.
    !
    !    Input/output, integer(ip) :: SEED, the random number seed.
    !
    !    Output, real(RP) A(N,N), the matrix.
    !
    implicit none
    integer(IP), intent(in)    :: nx
    integer(IP), intent(in)    :: ny
    integer(IP), intent(in)    :: n
    integer(IP), intent(inout) :: seed
    real(RP)   , intent(out)   :: a(n,n)

    integer(IP) :: i
    integer(IP) :: j
    integer(IP) :: kcol
    integer(IP) :: krow
    integer(IP) :: node(8)
    real(RP) :: rho

    a(1:n,1:n) = 0.0D+00

    do j = 1, ny
      do i = 1, nx

        node(1) = 3 * j * nx + 2 * j + 2 * i + 1
        node(2) = node(1) - 1
        node(3) = node(1) - 2
        node(4) = ( 3 * j - 1 ) * nx + 2 * j + i - 1
        node(5) = ( 3 * j - 3 ) * nx + 2 * j + 2 * i - 3
        node(6) = node(5) + 1
        node(7) = node(5) + 2
        node(8) = node(4) + 1

        rho = 100.0D+00 * r8_uniform_01 ( seed )

        do krow = 1, 8
          do kcol = 1, 8
            a(node(krow),node(kcol)) = a(node(krow),node(kcol)) &
            + rho * EM(krow,kcol)
          end do
        end do

      end do
    end do

    return
  end subroutine wathen_ge

  subroutine wathen_order ( nx, ny, n )

    !*****************************************************************************80
    !
    !! WATHEN_ORDER returns the order of the WATHEN matrix.
    !
    !  Discussion:
    !
    !    N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 January 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Nicholas Higham,
    !    Algorithm 694: A Collection of Test Matrices in MATLAB,
    !    ACM Transactions on Mathematical Software,
    !    Volume 17, Number 3, September 1991, pages 289-305.
    !
    !    Andrew Wathen,
    !    Realistic eigenvalue bounds for the Galerkin mass matrix,
    !    IMA Journal of Numerical Analysis,
    !    Volume 7, 1987, pages 449-457.
    !
    !  Parameters:
    !
    !    Input, integer(ip) :: NX, NY, values which determine the size of A.
    !
    !    Output, integer(ip) :: N, the order of the matrix,
    !    as determined by NX and NY.
    !
    implicit none


    integer(IP), intent(in)  :: nx
    integer(IP), intent(in)  :: ny
    integer(IP), intent(out) :: n
    n = 3 * nx * ny + 2 * nx + 2 * ny + 1
    return
  end subroutine wathen_order

  subroutine wathen_st ( nx, ny, nz_num, seed, row, col, a )

    !*****************************************************************************80
    !
    !! WATHEN_ST: Wathen matrix stored in sparse triplet (ST) format.
    !
    !  Discussion:
    !
    !    When dealing with sparse matrices in MATLAB, it can be much more efficient
    !    to work first with a triple of I, J, and X vectors, and only once
    !    they are complete, convert to MATLAB's sparse format.
    !
    !    The Wathen matrix is a finite element matrix which is sparse.
    !
    !    The entries of the matrix depend in part on a physical quantity
    !    related to density.  That density is here assigned random values between
    !    0 and 100.
    !
    !    The matrix order N is determined by the input quantities NX and NY,
    !    which would usually be the number of elements in the X and Y directions.
    !
    !    The value of N is
    !
    !      N = 3*NX*NY + 2*NX + 2*NY + 1,
    !
    !    The matrix is the consistent mass matrix for a regular NX by NY grid
    !    of 8 node serendipity elements.
    !
    !    The local element numbering is
    !
    !      3--2--1
    !      |     |
    !      4     8
    !      |     |
    !      5--6--7
    !
    !    Here is an illustration for NX = 3, NY = 2:
    !
    !     23-24-25-26-27-28-29
    !      |     |     |     |
    !     19    20    21    22
    !      |     |     |     |
    !     12-13-14-15-16-17-18
    !      |     |     |     |
    !      8     9    10    11
    !      |     |     |     |
    !      1--2--3--4--5--6--7
    !
    !    For this example, the total number of nodes is, as expected,
    !
    !      N = 3 * 3 * 2 + 2 * 2 + 2 * 3 + 1 = 29
    !
    !    The matrix is symmetric positive definite for any positive values of the
    !    density RHO(X,Y).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 July 2014
    !
    !  Author:
    !
    !    John Burkardt.
    !
    !  Reference:
    !
    !    Nicholas Higham,
    !    Algorithm 694: A Collection of Test Matrices in MATLAB,
    !    ACM Transactions on Mathematical Software,
    !    Volume 17, Number 3, September 1991, pages 289-305.
    !
    !    Andrew Wathen,
    !    Realistic eigenvalue bounds for the Galerkin mass matrix,
    !    IMA Journal of Numerical Analysis,
    !    Volume 7, Number 4, October 1987, pages 449-457.
    !
    !  Parameters:
    !
    !    Input, integer(ip) :: NX, NY, values which determine the size of
    !    the matrix.
    !
    !    Input, integer(ip) :: NZ_NUM, the number of values used to
    !    describe the matrix.
    !
    !    Input/output, integer(ip) :: SEED, the random number seed.
    !
    !    Output, integer(ip) :: ROW(NZ_NUM), COL(NZ_NUM), the row and
    !    column indices of the nonzero entries.
    !
    !    Output, real(RP) A(NZ_NUM), the nonzero entries of the matrix.
    !
    implicit none

    integer(IP), intent(in)    :: nx
    integer(IP), intent(in)    :: ny
    integer(IP), intent(in)    :: nz_num
    integer(IP), intent(inout) :: seed
    integer(IP), intent(out)   :: row(nz_num)
    integer(IP), intent(out)   :: col(nz_num)
    real(RP)   , intent(out)   :: a(nz_num)

    integer(IP) :: i
    integer(IP) :: j
    integer(IP) :: k
    integer(IP) :: kcol
    integer(IP) :: krow
    integer(IP) :: node(8)
    real(RP) :: rho

    row(1:nz_num) = 0
    col(1:nz_num) = 0
    a(1:nz_num) = 0.0D+00

    k = 0
    do j = 1, ny
      do i = 1, nx
        node(1) = 3 * j * nx + 2 * j + 2 * i + 1
        node(2) = node(1) - 1
        node(3) = node(1) - 2
        node(4) = ( 3 * j - 1 ) * nx + 2 * j + i - 1
        node(5) = ( 3 * j - 3 ) * nx + 2 * j + 2 * i - 3
        node(6) = node(5) + 1
        node(7) = node(5) + 2
        node(8) = node(4) + 1

        rho = 100.0D+00 * r8_uniform_01 ( seed )
        do krow = 1, 8
          do kcol = 1, 8
            k = k + 1
            row(k) = node(krow)
            col(k) = node(kcol)
            a(k) = rho * EM(krow,kcol)
          end do
        end do
      end do
    end do

    return
  end subroutine wathen_st

  subroutine wathen_st_size ( nx, ny, nz_num )

    !*****************************************************************************80
    !
    !! WATHEN_ST_SIZE: Size of Wathen matrix stored in sparse triplet format.
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
    !    John Burkardt.
    !
    !  Reference:
    !
    !    Nicholas Higham,
    !    Algorithm 694: A Collection of Test Matrices in MATLAB,
    !    ACM Transactions on Mathematical Software,
    !    Volume 17, Number 3, September 1991, pages 289-305.
    !
    !    Andrew Wathen,
    !    Realistic eigenvalue bounds for the Galerkin mass matrix,
    !    IMA Journal of Numerical Analysis,
    !    Volume 7, Number 4, October 1987, pages 449-457.
    !
    !  Parameters:
    !
    !    Input, integer NX, NY, values which determine the size of the matrix.
    !
    !    Output, integer NZ_NUM, the number of items of data used to describe
    !    the matrix.
    !
    implicit none
    integer(IP), intent(in)  :: nx
    integer(IP), intent(in)  :: ny
    integer(IP), intent(out) :: nz_num
    nz_num = nx * ny * 64
    return
  end subroutine wathen_st_size
end module wathen_problem_mod
