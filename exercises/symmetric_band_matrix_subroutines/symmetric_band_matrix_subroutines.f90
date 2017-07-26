!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Legacy code
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine wathen_pbu ( nx, ny, n, seed, a )

  !*****************************************************************************80
  !
  !! wathen_pbu  returns the Wathen matrix, using symmetric banded upper (PBU) storage.
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
  !    Output, real(RP) A(3*NX+5,N), the matrix.
  !
  implicit none
  integer(IP), intent(in)    :: nx
  integer(IP), intent(in)    :: ny
  integer(IP), intent(in)    :: n
  integer(IP), intent(inout) :: seed
  real(RP)   , intent(out)   :: a(3*nx+5,n)
  integer(IP) :: i
  integer(IP) :: ii
  integer(IP) :: j
  integer(IP) :: jj
  integer(IP) :: kcol
  integer(IP) :: krow
  integer(IP) :: mu
  integer(IP) :: node(8)
  real(RP) :: rho


  mu = 3 * nx + 4

  a(1:3*nx+5,1:n) = 0.0D+00

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
          if (ii <= jj) then
            a(ii-jj+mu+1,jj) = a(ii-jj+mu+1,jj) &
            + rho * EM(krow,kcol)
          end if
        end do
      end do

    end do
  end do
end subroutine wathen_pbu

subroutine cg_pbu ( n, mu, a, b, x )
  !*****************************************************************************80
  !
  !! PBU_CG uses the conjugate gradient method on an PBU system.
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
  !    The method is designed to reach the solution after N computational
  !    steps.  However, roundoff may introduce unacceptably large errors for
  !    some problems.  In such a case, calling the routine again, using
  !    the computed solution as the new starting estimate, should improve
  !    the results.
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
  !  Reference:
  !
  !    Frank Beckman,
  !    The Solution of Linear Equations by the Conjugate Gradient Method,
  !    in Mathematical Methods for Digital Computers,
  !    edited by John Ralston, Herbert Wilf,
  !    Wiley, 1967,
  !    ISBN: 0471706892,
  !    LC: QA76.5.R3.
  !
  !  Parameters:
  !
  !    Input, integer(IP) N, the order of the matrix.
  !    N must be positive.
  !
  !    Input, integer(IP) MU, the number of superdiagonals.
  !    MU must be at least 0, and no more than N-1.
  !
  !    Input, real(RP) A(MU+1,N), the PBU matrix.
  !
  !    Input, real(RP) B(N), the right hand side vector.
  !
  !    Input/output, real(RP) X(N).
  !    On input, an estimate for the solution, which may be 0.
  !    On output, the approximate solution vector.
  !
  implicit none

  integer(IP), intent(in)    :: n
  integer(IP), intent(in)    :: mu
  real(RP)   , intent(in)    :: a(mu+1,n)
  real(RP)   , intent(in)    :: b(n)
  real(RP)   , intent(inout) :: x(n)

  real(RP) :: alpha
  real(RP) :: ap(n)

  real(RP) :: beta
  integer(IP) :: it
  real(RP) :: p(n)
  real(RP) :: pap
  real(RP) :: pr
  real(RP) :: r(n)
  real(RP) :: rap

  !
  !  Initialize
  !    AP = A * x,
  !    R  = b - A * x,
  !    P  = b - A * x.
  !
  call mv_pbu ( n, n, mu, a, x, ap )

  r(1:n) = b(1:n) - ap(1:n)
  p(1:n) = b(1:n) - ap(1:n)
  !
  !  Do the N steps of the conjugate gradient method.
  !
  do it = 1, n
    !
    !  Compute the matrix*vector product AP=A*P.
    !
    call mv_pbu ( n, n, mu, a, p, ap )
    !
    !  Compute the dot products
    !    PAP = P*AP,
    !    PR  = P*R
    !  Set
    !    ALPHA = PR / PAP.
    !
    pap = dot_product ( p, ap )
    pr = dot_product ( p, r )

    if ( pap == 0.0D+00 ) then
      return
    end if

    alpha = pr / pap
    !
    !  Set
    !    X = X + ALPHA * P
    !    R = R - ALPHA * AP.
    !
    x(1:n) = x(1:n) + alpha * p(1:n)
    r(1:n) = r(1:n) - alpha * ap(1:n)
    !
    !  Compute the vector dot product
    !    RAP = R*AP
    !  Set
    !    BETA = - RAP / PAP.
    !
    rap = dot_product ( r, ap )

    beta = - rap / pap
    !
    !  Update the perturbation vector
    !    P = R + BETA * P.
    !
    p(1:n) = r(1:n) + beta * p(1:n)

  end do
  return
end subroutine cg_pbu

subroutine mv_pbu ( m, n, mu, a, x, b )

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

  integer(IP), intent(in)    :: m
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
