module iterative_solver_mod
  use types_mod
  implicit none
  private
  public :: mv_gb, mv_ge, mv_st, mv_pbu
  public :: cg_gb, cg_ge, cg_st, cg_pbu
contains

  subroutine cg_gb ( n, ml, mu, a, b, x )
    !*****************************************************************************80
    !
    !! CG_GB uses the conjugate gradient method for a general banded (GB) matrix.
    !
    !  Discussion:
    !
    !    The linear system has the form A*x=b, where A is a positive-definite
    !    symmetric matrix.
    !
    !    The method is designed to reach the solution to the linear system
    !      A * x = b
    !    after N computational steps.  However, roundoff may introduce
    !    unacceptably large errors for some problems.  In such a case,
    !    calling the routine a second time, using the current solution estimate
    !    as the new starting guess, should result in improved results.
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
    !    Input, integer(ip) :: N, the order of the matrix.
    !
    !    Input, integer(ip) :: ML, MU, the lower and upper bandwidths.
    !
    !    Input, real(RP) A(2*ML+MU+1,N), the band matrix.
    !
    !    Input, real(RP) B(N), the right hand side vector.
    !
    !    Input/output, real(RP) X(N).
    !    On input, an estimate for the solution, which may be 0.
    !    On output, the approximate solution vector.  
    !
    implicit none
    integer(IP), intent(in)    :: n
    integer(IP), intent(in)    :: ml
    integer(IP), intent(in)    :: mu
    real(RP)   , intent(in)    :: a(2*ml+mu+1,n)
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
    call mv_gb ( n, n, ml, mu, a, x, ap )

    r(1:n) = b(1:n) - ap(1:n)
    p(1:n) = b(1:n) - ap(1:n)
    !
    !  Do the N steps of the conjugate gradient method.
    !
    do it = 1, 2*n
       !
       !  Compute the matrix*vector product AP = A*P.
       !
       call mv_gb ( n, n, ml, mu, a, p, ap )
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
          exit
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
  end subroutine cg_gb
  subroutine cg_ge ( n, a, b, x )

    !*****************************************************************************80
    !
    !! CG_GE uses the conjugate gradient method for a general storage (GE) matrix.
    !
    !  Discussion:
    !
    !    The linear system has the form A*x=b, where A is a positive-definite
    !    symmetric matrix, stored as a full storage matrix.
    !
    !    The method is designed to reach the solution to the linear system
    !      A * x = b
    !    after N computational steps.  However, roundoff may introduce
    !    unacceptably large errors for some problems.  In such a case,
    !    calling the routine a second time, using the current solution estimate
    !    as the new starting guess, should result in improved results.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    01 June 2014
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
    !    Input, integer(ip) :: N, the order of the matrix.
    !
    !    Input, real(RP) A(N,N), the matrix.
    !
    !    Input, real(RP) B(N), the right hand side vector.
    !
    !    Input/output, real(RP) X(N).
    !    On input, an estimate for the solution, which may be 0.
    !    On output,  the approximate solution vector.  
    !
    implicit none

    integer(IP), intent(in)    :: n
    real(RP)   , intent(in)    :: a(n,n)
    real(RP)   , intent(in)    :: b(n)
    real(RP)   , intent(inout) :: x(n)
    
    real(RP) :: alpha
    real(RP) :: ap(n)
    
    real(RP)    :: beta
    integer(IP) :: it
    real(RP)    :: p(n)
    real(RP)    :: pap
    real(RP)    :: pr
    real(RP)    :: r(n)
    real(RP)    :: rap
    
    !
    !  Initialize
    !    AP = A * x,
    !    R  = b - A * x,
    !    P  = b - A * x.
    !
    ap = matmul ( a, x )

    r(1:n) = b(1:n) - ap(1:n)
    p(1:n) = b(1:n) - ap(1:n)
    !
    !  Do the N steps of the conjugate gradient method.
    !
    do it = 1, 2*n
       !
       !  Compute the matrix*vector product AP = A*P.
       !
       ap = matmul ( a, p )
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
       rap = dot_product ( r(1:n), ap(1:n) )

       beta = - rap / pap
       !
       !  Update the perturbation vector
       !    P = R + BETA * P.
       !
       p(1:n) = r(1:n) + beta * p(1:n)

    end do

    return
  end subroutine cg_ge
  subroutine cg_st ( n, nz_num, row, col, a, b, x )

    !*****************************************************************************80
    !
    !! CG_ST uses the conjugate gradient method for a sparse triplet (ST) matrix.
    !
    !  Discussion:
    !
    !    The linear system has the form A*x=b, where A is a positive-definite
    !    symmetric matrix, stored as a full storage matrix.
    !
    !    The method is designed to reach the solution to the linear system
    !      A * x = b
    !    after N computational steps.  However, roundoff may introduce
    !    unacceptably large errors for some problems.  In such a case,
    !    calling the routine a second time, using the current solution estimate
    !    as the new starting guess, should result in improved results.
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
    !    Input, integer(ip) :: N, the order of the matrix.
    !
    !    Input, integer(ip) :: NZ_NUM, the number of nonzeros.
    !
    !    Input, integer(ip) :: ROW(NZ_NUM), COL(NZ_NUM), the row and column 
    !    indices of the nonzero entries.
    !
    !    Input, real(RP) A(NZ_NUM), the nonzero entries.
    !
    !    Input, real(RP) B(N), the right hand side vector.
    !
    !    Input/output, real(RP) X(N).
    !    On input, an estimate for the solution, which may be 0.
    !    On output, the approximate solution vector.  
    !
    implicit none

    integer(IP), intent(in)    :: n
    integer(IP), intent(in)    :: nz_num
    integer(IP), intent(in)    :: row(nz_num)
    integer(IP), intent(in)    :: col(nz_num)
    real(RP)   , intent(in)    :: a(nz_num)
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
    call mv_st ( n, n, nz_num, row, col, a, x, ap )

    r(1:n) = b(1:n) - ap(1:n)
    p(1:n) = b(1:n) - ap(1:n)
    !
    !  Do the N steps of the conjugate gradient method.
    !
    do it = 1, 2*n
       !
       !  Compute the matrix*vector product AP = A*P.
       !
       call mv_st ( n, n, nz_num, row, col, a, p, ap )
       !
       !  Compute the dot products
       !    PAP = P*AP,
       !    PR  = P*R
       !  Set
       !    ALPHA = PR / PAP.
       !
       pap = dot_product ( p, ap )
       pr =  dot_product ( p, r )

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
  end subroutine cg_st
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
    !    Input, real(RP) A(NZ_NUM), the nonzero values in the matrix.
    !
    !    Input, real(RP) X(N), the vector to be multiplied.
    !
    !    Output, real(RP) B(M), the product A*X.
    !
    implicit none

    integer(IP), intent(in)  :: m
    integer(IP), intent(in)  :: n
    integer(IP), intent(in)  :: nz_num
    integer(IP), intent(in)  :: row(nz_num)
    integer(IP), intent(in)  :: col(nz_num)
    real(RP)   , intent(in)  :: a(nz_num)
    real(RP)   , intent(out) :: b(m)
    
    integer(IP) :: k
    real(RP) :: x(n)

    b(1:m) = 0.0D+00
    do k = 1, nz_num
       b(row(k)) = b(row(k)) + a(k) * x(col(k))
    end do

    return
  end subroutine mv_st
  
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


end module iterative_solver_mod
