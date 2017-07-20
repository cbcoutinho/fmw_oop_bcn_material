module iterative_solver
  use types
  implicit none
  private
  public :: mv_gb, mv_ge, mv_st
  public :: cg_gb, cg_ge, cg_st
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
    !    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the band matrix.
    !
    !    Input, real ( kind = 8 ) B(N), the right hand side vector.
    !
    !    Input/output, real ( kind = 8 ) X(N).
    !    On input, an estimate for the solution, which may be 0.
    !    On output, the approximate solution vector.  
    !
    implicit none

    integer(ip) :: ml
    integer(ip) :: mu
    integer(ip) :: n

    real ( kind = 8 ) a(2*ml+mu+1,n)
    real ( kind = 8 ) alpha
    real ( kind = 8 ) ap(n)
    real ( kind = 8 ) b(n)
    real ( kind = 8 ) beta
    integer(ip) :: it
    real ( kind = 8 ) p(n)
    real ( kind = 8 ) pap
    real ( kind = 8 ) pr
    real ( kind = 8 ) r(n)
    real ( kind = 8 ) rap
    real ( kind = 8 ) x(n)
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
    !    Input, real ( kind = 8 ) A(N,N), the matrix.
    !
    !    Input, real ( kind = 8 ) B(N), the right hand side vector.
    !
    !    Input/output, real ( kind = 8 ) X(N).
    !    On input, an estimate for the solution, which may be 0.
    !    On output,  the approximate solution vector.  
    !
    implicit none

    integer(ip) :: n

    real ( kind = 8 ) a(n,n)
    real ( kind = 8 ) alpha
    real ( kind = 8 ) ap(n)
    real ( kind = 8 ) b(n)
    real ( kind = 8 ) beta
    integer(ip) :: it
    real ( kind = 8 ) p(n)
    real ( kind = 8 ) pap
    real ( kind = 8 ) pr
    real ( kind = 8 ) r(n)
    real ( kind = 8 ) rap
    real ( kind = 8 ) x(n)
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
    !    Input, real ( kind = 8 ) A(NZ_NUM), the nonzero entries.
    !
    !    Input, real ( kind = 8 ) B(N), the right hand side vector.
    !
    !    Input/output, real ( kind = 8 ) X(N).
    !    On input, an estimate for the solution, which may be 0.
    !    On output, the approximate solution vector.  
    !
    implicit none

    integer(ip) :: n
    integer(ip) :: nz_num

    real ( kind = 8 ) a(nz_num)
    real ( kind = 8 ) alpha
    real ( kind = 8 ) ap(n)
    real ( kind = 8 ) b(n)
    real ( kind = 8 ) beta
    integer(ip) :: col(nz_num)
    integer(ip) :: it
    real ( kind = 8 ) p(n)
    real ( kind = 8 ) pap
    real ( kind = 8 ) pr
    real ( kind = 8 ) r(n)
    real ( kind = 8 ) rap
    integer(ip) :: row(nz_num)
    real ( kind = 8 ) x(n)
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
    !    Input, real ( kind = 8 ) A(2*ML+MU+1,N), the matrix.
    !
    !    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
    !
    !    Output, real ( kind = 8 ) B(M), the product A * x.
    !
    implicit none

    integer(ip) :: m
    integer(ip) :: ml
    integer(ip) :: mu
    integer(ip) :: n

    real ( kind = 8 ) a(2*ml+mu+1,n)
    real ( kind = 8 ) b(m)
    integer(ip) :: i
    integer(ip) :: j
    integer(ip) :: jhi
    integer(ip) :: jlo
    real ( kind = 8 ) x(n)

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
    !    Input, real ( kind = 8 ) A(M,N), the R8GE matrix.
    !
    !    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
    !
    !    Output, real ( kind = 8 ) B(M), the product A * x.
    !
    implicit none

    integer(ip) :: m
    integer(ip) :: n

    real ( kind = 8 ) a(m,n)
    real ( kind = 8 ) b(m)
    real ( kind = 8 ) x(n)

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

end module iterative_solver
