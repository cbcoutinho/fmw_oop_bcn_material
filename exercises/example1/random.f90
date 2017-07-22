module random_numbers_mod
  use types
  implicit none
  private
  
  integer(IP), parameter :: I4_HUGE = 2147483647

  public :: r8_uniform_01, r8vec_uniform_01
contains

  function r8_uniform_01 ( seed )

    !*****************************************************************************80
    !
    !! R8_UNIFORM_01 returns a unit pseudorandom R8.
    !
    !  Discussion:
    !
    !    An R8 is a real(RP) value.
    !
    !    For now, the input quantity SEED is an integer variable.
    !
    !    This routine implements the recursion
    !
    !      seed = 16807 * seed mod ( 2^31 - 1 )
    !      r8_uniform_01 = seed / ( 2^31 - 1 )
    !
    !    The integer arithmetic never requires more than 32 bits,
    !    including a sign bit.
    !
    !    If the initial seed is 12345, then the first three computations are
    !
    !      Input     Output      R8_UNIFORM_01
    !      SEED      SEED
    !
    !         12345   207482415  0.096616
    !     207482415  1790989824  0.833995
    !    1790989824  2035175616  0.947702
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 July 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Paul Bratley, Bennett Fox, Linus Schrage,
    !    A Guide to Simulation,
    !    Springer Verlag, pages 201-202, 1983.
    !
    !    Pierre L'Ecuyer,
    !    Random Number Generation,
    !    in Handbook of Simulation,
    !    edited by Jerry Banks,
    !    Wiley Interscience, page 95, 1998.
    !
    !    Bennett Fox,
    !    Algorithm 647:
    !    Implementation and Relative Efficiency of Quasirandom
    !    Sequence Generators,
    !    ACM Transactions on Mathematical Software,
    !    Volume 12, Number 4, pages 362-376, 1986.
    !
    !    Peter Lewis, Allen Goodman, James Miller
    !    A Pseudo-Random Number Generator for the System/360,
    !    IBM Systems Journal,
    !    Volume 8, pages 136-143, 1969.
    !
    !  Parameters:
    !
    !    Input/output, integer(ip) :: SEED, the "seed" value, which should
    !    NOT be 0. On output, SEED has been updated.
    !
    !    Output, real(RP) R8_UNIFORM_01, a new pseudorandom variate,
    !    strictly between 0 and 1.
    !
    implicit none
    integer(IP), intent(inout) :: seed
    real(RP)                   :: r8_uniform_01
    integer(IP) :: k
    
    if ( seed == 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
       write ( *, '(a)' ) '  Input value of SEED = 0.'
       stop 1
    end if

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
       seed = seed + I4_HUGE
    end if

    r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

    return
  end function r8_uniform_01
  subroutine r8mat_uniform_01 ( m, n, seed, r )

    !*****************************************************************************80
    !
    !! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
    !
    !  Discussion:
    !
    !    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    11 August 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Paul Bratley, Bennett Fox, Linus Schrage,
    !    A Guide to Simulation,
    !    Springer Verlag, pages 201-202, 1983.
    !
    !    Bennett Fox,
    !    Algorithm 647:
    !    Implementation and Relative Efficiency of Quasirandom
    !    Sequence Generators,
    !    ACM Transactions on Mathematical Software,
    !    Volume 12, Number 4, pages 362-376, 1986.
    !
    !    Peter Lewis, Allen Goodman, James Miller,
    !    A Pseudo-Random Number Generator for the System/360,
    !    IBM Systems Journal,
    !    Volume 8, pages 136-143, 1969.
    !
    !  Parameters:
    !
    !    Input, integer(ip) :: M, N, the number of rows and columns in
    !    the array.
    !
    !    Input/output, integer(ip) :: SEED, the "seed" value, which
    !    should NOT be 0.  On output, SEED has been updated.
    !
    !    Output, real(RP) R(M,N), the array of pseudorandom values.
    !
    implicit none

    integer(IP), intent(in)    :: m
    integer(IP), intent(in)    :: n
    integer(IP), intent(inout) :: seed
    real(RP)   , intent(out)   :: r(m,n)

    integer(IP) :: i
    integer(IP) :: j
    integer(IP) :: k
    

    do j = 1, n

       do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed < 0 ) then
             seed = seed + I4_HUGE
          end if

          r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

       end do
    end do

    return
  end subroutine r8mat_uniform_01
  subroutine r8vec_print ( n, a, title )

    !*****************************************************************************80
    !
    !! R8VEC_PRINT prints an R8VEC.
    !
    !  Discussion:
    !
    !    An R8VEC is a vector of R8's.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    22 August 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer(ip) :: N, the number of components of the vector.
    !
    !    Input, real(RP) A(N), the vector to be printed.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
    implicit none

    integer(IP)      , intent(in) :: n
    real(RP)         , intent(in) :: a(n)
    character (len=*), intent(in) :: title
    
    integer(IP) :: i
    
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '

    do i = 1, n
       write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
    end do

    return
  end subroutine r8vec_print
  subroutine r8vec_uniform_01 ( n, seed, r )

    !*****************************************************************************80
    !
    !! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
    !
    !  Discussion:
    !
    !    An R8VEC is a vector of R8's.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 July 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Paul Bratley, Bennett Fox, Linus Schrage,
    !    A Guide to Simulation,
    !    Springer Verlag, pages 201-202, 1983.
    !
    !    Bennett Fox,
    !    Algorithm 647:
    !    Implementation and Relative Efficiency of Quasirandom
    !    Sequence Generators,
    !    ACM Transactions on Mathematical Software,
    !    Volume 12, Number 4, pages 362-376, 1986.
    !
    !    Peter Lewis, Allen Goodman, James Miller
    !    A Pseudo-Random Number Generator for the System/360,
    !    IBM Systems Journal,
    !    Volume 8, pages 136-143, 1969.
    !
    !  Parameters:
    !
    !    Input, integer(ip) :: N, the number of entries in the vector.
    !
    !    Input/output, integer(ip) :: SEED, the "seed" value, which
    !    should NOT be 0.  On output, SEED has been updated.
    !
    !    Output, real(RP) R(N), the vector of pseudorandom values.
    !
    implicit none
    integer(IP), intent(in)    :: n
    integer(IP), intent(inout) :: seed
    real(RP)   , intent(out)   :: r(n)
    integer(IP) :: i
    integer(IP) :: k
    if ( seed == 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
       write ( *, '(a)' ) '  Input value of SEED = 0.'
       stop 1
    end if

    do i = 1, n

       k = seed / 127773

       seed = 16807 * ( seed - k * 127773 ) - k * 2836

       if ( seed < 0 ) then
          seed = seed + 2147483647
       end if

       r(i) = real ( seed, kind = 8 ) * 4.656612875D-10

    end do

    return
  end subroutine r8vec_uniform_01

end module random_numbers_mod
