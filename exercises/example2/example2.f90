program example2
  use types
  use random_numbers
  use wathen_problem_names
  use matrix_names
  use solver_names
  use dense_matrix_names
  use direct_solver_names
  implicit none
  class(matrix_t), allocatable :: A
  class(solver_t), allocatable :: S
  type(wathen_problem_t) :: problem
  
  integer (ip) seed
  real (rp), allocatable :: x1(:)
  real (rp), allocatable :: x2(:)
  real (rp), allocatable :: b(:)
  real (rp) e
  integer (ip) n
  integer (ip) nx
  integer (ip) ny

  ! Read from command line:
  nx = 10
  ny = 10

  ! Write information
  write ( *, '(a,i6)' ) '  Elements in X direction NX = ', nx
  write ( *, '(a,i6)' ) '  Elements in Y direction NY = ', ny
  write ( *, '(a,i6)' ) '  Number of elements = ', nx * ny

  ! Allocate matrix and solver
  allocate(dense_matrix_t  :: A)
  allocate(direct_solver_t :: S)

  ! Compute the number of unknown and create the matrix
  call problem%setup(nx,ny,n,A)
  write ( *, '(a,i6)' ) '  Number of nodes N = ', n

  ! Set up a random solution X.
  seed = 123456789
  allocate ( x1(1:n) )
  call r8vec_uniform_01 ( n, seed, x1 )

  ! Define problem
  seed = 123456789
  call problem%assembly(seed,A)

  ! Setup solver
  call S%create(A)

  ! Set up a RHS and solve
  allocate ( b(1:n) )
  call A%apply(x1,b)
  
  allocate ( x2(1:n) )
  call S%apply(A,b,x2)

  !  Compute the maximum solution error.
  e = maxval ( abs ( x1 - x2 ) )
  write ( *, '(a,g14.6)' ) '  Maximum solution error is ', e 

  !  Free memory.
  deallocate ( b )
  deallocate ( x1 )
  deallocate ( x2 )

end program example2
