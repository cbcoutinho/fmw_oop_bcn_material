#include "mcheck.i90"
program example2
  use types_mod
  use random_numbers_mod
  use wathen_problem_mod
  use matrix_mod
  use matrix_factory_mod
  use solver_mod
  use solver_factory_mod
  
  implicit none
  class(matrix_t), allocatable :: A
  class(solver_t), allocatable :: S
  type(wathen_problem_t) :: problem
  
  integer (ip) :: seed
  real (rp), allocatable :: x1(:)
  real (rp), allocatable :: x2(:)
  real (rp), allocatable :: b(:)
  real (rp) :: e
  integer (ip) :: n
  integer (ip) :: nx
  integer (ip) :: ny
  character(:), allocatable :: matrix_type, solver_type

  ! Read from command line:
  call read_and_check_command_line_parameters(  nx, ny, matrix_type, solver_type ) 

  ! Write information
  write ( *, '(a,i6)' ) '  Elements in X direction NX = ', nx
  write ( *, '(a,i6)' ) '  Elements in Y direction NY = ', ny
  write ( *, '(a,i6)' ) '  Number of elements = ', nx * ny
  write ( *, '(a,a)' ) '  Matrix type = ', matrix_type
  write ( *, '(a,a)' ) '  Solver type = ', solver_type
  
  ! Select the dynamic type of A and S by calling the appropiate
  ! factory method (single method factory OO design pattern)
  call matrix_factory(matrix_type, A)
  call solver_factory(solver_type, S)

  ! Compute the number of unknowns and create the matrix
  call problem%setup(nx,ny,n,A)
  write ( *, '(a,i6)' ) '  Number of nodes N = ', n

  ! Set up a random solution X.
  seed = 123456789
  allocate ( x1(1:n) )
  call r8vec_uniform_01 ( n, seed, x1 )

  ! Define problem
  seed = 123456789
  call problem%assembly(seed,A)

  ! Set up solver
  call S%create(A)

  ! Set up a RHS and solve
  allocate ( b(1:n) )
  call A%apply(x1,b)
  
  allocate ( x2(1:n) )
  call S%solve(b,x2)

  !  Compute the maximum solution error.
  e = maxval ( abs ( x1 - x2 ) )
  write ( *, '(a,g14.6)' ) '  Maximum solution error is ', e 

  !  Free memory
  call A%free()
  call S%free()
  deallocate(A)
  deallocate(S)
  deallocate ( b )
  deallocate ( x1 )
  deallocate ( x2 )
contains
  subroutine read_and_check_command_line_parameters(  NX, NY, matrix_type, solver_type ) 
     implicit none
     integer(IP)              , intent(out)   :: NX
     integer(IP)              , intent(out)   :: NY
     character(:), allocatable, intent(inout) :: matrix_type
     character(:), allocatable, intent(inout) :: solver_type
     
     character(:), parameter :: USAGE_ERROR_MSG       = "Usage: example2 NX NY matrix_type solver_type"
     character(:), parameter :: MATRIX_SOLVER_TYPE_ERROR_MSG = "matrix/solver type combination not allowed"
     character(:), parameter :: IO_ERROR_MSG = "Error while reading NX or NY, these should be integers!!!"
     CHARACTER(len=256) :: arg
     integer(IP) :: istat
     
     if (allocated(matrix_type)) deallocate(matrix_type)
     if (allocated(solver_type)) deallocate(solver_type)
     
     ! Recall that Fortran 2003 allows the retrieval of command line 
     ! arguments passed to the code
     CALL get_command_argument(1, arg)
     mcheck(len(trim(arg))>0, USAGE_ERROR_MSG)
     read(arg,*,iostat=istat) NX
     mcheck(istat==0, IO_ERROR_MSG)
     
     CALL get_command_argument(2, arg)
     mcheck(len(trim(arg))>0, USAGE_ERROR_MSG)
     read(arg,*,iostat=istat) NY
     mcheck(istat==0, IO_ERROR_MSG)
     
     CALL get_command_argument(3, arg)
     mcheck(len(trim(arg))>0, USAGE_ERROR_MSG)
     matrix_type = trim(arg)
     
     CALL get_command_argument(4, arg)
     mcheck(len(trim(arg))>0, USAGE_ERROR_MSG)
     solver_type = trim(arg)
     mcheck ( .not. (matrix_type == SPARSE_MATRIX .and. solver_type == DIR_SOLVE), MATRIX_SOLVER_TYPE_ERROR_MSG)     
  end subroutine read_and_check_command_line_parameters
end program example2
