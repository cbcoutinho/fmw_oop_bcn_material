#define mcheck(test,message) \
   if (.not.(test)) then ; \
      write(0,'(a,a,a,i10)') "Check failed in file ", __FILE__,", at line number", __LINE__  ; \
      write(0,'(a,a)') "Cause: ", message; \
      stop -1 ; \
   endif   

program example1
  use types
  use random_numbers_mod
  use wathen_problem_mod
  use direct_solver_mod
  use iterative_solver_mod
  implicit none
  
  character(:), parameter :: FULL_MATRIX   = "full_matrix"
  character(:), parameter :: BAND_MATRIX   = "band_matrix"
  character(:), parameter :: SPARSE_MATRIX = "sparse_matrix"
  character(:), parameter :: DIR_SOLVE     = "dir_solve"
  character(:), parameter :: CG_SOLVE      = "cg_solve"
  
  character(:), allocatable :: matrix_type, solver_type
  integer (IP) :: nz_num
  integer (IP), allocatable :: row(:)
  integer (IP), allocatable :: col(:)
  real (RP), allocatable :: a_st(:)
  integer(IP) :: lda
  integer(IP) :: ml
  integer(IP) :: mu
  integer(IP) :: md

  real (RP), allocatable :: a(:,:)
  real (RP), allocatable :: b(:)
  real (RP) :: e
  integer (IP) :: i
  integer (IP) :: info
  integer (IP), allocatable :: ipvt(:)
  integer (IP) :: job
  integer (IP) :: n
  integer (IP) :: nx
  integer (IP) :: ny
  integer (IP) :: seed
  real (RP), allocatable :: x1(:)
  real (RP), allocatable :: x2(:)

  ! Read from command line:
  call read_and_check_command_line_parameters(  NX, NY, matrix_type, solver_type ) 

  ! Write information
  write ( *, '(a,i6)' ) '  Elements in X direction NX = ', nx
  write ( *, '(a,i6)' ) '  Elements in Y direction NY = ', ny
  write ( *, '(a,i6)' ) '  Number of elements = ', nx * ny
  write ( *, '(a,a)' ) '  Matrix type = ', matrix_type
  write ( *, '(a,a)' ) '  Solver type = ', solver_type

  ! Compute the number of unknowns.
  call wathen_order ( nx, ny, n )
  write ( *, '(a,i6)' ) '  Number of nodes N = ', n

  ! Set up a random solution X.
  seed = 123456789
  allocate ( x1(1:n) )
  call r8vec_uniform_01 ( n, seed, x1 )

  seed = 123456789
  if(matrix_type == FULL_MATRIX) then
     ! Compute the matrix.
     allocate ( a(1:n,1:n) )
     call wathen_ge ( nx, ny, n, seed, a )
     !  Compute the corresponding right hand side B.
     allocate ( b(1:n) )
     b = matmul ( a, x1 )
     !  Solve the linear system.
     allocate ( x2(1:n) )
     if ( solver_type == DIR_SOLVE ) then
        allocate ( ipvt(1:n) )
        call dgefa ( a, n, n, ipvt, info ) ! Factorization
        x2 = b
        job = 0
        call dgesl ( a, n, n, ipvt, x2, job ) ! Back sub
        deallocate ( ipvt )
     else if( solver_type == CG_SOLVE ) then
        x2 = 1.0_rp
        call cg_ge ( n, a, b, x2 )
     end if
     deallocate ( a )
  else if(matrix_type == BAND_MATRIX) then
     call wathen_bandwidth ( nx, ny, ml, md, mu )
     ! Compute the matrix.
     allocate ( a(1:2*ml+mu+1,1:n) )
     call wathen_gb ( nx, ny, n, seed, a )
     !  Compute the corresponding right hand side B.
     allocate ( b(1:n) )
     call mv_gb ( n, n, ml, mu, a, x1, b )
     !  Solve the linear system.
     allocate ( x2(1:n) )
     if ( solver_type == DIR_SOLVE ) then
        lda = 2 * ml + mu + 1
        allocate ( ipvt(1:n) )
        call dgbfa ( a, lda, n, ml, mu, ipvt, info )
        x2 = b
        job = 0
        call dgbsl ( a, lda, n, ml, mu, ipvt, x2, job )
        deallocate ( ipvt )
     else
        x2 = 1.0_rp
        call cg_gb ( n, ml, mu, a, b, x2 )
     end if
     deallocate ( a )
  else if(matrix_type == SPARSE_MATRIX) then
     !  Compute the matrix size.
     call wathen_st_size ( nx, ny, nz_num )
     !  Compute the matrix.
     seed = 123456789
     allocate ( row(1:nz_num) )
     allocate ( col(1:nz_num) )
     allocate ( a_st(1:nz_num) )
     call wathen_st ( nx, ny, nz_num, seed, row, col, a_st )
     !  Compute the corresponding right hand side B.
     allocate ( b(1:n) )
     call mv_st ( n, n, nz_num, row, col, a_st, x1, b )
     !  Solve the linear system.
     allocate ( x2(1:n) )
     x2 = 1.0_rp
     call cg_st ( n, nz_num, row, col, a_st, b, x2 )
     deallocate ( a_st )
  end if

  !  Compute the maximum solution error.
  e = maxval ( abs ( x1 - x2 ) )
  write ( *, '(a,g14.6)' ) '  Maximum solution error is ', e 

  !  Free memory.
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
     
     character(:), parameter :: USAGE_ERROR_MSG       = "Usage: example1 NX NY matrix_type solver_type"
     character(:), parameter :: MATRIX_TYPE_ERROR_MSG = "matrix_type MUST BE either [full_matrix|band_matrix|sparse_matrix]"
     character(:), parameter :: SOLVER_TYPE_ERROR_MSG = "solver_type MUST BE either [dir_solve|cg_solve]"
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
     mcheck(trim(arg) == FULL_MATRIX .or. trim(arg) == BAND_MATRIX .or. trim(arg) == SPARSE_MATRIX, MATRIX_TYPE_ERROR_MSG)
     matrix_type = trim(arg)
     
     CALL get_command_argument(4, arg)
     mcheck(len(trim(arg))>0, USAGE_ERROR_MSG)
     mcheck(trim(arg) == DIR_SOLVE .or. trim(arg) == CG_SOLVE, SOLVER_TYPE_ERROR_MSG)
     solver_type = trim(arg)
     mcheck ( .not. (matrix_type == SPARSE_MATRIX .and. solver_type == DIR_SOLVE), MATRIX_SOLVER_TYPE_ERROR_MSG)     
   end subroutine read_and_check_command_line_parameters

end program example1
