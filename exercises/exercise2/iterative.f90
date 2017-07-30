module iterative_solver_mod
  use types_mod,  only: IP, RP
  use matrix_mod, only: matrix_t
  use solver_mod, only: solver_t

  implicit none
  private
  real(RP), parameter :: TOLERANCE = 1.0e-8_rp

  public :: cg_solver_t
  type, extends(solver_t) :: cg_solver_t
     private
     integer(IP) :: n = -1
     real(RP), allocatable :: p(:),r(:),Ap(:)
   contains
     procedure  :: create => cg_solver_create
     procedure  :: solve  => cg_solver_solve
     procedure  :: free   => cg_solver_free
  end type cg_solver_t


contains

  subroutine cg_solver_create(this,matrix)
    implicit none
    class(cg_solver_t), intent(inout) :: this
    class(matrix_t)   , intent(in)    :: matrix
    call this%free()
    this%n = matrix%get_size()
    allocate ( this%p(1:this%n) )
    allocate ( this%r(1:this%n) )
    allocate ( this%Ap(1:this%n) )
    call this%set_matrix(matrix)
  end subroutine cg_solver_create

  subroutine cg_solver_solve(this,rhs,x)
    implicit none
    class(cg_solver_t), intent(inout) :: this
    real(RP)          , intent(in)    :: rhs(:)
    real(RP)          , intent(inout) :: x(:)
    real(RP)          , allocatable   :: y(:)
    real(RP)    :: alpha, beta, rsqrd, pap
    integer(IP) :: it
    logical     :: converged
    class(matrix_t), pointer :: matrix

    matrix => this%get_matrix()

    call matrix%apply(x,this%Ap)
    this%r = rhs - this%Ap
    rsqrd  = dot_product ( this%r, this%r)
    converged = (sqrt(rsqrd) < TOLERANCE)
    this%p = this%r
    it = 1
    do while( (.not.converged) .or. (it<matrix%get_size()))
       call matrix%apply(this%p,this%Ap)
       pap    = dot_product ( this%p, this%Ap)
       if(pap==0.0_rp) then
          alpha = 0.0
       else
          alpha = rsqrd / pap
       end if
       x      = x      + alpha * this%p
       this%r = this%r - alpha * this%ap
       beta   = 1.0_rp / rsqrd
       rsqrd  = dot_product ( this%r, this%r)
       converged = (sqrt(rsqrd) < TOLERANCE)
       beta = beta * rsqrd
       this%p = this%r + beta * this%p
       it = it + 1
    end do
  end subroutine cg_solver_solve

  subroutine cg_solver_free(this)
    implicit none
    class(cg_solver_t), intent(inout) :: this
    this%n = -1
    if (allocated(this%p)) deallocate(this%p)
    if (allocated(this%r)) deallocate(this%r)
    if (allocated(this%Ap)) deallocate(this%Ap)
    call this%nullify_matrix()
  end subroutine cg_solver_free

end module iterative_solver_mod
