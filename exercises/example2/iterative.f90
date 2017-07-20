module iterative_solver_names
  use types
  use matrix_names
  use solver_names

  implicit none
  private
  real(rp), parameter :: tolerance =1.0e-8_rp
  
  type, extends(solver_t) :: cg_solver_t
     private
     integer(ip) :: n
     real(rp) , allocatable :: p(:),r(:),Ap(:)
   contains
     procedure  :: create => cg_solver_create
     procedure  :: apply  => cg_solver_apply
     procedure  :: free   => cg_solver_free
  end type cg_solver_t

  public :: cg_solver_t

contains

  subroutine cg_solver_create(this,matrix) 
    implicit none
    class(cg_solver_t), intent(inout) :: this
    class(matrix_t)   , intent(in)    :: matrix
    this%n = matrix%get_size()
    allocate ( this%p(1:this%n) )
    allocate ( this%r(1:this%n) )
    allocate ( this%Ap(1:this%n) )
  end subroutine cg_solver_create

  subroutine cg_solver_free(this) 
    implicit none
    class(cg_solver_t), intent(inout) :: this
    deallocate ( this%p )
    deallocate ( this%r )
    deallocate ( this%Ap )
  end subroutine cg_solver_free

  subroutine cg_solver_apply(this,matrix,rhs,x) 
    implicit none
    class(cg_solver_t), intent(inout) :: this
    class(matrix_t)   , intent(in)    :: matrix
    real(rp)          , intent(in)    :: rhs(:)
    real(rp)          , intent(inout) :: x(:)
    real(rp)          , allocatable   :: y(:)
    real(rp)    :: alpha, beta, rsqrd, pap
    integer(ip) :: it
    logical     :: converged

    call matrix%apply(x,this%Ap)
    this%r = rhs - this%Ap
    rsqrd  = dot_product ( this%r, this%r)
    converged = (sqrt(rsqrd) < tolerance) 
    this%p = this%r
    it = 1
    do while( (.not.converged) .or. (it<matrix%get_size()))
       call matrix%apply(this%p,this%Ap)
       rsqrd  = dot_product ( this%r, this%r)
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
       converged = (sqrt(rsqrd) < tolerance) 
       beta = beta * rsqrd
       this%p = this%r + beta * this%ap
       it = it + 1
    end do
    
  end subroutine cg_solver_apply
  
end module iterative_solver_names
