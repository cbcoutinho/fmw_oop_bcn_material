module direct_solver_names
  use types
  use matrix_names
  use solver_names
  implicit none
  private

  type, extends(solver_t) :: direct_solver_t
     private
     class(matrix_t), allocatable :: factors
     integer(ip)    , allocatable :: pivots(:)
     real(rp)       , allocatable :: work(:)
   contains
     procedure  :: create => direct_solver_create
     procedure  :: apply  => direct_solver_apply
  end type direct_solver_t

  public :: direct_solver_t

contains

  subroutine direct_solver_create(this,matrix) 
    implicit none
    class(direct_solver_t), intent(inout) :: this
    class(matrix_t)       , intent(in)    :: matrix
    allocate(this%pivots(matrix%get_size()))
    allocate(this%work(matrix%get_size()))
    call matrix%factorize(this%factors,this%pivots)
  end subroutine direct_solver_create

  subroutine direct_solver_apply(this,matrix,rhs,x) 
    implicit none
    class(direct_solver_t), intent(inout) :: this
    class(matrix_t)       , intent(in)    :: matrix
    real(rp)              , intent(in)    :: rhs(:)
    real(rp)              , intent(inout) :: x(:)
    real(rp)              , allocatable   :: y(:)
    if(this%factors%get_size()/=matrix%get_size()) then
       ! Write an error message
    else
       ! The algorithm will work
    end if
    call this%factors%backsolve(this%pivots,rhs,x)
    ! Check the solution
    call matrix%apply(x,this%work)
    if(maxval(abs(this%work-rhs))>=1.0e-8_rp) then
       ! write an error message and stop
    end if
  end subroutine direct_solver_apply
  
end module direct_solver_names
