#include "mcheck.i90"
module direct_solver_mod
  use types_mod,  only: IP, RP
  use matrix_mod, only: matrix_t
  use solver_mod, only: solver_t
  implicit none
  private

  public :: direct_solver_t
  type, extends(solver_t) :: direct_solver_t
     private
     class(matrix_t), allocatable :: factors
     integer(ip)    , allocatable :: pivots(:)
     real(rp)       , allocatable :: work(:)
   contains
     procedure  :: create => direct_solver_create
     procedure  :: solve  => direct_solver_solve
     procedure  :: free   => direct_solver_free
  end type direct_solver_t


contains

  subroutine direct_solver_create(this,matrix)
    implicit none
    class(direct_solver_t), intent(inout) :: this
    class(matrix_t)       , intent(in)    :: matrix
    allocate(this%pivots(matrix%get_size()))
    allocate(this%work(matrix%get_size()))
    call matrix%factorize(this%factors,this%pivots)
    call this%set_matrix(matrix)
  end subroutine direct_solver_create

  subroutine direct_solver_solve(this,rhs,x)
    implicit none
    class(direct_solver_t), intent(inout) :: this
    real(rp)              , intent(in)    :: rhs(:)
    real(rp)              , intent(inout) :: x(:)
    call this%factors%backsolve(this%pivots,rhs,x)
  end subroutine direct_solver_solve

  subroutine direct_solver_free(this)
    implicit none
    class(direct_solver_t), intent(inout) :: this
    if (allocated(this%factors)) then
      call this%factors%free()
      deallocate(this%factors)
    end if
    if (allocated(this%pivots)) deallocate(this%pivots)
    if (allocated(this%work)) deallocate(this%work)
    call this%nullify_matrix()
  end subroutine direct_solver_free

end module direct_solver_mod
