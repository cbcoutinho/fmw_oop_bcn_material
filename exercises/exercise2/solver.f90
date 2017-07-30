module solver_mod
  use types_mod,  only: IP, RP
  use matrix_mod, only: matrix_t
  implicit none
  private

  public :: solver_t
  type, abstract :: solver_t
     private
     class(matrix_t), pointer :: matrix => NULL()
   contains
     procedure, non_overridable :: set_matrix
     procedure, non_overridable :: get_matrix
     procedure, non_overridable :: nullify_matrix
     procedure(create_interface), deferred :: create
     procedure(solve_interface) , deferred :: solve
     procedure(free_interface)  , deferred :: free
  end type solver_t

  abstract interface
     subroutine create_interface(this,matrix)
       import :: solver_t, matrix_t
       implicit none
       class(solver_t), intent(inout) :: this
       class(matrix_t), intent(in)    :: matrix
     end subroutine create_interface
     subroutine solve_interface(this,rhs,x)
       import :: solver_t, matrix_t, rp
       implicit none
       class(solver_t), intent(inout) :: this
       real(rp)       , intent(in)    :: rhs(:)
       real(rp)       , intent(inout) :: x(:)
     end subroutine solve_interface
     subroutine free_interface(this)
       import :: solver_t
       implicit none
       class(solver_t), intent(inout) :: this
     end subroutine free_interface
  end interface

contains
  subroutine set_matrix(this,matrix)
    implicit none
    class(solver_t)        , intent(inout) :: this
    class(matrix_t), target, intent(in)    :: matrix
    this%matrix => matrix
  end subroutine set_matrix

  function get_matrix(this)
    implicit none
    class(solver_t), intent(inout) :: this
    class(matrix_t), pointer :: get_matrix
    get_matrix => this%matrix
  end function get_matrix

   subroutine nullify_matrix(this)
    implicit none
    class(solver_t), intent(inout) :: this
    nullify(this%matrix)
  end subroutine nullify_matrix
end module solver_mod
