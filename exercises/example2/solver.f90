module solver_names
  use types
  use matrix_names
  implicit none
  private

  type, abstract :: solver_t
     private
   contains
     procedure(create_interface), deferred :: create
     procedure(apply_interface) , deferred :: apply
  end type solver_t

  abstract interface
     subroutine create_interface(this,matrix) 
       import :: solver_t, matrix_t
       implicit none
       class(solver_t), intent(inout) :: this
       class(matrix_t), intent(in)    :: matrix
     end subroutine create_interface
     subroutine apply_interface(this,matrix,rhs,x) 
       import :: solver_t, matrix_t, rp
       implicit none
       class(solver_t), intent(inout) :: this
       class(matrix_t), intent(in)    :: matrix
       real(rp)       , intent(in)    :: rhs(:)
       real(rp)       , intent(inout) :: x(:)
     end subroutine apply_interface      
  end interface
  public :: solver_t

end module solver_names
