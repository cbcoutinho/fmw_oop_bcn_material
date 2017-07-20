module list_mod
  use node_list_mod
  implicit none
  private 
  
  type, abstract :: list_t
   private
   type(node_list_t), pointer :: head   => NULL()
   type(node_list_t), pointer :: tail   => NULL()
   type(node_list_t), pointer :: cursor => NULL()
  contains
   ! list_t construct and destruct TBPs
   procedure, non_overridable :: push_back_data
   procedure, non_overridable :: free
    
   ! list_t traversal TBPs
   procedure, non_overridable :: first
   procedure, non_overridable :: next
   procedure, non_overridable :: has_finished
   procedure, non_overridable :: get_current_data
   
   ! Deferred TBP to print list_t to stdout
   procedure(print_list), deferred :: print
  end type list_t

  abstract interface
   subroutine print_list(this)
    import :: list_t
    class(list_t), intent(inout) :: this
   end subroutine print_list
  end interface

  public :: list_t

contains
  subroutine push_back_data(this,data)
   implicit none
   class(list_t), intent(inout) :: this
   class(*)     , intent(in)    :: data
   if ( .not. associated(this%head)) then
      this%head => node_list_t(data)
      this%tail => this%head
   else
      call this%tail%set_next(node_list_t(data))
      this%tail => this%tail%get_next()
   end if
  end subroutine push_back_data

  subroutine free(this)
   implicit none
   class(list_t), intent(inout) :: this
   nullify(this%tail)
   call this%head%free()
   deallocate(this%head)
  end subroutine free

  subroutine first(this)
   implicit none
   class(list_t), intent(inout) :: this
   this%cursor => this%head
  end subroutine first
  
  subroutine next(this)
   implicit none
   class(list_t), intent(inout) :: this
   this%cursor => this%cursor%get_next()
  end subroutine next

  function has_finished(this)
   implicit none
   class(list_t), intent(in) :: this
   logical :: has_finished
   has_finished = .not. associated(this%cursor)
  end function has_finished
  
  function get_current_data(this)
   implicit none
   class(list_t), intent(in) :: this
   class(*), pointer :: get_current_data
   get_current_data => this%cursor%get_data()
  end function get_current_data
end module list_mod
