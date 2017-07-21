module node_list_mod
  implicit none  
  private 
  
  type node_list_t
     private 
     class(*)         , pointer :: data  => null()
     type(node_list_t), pointer :: next  => null()
   contains
     procedure, non_overridable :: get_data ! Return pointer to data
     procedure, non_overridable :: get_next  ! Return pointer to next
     procedure, non_overridable :: set_next  ! Set next pointer's target
     procedure, non_overridable :: free
  end type node_list_t
  
  interface node_list_t
     module procedure construct_node_list
  end interface node_list_t

  public :: node_list_t
contains

  function get_data(this)
    class(node_list_t), intent(in) :: this
    class(*), pointer :: get_data
    get_data => this%data
  end function get_data

  function get_next(this)
    class(node_list_t), intent(in) :: this
    type(node_list_t), pointer :: get_next
    get_next => this%next
  end function get_next

  subroutine set_next(this,next)
    class(node_list_t), intent(inout) :: this
    type(node_list_t), pointer :: next 
    this%next => next
  end subroutine set_next

  recursive subroutine free(this)
    class(node_list_t), intent(inout) :: this
    if ( associated(this%next) ) then
       call this%next%free()
       deallocate(this%next)
    end if
    deallocate(this%data)
  end subroutine free

  function construct_node_list(data)
    class(*), intent(in) :: data
    type(node_list_t), pointer :: construct_node_list
    allocate(construct_node_list)
    allocate(construct_node_list%data, source=data)
    nullify(construct_node_list%next)
  end function construct_node_list

end module node_list_mod
