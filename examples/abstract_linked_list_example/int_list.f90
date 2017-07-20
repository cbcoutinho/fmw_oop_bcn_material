module int_list_mod
 use list_mod
 implicit none
 private 
  
 type, extends(list_t) :: int_list_t
  private
 contains
  ! Child's TBPs to push back a new integer
  ! and get current's node list integer 
  procedure :: push_back
  procedure :: get_current

  ! TBP overriding parent's deferred method
  procedure :: print
 end type int_list_t

 public :: int_list_t

contains
 subroutine push_back(this,data)
  implicit none
  class(int_list_t), intent(inout) :: this
  integer          , intent(in)    :: data
  ! Call parent's unlimited polymorphic 
  ! variant of push_back
  call this%push_back_data(data)
 end subroutine push_back
  
 function get_current(this)
  implicit none
  class(int_list_t), intent(in) :: this
  integer :: get_current
  class(*), pointer :: data
  ! Call parent's unlimited polymorphic
  ! variant of get_current
  data => this%get_current_data()
  select type(data)
  type is(integer)
    get_current = data
  end select
 end function get_current

 subroutine print(this)
  implicit none
  class(int_list_t), intent(inout) :: this
  
  call this%first()
  do while (.not. this%has_finished() )
     write(*,*) this%get_current()
     call this%next()
  end do
 end subroutine print

end module int_list_mod
