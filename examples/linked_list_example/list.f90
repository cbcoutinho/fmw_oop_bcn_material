module list_mod
  use node_list_mod
  implicit none
  
  type list_t
     private
     type(node_list_t), pointer :: head => NULL()
     type(node_list_t), pointer :: tail => NULL()
  contains
    procedure, private :: push_back_integer
    procedure, private :: push_back_real
    procedure, private :: push_back_logical
    procedure, private :: push_back_data
    generic            :: push_back => push_back_integer, &
                                       push_back_real,    &
                                       push_back_logical  
    procedure          :: print
    procedure          :: free
  end type list_t

contains

  subroutine push_back_integer(this, data)
    implicit none
    class(list_t), intent(inout) :: this
    integer      , intent(in)    :: data
    call this%push_back_data(data)
  end subroutine push_back_integer

  subroutine push_back_real(this, data)
    implicit none
    class(list_t), intent(inout) :: this
    real         , intent(in)    :: data
    call this%push_back_data(data)
  end subroutine push_back_real

  subroutine push_back_logical(this, data)
    implicit none
    class(list_t), intent(inout) :: this
    logical      , intent(in)    :: data
    call this%push_back_data(data)
  end subroutine push_back_logical

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

  subroutine print(this)
    implicit none
    class(list_t), intent(in) :: this
    type(node_list_t), pointer :: current
    class(*), pointer :: data
    ! Print contents of link list nodes 
    ! to standard output
    current => this%head
    do while ( associated(current) )
       data => current%get_data()
       select type (data)
       type is (integer)
          write(*,*) data
       type is (logical)
          write(*,*) data
       type is (real)
          write(*,*) data
       end select
       current => current%get_next()
    end do
  end subroutine print

  subroutine free(this)
    implicit none
    class(list_t), intent(inout) :: this
    nullify(this%tail)
    call this%head%free()
    deallocate(this%head)
  end subroutine free

end module list_mod
