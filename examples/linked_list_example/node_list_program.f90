program node_list_program
  use node_list_mod
  implicit none
  type(node_list_t), pointer :: node
  class(*), pointer :: data
  type(node_list_t), pointer :: current
  integer :: i
  ! Create an integer linked list with 
  ! 10 nodes, and data 1, 2, ..., 10
  ! node becomes the head of the list
  node    => node_list_t(1)
  current => node
  do i=2, 10
     call current%set_next( node_list_t(i) )
     current => current%get_next()
  end do
  ! Print contents of link list nodes 
  ! to standard output
  current => node
  do while ( associated(current) )
     data => current%get_data()
     select type (data)
     type is (integer)
        write(*,*) data
     end select
     current => current%get_next()
  end do
  ! Free all dynamic memory 
  call node%free()
  deallocate(node)
end program node_list_program
