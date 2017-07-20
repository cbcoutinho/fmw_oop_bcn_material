program list_program
  use int_list_mod
  implicit none
  
  type(int_list_t) :: list
  ! type(list_t) :: list ! ILLEGAL !!!!

  ! Fill integer data type linked list
  call list%push_back(1)
  call list%push_back(2)
  call list%push_back(3)
  
  ! Print contents of link list nodes 
  ! to standard output used overrided
  ! version in int_list_t
  call list%print()
 
  ! Free all dynamic memory 
  call list%free()
end program list_program
