program list_program
  use list_mod
  implicit none
  
  type(list_t) :: list

  ! Create an heterogeneous data type linked list
  call list%push_back(3)
  call list%push_back(4.45)
  call list%push_back(.true.)
  
  ! Print contents of link list nodes 
  ! to standard output
  call list%print()
 
  ! Free all dynamic memory 
  call list%free()
end program list_program
