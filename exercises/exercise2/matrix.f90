module matrix_mod
  use types_mod
  implicit none
  private

  type, abstract :: matrix_t
     private
     integer(IP) :: n = -1
   contains
     procedure, non_overridable :: set_size => matrix_set_size
     procedure, non_overridable :: get_size => matrix_get_size
     procedure(create_interface)   , deferred :: create
     procedure(assembly_interface) , deferred :: assembly
     procedure(apply_interface)    , deferred :: apply
     procedure(factorize_interface), deferred :: factorize
     procedure(backsolve_interface), deferred :: backsolve
     procedure(free_interface)     , deferred :: free
  end type matrix_t

  abstract interface
     subroutine create_interface(this,n,ml,mu,nz)
       import :: matrix_t, ip
       class(matrix_t)      , intent(inout) :: this
       integer(ip)          , intent(in)    :: n
       integer(ip), optional, intent(in)    :: ml,mu,nz
     end subroutine create_interface
     subroutine assembly_interface(this,i,j,a) 
       import :: matrix_t, rp, ip
       implicit none
       class(matrix_t), intent(inout) :: this
       integer(ip)    , intent(in)    :: i,j
       real(rp)       , intent(in)    :: a
     end subroutine assembly_interface
     subroutine apply_interface(this,x,y) 
       import :: matrix_t, rp
       implicit none
       class(matrix_t), intent(in)    :: this
       real(rp)       , intent(in)    :: x(:)
       real(rp)       , intent(inout) :: y(:)
     end subroutine apply_interface      
     subroutine factorize_interface(this,factors,pivots) 
       import :: matrix_t, ip
       implicit none
       class(matrix_t), intent(in)    :: this
       class(matrix_t), allocatable, intent(inout) :: factors
       integer(ip)    , intent(inout) :: pivots(:)
     end subroutine factorize_interface
     subroutine backsolve_interface(this,pivots,rhs,x) 
       import :: matrix_t, ip, rp
       class(matrix_t), intent(in)    :: this
       integer(ip)    , intent(in)    :: pivots(:)
       real(rp)       , intent(in)    :: rhs(:)
       real(rp)       , intent(inout) :: x(:)
     end subroutine backsolve_interface
     subroutine free_interface(this) 
       import :: matrix_t, ip, rp
       class(matrix_t), intent(inout) :: this
     end subroutine free_interface
  end interface

  public :: matrix_t

contains

  subroutine matrix_set_size(this,n)
    class(matrix_t), intent(inout) :: this
    integer(IP)    , intent(in)    :: n
    this%n=n
  end subroutine matrix_set_size

  function matrix_get_size(this)
    class(matrix_t), intent(in) :: this
    integer(IP) :: matrix_get_size
    matrix_get_size = this%n
  end function matrix_get_size

end module matrix_mod
