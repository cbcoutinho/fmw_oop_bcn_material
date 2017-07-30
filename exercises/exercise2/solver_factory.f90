#include "mcheck.i90"
module solver_factory_mod
  use types_mod,            only: IP, RP
  use solver_mod,           only: solver_t
  use iterative_solver_mod, only: cg_solver_t
  use direct_solver_mod,    only: direct_solver_t

  implicit none
  private

  character(*), parameter :: DIR_SOLVE     = "dir_solve"
  character(*), parameter :: CG_SOLVE      = "cg_solve"
  character(*), parameter :: SOLVER_TYPE_ERROR_MSG = "solver_type MUST BE either [dir_solve|cg_solve]"

  public :: solver_factory
  public :: DIR_SOLVE, CG_SOLVE
contains

  subroutine solver_factory( solver_type, solver )
    implicit none
    character(len=*)            , intent(in)    :: solver_type
    class(solver_t), allocatable, intent(inout) :: solver

    if ( allocated(solver) ) then
      call solver%free(); deallocate(solver)
    end if

    mcheck(trim(solver_type) == DIR_SOLVE .or. trim(solver_type) == CG_SOLVE, SOLVER_TYPE_ERROR_MSG)

    select case ( trim(solver_type) )
    case (DIR_SOLVE)
      allocate( direct_solver_t :: solver )
    case (CG_SOLVE)
      allocate( cg_solver_t :: solver )
    end select
  end subroutine solver_factory

end module solver_factory_mod
