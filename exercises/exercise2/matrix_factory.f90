#include "mcheck.i90"
module matrix_factory_mod
  use types_mod
  use matrix_mod
  use dense_matrix_mod
  use band_matrix_mod
  use sparse_matrix_mod

  implicit none
  private

  character(*), parameter :: FULL_MATRIX   = "full_matrix"
  character(*), parameter :: BAND_MATRIX   = "band_matrix"
  character(*), parameter :: SYM_BAND_MATRIX   = "sym_band_matrix"
  character(*), parameter :: SPARSE_MATRIX = "sparse_matrix"
  character(*), parameter :: MATRIX_TYPE_ERROR_MSG = "matrix_type MUST BE either [full_matrix|band_matrix|sym_band_matrix|sparse_matrix]"

  public :: matrix_factory
  public :: FULL_MATRIX, BAND_MATRIX, SYM_BAND_MATRIX, SPARSE_MATRIX
contains

  subroutine matrix_factory( matrix_type, matrix )
    implicit none
    character(len=*)            , intent(in)    :: matrix_type
    class(matrix_t), allocatable, intent(inout) :: matrix

    if ( allocated(matrix) ) then
      call matrix%free(); deallocate(matrix)
    end if

    mcheck(trim(matrix_type) == FULL_MATRIX .or. trim(matrix_type) == BAND_MATRIX .or. trim(matrix_type) == SYM_BAND_MATRIX .or. trim(matrix_type) == SPARSE_MATRIX, MATRIX_TYPE_ERROR_MSG)

    select case ( trim(matrix_type) )
    case (FULL_MATRIX)
      allocate( dense_matrix_t :: matrix )
    case (BAND_MATRIX)
      allocate( band_matrix_t :: matrix )
    ! case (SYM_BAND_MATRIX)
    !   allocate( sym_band_matrix_t :: matrix )
    case (SPARSE_MATRIX)
      allocate( sparse_matrix_t :: matrix )
    end select
  end subroutine matrix_factory

end module matrix_factory_mod
