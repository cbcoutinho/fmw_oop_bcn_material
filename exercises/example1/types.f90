module types
  use, intrinsic :: iso_fortran_env, only: INT8, INT32, INT64, REAL32, REAL64
  implicit none
  integer, parameter :: ip   = INT32  ! Integer precision
  integer, parameter :: rp   = REAL64 ! Real precision
  integer, parameter :: drp  = REAL64 ! Double Real Precision
  integer, parameter :: dip  = INT64  ! Double Integer Precision
end module types
