module types
  use, intrinsic :: iso_fortran_env, only: INT8, INT32, INT64, REAL32, REAL64
  implicit none
  integer, parameter :: IP   = INT32  ! Integer precision
  integer, parameter :: RP   = REAL64 ! Real precision
  integer, parameter :: DRP  = REAL64 ! Double Real Precision
  integer, parameter :: DIP  = INT64  ! Double Integer Precision
end module types
