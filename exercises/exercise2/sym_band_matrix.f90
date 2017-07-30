#include "mcheck.i90"
module sym_band_matrix_mod
  use types_mod,  only: IP, RP
  use blas_mod,   only: idamax
  use matrix_mod, only: matrix_t
  use wathen_problem_mod, only: mv_pbu, dpbufa, dpbusl
  implicit none
  private

  public :: sym_band_matrix_t
  type, extends(matrix_t) :: sym_band_matrix_t
    private
    real(rp), allocatable :: a(:,:)
    integer(ip) :: ml  = -1
    integer(ip) :: mu  = -1
    integer(ip) :: lda = -1
  contains
    procedure :: create    => sym_band_matrix_create
    procedure :: assembly  => sym_band_matrix_assembly
    procedure :: apply     => sym_band_matrix_apply
    procedure :: factorize => sym_band_matrix_factorize
    procedure :: backsolve => sym_band_matrix_backsolve
    procedure :: free      => sym_band_matrix_free
  end type sym_band_matrix_t

contains

  subroutine sym_band_matrix_create(this,n,ml,mu,nz)
    class(sym_band_matrix_t) , intent(inout) :: this
    integer(ip)          , intent(in)    :: n
    integer(ip), optional, intent(in)    :: ml,mu,nz
    call this%free()
    call this%set_size(n)
    mcheck(present(ml),"ml dummy argument required by sym_band_matrix_create")
    mcheck(present(mu),"mu dummy argument required by sym_band_matrix_create")
    this%ml = ml
    this%mu = mu
    this%lda = mu+1
    allocate ( this%a(1:this%lda,1:n) )
    this%a = 0.0_rp
  end subroutine sym_band_matrix_create

  subroutine sym_band_matrix_assembly(this,i,j,a)
    implicit none
    class(sym_band_matrix_t), intent(inout) :: this
    integer(ip)         , intent(in)    :: i,j
    real(rp)            , intent(in)    :: a
    if (i <= j) then
      this%a(i-j+this%mu+1,j) = this%a(i-j+this%mu+1,j) + a
    endif
  end subroutine sym_band_matrix_assembly

  subroutine sym_band_matrix_apply(this,x,y)
    implicit none
    class(sym_band_matrix_t), intent(in)    :: this
    real(rp)            , intent(in)    :: x(:)
    real(rp)            , intent(inout) :: y(:)
    call mv_pbu ( this%get_size(), this%get_size(), this%mu, this%a, x, y)
  end subroutine sym_band_matrix_apply

  subroutine sym_band_matrix_factorize(this, factors, pivots)
    implicit none
    class(sym_band_matrix_t)        , intent(in)    :: this
    class(matrix_t), allocatable, intent(inout) :: factors
    integer(ip)                 , intent(inout) :: pivots(:)
    integer(ip) :: info, lda

    if (allocated(factors)) then
      call factors%free()
      deallocate(factors)
    end if

    allocate(sym_band_matrix_t :: factors)
    select type(factors)
    type is(sym_band_matrix_t)
      call factors%create(this%get_size(),this%ml,this%mu)
      factors%a = this%a
      call dpbufa ( this%get_size(), this%mu, factors%a, info )
      mcheck(info==0, 'Error in band factorization')
      class default
    end select
  end subroutine sym_band_matrix_factorize

  subroutine sym_band_matrix_backsolve(this,pivots,rhs,x)
    class(sym_band_matrix_t), intent(in)    :: this
    integer(ip)         , intent(in)    :: pivots(:)
    real(rp)            , intent(in)    :: rhs(:)
    real(rp)            , intent(inout) :: x(:)
    x = rhs
    call dpbusl ( this%get_size(), this%mu, this%a, x)
  end subroutine sym_band_matrix_backsolve

  subroutine sym_band_matrix_free(this)
    class(sym_band_matrix_t), intent(inout) :: this
    call this%set_size(-1)
    this%ml = -1
    this%mu = -1
    this%lda = -1
    if (allocated(this%a)) deallocate(this%a)
  end subroutine sym_band_matrix_free

end module sym_band_matrix_mod
