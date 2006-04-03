!
! $Id: param.f90,v 1.8 2001/09/30 10:07:20 jonas Exp $ 
! 

!
! Global parameters and access methods
!
module param
	implicit none
	integer, parameter, public :: ACCU=selected_real_kind(12,30)
	integer, parameter, public :: ALL_MF=42
	real(ACCU), parameter, public :: MU_TOLERANCE=1.d-4
	real(ACCU), parameter, public :: DGEN_TOLERANCE=.1d-9

	type rotor
		real(kind=ACCU) :: A, B, C
	end type rotor

	type potential
		real(kind=ACCU) :: Vj, Vk, Vm
	end type potential

	integer, private :: dim
	integer, private :: K_LIMIT=4

contains

function is_linear(mol)
	implicit none
	type(rotor), intent(in) :: mol
	logical :: is_linear

	if (mol%A > 0.0 .or. mol%C > 0.0 ) then
		is_linear=.false.
	else
		is_linear=.true.
	end if
end function is_linear

subroutine set_k_limit(k)
	implicit none
	integer, intent(in) :: k

	if ( k/=0 .and. k/=4 ) then
		stop 'set_k_limit(): Only K=0 or K=4 supported at the moment'
	end if
	K_LIMIT=k
end subroutine set_k_limit

function get_k_limit
	implicit none
	integer :: get_k_limit
	get_k_limit=K_LIMIT
end function get_k_limit

subroutine setdim(i)
	implicit none
	integer, intent(in) :: i
	dim=i
end subroutine setdim

function getdim
	implicit none
	integer :: getdim
	getdim=dim
end function getdim

end module param

