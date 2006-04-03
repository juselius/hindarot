!
! $Id: states.f90,v 1.2 2001/09/20 17:39:40 jonas Exp $
! 

module states
	use param
	implicit none

	type dstate
		real(ACCU) :: eval
		integer :: dgen
	end type dstate

	type dstates
		integer :: nstates
		type(dstate), dimension(:), pointer :: state
	end type

	public init_dstates, del_dstates, dstates, dstate

	private
contains

subroutine init_dstates(dvec)
	implicit none
	type(dstates), intent(out) :: dvec
	
	integer :: dim
	dim=getdim()
	allocate(dvec%state(dim))
end subroutine init_dstates

subroutine del_dstates(dvec)
	implicit none
	type(dstates), intent(out) :: dvec
	
	if (associated(dvec%state)) then
		deallocate(dvec%state)
	end if
end subroutine del_dstates

end module states

