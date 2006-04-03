!
! $Id: glue.f90,v 1.1 2001/09/25 16:32:23 jonas Exp $
!

!
! These two functions have to be in their own module to sort out silly
! module interdependencies...
!

module glue
	use param
	use states
	implicit none
	integer, public :: NPRINT=-1
	integer, public :: NMOM=-1

contains

! Returns number of states to _actually_ print
function n_print(dvec)
	implicit none
	type(dstates), intent(in) :: dvec
	integer :: n_print
	
	if ( NPRINT == -1 ) then
		n_print=dvec%nstates
	else if ( NPRINT > dvec%nstates ) then
		n_print=dvec%nstates
	else
		n_print=NPRINT
	end if
end function n_print

! Returns number of transition moments to _actually_ calculate
function n_moments(dvec)
	implicit none
	type(dstates), intent(in) :: dvec
	integer :: n_moments
	
	if ( NMOM == -1 ) then
		n_moments=n_print(dvec)
	else if ( NMOM > dvec%nstates ) then
		n_moments=dvec%nstates
	else
		n_moments=NMOM
	end if
end function n_moments

end module glue
