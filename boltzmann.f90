!
! $Id: boltzmann.f90,v 1.10 2001/09/30 10:07:20 jonas Exp $
!

module boltzmann
	use param
	use states
	implicit none
	
	real(ACCU), parameter, private :: k=0.695038769603
	real(ACCU), private :: E0=0.0, T=-1.0
	
contains

! returns partition function 
function partfunc(dvec) result(pf)
	implicit none
	type(dstates), intent(in) :: dvec
	real(ACCU) :: pf

	integer :: i
	real(ACCU) :: dgen, eval, beta
	
	if ( T > 0.0 ) then
		pf=0.0d0
		E0=dvec%state(1)%eval
		beta=k*T
		do i=1, dvec%nstates
			dgen=dble(dvec%state(i)%dgen)
			eval=dvec%state(i)%eval
!            pf=pf+dgen*exp((E0-eval)/beta)
			pf=pf+exp((E0-eval)/beta)
		end do
	else
		pf=1.0
	end if
end function partfunc

! returns statistical weight
function weight(state) result(wgt)
	implicit none
	type(dstate), intent(in) :: state
	real(ACCU) :: wgt

	real(ACCU) :: dgen, eval

	if ( T > 0.0 ) then
		dgen=dble(state%dgen)
		eval=state%eval
!        wgt=dgen*exp((E0-eval)/(k*T))
		wgt=exp((E0-eval)/(k*T))
	else
		wgt=1.0
	end if
end function weight

function weight2(eig, st) result(weight)
	implicit none
    real(ACCU), dimension(:), intent(in) :: eig
    integer(4), intent(in) :: st
	real(ACCU) :: weight

	if ( T > 0.0 ) then
		weight=exp((E0-eig(st))/(k*T))
	else
		weight=1.0
	end if
end function 

subroutine set_temperature(temp)
	implicit none
	real(ACCU), intent(in) :: temp
	T=temp
end subroutine set_temperature

function get_temperature()
	implicit none
	real(ACCU) :: get_temperature
	get_temperature=T
end function get_temperature

end module boltzmann
