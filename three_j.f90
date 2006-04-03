!
! $Id: three_j.f90,v 1.3 2001/09/28 06:14:35 jonas Exp $
!

module three_j
	use param
	implicit none

	public j3coef
	
	private
contains

! 3-j coefficients
function j3coef(a,b,c,al,be,ga)
	implicit none
	integer, intent(in) :: a, b, c, al, be, ga
	real (kind=ACCU) :: j3coef
	
	integer :: p, nyy, factor1, sigma
	logical :: flag
	real (kind=ACCU) :: factor2, factor3, factor4
	
	factor1=(-1.d0)**(a-b-ga)
	factor4=0.0d0
	sigma=a+b+c 
	
	if ((al+be) /= -ga) then
		j3coef=0.0d0
		return
	end if
	
	if (al == 0 .and. be == 0 .and. ga == 0 .and. mod(sigma,2) /= 0) then
		j3coef=0.0d0
		return
	endif
	
	flag=.false.
	do p=abs(a-b), (a+b)
		if (c == p) then
			flag=.true.
			exit
		end if
	end do
	 
	if ( .not. flag  ) then
		j3coef=0.0d0
		return
	endif

	factor2=sqrt(1.0d0*ker(a+b-c)*ker(a+c-b)*ker(b+c-a)/ker(a+b+c+1))
	factor3=sqrt(1.0d0*ker(a+al)*ker(a-al)*ker(b+be)*ker(b-be)*&
				   &ker(c+ga)*ker(c-ga))

	do nyy=0,15
		if (a-al-nyy < 0 .or. c-b+al+nyy < 0 .or. b+be-nyy < 0 .or. &
			 	& c-a-be+nyy < 0 .or. a+b-c-nyy < 0) then
			continue
		else
			factor4=factor4+((-1)**nyy)*(1.0d0/(ker(a-al-nyy)*ker(c-b+al+nyy)*&
			 		&ker(b+be-nyy)*ker(c-a-be+nyy)*ker(nyy)*ker(a+b-c-nyy)))
		end if
	end do

	j3coef=factor1*factor2*factor3*factor4
end function j3coef

! factorials
function ker(z) 
	implicit none
	integer, intent(in) :: z
	real (kind=ACCU) :: ker
	
	integer :: n
	
	if (z == 0) then
		ker=1
	else
		ker=1
		do n=1,z
			ker=ker*n
		end do
	end if
end function ker
end module three_j
