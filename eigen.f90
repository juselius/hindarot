!
! $Id: eigen.f90,v 1.10 2001/09/30 10:20:33 jonas Exp $
!

module eigen
	use param
	use states
	use externals
	implicit none
	public diag, sort_eig, degeneracy 
!	public diag2

	private
contains

! diag() does the diagonalization, and returns the eigenvalues in eigen
subroutine diag(hmat, eig)
	implicit none
	real(ACCU), dimension(:,:), intent(inout) :: hmat
	real(ACCU), dimension(:), intent(out) :: eig
	
	integer, parameter :: SORT_EV=1   ! sort eigenvalues
	integer, parameter :: EVAL_ONLY=0 ! want eigenvectors as well...
	integer :: i,j, dim
	real(ACCU), dimension(:), allocatable :: hpack, evmat

	dim=getdim()

	allocate(hpack(dim**2))
	allocate(evmat(dim**2))
	
	print '(a)', '*** entering packhm'
	call packhm(hmat, hpack, dim)
	print '(a)', '*** entering deigen'
	call deigen(hpack, evmat, dim, EVAL_ONLY, SORT_EV)
	hmat=reshape(evmat, (/dim,dim/))
	deallocate(evmat)

! unpack eigenvalues (black magick), see deigen.f for details
	do i=1,dim
	  	eig(i)=hpack(i*(i-1)/2+i)
	end do
	deallocate(hpack)
end subroutine diag 

!! LAPACK/BLAS version (slower than deigen!)
!! If you uncomment this routine, you must also enable linking against LAPACK
!! in the Makefile.
!subroutine diag2(hmat, eig)
!	implicit none
!	real(ACCU), dimension(:,:), intent(inout) :: hmat
!	real(ACCU), dimension(:), intent(out) :: eig
!	
!	integer, parameter :: SORT_EV=1   ! sort eigenvalues
!	integer, parameter :: EVAL_ONLY=0 ! want eigenvectors as well...
!	integer :: i, dim, info, lwork
!	real(ACCU), dimension(:), allocatable :: work
!
!	dim=getdim()
!	
!	lwork=ilaenv(1,'DSYTRD', 'U', 0, 0, 0, 0)
!	lwork=(lwork+2)*dim
!	allocate(work(lwork))
!
!	print '(a)', '*** entering dsyev'
!	call flush(6)
!	call dsyev('V', 'U', dim, hmat, dim, eig, work, lwork, info)
!	deallocate(work)
!
!end subroutine diag2 

! return a vector containing state degeneracies
subroutine degeneracy(eig, dvec)
	implicit none
	real(ACCU), dimension(:), intent(in) :: eig
	type(dstates), intent(out) :: dvec

	integer :: i, j, dgen, dim
	real(ACCU) :: cmp
	
	dim=getdim()

	dgen=1
	j=1
	cmp=eig(1)
	do i=2,dim
		if (abs(eig(i)-cmp) < DGEN_TOLERANCE) then
			dgen=dgen+1
		else
			dvec%state(j)%eval=cmp
			dvec%state(j)%dgen=dgen
			cmp=eig(i)
			dgen=1
			j=j+1
		endif
	end do
	dvec%state(j)%eval=cmp
	dvec%state(j)%dgen=dgen
	dvec%nstates=j
end subroutine degeneracy

! linear sort, returns a vector with sorted indeces into evec
! evec itself is untouched
! this is not used at the moment, deigen returns sorted eigenvalues...
subroutine sort_eig(evec, sv)
	implicit none
	real(ACCU), dimension(:), intent(in) :: evec
	integer, dimension(:), intent(out) :: sv

	integer :: dim, n, m, tmp

	dim=getdim()
	sv=(/(n, n=1,dim)/)

	do n=1, dim
		do m=n+1, dim
			if (evec(sv(m))-evec(sv(n)) < 0.0) then
				tmp=sv(n)
				sv(n)=sv(m)
				sv(m)=tmp
			end if
		end do
	end do
end subroutine sort_eig
end module eigen

