!
! $Id: basis.f90,v 1.10 2001/09/30 16:30:40 jonas Exp $
!

module basis
	use param
	implicit none

	integer, dimension(:,:), allocatable, target, private :: cbasis, kbasis
	integer, private :: active=1 ! 1=cbasis, 2=kbasis
	integer, private :: manifold=ALL_MF

	public nbasf, genbasis, delbasis, getbasis
	public restore_k_manifold, select_k_manifold
	public set_manifold, get_manifold

	private
contains

! Calculate the number of basis functions
function nbasf(jmax, mol) 
	implicit none
	integer, intent(in) :: jmax
	type(rotor), intent(in) :: mol
	integer :: nbasf

	integer :: j, k_limit
	logical :: linear

	linear=is_linear(mol)
	
	k_limit=get_k_limit()
	nbasf=0
	do j=0, jmax
		if ( linear ) then
			nbasf=nbasf+(2*j+1)
		else if (j > k_limit) then
			nbasf=nbasf+(2*j+1)*9
		else
			nbasf=nbasf+(2*j+1)**2
		end if
	end do
end function nbasf

! Create a 3-vector of spherical harmonic indeces.
! Use of basvec makes the rest of the code cleaner, and
! will at some point also enable a very convenient way of block diagonalizing
! the Hamilton matrix.
subroutine genbasis(jmax)
	implicit none
	integer, intent(in) :: jmax

	integer :: j, k, m 
	integer :: n, dim, k_limit

	dim=getdim()
	k_limit=get_k_limit()
	allocate(cbasis(dim,3))

	n=1
	do j = 0, jmax
		do k = -min(k_limit,j),min(k_limit,j)
			do m = -j,j
				cbasis(n,1:3)=(/j,k,m/)
				n=n+1
			end do 
		end do
	end do
	active=1
end subroutine genbasis

! k1 is lower limit, k2 upper limit 
function select_k_manifold(k1, k2)
	implicit none
	integer, intent(in) :: k1
	integer, optional, intent(in) :: k2
	integer :: select_k_manifold
	
	integer :: dim, i, j, ak1, ak2, bk
	
	ak1=abs(k1)
	if ( present(k2) ) then
		ak2=abs(k2)
		if ( ak2 < ak1 ) then
			stop 'select_k_manifold(): K2 > K1!'
		else if ( ak1+ak2 > get_k_limit() ) then
			stop 'select_k_manifold(): K1 + K2 > K_LIMIT!'
		end if
	else if ( ak1 > get_k_limit() ) then
		stop 'select_k_manifold(): K > K_LIMIT!'
	else
		ak2=ak1
	end if
	
	dim=getdim()
	if (.not.allocated(kbasis)) then
		allocate(kbasis(dim,3))
	end if
	
	j=0
	do i=1,dim
		bk=abs(cbasis(i,2))
		if (bk >= ak1 .and. bk <= ak2 ) then
			j=j+1
			kbasis(j,1:3)=cbasis(i,1:3)
		end if
	end do
	call setdim(j)
	active=2        ! set active basis to kbasis
	call set_manifold(ak1)
	select_k_manifold=j
end function select_k_manifold

function restore_k_manifold()
	implicit none
	integer :: restore_k_manifold

	integer :: dim(2)
	
	dim=shape(cbasis)
	call setdim(dim(1))
	active=1        ! set active basis to cbasis
	call set_manifold(ALL_MF)
	restore_k_manifold=dim(1)
end function restore_k_manifold

subroutine getbasis(bvec)
	implicit none
	integer, dimension(:,:), pointer :: bvec
	
	select case(active)
		case(1)
			bvec=>cbasis
		case(2)
			bvec=>kbasis
	end select
end subroutine getbasis

subroutine set_manifold(k)
	implicit none
	integer, intent(in) :: k
	
	if ( abs(k) > get_k_limit() .and. k /= ALL_MF ) then
		stop 'set_manifold(): Invalid manifold, |K| > 4'
	end if
	manifold=k
end subroutine set_manifold

function get_manifold()
	implicit none
	integer :: get_manifold

	get_manifold=manifold
end function get_manifold

! basvec destructor 
subroutine delbasis
	implicit none
	if (allocated(cbasis)) deallocate(cbasis)
	if (allocated(kbasis)) deallocate(kbasis)
end subroutine delbasis
end module basis
