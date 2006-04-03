!
! $Id: externals.f90,v 1.5 2001/10/01 13:23:18 jonas Exp $
!

!
! Interface blocks for external routines
!
module externals
	implicit none

	interface

!
! The arrays in deigen and packhm must be defined to scalar size due
! to an Alpha compiler bug... sucks.
!

		subroutine deigen(a, r, n, mv, mfkr)
			integer   :: n, mv, mfkr
			real(8), dimension(n) :: a, r
		end subroutine deigen

		subroutine packhm(h, hpack, dim)
			integer   :: dim
			real(8), dimension(dim,dim) :: h
			real(8), dimension(dim) :: hpack
		end subroutine packhm
		
		subroutine dsyev(c, uplo, dim, A, lda, eig, work, lwork, info)
			integer   :: dim, lda, lwork, info
			character :: c, uplo
			real(8), dimension(:) :: eig, work
			real(8), dimension(:,:) :: A
		end subroutine dsyev

		function ilaenv(a,b,c, n1, n2, n3, n4)
			integer :: a, n1, n2, n3, n4
			character :: b
			character(*) :: c
			integer ilaenv
		end function

	end interface
end module externals 

