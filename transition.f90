!
! $Id: transition.f90,v 1.15 2001/09/28 16:37:53 jonas Exp $
!

module transition
	use param
	use basis
	use three_j
	use states
	use boltzmann
	use glue
	implicit none
	
	public transmom, transition_moments, avg_transitions2, boltz_transitions
    public sum_dtrmom
	private
contains

! Calculate the transition moments.
subroutine transmom(evmat, trmat)
	implicit none
	real(ACCU), dimension(:,:), intent(in)  :: evmat
	real(ACCU), dimension(:,:), intent(out)  :: trmat

	integer ::  i, j, k, dim, nm
	real (kind=ACCU) :: bra, ket, mu, q
	real(ACCU), dimension(:,:), allocatable :: intmat, tmat
	
	dim=getdim()

	call flush(6)
	trmat=0.0d0
	
	! integral storage (x,y,z components)
	allocate(intmat(dim, dim))
	allocate(tmat(dim,dim))

	do k=-1,1 ! x-, y- and z-moments
		call muint(intmat, k)
        tmat=matmul(transpose(evmat), intmat)
        intmat=matmul(tmat, evmat)
        trmat=trmat+intmat**2
	end do
    
	deallocate(tmat)
	deallocate(intmat)
end subroutine 

subroutine transition_moments(evmat, trmat, dvec)
	implicit none
	real(ACCU), dimension(:,:), intent(in)  :: evmat
	real(ACCU), dimension(:,:), intent(out)  :: trmat
	type(dstates), intent(in) :: dvec

	integer ::  i, j, k, dim, nm
	real (kind=ACCU) :: bra, ket, mu, q
	real(ACCU), dimension(:,:), allocatable :: intmat
	real(ACCU), dimension(:), allocatable :: tvec

	nm=0
	do i=1,n_moments(dvec)
		nm=nm+dvec%state(i)%dgen
	end do
	
	dim=getdim()

	call flush(6)
	trmat=0.0d0
	
	! integral storage (x,y,z components)
	allocate(intmat(dim, dim))
	allocate(tvec(dim))

	do k=-1,1 ! x-, y- and z-moments
		call muint(intmat, k)
		do j=1,nm
			tvec=matmul(evmat(:,j), intmat)
			do i=1,j
				q=dot_product(tvec, evmat(:,i))
				trmat(i,j)=trmat(i,j)+q**2
			end do
		end do
	end do
	
	deallocate(tvec)
	deallocate(intmat)
end subroutine transition_moments

! transition moment integrals
subroutine muint(intmat, p)
	implicit none
	real(ACCU), dimension(:,:), intent(out) :: intmat
	integer, intent(in) :: p

	integer :: n, m, dim
	integer :: jg, je, kg, ke, mg, me
	real(ACCU) :: q, coef, bra, ket
	integer, dimension(:,:), pointer :: basvec

	dim=getdim()
	call getbasis(basvec)

	do n=1,dim
		jg=basvec(n,1)
		kg=basvec(n,2)
		mg=basvec(n,3)
		do m=1,dim
			je=basvec(m,1)
			ke=basvec(m,2)
			me=basvec(m,3)
			
			q=dble((2*jg+1)*(2*je+1))
			coef=sqrt(q)

			! <JMK|D10p|J'K'M'>
			bra=j3coef(1,je,jg,0,-ke,kg)
			ket=j3coef(1,je,jg,p,-me,mg)

			q=coef*bra*ket
			intmat(m,n)=q
		end do
	end do
end subroutine muint

subroutine boltz_transitions(trmat, dvec, eig)
	implicit none
	real(ACCU), dimension(:,:), intent(inout) :: trmat
    type(dstates), intent(in) :: dvec
	real(ACCU), dimension(:), intent(in) :: eig

	integer :: i, j, nst, ndim
	real(ACCU) :: Q

	nst=n_moments(dvec)
	Q=partfunc(dvec)
    ndim=getdim()

	do j=1,ndim
		do i=1,ndim
			trmat(i,j)=trmat(i,j)*weight2(eig, j)/Q
		end do
	end do
end subroutine boltz_transitions

subroutine sum_dtrmom(trmat, dvec, A)
	implicit none
	real(ACCU), dimension(:,:), intent(in) :: trmat
	type(dstates), intent(in) :: dvec
	real(ACCU), dimension(:,:), intent(out) :: A
	
	integer ::  i, j, nst
	integer :: pi, pj, di, dj
    real(ACCU) :: pf
	
	nst=dvec%nstates

	pf=partfunc(dvec)

	A=0.0d0

	pj=1
	do j=1,nst
		dj=dvec%state(j)%dgen
		pi=1
		do i=1,j
			di=dvec%state(i)%dgen
			A(i,j)=sum(trmat(pi:(di+pi-1), pj:(dj+pj-1)))
            A(j,i)=A(i,j)
			pi=pi+di
		end do
		pj=pj+dj
	end do
end subroutine 

subroutine avg_transitions2(trmat, A, dvec)
	implicit none
	real(ACCU), dimension(:,:), intent(in) :: trmat
	real(ACCU), dimension(:,:), intent(out) :: A
	type(dstates), intent(in) :: dvec
	
	integer ::  i, j, nst
	integer :: pi, pj, di, dj
	
	nst=n_moments(dvec)

	A(1:nst,1:nst)=0.0d0

	pj=1
	do j=1, nst
		dj=dvec%state(j)%dgen
		pi=1
		do i=1,j
			di=dvec%state(i)%dgen
			A(i,j)=sum(trmat(pi:(di+pi-1),pj:(dj+pj-1)))
			A(j,i)=A(i,j)
			pi=pi+di
		end do
		pj=pj+dj
	end do
end subroutine 

end module transition
