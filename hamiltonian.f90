!
! $Id: hamiltonian.f90,v 1.12 2001/10/01 16:11:25 jonas Exp $
! 

!
! Routines for setting up the Hamilton operator, and extracting 
! the eigenvalue spectrum. 
!
module hamiltonian
	use param
	use basis
	use three_j
	implicit none
	
	public hamilton
	private 

contains

! The physics...
subroutine hamilton(hmat, pot, mol)
	implicit none
	type(potential), intent(in) :: pot
	type(rotor), intent(in) :: mol
	real(ACCU), dimension(:,:), intent(out) :: hmat

!	integer(4) :: i,j, dimi

	print '(a)', '*** constructing hamiltonian for free rotor'
	call freerot(hmat, mol)

	print '(a)', '*** constructing hamiltonian for hindered rotor'
	call hindrot(hmat, pot)
	
!	open(20, file='h.mtx')
!	dimi=getdim()
!	do i=1,dimi
!		do j=1,dimi
!			write(20,'(f)',advance='no') hmat(i,j)
!		end do
!		write(20,*)
!	end do
!	close(20)

end subroutine hamilton

! set up Hamiltonian elements for hindered rotor
subroutine hindrot(hmat, pot)
	implicit none
	real(ACCU), dimension(:,:), intent(out) :: hmat
	type(potential), intent(in) :: pot

	integer :: jg, je, mg, me, kg, ke
	integer :: m, n, dim
	real (kind=ACCU) :: vja, vka, vma, vj1, vj2, vj3
	real (kind=ACCU) :: h1, h2, h3, qq
	real (kind=ACCU) :: c1, c2, c3
	integer, dimension(:,:), pointer :: basvec
	real(ACCU) :: vj, vk, vm

	vj=pot%vj
	vk=pot%vk
	vm=pot%vm

	c1=sqrt(5.0d0/14.0d0)
	c2=sqrt(1.0d0/90.0d0)
	c3=sqrt(1.0d0/10080.0d0)

	dim=getdim()
	call getbasis(basvec)

	do n=1,dim
		jg=basvec(n,1)
		kg=basvec(n,2)
		mg=basvec(n,3)
		do m=1,n
			je=basvec(m,1)
			ke=basvec(m,2)
			me=basvec(m,3)

			vj1=j3coef(4,je,jg, 0,-me,mg)
			vj2=j3coef(4,je,jg, 4,-me,mg)
			vj3=j3coef(4,je,jg,-4,-me,mg)

			vja=j3coef(4,je,jg,0,-ke,kg)
			vka=j3coef(4,je,jg,2,-ke,kg)
			vma=j3coef(4,je,jg,4,-ke,kg)
			
			qq=(-1.0d0)**(ke+me+1.d0)*(vj1+c1*(vj2+vj3))

			h1=vj*vja*qq
			h2=vk*vka*c2*qq
			h3=vm*vma*c3*qq

			qq=sqrt((2.0d0*jg+1.0d0)*(2.0d0*je+1.0d0))
			qq=hmat(n,m)+qq*(h1+h2+h3)

			hmat(n,m)=qq
			hmat(m,n)=qq
		end do
	end do
end subroutine hindrot

subroutine freerot(hmat, mol)
	implicit none
	real (kind=ACCU), dimension(:,:), intent(out) :: hmat
	type(rotor), intent(in) :: mol

	integer :: jg, je, mg, me, kg, ke
	integer :: m, n, dim
	real(kind=ACCU) :: c1, c2, c3, q
	integer, dimension(:,:), pointer :: basvec
	
	dim=getdim()
	call getbasis(basvec)

	c1=0.25d0*(mol%B-mol%A)
	c2=0.5d0*(mol%A+mol%B)
	c3=mol%C

	do n=1,dim
		jg=basvec(n,1)
		kg=basvec(n,2)
		mg=basvec(n,3)
		do m=1,n
			je=basvec(m,1)
			ke=basvec(m,2)
			me=basvec(m,3)

			if (jg == je .and. mg == me .and. kg == ke-2) then
				q=c1*sqrt(1.0*(jg*(jg+1)-kg*(kg+1))*(jg*(jg+1)-(kg+1)*(kg+2)))
			else if (jg == je .and. mg == me .and. kg == ke+2) then
				q=c1*sqrt(1.0*(jg*(jg+1)-kg*(kg-1))*(jg*(jg+1)-(kg-1)*(kg-2)))
			else if (jg == je .and. mg == me .and. kg == ke) then
				q=c2*(jg*(jg+1)-kg**2)+c3*kg**2
			else
				q=0.0d0
			end if
			hmat(n,m)=q
			hmat(m,n)=q
		end do
	end do
end subroutine freerot

!!!
!!! These routines are for testing purpouses only, not actually used...
!!!
subroutine hindrot_sym(hmat, pot)
	implicit none
	real(ACCU), dimension(:,:), intent(out) :: hmat
	type(potential), intent(in) :: pot

	integer :: jg, je, kg, ke, mg, me
	integer :: m, n, dim
	real (kind=ACCU) :: vja, vka, vj1, vj2, vj3
	real (kind=ACCU) :: h1, h2, qq
	real (kind=ACCU) :: c1, c2
	integer, dimension(:,:), pointer :: basvec
	real(ACCU) :: vj, vk

	vj=pot%vj
	vk=pot%vk

	c1=sqrt(5.0/14.0)
	c2=sqrt(2.0/5040.0)

	dim=getdim()
	call getbasis(basvec)

	do n=1,dim
		jg=basvec(n,1)
		kg=basvec(n,2)
		mg=basvec(n,3)
		do m=1,n 
			je=basvec(m,1)
			ke=basvec(m,2)
			me=basvec(m,3)

			vja=j3coef(4,je,jg,  0,-ke, kg)
			vka=j3coef(4,je,jg,  3,-ke, kg)

			vj1=j3coef(4,je,jg,  0,-me, mg)
			vj2=j3coef(4,je,jg,  4,-me, mg)
			vj3=j3coef(4,je,jg, -4,-me, mg)
			
			qq=(vj1+c1*(vj2+vj3))
			h1=vj*vja*qq
			h2=c2*vk*vka*qq
			
			qq=(-1.0)**(ke+me+1)*sqrt((2.0*jg+1.0)*(2.0*je+1.0))
			qq=hmat(n,m)+qq*(h1+h2)
			hmat(n,m)=qq
			hmat(m,n)=qq
		end do
	end do
end subroutine hindrot_sym

subroutine freerot_sym (hmat, mol)
	implicit none
	real (kind=ACCU), dimension(:,:), intent(out) :: hmat
	type(rotor), intent(in) :: mol

	integer :: m, n, dim
	integer :: j, k
	real(kind=ACCU) :: c1, c2, c3, q
	integer, dimension(:,:), pointer :: basvec
	
	dim=getdim()
	call getbasis(basvec)

	c1=mol%B       !B
	c2=(mol%C-c1)  !C-B

	hmat(:,:)=0.0

	do n=1,dim
		j=basvec(n,1)
		k=basvec(n,2)
		q=c1*j*(j+1)+c2*dble(k**2)
		hmat(n,n)=q
	end do
end subroutine freerot_sym 

end module hamiltonian

