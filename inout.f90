!
! $Id: inout.f90,v 1.37 2001/10/01 16:11:25 jonas Exp $
! 

!
! This module handles all of the input/output needed. It quires the parser for 
! keywords and sets global variables. It also handles printing/writing of
! files. Be warned, the code is a mess...
!

include 'string_m.h'

module inout_m
	use param
	use states
	use boltzmann
	use basis
	use transition
	use externals
	use glue
	use string_m
	implicit none

	include 'getkwf.h'
	
	integer, private, parameter :: eig_fd=10
	integer, private, parameter :: mom_fd=11
	integer, private, parameter :: spec_fd=12

	integer, private :: fileno=0
	integer, private :: jmax, loops=1, start_loops=1
	real(ACCU), private :: cm3
	real(ACCU), private :: vj_g, vk_g, vm_g

	character(80), private :: title
	character(80), private :: basename=''
	
	logical, public :: moments_on=.false.
	logical, public :: spectrum_on=.false.
	logical, public :: print_term=.false.
	logical, private :: number_files=.false.

	type s_array
		real(8), dimension(3) :: x
	end type

	
	public init_io, get_next_job, write_header, write_eig
	public write_transitions, s_array
	private
contains

! init calls the parser, and sets up keywords
subroutine init_io(mol, pot)
	implicit none
	type(rotor), intent(out) :: mol
	type(potential), intent(out) :: pot

	external getarg, iargc
	real(ACCU) :: kelvin=-1.0
	integer :: manifold=ALL_MF
	integer :: status=0, iargc
	integer :: spect, pterm
	character(80) :: arg

	status=iargc()
	if (status /= 1) then
		print *, 'usage: hindarot [infile]'
		call exit(1)
	end if
	call getarg(1, arg)
	call init_parse(trim(arg), status)

	if (status == 1) then
		print *, 'Invalid filename ', trim(arg)
		stop 
	else if (status == 2) then
		stop 'stop.'
	end if

	call set_active_section('HINDAROT')
	
	title=repeat(' ', 80)
	call getkw('title', title)
	call getkw('nprint', nprint)
	call getkw('Jmax', jmax)
	call getkw('file', basename)
	call getkw('manifold', manifold)
	call getkw('terminal', pterm)
	call set_manifold(manifold)
	if (pterm /= 0 ) then
		print_term=.true.
	end if
	
	call getkw('MOLECULE.A', mol%A)
	call getkw('MOLECULE.B', mol%B)
	call getkw('MOLECULE.C', mol%C)

	pot%Vj=0.0
	pot%Vk=0.0
	pot%Vm=0.0

	call getkw('POTENTIAL.Vj', pot%Vj)
	call getkw('POTENTIAL.Vk', pot%Vk)
	call getkw('POTENTIAL.Vm', pot%Vm)
	
	cm3=0.0
	nmom=0
	spect=0
	call getkw('MOMENTS', status)
	call getkw('MOMENTS.T', kelvin)
	call getkw('MOMENTS.nmom', nmom)
	call getkw('MOMENTS.spectrum', spect)
	call getkw('MOMENTS.cm1', cm3)
	if ( status == 1 ) then
		moments_on=.true.
		call set_temperature(kelvin)
		if ( spect /= 0 ) then
			spectrum_on=.true.
		end if
	end if

	loops=1
	call getkw('LOOP', loops)
	if (loops > 1) then
		number_files=.true.
		start_loops=loops
	end if

	if (abs(mol%A) < 1.d-9 .and. abs(mol%C) < 1.d-9 ) then 
		! linear molecule
		call set_k_limit(0)
	else
		call set_k_limit(4)
	end if

	call end_parse()

end subroutine init_io

function get_next_job(j, mol, pot)
	implicit none
	integer, intent(out) :: j
	type(rotor), intent(in) :: mol
	type(potential), intent(inout) :: pot
	logical :: get_next_job

	integer :: jinc,  Kinc, K, status
	real(ACCU) :: vj, vk, vm, Tinc, T
	
	jinc=0
	vj=0.d0; vk=0.d0; vm=0.d0
	Tinc=0.d0
	Kinc=0

	if ( loops /= start_loops ) then
		call getkw('LOOP.Jmax',jinc)
		call getkw('LOOP.Vk',vk)
		call getkw('LOOP.Vj',vj)
		call getkw('LOOP.Vm',vm)
		call getkw('LOOP.T',Tinc)
		call getkw('LOOP.manifold',Kinc)
		vj_g=vj
		vk_g=vk
		vm_g=vm
	end if
	K=get_manifold()
	T=get_temperature()
	
	get_next_job=.false.
	if ( loops > 0) then
		loops=loops-1
		jmax=jmax+jinc
		pot%vj = pot%vj + vj
		pot%vk = pot%vk + vk
		pot%vm = pot%vm + vm
		K=K+Kinc
		call set_manifold(K)
		T=T+Tinc
		call set_temperature(T)
		if (number_files) fileno=fileno+1
		get_next_job=.true.
		vj_g=pot%vj
		vk_g=pot%vk
		vm_g=pot%vm
	end if

	j=jmax
end function get_next_job

! Construct numbered file names
subroutine set_file_name(new, base, ext)
	implicit none
	character(*), intent(out) :: new
	character(*), intent(in) :: base
	character(*), intent(in) :: ext
	
	character(3) :: tmp

	new=repeat(' ',80)
	if ( loops == 0 .and. fileno == 0 ) then
		new=trim(base)//'.'//trim(ext)
		return
	end if
	
	if ( fileno < 10 ) then
		write(tmp, '(i1)') fileno
	else
		write(tmp, '(i2)') fileno
	end if
	
	new=trim(base)//'_'//trim(tmp)//'.'//trim(ext)
end subroutine set_file_name

function open_file(fd, ext, avail, append)
	implicit none
	integer, intent(in) :: fd
	logical, intent(in) :: avail
	character(*), intent(in) :: ext
	integer, optional, intent(in) :: append
	logical :: open_file

	character(80) :: fname

	if (.not. avail) then
		open_file=.false.
		return
	end if

	if ( trim(basename) /= '') then
		call set_file_name(fname, basename, ext)
		if (present(append)) then
			open(fd, file=fname, status='old', position='append')
		else
			open(fd, file=fname, status='replace')
		end if	
		open_file=.true.
	else
		open_file=.false.
	end if
end function

! print an ugly header
subroutine write_header(mol, pot)
	implicit none
	type(rotor), intent(in) :: mol
	type(potential), intent(in) :: pot

	call dump_header(6, mol, pot)

	if ( open_file(eig_fd, 'eig', .true.) ) then
		call dump_header(eig_fd, mol, pot)
		close(eig_fd)
	end if

end subroutine write_header

subroutine dump_header(fd, mol, pot)
	implicit none
	integer, intent(in) :: fd 
	type(rotor), intent(in) :: mol
	type(potential), intent(in) :: pot

	integer :: i, dim
	character :: c

	c='#'
	if ( fd == 6 ) c='*'

	dim=getdim()
	
	write (fd,88)  c,' '
	write (fd,66)  c,'   Title: ', trim(title)
	write (fd,88)  c,' '
	write (fd,88)  c,'   Rotational constants:'
	write (fd,88)  c,'   ----------------------'
	write (fd,99)  c,'   A =', mol%A
	write (fd,99)  c,'   B =', mol%B
	write (fd,99)  c,'   C =', mol%C
	write (fd,88)  c,' '
	write (fd,77)  c,'   Jmax = ', jmax
	write (fd,100) c,'   Vj = ', pot%vj
	write (fd,100) c,'   Vk = ', pot%vk
	write (fd,100) c,'   Vm = ', pot%vm
	write (fd,66)  c,' '
	i=get_manifold()
	if (i /= ALL_MF) then
		write (fd,'(a1,a,i2)') c,'   Selected manifold K =', i
	end if
 	write (fd,'(a1,a,i5)') c,'   The Number of Basis Functions is', dim
	write (fd,*)
	call flush(fd)
66  format(a1,2a)
77  format(a1,a,i3)
88  format(a1,a)
99 	format(a1,a, f11.7, a, f11.7, a, f11.7)
100 format(a1,a, f7.3, a, f7.3, a, f7.3)
end subroutine dump_header

! Print the eigenvalues and their degenerasy
subroutine write_eig(dvec)
	implicit none
	type(dstates), intent(in) :: dvec

	if ( print_term ) call dump_eig(6, dvec)
	
	if ( open_file(eig_fd, 'eig', .true., 1) ) then
		call dump_eig(eig_fd, dvec)
		close(eig_fd)
	end if
end subroutine write_eig

subroutine dump_eig(fd, dvec)
	implicit none
	integer, intent(in) :: fd
	type(dstates), intent(in) :: dvec

	integer :: i, dim, np
	real(ACCU) :: cmp, b
	character :: c

    call getkw('MOLECULE.B', B)
	c='#'
	if ( fd == 6 ) c=' '
	np=n_print(dvec)
	dim=getdim()
	
	write(fd,88) c,' State', 'Eigenvalue', 'g'
	write(fd,'(2a1,a)') c,' ', repeat('-',  26)
	do i=1,dvec%nstates
		if ( i > np ) return
			write(fd, fmt=99) i, dvec%state(i)%eval, dvec%state(i)%dgen
!            write(fd, fmt=99) i, dvec%state(i)%eval/B, dvec%state(i)%dgen
	end do
	call flush(6)
88  format(a1,a6, a15, a5)
99	format(i7,f15.8,i5)
end subroutine dump_eig

subroutine write_transitions(trmat, dvec, eig)
	implicit none
	real(ACCU), dimension(:,:), intent(inout) :: trmat
	type(dstates), intent(in) :: dvec
	real(ACCU), dimension(:), intent(in) :: eig

	if ( print_term ) call dump_trmom(6, trmat, dvec, eig)
	call flush(6)
	
	if ( open_file(mom_fd, 'mom', moments_on) ) then
		call dump_trmom(mom_fd, trmat, dvec, eig)
		close(mom_fd)
	end if

end subroutine 

subroutine dump_trmom(fd, trmat, dvec, eig)
	implicit none
	real(ACCU), dimension(:,:), intent(inout) :: trmat
	type(dstates), intent(in) :: dvec
	real(ACCU), dimension(:), intent(in) :: eig
	integer, intent(in) :: fd
	
	integer :: i, j, k, nst
	real(ACCU) :: kelvin
	real(ACCU), dimension(:,:), allocatable :: A
	character :: c	

	c='#'
	if ( fd == 6 ) c=' '

	kelvin=get_temperature()

    nst=dvec%nstates !n_moments(dvec)
!    nst=getdim()
	allocate(A(nst, nst))

!    call boltz_transitions(trmat,dvec,eig)
!    do i=1,nst
!        do j=1,nst
!            if (trmat(j,i) >1.d-5) print *, trmat(j,i), j, i
!        end do
!    end do
	print *, '*** Summing degenerate transitions'
	print *
    call sum_dtrmom(trmat, dvec, A)

	write (fd, '(a1,a)') c,' Transition moments between states:'
	write (fd, '(a1,a)') c,' -----------------------------------'
	write (fd, '(a,$)') '  '
	write(fd, '(i8,$)') (i,i=1,nst)
	write(fd, *)
	
	do i=1,nst
		write (fd, '(i2,$)') i
!        write (fd, '(a8,$)') ('        ', k=2,i)
		write (fd, '(f8.4,$)') (A(i,j), j=1,i)
		write(fd, '(1a)') ' '
	end do
	
	write(fd, *)
	write(fd,'(a1,a)') c,' Transition moments:'
	write(fd,'(a1,a)') c,' -------------------------'

	write(fd,'(a1,a)') c,' From    To        mu'
	write(fd,'(a1,a)') c,' -------------------------'
	
	do j=1,nst
		do i=1,j
			if (A(i,j) > MU_TOLERANCE ) then
				write(fd,80) i, '  -', j, A(i,j)
			end if
		end do
	end do

	if ( spectrum_on ) then
		if (fd == 6) then
			print *
			if ( print_term ) call dump_spectrum(fd, dvec, A, eig)
		else if ( open_file(spec_fd, 'spc', spectrum_on) ) then
			call dump_spectrum(spec_fd, dvec, A, eig)
			close(fd)
		end if
	end if

	deallocate(A)
80	format(i5,a3,i4,f14.9)
end subroutine dump_trmom

subroutine write_spectrum(dvec, A, eig)
	implicit none
	type(dstates), intent(in) :: dvec
	real(ACCU), dimension(:,:), intent(in) :: A
	real(ACCU), dimension(:), intent(in) :: eig

	character(80) :: fname

	call dump_spectrum(6, dvec, A, eig)
	call flush(6)
	if ( open_file(spec_fd, 'spc', spectrum_on) ) then
		call dump_spectrum(spec_fd, dvec, A, eig)
		close(spec_fd)
	end if
end subroutine write_spectrum

! make a stick spectrum of transitions
subroutine dump_spectrum(fd, dvec, A, eig)
	implicit none
	integer, intent(in) :: fd
	type(dstates), intent(in) :: dvec
	real(ACCU), dimension(:,:), intent(in) :: A
	real(ACCU), dimension(:), intent(in) :: eig

	real(ACCU) :: kelvin
	character :: c	

	print *, '*** Calculating spectrum'
	c='#'
	if ( fd == 6 ) c=' '
	kelvin=get_temperature()

	write(fd, '(a1,a)') c,' Spectrum:'
	write(fd, '(a1,a)') c,    ' -----------------------------------------'
	if ( kelvin > 0.0 ) then
		write(fd, '(a1,a,f5.2,1a)') c,'      cm^-1          I          I('&
		&,kelvin,')'
		write(fd, '(a1,a)') c,' -----------------------------------------'
	else
		write(fd, '(a1,a)') c,'      cm^-1       Intensity'
		write(fd, '(a1,a)') c,' ---------------------------'
	end if
	call branches(fd, dvec, A, eig)
end subroutine dump_spectrum

subroutine branches(fd, dvec, A, eig)
	implicit none
	integer, intent(in) :: fd
	type(dstates), intent(in) :: dvec
	real(ACCU), dimension(:,:), intent(in) :: A
	real(ACCU), dimension(:), intent(in) :: eig

	real(ACCU), dimension(:,:), allocatable :: specvec

	integer :: nst, i, j, k, sdim
	real(ACCU) :: q, kelvin, pf
	nst=n_moments(dvec)
	kelvin=get_temperature()
	
	allocate(specvec(nst*nst,3))
		
	pf=partfunc(dvec)
	k=0
	do j=1,nst
		do i=1,nst
			if (A(i,j) > 1.d-6) then
				k=k+1
				q=dvec%state(j)%eval - dvec%state(i)%eval
				specvec(k,1)=q+cm3
				specvec(k,2)=A(i,j)
!                if (A(i,j) >1.d-4) specvec(k,3)=weight(dvec%state(i))/pf
				specvec(k,3)=A(i,j)*weight(dvec%state(i))/pf
!                specvec(k,3)=weight(dvec%state(i))/pf
			end if
		end do
	end do

	sdim=k

!    print *, '*** Summing transitions'
	call sum_spec(specvec, sdim)
	
	call plot_header(99)
	do i=1,sdim
		if (kelvin > 0.0) then
			q=specvec(i,3)!*specvec(i,3)
			if (q > 1.d-6) then
!                write(fd,99) specvec(i,1:2), q
				write(fd,99) specvec(i,1:3)
			end if
		else if (specvec(i,2) > MU_TOLERANCE) then
			write(fd,88) specvec(i,1:2)
		end if
	end do

	print *, '*** Constructing line spectrum'
	call gauss_spec(dvec, specvec, sdim)

	deallocate(specvec)

88 	format(2f14.7)
99 	format(3f14.7)
end subroutine branches

subroutine gauss_spec(dvec, spec, dim)
	type(dstates), intent(in) :: dvec
	real(8), dimension(:,:), intent(in) :: spec
	integer(4) :: dim
	
	real(8) :: head, step, cm1, xp, y, q
	integer(4) :: i, j, npts

!    head=-dvec%state(dvec%nstates)%eval
	head=spec(1,1)
	step=5.d-2
	call getkw('MOMENTS.linestep',step)
	npts=-head/step*2.d0
	xp=5.d0
	call getkw('MOMENTS.gauss_exp',xp)

	if (.not.open_file(99, 'lns', .true.)) return
	call plot_header(99)
	do i=0,npts
		cm1=head+i*step
		y=0.d0
		do j=1,dim
			q=cm1-spec(j,1)
			if (abs(q) > 2.5d0) cycle
			y=y+spec(j,3)*exp(-xp*q**2)
		end do
		if (y > 1.d-3) write(99, *) cm1, y
	end do
	close(99)
end subroutine

! Sum all transitions with the same wave numbers.
! The spectra must be sorted for this to work.
subroutine sum_spec(svec, dim)
	implicit none
	real(ACCU), dimension(:,:), intent(inout) :: svec
	integer, intent(inout) :: dim

	interface
		function kusse(a,b) result(p)
			type s_array
				real(8), dimension(3) :: x
			end type
			type(s_array) :: a, b
			integer(4) :: p
		end function
	end interface

	integer :: n, m, k
	real(ACCU), dimension(:,:), allocatable :: foo

	type(s_array), dimension(:), allocatable :: sarray

	allocate(sarray(dim))
	do n=1,dim
		sarray(n)%x=svec(n,:)
	end do

	print *, '*** Sorting spectrum', dim
	call qsort(sarray,dim,3*8,kusse)
	do n=1,dim
		svec(n,:)=sarray(n)%x
	end do

	deallocate(sarray)
	
!    allocate(foo(dim,3))
!    foo=0.0
!    k=1
!    m=1
!    do n=1,dim
!        if (abs(svec(m,1)-svec(n,1)) < 1.d-4) then
!            foo(k,1)=svec(m,1)
!            foo(k,2)=foo(k,2)+svec(n,2)
!            foo(k,3)=foo(k,3)+svec(n,3)
!        else
!            k=k+1
!            foo(k,:)=svec(n,:)
!            m=n
!        end if
!    end do
!    dim=k
!    svec=0.0
!    svec(1:dim,:)=foo(1:dim,:)
!    deallocate(foo)
end subroutine sum_spec

subroutine sort_spec(svec, dim)
	implicit none
	real(ACCU), dimension(:,:), intent(inout) :: svec
	integer, intent(in) :: dim
	
	integer :: n, m
	real(ACCU) :: tmp(3)
	
	do n=1, dim
		do m=n+1, dim
			if (svec(m,1)-svec(n,1) < 0.0) then
				tmp=svec(n,:)
				svec(n,:)=svec(m,:)
				svec(m,:)=tmp
			end if
		end do
	end do
end subroutine sort_spec

subroutine plot_header(fd)
	integer(4), intent(in) :: fd

	integer(4) :: K, J
	real(8) :: T, gexp
	K=get_manifold()
	T=get_temperature()

	call getkw('Jmax', J)
	call getkw('MOMENTS.gauss_exp', gexp)

	write(fd,100) J, K, T, vj_g, vk_g, vm_g, gexp

100 format('# J=',i2,' K=',i1, ' T=', f6.2,' Vj=',f6.2,' Vk=',f6.2,' Vm=', &
	f6.2,' exp=',f6.2)
end subroutine

end module inout_m

function kusse(a,b) result(p)
	use inout_m
	type(s_array) :: a, b
	integer(4) :: p

	if (a%x(1) < b%x(1)) then
		p=-1
	else if (a%x(1) > b%x(1)) then
		p=1
	else
		p=0
	end if
end function

