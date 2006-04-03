! Program HINDAROT
! Copyright (C) 2001 Jonas Juselius
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!
! Written by Jonas Juselius <jonas@iki.fi> Sep 2001
! Original implementation by Mika Petterson 
! Diagonalization routines by Dage Sundholm
!
! $Id: hindarot.f90,v 1.26 2001/09/26 15:52:52 jonas Exp $
!
! Be ware! Linear molecules: A=B, C=0 and manifold=0! 
! If A=0, C=0 and B>0 you get the wrong result!!!
!

program hindarot
	use param
	use basis
	use states
	use eigen
	use hamiltonian
	use inout_m
	use transition

	implicit none
	intrinsic abs, mod, sqrt, min

	integer :: jmax, dim, kmf, i, j
	real (kind=ACCU), dimension(:,:), allocatable :: hmat, trmat
	real (kind=ACCU), dimension(:), allocatable :: eig
	type(dstates) :: dvec
	type(potential) :: pot
	type(rotor) :: mol
	
	call init_io(mol, pot)

	do while ( get_next_job(jmax, mol, pot) /= 0 )

		dim=nbasf(jmax, mol)
		call setdim(dim)

		call genbasis(jmax)
		
		kmf=get_manifold()
		if ( kmf /= ALL_MF ) then
			dim=select_k_manifold(kmf)
		end if

		call write_header(mol, pot)

		allocate(hmat(dim, dim))
		allocate(eig(dim)) 

		call hamilton(hmat, pot, mol)
		call diag(hmat, eig)

		call init_dstates(dvec)
		call degeneracy(eig, dvec)
		
		print *
		call write_eig(dvec)
		
		if ( moments_on ) then
			allocate(trmat(dim, dim))
			print *
			print '(a)', '*** calculating transition moments'
            call transmom(hmat, trmat)

			print *
			call write_transitions(trmat, dvec, eig)
			deallocate(trmat)
		end if

		call delbasis
		call del_dstates(dvec)
		deallocate(eig)
		deallocate(hmat)
		call flush(6)
	end do

	print *
	print *, ' *** IT''S A MIRACLE!!! ***'
	call flush(6)
end program hindarot

