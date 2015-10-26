

	module mod_output
	
	INTEGER :: AtomFile=10, EnergyFile=11, TempFile=17
	
	contains
	
	subroutine initOutput(atomFileName, energyFileName, tempFileName)
	character *80 atomFileName, energyFileName, tempFileName	
	
		open(unit=AtomFile, file=atomFileName, status='unknown')
		open(unit=EnergyFile, file=energyFileName, status='unknown')
		open(unit=TempFile, file=tempFileName, status='unknown')

		
	end subroutine initOutput

	
	
	subroutine writeData(steps, atomCoord,
     $		velocity, atomDispl, atomMass, isDamped)
	use mod_global
	implicit none
	double precision atomCoord(ndf, *), velocity(ndf, *), atomMass,
     $		atomDispl(ndf,*)
	double precision cellSize, systemEnergy, sysEnergy75, temp75,
     $		temp55, x0, y0, x1, y1, x, y, e0,e1,e2,size,tempwhole,
     $		sysEnergy55, tempIntermediate,temp45,temp35,temp25,
     $		temp15

	integer i0, i1, i2, i3, i, steps, isDamped(*), numAtoms200,
     $		numAtoms75, numAtoms55, numAtoms,numAtoms45,
     $		numAtoms35,numAtoms25,numAtoms15
     
! 	size = 55
! 	rc = findCornerAtoms(atomCoord, atomID,
!      $		i0, i1, i2, i3, size)
  	
C---	Print Data	
	cellSize = 200.0
	systemEnergy= getSystemEnergy(atomCoord, velocity, 
     $	 atomMass, cellSize)

	cellSize = 75.0
	sysEnergy75 = getSystemEnergy(atomCoord, velocity, 
     $	 atomMass, cellSize)

	cellSize = 55.0
	sysEnergy55 = getSystemEnergy(atomCoord, velocity, 
     $	 atomMass, cellSize)

    	size = 200.0
     	tempWhole = getTempRegion(atomCoord,
     $		isDamped, velocity, atomMass, size, numAtoms200) 

	
	size = 55.0
    	temp55 = getTempRegion(atomCoord,
     $		isDamped, velocity, atomMass, size, numAtoms55) 
        size = 45.0
        temp45 = getTempRegion(atomCoord,
     $          isDamped, velocity, atomMass, size, numAtoms45)
        size = 35.0
        temp35 = getTempRegion(atomCoord,
     $          isDamped, velocity, atomMass, size, numAtoms35)
        size = 25.0
        temp25 = getTempRegion(atomCoord,
     $          isDamped, velocity, atomMass, size, numAtoms25)
        size = 15.0
        temp15 = getTempRegion(atomCoord,
     $          isDamped, velocity, atomMass, size, numAtoms15)
          	
	numAtoms = numAtoms200 - numAtoms55
	tempIntermediate = (tempWhole * numAtoms200 - 
     $		temp55 * numAtoms55)/dfloat(numAtoms)
C	write(*,*)numAtoms55, numAtoms45,numAtoms35,
C     $            numAtoms25,numAtoms15,numAtoms200
C	stop

	write(EnergyFile,9)steps, systemEnergy, sysEnergy75, sysEnergy55
9	format('steps: ', 1x, i5, 2x, 'energy:', 3(1x, 1pe16.9))

	write (TempFile,11)steps,tempWhole,temp55,temp45,temp35,temp15
11	format('steps: ', 1x, i5, 1x, 'temp: ', 1x, 5(f12.4, 1x))

	
C	if (steps .ge. 5000 .and. mod(steps, 1000) .eq. 0) then
C		size = 75.0
C		call writeAtomConfig(atomCoord, atomDispl, size)
C	endif

								
	return
	end subroutine writeData
		

	
	subroutine writeAtomConfig(atomCoord, atomDispl, size) 
	use mod_global
	use mod_dynamo
	implicit none	
	double precision atomCoord(ndf, *), atomDispl(ndf, *)
	double precision, pointer::initCoord(:,:)
	double precision dx, dy, dz, x, y, z, size
	integer iAtom, atomType
	logical firstTime
	data firstTime /.true./ 

	if (firstTime .eq. .true.) then
		allocate(initCoord(ndf, numnp))

		do 5 iAtom = 1, numnp
			initCoord(1, iAtom) = atomCoord(1, iAtom)
			initCoord(2, iAtom) = atomCoord(2, iAtom)
			initCoord(3, iAtom) = atomCoord(3, iAtom)
5		continue

		firstTime = .false.
	endif

		
	do 10 iAtom = 1, numnp
		atomType = isRelaxed(iAtom)
		
		x = atomCoord(1, iAtom)
		y = atomCoord(2, iAtom)
		z = atomCoord(3, iAtom)

		dx = atomDispl(1, iAtom)
		dy = atomDispl(2, iAtom)
		dz = atomDispl(3, iAtom)


		if (atomType .eq. indexContinuum) goto 10
		if (dabs(x) .gt. size) goto 10
		if (dabs(y) .gt. size) goto 10
		write (AtomFile,9)  x, y, z, dx, dy, dz, isRelaxed(iAtom)
10	continue
	write (AtomFile,*)

9		format(3(f16.9,2x), 3(1pe16.9,2x), i3)	

	return
	end subroutine writeAtomConfig


	
	
C	Finds the energy of an annular rectangular box with dimensions
C	xmin, ymin, xmax, ymax.
	double precision function getEnergyRegion(atomCoord, velocity, 
     $	 atomMass, xmin, xmax, ymin, ymax)
	use mod_global
	implicit none
	double precision atomCoord(ndf, *), velocity(ndf,*),
     $	atomMass, ke, pe, v, size, x, y, xmin, xmax, 
     $	ymin, ymax
     
	integer iAtom, j, numAtoms

	ke = 0.0		!kinetic energy
	pe = 0.0		!potential energy
	boltzmannConst = 8.62906e-5	! eV/K
	numAtoms = 0

	do 10 iAtom = 1, numnp

		if (isRelaxed(iAtom) .lt. 1) then
			   goto 10
		endif	

		x = atomCoord(1, iAtom)
		y = atomCoord(2, iAtom)

		if ( (dabs(x) .lt. xmin) .and. (dabs(y) .lt. ymin)) goto 10 
		if ( (dabs(x) .gt. xmax) .or. (dabs(y) .gt. ymax)) goto 10
		
		do 5 j = 1, ndf
			v = velocity(j, iAtom)		
			ke = ke + v*v 
5		continue

C		write(6,9) x, y, xmin, xmax, ymin, ymax
9		format(6(2x, f11.4))
		pe = pe + energy(iAtom)
		numAtoms = numAtoms + 1
10	continue

	ke = 0.5 * atomMass * ke		! gives kinetic energy in eV 
	getEnergyRegion = pe + ke	! gives total energy
C	print*, 'numatoms: ', numAtoms, ' size: ', xmin, xmax, ymin, ymax
	return
	end function getEnergyRegion

	
	
	
	double precision function getTempRegion(atomCoord,
     $		isDamped, velocity, atomMass, size, numAtoms) 

	use mod_global
	implicit none
	double precision atomCoord(ndf, *),  atomMass,
     $			velocity(ndf, *)
	double precision  kineticEnergy, v
	double precision  temperature, size, x, y
	integer iAtom, numAtoms, j, isDamped(*)


	kineticEnergy = 0.0
	numAtoms = 0
	do 10 iAtom = 1, numnp

		if (isRelaxed(iAtom) .lt. 1) goto 10

		x = atomCoord(1, iAtom)
		Y = atomCoord(2, iAtom)

		if (abs(x) .gt. size .or. abs(y) .gt. size) goto 10

		numAtoms = numAtoms + 1
		do 20 j = 1, ndf
			v = velocity(j, iAtom)
			kineticEnergy = kineticEnergy + v *v
20		continue

10	continue
	

	kineticEnergy = 0.5 * kineticEnergy * atomMass
	temperature = kineticEnergy/(1.5 * boltzmannConst * numAtoms)

	getTempRegion = temperature
	return
	end function getTempRegion


	
	double precision function getSystemEnergy(atomCoord, velocity, 
     $	 atomMass, size)
	use mod_global
	implicit none
	double precision atomCoord(ndf, *), velocity(ndf,*),
     $	atomMass, ke, pe, v, size, x, y
	integer iAtom, j, numAtoms

	ke = 0.0		!kinetic energy
	pe = 0.0		!potential energy
	numAtoms = 0

	do 10 iAtom = 1, numnp

		if (isRelaxed(iAtom) .lt. 1) then
			   goto 10
		endif	

		x = atomCoord(1, iAtom)
		y = atomCoord(2, iAtom)

		if (abs(x) .gt. size .or. abs(y) .gt. size) goto 10

		do 5 j = 1, ndf
			v = velocity(j, iAtom)		
			ke = ke + v*v 
5		continue

	     pe = pe + energy(iAtom)
		numAtoms = numAtoms + 1
10	continue

	ke = 0.5 * atomMass * ke		! gives kinetic energy in eV 
	getSystemEnergy = pe + ke	! gives potential energy
	return
	end function getSystemEnergy
	
	
	
	integer function findCornerAtoms(atomCoord, atomID,
     $		i0, i1, i2, i3, size)

	use mod_global
	implicit none
	double precision atomCoord(ndf, *)
	double precision xmin, xmax, ymin, ymax, x, y, size
	integer iAtom, j, flag, atomID(ndf, *), nAtoms
	integer i0, i1, i2, i3
	double precision d0, d1, d2, d3, d
	double precision distance


	d0 = 10.0
	d1 = 10.0
	d2 = 10.0
	d3 = 10.0

	nAtoms = 0
	do 10 iAtom = 1, numnp
		
		if (isRelaxed(iAtom) .lt. 1) goto 10

		nAtoms = nAtoms + 1

		x = atomCoord(1, iAtom)
		y = atomCoord(2, iAtom)

		d= distance(x, y, -size, -size)
		if (d .lt. d0) then
		   d0 = d
		   i0 = iAtom
		endif

		d= distance(x, y, size, -size)
		if (d .lt. d1) then
		   d1 = d
		   i1 = iAtom
		endif

		d= distance(x, y, -size, size)
		if (d .lt. d2) then
		   d2 = d
		   i2 = iAtom
		endif

		d= distance(x, y, size, size)
		if (d .lt. d3) then
		   d3 = d
		   i3 = iAtom
		endif
10	continue

	findCornerAtoms = 0 
! 	print*, 'Corner Atoms:  ', i0, i1, i2, i3
! 	print*, 'distances: ', d0, d1, d2, d3

	return
	end function findCornerAtoms
	
	
	subroutine writeNeighborList(neighborFileName)
	use mod_global
	use mod_dynamo
	implicit none
	integer iAtom, lowIndex, highIndex, range, j
	character *80 neighborFileName
	
	open(unit=14, file=neighborFileName, status='unknown')
!	Write the list of indices	
	write(14, *)nnindx(0)
	do 10 iAtom = 1, numnp
C		if (isRelaxed(iAtom) .eq. indexContinuum) goto 10
		write(14, *) iAtom, nnindx(iAtom)
10	continue

!	Write the list of neighbors
	do 20 iAtom = 1, numnp
		if (isRelaxed(iAtom) .eq. indexContinuum) goto 20
		
cc		lowIndex = nnindx(iAtom - 1)
		lowIndex = 0
		highindex = nnindx(iAtom)
		range = highIndex - lowIndex
		
		do 30 j = 1, range
			write(14, *) iAtom, nnlst(iAtom,j+lowIndex)
30		continue	
		write(14,*)
20	continue	
	close(14)
	
	
	return
	end subroutine writeNeighborList
	

C	Read the neighborlist
	subroutine readNeighborList(neighborFileName)
	use mod_global
	use mod_dynamo
	implicit none
	integer iAtom, lowIndex, highIndex, range, j, k
	character *80 neighborFileName
	
	open(unit=14, file=neighborFileName, status='old')
	
!	Read the list of indices	
	read(14, *)nnindx(0)
	do 10 iAtom = 1, numnp
C		if (isRelaxed(iAtom) .eq. indexContinuum) goto 10
		read(14, *) j, nnindx(iAtom)
		if (j .ne. iAtom) then
			print*, 'Stopping.. Error in reading index list'
			stop
		endif
10	continue



!	Read the list of neighbors
	do 20 iAtom = 1, numnp
		if (isRelaxed(iAtom) .eq. indexContinuum) goto 20
		
cc		lowIndex = nnindx(iAtom - 1)
		lowIndex = 0
		highindex = nnindx(iAtom)
		range = highIndex - lowIndex
		

		do 30 j = 1, range
			read(14, *) k, nnlst(k,j+lowIndex)
			
		if (k .ne. iAtom) then
			print*, 'Stopping.. Error in reading list of neighbors'
			print*, k, j, iAtom, range
			stop
		endif
		
! 		print*, iAtom, k, nnlst(iAtom, j+lowIndex)
30		continue
! 		print*
		read(14,*)
20	continue	
	close(14)	
	
	return
	end subroutine readNeighborList
	
		
	subroutine closeOutput()
		close(AtomFile)
		close(EnergyFile)
		close(TempFile)
	end subroutine closeOutput
	
	
!	subroutine writeRestartFile(atomCoord, velocity, acceleration)
	
	end module mod_output
	
