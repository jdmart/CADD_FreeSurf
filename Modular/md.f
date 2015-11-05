!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Find forces on atoms
	subroutine getEnergiesAndForces (id, atomCoord,ix, f, atomDispl,
     &		aveDispl, atomForce, atomMass, systemEnergy, 
     &		MoveAtoms,MoveDisl, FullField, solveFEM,strainE0,iFem)

	use mod_global
        use mod_timming
	implicit none


	double precision atomDispl(ndf,*), atomCoord(ndf,*),f(ndf,*), 
     &		atomForce(ndf,*), systemEnergy, atomMass, aveDispl(ndf,*)
	integer  id(ndf,*), ix(nen1,*)
	logical MoveAtoms,MoveDisl, FullField, solveFEM, ChangeTime
	integer iAtom, j, i, iFem, jdof
	double precision displ_old(3), displ_new(3)
        double precision strainE0


        call cpu_time(ct2)
	if (solveFEM .eq. .true.) then
!!	   Solve FEM
           call vafuncMD(id, atomCoord, ix, f, aveDispl, 
     &	    atomForce, systemEnergy, MoveAtoms, MoveDisl, 
     &      FullField, strainE0,iFem, Moved)

!!		Get Forces and displacements, specifically on PAD atoms
!JM   loop adjustment made to include indenter atoms
!JM   and BC's applied to them. numnp vs numnpp1
!	   do iatom = 1, numnp
	   do iatom = 1, numnpp1
		atomForce(1:ndf, iatom) =
     &		 -atomForce(1:ndf, iatom)
   	    if (isRelaxed(iatom) .eq. indexContinuum
     &		 .or. isRelaxed(iatom) .eq. indexPad) then
 		    atomDispl(1:ndf, iatom) =
     &		     aveDispl(1:ndf, iAtom)

	    endif

!!JM Hardcode fix to add displacements to BC atoms
   	     if (isRelaxed(iatom) .eq. indexAtom
     &		 .or. isRelaxed(iatom) .eq. indexInterface) then
           do jdof = 1,ndf
              if (id(jdof,iatom) .eq. 1) then

                atomDispl(jdof, iatom) = aveDispl(jdof, iAtom)    	      

              endif

	       enddo

	     endif

	    enddo

!!	   Zero out forces for fixed nodes
	   do i=1,ndf
	      do j=1,numnp
                 if (idtemp(i,j)) then
                    atomForce(i,j)=0
                 endif
            enddo
	    enddo	    
	endif
        call cpu_time(ct3)
        ct7=ct7+ct3-ct2
	
!!	Zero out forces on all MD atoms in the atomistic region	
!	call InitialiseEnergy(.true.,.true.,id,atomForce,f)
	do iAtom = 1, numnp
	  if (isRelaxed(iAtom) .eq. indexAtom .or.
     &	   isRelaxed(iAtom) .eq. indexInterface)	then
     		atomForce(1:3, iAtom) = 0.d0
	  endif	
        enddo
		
!!	Find forces from EAM
        call cpu_time(ct2)
	call ProcessClump(id,atomCoord,ix,f,atomDispl,
     &	  atomForce,.true.,.true.)
        call cpu_time(ct3)
        ct4=ct4+ct3-ct2

!!JM Hard Code to 0 out forces on BC atoms
       do i=1,ndf
!JM   loop adjustment made to include indenter atoms
!JM   and BC's applied to them. numnp vs numnpp1
!          do j=1,numnp
          do j=1,numnpp1
             if (idtemp(i,j)) then
	            if (isRelaxed(iAtom) .eq. indexAtom .or.
     &	          isRelaxed(iAtom) .eq. indexInterface)	then
                     atomForce(i,j)=0
                endif
             endif
          enddo
       enddo

	return
	end






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--	Move atoms according to a velocity verlet scheme
CC--Jun Song Comment: Add new variables and embed
CC--getEnergiesAndForces subroutine 
	subroutine velocityVerlet(atomCoord, atomDispl, atomForce,
     &		atomID, isDamped, oldAcceleration, acceleration, 
     &		oldVelocity, velocity, timestep, atomMass, 
     &		langevinCoeff, requiredTemp, currentTemp,
     &		simulationCell, rampLevel, Thermostat, ix, f, 
     &		aveDispl, systemEnergy, MoveAtoms, MoveDisl,
     & 		Fullfield, solveFEM, strainE0,iFem)	 

	use mod_global
	use mod_grain
        use mod_poten
	use mod_common
	implicit none
	
!	type(graintype), dimension(:), pointer :: grains
	
	integer n,iprint,maxfn,NumDis, atomID(ndf, *)
	double precision atomDispl(ndf,*), atomForce(ndf,*), 
     &		atomCoord(ndf,*), rampLevel(*),
     &		oldAcceleration(ndf,*), acceleration(ndf,*),
     &		oldVelocity(ndf,*), velocity(ndf,*), timestep, atomMass,
     &		maxDispl, displ, langevinCoeff, dampCoeff, interiorCoeff,
     &		currentTemp, requiredTemp, xmin, xmax, ymin,
     &		ymax, damped_width, zetaDot, zeta
     
CC--Jun Song: new local variables
	double precision f(ndf,*), systemEnergy, aveDispl(ndf,*),
     & 		strainE0, newdispl(ndf), currentCoord(ndf), timestepHH
	integer ix(nen1,*)
	logical MoveAtoms, MoveDisl, Fullfield, solveFEM
CC--Mod Ends

	logical ConvOnFn,AddedSlip,CheckSlip,LostSlip,dislcheck
	logical Printing, debug, plot
	Character(len=16) :: Damping_Mode
	TYPE(region) simulationCell
	TYPE(MD_Thermostat) Thermostat
	common/debugger/debug
	data zeta /0.d0/
CC--initialize zeta, only called once!!
	
C--	Local variables	
	integer iAtom, j, nAtoms, isDamped(*), location
CC--constrain motion to be 2D
	integer ndf2D
        integer selectThermostat
	double precision  Tatom, TInterior, rampValue,
     &		 meanSqDispl, position, size, rmsForce, TDampRegion, 
     &		ramp, stadiumTemp, neighborTemp, Tdamp,randForce,
     &          ZBQLU01
C--	Functions
	double precision getKineticTemp, getNearNeighborTemp,
     &		 getRamp, BerendsenDampCoeff, getTemperature,
     &		NoseHooverCoeff

CC--Jun Song: half the timestep
CC--dampfactor for thermstat, trsfactor for T rescale
CC--store damp coefficients for each atom to use in final int
CC--hack: Htrsfactor used to stablize H initially (maintain T_H<10000)
CC--Global NHrescale controls # of steps do Htrsfactor
cc--Store the intemediate V for final integration for H atom
CC--Store the random force in initial integration for final integration
	double precision dthalf,dampfactor(i_final),NTdamp(i_final)
	double precision NdampCoeff(i_final), trsfactor,deltaTemp
	double precision :: Htrsfactor=1.0d0
	double precision TempOldVelocity(ndf,i_final)
	double precision RandFc(ndf,i_final)
CC--Mod Ends
	
	logical useNoseHoover, useLangevin, useMarder, useNve
	integer iFem
	
	dampfactor=0d0
	deltaTemp=0d0
	nAtoms = 0
	maxDispl = 0.5d0
CCJSTest: Constrain motion to be 2D
        ndf2D=ndf-1

!--	timestep for H atom
	timestepHH=timestep1/indextimeH


	call getBoxTemp(atomCoord, Velocity, atomMass, simulationCell,
     &		TInterior, stadiumTemp, currentTemp) 
 	
        selectThermostat=0
     	!! Check for NoseHoover
	useNoseHoover = .false.
     	location = Index(Thermostat%Type, 'NoseHoover')
	if (location .gt. 0) useNoseHoover = .true.
       if (location .gt. 0) selectThermostat=1

       location=0 
   	!! Check for Langevin
	useLangevin = .false.
     	location = Index(Thermostat%Type, 'Langevin')
	if (location .gt. 0) useLangevin = .true.
       if (location .gt. 0) selectThermostat=2

        location=0 
   	!! Check for Marder
	useMarder = .false.
     	location = Index(Thermostat%Type, 'Marder')
	if (location .gt. 0) useMarder = .true.
       if (location .gt. 0) selectThermostat=3
       
CC--Jun Song: Adding new thermostat-T rescaling
CC--UseRescale is global parameter, set in dosteps
       if (UseRescale) selectThermostat=4
       
        location=0 
   	!! Check for NVESTAT
	useNve = .false.
     	location = Index(Thermostat%Type, 'Nvestat')
	if (location .gt. 0) useNve = .true.
        if (location .gt. 0) selectThermostat=5

CC--Mod Ends

c       print*, 'Select Thermostat', selectThermostat
		     
	Damping_Mode = Thermostat%Damping_Mode

ccccccccccccccccccccccccccccccJun Songcccccccccccccccccccccccccccccccc
CC--Replace useNoseHoover with selectThermostat.eq.1
CC--Thus do not update coefficient if doing T rescaling
CC--Use timestepHH because zeta is updated when both H and Ni move
CC--and do MD first for H for SimStep=1 to indextimeH-1 	     
	if (selectThermostat.eq.1) then  
     		zetaDot = NoseHooverCoeff(requiredTemp, currentTemp)
		zeta = zeta + zetaDot * timestepHH
c		print*, 'Nose Hoover'
	endif

CC--Jun Song: opening files to output T and disp for H atoms	
	if(SimStep .eq. 1) then
       	  open(unit=9191,file='TempHydrogen.dat',status='unknown')
       	  open(unit=9192,file='DisplHydrogen.dat',status='unknown')
       	endif
CC--Comment End

	deltaTemp=ABS(stadiumTemp-requiredTemp)
!-	JS: rescale T if exceeds Twindow 
	trsfactor=1.0d0
	if(deltaTemp .gt. Twindow) then
	   trsfactor=dsqrt(requiredTemp/currentTemp)
	endif
	

CC--Jun Song: Initial integration, update V1 and Disp
      do 101 iAtom = i_initial, i_final
	  if(atomSpecie(iAtom).eq.2)then
	    timestep=timestepHH
	  else
            timestep=timestep1
	  endif
	  
!-	used to store dampfactor for second integral	  
	 dampfactor(iAtom)=0.0d0
	 NdampCoeff(iAtom)=0.0d0
	 
	 dthalf=0.5*timestep
		
!-	Skip near continuum nodes and pad atoms
         atomMass=amass(atomSpecie(iAtom))*1.0d-24 
	  if (isRelaxed(iAtom) .eq. indexContinuum) goto 101
	  if (isRelaxed(iAtom) .eq. indexPad) goto 101

	  neighborTemp = getNearNeighborTemp(iAtom, atomCoord,
     &	  atomDispl,  isDamped, oldVelocity, atomMass)   

!--	output T and Disp for all H atoms
	if(mod(SimStep,10).eq.1) then
         if (atomSpecie(iAtom).eq.2) then 
          Tatom = getKineticTemp(atomCoord, 
     $	    oldVelocity, iAtom, atomMass)
          write(9191,*) atomMass/1.0365*1.0d+28,Tatom,iAtom
          write(9192,'(3f16.11,i8)') atomDispl(1,iAtom), 
     $      atomDispl(2,iAtom),atomDispl(3,iAtom),iAtom

         endif
        endif

!-	  Find the damping coeff
	dampCoeff = 0.d0
 	if (isDamped(iAtom) .eq. 0) then	! Interior atom
   	   dampCoeff = 0.d0
      	else	       
	  Tatom = getKineticTemp(atomCoord, 
     $	    oldVelocity, iAtom, atomMass)
     
!-	Initialize Tdamp to be stadiumTemp at beginning
	Tdamp = stadiumTemp

!- 	Update Tdamp and coefficient only when not rescaling	
	if(selectThermostat.ne.4) then
 	    if (useNoseHoover) then
		dampCoeff = zeta
	    else
		rampValue = rampLevel(iAtom)
		location = index(Damping_mode,'Neighbor')
		if ( location .gt. 0 ) then
		   Tdamp = neighborTemp
		endif
					
		location=index(Damping_mode,'Stadium')	
		if (location .gt. 0) then
		   Tdamp = stadiumTemp
		endif
					
		location=index(Damping_mode,'Atom')
		if (location .gt. 0) then
		    Tdamp = Tatom
		endif

		location=index(Damping_mode,'OFF')
		if (location .gt. 0) then
		    Tdamp = requiredTemp
		endif
!-	damp parameter for Marder thermostat		
 		dampCoeff = BerendsenDampCoeff(requiredTemp, 
     $			Tdamp, langevinCoeff, rampValue)
	    endif

          endif
       endif


CC     Initial newdispl vector
	do j =1, ndf
	     newdispl(j)=0.0d0
cc	     currentCoord(j)=0.0d0
	enddo
        

CC--Updating for different thermostat
        if (selectThermostat.eq.1)then                                 !Nose Hoover
	   dampCoeff = zeta
	   NdampCoeff(iAtom) = dampCoeff
	   dampfactor(iAtom) = EXP(-dthalf*dampCoeff)
!-	  Repeat for each dimension
	   do j =1, ndf2D
!-           Intemediate velocity (same as lammps nvt)
	     velocity(j, iAtom) = oldVelocity(j, iAtom)*
     &	     	dampfactor(iAtom)+
     &		dthalf*atomForce(j,iAtom) / atomMass
!-	     update displacement
	     newdispl(j) = velocity(j,iAtom) * timestep
	   enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc--JS: override to do NVE update for H atom if HNVEFlag==1	   
	   if((atomSpecie(iAtom).eq.2).and.(HNVEFlag.eq.1)) then
	     do j=1,ndf2D 
	       velocity(j, iAtom) = oldVelocity(j, iAtom)
     &		+dthalf*atomForce(j,iAtom) / atomMass
cc--Store the intemediate V for final integration
	       TempOldVelocity(j,iAtom)=velocity(j, iAtom)	 
	       newdispl(j) = velocity(j,iAtom) * timestep

	     enddo
	   endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	       
        endif
        
        if (selectThermostat.eq.3)then                                   !Marder
!-	store Tdamp for final integration	
          NTdamp(iAtom)=Tdamp        
          
	  do j=1, ndf2D
CC--Intemediate V (same as lammps nvt)
	    velocity(j, iAtom) = oldVelocity(j, iAtom)+
     &      +dthalf*(atomForce(j,iAtom)+DampForce(j,iAtom))/atomMass
cc--Updating Disp
            newdispl(j) = velocity(j,iAtom) * timestep
	  enddo	  
	  
	endif
	

        if (selectThermostat.eq.2) then                                    !Langevin
CC--JS: Real LangevinCoeff!
	  dampCoeff=langevinCoeff*LVscaleRatio
	  NdampCoeff(iAtom) = dampCoeff
!-	store Tdamp for final integration	
          NTdamp(iAtom)=Tdamp
                
	  do j=1, ndf2D
CC--Intemediate V
	    velocity(j, iAtom) = oldVelocity(j, iAtom)+
     &      +dthalf*(atomForce(j,iAtom)+DampForce(j,iAtom))/atomMass
cc--Updating Disp
            newdispl(j) = velocity(j,iAtom) * timestep
	  enddo

	endif
	
	
	if(selectThermostat.eq.4) then                                    !Trescale
	     if(stadiumTemp.lt.1E-10) then
	     	write(*,*) 'No rescale with T=0! Exit...'
	     	stop
	     endif
!-	rescale T is exceeds window
	  do j=1, ndf2D
	     if(deltaTemp .gt. Twindow) then
	       velocity(j, iAtom) = oldVelocity(j, iAtom)*
     &		trsfactor
	     endif
	        	  	  
     	  	  
	     velocity(j, iAtom) = velocity(j, iAtom)
     &		+dthalf*atomForce(j,iAtom) / atomMass
	     newdispl(j) = velocity(j,iAtom) * timestep	 
	  enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc--JS: override to do NVE update for H atom if HNVEFlag==1	   
	   if((atomSpecie(iAtom).eq.2).and.(HNVEFlag.eq.1)) then
	     do j=1,ndf2D 
	       velocity(j, iAtom) = oldVelocity(j, iAtom)
     &		+dthalf*atomForce(j,iAtom) / atomMass
cc--Store the intemediate V for final integration
	       TempOldVelocity(j,iAtom)=velocity(j, iAtom)
	       newdispl(j) = velocity(j,iAtom) * timestep	 
	     enddo
	   endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	endif  
	
	if(selectThermostat.eq.5) then                                    !NVESTAT
!-	NVE velocity & displacement update 
	  do j=1, ndf2D
	     velocity(j, iAtom) = oldVelocity(j, iAtom)
     &		+dthalf*atomForce(j,iAtom) / atomMass
	     newdispl(j) = velocity(j,iAtom) * timestep	 
	  enddo
	endif      


ccccccccccccccccccccccc*H-Stablizer*cccccccccccccccccccccccc
CC--Jun Song: hack for H atom
cc--avoiding Temperature of H being out of control
	if(SimStep .lt. NHrescale) then
	   if(atomSpecie(iAtom) .eq. 2) then 
	     Tatom = getKineticTemp(atomCoord, 
     $	    Velocity, iAtom, atomMass)
     	     if(Tatom>MaxHTemp) then
CC--	The scaling factor of H atom
	       write(*,*)"******H Stablizer******"
     	       Htrsfactor=dsqrt(requiredTemp/Tatom) 
     	       do j=1, ndf2D
     	         velocity(j,iAtom)=velocity(j,iAtom)
     &	     *Htrsfactor
     	       enddo	
     	     endif
     	   endif
     	endif
ccccccccccccccccccccccc*H-Stablizer*ccccccccccccccccccccccccc  	  
	
!-	JS: update atomDisp
	do 30 j=1, ndf2D
	    if (newdispl(j) .gt. maxDispl) then
!-	JS: warning message
	    	write(*,*) 'Warning! DeltaDispl bigger than maxDispl!'
		newdispl(j) = maxDispl
	    else if (newdispl(j) .lt. -maxDispl) then
!-	JS: warning message
	    	write(*,*) 'Warning! DeltaDispl less than -maxDispl!'
		newdispl(j) = -maxDispl
	    endif
!--	    if iAtom is not allowed to move in the j'th direction	
	    if (atomID(j, iAtom) .eq. 1) goto 30

!--	    Increment atom displacements
	    atomDispl(j, iAtom) = atomDispl(j, iAtom) + newdispl(j)
			
	    rmsForce = rmsForce + atomForce(j, iAtom)**2
	    meanSqDispl = meanSqDispl + newdispl(j)**2
cc	    currentCoord(j)=atomCoord(j,iAtom)+atomDispl(j,iAtom)
	    
30	continue

!--Jun Song: put atoms in the simulation box 
cc	    if(currentCoord(3) .gt. z_length) then
cc	    	atomDispl(3, iAtom)=atomDispl(3,iAtom)-z_length
cc	    endif
cc	    if(currentCoord(3). lt. 0.0d0) then
cc	    	atomDispl(3, iAtom)=atomDispl(3,iAtom)+z_length
cc	    endif
	    
          nAtoms = nAtoms + 1
101    continue



!-	JS: update E and F	     
	call getEnergiesAndForces (atomID, atomCoord,ix, f, 
     &       atomDispl, aveDispl, atomForce, atomMass,  
     &       systemEnergy, MoveAtoms, MoveDisl, FullField, 
     &       solveFEM, strainE0,iFem)


!-    JS: Final integration, update V2
      do 102 iAtom = i_initial, i_final
	  if(atomSpecie(iAtom).eq.2)then
	  timestep=timestepHH
	  else
          timestep=timestep1
	  endif

!-    Again, need to use the right timestep for different specie	  
	  dthalf=0.5*timestep
		
!-	  Skip near continuum nodes and pad atoms
          atomMass=amass(atomSpecie(iAtom))*1.0d-24 
	  if (isRelaxed(iAtom) .eq. indexContinuum) goto 102
	  if (isRelaxed(iAtom) .eq. indexPad) goto 102

CC--Updating V2 for different thermostat
        if (selectThermostat.eq.1)then                                 !Nose Hoover
	do j=1, ndf2D
	  velocity(j, iAtom) = velocity(j, iAtom)*
     &	  	dampfactor(iAtom)+dthalf*dampfactor(iAtom)*
     &	  	atomForce(j,iAtom) / atomMass          
	enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc--JS: override to do NVE update for H atom if HNVEFlag==1	   
	   if((atomSpecie(iAtom).eq.2).and.(HNVEFlag.eq.1)) then
	     do j=1,ndf2D 
	       velocity(j, iAtom) = TempOldVelocity(j,iAtom)
     &		+dthalf*atomForce(j,iAtom) / atomMass	 
	     enddo
	   endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	  
	endif


        if (selectThermostat.eq.3)then                                 !Marder

	  if((isDamped(iAtom).gt.0))then
 
	    dampCoeff=langevinCoeff*rampValue
	    do j=1, ndf2D
	      randForce=(-1+2*ZBQLU01(0.0d0))*dsqrt(6*dampCoeff*
     &		atomMass*boltzmannConst*NTdamp(iAtom)/timestep)
cc--JS: override to do NVE update for H atom if HNVEFlag==1	   
	      if((atomSpecie(iAtom).eq.2).and.(HNVEFlag.eq.1)) then
	        randForce=0.0d0
	        dampCoeff=0.0d0
	      endif
cc--JS: update Dampforce for damped atoms	      
	      DampForce(j, iAtom)=randForce
     &	      -dampCoeff*velocity(j,iAtom)*atomMass         
           enddo
            
          endif
cc--JS: final integration: update velocities
          do j=1, ndf2D
           velocity(j, iAtom)=velocity(j, iAtom) 
     &      +dthalf*(atomForce(j,iAtom)+DampForce(j,iAtom))/atomMass
          enddo
                             	  
	endif


        if (selectThermostat.eq.2) then                                    !Langevin
	  dampCoeff=NdampCoeff(iAtom)
	  do j=1, ndf2D
	    randForce=(-1+2*ZBQLU01(0.0d0))*dsqrt(6*dampCoeff*
     $        atomMass*boltzmannConst*NTdamp(iAtom)/timestep)
cc--JS: override to do NVE update for H atom if HNVEFlag==1	   
	    if((atomSpecie(iAtom).eq.2).and.(HNVEFlag.eq.1)) then
	      randForce=0.0d0
	      dampCoeff=0.0d0
	    endif
cc--JS: update Dampforce for damped atoms	      
	    DampForce(j, iAtom)=randForce
     &	      -dampCoeff*velocity(j,iAtom)*atomMass
cc--JS: final update V     
            velocity(j, iAtom)=velocity(j, iAtom) 
     &      +dthalf*(atomForce(j,iAtom)+DampForce(j,iAtom))/atomMass		                           
	  enddo
	  	  
	 endif	

	 
	if(selectThermostat.eq.4) then	                                   !Trescale
	  do j=1, ndf2D
	     velocity(j, iAtom) = velocity(j, iAtom)
     &		+dthalf*atomForce(j,iAtom)/atomMass 
	  enddo 
     	endif

	if(selectThermostat.eq.5) then	                                   !NVESTAT
	  do j=1, ndf2D
	     velocity(j, iAtom) = velocity(j, iAtom)
     &		+dthalf*atomForce(j,iAtom)/atomMass 
	  enddo 
     	endif

	  
102   continue


	nAtoms = 0
	do iAtom = 1, numnp
	if (isRelaxed(iAtom) .ge. 1) then
		do j = 1, ndf2D
			oldAcceleration(j,iAtom) = acceleration(j,iAtom)
			oldVelocity(j,iAtom) = velocity(j,iAtom)
		enddo
		nAtoms = nAtoms + 1
	endif
	enddo

        return
	end
     	  
	


	integer function getSystemSize(atomCoord, atomID, simulationCell)

	use mod_global
	use mod_common

	implicit none
	double precision atomCoord(ndf, *)
	double precision xmin, xmax, ymin, ymax, x, y, size
	integer iAtom, j, flag, atomID(ndf, *), nAtoms
	TYPE(region) simulationCell

	xmax = -1.0d10
	xmin = 1.0d10
	ymax = xmax
	ymin = xmin

	nAtoms = 0
	do 10 iAtom = 1, numnp
		
		if (isRelaxed(iAtom) .lt. 1) goto 10


		nAtoms = nAtoms + 1

		x = atomCoord(1, iAtom)
		y = atomCoord(2, iAtom)

		if (x .gt. xmax) then
			xmax = x
		else if (x .lt. xmin) then
			xmin = x
		endif

		if (y .gt. ymax) then
			ymax = y
		else if (y .lt. ymin) then
			ymin = y
		endif

10	continue

	getSystemSize = 0 

	simulationCell%xmin = xmin
	simulationCell%xmax = xmax
	simulationCell%ymin = ymin
	simulationCell%ymax = ymax
	
	print*, 'nAtoms: ', nAtoms
	return
	end




	double precision function distance(x0, y0, x1, y1)
	implicit none
	double precision x0, y0, x1, y1

	distance = sqrt( (x1-x0)**2 + (y1-y0)**2)

	return
	end

	
	
	
	integer function findMDAtoms(atomCoord, isDamped, simulationCell)
	use mod_global
	use mod_common
	implicit none
	double precision atomCoord(ndf, *), damped_width
	double precision xmin, xmax, ymin, ymax 
	double precision x, y, rampValue, getRamp
	integer iAtom, isDamped(*), numAtoms, numDampedAtoms
	TYPE(region) simulationCell

	numAtoms = 0
	numDampedAtoms = 0
	do 5 iAtom = 1, numnp
		isDamped(iAtom) = -1
5	continue

	do 10 iAtom = 1, numnp

		if (isRelaxed(iAtom) .eq. indexContinuum) goto 10
		if (isRelaxed(iAtom) .eq. indexPad) goto 10
		
		x = atomCoord(1, iAtom)
		y = atomCoord(2, iAtom)
		rampValue = getRamp(atomCoord, iAtom, simulationCell)
     
     		if (rampValue .ge. 1.d-5) then
			isDamped(iAtom) = 1
			numDampedAtoms = numDampedAtoms + 1
		else
			isDamped(iAtom) = 0
		endif
		
		numAtoms = numAtoms + 1

10	continue	

	print*, 'Number of Damped Atoms: ', numDampedAtoms
	findMDAtoms = numAtoms
	return
	end


	double precision function getTemperature(atomCoord,
     &		isDamped, velocity, atomMass) 

	use mod_global
       use mod_poten
	implicit none
	double precision atomCoord(ndf, *), isDamped(*), atomMass,
     &			velocity(ndf, *)
	double precision kineticEnergy, v
	double precision  temperature
	integer iAtom, numAtoms, j


	kineticEnergy = 0.d0
	numAtoms = 0
	do 10 iAtom = 1, numnp
        atomMass=amass(atomSpecie(iAtom))*1.0d-24
		if (isRelaxed(iAtom) .lt. 1) goto 10

C		if (isDamped(iAtom) .eq. -1) goto 10
C		if (isDamped(iAtom) .eq. 0) goto 10

		numAtoms = numAtoms + 1
		do 20 j = 1, ndf
			v = velocity(j, iAtom)
			kineticEnergy = kineticEnergy + 
     &				v *v
! 			print*, iAtom, v , isRelaxed(iAtom)
		
20		continue

10	continue
	

	kineticEnergy = 0.5 * kineticEnergy * atomMass

	temperature = kineticEnergy/(1.5 * boltzmannConst *
     &			numAtoms)

c	print*, 'getBathTemperature: ', kineticEnergy, numDampedAtoms
	getTemperature = temperature

	return
	end

	
	
	
	subroutine getBoxTemp(atomCoord,
     &		velocity, atomMass, simulationCell,
     &		interiorTemp, stadiumTemp, currentTemp) 

	use mod_global
	use mod_poten
	use mod_common

	implicit none
	double precision atomCoord(ndf, *),  atomMass,
     &			velocity(ndf, *)
	double precision  keInterior, keStadium, v, keAtom, keTotal
	double precision  temperature, xmin, xmax, ymin, ymax
	double precision damped_width, rampValue, getRamp
	double precision interiorTemp, stadiumTemp, currentTemp
	integer iAtom, numIntAtoms, numStadiumAtoms, j, numAtoms
	TYPE(region) simulationCell

	keInterior = 0.d0
	keStadium = 0.d0
	
	numIntAtoms = 0
	numStadiumAtoms = 0
	numAtoms = 0
	
	do 10 iAtom = 1, numnp

		if (isRelaxed(iAtom) .lt. 1) goto 10
                atomMass=amass(atomSpecie(iAtom))*1.0d-24 
		rampValue = getRamp(atomCoord, iAtom, simulationCell)
     
     		keAtom = 0.0
		do 20 j = 1, ndf
			v = velocity(j, iAtom)
			keAtom = keAtom + v *v* atomMass
20		continue
     
		if (rampValue .gt. 1.d-5) then
			numStadiumAtoms = numStadiumAtoms + 1
			keStadium = keStadium + keAtom
		else
			numIntAtoms = numIntAtoms + 1
			keInterior = keInterior + keAtom			
		endif
		
		keTotal = keTotal + keAtom
		numAtoms = numAtoms + 1
		
10	continue
	

	keStadium = 0.5 * keStadium 
	keInterior = 0.5 * keInterior 
	keTotal = 0.5 * keTotal 
	
	interiorTemp = keInterior/(1.5 * boltzmannConst *
     &			numIntAtoms)     
	stadiumTemp = keStadium/(1.5 * boltzmannConst *
     &			numStadiumAtoms)
	currentTemp = keTotal/(1.5 * boltzmannConst * 
     &			numAtoms)     

     	return
	end


	

			
	
C	get the temperature of the surrounding region
	double precision function getNearNeighborTemp(iAtom, atomCoord,
     &		atomDispl,  isDamped, velocity, atomMass) 
	use mod_global
	use mod_dynamo
	use mod_poten
	implicit none

	double precision atomCoord(ndf, *), atomMass,
     &			velocity(ndf, *), atomDispl(ndf, *)
	integer isDamped(*), iAtom
	
C--	Local variables	
	double precision kineticEnergy, v, rcut, weight
	double precision x1, y1, z1, x2, y2, z2, dist,  value
	integer totNeighbors, jAtom, nAtoms, i, j, k
	
C--	Functions	
	double precision  CutoffRadius	



	rcut = CutoffRadius(1)
     	x1 = atomCoord(1, iAtom)
	y1 = atomCoord(2, iAtom)
	z1 = atomCoord(3, iAtom)


C	Find the kinetic energy of the surrounding atoms
	nAtoms = 0
	kineticEnergy = 0.d0	
	weight = 0.d0
	do 10 k = 1, NumNeighbors(iAtom)

		jAtom = NeighborList(k, iAtom)
		if (isRelaxed(jAtom) .eq. indexContinuum) goto 10
		if (isRelaxed(jAtom) .eq. indexPad) goto 10
		
! 		if (jAtom .eq. iAtom) then
! C			print*, 'Warning; jAtom = iAtom'
! 		endif

		x2 = atomCoord(1, jAtom)
		y2 = atomCoord(2, jAtom)
		z2 = atomCoord(3, jAtom)
              atomMass=amass(atomSpecie(jAtom))*1.0d-24 
		dist = sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
		do 20 j = 1, ndf
			v = velocity(j, jAtom)
			value = v*v* atomMass
			kineticEnergy = kineticEnergy + value
20		continue
		nAtoms = nAtoms + 1
10	continue
	
C	Add the KE of the atom itself
       atomMass=amass(atomSpecie(iAtom))*1.0d-24 
	do 30 j = 1, ndf
		v = velocity(j, iAtom)
		kineticEnergy = kineticEnergy + v *v* atomMass
30	continue
	nAtoms = nAtoms+ 1
	
	
C	Get the temperature
	kineticEnergy = 0.5 * kineticEnergy 	
	getNearNeighborTemp = kineticEnergy/(1.5 * boltzmannConst *
     &			nAtoms)

!      	write(6,9) iAtom, nAtoms
! 9	format('iAtom:', 1x, i5, 2x, 'neighboring atoms:', 1x, i5)	
	return
	end
	
	

	
	
	subroutine setAtomVelocity(atomCoord, isDamped, velocity, 
     &	 atomMass, requiredTemp, rseed)

	use mod_global
        use mod_parallel
	use mod_poten
	implicit none
	double precision atomCoord(ndf, *), velocity(ndf,*),
     &	atomMass,  kineticEnergy,
     &	requiredTemp, v, tmpValue, dfVelocity, aveVelocity(ndf),
     &  totmass 
	integer iAtom, j, isDamped(*)
	integer numAtoms, rseed
	double precision  ZBQLU01, ZBQLUAB,random, getRandomNumber 
	
	numAtoms = 0

	do 5 j = 1,3
		aveVelocity(j) = 0.d0
5	continue

	totmass=0.d0

C--	Initialize random number generator
!	call ZBQLINI(2478945)
        call ZBQLINI(rank+rseed)
        print*,rank,'random seed = ',rank+rseed

	kineticEnergy = boltzmannConst * requiredTemp
cc	v = sqrt(kineticEnergy/atomMass)
	
	do 10 iAtom = 1, numnp
           atomMass=amass(atomSpecie(iAtom))*1.0d-24
	   if (isRelaxed(iAtom) .eq. indexContinuum) goto 10
	   if (isRelaxed(iAtom) .eq. indexPad) goto 10

	   v = sqrt(kineticEnergy/atomMass)			
		do 20 j = 1, ndf 
		  random = getRandomNumber()
		  dfVelocity = v * random
		  velocity(j, iAtom) = dfVelocity					
		  aveVelocity(j)=aveVelocity(j)+dfVelocity*atomMass
		  totmass=totmass+atomMass
20		continue
		
		numAtoms = numAtoms + 1

10	continue
	

	do 30 j = 1, ndf
		aveVelocity(j) = aveVelocity(j)/totmass
		write(6,29)j, aveVelocity(j)
29		format('Average velocity[', i1, ']:',  2x, 1pe11.4)	
30	continue		
	

C	Subract Average velocities
	do 40 iAtom = 1, numnp
	
		if (isRelaxed(iAtom) .eq. indexContinuum) goto 40
		if (isRelaxed(iAtom) .eq. indexPad) goto 40
	
     		do 45 j = 1, ndf				
			velocity(j, iAtom) = velocity(j, iAtom) - aveVelocity(j)
45		continue
	
40	continue	

	print*, 'setAtomVelocity: Total MD Atoms: ', numAtoms
	return
	end




	
	double precision function getRandomNumber() 
	implicit none
	double precision random,  ZBQLU01
	
		random = -1.d0 + 2.d0 *ZBQLU01(0.0D0)
		
		if (random .gt. 0) then
			random = 1.d0
		else if (random .lt. 0) then
			random = -1.d0
		else if (dabs(random) .lt. 1.0d-6) then
			random = 0.d0
		endif
		
		getRandomNumber = random				
	return
	end
	

		

	
	
	double precision function getKineticTemp(atomCoord, velocity, 
     &	iAtom, atomMass)
	use mod_global
       use mod_poten
	implicit none
	double precision temp, ke, atomMass,
     &	atomCoord(ndf, *), velocity(ndf, *), v, x, y
	integer j, iAtom
	

	ke = 0.d0
	do j = 1, ndf
		v = velocity(j, iAtom)
		ke = ke + v*v
	enddo	
       atomMass=amass(atomSpecie(iAtom))*1.0d-24
	ke = ke * 0.5 * atomMass
	temp = ke/(1.5 * boltzmannConst)

	getKineticTemp = temp
	
	return
	end
	
	
	


	double precision function getRamp(atomCoord, iAtom, simulationCell)
	use mod_global
	use mod_common
	implicit none
	double precision atomCoord(ndf, *), xmin, xmax, ymin, ymax
	integer iAtom
	double precision v1, v2, v3, v4, x, y, damp_width, minValue
	TYPE(region) simulationCell
	
	
	x = atomCoord(1, iAtom)
	y = atomCoord(2, iAtom)

	v1 = abs(x - simulationCell%xmin)
	v2 = abs(x - simulationCell%xmax)
	v3 = abs(y - simulationCell%ymin)
	v4 = abs(y - simulationCell%ymax)


	getRamp = minValue(v1, v2, v3, v4)
	if (getRamp .gt. simulationCell%damped_width) then
		getRamp = 0.d0
	else
		getRamp =  1.d0 - getRamp/simulationCell%damped_width
	endif
C	getRamp = 1.0

c	print*, 'getRamp: ', getRamp

	return
	end


	
	

	double precision function BerendsenDampCoeff(requiredTemp, kineticTemp, 
     &	langevinCoeff, rampValue)
		use mod_global
		implicit none
		double precision  langevinCoeff, requiredTemp, kineticTemp,  
     &		 damped_width, deltaT, T0, T, rampValue
				

	deltaT = 1.d0/20.d0 
	T = kineticTemp
	T0 = requiredTemp
		
	BerendsenDampCoeff  = langevinCoeff *
     &				(T- T0)/
     &				sqrt(T**2 + (deltaT)**2) *  rampValue

  	
	return
	end
	
	
	double precision function NoseHooverCoeff(requiredTemp, kineticTemp)
	use mod_global
	implicit none
	double precision  requiredTemp, kineticTemp, zetaDot, OmegaE, Q,
     &		prefactor, unitpico
	integer numDampedAtoms

cc	OmegaE = 8.0e14		!Einstein Frequency s^-1	
cc	prefactor = 0.005
cc	Q = boltzmannConst * requiredTemp/OmegaE**2
cc	zetaDot = 3.d0/Q * boltzmannConst *(kineticTemp - requiredTemp)
cc 	NoseHooverCoeff = zetaDot * prefactor

	unitpico=1.0e-12
	if(NHDampCoeff .lt. 0.0d0) then
	  write(*,*)"Nose Hoover dampcoeff less than 0!!!"
	  stop
	endif

	NoseHooverCoeff=(kineticTemp/requiredTemp-1.0d0)/unitpico**2
 	NoseHooverCoeff=NoseHooverCoeff/NHDampCoeff**2		
	
	return
	end
	
	
C	minimum of four real numbers, modified to calculate the minimum of first
C	two
	double precision function minValue(a, b, c, d)
	implicit none
	double precision a, b, c, d, value(4)
	double precision min_val
	integer i

	value(1) = a
	value(2) = b
	value(3) = c
	value(4) = d
	
	do i=1,4
		if (i .eq. 1) then
			min_val = value(1)
		else
			if (value(i) .lt. min_val) then
				min_val = value(i)
			endif
		endif	
	enddo

	minValue = min_val
	return
	end


	


	
	subroutine BoxMuller(atomCoord, isDamped, velocity, 
     &	 atomMass, requiredTemp)
     	use mod_global
	implicit none
	
	double precision atomCoord(ndf, *), velocity(ndf,*),
     &	aveVelocity(ndf), atomMass,  requiredTemp, v, dfVelocity 
	integer iAtom, j, isDamped(*), numAtoms
	double precision  ZBQLU01, r1, r2, sigma, random, pi
	
	
	
	pi = dacos(-1.d0)
	
	do 5 j = 1,3
		aveVelocity(j) = 0.d0
5	continue


C	Initialize random number generator with a chosen random seed
	call ZBQLINI(2478945)
	
C	The variance in velocities	
	sigma = sqrt(boltzmannConst * requiredTemp/atomMass)

	do 10 iAtom = 1, numnp

		if (isRelaxed(iAtom) .eq. indexContinuum) goto 10
		if (isRelaxed(iAtom) .eq. indexPad) goto 10
			
		do 20 j = 1, ndf 
		
C	Get two random numbers	
			r1 = ZBQLU01(0.0D0)
			r2 = ZBQLU01(0.0D0)
C	Derive the random number used			
			random = sqrt(-2.0 * log(r1)) * cos(2.0 * pi * r2)
			
			dfVelocity = sigma * random
			velocity(j, iAtom) = dfVelocity
			aveVelocity(j) = aveVelocity(j)  + dfVelocity 
20		continue
		
		numAtoms = numAtoms + 1

10	continue
	

	do 30 j = 1, ndf
		aveVelocity(j) = aveVelocity(j)/dfloat(numAtoms)
		write(6,29)j, aveVelocity(j)
29		format('Average velocity[', i1, ']:',  2x, 1pe11.4)	
30	continue		

		
C	Subract Average velocities
	do 40 iAtom = 1, numnp
	
		if (isRelaxed(iAtom) .eq. indexContinuum) goto 40
		if (isRelaxed(iAtom) .eq. indexPad) goto 40
	
     		do 45 j = 1, ndf				
			velocity(j, iAtom) = velocity(j, iAtom) - aveVelocity(j)
45		continue
	
40	continue	

	print*, 'BoxMuller: Total MD Atoms: ', numAtoms	
     	return
	end

	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! New subroutine for calculate the average stresses !!!!!
!!!!! Based on the subroutine hstress - Jun Song comment!!!!!

      subroutine  Nhstress(tempAvgStress,Pnump,PIsRelaxed,Pvirst)
      implicit none
      integer i, j, k, Pnatoms
      integer Pnump
      integer PIsRelaxed(Pnump)
      double precision Pvirst(3,3,Pnump)
      double precision Phydro_stress
      double precision Pavg_stress(3,3)
      double precision Ptstress
      double precision tempAvgStress(3,3)
      
  

      Phydro_stress=0.0
      Pavg_stress(1:3,1:3)=0.0
      Pnatoms=0
      do i=1,Pnump
        if (PIsRelaxed(i).eq.1) then
          Pnatoms=Pnatoms+1
          do j=1,3
            Phydro_stress=Phydro_stress+Pvirst(j,j,i) 
            do k=1,3
              Pavg_stress(j,k)=Pavg_stress(j,k)+Pvirst(j,k,i)
            enddo
          enddo
        endif
      enddo

      Pavg_stress(1:3,1:3)=Pavg_stress(1:3,1:3)/Pnatoms
      Phydro_stress=Phydro_stress*1.0/3.0/Pnatoms
     
      Ptstress=0.5*(Pavg_stress(1,1)+Pavg_stress(3,3))
      
  
      tempAvgStress(1:3,1:3)=Pavg_stress(1:3,1:3)
      	  
      return
      end 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Jun Song comment end			!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The subroutine called by ma06 to do MD steps in cgma05.f!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
! ---- include timeplot for plot view

	subroutine dosteps(n, atomDispl, atomForce,db,fl,epps,
     &	  iprint, dsmax, rseed, dfn, atomID, atomCoord,ix,f,itx,
     &	  CheckSlip,AddedSlip,LostSlip,MoveDisl, MoveAtoms) 

	use mod_global
	use mod_poten

	use mod_dynamo
	use mod_output
        use mod_material
        use mod_parallel

        use mod_timming
	use mod_common	
	
	use mod_file

	implicit none
	include '../Disl/disl_parameters.par'
	integer maxAtoms
	parameter (maxAtoms = 100000)
	
	integer n,ix(*),atomID(ndf,*),iprint,itx(*),rseed,NumDis, numH
	double precision atomDispl(ndf, numnp),db(n), atomForce(3,*),fl,
     &		epps,dsmax,dfn,atomCoord(ndf, *),f(*), mass
	logical ConvOnFn,AddedSlip,CheckSlip,MoveDisl,LostSlip,
     &     	MoveAtoms, dislcheck
	logical dislpass
	logical Printing, debug, plot, FullField, solveFEM, reStartAveraging
	common/debugger/debug
	common /MD/ oldAcceleration, oldVelocity, acceleration, 
     &		velocity, aveDispl, rampLevel,numMDAtoms

	
!--	Local variables
	integer isDamped(maxAtoms)
	integer Nsteps, i, iAtom, j, rc,intAtom1,intAtom2,intAtom3,intAtom4
        integer ii,jj
	integer numMDAtoms, FEMSteps, ndis_checked
	integer MDSteps, nnsteps, MaxMDSteps
     	integer FEMStepMin, FEMStepRange, FEMStepCounter, iStep,  
     &		readNeighbors
!--	Indexes for H atoms, using the fact that they are adjacent
     	integer :: HIndex_Init=0,HIndex_Final=-1,Indexmin,IndexMax
     	integer :: NumInterstitial=0

	double precision atomMass, damped_width, 
     &		langevinCoeff, LCNHDampCoeff, requiredTemp, 
     &		currentTemp, picoSecond, systemEnergy,  timestep,
     &		TInterior, stadiumTemp,atomTemp	
	
	double precision
     &		oldAcceleration(3,maxAtoms), oldVelocity(3,maxAtoms), 
     &		acceleration(3,maxAtoms), velocity(3,maxAtoms), 
     &		aveDispl(3,maxAtoms), rampLevel(maxAtoms)
     
	double precision xmin, xmax, ymin, ymax, plottime
	character *80 atomFileName, energyFileName, tempFileName, 
     &		neighborFileName

        logical newMD, finishMD
        
        !!! JS comment For output Avgstress !!!
        double precision tempAvgStress(3,3)
	integer iFem

	character*80 filename
	integer logic, Nstepsorig, npass
	double precision :: dtol
		
	Type(region) SimulationCell
	Type(MD_Thermostat) Thermostat
	

!--	Functions
	integer getSystemSize, findMDAtoms
	double precision getTemperature, random, ZBQLU01, getRamp,strainE0
	data MDSteps /0/
	data nnsteps /0/
	data intAtom1 /0/
	data intAtom2 /0/
	data intAtom3 /0/
	data intAtom4 /0/
        data FEMStepCounter /0/
        data ndis_checked /0/

	npass = 0

        if (maxAtoms.lt.numnp) then
          print*,'increase size of maxAtoms in md.f'
          stop
        endif
	
	plot = .false.
	reStartAveraging = .false.
	
	call lostslipinit(LostSlip)
	
	picoSecond = 1.0d-12	! s
	atomMass = amass(1) * picoSecond**2	! eV s^2/A^2
	timestep = 1.0d-15			! seconds
	FEMSteps = 1
	FullField = .true.
	FEMStepMin = 5
	FEMStepRange = 0
	
	dtol = 1.d-3

	print*, 'Entering dosteps'
		
!	Read Data
CC--Jun Song NumMDRescale is # of MD steps 
CC--doing temperature rescaling
CC--Reading neighborlist update parameter			
	open(unit=200, file='md.inp', status='old')
	read(200,*) damped_width 
	read(200,*) langevinCoeff, LVscaleRatio
	read(200,*) LCNHDampCoeff
	read(200,*) requiredTemp
        read(200,*) timestep1
        read(200,*) indextimeH
	read(200,*) FEMStepMin
	read(200,*) Nsteps 
	read(200,*) Thermostat%Type
	read(200,*) Thermostat%Damping_Mode
	read(200,*) NumMDRescale
	read(200,*) Twindow
	read(200,*) NHrescale, MaxHTemp, HNVEFlag 
	read(200,*) MaxMDSteps
	close(200)


!! Jun Song comments: output step, energy and temperature per steps 
!!(specified when writing to the file)
	open(5800,file='MDlog.CADD',status='unknown')
!! Jun Song comment: Initialization

!	Write Data	
	write(6,4) 'damped_width: ', damped_width
	write(*,*) 'Base langevinCoeff: ', langevinCoeff
	write(*,*) 'LangevinCoeff scale ratio', LVscaleRatio
CC--Jun Song: Care with Langevin Coefficient!!
	write(*,*) 'Langevin Coefficient', langevinCoeff*LVscaleRatio
	write(6,4) 'requiredTemp: ', requiredTemp
        write(6,4) 'timestep1: ', timestep1
        write(6,4) 'timestepH: ', timestep1/indextimeH
	write(6,6) 'FEMStepMin: ', FEMStepMin
	write(6,6) 'Nsteps: ', Nsteps
	write(6,4) 'atomic mass: ', atomMass
	write(6,8) 'Thermostat: ', Thermostat%Type
	write(6,8) 'Damping Mode: ', Thermostat%Damping_Mode
	write(*,*) '# T Rescale MD steps: ', NumMDRescale
	write(*,*) 'The windown for T rescale ', Twindow
	write(*,*) '**********For H atom only**********'
	write(*,*) 'H stablizer Steps',NHrescale,'for T>',MaxHTemp
	write(*,*) 'Exclude H from Thermostat? ', HNVEFlag
4	format(a16, 2x, 1pe11.4)		
6	format(a16, 2x, i5)
8	format(a16, 2x, a16)

	finishMD = .false.
	newMD = .false.
        timestep=timestep1
       
CC--Jun Song: set SysTemp and CUTFACT to user defined value       
CC--	set Nose Hoover damping coefficient as from input
	NHDampCoeff=LCNHDampCoeff
	SysTemp=requiredTemp
cc--	Stop if CUTFACT value if too small or too big
	if(CUTFACT .lt. 1.0d0 .or. CUTFACT .gt. 2.0d0) then
	  write(*,*) "*****************************************"
	  write(*,*) "Neigh-cutoff not suitable! Set to default"
	  write(*,*) "*****************************************"
	  stop
	endif
       
	print*, 'MDSteps: ', MDSteps
	if (MDSteps .eq. 0) newMD = .true.
cc	if (MDSteps .eq. MAXMDSteps) finishMD = .true.

	UseRescale=.false.
	if (MDSteps .lt. NumMDRescale) then
		UseRescale=.true.
	endif
	
	if (finishMD) goto 500	! Finish MD simulation

CC--Here comes inilization if newMD	
	if (newMD) then
	      SimStep=0
              allocate(avevirst(3,3,numnp))


!      Find the atom# of two interface atom along crack line
            do iAtom=1,numnp
             if(isRelaxed(iAtom).eq.indexInterface.and.
     &          atomCoord(2,iAtom).eq.0)then
                if(atomCoord(1,iAtom).lt.0)then
                   intAtom1=iAtom
                endif
                if(atomCoord(1,iAtom).gt.0)then
                   intAtom2=iAtom
                endif
             endif
             if(isRelaxed(iAtom).eq.indexInterface.and.
     &          atomCoord(1,iAtom).eq.0)then
                if(atomCoord(2,iAtom).lt.0)then
                   intAtom3=iAtom
                endif
                if(atomCoord(2,iAtom).gt.0)then
                   intAtom4=iAtom
                endif
             endif
           enddo

!         Get System Size information
          rc = getSystemSize(atomCoord, atomID, SimulationCell)
          SimulationCell%damped_width = damped_width

          numMDAtoms = findMDAtoms(atomCoord, isDamped, SimulationCell)


         print*,'Simulation Cell Dimensions:'
         write(6,9) 'xmin:', SimulationCell%xmin, 
     &		'xmax:', SimulationCell%xmax
	write(6,9) 'ymin:', SimulationCell%ymin,
     &		'ymax:', SimulationCell%ymax
9	format(2(a6, 2x, f11.4,2x))
	print*, 'numMDAtoms: ', numMDAtoms


!    --Initialize data
		do iAtom = 1, numnp
			if (isRelaxed(iAtom) .eq. indexAtom .or.
     &			isRelaxed(iAtom) .eq. indexInterface ) then		
				velocity(1:ndf, iAtom) = 0.d0
				oldVelocity(1:ndf, iAtom) = 0.d0
				oldAcceleration(1:ndf, iAtom) = 0.d0
				aveDispl(1:ndf, iAtom) = 0.d0
				rampLevel(iAtom) = 0.d0
				DampForce(1:ndf, iAtom)=0.d0
			endif
		enddo	


!!	Set atom velocities 
CC	Only do this if newMD
		call setAtomVelocity(atomCoord, isDamped, velocity, 
     &	 	atomMass, requiredTemp * 2.d0,rseed)
	
	
!!	Get Ramp Value
		do iAtom = 1, numnp
			if (isRelaxed(iAtom) .eq. indexAtom .or. 
     &			isRelaxed(iAtom) .eq. indexInterface ) then
			rampLevel(iAtom) = getRamp(atomCoord, iAtom, simulationCell)
			endif
		enddo

	endif	
!     end of if new md loop
CC--New MD initialization ends


!	Assign oldVelocity = velocity
        oldVelocity(1:ndf, 1:numnp) = velocity(1:ndf, 1:numnp)
        call getBoxTemp(atomCoord, velocity, atomMass, 
     &     simulationCell, TInterior, stadiumTemp, 
     &     currentTemp) 

        print*, 'Total Temperature:', currentTemp,'iT ',TInterior
     &     ,'Stemp ',stadiumTemp 

!	Main Loop
!--	get the min and max of H indexes
        Indexmin=numnp
        Indexmax=1
        
        FEMSteps=FEMStepMin
        do iAtom = 1, numnp
          if(atomSpecie(iAtom).eq.2)then
	    numH=iAtom
            write(*,*)"*********numH*********",numH
            if(iAtom .gt. Indexmax) Indexmax=iAtom
            if(iAtom .lt. Indexmin) Indexmin=iAtom
	  endif
        enddo

        
!--	when there is H atom. Set HIndex_Init and Final
	if(Indexmax .ge. Indexmin) then
	  HIndex_Init=Indexmin
	  HIndex_Final=Indexmax
	  NumInterstitial=Indexmax-Indexmin+1
!	Output the min and max of H indexes
	  write(*,*) "******HIndex_Init is ", HIndex_Init
	  write(*,*) "******HIndex_Final is ", HIndex_Final
	endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CC--Jun Song: No multiple timestep if no interstitial!
	if((NumInterstitial.le.0).and.(IndextimeH.gt.1)) then
	  write(*,*)"Only use mult-timestep with >0 intersitial!"
	  stop
	endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    

	iFem = 0
	dislpass = .false.
	Nstepsorig = Nsteps
	do 10 iStep = 0, Nsteps
CC--Increment SimStep for each simulation step
	    SimStep=SimStep+1
c$$$	    if (iStep == Nsteps-f) then 
c$$$	       MoveDisl = .true.
c$$$	    else
c$$$	       MoveDisl = .false.
c$$$	    endif
	    if(mod(SimStep, indextimeH).eq.0)then
!	    write(*,*)"iStep_: ",iStep
	       i_initial=1
	       i_final=numnp
	    else 
	       i_initial=HIndex_Init
	       i_final=HIndex_Final
	    endif

 
          if(mod(iStep,101).eq.0)then
             write(*,*)"Simlation Step: ",iStep
          endif


          if (nnsteps .eq. 0) then
		!! At the very first step, aveDispl = atomDispl for all nodes
	      aveDispl(1:ndf, 1:numnp) = atomDispl(1:ndf, 1:numnp)
	      avevirst(1:3,1:3,1:numnp) = 0.d0
	  endif
	
	   !!Find Average positions of free and interface atoms
	  if (mod(FEMStepCounter, FEMSteps) .eq. 0) then	
     	     iFem = iFem + 1
	     solveFEM = .true.
	     restartAveraging = .true.
	     print *, 'Ifem = ', iFem
	  else
c$$$	     if (dislpass) then 
c$$$		solveFEM = .true. 
c$$$		restartAveraging = .true.
c$$$		iFem=1
c$$$	     else
		solveFEM = .false.
c$$$	     end if
	  endif
c$$$	  if (iStep == NSteps-FEM) then 
c$$$	     solveFEM = .true.
c$$$	  endif
          !! Get Energies and Forces on MD atoms
CC--Jun Song: Get forces and energies for first runs
          if(newMD) then	
          call getEnergiesAndForces (atomID, atomCoord,ix, f, 
     &       atomDispl, aveDispl, atomForce, atomMass,  
     &       systemEnergy, MoveAtoms, MoveDisl, FullField, 
     &       solveFEM,strainE0,iFem)
     	  newMD=.false.
     	  endif

!         !!Get Current temperature of the MD region
          currentTemp = getTemperature(atomCoord,
     &       isDamped, velocity, atomMass) 

!         !! Increment atom positions
          call velocityVerlet(atomCoord, atomDispl, atomForce,
     &		atomID, isDamped, oldAcceleration, acceleration,  
     &		oldVelocity, velocity, timestep, atomMass, 
     &		langevinCoeff, requiredTemp, currentTemp,
     &		SimulationCell, rampLevel, Thermostat, ix, f, 
     &		aveDispl, systemEnergy, MoveAtoms, MoveDisl,
     & 		Fullfield, solveFEM, strainE0,iFem)
     

          !! Restart averaging displacements for MD atoms
 	  if (reStartAveraging .eq. .true.) then
	      do iAtom = 1, numnp	
		if (isRelaxed(iAtom) .eq. indexAtom
     &              .or. isRelaxed(iAtom) .eq. indexInterface) then
			aveDispl(1:ndf, iAtom) = 0.d0
		endif	
	      enddo
			
	      random = ZBQLU01(0.0D0) 
	      FEMSteps = FEMStepMin + int(random  * FEMStepRange)
	      if (FEMSteps.gt.Nsteps)then
	          FEMSteps = Nsteps+1
	      endif
	      !print*, 'FEMSteps: ', FEMSteps
	      FEMStepCounter = 0
	      reStartAveraging = .false.
          endif

	  !!Calculate average displacements for interface and free atoms
     	  do 30 iAtom = 1, numnp
	     if (isRelaxed(iAtom) .eq. indexAtom
     &	      .or. isRelaxed(iAtom) .eq. indexInterface) then	
     	         aveDispl(1:ndf, iAtom) =  
     &		 aveDispl(1:ndf, iAtom) + 
     &		 atomDispl(1:ndf, iAtom)/FEMSteps
	     endif	
30	  continue

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CC--Jun Song Test
cc	  do iAtom=1, numnp
cc	     if(atomSpecie(iAtom).eq.2) then
cc	       write(*,*) "# of Steps and H neighbors", SimStep,
cc     $	       NumNeighbors(iAtom)
cc	     endif
cc	  enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          !! calculate average virial stresses
	  if (iStep.ge.Nsteps/2)then
	     do iAtom = 1,numnp
            atomMass=amass(atomSpecie(iAtom))*1.0d-24
	       if (isRelaxed(iAtom) .eq. indexAtom
     &          .or. isRelaxed(iAtom) .eq. indexInterface) then

                 !! add kinetic term
                 do ii=1,3
                  do jj=1,3
                   virst(ii,jj,iAtom)=virst(ii,jj,iAtom)-
     &               velocity(ii,iAtom)*velocity(jj,iAtom)*atomMass
     &               /((material(1)%a0**3)/4.d0)
                  enddo
                 enddo

	         avevirst(1:3,1:3,iAtom)=(avevirst(1:3,1:3,iAtom)*
     &             (iStep-Nsteps/2)+virst(1:3,1:3,iAtom))
     &             /(iStep-Nsteps/2+1)
	       endif
	     enddo
	  endif

          FEMStepCounter = FEMStepCounter + 1
!	  print *, 'FEMSTepCounter =', FEMStepCounter
          !!Check for emitted Dislocations
          plottime = nnsteps*timestep1*1.d12

          call cpu_time(ct2)
!	  print *, 'Calling dislcheck'
	  if (dislpass) then 
	     filename='out/atom_pass.cfg'
	     call iofile(filename,'formatted  ',logic,.false.)
	     call dump_atom(atomCoord, atomDispl, logic)
	     close(logic)

	     filename='out/atom_pass.vtk'
	     call iofile(filename,'formatted  ',logic,.false.) 
	     call dump_mesh(atomCoord, atomDispl, ix,logic)
	     close(logic)

	  end if
	  dislpass = .false. 
          if (dislcheck(CheckSlip,LostSlip,AddedSlip,MoveDisl,ix,
     &      atomCoord, atomDispl,itx,IsRelaxed,numnp,ndf,
     &      nxdm,numel,nen1,newmesh,plottime,dislpass, npass)) then
	     if (dislpass) then 
		print *,'Dislocation removed from atomistics'
	     end if
c$$$	     if (npass > 1) then 
c$$$		Nsteps = NstepsOrig*2
c$$$		npass = 0
c$$$	     end if

             write(*,*)'time = ', plottime,'ps'
             write(*,*)'Disl checked!'
             MDSteps=MDSteps+1
             write(*,*)'Disl checked by rank',rank
             ndis_checked=ndis_checked+1
!             if (ndis_checked.eq.1) then
!               call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
!               stop
!                goto 500
!             endif
          endif
          call cpu_time(ct3)
          ct5=ct5+ct3-ct2

          call getBoxTemp(atomCoord, velocity, atomMass, 
     &		simulationCell, TInterior, stadiumTemp, 
     &		currentTemp) 

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CC--Jun Song: output step, energy, temperature, stress per steps
          tempAvgStress(1:3,1:3)=0.0
	  if(mod(iStep,indextimeH).eq.0) then
!! Jun Song:  get the average stress of the system
CC	     call Nhstress(tempAvgStress,numnp,IsRelaxed,avevirst)
!! Can output average stresses if needed  
	     write(5800,5801) SimStep, currentTemp, stadiumTemp,  
     $       TInterior, total_energyMD, energyH
	   endif

5801	  format(i8,3f10.3,2f12.4)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
	  nnsteps=nnsteps+1
c$$$	  if (Moved) then 
c$$$c       Output the new_atom config after moving atom_displacements
c$$$	     filename='out/moved_atoms.cfg'
c$$$	     call iofile(filename,'formatted  ',logic,.false.)
c$$$	     call dump_atom(atomCoord, atomDispl, logic)
c$$$	     close(logic)
c$$$	     Moved = .false. 
c$$$	  end if

	  if (MoveMesh) then 
c$$$	     if (Moved) then 
c$$$		Moved = .false. 
c$$$	     end if
	     if (iStep .eq. 1) then 
		call get_crack_tip(atomCoord, atomDispl)
		if (xtip(1) > 0.0d0) then 
		   if (abs(xtip(1)-xtip_init(1)) > (x_move_mesh)) then 
		      do i = 1,2
			 xtip_actual(i) = abs(xtip(i)-xtip_init(i))
			 if (xtip(i) > xtip_init(i)) then 
			    x_tip_dir(i) = 1.0d0
			 else
			    x_tip_dir(i) = -1.0d0
			 end if
		      end do
		      print *, 'xtip_actual = ', xtip_actual, x_tip_dir
		      if (x_tip_dir(1) > 0.0d0) then 
			 if (int(xtip_actual(1)/x_move_mesh)>0) then 
			    call move_atomistic_crack(atomCoord, ix,
	1			 atomDispl)
			    Moved = .true. 
			 end if
		      end if
c$$$		   return
		   end if
		end if
	     end if
	  end if
c$$$	  if (ndisl > 4) then 
c$$$	     if (mod(iStep,10) .eq. 0) then 
c$$$		filename='out/atom_pass2.cfg'
c$$$		call iofile(filename,'formatted  ',logic,.false.)
c$$$		call dump_atom(atomCoord, atomDispl, logic)
c$$$		close(logic)
c$$$	     end if
c$$$	  end if
10      continue

	MDSteps = MDSteps + 1
 
 	print*,'total_energyMD',total_energyMD
	if (finishMD .eq. .false.) return
	
	! Finish MD
        stop

500    continue


	return
	end


