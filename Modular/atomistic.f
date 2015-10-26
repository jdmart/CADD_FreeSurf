c**---------------------------------------------------------------------
c** atomistic.f : atomistic routines
c--

c**-----------------------------------------------------------------------
c** ProcessClump : process clump
c**
c**      Non-Obvious Parameters :
c**          dr (out) : out-of-balance forces
c**            fls (in)  : true = compute forces
c**            flw (in)  : true = compute energy
c**      Algorithm :
c**              Loop through each atom
c**                   And add its contribution
c**              end do
c--

      subroutine ProcessClump(id,x,ix,f,b,dr,fls,flw)
      use mod_global
      use mod_poten
      use mod_dynamo
      use mod_file
      implicit none

c--Variables Transferred

      double precision  b(ndf,*),x(nxdm,*),f(ndf,*), dr(ndf,*)
      double precision total_energy,xdef, ydef, zdef, xini, yini, zini
      double precision bb1, bb2, bb3, aa1, aa2, aa3
      integer  id(ndf,*), ix(nen1,*), no_MDatoms

      logical fls, flw,NeedList


c--Local Variables
      integer irep, iel, i, j, inode, idf
      integer logic
      character*80 filename
      logical Stressflag
      common/flag/Stressflag
c
      if(.not.Stressflag) then
         allocate(virst(3,3,numnp))
         Stressflag = .true.
      endif
      virst(1:3,1:3,1:numnp)=0.d0


c     chkdis checks whether neighbor lists need updating
c
      call chkdis(b,x)
       total_energy=0.d0
       no_MDatoms=0
       total_energyMD=0.d0 
       energyH=0.d0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CC--Jun Song Screen Print if updating neighborlist
      if(newlst .eq. 1) then
        write(*,*)"ccccccccccccccccccccccccccccccc"
        write(*,*)"Updating NeighborList at ", SimStep
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     
      do irep = 1, numnp
         NeedList=IsRelaxed(irep).ne.0
c         NeedList=IsRelaxed(irep).ge.1
         call gneigh(irep,b,x,NeedList)
         if(NeedList) then
c
c     NonLocal adds the energy and forces due to the E_i term (energy of
c     atom i) in the total energy functional.
c


            call NonLocal(id,x,ix,f,b,dr,fls,flw,irep)

        if (IsRelaxed(irep).eq.1) then
	total_energyMD=total_energyMD +energy(irep)

CC--MS: Get H displacements and E
      if (atomSpecie(irep).eq.2) then
        bb1=b(1,irep)
        bb2=b(2,irep)
        bb3=b(3,irep)
        aa1=x(1,irep)
        aa2=x(2,irep)
        aa3=x(3,irep)
        energyH=energyH+energy(irep)
        end if 
       no_MDatoms=no_MDatoms+1
         end if
        end if
      end do
c       print*, h_no
c          xdef =(x(1,h_no)+ b(1,h_no))/95.5171
c           ydef =(x(2,h_no)+ b(2,h_no))/98.2759
c           zdef = (x(3,h_no) + b(3,h_no))/4.999  
c          xini =(x(1,h_no))/95.5171
c           yini =(x(2,h_no))/98.2759
c           zini = (x(3,h_no))/4.999  
c       print*,'poz_H',xdef,  ydef, zdef
c       print*,'poz_H1',xini, xini,xini
c       print*,'poz_b',bb1, bb2,bb3
c       print*,'poz_a',aa1, aa2,aa3
c      print*,'total_energy', total_energy
c      print*,'total_energy',total_energy
c      
c	print*,'total_energyMD',total_energyMD
c	print*,'no of MD atoms', no_MDatoms
c	stop

c      print*, 'energy of H atom', energyH
c
c     anything involving "numnpp1" is hard-wired to the brinell indenter
c     and is written so that it will be "turned off" if numnpp1=-1,
c     which is the initial value.
c
      if(numnpp1.lt.numnp) return
      call gneigh(numnpp1,b,x,.true.)
      call NonLocal(id,x,ix,f,b,dr,fls,flw,numnpp1)
      return
      end


	subroutine ProcessClumpWrapper(id,x,ix,f,b,dr,Eatom)
      use mod_global
      use mod_poten
      use mod_dynamo
      use mod_file
  	implicit none

      double precision  b(ndf,*),x(nxdm,*),f(ndf,*), dr(ndf,*), Eatom(*)
      integer  id(ndf,*), ix(nen1,*), iatom
      logical fls, flw

	 fls = .true.
	 flw = .true.

	 call ProcessClump(id,x,ix,f,b,dr,fls,flw)
	 do iatom = 1, numnp
		Eatom(iatom) = energy(iatom)
	 enddo
 

	return
	end


