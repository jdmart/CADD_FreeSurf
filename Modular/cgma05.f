      
      
      subroutine vafuncMD(id,x,ix,f,b,dr,totener2,MoveAtoms,
     &  MoveDisl, FullField,strainE0, iFem)
      use mod_global
      implicit none
c     
      double precision b(ndf,*),f(ndf,*),dr(ndf,*),x(nxdm,*),totener2
      integer id(ndf,*),ix(4,*),idf,iFem
      logical MoveAtoms,MoveDisl, FullField
c     
c---- Calculate out-of-balance force residual, strain energy, external
c     work,
c---- and total potential energy associated with mesh
c---- 
c     
      integer i,j
      double precision ftot,sed,work,etot,edd,StrainEnergy
      double precision totener,strener,strainE0
      common/enerdat/ totener,strener
c     
c     Calculate energy and force residual
c     
c     
c     apply boundary conditions
c     
      do i=1,ndf
!JM   loop adjustment made to include indenter atoms
!JM   and BC's applied to them. numnp vs numnpp1
!         do j=1,numnp
         do j=1,numnpp1
            if (id(i,j).eq.1) b(i,j) = time*f(i,j)
         end do
      enddo
c     
c     
c     

      edd=0
      call InitialiseEnergy(.true.,.true.,id,dr,f)
!      print *, 'Ifem in VafuncMD is', iFem
      call fem_move_pad(x,b,ix,f,time,z_length,id,IsRelaxed,edd,
     &     dr, FullField, MoveDisl,strainE0,numel,avevirst,SysTemp,iFem,
     $     Moved)
c     
c     initialize
c
CCCCC
	if(FullField) then
	  do i = 1, numnp
	    if(IsRelaxed(i).eq.2) dr(1:ndf,i)=0.0d0
	  enddo
	else
          call InitialiseEnergy(.true.,.true.,id,dr,f)
	endif
c     
c     atomistic energy and forces.
c
      return
      end


********************************************************************
      subroutine ma06(id,x,ix,f,b,dr,db,input,itx)
      use mod_global
      use mod_dislocation
      use mod_timming
      implicit none
c
      double precision x(nxdm,numnp),b(ndf,numnp),dr,db,f
      integer id(ndf,numnp),ix(nen1,numel),itx
      logical AddedSlip,LostSlip, MoveDisl
      character*80 input,filename
c
      integer lower,upper,next,maxfn,idum,iprint,n,i,j,iat,ii(3),k
     $     ,maxfn2,natoms,iprint2,rseed
      double precision dum,dsmax,dsmax2,dfn,tolm,ener,LENGTHTOL,xc(2)
     $     ,xx(2),CutoffR2,rcutsq,told
      parameter(LENGTHTOL=1.d-4)
      integer minflag
      common/cgflag/minflag
      Character(len=16) :: MD_Status

      logical debug, md
      common/debugger/debug



c     calling format:   ma06,,rseed,,
c     rseed:   integer, seed for random # generator
c
 10   continue

      call cpu_time(ct1)

      lower = 4
      upper = next(lower,input)
      call freein(input,lower,upper,rseed,dum,1)
      lower = upper
      upper = next(lower,input)
      call freein(input,lower,upper,idum,dsmax,2)
      lower = upper
      upper = next(lower,input)
      call freein(input,lower,upper,iprint,dum,1)

      idtemp=.false.
      do i=1,numnp
         if(IsRelaxed(i).le.0) then
            idtemp(1:ndf,i)=.true.
         else
            do j=1,ndf
               idtemp(j,i)=(id(j,i).eq.1)
            enddo
         endif
      enddo
      debug=.false.

                call dosteps(n,b,dr,db,ener,tolm,iprint,
     $          dsmax,rseed,dfn,id,x
     $          ,ix,f,itx,.true.,AddedSlip,LostSlip,
     $          .true.,.true.)



        call cpu_time(ct3)
        ct6=ct6+ct3-ct1
        print*,' '
        print*,'****amount of cpu time devoted to various tasks '
        print*,' in detection band = ',ct5
        print*,' in md = ',ct4
        print*,' in fem = ',ct7
        print*,' in ma06 = ',ct6
        print*,' '
        if (Moved) then 
           print *, 'Crack tip moved to ', xtip_init(1)
           print *, 'Repeating calculation without increasing load'
c$$$           Moved = .false. 
c$$$           goto 10
        end if

      return
      end

