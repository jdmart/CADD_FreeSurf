c**---------------------------------------------------------------
c**    PROGRAM FEAP:
c**      executes the subprograms "feap", "macros" or "stop"
c**
c**   Sub-Programs:
c**      feap: reads global settings, allocates storage and calls pmesh
c**            where material, grain, model and constitutive info is read.
c**      macros: calls the solution macros through the control
c**              subroutine pmacr1
c--
c
      program feap
      use mod_main
      use mod_global
      use mod_boundary
      use mod_dynamo
      use mod_parallel
      implicit none
      character*4 wd(7)
      character*80 input
      data wd/'feap','macr','stop','xxxx','xxxx','xxxx','xxxx'/
c

      integer neqad,l
      integer i
      double precision rcutsq,CutoffR2
c
!      call MPI_INIT(ierr)
!      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
!      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
       nprocs=1
       rank=0
       ierr=0

c
c---- variable definition
c
c---- numnp  = number of nodal points
c---- numel  = number of elements
c---- ndm    = dimension of ambient space
c---- nxdm   = dimension of x array (=3)
c---- ndf    = number of degrees of freedom per node
c---- nen    = number of nodes per element
c---- nen1    = number of nodes per element + 1
c---- nad    = added size to element matrices in excess of ndf*nen
c---- nsdm   = dimension of stress array
c---- nquad  = number of quadrature points per element
c---- neq    = number of equations of system
c---- nq     = dimension of element internal variable array
c---- nshp   = dimension of element shape function array
c---- nxsj   = dimension of element jacobian array
c
c
c---- write banner
      write(6,'('' f i n i t e   e l e m e n t   a n a l y s i s'',
     1  ''   p r o g r a m'')')
c
c---- read a card and compare first 4 columns with macro list
1     read(5,'(a80)') input
      if (input(1:4).eq.wd(1)) go to 100
      if (input(1:4).eq.wd(2)) go to 200
      if (input(1:4).eq.wd(3)) go to 999
      if (input(1:4).eq.wd(4)) go to 150
      if (input(1:4).eq.wd(5)) go to 125
      if (input(1:4).eq.wd(6)) go to 300
      if (input(1:4).eq.wd(7)) go to 400
      go to 1
c
c---- macro 'feap'
c---- initialize variables
 100  newmesh=.true.
c for neighbor finding:
      nmeth=2
      newlst=1
      ngtlst=0
      perub(1:2)=1.e30
      perlb(1:2)=-1.e30
      perlen(1:2)=2e30
c
      DD_set=.false.
      nce=0;ncb=0;NCEMAX=200
      allocate(elist(2,NCEMAX))
c---- read and print control information
      head = input(5:)
      call GlobalSettings
c
      write(6,'(/'' '',a80//
     1 5x,''number of nodal points (max) ='',i6/
     2 5x,''number of elements     (max) ='',i6/
     4 5x,''dimension of coordinate space='',i6/
     5 5x,''degrees of freedom/node      ='',i6/
     6 5x,''nodes/element (maximum)      ='',i6/
     7 5x,''dimension of stress array    ='',i6/
     9 5x,''number of quad pts/element   ='',i6/
     o 5x,''number of local nodes/element='',i6)')
     1 head,maxnp,maxel,ndm,ndf,nen,nsdm,nquad,nad

      write(*,*)

      neq  = maxnp*ndf
      neqad = (maxnp + maxel*nad)*ndf

      allocate(dr(ndf*maxnp))
      allocate(x(nxdm*(maxnp+maxel*nad)))
      allocate(Isrelaxed(maxnp))
      allocate(idtemp(ndf,maxnp))
      allocate(energy(maxnp))
      allocate(f(ndf*maxnp))
      allocate(db(neqad))
      allocate(id(ndf*maxnp))
      allocate(ix(nen1*maxel))
      allocate(itx(3,maxel))
      allocate(b(neqad))

      allocate(b0(ndf,maxnp+maxel*nad))
CC--Jun Song: allocate neighflag for each atom
c     c--DampForce for Marder and langevin Thermostats
      allocate(UpdateNeigh(maxnp))
      allocate(DampForce(ndf,maxnp))

c
c dynamo variables
c
      allocate(dis(3,maxnp),rdyn(maxnp),knbr(maxnp))
      allocate(rold(3,maxnp))
      allocate(nnindx(0:maxnp),nnlst(maxnp,neimax))
cc--JS      nnindx(0)=0
cc--JS      nnindx(1)=0
      nnindx=0

c     Initializing dynamical allocated variables
      IsRelaxed=1
      idtemp=.false.
      energy=0.d0
      dr = 0.d0
      x = 0.d0
      f = 0.d0
      b = 0.d0
      b0 = 0.d0
      db = 0.d0
      id = 0
      ix = 0
      UpdateNeigh=0

c     Initialize crack tip position to 0,0
c     If getc macro is not called then crack tip is assumed to be 0,0
c     if getc is called it will initialize the crack tip position 
      xtip_init = 0.0d0
c---- call mesh input subroutine to read and print all mesh data


      do i=1,nprocs
       if(i.eq.rank+1) then
        call pmesh(id,x,ix,f,b,dr,itx)
        call chkper
        neq = numnp*ndf
        rcutsq=CutoffR2(1)
        if(cutfact.le.1.d0) then
         dradn = 0.1*dsqrt(rcutsq)
        else
         dradn = (CUTFACT-1.d0)*dsqrt(rcutsq)
        endif
        rctsqn = (sqrt(rcutsq) + dradn)**2
        call EchoSettings
        call latticecheck(x)
        call StatusCalc(x,ix,.false.)
        if(nqc.ne.0) call GetDetectionBand(ix,nen1,numel,x,num2Dnode,
     $     nxdm,itx,IsRelaxed)
        call PlotEsi(x,ix)
        call PlotEsi_vtk(x,ix)
        
       endif
!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      enddo
        go to 1
c
c---- macro 'xxxx'
125   go to 1
c
c---- macro 'xxxx'
150   continue
      go to 1
c
c---- macro 'macr'
c---- set up macro program for execution
c---- call macro solution module for establishing solution algorithm
 200  call pmacr(id,x,ix,f,b,dr,db,itx)
      go to 1
c
c---- macro 'xxxx'
300   continue
      go to 1
c
c---- macro 'xxxx'
400   continue
      go to 1
c
c---- macro 'stop'
 999  do l=7,99
         close(l)
      enddo
!      call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
!      call MPI_FINALIZE(ierr)
       stop
      end


