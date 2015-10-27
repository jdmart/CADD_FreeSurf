!*******************************************************************
!**   
!     <**  MODULE mod_global: contains routines which set all relevent
!     global
!     !**                        flags and dimensions at the start of a
!     run.
!     !**
!     !**  Variable Definitions:
!     !**  ---------------------
!**   
!     !**  NOTE: variables marked by (????) are ones that I believe
!     should be
!**   eliminated eventually (R.E. Miller (jul 02))
!**   
!     <**  Integer Variables:
!     <**    nxdm                 -  Dimension of the nodal coord array
!     x()
!     <**    ndm                  -  Number of spatial dimensions
!     <**    ndf                  -  Number of degrees of freedom per
!     node
!     <**    nen                  -  Number of nodes per element
!     <**    nen1                 -  Number of nodes per element plus 1
!     <**    nsdm                 -  Number of stress components
!     <**    numnp                -  current number of nodes in the mesh
!     <**    numel                -  current number of elements in the
!     mesh
!     <**    neq                  -  number of equations (ndf*numnp)
!     <**    nquad         (????) -  number of quadrature points
!     <**    nad           (????) -  number of additional degrees of
!     freedom per element
!     <**    nshpdm        (????) -
!     <**    nstad         (????) -
!     <**    nst0          (????) -
!     <**    nentot        (????) -
!     <**    maxnp                -  SIZE of allocated storage for nodes
!     <**    maxel                -  SIZE of allocated storage for
!     elements
!     <**    nxsj                 -  dimension of shape function array,
!     xsj
!     <**    nshp                 -  number of shape functions, 
!     <**    nstr                 -  dimension of stress and strain
!     arrays, str,eps
!     <**    nst           (????) -
!     <**   CUTFACT:
!     !**   EFFECTIVE CUTOFF FACTOR (in multiples of rcut gives the
!     radius in
!     !**   the deformed configuration contributing to an atom's
!     neighbor lis).
!     !**   The larger this factor the fewer updates to the neighbor
!     lists will
!     !**   be necessary during relaxation, but more atoms will be
!     sampled at
!     !**   each step.
!**   
!     <**   New CADD variables:
!**   
!     <**   nqc       -        number of atoms in the atomistic region
!     <**   nspring   -        number of interface atoms
!     <**   numnpp1   -        if its -1, all reference to the brinell
!**   indenter is ignored.  Otherwise it should be
!**   set to numnp+1, the location where the brinell indenter is stored
!**   in
!**   the atom list
!     <**   z_length  - periodic length in the z direction
!     <**   DD_set    - flag to show that leo's stuff is initialized
!     <**   newmesh   - flag to tell the detection routine to re
!     -initialize
!     <**   b0        - stores the displacements at the start of the
!     time
!     !**               step, before any relaxation starts.
!     <**   INDRAD    - brinell indenter radius
!     <**   INDRADSQ  - INDRAD**2
!     <**   energy    - site energy of each atom
!     <**   IsRelaxed - identifies each atom type:  -1 pad
!     !**                                            0 continuum
!     !**                                            1 atom
!     !**                                            2 interface
!     <**   idtemp    - used to fix (set force to zero) constrained
!     degrees of freedom.
!**   
!**   
!**   
!**   
!**   
!**   Contains Routines:
!**   ------------------
!**   
!**   GlobalSettings -- Reads in global feap settings and sets others to
!**   default values
!**   EchoSettings   -- Echoes all settings.
!**   
!**************************************************************

      module mod_global

!     * Variable Definintions
      integer nxdm,ndm,ndf,nen,nen1,nsdm,numnp,numel,neq,nquad,nad
     $     ,nshpdm ,nstad,nst0,nentot,maxnp,maxel,nxsj,nshp,nstr,nst
     $     ,niter,nqc,nspring
      integer :: numnpp1=-1
      double precision z_length,CONVERGETOL,rnmax
      logical DD_set,newmesh
      double precision, pointer:: b0(:,:)
      double precision :: CUTFACT,INDRAD=-1.d0,INDRADSQ
      double precision, allocatable:: energy(:)
!     double precision, allocatable:: amass(:)
      integer, allocatable:: IsRelaxed(:)
      logical, allocatable :: idtemp(:,:)
      integer i_initial, i_final, indextimeH
      integer :: NumTotH=0 
!     Jun Song Modification
!     NHrescale means how many steps to do H stablizer
!     MaxTemp, do H stablizer if exceed this Temp
!     NVEFlag determines whether to include H in thermostat
!     Number of MD steps doing temp rescaling
      integer NumMDRescale,Debugflag  
      logical UseRescale
      integer NHrescale, HNVEFlag
      integer numperiodz,num2Dnode
      double precision Twindow, MaxHTemp
!     *Qu modification begins related to virial stress
      integer, allocatable :: imaterial(:) !!, boundary(:)
      double precision, allocatable :: virst(:,:,:), rssgb(:)
      double precision, allocatable :: avevirst(:,:,:)
      double precision, allocatable :: rsatomstress(:,:,:)
      double precision, allocatable :: rsatomstressf(:,:,:)
      integer, allocatable :: atomSpecie(:)
!     Jun Song: neighbor flag-determine if updating an atom's neighbor
!     DampForce store the damp forces for Marder&langevin, allocate in
!     pmain.f
      integer, allocatable :: UpdateNeigh(:)
      double precision, allocatable :: DampForce(:,:)
!     Jun Song: total_energyMD is total energy for system
!     energyH is energy for H atom
!     timestep1 is to implement 2 timesteps
!     SimStep is current # of simlation steps
      double precision total_energyMD,energyH, timestep1
      integer dim, SimStep
      double precision perthick
!     Global parameter for Tempset-set in dosteps. Default 0.01
!     Global parameter for NoseHoover, from input
      double precision :: SysTemp = 0.01d0, NHDampCoeff
!     JS-Scale Rationfor Langevin (e.g., Use the coefficient in Marder
!     as base
!     If 0.5, use half the coefficient in Marder)
      double precision LVscaleRatio
!     *Qu modification ends
      character*80 head
      double precision timeol,time,dt, boltzmannConst
      integer indexPad, indexInterface, indexContinuum, indexAtom
      integer MAXNEIGHBORS, iHnumber
      logical initPlot
!     ------Hack -SC
!     This is set in mesh.f to lattice spacing in x. 
!      double precision, parameter :: x_move_mesh = 20.0     
      double precision :: x_move_mesh
      logical :: MoveMesh, Moved
!     Current Crack tip position and other crack tip related parameters
      double precision :: xtip(2),  xtip_old(2), xtip_init(2)
      double precision :: xtip_actual(2), x_tip_dir(2), crack_motion(2)
!     -----End Hack -SC
!     Since the xtip is located in the atomistic region and is always an
!     atom 
!     position in the reference coordinate it will always be a multiple
!     of the 
!     lattice spacing in x
      double precision :: pad_width
!     Adjustable parameter read from input file to control the width of
!     the pad atom region. Needs to be bigger to accomodate more
!     dislocations on the same slip plane entering the atomistic region.
!     Dislocations can be in the pad region provided the template method
!     is used. 

!     As the pad region gets bigger, there is a jump in the dislocation
!     passing algorithm that makes the dislocation move past the pad
!     atoms going either way Atom-> continuum or Continuum-> atom. 

!     Template method not used here .... need implementation 


      contains

!------------------------------------------------------------------
! GlobalSettings -- Reads in global feap settings and sets others to
!                   default values
!
!      Author :
!            R. Miller (01/13/98)
!
!      Revisions :
!            R. Miller (10/10/98)
!            S. Chakravarthy (2012-2014)
!---
      subroutine GlobalSettings

      implicit none
!     
!**   Store FEAP values
!     
      read(5,*) maxnp,maxel,ndm,ndf,nen,nsdm,nquad,nad
      read(5,*) numperiodz
      nxdm=3
      nen1=nen+1
!     
      nsdm = max0(nsdm,1)
      nquad= max0(nquad,1)
      nad  = max0(nad,0)
!     
      nen1 = nen + 1
      nst0 = nen*ndf
      nstad = nad*ndf
      nst  = nst0 + nstad
      nentot = nen + nad
      nstr = nquad*nsdm
      nshpdm = (ndm+1)*(nen+nad)
      nshp = nquad*nshpdm
      nxsj = nquad

      indexContinuum = 0
      indexPad = -1
      indexInterface = 2
      indexAtom = 1
      initPlot = .false.
      
!**   default
      CUTFACT=1.4d0
      MAXNEIGHBORS=100
      
!**   physical constants
      boltzmannConst = 8.629064e-05 ! eV/K
      end subroutine globalsettings
!**   ---------------------------------------------------------------
!**   EchoSetting  :  Echo all global settings of note to std out.
!**   
      subroutine EchoSettings

      implicit none
!**   echo
      write (*,*)
      write(*,*) '*** Global Parameter Settings ***'
      write(*,*)
      write(*,300) 'Number of spatial dimensions: ',ndm
      write(*,300) 'Number of degrees of freedom per node: ',ndf
      write(*,300) 'Number of nodes per element: ',nen
      write(*,300) 'Number of stress components: ',nsdm
      write(*,*)
      write(*,400) 'Effective cutoff factor: ', CUTFACT
      write(*,*)
 300  format(5x,a50,i7)
 400  format(5x,a50,g14.5)
      end subroutine echosettings

      End module mod_global


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!module with mpi variables!!!!!!!!!!!!!!!!!
      module mod_parallel

!     include '/opt/hpmpi/include/mpif.h'
      integer nprocs,rank,ierr

      end module mod_parallel




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!module with timming variables!!!!!!!!!!!!!!!!!
      module mod_timming
!     
!     ct1 - time at call to ma06
!     ct2 - local start time
!     ct3 - local end time
!     ct4 - MD time 
!     ct5 - detection band time
!     ct6 - time in ma06
!     ct7 - time in fem

      real ct1, ct2, ct3, ct4, ct5, ct6, ct7, ct8

      end module mod_timming
