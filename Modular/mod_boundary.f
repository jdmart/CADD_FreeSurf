      module mod_boundary
! elist is the model boundary, see explanation in CONTRI.F
! dbpoly contains the polygons where you want to put one or more
! detection bands.  Just allocate the necessary storage, define the
! polygon (CCW) and the code will find the elements necessary to make
! the band.  If the band should go right up to the free surface (like
! when the atomistic region is on the surface) define the polygon to
! extend out into free space and then come back into the mesh at the
! right place where the detection band should start again.  
!
!  ndbpoly: number of detection band polygons.  Usually 1 for each
!           atomistic region
!  ndbvtx(ndbpoly): number of vertices for each det band polygons
!  dbpoly(2,ndbpoly):  coordinates of the vertices defining the det 
!                      band polygons
!  eldb(:,ndbpoly):  element global # of the detection band element 
!     on idb
!  nelidb(ndbpoly):  # of elements on idb ring
!  dbbound(ndb
      integer nce,ncb,NCEMAX,ndbpoly
      integer, allocatable :: ndbvtx(:)
      integer, allocatable :: eldb(:,:)
      integer, allocatable :: nelidb(:)
      integer, allocatable::elidb(:)
      integer, pointer:: elist(:,:)
      double precision, allocatable :: dbpoly(:,:,:)
      logical, allocatable :: dbbound(:,:), dbboundnear(:)
      integer,parameter :: dbNmax=300

      contains
!**
!**---------------------------------------------------------------
!**  IncreaseElist : Increase storage for constrained edges
!**
      subroutine IncreaseElist(nadd)
      implicit none
      integer nadd
      integer, pointer:: e(:,:)
      allocate(e(2,NCEMAX))
      e(1:2,1:NCEMAX)=elist(1:2,1:NCEMAX)
      deallocate(elist)
      allocate(elist(2,NCEMAX+nadd))
      elist(1:2,1:NCEMAX)=e(1:2,1:NCEMAX)
      NCEMAX=NCEMAX+nadd
      deallocate(e)
      end subroutine IncreaseElist

      subroutine allocate_db
      implicit none
      	allocate(ndbvtx(ndbpoly)) 
	allocate(dbpoly(2,4,ndbpoly))
        allocate(dbboundnear(ndbpoly))
        allocate(nelidb(ndbpoly))
        allocate(eldb(dbNmax, ndbpoly))
        allocate(dbbound(dbNmax, ndbpoly))

        dbboundnear(1:ndbpoly) = .false. 
	ndbvtx(1:ndbpoly)=4
        dbbound = .false. 
      end subroutine allocate_db

      end module mod_boundary
