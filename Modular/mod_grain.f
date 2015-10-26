!****************************************************************
!**   
!**   MODULE mod_grain : contains definition of model grains and related
!**   routines.
!**   
!**   
!**   Variable Definitions:
!**   ---------------------
!**   
!**   integer variables:
!**   ngrains          - number of grains
!**   type(graintype) variables:
!**   grains(ngrains) -  All data relevant to the grain structure
!**   ncell -  number of atoms in the orthogonal repeating cell
!**   cell(3,ncell) -  coordinates of repeating cell atoms
!**   dcell(3) -  length of periodic cell in xyz directions
!**   
!**   matgrain -  material in the grain
!**   xlatvect(3,3) -  lattice vectors defining natural coordinates of
!**   the grain
!**   rotation -  rotation (degrees) of the natural coordinates defining
!**   the grain
!**   rotmat -  rotation matrix elements
!**   numvrts -  number of vertices defining the grain
!**   grainarr -  polygon defining grain
!**   rotgrain -  polygon defining grain rotated to grain's natural
!**   coordinate system
!**   refatom -  reference atom coordinates
!**   fudgevec -  fudge vector
!**   rfudgvec -  fudge vector rotated to grain's natural coordinate
!**   system
!**   
!**   Contains Routines:
!**   ------------------
!**   
!**   ReadGrainData        - Reads in the number of grains and their
!**   properties
!**   OutputGrainData      - Echoes all grain data, either read in or
!**   computed
!**   ProcessGrains        - Initializes grain data based on that which
!**   is read in
!**   ProcessLatVect       - given 2 lattice vectors, computes third,
!**   and computes rotmat
!**   ProcessGrainGeometry - Computes rotgrain and rfudgvec
!**   
!***************************************************************

      Module mod_grain

!     * Type Defintions
      type graintype
!     stuff for repeating cell
      integer ncell
      double precision, dimension (:,:), pointer :: cell
      double precision, dimension (3)                :: dcell
!     stuff for grain structure and energy
      double precision rotation,xlatvect(3,3),rotmat(2),refatom(3)
      integer matgrain,numvrts
      double precision, dimension (:), pointer   ::  fudgevec,rfudgvec
      double precision, dimension (:,:), pointer ::  grainarr,rotgrain
      end type graintype

!     * Variable Definintions
      integer ngrains
      type(graintype), dimension(:), pointer :: grains

      contains
!------------------------------------------------------------------
!     ReadGrainData -- Read in th grain data from .geo file
!     
!     Passed Parameters :
!     geomfile  (in):  name of geometry file
!     
!     Module Parameters :
!     ngrains         (out) :defined above
!     grains(ngrains) (out)
!     
!     Algorithm :
!     obvious
!     
!     Notes :
!     may not be general for 3D
!     
!     Author :
!     R. Miller (01/17/98)
!     
!     Revisions :
!     none
!     
!--   
      subroutine ReadGrainData(geomfile,key)
      use mod_global
      use mod_file
      use mod_material
      implicit none

!--   Variables transferred
      character(len=80) :: geomfile
      character(len=4) ::key
!--   Local variables
      integer i, j, k, iunit
!     
!--   Check that the materials are ready
      if(nmaterials.lt.1) then
         write(*,*) '***ERROR: Material definitions must 
     $come before grain definitions'
         stop
      endif
!--   Open the specified file
      if(key.eq.'dire') then
         iunit=5
      else
         if(FileExists(geomfile,.false.)) then
            call iofile(geomfile,'formatted  ',iunit,.true.)
         else
            write(*,*) '**ERROR: Grain geometry file not found.**'
            write(*,*) '         Filename:',geomfile
            stop
         endif
      endif

!--   Readin data
      read(iunit,*)ngrains
      allocate(grains(ngrains))
      do i=1,ngrains
         allocate(grains(i)%fudgevec(ndm))
         allocate(grains(i)%rfudgvec(ndm))
      enddo

      do i = 1, ngrains
         read(iunit,*) grains(i)%matgrain
         if(grains(i)%matgrain.lt.1.or.grains(i)%matgrain.gt.nmaterials)
     $        then
            write(*,*) '***ERROR: An undefined material type '
     $           ,grains(i)%matgrain
            write(*,*) '          has been assigned to grain number ',i
            stop
         endif
         read(iunit,*) (grains(i)%refatom(j),j=1,nxdm)
         read(iunit,*) (grains(i)%fudgevec(j),j=1,ndm)
         read(iunit,*) ((grains(i)%xlatvect(j,k),k=1,3),j=1,2)
         read(iunit,*) grains(i)%rotation
         read(iunit,*) grains(i)%numvrts

         allocate(grains(i)%grainarr(ndm,grains(i)%numvrts))
         allocate(grains(i)%rotgrain(ndm,grains(i)%numvrts))

         read(iunit,*)((grains(i)%grainarr(k,j),k=1,ndm),
     $    j=1,grains(i)%numvrts)
      enddo
!--   Close file
      if(iunit.ne.5) close(iunit)

      end subroutine ReadGrainData
!------------------------------------------------------------------
! OutputGrainData : Print out all the grain data
!
!      Passed Parameters :
!            none
!
!      Module Parameters :
!            ngrains         (out) :defined above
!            grains(ngrains) (out)
!
!      Algorithm :
!            obvious
!
!      Notes :
!
!
!      Author :
!            R. Miller (01/17/98)
!
!      Revisions :
!              none
!
!--
      subroutine OutputGrainData()
      use mod_global
      implicit none
!--   Local variables
      integer i, j, k
      write(*,'(a)') 'GRAIN INFORMATION'
      write(6,9000)ngrains
      do i = 1, ngrains
         write(6,9001) i
         write(6,9002) grains(i)%matgrain
         write(6,9007)(grains(i)%refatom(j),j=1,nxdm)
         write(6,9003)((grains(i)%xlatvect(j,k), k =1,3),j=1,3)
         write(6,9004) grains(i)%rotation
         write(6,9005) grains(i)%numvrts
         write(6,9006)((grains(i)%grainarr(k,j),k=1,ndm),j=1
     $        ,grains(i)%numvrts)
         write(6,9008)
         write(6,9009) (grains(i)%dcell(j),j=1,3)
         write(6,9010) grains(i)%ncell
         write(6,9011)((grains(i)%cell(k,j),k=1,3),j=1,grains(i)%ncell)
      enddo

!--   Format statements
 9000 format('    Number of Grains   = ',i3)
 9001 format('    ================================================
     $ == Grain ', i4)
 9002 format('        Material ',i4)
 9003 format('        Lattice Vectors',3(/8x,3f10.5))
 9004 format('        Rotation (degrees) ', f10.5)
 9005 format('        Number of vertices = ', i4)
 9006 format(8x,2e15.5)
 9007 format('        Reference Atom at '/8x,3f10.5)
 9008 format('        Cell Structure ')
 9009 format('         Cell dimensions',/8x,3f10.5)
 9010 format('         Cell Bravais Lattice sites :', i4)
 9011 format(8x,3e15.5)

      return
      end subroutine OutputGrainData
!------------------------------------------------------------------
! ProcessGrains -- Given the grain data that was read in, compute other
!                  grain related data.
!
!      Passed Parameters :
!            none
!
!      Module Parameters :
!            none
!
!      Algorithm :
!            obvious
!
!      Notes :
!          GrainCry call should soon be obsolete
!
!      Author :
!            R. Miller (01/17/98)
!
!      Revisions :
!              none
!
!--
      subroutine ProcessGrains
      implicit none
      integer i
!--   Process Lattice Vectors
      call ProcessLatVect()

!--   Check grain geometry
!--   Can we relax the requirement of convex grains?  I think so.
      call ProcessGrainGeometry()

!--   Make Cell Structure and representative crystallite.
      do i = 1, ngrains
         call GetCellData(grains(i))
      end do
      end subroutine ProcessGrains

!------------------------------------------------------------------
! ProcessLatVect : Form three unit lattice vectors
!                  for each grain.
!
!      Passed Parameters :
!            none
!
!      Module Parameters :
!            grains%xlatvect  (in/out)
!            grains%rotmat    (out)
!
!      Algorithm :
!        For each grain do
!           1)  Find third lattice vector by cross product
!           2)  Normalise the three vectors
!           3)  Check if determinent of resulting matrix is unity
!           4)  Make cosine and sine of rotation angle
!        end do
!
!      Notes :
!            the only variables in the grain data that are changed
!            are xlatvect and rotmat
!
!      Author :
!            R. Miller (01/17/98)
!
!      Revisions :
!              none
!
!--
      subroutine ProcessLatVect()
      implicit none

!--   Local Variables
      integer i, j, k
      double precision detq, qmag, q(3,3)
      double precision, parameter ::eps = 1.D-09, Pi=3.1415926535898


!--   Loop over all grains
      do i = 1, ngrains

!Load in lattice vectors to local arrays
         do j = 1,2
            do k = 1,3
               q(j,k) = grains(i)%xlatvect(j,k)
            end do
         end do

!Find third vector by cross product
         q(3,1)=q(1,2)*q(2,3)-q(1,3)*q(2,2)
         q(3,2)=q(1,3)*q(2,1)-q(1,1)*q(2,3)
         q(3,3)=q(1,1)*q(2,2)-q(1,2)*q(2,1)

!Normalise the  vectors
         do j=1,3
            qmag=dsqrt(1.d0/(q(j,1)*q(j,1)+q(j,2)*q(j,2)+q(j,3)*q(j,3)))
            do k=1,3
               q(j,k)=q(j,k)*qmag
            enddo
         enddo


!Check if determinant is unity
         detq = q(1,1)*(q(2,2)*q(3,3)-q(2,3)*q(3,2))-   
     $    q(1,2)*(q(2,1)*q(3,3)-q(2,3)*q(3,1))+   
     $    q(1,3)*(q(2,1)*q(3,2)-q(2,2)*q(3,1))
         if (dabs(detq-1.d0).gt.eps) then
            print *
            print *,'*** ERROR: Lattice Vector Error in grain ', i
            print *,'Determinant = ',detq
            print *
            stop
         endif

!Put back the lattice vector into grainstrc
         do j = 1,3
            do k = 1,3
               grains(i)%xlatvect(j,k) = q(j,k)
            end do
         end do

!Make cosine and sine of the rotation angle
         grains(i)%rotmat(1) = cos(Pi*grains(i)%rotation/180.0)
         grains(i)%rotmat(2) = sin(Pi*grains(i)%rotation/180.0)
      end do
      return
      end subroutine ProcessLatVect
!------------------------------------------------------------------
! ProcessGrainGeometry -- perform rotations on grain data
!
!      Passed Parameters :
!            none
!
!      Module Parameters :
!            ngrains        (in)
!            grains%refatom (in)
!            grains%numvrts (in)
!            grains%rotmat  (in)
!            grains%grainarr(in)
!            grains%rotgrain(out)
!            grains%fudgevec(in)
!            grains%rfudgvec(out)
!
!      Algorithm :
!       Algorithm :
!          1) check if reference atoms are in the interior of grain
!          1) compute rotgrain and rfudgvec
!
!      Notes :
!            notes
!
!      Author :
!            R. Miller (01/17/98)
!
!      Revisions :
!              none
!
!--
      subroutine ProcessGrainGeometry
      use mod_global
      implicit none
!--   Local variables
      integer i, j, k, nvt
      logical PointInGrain
      double precision u(ndm)

      do i = 1, ngrains
         if(.not.PointInGrain(grains(i)%refatom,i)) then
            write(*,*) '***ERROR: The reference atom of grain',i
            write(*,*) '          is outside of the grain polygon.'
            stop
         endif
         nvt = grains(i)%numvrts
         do j = 1, nvt
            do k = 1, ndm
               u(k) = grains(i)%grainarr(k,j) - grains(i)%refatom(k)
            end do
            grains(i)%rotgrain(1,j) =  u(1)*grains(i)%rotmat(1) + u(2)
     $           *grains(i)%rotmat(2)
            grains(i)%rotgrain(2,j) = -u(1)*grains(i)%rotmat(2) + u(2)
     $           *grains(i)%rotmat(1)
         end do
         grains(i)%rfudgvec(1) =   grains(i)%fudgevec(1)
     $        *grains(i)%rotmat(1) + 
     $    grains(i)%fudgevec(2)*grains(i)%rotmat(2)
         grains(i)%rfudgvec(2) =  -grains(i)%fudgevec(1)
     $        *grains(i)%rotmat(2) + grains(i)%fudgevec(2)
     $        *grains(i)%rotmat(1)
      enddo
      end subroutine ProcessGrainGeometry
!************************************************************************
      subroutine GetCellData(grain)
      use mod_global
      use mod_material
      use mod_cluster
      use mod_dynamo
      implicit none

!** Transferred Variables **!
      type(graintype) grain

!** Parameter Declarations **!
      double precision,parameter::  tol=1.d-6

!** Local Variables **!
      double precision tmp(3),brot(3,3),rclust,rtol,cx,cy,cz
      double precision, pointer ::cell2(:,:)
      integer icry,i,j,k,maxcell,nlast,nsort
      logical xzero,yzero,zzero,break
      type(bravaismat) matl
      type(cluster) clust
!     c
!     C     Analyze rotated crystal structure
!     c
!     c Find the repeat cell of bravais sites in the natural coord
!     system of the grain
!     c by building a temporary cluster of bravais sites and analysing
!     it.  Note that
!     c if the assumed size of the cluster is too small, it increases it
!     and tries
!     c again.
!     c
!     ctry something simple: find atoms closest to origin along
!     c each axis and take this as the unit cell.

! Rotate BL vectors into slip c.s. (note sexy f90 function call)
      brot=matmul(grain%xlatvect,material(grain%matgrain)%bvec)
      rclust=sqrt(dot_product(brot(1:3,1),brot(1:3,1)))
      rtol=tol*rclust
      rclust=rclust*5.d0
      maxcell=100
!     c
 1    continue
      allocate(cell2(3,maxcell))
      matl%bvec=brot
      matl%nbasis=1
      allocate(matl%basis(3,1))
      matl%basis(1:3,1)=0.0d0
      allocate(matl%ispec(1))
      matl%ispec(1)=1
      matl%volume = material(grain%matgrain)%volume

      call buildcluster(matl,rclust,clust)

      deallocate(matl%basis)
      deallocate(matl%ispec)
      cx=1.d30
      cy=1.d30
      cz=1.d30
      do i=1,clust%natoms
         if ((clust%x(1,i).lt.-rtol).or.(clust%x(2,i).lt.-rtol) 
     $    .or.(clust%x(3,i).lt.-rtol)) go to 10
         xzero=dabs(clust%x(1,i)).lt.rtol
         yzero=dabs(clust%x(2,i)).lt.rtol
         zzero=dabs(clust%x(3,i)).lt.rtol
         if (xzero.and.yzero.and.(.not.zzero)) then
            if (clust%x(3,i).lt.cz) cz=clust%x(3,i)
         else if (xzero.and.zzero.and.(.not.yzero)) then
            if (clust%x(2,i).lt.cy) cy=clust%x(2,i)
         else if (yzero.and.zzero.and.(.not.xzero)) then
            if (clust%x(1,i).lt.cx) cx=clust%x(1,i)
         endif
 10   enddo
      if (cx.gt.1.d29) then
         rclust=2*rclust
         print *,'***WARNING: Repeating pattern not located in x-dir''n'
         print *,'          Increased computed cluster size:',rclust
         go to 1
      else if (cy.gt.1.d29) then
         rclust=2*rclust
         print *,'***WARNING: Repeating pattern not located in y-dir''n'
         print *,'          Increased computed cluster size:',rclust
         go to 1
      else if (cz.gt.1.d29) then
         rclust=2*rclust
         print *,'***WARNING: Repeating pattern not located in z-dir''n'
         print *,'          Increased computed cluster size:',rclust
         go to 1
      endif
      grain%ncell=0
      do i=1,clust%natoms
         xzero=(clust%x(1,i).lt.cx-rtol).and.(clust%x(1,i).gt.-rtol)
         yzero=(clust%x(2,i).lt.cy-rtol).and.(clust%x(2,i).gt.-rtol)
         zzero=(clust%x(3,i).lt.cz-rtol).and.(clust%x(3,i).gt.-rtol)
         if (xzero.and.yzero.and.zzero) then
            grain%ncell=grain%ncell+1
            if (grain%ncell.gt.maxcell) then
               maxcell=maxcell*2
               print *,'***WARNING: Insufficient storage in cell matrix'
               print *,'            increased maxcell to:',maxcell
               go to 1
            endif
            do j=1,3
               cell2(j,grain%ncell)=clust%x(j,i)
            enddo
         endif
      enddo
      call deallocate_cluster(clust)
      allocate(grain%cell(3,grain%ncell))
      grain%cell(1:3,1:grain%ncell)=cell2(1:3,1:grain%ncell)
!     c
!     c sort into ascending y-layers.  In each layer sort into ascending
!     c x coord.  If more than one with same x coord, sort into
!     ascending z
!     c coord.
!     c
!     c  First sort, y coord as key:
!     c
      call sort(grain%cell,3,grain%ncell,2,grain%ncell)
!     c
!     c  Second sort, within each y-layer sort by x coord
!     c
      nlast=0
      do i=1,grain%ncell
         break=(i.eq.grain%ncell)
         if (.not.break) break=(abs(grain%cell(2,i+1)-grain%cell(2
     $        ,i)).gt.rtol)
         if (break) then
            nsort=i-nlast
            do j=1,nsort
               do k=1,3
                  cell2(k,j)=grain%cell(k,nlast+j)
               enddo
            enddo
            call sort(cell2,3,maxcell,1,nsort)
            do j=1,nsort
               do k=1,3
                  grain%cell(k,nlast+j)=cell2(k,j)
               enddo
            enddo
            nlast=i
         endif
      enddo
!     c
!     c third sort, within each line of atoms with the same x and y,
!     c sort by z-coord.
!     c
      nlast=0
      do i=1,grain%ncell
         break=(i.eq.grain%ncell)
         if (.not.break) break=(abs(grain%cell(2,i+1)-grain%cell(2
     $        ,i)).gt.rtol) 
     $    .or.(abs(grain%cell(1,i+1)-grain%cell(1,i)).gt.rtol)
         if (break) then
            nsort=i-nlast
            do j=1,nsort
               do k=1,3
                  cell2(k,j)=grain%cell(k,nlast+j)
               enddo
            enddo
            call sort(cell2,3,maxcell,3,nsort)
            do j=1,nsort
               do k=1,3
                  grain%cell(k,nlast+j)=cell2(k,j)
               enddo
            enddo
            nlast=i
         endif
      enddo
      deallocate(cell2)
!     c
!     c at this point, cell is sorted.
!     c
!     c     Compute projected bravais site area in the global xy plane
      grain%dcell(1)=cx
      grain%dcell(2)=cy
      grain%dcell(3)=cz
      z_length=numperiodz*cz
      perlen(3)=numperiodz*cz
      perub(3)=numperiodz*cz
      perlb(3)=0.d0
      end subroutine getcelldata

      SUBROUTINE SORT(RA,MRA,NRA,M,N)
!     C     Based on a Numerical Recipes routine.
      IMPLICIT NONE
      integer, intent(in) :: nra,mra,n,m
      double precision,intent(inout)::  RA(MRA,NRA)

      integer, PARAMETER ::MAXM=10
      double precision RRA(MAXM)
      integer ir,i,j,k,L
      IF (N.LE.1) RETURN
      IF (MRA.GT.MAXM) THEN
         PRINT *,'***ERROR: MRA too large in SORT'
         STOP
      ENDIF
      L=N/2+1
      IR=N
 10   CONTINUE
      IF(L.GT.1)THEN
         L=L-1
         DO K=1,MRA
            RRA(K)=RA(K,L)
         ENDDO
      ELSE
         DO K=1,MRA
            RRA(K)=RA(K,IR)
         ENDDO
         DO K=1,MRA
            RA(K,IR)=RA(K,1)
         ENDDO
         IR=IR-1
         IF(IR.EQ.1)THEN
            DO K=1,MRA
               RA(K,1)=RRA(K)
            ENDDO
            RETURN
         ENDIF
      ENDIF
      I=L
      J=L+L
 20   IF(J.LE.IR)THEN
         IF(J.LT.IR)THEN
            IF(RA(M,J).LT.RA(M,J+1))J=J+1
         ENDIF
         IF(RRA(M).LT.RA(M,J))THEN
            DO K=1,MRA
               RA(K,I)=RA(K,J)
            ENDDO
            I=J
            J=J+J
         ELSE
            J=IR+1
         ENDIF
         GO TO 20
      ENDIF
      DO K=1,MRA
         RA(K,I)=RRA(K)
      ENDDO
      GO TO 10
      END subroutine sort


      end module mod_grain


