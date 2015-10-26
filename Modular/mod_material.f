!*************************************************************
!**   
!**   MODULE mod_material : contains definition of model materials and
!**   related
!**   routines.
!**   
!**   
!**   Variable Definitions:
!**   ---------------------
!**   
!**   Integer Variables:
!**   nmaterials           -  Number of materials
!**   nbasismax            -  Maximum number of basis atoms
!**   
!**   Type(bravaismat) Variables:
!**   material(nmaterials) -  All data relevant to the material
!**   structure:
!**   %name -  Materials name (character string)
!**   %bvec(:,i) -  Components of Bravais vector i
!**   %nbasis -  Number of atoms per Bravais site
!**   %basis(:,i) -  Coordinates of basis atom i
!**   %ispec(nbasis) -  Atomic species of basis atoms (use periodic
!**   table numbers)
!**   %bmag(:) -  Magnitude of Bravais vectors
!**   %volume -  Bravais cell volume
!**   %structure -  3 chr string identifying the structure
!**   (fcc, hcp, hex, etc)cubic lattice constant
!**   %a0     -  cubic lattice constant
!**   %cc     -  cubic elastic constants
!**   
!**   Contains Routines:
!**   ------------------
!**   
!**   ReadMaterials   - Reads in the number of materials and their
!**   properties
!**   OutputMaterials - Prints out information on the materials read in
!**   
!**********************************************

      Module mod_material

!     * Type Defintions
      type bravaismat
      character(len=20) :: name
      double precision bvec(3,3)
      integer nbasis
      double precision, dimension(:,:), pointer :: basis
      integer, dimension(:), pointer :: ispec
      double precision bmag(3)
      double precision volume,a0,cc(6,6)
      character*3 structure
      end type bravaismat

!     * Variable Definintions
      integer nmaterials,nbasismax
      type(bravaismat), dimension(:), pointer :: material

       contains

!------------------------------------------------------------------
! ReadMaterials -- Read in relevant data for all model materials
!                  and store in the structured array material
!
!      Passed Parameters :
!                      none
!
!      Module Parameters :
!          nmaterials (out) : number of materials
!          nbasismax  (out) : maximum number of basis atoms
!         material(:) (out) : material data (see defn above)
!
!      Algorithm :
!            Read in number of materials and then for each material
!            read in all relevant data as defined above.
!
!      Notes :
!            Does not currently allow for non-stoichiometric structures
!
!      Author :
!            E.B.Tadmor (12/31/97)
!
!      Revisions :
!
!--
      subroutine ReadMaterials(ndf,key)

      use mod_file
      implicit none
      integer ndf
      character*4 key

!** Local Variables **!
      integer i,j,iunit,icc,jcc,ncc
      double precision det33,ctmp
      logical error

      character(len=80) datafile

      read *,nmaterials
      print*, 'number of materials', nmaterials
      if (nmaterials.lt.1) then
         print *,'***ERROR: Illegal number of materials specified.'
         stop
      endif
      allocate(material(nmaterials))
      nbasismax=0
      do i=1,nmaterials
         if(key.eq.'dire') then
            iunit=5
         else
            read (*,'(a)') datafile
            if (.not.FileExists(datafile,.true.)) stop
            call iofile(datafile,'formatted  ',iunit,.true.)
         endif
         read(iunit,*) material(i)%name
         read(iunit,*) (material(i)%bvec(:,j),j=1,3)
         read(iunit,*) material(i)%nbasis
         allocate(material(i)%basis(3,material(i)%nbasis))
         material(i)%basis = 0.
         allocate(material(i)%ispec(material(i)%nbasis))
         material(i)%ispec = 0
         read(iunit,*) (material(i)%basis(:,j),material(i)%ispec(j),j=1
     $        ,material(i)%nbasis)
         read(iunit,'(a3)') material(i)%structure
         read(iunit,*) material(i)%a0
         read(iunit,*) ncc
         material(i)%cc=0.d0
         do j=1,ncc
            read(iunit,*) icc,jcc,ctmp
            material(i)%cc(icc,jcc)=ctmp
            material(i)%cc(jcc,icc)=ctmp
         enddo
         do j=1,3
            material(i)%bmag(j) = dsqrt(dot_product(material(i)%bvec(:
     $           ,j),material(i)%bvec(:,j)))
         enddo
         material(i)%volume = det33(material(i)%bvec)
         if(material(i)%volume.le.0.d0) then
            write(*,*) '** ERROR: unit cell volume is negative
     $       or zero for material',i
            write(*,*) 'Make sure Bravais lattice is right-handed'
            stop
         endif
         if (nbasismax.lt.material(i)%nbasis) nbasismax =
     $        material(i)%nbasis
         if(key.ne.'dire') close(iunit)
      enddo

      error=.false.
      if(ndf.lt.3*nbasismax) then
         write(*,*) '***ERROR: ndf is too small to accomodate the'
         write(*,*) '          defined materials.'
         write(*,*) '          Increase ndf to:',3*nbasismax
         error=.true.
      endif
      if(error) stop
      return
      end subroutine ReadMaterials

   !------------------------------------------------------------------
   ! OutputMaterials -- Prints information on loaded materials in an
   !                    orderly fashion
   !
   !      Passed Parameters :
   !                      none
   !
   !      Module Parameters :
   !          nmaterials (out) : number of materials
   !         material(:) (out) : material data (see defn above)
   !
   !      Algorithm :
   !            self explanatory
   !
   !      Notes :
   !          none
   !
   !      Author :
   !            E.B.Tadmor (12/31/97)
   !
   !      Revisions :
   !              none
   !
   !--
      subroutine OutputMaterials()

      implicit none

      !** Local Variables **!
      integer i,j

      print *
      print '(a)','MATERIAL INFORMATION'
      do i=1,nmaterials
         print *
         print '("Mat #",i2," : ",A)',i,material(i)%name
         print '(a)','Bravais lattice vectors:'
         print '("a",i1," = ",3f10.5)',(j,material(i)%bvec(:,j),j=1,3)
         print '(a)','Basis atoms coordinates and species:'
         print '("Atom #",i2," : ",3f10.5,5x,i3)',(j,material(i)%basis(:
     $        ,j),material(i)%ispec(j),j=1,material(i)%nbasis)
      enddo
      print *

      end subroutine OutputMaterials

      End module mod_material
