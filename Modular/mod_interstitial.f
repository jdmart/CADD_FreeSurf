!
!
!
!

      Module mod_interstitial

!* Type Defintions
      type interstitial_atoms
      integer          natoms, maxAtoms
      double precision, dimension(:,:), pointer :: r
      integer, dimension(:), pointer :: element
      end type interstitial_atoms

      Contains

!- Initializes the collection of interstitials
      subroutine Init(interstitial)
      implicit none

      type(interstitial_atoms) interstitial
      integer MAXATOMS
      
      MAXATOMS=100

      allocate(interstitial%r(3,MAXATOMS))
      allocate(interstitial%element(MAXATOMS))
      interstitial%maxAtoms = MAXATOMS

      end subroutine Init


!     - Reads in interstitial data
      subroutine ReadInterstitialData(interstitial)
      implicit none

      type(interstitial_atoms) interstitial
      character (len=80) :: input_file
      double precision x, y, z
      integer ielement, iAtoms, n

      input_file = 'interstitial.dat'
      open(unit=10, file = input_file, status = 'old')

      n = 0
      do iAtoms	= 1, interstitial%maxAtoms
	read(10,*,end=10) x, y, z, ielement
	n = n+1
	interstitial%r(1,n) = x	
	interstitial%r(2,n) = y	
	interstitial%r(3,n) = z	
	interstitial%element(n) = ielement	
      end do
10	continue
close(10)

        interstitial%natoms = n

        end subroutine ReadInterstitialData


      end module mod_interstitial

