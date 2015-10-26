!****************************************************************
!     Defines all the basic variables in pmain to use this module
!     These variables are 
!     x - reference configuration of atoms/nodes
!     b - atom/node displacements
!     f - boundary conditions (id=1, displacement, id=0, force)
!     id - b.c. flag (0=free, 1=fixed)
!     ix - ix(1:3,i) nodes making up element i. 
!          ix(4,i): -1,-2,-3: element is in the detection band, 
!                             abs(ix(4,i)) is the "entry side"
!                   1:        element is in the atomistic region
!                   0:        element is in the continuum region
!     itx - adjacency matrix: itx(j,i) is the number of the element
!           adjacent to side j of element i.
!     dr  - forces
!     db  - displament increment/step.
!

      module mod_main
      double precision, pointer:: dr(:),x(:),f(:),b(:),db(:)
      integer, pointer:: id(:),ix(:),itx(:,:)
      end module mod_main
!****************************************************************
