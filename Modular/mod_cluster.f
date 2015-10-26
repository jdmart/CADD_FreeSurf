!******************************************************************
!**   
!**   MODULE mod_cluster: routines pertaining to building, storing, and
!**   manipulating clusters.
!**   
!**   Variable TYPES Defined:
!**   ---------------------
!**   
!**   type(cluster) variables:
!**   %natoms          - number of atoms in the cluster
!**   %x(3,natoms)     - coordinates of the atoms
!**   %spec(natoms)    - species of each atom
!**   %lattice(natoms) - lattice from which the atom is derived (for
!**   complex lattices)
!**   
!**   Contains Routines:
!**   ------------------
!**   
!**   BuildCluster         - Builds a cluster of given radius from a
!**   given material
!**   deallocate_cluster   - Deallocates a type(cluster) variable
!**   SortBravais          - pre-processing for buildcluster
!**   
!********************************************************************

      Module mod_cluster

!     * Type Defintions
      type cluster
      integer                                      natoms
      double precision, dimension(:,:), pointer :: x,rlat
      integer, dimension(:),pointer             :: spec,lattice
      end type cluster

      contains

!------------------------------------------------------------------
! BuildCluster -- Build a cluster of atoms
!
!      Passed Parameters :
!            mat  (in) : type(bravaismat), contains the material
!                        from which to build the cluster
!            rcut (in) : radius of cluste
!            atoms(out): type(cluster), contains the built cluster
!
!      Algorithm :
!         Span the directions of the three Bravais vectors,
!         but try to be clever about minimizing the range of
!         each loop (important for non-orthogonal bravais lattice
!         vectors).
!
!         Add to the cluster all atoms associated with bravais
!         sites less than rcut from the origin, making the origin
!         site atoms first in the list.
!
!      Notes :
!         Needs to be rigourously tested for speed and optimized
!         if possible - this will be called for every call to
!     "cauchyborn"
!
!      Author :
!            R. Miller (01/17/98)
!
!      Revisions :
!              none
!
!--
      subroutine BuildCluster(mat,rcut,atoms)
      use mod_material
      implicit none

!     variables that are passed
      type(bravaismat) mat
      double precision rcut
      type(cluster) atoms

!--   Local Variables
!     
      type(cluster) at
      integer ii(3),i1,i2,i3,is1,is2,is3,if1,if2,if3,ic1,ic2,ic3,nguess
     $     ,i,j,ia,ib,ib2
      double precision, parameter :: PI43=4.18879d0
      double precision latlen(3),org(2),per3(2),bperp,rlen2,r(3),rb(3)
     $     ,per2,cperp,rcut2,vec(3),rtol
      logical notzero1,notzero2,notzero3
!     

      rcut2  = rcut*rcut
!     
!     estimate cluster size using (volume of sphere)/(volume per atom),
!     with
!     a safty factor of 2
!     
      nguess=2*nint(mat%nbasis*(PI43*rcut*rcut2)/(mat%volume))
      allocate (at%x(3,nguess))
      allocate (at%rlat(mat%nbasis,nguess))
      allocate (at%spec(nguess))
      allocate (at%lattice(nguess))
!     C@
      at%x = 0.
      at%rlat = 0.
      at%spec = 0
      at%lattice = 0
!     
!     compute lengths of bravais lattice vectors
!     
 1    do j=1,3
         latlen(j)=dot_product(mat%bvec(1:3,j),mat%bvec(1:3,j))
         latlen(j)=sqrt(latlen(j))
      enddo
!     
!     to optimize speed, sort (using sorting map ii) the bvec.
!     Call the bvec vectors "a" "b" and "c".  "a" is chosen to be
!     the shortest of the three vectors.  Of the remaining 2, "b" is
!     the one that is "least perpendicular" to "a".  Finally, "c" is the
!     remaining vector.  map these using ii(1)=a, ii(2)=b, ii(3)=c.
!     
!     bperp is the length of b perp. to a.  Similarly cperp.
!     "per3" is the projection of c in the plane defined by a and b
!     "per2" is the projection of b along a.
!     
      call SortBravais(mat%bvec,latlen,ii,bperp,cperp,per3,per2)
!     
!     Put in the center site atoms (bravais lattice (0,0,0))
!     now, so that they appear first in list.
!     
      at%natoms=0
      do ib=1,mat%nbasis
         do j=1,3
            rb(j)= mat%basis(ii(1),ib)*mat%bvec(j,ii(1))    
     $           + mat%basis(ii(2),ib)*mat%bvec(j,ii(2))    
     $           + mat%basis(ii(3),ib)*mat%bvec(j,ii(3))
         enddo
         at%natoms=at%natoms+1
         at%x(1:3,at%natoms)=rb
         at%spec(at%natoms)=mat%ispec(ib)
         at%lattice(at%natoms)=ib
      enddo
!     update the rlat
      do ia=1,at%natoms
         do ib=1,at%natoms
            if(ia.eq.ib) then
               at%rlat(ib,ia)=0
            else
               vec=at%x(1:3,ia)-at%x(1:3,ib)
               at%rlat(ib,ia)=sqrt(Dot_Product(vec,vec))
            endif
         enddo
      enddo

!     
!     decide how big to make the loop over the "c" direction
!     
      if3=int(rcut/cperp)
      is3=int(-rcut/cperp)
      if (if3.gt.is3) then
         ic3=1
      else
         ic3=-1
      endif
!     
!     loop in the c direction
!     
      do i3=is3,if3,ic3
         notzero3=i3.ne.0
         org(2)=i3*per3(2)
!     
!     decide how big to make the loop over the "b" direction
!     
         if2=int((rcut-org(2))/bperp)
         is2=int((-rcut-org(2))/bperp)
         if (if2.gt.is2) then
            ic2=1
         else
            ic2=-1
         endif
!     
!     loop in the b direction
!     
         do i2=is2,if2,ic2
            notzero2=i2.ne.0
            org(1)=i3*per3(1)+i2*per2
!     
!     decide how big to make the loop over the "a" direction
!     
            if1=int((rcut-org(1))/latlen(ii(1)))
            is1=int((-rcut-org(1))/latlen(ii(1)))
            if (if1.gt.is1) then
               ic1=1
            else
               ic1=-1
            endif
!     
!     loop in the a direction
!     
            do i1=is1,if1,ic1
               notzero1=i1.ne.0
               if(notzero1.or.notzero2.or.notzero3) then
                  rlen2=0.
                  do j=1,3
                     r(j) = i1*mat%bvec(j,ii(1))    
     $                    + i2*mat%bvec(j,ii(2)) 
     $                    + i3*mat%bvec(j,ii(3))
                     rlen2 = rlen2 + r(j)*r(j)
                  enddo
                  if(rlen2.le.rcut2) then
!     
!     add all atoms associated with this BL site
!     
                     do ib=1,mat%nbasis
                        do j=1,3
                           rb(j) = r(j) + mat%basis(ii(1),ib)*mat%bvec(j
     $                          ,ii(1))    
     $                          + mat%basis(ii(2),ib)*mat%bvec(j,ii(2))
     $                          + mat%basis(ii(3),ib)*mat%bvec(j,ii(3))
                        enddo
!     add an atom
                        at%natoms=at%natoms+1
!     
!     make sure the assumed array size is still okay.  If not, start
!     again.
!     
                        if(at%natoms.gt.nguess) then
                           write(*,*) 
     $           '***WARNING: BuildCluster had to reallocate!'
                           nguess=nguess*2
                           allocate (at%x(3,nguess))
                           allocate (at%rlat(mat%nbasis,nguess))
                           allocate (at%spec(nguess))
                           allocate (at%lattice(nguess))
!C@
                           at%x = 0.
                           at%rlat = 0.
                           at%spec = 0
                           at%lattice = 0

                           go to 1
                        endif
                        at%x(1:3,at%natoms)=rb
                        at%spec(at%natoms)=mat%ispec(ib)
                        at%lattice(at%natoms)=ib
! compute rlat
                        do ib2=1,mat%nbasis
                           vec=at%x(1:3,at%natoms)-at%x(1:3,ib2)
                           at%rlat(ib2,at%natoms)=sqrt(Dot_Product(vec
     $                          ,vec))
                        enddo
                     enddo
                  endif
               endif
            enddo
         enddo
      enddo
!
! Now allocate the array to be returned, and copy the data from
!  the temporary storage
!
      atoms%natoms=at%natoms
      allocate (atoms%x(3,at%natoms))
!C@
      atoms%x = 0.
      atoms%x(1:3,1:at%natoms)=at%x(1:3,1:at%natoms)
      allocate (atoms%rlat(mat%nbasis,at%natoms))
      atoms%rlat = 0.
      atoms%rlat(1:mat%nbasis,1:at%natoms)=at%rlat(1:mat%nbasis
     $     ,1:at%natoms)
      allocate (atoms%spec(at%natoms))
      atoms%spec = 0
      atoms%spec(1:at%natoms)=at%spec(1:at%natoms)
      allocate (atoms%lattice(at%natoms))
      atoms%lattice = 0
      atoms%lattice(1:at%natoms)=at%lattice(1:at%natoms)
!
! deallocate the temporary storage
!
      call deallocate_cluster(at)
      end subroutine BuildCluster

   !------------------------------------------------------------------
   ! deallocate_cluster -- deallocate a type(cluster) variable
   !
   !      Passed Parameters :
   !            c  (in) : type(cluster) to be deallocated
   !
   !      Notes :
   !         It is important to deallocate all temporary clusters after
   !         use!
   !
   !      Author :
   !            R. Miller (01/17/98)
   !
   !      Revisions :
   !              none
   !
   !--

      subroutine deallocate_cluster(c)
      implicit none
      type(cluster) c
      deallocate (c%x)
      deallocate (c%rlat)
      deallocate (c%spec)
      deallocate (c%lattice)
      end subroutine deallocate_cluster

   !------------------------------------------------------------------
   ! SortBravais -- Pre-processing for buildcluster
   !
   !      Passed Parameters :
   !            lat   (in): bravais lattice vectors
   !            latlen(in): lengths of lat
   !            ii   (out): sort map for lat a=1,b=2,c=3
   !            bperp(out): component of b perp. to a
   !            bperp(out): component of c perp. to a
   !            per3 (out): projection of c on a-b plane
   !            per2 (out): projection of b on a
   !
   !      Algorithm :
   !         make the shortest vector "a", make "b" the one
   !         that is least perpendicular to "a" and "c" is the leftover
   !
   !      Notes :
   !         Needs to be rigourously tested for speed and optimized
   !         if possible - this will be called for every call to "cauchyborn"
   !
   !      Author :
   !            R. Miller (01/17/98)
   !
   !      Revisions :
   !              none
   !
   !--

      subroutine SortBravais(lat,latlen,ii,bperp,cperp,per3,per2)
      implicit none
      integer ii(3),j,i,itemp
      double precision lat(3,3),bperp,per3(2),latlen(3),t3(3),per2,cperp
     $     ,t2(3),bperp2,bperp3,temp2,temp3,per22,per23
!
! make the shortest vector "a"
!
      ii(1)=1
      do i=2,3
         if(latlen(i).lt.latlen(ii(1))) ii(1)=i
      enddo
      ii(2)=mod(ii(1),3)+1
      ii(3)=mod(ii(2),3)+1
!
! bperp is the length of the component of b perp. to a.
! bperp2 is using vector 2 as b,
! bperp3 is using vector 3 as b.
!
! choose the shorter of bperp2,bperp3 to be bperp
!
      per22=dot_product(lat(1:3,ii(1)),lat(1:3,ii(2)))/latlen(ii(1))
      temp2=per22/latlen(ii(1))
      per23=dot_product(lat(1:3,ii(1)),lat(1:3,ii(3)))/latlen(ii(1))
      temp3=per23/latlen(ii(1))
      bperp2=0.
      bperp3=0.
      do j=1,3
         t2(j)=lat(j,ii(2))-temp2*lat(j,ii(1))
         bperp2=bperp2+t2(j)**2
         t3(j)=lat(j,ii(3))-temp3*lat(j,ii(1))
         bperp3=bperp3+t3(j)**2
      enddo
      if(bperp2.lt.bperp3) then
         bperp=dsqrt(bperp2)
         per2=per22
         per3(1)=dot_product(lat(1:3,ii(3)),lat(1:3,ii(1)))
     $        /latlen(ii(1))
         per3(2)=dot_product(lat(1:3,ii(3)),t2)/bperp
      else
         bperp=dsqrt(bperp3)
         per2=per23
         itemp=ii(2)
         ii(2)=ii(3)
         ii(3)=itemp
         per3(1)=dot_product(lat(1:3,ii(3)),lat(1:3,ii(1)))
     $        /latlen(ii(1))
         per3(2)=dot_product(lat(1:3,ii(3)),t3)/bperp
      endif
!
      call Cross_Product(lat(1:3,ii(1)),lat(1:3,ii(2)),t3)
      cperp=dot_product(t3,lat(1:3,ii(3)))/sqrt(dot_product(t3,t3))
      end subroutine SortBravais
      end module mod_cluster
!


