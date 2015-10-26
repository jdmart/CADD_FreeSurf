c//   Qu modified on 01/11/2005 to check if the dected dislocation
c//   related to 
c//   thermal fluctuation. One subroutine findNearburgers() is added
c//   into.
c     
c//   Qu modified on 05/18/2006 to use detection band rings in the
c//   atomistic region.
C//   when dislocation is detected, the corresponding ring is removed,
C//   which makes
C//   the # of the detection band rings decreased.
c     
c     this file contains routines that manage the passing of
c     dislocations to and from the atomistic region
c     
c     global module variables:
c     
c     nburger - number of possible burgers vectors for this crystal
c     structure.  Note that this only includes those whose slip plane
c     normals or in the x-y plane in 2D CADD.
c     nslip - number of detection band elements
c     utilde - current u~ field superimposed of all the existing DD's
c     epsloc(1:3,1:3,i) - strain in detection band (DB) element i
c     amat(1:3,1:2,i) - matrix of shape function derivatives in DB
c     element i
c     enorm - for a given DB element, contains the L2 norm of epsloc
c     -epslib
c     imap(i) - global element number of DB element i
c     iburg(i) - the burgers vector type found in DB element i
c     possible(1:nburger,1:nslip) - .true. if a particular burgers
c     vector is possible in a particular DB element.  This means that
c     the burgers vector is not parallel to the "entry side" of the
c     element.  The entry side is the side of the element facing the
c     atomistics - into which we expect dislocations to enter the DB.
c     newslip - flag that signals when new DD's have been added,
c     necessitating an update of utilde.
c     burg(1:3,i) - library entry of the burger's vector for dislocation
c     i
c     normal(1:3,i) - library entry of the slip plane normal for
c     dislocation i
c     epslib(1:3,1:3,i) - eigenstrain for dislocation i
c     flib(1:3,1:3,i) - deformation gradient for dislocation i
c     
      module mod_dislocation
      integer nburger,nslip
      double precision, allocatable:: utilde(:,:),epsloc(:,:,:),amat(:,:
     $     ,:),enorm(:)
      integer, allocatable:: imap(:),iburg(:)
      logical, allocatable:: possible(:,:)
      logical newslip
      double precision, pointer:: burg(:,:),normal(:,:)
      double precision, pointer:: epslib(:,:,:),flib(:,:,:)
      end module mod_dislocation
c**********************************************************************
c     
c     NewDisocation: provides a means to insert continuum dislocations
c     "by
c     hand" using the macro "newd":
c     
c     newd,direct,xpos,ypos,bx,by,bz,itheta_e,theta_s
c     
c     where:
c     xpos,ypos - dislocation coordinates
c     bx,by,bz  - burgers vector
c     itheta_e  - angle of branch cut of edge field, in multiples
c     of PI, measured from the direction of b
c     theta_s   - angle of branch cut of the screw field, measured
c     from pos. x axis.
c     
c     or:
c     
c     newd,file,filename
c     
c     to put multiple dislocations in a separate file.
c     
      subroutine NewDislocation(input,x,b,IsRelaxed,numnp,nxdm,ndf)
      use mod_file
      implicit none
c     
c     transferred variables
c     
      integer ndf,numnp,nxdm
      integer IsRelaxed(numnp),x(nxdm,numnp),b(ndf,numnp)
      character input*80
c     
c     local variables
c     
      double precision disx(2),disb(3),theta_e,theta_s,PI
c---  PI is deliberately less than full 3.14159 to be PI-epsilon
      parameter (PI=3.1415)

      integer upper,lower,next,idum,logic,nnewdis,inew,itheta_e
      character*80 filename
      character*4 key
      key=input(1:4)

      lower=4
      upper = next(lower,input)
      if(key.eq.'dire') then
         call freein(input,lower,upper,idum,disx(1),2)
         lower=upper
         upper=next(lower,input)
         call freein(input,lower,upper,idum,disx(2),2)
         lower=upper
         upper=next(lower,input)
         call freein(input,lower,upper,idum,disb(1),2)
         lower=upper
         upper=next(lower,input)
         call freein(input,lower,upper,idum,disb(2),2)
         lower=upper
         upper=next(lower,input)
         call freein(input,lower,upper,idum,disb(3),2)
         lower=upper
         upper=next(lower,input)
         call freein(input,lower,upper,itheta_e,idum,1)
         lower=upper
         upper=next(lower,input)
         call freein(input,lower,upper,idum,theta_s,2)
         nnewdis=1
      else
         filename=input(1:(upper-1))
         if(.not.FileExists(filename,.false.)) then
            write(*,*) '** WARNING: no dislocation file found'
            return
         endif
         call iofile(filename,'formatted  ',logic,.true.)
         read(logic,*) nnewdis
         read(logic,*) disx,disb,itheta_e,theta_s
      endif
      do inew=1,nnewdis
         if(inew.gt.1) read(logic,*) disx,disb,itheta_e,theta_s
         theta_e=PI*mod(itheta_e,2)

         call disl_pass(disx,disx,disb,theta_e,theta_s,x,b,IsRelaxed
     $        ,numnp,.false.,.true.)
         write(*,*) '    New Dislocation Added at:',disx
         write(*,*) '         with Burgers vector:',disb
         write(*,*) '                     theta_e:',theta_e
         write(*,*) '                 and theta_s:',theta_s
      enddo
      if(key.ne.'dire') close(logic)
      end subroutine NewDislocation
c**********************************************************************
c     
c     dislcheck:
c     main point of contact for checking if dislocations want to pass
c     either way.  
c     
      logical function dislcheck(CheckSlip,LostSlip,AddedSlip,MoveDisl
     $     ,ix,x,b,itx,IsRelaxed,numnp,ndf,nxdm,numel,nen1,newmesh
     $     ,plottime, dislpass, npass)
      implicit none
      integer npass 
      logical CheckSlip,LostSlip,AddedSlip,MoveDisl,newmesh
      integer numnp,ndf,nxdm,numel,nen1
      integer ix(nen1,numel),itx(3,numel),IsRelaxed(numnp)
      double precision x(nxdm,numnp),b(ndf,numnp),plottime
      logical dislpass
      dislcheck=.false.
c     
c     if we are checking for dislocation passings, proceed, otherwise
c     just return
c     
      if(CheckSlip) then
c     
c     only check for dislocations leaving the continuum if the
c     dislocations are mobile.  Lostslipcheck looks for dislocations
c     that want to go from continuum to atomistic
c     
         if(MoveDisl) then
            call lostslipcheck(LostSlip,ix,x,b, npass)
            if(LostSlip) then
               dislcheck=.true.
               return
            endif
         endif
c     
c     check detection band for dislocations that want to go from
c     atomistic to continumm
c     
         if (.not. LostSlip) then 
            call slipcheck(x,b,ix,itx,IsRelaxed,numnp,ndf,nxdm,numel
     $           ,nen1,newmesh,AddedSlip,plottime, dislpass)
            LostSlip = .false. 
         
            if(AddedSlip) then
               dislcheck=.true.
               return
            endif
         end if
      endif
      end function dislcheck
c*********************************************************************
c     
c     slipcheck:  checks detection band for dislocations leaving the
c     atomistic region.
c     
      subroutine slipcheck(x,b,ix,itx,IsRelaxed,numnp,ndf,nxdm,numel
     $     ,nen1,newmesh,AddedSlip,plottime, dislpass)
      use mod_dislocation
      use mod_file
      use mod_boundary
      use mod_parallel
      implicit none
      Logical  newmesh,AddedSlip
      integer numnp,ndf,numel,nen1,nxdm,i,n1,n2,n3
      double precision x(nxdm,numnp),b(ndf,numnp),det,plottime
      integer ix(nen1,numel),IsRelaxed(numnp),itx(3,*)
      logical dislpass

      integer j1,j2,j3,i1,i2,i3,node1,node2,j,k1,k2,iel,k,kmod
      integer itotal,itheta,kp1,kp2,ndisl,numnew,ifactor,idb
      integer, volatile :: idbCopy, jcopy
      logical found(3),InContinuum,InAtoms
      character*80 filename
      double precision dvec(2,3),delp(2),delm(2),xav(2),x3(2)
     $     ,theta_e,theta_s,xd(3),x0(3),xi(3),cross,cvec(2),bvec(3)
      double precision LENGTHTOL,LENGTHTOL2,PI
c     Qu modification begins
c     
      integer neleNear,neleNext,nsideNext,nsideRight,nsideLeft,ibNext
     &     ,nsideNearR,nsideNearL,nside,kk
      double precision amatNext(3,2),epslocNext(3,3),coordSide(2)
     &     ,coordRight(2),coordLeft(2),vecRight(2),vecLeft(2)
     &     ,vecLength,dotRight,dotLeft
      logical possibleNext(nburger)
      logical, save, allocatable :: examined(:,:)
      integer, save :: oidb
      double precision :: s_dis

      integer, volatile :: islp

c     Qu modification ends
c     
c---  PI is deliberately less than full 3.14159 to be PI-epsilon
      parameter (LENGTHTOL=1.e-4, LENGTHTOL2=1.e-8, PI=3.1415)
c     
c     allocate and initialize data structure
c     
      AddedSlip=.false.
      if(newmesh) then
         oidb=0
C     write(*,*) '** First Entry to slipcheck: initializing'
         newmesh=.false.
         newslip=.true.
         if(allocated(epsloc)) then
            deallocate(epsloc,imap,utilde,amat,
     &           iburg,possible,enorm,elidb,examined)
         endif
c--   count the slip detection elements
c     Qu modified detection band rings starts
         nslip=0
         do idb=ndbpoly,1,-1
            do i=1,nelidb(idb)
               nslip=nslip+1
            enddo
         enddo
         print *, 'Number of detection band slip elements =',nslip
c     do i=1,numel
c     if(ix(nen1,i).lt.0) then
c     nslip=nslip+1
c     endif
c     enddo
c--   allocate
         if (.not. allocated(epsloc)) then 
         allocate(epsloc(3,3,nslip),imap(nslip),utilde(3,numnp),amat(3,2
     $        ,nslip),iburg(nslip),possible(nburger,nslip),enorm(nburger
     $        ),elidb(nslip),examined(nburger,numel))
         else
            epsloc = 0.0d0
            imap = 0
            utilde = 0.0d0
            amat = 0.0d0
            iburg = 0
            enorm = 0.0d0
            elidb = 0
            examined = 0
         end if

c--   compute imap
         nslip=0
         do idb=ndbpoly,1,-1
            do i=1,nelidb(idb)
               nslip=nslip+1
               imap(nslip)=eldb(i,idb)
               elidb(nslip)=idb
            enddo
         enddo
c     do i=1,numel
c     if(ix(nen1,i).lt.0) then
c     nslip=nslip+1
c     imap(nslip)=i
c     endif
c     enddo
c     Qu modified detection band rings ends
c     
c--   pre-process the slip elements (compute amat for each)
c     
         do i=1,nslip
            n1=ix(1,imap(i))
            n2=ix(2,imap(i))
            n3=ix(3,imap(i))
            det =(x(1,n1)-x(1,n3))*(x(2,n2)-x(2,n3)) -
     &           (x(1,n2)-x(1,n3))*(x(2,n1)-x(2,n3))
            amat(1,1,i)=(x(2,n2)-x(2,n3))/det
            amat(2,1,i)=(x(2,n3)-x(2,n1))/det
            amat(3,1,i)=(x(2,n1)-x(2,n2))/det
            amat(1,2,i)=(x(1,n3)-x(1,n2))/det
            amat(2,2,i)=(x(1,n1)-x(1,n3))/det
            amat(3,2,i)=(x(1,n2)-x(1,n1))/det
c     
c     only allow burgers vectors that are not parallel to entry side of
c     the element.
c     
            possible(nburger,i)=.true.
            k=abs(ix(nen1,imap(i)))
            if(k.ne.0) then
               kp1=mod(k,3)+1
               k=ix(k,imap(i))
               kp1=ix(kp1,imap(i))
               do j=1,nburger-1
                  cross=burg(1,j)*(x(2,k)-x(2,kp1))-burg(2,j)*(x(1,k)
     $                 -x(1,kp1))
                  possible(j,i)=(abs(cross).gt.LENGTHTOL)
c     dw hack
c$$$  if (j.ne.10) possible(j,i)=.false.
c     end hack
               enddo
            else
               possible(1:nburger-1,i)=.false.
            endif
         enddo
      endif
c     
c     subtract all previously passed dislocations for computation of DB
c     elemental strains
c     

!     print *, 'In slipcheck and adding displ_contribution'
 300  if(newslip) then
         utilde=0
         do i=1,numnp
            call disl_displ(x(1:3,i),utilde(1:3,i))
         enddo
         newslip=.false.
      endif
      b(1:3,1:numnp)=b(1:3,1:numnp)-utilde(1:3,1:numnp)
c     
c     get current strain in each DB element
c     
      do i=1,nslip
         n1=ix(1,imap(i))
         n2=ix(2,imap(i))
         n3=ix(3,imap(i))
        call GetElementStrain(b(1,n1),b(1,n2),b(1,n3),amat(1:3,1:2,i)
     $        ,epsloc(1:3,1:3,i))
      enddo
c     
c     find nearest possible slip vector for each element
c     
      ndisl=0
      do i=1,nslip 
!     if(elidb(nslip).ne.1) then 
         call findburgers(epsloc(1:3,1:3,i),possible(1:nburger,i)
     $        ,iburg(i),ndisl,enorm)
!     endif
      enddo
c     **** Finished finding Burger's vector in elements **
      do i=1,nslip
         do j=1,nburger
            examined(j,imap(i))=.false.
         enddo
      enddo

c$$$c     Qu's modification begins
c$$$c     check if the detected slip is related to thermal fluctuation
c$$$      if(ndisl.gt.0) then
c$$$         do i=1,nslip
c$$$
c$$$c     Dw's mod begin
c$$$            if(examined(iburg(i),imap(i))) then
c$$$               ndisl=ndisl-1
c$$$               iburg(i)=nburger
c$$$               goto 21
c$$$            endif
c$$$c$$$            print *, 'XXXX Dislocations in Atomistic = ', ndisl
c$$$
c$$$c     Dw's mod end
c$$$
c$$$            if(iburg(i).ne.nburger) then
c$$$               call checkburgers(epsloc(1:3,1:3,i),possible(1:nburger,i)
c$$$     &              ,iburg(i),x,b,ix,imap(i),enorm)
c$$$c     find the neighbor element of the detection band element i
c$$$               neleNear=itx(abs(ix(nen1,imap(i))),imap(i))
c$$$c--   pre-process the element (compute amat)
c$$$c     
c$$$               if(neleNear.ne.0) then
c$$$                  possibleNext(1:nburger)=possible(1:nburger,i)
c$$$                  k=abs(ix(nen1,imap(i)))
c$$$c     vector from centroid to midside of edge k
c$$$                  kp1=mod(k,3)+1
c$$$                  kp2=mod(kp1,3)+1
c$$$                  k=ix(k,imap(i))
c$$$                  kp1=ix(kp1,imap(i))
c$$$                  kp2=ix(kp2,imap(i))
c$$$                  cvec=(x(1:2,k)+x(1:2,kp1)-2*x(1:2,kp2))/6.d0
c$$$                  cross=normal(1,iburg(i))*cvec(2)-normal(2,iburg(i))
c$$$     $                 *cvec(1)
c$$$                  if(cross.lt.0.d0) then
c$$$                     bvec=-burg(1:3,iburg(i))
c$$$                  else
c$$$                     bvec=burg(1:3,iburg(i))
c$$$                  endif
c$$$c     
c$$$                  do j=1,3
c$$$                     if(itx(j,neleNear).eq.imap(i))then
c$$$                        nside=j
c$$$                     endif
c$$$                  enddo
c$$$                  nsideRight=mod(nside,3)+1
c$$$                  nsideLeft=mod(nsideRight,3)+1
c$$$
c$$$                  n1=ix(nside,neleNear)
c$$$                  n2=ix(nsideRight,neleNear)
c$$$                  n3=ix(nsideLeft,neleNear)
c$$$c     find the coord of the mid point of nside
c$$$                  coordSide(1:2)=(x(1:2,n1)+x(1:2,n2))/2.d0
c$$$c     find the coord of the mid point of nsideRight
c$$$                  coordRight(1:2)=(x(1:2,n2)+x(1:2,n3))/2.d0
c$$$c     find the coord of the mid point of nsideLeft
c$$$                  coordLeft(1:2)=(x(1:2,n3)+x(1:2,n1))/2.d0
c$$$c     
c$$$c     find vector from mid isideRight to mid j
c$$$                  vecRight(1:2)=coordRight(1:2)-coordSide(1:2)
c$$$                  vecLength=dot_product(vecRight(1:2),vecRight(1:2))
c$$$                  vecRight(1:2)=vecRight(1:2)/vecLength
c$$$c     find vector from mid isideLeft to mid j
c$$$                  vecLeft(1:2)=coordLeft(1:2)-coordSide(1:2)
c$$$                  vecLength=dot_product(vecLeft(1:2),vecLeft(1:2))
c$$$                  vecLeft(1:2)=vecLeft(1:2)/vecLength
c$$$                  
c$$$                  dotRight=dot_product(bvec(1:2),vecRight(1:2))
c$$$                  dotLeft=dot_product(bvec(1:2),vecLeft(1:2))
c$$$                  if(dotRight.gt.dotLeft)then
c$$$                     nsideNext=nsideRight
c$$$                  else
c$$$                     nsideNext=nsideLeft
c$$$                  endif
c$$$                  neleNext=neleNear
c$$$ 10               n1=ix(1,neleNext)
c$$$                  n2=ix(2,neleNext)
c$$$                  n3=ix(3,neleNext)
c$$$
c$$$                  det =(x(1,n1)-x(1,n3))*(x(2,n2)-x(2,n3)) -
c$$$     &                 (x(1,n2)-x(1,n3))*(x(2,n1)-x(2,n3))
c$$$                  amatNext(1,1)=(x(2,n2)-x(2,n3))/det
c$$$                  amatNext(2,1)=(x(2,n3)-x(2,n1))/det
c$$$                  amatNext(3,1)=(x(2,n1)-x(2,n2))/det
c$$$                  amatNext(1,2)=(x(1,n3)-x(1,n2))/det
c$$$                  amatNext(2,2)=(x(1,n1)-x(1,n3))/det
c$$$                  amatNext(3,2)=(x(1,n2)-x(1,n1))/det
c$$$c     
c$$$c     only allow burgers vectors that are not parallel to entry side of
c$$$c     the element.
c$$$c     
c$$$                  possibleNext(nburger)=.true.
c$$$                  kk=abs(ix(nen1,neleNext))
c$$$                  if(kk.ne.0)then
c$$$                     k=ix(nsideNext,neleNext)	           
c$$$                     kp1=mod(nsideNext,3)+1
c$$$                     kp1=ix(kp1,neleNext)
c$$$                     do j=1,nburger-1
c$$$                        cross=burg(1,j)*(x(2,k)-x(2,kp1))-burg(2,j)
c$$$     &                       *(x(1,k)-x(1,kp1))
c$$$                        possibleNext(j)=(abs(cross).gt.LENGTHTOL)
c$$$                     enddo
c$$$                  else
c$$$                     possibleNext(1:nburger-1)=.false.
c$$$                  endif
c$$$c     
c$$$c     get current strain in each DB element
c$$$c     
c$$$                  n1=ix(1,neleNext)
c$$$                  n2=ix(2,neleNext)
c$$$                  n3=ix(3,neleNext)
c$$$                  call GetElementStrain(b(1,n1),b(1,n2),b(1,n3),
c$$$     &                 amatNext(1:3,1:2),epslocNext(1:3,1:3))
c$$$     $                 
c$$$c     find nearest possible slip vector for each element
c$$$                  call findNextburgers(epslocNext,possibleNext,ibNext)
c$$$                  if(ibNext.ne.nburger)then
c$$$                     if(itx(nsideNext,neleNext).eq.0)then
c$$$	                goto 20
c$$$                     endif
c$$$                     if(ibNext.eq.iburg(i)) then
c$$$                        examined(iburg(i),neleNext)=.true.
c$$$                     endif
c$$$                     neleNear=neleNext
c$$$                     nsideNearR=nsideRight
c$$$                     nsideNearL=nsideLeft
c$$$                     neleNext=itx(nsideNext,neleNear)
c$$$                     do j=1,3
c$$$                        if(itx(j,neleNext).eq.neleNear)then
c$$$                           nside=j
c$$$                        endif
c$$$                     enddo
c$$$                     nsideRight=mod(nside,3)+1
c$$$                     nsideLeft=mod(nsideRight,3)+1
c$$$                     if(nsideNext.eq.nsideNearR)then
c$$$                        nsideNext=nsideLeft
c$$$                     else
c$$$                        nsideNext=nsideRight
c$$$                     endif
c$$$                     goto 10
c$$$                  else
c$$$                     ndisl=ndisl-1
c$$$                     iburg(i)=nburger
c$$$                  endif
c$$$               else
c$$$                  write(*,*)'!Warning: No past Dislocation path'
c$$$               endif
c$$$ 20            continue
c$$$            endif
c$$$ 21         continue
c$$$c$$$            print *, 'NSLIP = ', nslip, nburger
c$$$	 enddo
c$$$      end if
c$$$c     Qu's modification ends

c     Qu modified detection band rings starts
!     if(ndisl.gt.1)then
!     write(*,*)'ERROR----, more than one disl'
c     stop
!     endif
      if(ndisl.gt.0)then
!     find outermost detection ring that is triggered
         j=0
         do i=1,nslip
            if(j.eq.0.and.iburg(i).ne.nburger) then
               idb=elidb(i)
               j=1
               idbCopy = idb
               jCopy = j
            endif
         enddo

         if(idb.le.ndbpoly)then
            do i=1,nslip
               if(iburg(i).ne.nburger.and.idb.ne.oidb) then
                  if(ndisl.gt.1)then
                     do j=1,nslip
                        if(j.ne.i) then
                           if(elidb(j).eq.3.and.iburg(j).ne.nburger)
     $                          then
!     AddedSlip=.true.
                              write(*,*)'!!!!!more than one disl on '
     $                             ,rank
!     ndisl=0
                           endif
                        endif
                     enddo
                     write(*,*)'ERROR----, more than one disl'
c     stop
                  endif
c     
c     recompute enorms
c     
                  call findburgers(epsloc(1:3,1:3,i),possible(1:nburger
     $                 ,i),iburg(i),ndisl,enorm)
c     
c     resolve degeneracies - sometimes 2 or more dislocations produce
c     the
c     same strain matrix.  checkburgers chooses the b that
c     produces the best fit with the rotation of the element
c     
                  call checkburgers(epsloc(1:3,1:3,i),possible(1:nburger
     $                 ,i),iburg(i),x,b,ix,imap(i),enorm)



c$$$                  write(*,*)
c$$$                  write(*,*) 'slip found in element',imap(i),abs(ix(nen1
c$$$     $                 ,imap(i))), elidb(i), ndbpoly, ' :'m
c$$$

!     write(*,*) x(1:2,ix(1,imap(i))),b(1:3,ix(1,imap(i)))
!     write(*,*) x(1:2,ix(2,imap(i))),b(1:3,ix(2,imap(i)))
!     write(*,*) x(1:2,ix(3,imap(i))),b(1:3,ix(3,imap(i)))
!     write(*,*)
!     write(*,*) 'strain matrix:'
!     write(*,'(3e15.6)') (epsloc(j,1:3,i),j=1,3)

c$$$                  write(*,*) 'burgers vector:',iburg(i)
c$$$                  write(*,*) burg(1:3,iburg(i))
c$$$                  x0=0.d0
c$$$                  do k=1,3
c$$$                     x0(1:2)=x0(1:2)+x(1:2,ix(k,imap(i)))
c$$$                  enddo
c$$$                  x0(1:2)=x0(1:2)/3.d0
c$$$
c$$$                  write(*,*) 'dislocation at ',x0(1:2)
c$$$                  write(*,*)'time = ', plottime,'ps'
c$$$                  x0(3)=dsqrt(x0(2)*x0(2)+x0(1)*x0(1))
c$$$                  write(*,107) rank,iburg(i),plottime,x0(3),x0(1),x0(2)
c$$$ 107              format (I4,I4,' disdata ',4f10.3)
c$$$                  write(*,*)
                  oidb=idb

               endif
            enddo
!     write(*,*) 'dislocation at ',x0
!     write(*,*)'time = ', plottime,'ps'
!     write(*,*)
!     nslip=nslip-nelidb(idb)
!            ndisl=0
         endif
      endif
      if (ndisl .gt. 0) then 
         print *, 'Outermost Detection Band', idb
      end if

c     Qu modified detection band rings ends
c     
c     put in the new dislocations
c     
      if(ndisl.gt.0) then
         write(*,*) '   ***** New Slip Detected *****   '
         numnew=0
c$$$         filename='out/detection.plt'
c$$$         call plottrigger(x,b,ix,numel,numnp,nen1,ndf,nxdm,iburg,imap
c$$$     $        ,nslip,utilde,filename,nburger)
c$$$         call plotdisp(x,b,ix,numel,numnp,nen1,ndf,nxdm,IsRelaxed
c$$$     $        ,filename)
         newslip=.true.
         AddedSlip=.true.
c     
c     assume core is at center of each slipped element.
c     
c     find x0, xd, theta
c     


!     dw hack
!     goto 505

         do i=1,nslip
            if(iburg(i).ne.nburger) then
c     
c     recompute enorms
c     
               call findburgers(epsloc(1:3,1:3,i),possible(1:nburger,i)
     $              ,iburg(i),ndisl,enorm)
c     
c     resolve degeneracies - sometimes 2 or more dislocations produce
c     the
c     same strain matrix.  checkburgers chooses the b that
c     produces the best fit with the rotation of the element
c     
               call checkburgers(epsloc(1:3,1:3,i),possible(1:nburger,i)
     $              ,iburg(i),x,b,ix,imap(i),enorm)
               write(*,*) 'slip found in element',imap(i),abs(ix(nen1
     $              ,imap(i))), elidb(i), ndbpoly,' :'
c$$$               write(*,*) x(1:2,ix(1,imap(i))),b(1:3,ix(1,imap(i)))
c$$$               write(*,*) x(1:2,ix(2,imap(i))),b(1:3,ix(2,imap(i)))
c$$$               write(*,*) x(1:2,ix(3,imap(i))),b(1:3,ix(3,imap(i)))
c$$$               write(*,*)
c$$$               write(*,*) 'strain matrix:'
c$$$               write(*,'(3e15.6)') (epsloc(j,1:3,i),j=1,3)
c$$$               write(*,*)
               write(*,*) 'burgers vector:',iburg(i),burg(1:3,iburg(i))

c     
c     process the dislocation:  determine its location, the direction of
c     its
c     branch cut
c     x0: initial location (centroid of DB element)
c     xd: place in the continuum to which disl will be moved.
c     xi: location of the "image" dislocation as far from any continuum
c     region as possible but still on the slip plane.
c     
c--   ix(nen1,i) stores the side of the element into which dislocations
c--   may
c     pass
               k=abs(ix(nen1,imap(i)))
               if(k.ne.0) then
c     vector from centroid to midside of edge k
c     Check if dislocation is in the last detection band polygon
                  idb = elidb(i)
                  idbCopy = idb
c     Pass dislocations only if the stored Outer polygon is detected
                  kp1=mod(k,3)+1
                  kp2=mod(kp1,3)+1
                  k=ix(k,imap(i))
                  kp1=ix(kp1,imap(i))
                  kp2=ix(kp2,imap(i))
                  cvec=(x(1:2,k)+x(1:2,kp1)-2*x(1:2,kp2))/6.d0
                  cross=normal(1,iburg(i))*cvec(2)-normal(2,iburg(i))
     $                 *cvec(1)
                  if(cross.lt.0.d0) then
                     bvec=-burg(1:3,iburg(i))
                  else
                     bvec=burg(1:3,iburg(i))
                  endif
                  if(dot_product(bvec(1:2),cvec(1:2)).lt.0.d0) then
                     itheta=1
                     theta_s=datan2(-bvec(2),-bvec(1))
                  else
                     itheta=0
                     theta_s=datan2(bvec(2),bvec(1))
                  endif
                  x0=0.d0
                  do k=1,3
                     x0(1:2)=x0(1:2)+x(1:2,ix(k,imap(i)))
                  enddo
                  x0(1:2)=x0(1:2)/3.d0
c     
c     compute location for the continuum core.  FindEntryPoint finds the
c     place where the slip plane intersects the atom/continuum
c     interface.
c     Then we move it 1 burgers vector further along to get it out of
c     the continuum detection region.
c     
                  numnew=numnew+1
c     
c     Check of Detection band polygon is a boundary (dbboundnear))
c     pass it ony if it is in this polygon
                  if (dbboundnear(idb)) then 
                     if(itheta.eq.0) then
                        ifactor=-1
                     else
                        ifactor=1
                     endif
                     islp = 0
                     s_dis = 0.0d0
                     call FindEntryPoint(bvec,x0,xd, ifactor)
                     call FindSlipPlane(bvec, x0, ifactor, islp, theta_s
     $                    , s_dis, xd, xi)
c     
c--   choose location for the "image" dislocation as far from any
c--   detection bands as possible.
c     
c$$$                     call FindImageLocation(xi,ifactor,x0,bvec(1:3),ix,x
c$$$     $                    ,nxdm,numnp,numel,nen1)
                     theta_e=itheta*PI
 3                   write(*,*) 'dislocation at ',x0
                     write(*,*) 'passed to ',xd
                     write(*,*) 'Image Location', xi
                     write(*,*) 'with b=',bvec(1:3)
                     write(*,*) 'and theta_e=',theta_e
                     write(*,*) 'and theta_s=',theta_s

                     call disl_pass(x0,xd,bvec(1:3),theta_e,theta_s,x,b
     $                    ,IsRelaxed,numnp,.true.,.true.,0, islp, s_dis)
c$$$                     itheta=mod(itheta+1,2)
c$$$                     theta_e=itheta*PI
c$$$                     call disl_pass(xi,xi,-bvec(1:3),theta_e,theta_s,x
c$$$     $                    ,b,IsRelaxed,numnp,.true.,.true.)
                     call disl_print(0)
 1000                format(5e15.6)
                     dislpass = .true. 
                  end if
               else
                  write(*,*)
     $                 'warning: slip in an element on the interior of'
                  write(*,*) '       the detection band'
               endif
 500           continue
            endif
         enddo
         if(numnew.eq.0) stop 'ERROR: couldn''t resolve dislocations'
         filename='out/detection.plt'
         call plotdisp(x,b,ix,numel,numnp,nen1,ndf,nxdm,IsRelaxed
     $        ,filename)
      endif

 505  continue
c     
c     restore b-vector
c     
      b(1:3,1:numnp)=b(1:3,1:numnp)+utilde(1:3,1:numnp)

      end subroutine slipcheck
c*******************************************************************
c     
c     GetElementStrain:  compute strain in an element
c     
      subroutine GetElementStrain(b1,b2,b3,amat,eps)
      implicit none
      double precision amat(3,2),b1(3),b2(3),b3(3),u(3,3),eps(3,3),ua(3
     $     ,3)
      u(1:3,1)=b1
      u(1:3,2)=b2
      u(1:3,3)=b3
      ua(1:3,1:2)=matmul(u,amat)
      ua(1:3,3)=0.d0
      eps=matmul(transpose(ua),ua)
      eps=eps+ua+transpose(ua)
      end subroutine GetElementStrain
c*********************************************************************
c     
c     checkburgers:
c     resolve degeneracies in the dislocation strain library.
c     
      subroutine checkburgers(eps,p,ib,x,b,ix,iel,test)
      use mod_global
      use mod_dislocation
      implicit none
      double precision eps(3,3),x(nxdm,*),b(ndf,*),test(*)
      integer iel,ix(nen1,*),ib,ibest,i
      logical p(*)
c     
      double precision dx1(3),dx2(3),dy1(3),dy2(3),norm,del(3),normmin
     $     ,testmin,dz1(3),dz2(3)
      dx1=x(1:3,ix(2,iel))-x(1:3,ix(1,iel))+b(1:3,ix(2,iel))-b(1:3,ix(1
     $     ,iel))
      dx2=x(1:3,ix(3,iel))-x(1:3,ix(1,iel))+b(1:3,ix(3,iel))-b(1:3,ix(1
     $     ,iel))
      dy1=x(1:3,ix(2,iel))-x(1:3,ix(1,iel))
      dy2=x(1:3,ix(3,iel))-x(1:3,ix(1,iel))
      normmin=1.e30
      testmin=1.001*test(ib)
      do i=1,nburger-1
         if(p(i).and.test(i).lt.testmin) then
            dz1=matmul(flib(1:3,1:3,i),dy1)
            dz2=matmul(flib(1:3,1:3,i),dy2)
            del=dx1-dz1
            norm=dot_product(del,del)
            del=dx2-dz2
            norm=norm+dot_product(del,del)
            if(norm.lt.normmin) then
               normmin=norm
               ibest=i
            endif
         endif
      enddo
      if(ibest.ne.ib) then
c//   Qu modification begins
c     write(*,*) '**NOTICE: burgers vector changed in checkburgers'
c//   Qu modification ends
         ib=ibest
      endif
      end subroutine checkburgers
c*********************************************************************
c     
c     findburgers:
c     determine what burgers vector lies in a DB element.
c     
      subroutine findburgers(eps,p,ib,ndisl,test)
      use mod_dislocation
      implicit none
      double precision eps(3,3),test(*),xmin,del(3,3)
      integer ib,ndisl,i,k,j
      logical p(nburger)
      xmin=1.e30
      do i=1,nburger
         if(p(i)) then
            del=eps-epslib(1:3,1:3,i)
            test(i)=0
            do j=1,3
               do k=j,3
                  test(i)=test(i)+(del(j,k))**2
               enddo
            enddo
            if(test(i).lt.xmin) then
               ib=i
               xmin=test(i)
            endif
         else
            test(i)=1.e30
         endif
      enddo
      if(ib.eq.nburger) return
      ndisl=ndisl+1
      end subroutine findburgers
c*********************************************************************
c     Qu modification begins
c     
c     findNearburgers:
c     determine what burgers vector lies in a DB element.
c     
      subroutine findNextburgers(eps,p,ib)
      use mod_dislocation
      implicit none
      double precision eps(3,3),test(nburger),xmin,del(3,3)
      integer ib,ndislNear,i,k,j
      logical p(nburger)
      xmin=1.e30
      do i=1,nburger
         if(p(i)) then
            del=eps-epslib(1:3,1:3,i)
            test(i)=0
            do j=1,3
               do k=j,3
                  test(i)=test(i)+(del(j,k))**2
               enddo
            enddo
            if(test(i).lt.xmin) then
               ib=i
               xmin=test(i)
            endif
         else
            test(i)=1.e30
         endif
      enddo
      end subroutine findNextburgers
c**********************************************************************
c     Qu modification ends
c     
c     rotateburgers:
c     for a given crystal structure and orientation, find all the
c     possible burgers vectors and rotate them to the global coord
c     system.  Also build the strain and defm gradient libraries.
c     
      subroutine rotateburgers(q,rotmat,a0,struct)
      use mod_dislocation
      implicit none
      double precision q(3,3),rotmat(2),q2(3,3),a6,a62,x(3),a0,b(3),b1(3
     $     ),m(3),bx,bz,rt6,rt3,d,tol,rt2,b2(3)
      parameter(tol=1.e-6)
      integer i1,i2,i3,i,ib,ip,ipp,is,j,ip2,ip3,is2
      character*3 struct
      rt6=dsqrt(6.d0)
      rt3=dsqrt(3.d0)
      rt2=dsqrt(2.d0)
      if(struct.eq.'fcc') then
         nburger=25
      else if(struct.eq.'hex') then
         nburger=7
      else if(struct.eq.'bcc') then
         nburger=25
      else
         write(*,*) 'structure type:',struct
         stop 'not recognized in rotateburgers'
      endif
      allocate(epslib(3,3,nburger),flib(3,3,nburger),burg(4,nburger)
     $     ,normal(3,nburger))
      q2=0.d0
      q2(1,1)=rotmat(1)
      q2(2,2)=rotmat(1)
      q2(1,2)=-rotmat(2)
      q2(2,1)=rotmat(2)
      q2(3,3)=1.d0
      q2=matmul(q2,q)
c     
      nburger=0
      if(struct.eq.'fcc') then
         a6=a0/6.d0
         a62=2*a6
         d=a0/rt3
         do ip=0,3
            m(1:3)=1.d0/rt3
            if(ip.ne.0) m(ip)=-m(ip)
            m=matmul(q2,m)
            if(abs(m(3)).lt.tol) then
               do ib=1,3
                  call getb(struct,ip,ib,b1,a6,a62)
                  do is=-1,1,2
                     nburger=nburger+1
                     normal(1:3,nburger)=m
                     burg(1:3,nburger)=is*matmul(q2,b1)
                     burg(4,nburger)=sqrt(dot_product(b1,b1))
                     call getstrain(burg(1:3,nburger),d,m,epslib(1:3,1:3
     $                    ,nburger),flib(1:3,1:3,nburger))
                  enddo
               enddo
            endif
         enddo
      else if(struct.eq.'bcc') then
         a6=rt3*a0/2.d0
         d=a0/rt2
         do ip=1,3
            do is2=-1,1,2
               m(1:3)=1.d0/rt2
               m(ip)=0
               ip2=mod(ip,3)+1
               ip3=mod(ip2,3)+1
               m(ip2)=m(ip2)*is2
               b1(1:3)=a0/2.d0
               b2(1:3)=a0/2.d0
               b1(ip2)=-b1(ip2)*(m(ip2)/abs(m(ip2))) !sgn fcn
               b1(ip3)= b1(ip3)*(m(ip3)/abs(m(ip3))) !sgn fcn
               b2(ip2)=-b1(ip2)
               b2(ip3)=-b1(ip3)
               m=matmul(q2,m)
               if(abs(m(3)).lt.tol) then
                  do is=-1,1,2
                     nburger=nburger+1
                     normal(1:3,nburger)=m
                     burg(1:3,nburger)=is*matmul(q2,b1)
                     burg(4,nburger)=sqrt(dot_product(b1,b1))
                     call getstrain(burg(1:3,nburger),d,m,epslib(1:3
     $                    ,1:3,nburger),flib(1:3,1:3,nburger))
                     nburger=nburger+1
                     normal(1:3,nburger)=m
                     burg(1:3,nburger)=is*matmul(q2,b2)
                     burg(4,nburger)=sqrt(dot_product(b2,b2))
                     call getstrain(burg(1:3,nburger),d,m,epslib(1:3
     $                    ,1:3,nburger),flib(1:3,1:3,nburger))
                  enddo
               endif
            enddo
         enddo
C--   Qu modification begins
         d=a0/rt6
         do ip=1,3
            ip2=mod(ip,3)+1
            ip3=mod(ip2,3)+1
            do ipp=0,3
               m(ip)=2.d0/rt6
               m(ip2)=1.d0/rt6
               m(ip3)=1.d0/rt6
               if (ipp.ne.0) m(ipp)=-m(ipp)
               b1(1:3)=a0/2.d0
               do is=1,3
                  b1(is)=b1(is)*(m(is)/abs(m(is)))
               enddo
               b1(ip)=-b1(ip)
               m=matmul(q2,m)
               if (abs(m(3)).lt.tol)then
                  do is=-1,1,2
                     nburger=nburger+1
                     normal(1:3,nburger)=m
                     burg(1:3,nburger)=is*matmul(q2,b1)
                     burg(4,nburger)=sqrt(dot_product(b1,b1))
                     call getstrain(burg(1:3,nburger),d,m,epslib(1:3
     $                    ,1:3,nburger),flib(1:3,1:3,nburger))
                  enddo
               endif
            enddo
         enddo
C--   Qu modification ends
      else if (struct.eq.'hex') then
         a6=a0/2.d0
         a62=rt3*a6
         d=a62
         do ip=0,2
            m(1:3)=0
            if(ip.eq.0) then
               m(2)=1
            else if(ip.eq.1) then
               m(1)=-rt3/2.d0
               m(2)=0.5d0
            else
               m(1)=-rt3/2.d0
               m(2)=-0.5d0
            endif
            m=matmul(q2,m)
            if(abs(m(3)).lt.tol) then
               call getb(struct,ip,1,b1,a6,a62)
               do is=-1,1,2
                  nburger=nburger+1
                  normal(1:3,nburger)=m
                  burg(1:3,nburger)=is*matmul(q2,b1)
                  burg(4,nburger)=sqrt(dot_product(b1,b1))
                  write(*,'(i3,3e15.6,2x,4e15.6)') nburger,b1
     $                 ,burg(1:4,nburger)
                  call getstrain(burg(1:3,nburger),d,m,epslib(1:3,1:3
     $                 ,nburger),flib(1:3,1:3,nburger))
                  write(*,'(10x,3e15.6)') (epslib(j,1:3,nburger),j=1
     $                 ,3)
               enddo
            endif
         enddo
      endif
      nburger=nburger+1
      epslib(1:3,1:3,nburger)=0.d0
      flib(1:3,1:3,nburger)=0.d0
      do i=1,3
         flib(i,i,nburger)=1.d0
      enddo
      b1(1:3)=0.d0
      normal(1:3,nburger)=0.d0
      burg(1:4,nburger)=0.d0
      write(*,*) '------burgers vectors-------'
      do i=1,nburger
         write(*,'(i3,3e15.6,5x,4e15.6)') i,normal(1:3,i),burg(1:4,i)
         write(*,'(2(10x,3e15.6))')(epslib(j,1:3,i),flib(j,1:3,i),j=1,3)
      enddo
      write(*,*) '-----------------------'
      end subroutine rotateburgers
c*********************************************************************
c     
c     getstrain:
c     for a given dislocation find the eigenstrain and defm gradient.
c     
      subroutine getstrain(b,d,m,eps,f)
      implicit none
      double precision b(3),d,m(3),eps(3,3),f(3,3)
      integer i,j
c     
c     store dudx in f:
c     
      do i=1,3
         do j=1,3
            f(i,j)=b(i)*m(j)/d
         enddo
      enddo
c     
c     get strain:
c     
      eps=matmul(transpose(f),f)
      eps=(eps+f+transpose(f))
c     
c     add delta to f:
c     
      do i=1,3
         f(i,i)=f(i,i)+1.d0
      enddo
      end subroutine getstrain
c*********************************************************************
c     
c     getb:
c     just a lookup table of burgers vectors for different crystal
c     structures.
c     
      subroutine getb(struct,ip,ib,b,a6,a62)
      implicit none
      integer ip,ib
      double precision b(3),a6,a62
      character*3 struct
      if(struct.eq.'fcc') then
         if(ip.eq.0) then
            if(ib.eq.1) then
               b(1)=-a6
               b(2)= a62
               b(3)=-a6
            else if(ib.eq.2) then
               b(1)= a62
               b(2)=-a6
               b(3)=-a6
            else
               b(1)=-a6
               b(2)=-a6
               b(3)= a62
            endif
         else if(ip.eq.1) then
            if(ib.eq.1) then
               b(1)=-a6
               b(2)= a6
               b(3)=-a62
            else if(ib.eq.2) then
               b(1)= a62
               b(2)= a6
               b(3)= a6
            else
               b(1)=-a6
               b(2)=-a62
               b(3)= a6
            endif
         else if(ip.eq.2) then
            if(ib.eq.1) then
               b(1)= a6
               b(2)=-a6
               b(3)=-a62
            else if(ib.eq.2) then
               b(1)=-a62
               b(2)=-a6
               b(3)= a6
            else
               b(1)= a6
               b(2)= a62
               b(3)= a6
            endif
         else
            if(ib.eq.1) then
               b(1)=-a62
               b(2)= a6
               b(3)=-a6
            else if(ib.eq.2) then
               b(1)= a6
               b(2)= a6
               b(3)= a62
            else
               b(1)= a6
               b(2)=-a62
               b(3)=-a6
            endif
         endif
      else if (struct.eq.'bcc') then
         if (ip.eq.-3) then

         else if (ip.eq.-2) then
         else if (ip.eq.-1) then
         else if (ip.eq.1) then
         else if (ip.eq.2) then
         else if (ip.eq.3) then
         endif
      else if (struct.eq.'hex') then
         if(ip.eq.0) then
            b(1)=2.0*a6
            b(2)=0.d0
            b(3)=0.d0
         else if (ip.eq.1) then
            b(1)=a6
            b(2)=a62
            b(3)=0.d0
         else
            b(1)=-a6
            b(2)=a62
            b(3)=0.d0
         endif
      endif
      end subroutine getb
c*********************************************************************
c     
c     plottrigger:
c     mostly for debugging, plots the DB element that triggered a
c     dislocation pass.
c     
      subroutine plottrigger(x,b,ix,numel,numnp,nen1,ndf,nxdm,iburg,imap
     $     ,nslip,utilde,filename,nburger)
      use mod_file
      implicit none
      integer numel,numnp,nen1,ndf,nxdm,nslip,idisfile
      integer imap(nslip),ix(nen1,numel),iburg(nslip),nburger
      double precision x(nxdm,numnp),b(ndf,numnp)
     $     ,utilde(ndf,numnp)
      integer i,j,k,im,nplot
      character*80 filename
      call iofile(filename,'formatted  ',idisfile,.false.)
      nplot=0
      do i=1,nslip
         if(iburg(i).ne.nburger) nplot=nplot+1
      enddo
      write(idisfile,123) 'zone, f=fepoint, et=triangle, n=',3*nplot
     $     ,', e=',nplot
 123  format(a,i10,a,i10)
      do i=1,nslip
         do j=1,3
            if(iburg(i).ne.nburger) then
               do k=1,3
                  im=ix(k,imap(i))
                  write(idisfile,124) x(1:3,im)+b(1:3,im),iburg(i)
     $                 ,iburg(i),iburg(i)
 124              format(3e15.6,3i4)
               enddo
               go to 2
            endif
         enddo
 2       continue
      enddo
      do i=1,nplot
         write(idisfile,125) 3*i-2,3*i-1,3*i
 125     format(3i6)
      enddo
      call flush(idisfile)
      end subroutine plottrigger
c*********************************************************************
c     
c     plotdisp:
c     mostly for debugging, plots the current state of affairs in the
c     atomistic region.
c     
      subroutine plotdisp(x,b,ix,numel,numnp,nen1,ndf,nxdm,IsRelaxed
     $     ,filename)
      use mod_file
      implicit none
      integer numel,numnp,nen1,ndf,nxdm,i
      integer ix(nen1,numel),numel1,IsRelaxed(numnp),idisfile
      double precision x(nxdm,numnp),b(ndf,numnp)
      character*80 filename
      call iofile(filename,'formatted  ',idisfile,.false.)
      numel1=0
      do i=1,numel
         if(ix(nen1,i).ne.0) numel1=numel1+1
      enddo
      write(idisfile,123) 'zone, f=fepoint, et=triangle, n=',numnp
     $     ,', e=',numel1
 123  format(a,i10,a,i10)
      do i=1,numnp
         if(IsRelaxed(i).eq.0) then
            write(idisfile,127) 1000,1000,0,0,0,0
         else
            write(idisfile,124) x(1:3,i)+b(1:3,i),b(1:3,i)
         endif
 127     format(2i5,4i2)
 124     format(6e15.6)
      enddo
      do i=1,numel
         if(ix(nen1,i).ne.0) write(idisfile,125) ix(1:3,i)
 125     format(3i6)
      enddo
      call flush(idisfile)
      end subroutine plotdisp
c*********************************************************************
c     
c     InContinuum: check if the point xd is in one of the the continuum
c     region elements
c     
      logical function InContinuum(xd,ix,x,nxdm,numnp,numel,nen1,
     $     InAtoms)
      implicit none
      integer nxdm,numnp,numel,nen1
      integer ix(nen1,numel),i
      double precision xd(2),x(nxdm,numnp),s(3)
      logical intri,ontri,in,InAtoms
      InContinuum=.false.
      InAtoms=.false.
      do i=1,numel
         in=intri(x(1:2,ix(1,i)),x(1:2,ix(2,i)),x(1:2,ix(3,i)),xd,s
     $        ,ontri)
         if(in.or.ontri) then
            if(ix(nen1,i).eq.0) then
               InContinuum=.true.
            else
               InAtoms=.true.
            endif
            return
         endif
      enddo
      return
      end function InContinuum
c*********************************************************************
c     
c     FindImageLocation:  figure out where to put the image dislocation
c     so that it is far away from any continuum regions.
c     
      subroutine FindImageLocation(xi,ifactor,x0,bmat,ix,x,nxdm,numnp
     $     ,numel,nen1)
      implicit none
      integer ifactor,numel,nen1,numnp,nxdm
      integer ix(nen1,numel)
      double precision x(nxdm,numnp)
      double precision bmat(3),x0(2),xi(3),xnew(2)
      logical ic,ia,InContinuum,towardscontinuum
      integer istep
c     
c     make sure that x0 is in the atomistic region, as it should be.
c     
      ic=InContinuum(x0,ix,x,nxdm,numnp,numel,nen1,ia)
      if(.not.ia) stop 'ERROR: bug in findimagelocation'
c     
c     search for the "other end" of the slip plane, in initial steps of
c     16b, quit refining when step is down to 1b.  March away from the
c     place where the dislocation is passing across the interface until
c     you are in the continuum (or free space, which is even better),
c     then march back with smaller steps until you narrow down the
c     "other end" of the slip plane.
c     
      istep=-16
      xnew=x0(1:2)
 1    continue
      towardscontinuum=istep.lt.0
      if(abs(istep).lt.1) go to 10
 2    xnew=xnew+ifactor*istep*bmat(1:2)
      ic=InContinuum(xnew,ix,x,nxdm,numnp,numel,nen1,ia)
      if(towardscontinuum) then
         if(ic.or.(.not.ic.and..not.ia)) then
            istep=istep/(-2)
            go to 1
         endif
      else
         if(ia) then
            istep=istep/(-2)
            go to 1
         endif
      endif
      go to 2
 10   continue
c     
c     at this point, we have either found the point where the slip plane
c     re-enters the continuum or we have found a free surface
c     
c--   re-entry: put the image half-way between the two points where slip
c     plane meets continuum.
c--   free:surface: put image just outside the mesh near the free
c--   surface.
c     
      if(ic) then
         xnew=0.5d0*(xnew+x0(1:2))
      endif
      xi(1:2)=xnew
      xi(3)=0.d0
      end subroutine FindImageLocation
c*********************************************************************
c     
c     GetDetectionBand:  Takes the user-defined DB polygons (in
c     mod_boundary) and finds the elements that make up the DB.  One
c     should check "esi.tec" to make sure the DB makes sense.  The DB
c     must only hit atomistic elements (ix(4,i)=1) or this routine will
c     complain.
c     
      subroutine GetDetectionBand(ix,nen1,numel,x,numnp,nxdm,itx
     $     ,IsRelaxed)
      use mod_boundary
 1    implicit none
      integer nen1,numel,numnp,nxdm
      integer ix(nen1,numel),itx(3,numel),IsRelaxed(numnp)
      double precision x(nxdm,numnp)
      include '../Disl/disl_parameters.par'

      integer next,NMAX,ifound,iel,j,i,n1,n2,jp1,iring
      parameter(NMAX=300)
      logical found
      double precision d1,dmin,vec(2)
      integer, allocatable:: nn(:),neigh(:,:)
      integer ip(2),idb
      allocate(nn(numnp),neigh(dbnmax,numnp))
c$$$      allocate(eldb(NMAX,ndbpoly))
c$$$      allocate(dbbound(NMAX,ndbpoly))
c$$$      allocate(nelidb(ndbpoly))
c$$$      dbbound = .false.
c     
c     using ix and itx, compute the possible paths between nodes
c     
      write(*,*) ' *** computing detection band'
      nn=0
      do iel=numel,1,-1
         do j=1,3
            if(itx(j,iel).lt.iel) then
               jp1=mod(j,3)+1
               n1=ix(j,iel)
               n2=ix(jp1,iel)
               nn(n1)=nn(n1)+1
               if(nn(n1).gt.NMAX) stop 'increase NMAX'
               neigh(nn(n1),n1)=n2
               nn(n2)=nn(n2)+1
               if(nn(n2).gt.NMAX) stop 'increase NMAX'
               neigh(nn(n2),n2)=n1
            endif
         enddo
      enddo
      do idb=1, ndbpoly
         iring=0
c     
c     find nearest node to the first vertex of the dbpoly
c     
         dmin=1.e20
         do i=1,numnp
            if(IsRelaxed(i).eq.1) then
               vec=x(1:2,i)-dbpoly(1:2,1,idb)
               d1=dot_product(vec,vec)
               if(d1.lt.dmin) then
                  dmin=d1
                  ip(1)=i
               endif
            endif
         enddo
c     
c     error:
c     
         if(nn(ip(1)).eq.0) then
            write(*,*) ip(1),x(1:2,ip(1))
            do i=1,ndbvtx(idb)
               write(*,*) dbpoly(1:2,i,idb)
            enddo
            do i=1,numel
               do j=1,3
                  if(ix(j,i).eq.ip(1)) write(*,*) ix(j,i)
               enddo
            enddo
            stop 'ERROR: detection band hits atom pad'
         endif
c     
c     choose minimum path to next vertex by moving from one node to the
c     next, always choosing the path that minimizes distance to the
c     next DB polygon vertex.
c     
         next=2
 10      vec=x(1:2,ip(1))-dbpoly(1:2,next,idb)
         dmin=dot_product(vec,vec)
         found=.false.
         do i=1,nn(ip(1))
            vec=x(1:2,neigh(i,ip(1)))-dbpoly(1:2,next,idb)
            d1=dot_product(vec,vec)
            if(d1.lt.dmin) then
               found=.true.
               ifound=neigh(i,ip(1))
               dmin=d1
            endif
         enddo
         if(found) then
            ip(2)=ifound
            do iel=1,numel
               do j=1,3
                  jp1=mod(j,3)+1
                  if(ix(j,iel).eq.ip(2).and.ix(jp1,iel).eq.ip(1)) then
c$$$  print *, 'Detection band element = ', iel, -j
                     ix(nen1,iel)=-j
c     Qu modified detection band ring starts
                     iring=iring+1
                     if (iring.gt.NMAX) then
                        write(*,*)'error---NMX should be increased'
                        stop
                     endif
                     eldb(iring,idb)=iel
                     dbbound(iring,idb) = dbboundnear(idb)
c     Qu modified detection band ring ends
                     go to 20
                  endif
               enddo
            enddo
 20         ip(1)=ip(2)
            go to 10
         else 
            next=mod(next,ndbvtx(idb))+1
            if(next.ne.1)then
               go to 10
            endif
         endif
         nelidb(idb)=iring
         print *, 'DDDD', nelidb(idb), idb, iring
      enddo
C     Find the boundary detection band elements
c     Found by distance to atom boundary such that the minimum
c       DB elements at minimum distance to atom boundary are 
c       chosen to be the boundary elements
      deallocate(nn,neigh)
      end subroutine GetDetectionBand
c*********************************************************************
      subroutine lostslipinit(LostSlip)
      implicit none
      include '../Disl/disl_parameters.par'
      logical LostSlip
      LostSlip=.false.
      r_old(1:3,1:ndisl)=r_disl(1:3,1:ndisl)
      end subroutine lostslipinit
c*********************************************************************
c     
c     lostslipcheck:
c     checks for dislocations that want to leave the continuum.  Each
c     dislocation has a range in which it can live (disl_range) that
c     spans the continuum region but for a little gap near the
c     atom/continuum interface and new free surfaces.
c     
      subroutine lostslipcheck(LostSlip,ix,x,b, npass)
      implicit none
      include '../Disl/disl_parameters.par'
      double precision x,b,r, r1, r2
      integer ix,i, npass
      logical LostSlip, pass
      LostSlip=.false.
      npass = 0
!print *, 'In lostslipcheck'
      do i=1,ndisl
         if(elem_disl(i).ne.0) then
            pass = .false.
            r1=r_disl(1,i)
            r2 = r_disl(2,i)
            r = sqrt(r1**2 + r2**2)
c$$$            print *, 'In lostslipcheck', i, r1, disl_range(1,i),
c$$$     $           r2, disl_range(2,i)
            if (abs(r1) .lt. abs(disl_range(1,i))) then 
               if (abs(r2) .lt. abs(disl_range(2,i))) then 
                  pass = .true. 
               endif
               npass = npass + 1
            endif
            if (pass) then 
               print *,'Entering PasstoAtomistic',i,r, disl_range(
     $              disl_index(i),i), elem_disl(i), r_old(1:2,i),
     $              r_disl(1:2,i)
!           if(r.lt.disl_range(1,i).or.r.gt.disl_range(2,i)) then
               call PassToAtomistic(r_disl(1,i),r_old(1,i),burgers(1,i)
     $              ,theta_e(i),theta_s(i),ix,x,b,LostSlip,ndisl_dd(i))
               elem_disl(i)=0
               LostSlip = .true. 
c               exit
            endif
         endif
      enddo
      r_old(1:3,1:ndisl)=r_disl(1:3,1:ndisl)
      end subroutine lostslipcheck
c*********************************************************************
c     
c     PassToAtomistic.  As the name suggests, come here when we have
c     found a DD that wants to leave the continuum.  
c     
      subroutine PassToAtomistic(r,rold,burgers,te,ts,ix,x,b
     $     ,LostSlip,idis_slip)
      use mod_dislocation
      use mod_global
      use mod_boundary
      implicit none
      double precision r(3),rold(3),burgers(3),te,ts,vec(2),xi(3),xd(3)
     $     ,x(nxdm,numnp),s(3),b(ndf,numnp)
      integer ifactor,ix(nen1,numel),idb,iel,idis_slip
      integer, static:: NSHIFT=1
      character*80 filename
      logical PointInPoly,LostSlip,notin,intri,ontri,ntin
      do iel=1,numel
         if(intri(x(1:2,ix(1,iel)),x(1:2,ix(2,iel)),x(1:2,ix(3,iel)),r,s
     $        ,ontri)) go to 10
      enddo
c     
c     dislocation is outside of any element, assume it passed to free
c     space.
c     
      return
c     
c     figure out where to pass the dislocation
c     
 10   continue
c     
c     use last known location of the dislocation to figure out which way
c     it is moving relative to b
c     
      vec(1:2)=r(1:2)-rold(1:2)
      if(dot_product(vec(1:2),burgers(1:2)).lt.0.d0) then
         ifactor=-1
      else
         ifactor=1
      endif
      xd(1:2)=r(1:2)+ifactor*burgers(1:2)
c     
c     if disl is still in the continuum, march along by 'b' until it is
c     either in free space or the atomistics
c     
 400  if(ix(nen1,iel).ne.0) go to 450
c     
c     disl is still in the continuum - nudge it along by b:
c     
      r(1:2)=r(1:2)+ifactor*burgers(1:2)      
      do iel=1,numel
         if(intri(x(1:2,ix(1,iel)),x(1:2,ix(2,iel)),x(1:2,ix(3,iel)),r,s
     $        ,ontri)) go to 400
      enddo
c     
c     disl is in free space, return without passing to atomistics0
c     
      return
c     
c     disl is definitely moving to atomistic region.  march along until
c     you are inside a detection band polygon.
c     
 450  continue
      print *, 'dislocation in atomistic region', xd
      notin=.true.
      do idb=1,ndbpoly
c     Make sure that the dislocation is just outside the outermost
c      Detection band polygon
         if (dbboundnear(idb)) then 
            notin=notin.and.(.not.PointInPoly(xd,ndbvtx(idb),dbpoly(1,1
     $           ,idb))) 
            print *, 'Checking dislocation in ', idb
         end if
      enddo
      if(notin) then
         xd(1:2)=xd(1:2)+ifactor*burgers(1:2)
         go to 450
      endif
      xd(1:2)=xd(1:2)+NSHIFT*ifactor*burgers(1:2)
      print *, 'Disl. Position in atomistics', ifactor, xd
      ifactor=-ifactor
      call FindImageLocation(xi,ifactor,xd,burgers,ix,x,nxdm,numnp,numel
     $     ,nen1)
c     
c     puts the discrete dislocation as far from the detection bands as
c     possible
c     
      r=xi
c     
c     adjusts the atomistic displacements so that the new atomistic core
c     is just inside the detection band
c     

      print *, 'Dislocation initially at', rold(1:2)
      call disl_pass(rold,xd,burgers,te,ts,x,b,IsRelaxed,numnp
     $     ,.true.,.false., idis_slip)
      newslip=.true.
      LostSlip=.true.
      write(*,*) 'Dislocation passed to Atomistics'
      write(*,*),'With burgers vector', burgers(1:2)
      write(*,*) 'initially at:',rold(1:2)
      write(*,*) 'moved to:    ',r(1:2)
      write(*,*) 'atomistic core at:    ',xd(1:2)
      filename='out/final.plt'
      call plotdisp(x,b,ix,numel,numnp,nen1,ndf,nxdm,IsRelaxed
     $     ,filename)
      end subroutine PassToAtomistic
c*********************************************************************
c     
c     disl_restart:
c     write a restart file during minimization, just in case something
c     goes wrong.
c     
      subroutine disl_restart(id,x,ix,itx,f,b)
      use mod_global
      use mod_file
      implicit none
      include '../Disl/disl_parameters.par'
      integer id(*),ix(*),itx(*),logic
      character key*4, filename*80
      double precision x(*),f(*),b(*)
      key='writ'
      filename='disloc.res'
      call iofile(filename,'unformatted',logic,.false.)
      call rest(id,x,ix,itx,f,b,logic,key)
      close(logic)
      write(*,*) '** updating disloc.res ',ndisl
      end subroutine disl_restart
