      subroutine pmacr(id,x,ix,f,b,dr,db,itx)
      use mod_file
      use mod_global
      use mod_boundary
      use mod_dynamo
      use mod_grain
      use mod_material
      implicit none
c
      double precision b(ndf,*),x(nxdm,*),f(*),dr(ndf,*),db(*)
      integer id(*),ix(nen1,*)
      integer itx(3,*)

c
c---- macro instruction interpreter
c---- controls solution algorithms
c
      integer nwd
      parameter (nwd=56)
      logical pcomp,check,boundary,fladapt,statrue
      character*80 input,filename, ctemp
      character*4 key
      real*4 wd(nwd)
      integer icrit,isign,nodecrit,idfcrit
      double precision arcrate
      common /cstep/ boundary,icrit,arcrate,isign,nodecrit,idfcrit

      integer lvs(9),lve(9),neqmax,i,j,k,l,lv,lx,ll,lmax,ii,lower,upper
     $     ,next,ldif,i1,i2,loops,idum,n1,n2,n,logic,l0,numelsve
     $     ,loopshow,nloop,ngt0
      double precision dum,bnew,varstp,dtstep,tmax,dttol,dt0,qrot(3,3)
     $     ,cc2(6,6),StrainEnergy, dtstp_orig
      logical flag02
      real*4, pointer:: ct(:,:),ctpass(:,:)
      integer, allocatable :: itemp1(:),itemp2(:)
      data wd/'tole','dtim','dump','stre','xxxx','xxxx','loop','next',
     1        'xxxx','chec','time','xxxx','chti','xxxx','xxxx','mesh',
     2        'xxxx','end ','pdel','xxxx','xxxx','xxxx','prec','ddse',
     3        'xxxx','newd','ma05','ma06','getc','xxxx','xxxx','xxxx',
     4        'xxxx','xxxx','xxxx','rest','xxxx','ma01','ma02','ma03',
     5        'xxxx','xxxx','clea','xxxx','stat','xxxx','xxxx','xxxx',
     6        'xxxx','xxxx','xxxx','xxxx','xxxx','xxxx','xxxx','xxxx'/
c
c---- set initial values of parameters
      neqmax = maxnp*ndf
      boundary  = .true.
      icrit     = 0
      niter     = 0
      rnmax     = 0.d0
      timeol      = 0.d0
      time      = 0.d0
      flag02=.true.
      CONVERGETOL = 1.d-9
c---- read macro program
      write(6,*)
      write(6,'('' '',a80//5x,''macro instructions''//1x,
     1 ''nesting  statement'',10x,''variables'')') head
      write(6,'(1x,''-------  ---------'',10x,''---------'')')
      nloop=0
      loopshow=0
      ll = 1
      lmax = 16
      allocate(ct(20,lmax),ctpass(20,lmax))
      ct(1,1) = wd(7)
      ct(3,1) = 1.0
100   ll = ll + 1
      if (ll.lt.lmax) go to 110
      lmax = lmax + 16
      ctpass=ct
      deallocate(ct)
      allocate(ct(20,lmax))
      ct(1:20,1:lmax)=ctpass
      deallocate(ctpass)
      allocate(ctpass(20,lmax))
110   continue
      read(5,'(a80)') input
c---- added to allow indentation with spaces in input file
      i=1
      do while ((input(i:i).eq.' ').and.(i.le.80))
         i=i+1
      enddo
      if(i.le.80) then
         if (input(i:i).eq.'%') go to 110
      endif
      upper = i-1
      do 115 ii = 1,2
      lower = upper
      upper = next(lower,input)
      ldif = upper - lower - 1
      if (ldif.eq.0) then
      ct(ii,ll) = '    '
      else
      call equal(ct(ii,ll),input,lower+1,upper-1)
      end if
115   continue
      if (pcomp(ct(1,ll),wd(7))) then
c-- added to show loop nesting
         nloop=nloop+1
         loopshow=loopshow+nloop*10**(nloop-1)
         i1 = upper
         i2 = next(i1,input)
         call freein(input,i1,i2,loops,dum,1)
         ct(3,ll) = loops
      else
      call equal(ct(3,ll),input,upper+1,80)
      end if
      if((loopshow.gt.0)) then
         write(6,'(i8,2x,a4,1x,a4,9x,a40)') loopshow,(ct(j,ll),j=1,2)
     $        ,input(upper+1:)
      else
         write(6,'(10x,a4,1x,a4,9x,a40)') (ct(j,ll),j=1,2)
     $        ,input(upper+1:)
      endif
c--- added to show loop nesting
      if (pcomp(ct(1,ll),wd(8))) then
         loopshow=loopshow-nloop*10**(nloop-1)
         nloop=nloop-1
      endif
      if (.not.pcomp(ct(1,ll),wd(18))) go to 100
200   ct(1,ll) = wd(8)
c---- set loop markers
      lx = ll - 1
      do 230 l = 1,lx
      if (.not.pcomp(ct(1,l),wd(7))) go to 230
      j = 1
      k = l + 1
      do 210 i = k,ll
      if (pcomp(ct(1,i),wd(7))) j = j + 1
      if (j.gt.9) then
      write(6,'(/'' **error** loops nested deeper than 8'')')
      stop
      end if
      if (pcomp(ct(1,i),wd(8))) j = j - 1
210   if (j.eq.0) go to 220
      write(6,'(/'' **error** unbalanced loop/next macros'')')
      stop
220   ct(4,i) = l
      ct(4,l) = i
230   continue
      j = 0
      do 240 l = 1,ll
      if (pcomp(ct(1,l),wd(7))) j = j + 1
240   if (pcomp(ct(1,l),wd(8))) j = j - 1
      if (j.ne.0) then
      write(6,'(/'' **error** unbalanced loop/next macros'')')
      stop
      end if
c---- execute macro instruction program
      lv = 0
      l = 1
299   do 300 j = 1,nwd
300   if (pcomp(ct(1,l),wd(j))) go to 310
      go to 330
310   i = l - 1
      if (l.ne.1.and.l.ne.ll) then
          write(6,'(2x,''**macro instruction'',i4,'' ** '',
     1     2(a4,2x))') i,(ct(k,l),k = 1,2)
          call flush(6)
      endif
      if (j.le.24) then
      go to (1,2,3,4,5,6,7,8,9,10,11,12,
     1 13,14,15,16,17,18,19,20,21,22,23,24),j
      else if(j.le.48) then
      go to (25,26,27,28,29,31,32,33,34,35,36,37,
     1 38,39,41,42,43,44,45,46,47,48,49,50),(j-24)
      else
      go to (51,52,53,54,55,56,57,58),(j-48)
      end if
c
c---- macro 'tole'
c---- set solution tolerance
1     continue
      lower = 0
      upper = next(lower,ct(3,l))
      call freein(ct(3,l),lower,upper,idum,CONVERGETOL,2)
      write(6,'(2x,''**tolerance set to '',e15.5)') CONVERGETOL
      go to 330
c
c---- macro 'dtim'
c---- set time increment
2     continue
      lower = 0
      upper = next(lower,ct(3,l))
      call freein(ct(3,l),lower,upper,idum,dt,2)
      write(6,'(2x,''**time step set to '',e15.5)') dt
      go to 330
c
c---- macro 'dump'
3     continue
      call leodump(numnp,numel,ndf,nxdm,nen1,x,ix,id,IsRelaxed,f,itx)
      go to 330
c
c---- macro 'stre'
4     continue
      call hstress(numnp,IsRelaxed,avevirst)
      go to 330
c
c---- macro 'xxxx'
5     continue
      go to 330
c
c---- macro 'xxxx'
6     continue
      go to 330
c
c---- macro 'loop'
c---- set loop start indicators
7     lv = lv + 1
      lx = ct(4,l)
      lvs(lv) = l
      lve(lv) = lx
      ct(3,lx) = 1.
      go to 330
c
c---- macro 'next'
c---- loop terminator control
8     n = ct(4,l)
      ct(3,l) = ct(3,l) + 1.0
      if(ct(3,l).gt.ct(3,n)) lv = lv - 1
      if(ct(3,l).le.ct(3,n)) l = n
      go to 330
c
c---- macro 'xxxx'
9     continue
      go to 330
c
c---- macro 'chec'
c
   10 continue
      write(*,*) " chec no longer in this version of CADD "
      stop       
c      call derivcheck(id,b,x,ix,f,dr,itx)
      go to 330
c
c---- macro 'time'
c---- increment time
11    continue       
      timeol=time
      upper = 0
      dtstep = dt         
      if(dtstep.ne.0.d0) b0(1:ndf,1:numnp)=b(1:ndf,1:numnp)
      if (Moved) then 
         dtstep = 0.0d0
         if (xtip_init(1) < 0.0d0) then 
            call get_crack_tip(x,b)
            xtip_init(1:2)=xtip(1:2)
         end if
         print *, 'Moving mesh, repeat calculation without load change'
         Moved = .false. 
c     Reset moved flag to .false.
      end if
      time = time + dtstep
      db(1:maxnp*ndf)=0.d0
      write(6,'(2x,''**time set to '',e15.5)') time
      lower = upper
      upper = next(lower,ct(3,l))
      call freein(ct(3,l),lower,upper,idum,tmax,2)
      if (tmax.gt.0.0) then
      if (time.gt.tmax) then
      write(6,'(2x,''**maximum time attained**'')')
      write(6,'(2x,''**end of macro execution**'')')
      return
      end if
      end if
      rnmax = 0.0d0
      niter = 0
      go to 330
c
c---- macro 'xxxx'
   12 continue
      go to 330
c
c---- macro 'chti'
   13 continue
      call checktime
      go to 330
c
c---- macro 'xxxx'
14    continue
      go to 330
c
c---- macro 'xxxx'
15    continue
      go to 330
      
c
c---- macro 'mesh'
c---- modify mesh data
16    call pmesh(id,x,ix,f,b,dr,itx)
      go to 330
      
c
c---- macro 'xxxx'
c----
 17   continue
      go to 330
c
c---- macro 'xxxx'
   18 continue
      go to 330
c
c---- macro 'pdel'
   19 continue
      lower = 0
      upper = next(lower,ct(3,l))
      call equal(filename,ct(3,l),lower+1,upper-1)
      call iofile(filename,'formatted  ',logic,.true.)
      call pdel(f,dr,id,logic,ndf,x,time)
      go to 330
c
c---- macro 'xxxx'
20    continue
      go to 330
c
c---- macro 'xxxx'
21    continue
      go to 330
c
c---- macro 'xxxx'
22    continue
      go to 330
c
c---- macro 'prec'
c---- initializes the displacement to the K-field
23    continue
      call precrack(id, x,b,f, ct(2,l))
      go to 330
c
c---- macro 'ddse'
c---- set up the discrete dislocations solving routine
24    continue
      if(DD_set) return
      write(*,*)
      write(*,*) 'setting up the discrete dislocation solver'
      write(*,*)
      call bandnl(id,x,ix,f,b)
      qrot=0.d0
      if(ngrains.ne.1) stop 'hardwired here for 1 grain'
      qrot(1,1)=grains(1)%rotmat(1)
      qrot(2,2)=grains(1)%rotmat(1)
      qrot(1,2)=-grains(1)%rotmat(2)
      qrot(2,1)=grains(1)%rotmat(2)
      qrot(3,3)=1.d0
      qrot=matmul(qrot,grains(1)%xlatvect(1:3,1:3))
      call rotatevoigt(qrot,material(1)%cc,cc2)
      write(*,*)
      write(*,*) 'Rotated elastic matrix:'
      do i=1,6
         write(*,'(6e13.4)') cc2(i,1:6)
      enddo
      allocate(itemp1(numnp),itemp2(numel))
      call fem_setup(numnp, numel, x, id, IsRelaxed, ix, itx, cc2,
     $     itemp1, itemp2)

      deallocate(itemp1,itemp2)
      call disl_setup
      DD_set=.true.
      go to 330
c
c---- macro 'xxxx'
25    continue
      go to 330
c
c---- macro 'newd'
c---  add a new discrete dislocation
26    continue
      print *, 'Entering new dislocation'
      call NewDislocation(ct(2,l),x,b,IsRelaxed,numnp,nxdm,ndf)
      go to 330
c
c---- macro 'ma05'
c---- user supplied macro
27    continue
      if(.not.DD_set) stop 'ERROR: must initialize D.D.'
      ngt0=ngtlst

      print*,'did you mean ma06?'
      stop

      write(*,*) 'NUMBER NEIGHBOR UPDATES:',ngtlst-ngt0
      go to 330
c
c---- macro 'xxxx'
c---- user supplied macro to call MD/DD routine
28    continue

      if(.not.DD_set) stop 'ERROR: must initialize D.D.'
      ngt0=ngtlst

      call ma06(id,x,ix,f,b,dr,db,ct(2,l),itx)
      write(*,*) 'NUMBER NEIGHBOR UPDATES:',ngtlst-ngt0
      call get_crack_tip(x,b)
c$$$      if (abs(xtip(1)-xtip_init(1)) > x_move_mesh) then 
c$$$         do j = 1, 2
c$$$            xtip_actual(j)=abs(xtip(j)-xtip_init(j))
c$$$            if (xtip(j) > xtip_init(j)) then 
c$$$               x_tip_dir(j) = 1.0d0
c$$$            else
c$$$               x_tip_dir(j) = -1.0d0
c$$$            end if
c$$$         end do
c$$$         call move_atomistic_crack(x,ix,b)
c$$$      end if

      go to 330
c
c---- macro 'getc'
c---- macro to set the initial crack tip position
 29   continue
      call get_crack_tip(x,b)
      xtip_init(1:2)=xtip(1:2)
      crack_motion = 0.0
      MoveMesh=.true.
      go to 330
c
c---- macro 'xxxx'
31    continue
      go to 330
c
c---- macro 'xxxx'
   32 continue
      go to 330
c
c---- macro 'xxxx'
   33 continue
      go to 330
c
c---- macro 'xxxx'
 34   continue
      go to 330
c
c---- macro 'xxxx'
   35 continue
      go to 330
c
c---- macro 'xxxx'
   36 continue
      go to 330
c
c---- macro 'rest'
c---- read/write restart files
   37 continue
      write(*,*) " (rest macro) not in this version! "
      stop
c
c---- macro 'xxxx'
   38 continue
      go to 330
c
c---- macro 'ma01'
c---- user supplied macro - boundary conditions/loading
   39 call ma01(id,x,ix,f,b,dr,db,ct(2,l))
      go to 330
c
c---- macro 'ma02'
c---- user supplied macro - tecplot files
   41 continue
       call ma02(id,x,ix,f,b,dr,db,ct(2,l),flag02)
      go to 330
c
c---- macro 'ma03'
c---- user supplied macro - energy calculation
   42 continue
      write(*,*) "ma03 not in this version of CADD"
      stop

      go to 330
c
c---- macro 'xxxx'
   43 continue
      go to 330
c
c---- macro 'xxxx'
   44 continue
      go to 330
c
c---- macro 'clea'
c---- Reinitialize variables within same time step for a different
c---- solution procedure.
   45 continue
      b(1:ndf,1:numnp)=b0(1:ndf,1:numnp)
      db(1:maxnp*ndf)=0.d0
      write(6,*) ' **Displacement increments initialized'
      rnmax = 0.0d0
      niter = 0
      go to 330
c
c---- macro 'xxxx'
   46 continue
      go to 330
c
c---- macro 'stat'
c---- Recompute nonlocal element status
  47  continue
      call latticecheck(x)
      call StatusCalc(x,ix,.false.)
      go to 330
c
c---- macro 'xxxx'
   48 continue
      go to 330
c
c---- macro 'xxxx'
   49 continue
      go to 330
c
c---- macro 'xxxx'
   50 continue
      go to 330
c
c---- macro 'xxxx'
 51   continue
      go to 330
c
c---- macro 'xxxx'
 52   continue
      go to 330
c
c---- macro 'xxxx'
 53   continue
      go to 330
c
c---- macro 'xxxx'
 54   continue
      go to 330
c
c---- macro 'xxxx'
 55   continue
      go to 330
c
c---- macro 'xxxx'
 56   continue
      go to 330
c
c---- macro 'xxxx'
 57   continue
      go to 330
c
c---- macro 'xxxx'
 58   continue
      go to 330
c
330   l = l + 1
      if (l.le.ll) go to 299
      write(6,'(2x,''**end of macro execution**'')')
      return
      end


      subroutine equal(char1,char2,l1,l2)
      implicit double precision (a-h,o-z)
c
      character*80 char1
      character*80 char2
c
      char1 = char2(l1:l2)
      return
      end


      subroutine keyequal(char1,char2,l1,l2)
      implicit double precision (a-h,o-z)
c
      character*4 char1
      character*80 char2
c
      char1 = char2(l1:l2)
      return
      end



************************************************************************
      subroutine rotatevoigt(q,c,c2)
      implicit none
      double precision q(3,3),c(6,6),c2(6,6),cc(3,3,3,3),cc2(3,3,3,3)
      integer i1,i2,i3,i4,j1,j2,j3,j4
      call unvoigt(c,cc)
      do i1=1,3
         do i2=1,3
            do i3=1,3
               do i4=1,3
                  cc2(i1,i2,i3,i4)=0.
                  do j1=1,3
                     do j2=1,3
                        do j3=1,3
                           do j4=1,3
                              cc2(i1,i2,i3,i4)=cc2(i1,i2,i3,i4) + q(i1
     $                             ,j1)*q(i2,j2)*q(i3,j3)*q(i4,j4)*cc(j1
     $                             ,j2,j3,j4)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      call voigt(cc2,c2)
      end subroutine rotatevoigt
************************************************************************
      subroutine voigt(c,cv)
      implicit none
      double precision c(3,3,3,3),cv(6,6)
      integer i,j
      cv(1,1)=c(1,1,1,1)
      cv(1,2)=c(1,1,2,2)
      cv(1,3)=c(1,1,3,3)
      cv(1,4)=c(1,1,2,3)
      cv(1,5)=c(1,1,1,3)
      cv(1,6)=c(1,1,1,2)
      cv(2,2)=c(2,2,2,2)
      cv(2,3)=c(2,2,3,3)
      cv(2,4)=c(2,2,2,3)
      cv(2,5)=c(2,2,1,3)
      cv(2,6)=c(2,2,1,2)
      cv(3,3)=c(3,3,3,3)
      cv(3,4)=c(3,3,2,3)
      cv(3,5)=c(3,3,1,3)
      cv(3,6)=c(3,3,1,2)
      cv(4,4)=c(2,3,2,3)
      cv(4,5)=c(2,3,1,3)
      cv(4,6)=c(2,3,1,2)
      cv(5,5)=c(1,3,1,3)
      cv(5,6)=c(1,3,1,2)
      cv(6,6)=c(1,2,1,2)
      do i=1,6
         do j=i,6
            cv(j,i)=cv(i,j)
         enddo
      enddo
      end subroutine voigt
************************************************************************
      subroutine unvoigt(cv,c)
      implicit none
      double precision c(3,3,3,3),cv(6,6)
      integer i1,i2,j1,j2,j3,j4
      do i1=1,6
         if (i1.eq.1) then
            j1=1
            j2=1
         else if (i1.eq.2) then
            j1=2
            j2=2
         else if (i1.eq.3) then
            j1=3
            j2=3
         else if (i1.eq.4) then
            j1=2
            j2=3
         else if (i1.eq.5) then
            j1=1
            j2=3
         else if (i1.eq.6) then
            j1=1
            j2=2
         endif
         do i2=i1,6
            if (i2.eq.1) then
               j3=1
               j4=1
            else if (i2.eq.2) then
               j3=2
               j4=2
            else if (i2.eq.3) then
               j3=3
               j4=3
            else if (i2.eq.4) then
               j3=2
               j4=3
            else if (i2.eq.5) then
               j3=1
               j4=3
            else if (i2.eq.6) then
               j3=1
               j4=2
            endif
            c(j1,j2,j3,j4)=cv(i1,i2)
            c(j2,j1,j3,j4)=cv(i1,i2)
            c(j1,j2,j4,j3)=cv(i1,i2)
            c(j2,j1,j4,j3)=cv(i1,i2)
            c(j3,j4,j1,j2)=cv(i1,i2)
            c(j3,j4,j2,j1)=cv(i1,i2)
            c(j4,j3,j1,j2)=cv(i1,i2)
            c(j4,j3,j2,j1)=cv(i1,i2)
         enddo
      enddo
      end subroutine unvoigt
************************************************************************
      subroutine leodump(numnp,numel,ndf,nxdm,nen1,x,ix,id,IsRelaxed,f
     $     ,itx)
      implicit none
      integer numnp,numel,ndf,nxdm,nen1,i
      integer ix(nen1,numel),id(ndf,numnp),IsRelaxed(numnp),itx(3,numel)
      double precision x(nxdm,numnp),b(ndf,numnp),f(ndf,numnp)
      write(*,*) numnp,numel
      write(*,*) 'x:'
      do i=1,numnp
         write(*,*) x(1:3,i)
      enddo
      write(*,*) 'f:'
      do i=1,numnp
         write(*,*) f(1:3,i)
      enddo
      write(*,*) 'id:'
      do i=1,numnp
         write(*,*) id(1:3,i)
      enddo
      write(*,*) 'IsRelaxed:'
      do i=1,numnp
         write(*,*) IsRelaxed(i)
      enddo
      write(*,*) 'ix:'
      do i=1,numel
         write(*,*) ix(1:4,i)
      enddo
      write(*,*) 'itx:'
      do i=1,numel
         write(*,*) itx(1:3,i)
      enddo


      end subroutine leodump


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine  hstress(numnp,IsRelaxed,virst)
      implicit none
      integer i, j, k, natoms
      integer numnp
      integer IsRelaxed(numnp)
      double precision virst(3,3,numnp)
      double precision hydro_stress
      double precision avg_stress(3,3)
      double precision tstress

      hydro_stress=0.0
      avg_stress(1:3,1:3)=0.0
      natoms=0
      do i=1,numnp
        if (IsRelaxed(i).eq.1) then
          natoms=natoms+1
          do j=1,3
            hydro_stress=hydro_stress+virst(j,j,i) 
            do k=1,3
              avg_stress(j,k)=avg_stress(j,k)+virst(j,k,i)
            enddo
          enddo
        endif
      enddo

      avg_stress(1:3,1:3)=avg_stress(1:3,1:3)/natoms
      hydro_stress=hydro_stress*1.0/3.0/natoms
      print*,' '
      print*,' **** Average Stresses in Atomistic Region (eV/A^3) ****'
c      print*, avg_stress(1:3,1:3) 
      print*,'hydrostatic stress = ',hydro_stress
      tstress=0.5*(avg_stress(1,1)+avg_stress(3,3))
      print*,'s(1,1) s(2,2) s(3,3) '
      print*,'s(2,3) s(1,3) s(1,2) '
      print*,'strse=',avg_stress(1,1),avg_stress(2,2),avg_stress(3,3)
      print*,'strss=',avg_stress(2,3),avg_stress(1,3),avg_stress(1,2)
      print*,' '
      
      end subroutine hstress
