
      subroutine bandnl(id,x,ix,f,b)
      use mod_global
      implicit none
c
c---- optimize  bandwidth/profile
c
      integer ngra,idpth,ideg
      common /gra/ ngra,idpth,ideg
      double precision  x(nxdm,*),f(ndf,*),b(ndf,*)
      integer ix(nen1,*),id(ndf,*)
      integer, pointer:: ndstk(:,:),ndeg(:),iold(:),lvl(:),renum(:)
     $     ,ccstor(:),lvls1(:),lvls2(:)
      logical deglimit
      integer maxdeg,i,ibw2,ipf2

      print *,'** Bandwidth Optimization **'

      maxdeg=600
      allocate(ndstk(maxdeg,numnp),ndeg(numnp),iold(numnp)
     $     ,renum(numnp+1),lvl(numnp),lvls1(numnp),lvls2(numnp)
     $     ,ccstor(numnp))
      ngra  = numnp
c
      renum=0
      lvl=0
      lvls1=0
      lvls2=0
      ccstor=0
      do i = 1,numnp
         iold(i)=i
      end do
c
 1    call setcongr(ix,ndstk,ndeg,maxdeg,deglimit,b,x)
      if(deglimit) then
         maxdeg=maxdeg*2
         deallocate(ndstk)
         allocate(ndstk(maxdeg,numnp))
         write(*,*) '**WARNING: maxdeg increased to ',maxdeg
         write(*,*) '           in bandwidth optimization'
         go to 1
      endif
      call reduce(ndstk,maxdeg,iold,renum,ndeg,lvl,lvls1,lvls2,ccstor
     $     ,ibw2,ipf2)
      call swapallgr(renum,id,x,ix,f,b)
      deallocate(ndstk,ndeg,iold,renum,lvl,lvls1,lvls2,ccstor)
      return
      end

c---VBS------------------------------
c-------------------
      subroutine setcongr(ix,ndstk,ndeg,maxdeg,deglimit,b,x)
      use mod_global
      use mod_dynamo
      implicit none
c
      integer maxdeg
      integer  ix(nen1,*),ndstk(maxdeg,*),ndeg(*)
      double precision b(ndf,*),x(nxdm,*)
      logical deglimit,NeedList
c
      integer ngra,idpth,ideg
      common /gra/ ngra,idpth,ideg
      integer  nloclist(40)
      integer node,j,irep,jpoint,nnlist,ilist,jlist
     &     ,ii,l,iel,natms,igrain,i
      deglimit=.false.

c
      do node = 1,numnp
         ndeg(node) = 0
         do j = 1, maxdeg
            ndstk(j,node) = 0
         enddo
      enddo
c
      !Now Handle Local Elements

      do iel = 1,numel
         nnlist = 3
         do ilist = 1,nnlist
            node = ix(ilist,iel)
            do jlist = 1,nnlist
               ii = ix(jlist,iel)
               if (ii.eq.node) goto 41
               if (ndeg(node).gt.0) then
                  do l = 1,ndeg(node)
                     if (ndstk(l,node).eq.ii) goto 41
                  enddo
               end if
               if (ndeg(node).lt.maxdeg) then
                  ndeg(node) = ndeg(node) + 1
                  ndstk(ndeg(node),node) = ii
               else
                  deglimit=.true.
                  return
               end if
 41            continue
            enddo
         enddo
      enddo

      ideg = 0
      do i = 1,numnp
         if (ndeg(i).gt.ideg) ideg = ndeg(i)
         ndeg(i) = 0
      enddo
      print *,'Maximum node degree = ',ideg
      return
      end


c-------------------VBS---------------------
c-------------------------------------------
      subroutine swapallgr(num,id,x,ix,f,b)
      use mod_grain
      use mod_global
      use mod_dynamo
      use mod_boundary
      implicit none
c
      integer num(*),id(ndf,*),ix(nen1,*)
      double precision x(nxdm,*),f(ndf,*),b(ndf,*)
c
      integer i,j,irep,jpoint,nnlist,jlist,ii
      integer, dimension(:), pointer :: itemp
      double precision, dimension(:), pointer :: temp

c
c  force a new neighbor list
c
      newlst=1
c
c     allocate
c
      allocate(temp(numnp),itemp(numnp))
c
      do j = 1,ndf
         do i = 1,numnp
            itemp(num(i)) = id(j,i)
         end do
         do i = 1, numnp
            id(j,i) = itemp(i)
         end do
      end do
c
      do j = 1,nxdm
         do i = 1,numnp
            temp(num(i)) = x(j,i)
         end do
         do i = 1, numnp
            x(j,i) = temp(i)
         end do
      end do
c
      do j = 1,ndf
         do i = 1,numnp
            temp(num(i)) = f(j,i)
         end do
         do i = 1, numnp
            f(j,i) = temp(i)
         end do
      end do
c
      do j = 1,ndf
         do i = 1,numnp
            temp(num(i)) = b(j,i)
         end do
         do i = 1, numnp
            b(j,i) = temp(i)
         end do
      end do
c
      do j = 1,nen
         do i = 1,numel
            ix(j,i) = num(ix(j,i))
         end do
      end do
c
c
c for constrained delaunay
c
      if (nce.gt.0) then
         do i=1,nce
            do j=1,2
               ii = elist(j,i)
               elist(j,i) = num(ii)
            enddo
         enddo
      endif

      do i = 1,numnp
         itemp(num(i)) = IsRelaxed(i)
      enddo
      do i = 1,numnp
         IsRelaxed(i)= itemp(i)
      enddo
      do i=1,numnp
      itemp(num(i))= atomSpecie(i)
      enddo
      do i=1,numnp
      atomSpecie(i)=itemp(i)
      enddo
c
      do i = 1,numnp
         temp(num(i)) = energy(i)
      enddo
      do i = 1,numnp
         energy(i)= temp(i)
      enddo
c
      deallocate(temp,itemp)
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dgree(ndstk, nr, ndeg, iold, ibw1, ipf1)
c
c
c  dgree computes the degree of each node in ndstk and stores
c  it in the array ndeg.  The bandwidth and profile for the original
c  or input renumbering of the graph is computed also.
c  Use integer*2 ndstk  with an ibm 360 or 370.
c
c
      integer ndstk
c
      common /gra/ n, idpth, ideg
c
      dimension ndstk(nr,*), ndeg(*), iold(*)
c
      ibw1 = 0
      ipf1 = 0
c
      do 40 i = 1, n
        ndeg(i) = 0
        irw     = 0
        do 20 j = 1, ideg
          itst = ndstk(j,i)
          if (itst) 30, 30, 10
   10     ndeg(i) = ndeg(i) + 1
          idif    = iold(i) - iold(itst)
          if (irw.lt.idif) irw = idif
   20   continue
   30   ipf1 = ipf1 + irw
        if (irw.gt.ibw1) ibw1 = irw
   40 continue
c
c
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fndiam(snd1, snd2, ndstk, nr, ndeg, lvl, lvls1,
     * lvls2, iwk, idflt)
c
c
c  fndiam is the control procedure for finding the pseudo-diameter of
c  ndstk as well as the level structure from each end
c  snd1-        on input this is the node number of the first
c               attempt at finding a diameter.  on output it
c               contains the actual number used.
c  snd2-        on output contains other end of diameter
c  lvls1-       array containing level structure with snd1 as root
c  lvls2-       array containing level structure with snd2 as root
c  idflt-       flag used in picking final level structure, set
c               =1 if width of lvls1 .le. width of lvls2, otherwise =2
c  lvl,iwk-     working storage
c  Use integer*2 ndstk  with an ibm 360 or 370.
c
      integer ndstk
      integer flag, snd, snd1, snd2
c
      common /gra/ n, idpth, ideg
c
c  It is assumed that the last level has at most 'maxlvl' nodes.
c
      integer maxlvl
      parameter (maxlvl=2000)
      common / cc / ndlst(maxlvl)
      dimension ndstk(nr,*), ndeg(*), lvl(*), lvls1(*), lvls2(*),
     +          iwk(*)
c
      flag = 0
      mtw2 = n
      snd  = snd1
c
c  Zero lvl to indicate all nodes are available to tree.
c
   10 do 20 i = 1, n
        lvl(i) = 0
   20 continue
      lvln = 1
c
c  Drop a tree from snd.
c
      call tree(snd, ndstk, nr, lvl, iwk, ndeg, lvlwth, lvlbot,
     *          lvln, maxlw, mtw2)
      if (flag.ge.1) go to 50
c
      flag  = 1
   30 idpth = lvln - 1
      mtw1  = maxlw
c
c  Copy level structure into lvls1.
c
      do 40 i = 1, n
        lvls1(i) = lvl(i)
   40 continue
      ndxn = 1
      ndxl = 0
      mtw2 = n
c
c  Sort last level by degree  and store in ndlst.
c
      call sortdg(ndlst, iwk(lvlbot), ndxl, lvlwth, ndeg)
      snd = ndlst(1)
      go to 10
   50 if (idpth.ge.lvln-1) go to 60
c
c  Start again with new starting node.
c
      snd1 = snd
      go to 30
   60 if (maxlw.ge.mtw2) go to 80
      mtw2 = maxlw
      snd2 = snd
c
c  Store narrowest reverse level structure in lvls2.
c
      do 70 i = 1, n
        lvls2(i) = lvl(i)
   70 continue

   80 if (ndxn.eq.ndxl) go to 90
c
c  Try next node in ndlst.
c
      ndxn = ndxn + 1
      if (ndxn.gt.maxlvl) then
         print *,'***ERROR: Insufficient storage in fndiam for ndlst()'
         stop
      endif
      snd  = ndlst(ndxn)
      go to 10
   90 idflt = 1
      if (mtw2.le.mtw1) idflt = 2
c
c
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine number(snd, num, ndstk, lvls2, ndeg, renum, lvlst,
     +                  lstpt, nr, nflg, ibw2, ipf2, ipfa, isdir)
c
c  Number produces the numbering of the graph for min bandwidth
c  snd-         on input the node to begin numbering on
c  num-         on input and output, the next available number
c  lvls2-       the level structure to be used in numbering
c  renum-       the array used to store the new numbering
c  lvlst-       on output contains level structure
c  lstpt(i)-    on output, index into lvlst to first node in ith lvl
c               lstpt(i+1) - lstpt(i) = number of nodes in ith lvl
c  nflg-        =+1 if snd is forward end of pseudo-diam
c               =-1 if snd is reverse end of pseudo-diam
c  ibw2-        bandwidth of new numbering computed by number
c  ipf2-        profile of new numbering computed by number
c  ipfa-        working storage used to compute profile and bandwidth
c  isdir-       indicates step direction used in numbering(+1 or -1)
c  Use integer*2 ndstk  with an ibm 360 or 370.
c
      integer ndstk
      integer snd, stka, stkb, stkc, stkd, xa, xb, xc, xd, cx, end,
     +        renum, test
c
      common /gra/ n, idpth, ideg
c
c  The storage in common blocks cc and lvlw is now free and can
c  be used for stacks.
c
      integer maxlvl
      parameter (maxlvl=2000)
      common / lvlw /   stka(maxlvl), stkb(maxlvl), stkc(maxlvl)
      common / cc   /   stkd(maxlvl)
c
      dimension ipfa(*)
      dimension ndstk(nr,*), lvls2(*), ndeg(*), renum(*),
     +          lvlst(*),    lstpt(*)
c
c  Set up lvlst and lstpt from lvls2.
c
      do 10 i = 1, n
        ipfa(i) = 0
   10 continue
      nstpt = 1
      do 30 i = 1, idpth
        lstpt(i) = nstpt
        do 20 j = 1, n
          if (lvls2(j).ne.i) go to 20
          lvlst(nstpt) = j
          nstpt = nstpt + 1
   20   continue
   30 continue
      lstpt(idpth+1) = nstpt
c
c  stka, stkb, stkc and stkd are stacks with pointers
c  xa,xb,xc, and xd.  cx is a special pointer into stkc which
c  indicates the particular node being processed.
c  lvln keeps track of the level we are working at.
c  initially stkc contains only the initial node, snd.
c
      lvln = 0
      if (nflg.lt.0) lvln = idpth + 1
      xc = 1
      stkc(xc) = snd
   40 cx = 1
      xd = 0
      lvln = lvln + nflg
      lst  = lstpt(lvln)
      lnd  = lstpt(lvln+1) - 1
c
c  Begin processing node stkc(cx).
c
   50 ipro        = stkc(cx)
      renum(ipro) = num
      num         = num + isdir
      end         = ndeg(ipro)
      xa          = 0
      xb          = 0
c
c  Check all adjacent nodes.
c
      do 80 i = 1, end
        test = ndstk(i,ipro)
        inx  = renum(test)
c
c  Only nodes not numbered or already on a stack are added.
c
        if (inx.eq.0) go to 60
        if (inx.lt.0) go to 80
c
c  Do preliminary bandwidth and profile calculations.
c
        nbw = (renum(ipro)-inx)*isdir
        if (isdir.gt.0) inx = renum(ipro)
        if (ipfa(inx).lt.nbw) ipfa(inx) = nbw
        go to 80
   60   renum(test) = -1
c
c  Put nodes on same level on stka, all others on stkb.
c
        if (lvls2(test).eq.lvls2(ipro)) go to 70
        xb = xb + 1
        if (xb.gt.maxlvl) then
           print *,'***ERROR: Insufficient storage in subroutine number'
           print *,'           for stkb()'
           stop
        endif
        stkb(xb) = test
        go to 80
   70   xa = xa + 1
        if (xa.gt.maxlvl) then
           print *,'***ERROR: Insufficient storage in subroutine number'
           print *,'           for stka()'
           stop
        endif
        stka(xa) = test
   80 continue
c
c  Sort stka and stkb into increasing degree and add stka to stkc
c  and stkb to stkd.
c
      if (xa.eq.0) go to 100
      if (xa.eq.1) go to 90
      call sortdg(stkc, stka, xc, xa, ndeg)
      go to 100
   90 xc = xc + 1
      if (xc.gt.maxlvl) then
         print *,'***ERROR: Insufficient storage in subroutine number'
         print *,'           for stkc()'
         stop
      endif
      stkc(xc) = stka(xa)
  100 if (xb.eq.0) go to 120
      if (xb.eq.1) go to 110
      call sortdg(stkd, stkb, xd, xb, ndeg)
      go to 120
  110 xd = xd + 1
      if (xd.gt.maxlvl) then
         print *,'***ERROR: Insufficient storage in subroutine number'
         print *,'           for stkd()'
         stop
      endif
      stkd(xd) = stkb(xb)
c
c  Be sure to process all nodes in stkc.
c
  120 cx = cx + 1
      if (cx.gt.maxlvl) then
         print *,'***ERROR: Insufficient storage in subroutine number'
         print *,'           for stkc() (index cx)'
         stop
      endif
      if (xc.ge.cx) go to 50
c
c  When stkc is exhausted look for min degree node in same level
c  which has not been processed.
c
      max = ideg + 1
      snd = n + 1
      do 130 i = lst, lnd
        test = lvlst(i)
        if (renum(test).ne.0) go to 130
        if (ndeg(test).ge.max) go to 130
        renum(snd) = 0
        renum(test) = -1
        max = ndeg(test)
        snd = test
  130 continue
      if (snd.eq.n+1) go to 140
      xc = xc + 1
      if (xc.gt.maxlvl) then
         print *,'***ERROR: Insufficient storage in subroutine number'
         print *,'           for stkc()'
         stop
      endif
      stkc(xc) = snd
      go to 50
c
c  If stkd is empty we are done, otherwise copy stkd onto stkc
c  and begin processing new stkc.
c
  140 if (xd.eq.0) go to 160
      if (xd.gt.maxlvl) then
         print *,'***ERROR: Insufficient storage in subroutine number'
         print *,'           for stkc()'
         stop
      endif
      do 150 i = 1, xd
        stkc(i) = stkd(i)
  150 continue
      xc = xd
      go to 40
c
c  Do final bandwidth and profile calculations.
c
  160 do 170 i = 1, n

        if (ipfa(i).gt.ibw2) ibw2 = ipfa(i)
        ipf2 = ipf2 + ipfa(i)
  170 continue
c
c
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine piklvl(lvls1, lvls2, ccstor, idflt, isdir)
c
c piklvl chooses the level structure  used in numbering graph
c lvls1-       on input contains forward leveling info
c lvls2-       on input contains reverse leveling info
c              on output the final level structure chosen
c ccstor-      on input contains connected component info
c idflt-       on input =1 if wdth lvls1.le.wdth lvls2, =2 otherwise
c nhigh        keeps track of level widths for high numbering
c nlow-        keeps track of level widths for low numbering
c nacum-       keeps track of level widths for chosen level structure
c xc-          number of connected components
c size(i)-     size of ith connected component
c stpt(i)-     index into ccstore of 1st node in ith con compt
c isdir-       flag which indicates which way the largest connected
c              component fell.  =+1 if low and -1 if high
c
      integer ccstor, size, stpt, xc, end
c
      common /gra/ n, idpth, ideg
c
c  It is assumed that the graph has at most 'maxcmp' (was 50) components and
c  that there are at most 'maxlvl' (was 100) levels.
c
      integer maxlvl
      parameter (maxlvl=2000,maxcmp=400) ! Note: maxcmp must be less than
                                         !       (maxlvl-1)/2
      common / lvlw / nhigh(maxlvl), nlow(maxlvl), nacum(maxlvl)
      common / cc   / xc, size(maxcmp), stpt(maxcmp)
c     common / ccc  / xc, size(maxcmp), stpt(maxcmp)
c
      dimension lvls1(*), lvls2(*), ccstor(*)
c
c  For each connected component do.
c
      if (xc.gt.maxcmp) then
         print *,'***ERROR: Insufficient storage in piklvl for'
         print *,'          stpt() and size()'
         stop
      endif
      do 80 i = 1, xc
        j   = stpt(i)
        end = size(i) + j - 1
c
c  Set nhigh and nlow equal to nacum.
c
        if (idpth.gt.maxlvl) then
           print *,'***ERROR: Insufficient storage in piklvl for'
           print *,'          nacum(), nhigh() and nlow()'
           stop
        endif
        do 10 k = 1, idpth
          nhigh(k) = nacum(k)
          nlow(k)  = nacum(k)
   10   continue
c
c  Update nhigh and nlow for each node in connected component.
c
        do 20 k = j, end
          inode        = ccstor(k)
          lvlnh        = lvls1(inode)
          if (lvlnh.gt.maxlvl) then
             print *,'***ERROR: pointer in lvls1 out-of-bounds'
             stop
          endif
          nhigh(lvlnh) = nhigh(lvlnh) + 1
          lvlnl        = lvls2(inode)
          if (lvlnl.gt.maxlvl) then
             print *,'***ERROR: pointer in lvls2 out-of-bounds'
             stop
          endif
          nlow(lvlnl)  = nlow(lvlnl) + 1
   20   continue
        max1 = 0
        max2 = 0
c
c  Set max1=largest new number in nhigh.
c  Set max2=largest new number in nlow.
c
        do 30 k = 1, idpth
          if (2*nacum(k).eq.nlow(k)+nhigh(k)) go to 30
          if (nhigh(k).gt.max1) max1 = nhigh(k)
          if (nlow(k).gt.max2) max2 = nlow(k)
   30   continue
c
c  Set it= number of level structure to be used.
c
        it = 1
        if (max1.gt.max2) it = 2
        if (max1.eq.max2) it = idflt
        if (it.eq.2) go to 60
        if (i.eq.1) isdir = -1
c
c  Copy lvls1 into lvls2 for each node in connected component.
c
        do 40 k = j, end
          inode        = ccstor(k)
          lvls2(inode) = lvls1(inode)
   40   continue
c
c  Update nacum to be the same as nhigh.
c
        do 50 k = 1, idpth
          nacum(k) = nhigh(k)
   50   continue
        go to 80
c
c  Update nacum to be the same as nlow.
c
   60   do 70 k = 1, idpth
          nacum(k) = nlow(k)
   70   continue
   80 continue
c
c
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine reduce(ndstk, nr, iold, renum, ndeg, lvl, lvls1,
     * lvls2, ccstor, ibw2, ipf2)
c
c  Subroutine reduce determines a row and column permutation which,
c  when applied to a given sparse matrix, produces a permuted
c  matrix with a smaller bandwidth and profile.
c  The input array is a connection table which represents the
c  indices of the nonzero elements of the matrix, a.  The algo-
c  rithm is described in terms of the adjacency graph which
c  has the characteristic that there is an edge (connection)
c  between nodes i and j if a(i,j) .ne. 0 and i .ne. j.
c
c  Dimensioning information:
c  The following integer arrays must be dimensioned in the calling routine.
c    ndstk(nr,d1)        nr is .ge. maximum degree of all nodes.
c    iold(d2)            d2 and d1 are .ge. the total number of
c    renum(d2+1)         nodes in the graph.
c    ndeg(d2)            storage requirements can be significantly
c    lvl(d2)             decreased for ibm 360 and 370 computers
c    lvls1(d2)           by replacing integer ndstk by
c    lvls2(d2)           integer*2 ndstk in subroutines reduce,
c    ccstor(d2)          dgree, fndiam, tree and number.
c
c  Common information:
c  The following common block must be in the calling routine.
c
c    common / gra / n, idpth, ideg
c
c  Explanation of input variables:
c    ndstk-     connection table representing graph.
c               ndstk(i,j)=node number of ith connection to node
c               number-j.  a connection of a node to itself is not
c               listed.  extra positions must have zero fill.
c    nr-        row dimension assigned ndstk in calling program.
c               this is the maximum number of node point connections
c               allowed in the graph.
c    iold(i)-   numbering of ith node upon input.
c               if no numbering exists then iold(i)=i.
c    n-         number of nodes in graph (equal to order of matrix).
c    ideg-      maximum degree of any node in the graph.
c
c  Explanation of output variables:
c    renum(i)-  the new number for the ith node.
c    ndeg(i)-   the degree of the ith node.
c    ibw2-      the bandwidth after renumbering.
c    ipf2-      the profile after renumbering.
c    idpth-     number of levels in reduce level structure.
c  The following only have meaning if the graph was connected--
c    lvl(i)-    index into lvls1 to the first node in level i.
c               lvl(i+1)-lvl(i)= number of nodes in ith level
c    lvls1(i)-  node numbers listed by level.
c    lvls2(i)-  the level assigned to node i by reduce.
c
c  Working storage variable:
c    ccstor
c
c  Local storage:
c    common/cc/-subroutines reduce, sort2 and piklvl assume that
c               the graph has at most 100 (was 50) connected components.
c               subroutine fndiam assumes that there are at most
c               500 (was 100) nodes in the last level.
c    common/lvlw/-subroutines setup and piklvl assume that there
c               are at most 500 (was 100) levels.
c
c  Use integer*2 ndstk  with an ibm 360 or 370.
c
      integer ndstk
      integer stnode, rvnode, renum, xc, sort2, stnum, ccstor,
     +        size, stpt, sbnum
c
      common /gra/ n, idpth, ideg
c
c  It is assumed that the graph has at most 'maxcmp' (was 50) connected
c  components.
c
      parameter (maxlvl=2000,maxcmp=400) ! Note: maxcmp must be less than
                                         !       (maxlvl-1)/2
      common / cc   / xc, size(maxcmp), stpt(maxcmp)
c     common / ccc  / xc, size(maxcmp), stpt(maxcmp)
      common / lvlw / nhigh(maxlvl), nlow(maxlvl), nacum(maxlvl)
c
      dimension ccstor(*), iold(*)
      dimension ndstk(nr,*), lvl(*), lvls1(*), lvls2(*), renum(*),
     +          ndeg(*)
c
      ibw2 = 0
      ipf2 = 0
c
c  Set renum(i)=0 for all i to indicate node i is unnumbered.
c
      do 10 i = 1, n
        renum(i) = 0
   10 continue
c
c  Compute degree of each node and original bandwidth and profile.
c
      call dgree(ndstk, nr, ndeg, iold, ibw1, ipf1)
c
c  Display the original bandwidth and profile.
c
c     write(6,'(x,''..original bandwidth  = '',i10,
c    +            /,x,''..original profile    = '',i10)')
c    +            ibw1, ipf1
c
c  sbnum = low end of available numbers for renumbering
c  stnum = high end of available numbers for renumbering
c
      sbnum = 1
      stnum = n
c
c  Number the nodes of degree zero.
c
      do 20 i = 1, n
        if (ndeg(i).gt.0) go to 20
        renum(i) = stnum
        stnum    = stnum - 1
   20 continue
c
c  Find an unnumbered node of min degree to start on.
c
   30 lowdg = ideg + 1
      nflg  = 1
      isdir = 1
      do 40 i = 1, n
        if (ndeg(i).lt.lowdg) then
         if (renum(i).le.0) then
          lowdg  = ndeg(i)
          stnode = i
         endif
        endif
   40 continue
c
c  Find pseudo-diameter and associated level structures.
c  stnode and rvnode are the ends of the diam and lvls1 and lvls2
c  are the respective level structures.
c
      call fndiam(stnode, rvnode, ndstk, nr, ndeg, lvl, lvls1,
     * lvls2, ccstor, idflt)
      if (ndeg(stnode).le.ndeg(rvnode)) go to 50
c
c  nflg indicates the end to begin numbering on
c
      nflg   = -1
      stnode = rvnode
   50 call setup(lvl, lvls1, lvls2)
c
c  Find all the connected components  (xc counts them).
c
      xc    = 0
      lroot = 1
      lvln  = 1
      do 60 i = 1, n
        if (lvl(i).ne.0) go to 60
        xc       = xc + 1
        if (xc.gt.maxcmp) then
           print *,'***ERROR: Insufficient storage in reduce for stpt()'
           stop
        endif
        stpt(xc) = lroot
        call tree(i, ndstk, nr, lvl, ccstor, ndeg, lvlwth, lvlbot,
     *   lvln, maxlw, n)
        size(xc) = lvlbot + lvlwth - lroot
        lroot    = lvlbot + lvlwth
        lvln     = lroot
   60 continue
      if (sort2(dmy).eq.0) go to 70
      call piklvl(lvls1, lvls2, ccstor, idflt, isdir)
c
c  On return from piklvl, isdir indicates the direction the largest
c  component fell.  isdir is modified now to indicate the numbering
c  direction.  num is set to the proper value for this direction.
c
   70 isdir = isdir*nflg
      num   = sbnum
      if (isdir.lt.0) num = stnum
      call number(stnode, num, ndstk, lvls2, ndeg, renum, lvls1,
     * lvl, nr, nflg, ibw2, ipf2, ccstor, isdir)
c
c  Update stnum or sbnum after numbering.
c
      if (isdir.lt.0) stnum = num
      if (isdir.gt.0) sbnum = num
      if (sbnum.le.stnum) go to 30
c
c     write(6,'(x,''..minimized bandwidth = '',i10,
c    +            /,x,''..minimized profile   = '',i10)')
c    +            ibw2, ipf2
c
c  If the original bandwidth is the same as the minimum bandwidth
c  keep the old labeling and return.
c
      if (ibw2.le.ibw1) then
c      write(6,'(x,''..original bandwidth is minimized'')')
       return
      endif
c
c  If original numbering is better than new one, set up to return it.
c
      do 80 i = 1, n
        renum(i) = iold(i)
   80 continue
      ibw2 = ibw1
      ipf2 = ipf1
c
c
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine setup(lvl, lvls1, lvls2)
c
c Setup computes the reverse leveling info from lvls2 and stores
c it into lvls2.  nacum(i) is initialized to nodes/ith level for nodes
c on the pseudo-diameter of the graph.  lvl is initialized to non-
c zero for nodes on the pseudo-diam and nodes in a different
c component of the graph.
c
      common / gra / n, idpth, ideg
c
c It is assumed that there are at most 'maxlvl' (was 100) levels.
c
      parameter (maxlvl=2000)
      common / lvlw / nhigh(maxlvl), nlow(maxlvl), nacum(maxlvl)
c
      dimension lvl(*), lvls1(*), lvls2(*)
c
      if (idpth.gt.maxlvl) then
         print *,'***ERROR: Insufficient storage in setup for nacum()'
         stop
      endif
      do 10 i = 1, idpth
        nacum(i) = 0
   10 continue
c
      do 30 i = 1, n
        lvl(i) = 1
        lvls2(i) = idpth + 1 - lvls2(i)
        itemp = lvls2(i)
        if (itemp.gt.idpth) go to 30
        if (itemp.ne.lvls1(i)) go to 20
        nacum(itemp) = nacum(itemp) + 1
        go to 30
   20   lvl(i) = 0
   30 continue
c
c
      return
      end


      integer function sort2(dmy)
c
c Sort2 sorts size and stpt into descending order according to
c values of size. xc=number of entries in each array.
c
      integer temp, ccstor, size, stpt, xc
c
c It is assumed that the graph has at most 'maxcmp' (was 50) connected
c components.
c
      parameter (maxlvl=2000,maxcmp=400) ! Note: maxcmp must be less than
                                         !       (maxlvl-1)/2
      common / cc   / xc, size(maxcmp), stpt(maxcmp)
c     common / ccc  / xc, size(maxcmp), stpt(maxcmp)
c
      sort2 = 0
      if (xc.eq.0) return
      if (xc.gt.maxcmp) then
         print *,'***ERROR: Insufficient storage in sort2 for size()'
         print *,'          and stpt()'
         stop
      endif
c
      sort2 = 1
      ind   = xc
   10 itest = 0
      ind   = ind - 1
      if (ind.lt.1) return
c
      do 20 i = 1, ind
        j = i + 1
        if (size(i).ge.size(j)) go to 20
        itest   = 1
        temp    = size(i)
        size(i) = size(j)
        size(j) = temp
        temp    = stpt(i)
        stpt(i) = stpt(j)
        stpt(j) = temp
   20 continue
c
      if (itest.eq.1) go to 10
c
c
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sortdg(stk1, stk2, x1, x2, ndeg)
c sortdg sorts stk2 by degree of the node and adds it to the end
c of stk1 in order of lowest to highest degree.  x1 and x2 are the
c number of nodes in stk1 and stk2 respectively.
      integer x1, x2, stk1, stk2, temp
      common /gra/ n, idpth, ideg
      dimension ndeg(1), stk1(1), stk2(1)
      ind = x2
   10 itest = 0
      ind = ind - 1
      if (ind.lt.1) go to 30
      do 20 i=1,ind
        j = i + 1
        istk2 = stk2(i)
        jstk2 = stk2(j)
        if (ndeg(istk2).le.ndeg(jstk2)) go to 20
        itest = 1
        temp = stk2(i)
        stk2(i) = stk2(j)
        stk2(j) = temp
   20 continue
      if (itest.eq.1) go to 10
   30 do 40 i=1,x2
        x1 = x1 + 1
        stk1(x1) = stk2(i)
   40 continue
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tree(iroot, ndstk, nr, lvl, iwk, ndeg, lvlwth,
     +                lvlbot, lvln, maxlw, ibort)
c
c  tree drops a tree in ndstk from iroot.
c  lvl-         array indicating available nodes in ndstk with zero
c               entries. tree enters level numbers assigned
c               during execution of this procedure
c  iwk-         on output contains node numbers used in tree
c               arranged by levels (iwk(lvln) contains iroot
c               and iwk(lvlbot+lvlwth-1) contains last node entered)
c  lvlwth-      on output contains width of last level
c  lvlbot-      on output contains index into iwk of first
c               node in last level
c  maxlw-       on output contains the maximum level width
c  lvln-        on input the first available location in iwk
c               usually one but if iwk is used to store previous
c               connected components, lvln is next available location.
c               on output the total number of levels + 1
c  ibort-       input param which triggers early return if
c               maxlw becomes .ge. ibort
c use integer*2 ndstk  with an ibm 360 or 370.
c
      integer ndstk
c
      dimension ndstk(nr,*), lvl(*), iwk(*), ndeg(*)
c
      maxlw  = 0
      itop   = lvln
      inow   = lvln
      lvlbot = lvln
      lvltop = lvln + 1
      lvln   = 1
      lvl(iroot) = 1
      iwk(itop)  = iroot
   10 lvln       = lvln + 1
   20 iwknow     = iwk(inow)
      ndrow      = ndeg(iwknow)
c
      do 30 j = 1, ndrow
        itest = ndstk(j,iwknow)
        if (lvl(itest).ne.0) go to 30
        lvl(itest) = lvln
        itop       = itop + 1
        iwk(itop)  = itest
   30 continue
c
      inow = inow + 1
      if (inow.lt.lvltop) go to 20
      lvlwth = lvltop - lvlbot
      if (maxlw.lt.lvlwth) maxlw = lvlwth
      if (maxlw.ge.ibort) return
      if (itop.lt.lvltop) return
      lvlbot = inow
      lvltop = itop + 1
      go to 10
c
c
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
