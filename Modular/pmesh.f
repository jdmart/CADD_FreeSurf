c**---------------------------------------------------------------
c**    pmesh  :  control routine for the macro "feap".
c**
c**   Non-Obvious Parameters : NONE
c**
c--
      subroutine pmesh(id,x,ix,f,b,dr,itx)
      use mod_global
      use mod_poten
      use mod_grain
      use mod_material
      use mod_boundary
      implicit none
c
c---- data input routine for mesh description
c
      integer id(ndf,*),ix(nen1,*),itx(3,*)
      double precision b(ndf,*),f(ndf,*),x(nxdm,*),dr(ndf,*)
c
      logical pcomp,init,die
      character*4 va(2),xm(3),xv(3),wd(25),cd(3),te(3),fd(3),u0(3),cc(2)
      character*4 bl
      character*80 geomfile,input

      integer list,n,i,ii,lower,upper,next,ldif


      data wd/'xxxx','xxxx','mate','xxxx','xxxx','xxxx','end ','xxxx',
     1        'xxxx','xxxx','xxxx','xxxx','xxxx','xxxx','xxxx','xxxx',
     2        'mp01','xxxx','xxxx','xxxx','cons','grai','xxxx','xxxx',
     3        'xxxx'/
      data bl/'blan'/,list/25/
      data va/' val','ue  '/
      data xm/' mas','ses ','    '/,xv/' vel','ocit','ies '/,
     1     cd/' coo','rdin','ates'/,te/' tem','pera','ture'/,
     2     fd/' for','ce/d','ispl'/,u0/' dis','pl. ','    '/
c
c
c---- initialize arrays
      init = .true.
      do 101 n = 1,numnp
      do 101 i = 1,nxdm
      if (x(i,n).ne.0.) init = .false.
101   continue
      if (init) then
      do 102 n = 1,numnp
102   x(1,n) = 0.0
      end if
c---- read macro cards
10    continue
      read(5,'(a80)') input
      upper = 0
      do ii = 1,2
         lower = upper
         upper = next(lower,input)
         ldif = upper - lower - 1
         if (ldif.eq.0) then
            cc(ii) = '    '
         else
            call keyequal(cc(ii),input,lower+1,upper-1)
         end if
      enddo
      do 20 i = 1,list
20    if(pcomp(cc(1),wd(i))) go to 30
      go to 10
c---- process macros
 30   go to (1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19
     @     ,21,22,23,24,25,26,27),i
c
c---- macro 'xxxx'
1     continue
      go to 10
c
c---- macro 'xxxx'
2     continue
      go to 10
c
c---- macro 'mate'
c---- material data input
3     continue
c--Read Material Property table (Ellad's routine)
      write(*,*) '** Reading Material Data'
      write(*,*)
      call ReadMaterials(ndf,cc(2))
      call OutputMaterials()
c
c     Verify user input is correct for this element
c
      if (ndm.ne.2) then
         print *,'***ERROR: Incorrect spatial dimension (ndm=2)'
         stop
      endif
      if (nen.ne.3) then
         print *,'***ERROR: Incorrect #nodes per element (nen=3)'
         stop
      endif
      if (nsdm.ne.6) then
         print *,'***ERROR: Incorrect stress dimension (nsdm=6)'
         stop
      endif
      if (nquad.ne.1) then
         print *,'***ERROR: Incorrect #quadrature points (nquad=1)'
         stop
      endif
      if (nad.ne.0) then
         print *,'***ERROR: Incorrect #internal nodes (nad=0)'
         stop
      endif
      go to 10
c
c---- macro 'xxxx'
    4 continue
      go to 10
c
c---- macro 'xxxx'
c---- define the model boundaries
5     continue
      go to 10
c
c---- macro 'xxxx'
6     continue
      go to 10
c
c---- macro 'end '
c---- terminate mesh input
7     continue
c
c check that some basic stuff is properly defined
c
      die=.false.
      if (nmaterials.eq.0) then
         write(*,*) '**ERROR: no materials are defined'
         die=.true.
      endif
      if (ngrains.eq.0) then
         write(*,*) '**ERROR: no grains are defined'
         die=.true.
      endif
      do i=1,103
         if(MapSpecies(i).ne.0) go to 99
      enddo
         write(*,*) '**ERROR: no constitutive information found'
         die=.true.
 99   continue
      if(die) stop
      write(6,'(/2x,''**end of mesh definition**''/)')
      return
c
c---- macro 'xxxx'
 8    continue
      go to 10
c
c---- macro 'xxxx'
 9    continue
      go to 10
c
c---- macro 'xxxx'
11    continue
      go to 10
c
c---- macro 'xxxx'
12    go to 10
c
c---- macro 'xxxx'
13    continue
      go to 10
c
c---- macro 'xxxx'
   14 continue
      go to 10
c
c---- macro 'xxxx'
   15 continue
      go to 10
c
c---- macro 'xxxx'
   16 continue
      go to 10
c
c---- macro 'xxxx'
   17 continue
      go to 10
c
c---- macro 'mp01'
c---- user supplied macro
   18 continue
      call mp01(id,x,ix,f,b,itx)
      go to 10
c
c---- macro 'xxxx'
   19 continue
      go to 10
c
c---- macro 'xxxx'
   21 continue
      go to 10
c
c---- macro 'xxxx'
   22 continue
      go to 10
c
c---- macro 'cons'
c---- read in potential specific constitutive input
 23   continue
      call ReadConstitutive()
      go to 10
c
c---- macro 'grai'
c---- user supplied macro
 24   continue
c--Read in filename
      if(cc(2).eq.'dire') then
         geomfile=' '
      else
         read(5,9000)geomfile
 9000    format(a)
      endif

c--Read data form  file
      call ReadGrainData(geomfile,cc(2))
c
      call ProcessGrains
c--Print out grain structure
      call OutputGrainData()
      if(ngrains.gt.1) stop 'hardwired for 1 grain'
      call rotateburgers(grains(1)%xlatvect(1:3,1:3),grains(1
     $     )%rotmat(1:2),material(grains(1)%matgrain)%a0
     $     ,material(grains(1)%matgrain)%structure)
      go to 10
c
c---- macro 'xxxx'
 25   continue
      go to 10
c
c---- macro 'xxxx'
 26   continue
      go to 10
c
c---- macro 'xxxx'
 27   continue
      go to 10
c
2003  format(/5x,'material set',i3,' for element type',i2)
2004  format(/' ',a80//5x,'material properties')
2006  format(i10,9e13.3)
      end


