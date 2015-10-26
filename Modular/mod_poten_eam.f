!****************************************************************************
!**
!**  MODULE mod_potential : Embedded Atom Method (EAM) Potential Module
!**
!**
!**  Variable Definitions:
!**  ---------------------
!**
!**  Dynamo EAM variables
!**
!**  Contains Required Routines:
!**  ---------------------------
!**
!**  ReadConstitutive - Reads in potential specific data
!**  Nonlocal()       - Computes energy, stress and moduli using nonlocal
!**                     limit of formulation
!**
!****************************************************************************

      Module mod_poten

c--This file is a modified form of "dyn87.inc" which is the include file
c  for dynamo v8.7.  This is to store the interatomic potentials as done
c  in dynamo

      integer, parameter  :: ngrid=1000
      integer, parameter :: ngridar=1000
      integer, parameter :: nelmax=3
      integer ntypes,nrhoar,nrar,nrho,nr
      integer ielement(nelmax),netype(nelmax)
      double precision amass(nelmax)
      double precision frho(ngrid,nelmax),z2r(ngrid,nelmax,nelmax),
     .  rhor(ngrid,nelmax),drho,drad,rcutsq,sqrtrcutsq
      double precision frhoar(ngridar,nelmax),frhoar1(ngridar,nelmax),
     .  frhoar2(ngridar,nelmax),frhoar3(ngridar,nelmax),
     .  frhoar4(ngridar,nelmax),frhoar5(ngridar,nelmax),
     .  frhoar6(ngridar,nelmax),frhoar7(ngridar,nelmax),
     .  rhorar(ngridar,nelmax),
     .  rhorar1(ngridar,nelmax),rhorar2(ngridar,nelmax),
     .  rhorar3(ngridar,nelmax),rhorar4(ngridar,nelmax),
     .  rhorar5(ngridar,nelmax),rhorar6(ngridar,nelmax),
     .  rhorar7(ngridar,nelmax),
     .  z2rar(ngridar,nelmax,nelmax),z2rar1(ngridar,nelmax,nelmax),
     .  z2rar2(ngridar,nelmax,nelmax),z2rar3(ngridar,nelmax,nelmax),
     .  z2rar4(ngridar,nelmax,nelmax),z2rar5(ngridar,nelmax,nelmax),
     .  z2rar6(ngridar,nelmax,nelmax),z2rar7(ngridar,nelmax,nelmax),
     .  drhoar,drar
      integer MapSpecies(103)

      double precision CON1
      parameter (CON1=100.d0)

      Contains

                   !!! REQUIRED ROUTINES !!!

c**------------------------------------------------------------------
c** ReadConstitutive : Read in potential specific data
c**
c--
      subroutine ReadConstitutive()

      use mod_grain
      use mod_material
      implicit none

      !** Local Variables **!
      integer i,imat,ibasis

      ! Make sure Grains have been defined
      if(ngrains.lt.1) then
         write(*,*) '***ERROR: Grain definitions must come before'
         write(*,*) '          constitutive information'
         write(*,*)
         stop
      endif

      ! Initialize species mapping vector
      MapSpecies=0
      ! Read Potentials
      write(*,*)
      write(*,*)
      write(*,*) '-------------here comes some dynamo output ------'
      write(*,*)
      call ReadPoten()
      write(*,*)
      write(*,*) '-------------back to our output -------------'
      write(*,*)
      write(*,*)

      ! Construct species mapping vector
      do i=1,ntypes
         MapSpecies(ielement(i))=i
      enddo

      ! Check that all species in loaded materials are handled by pot
      do imat=1,nmaterials
         do ibasis=1,material(imat)%nbasis
            if(MapSpecies(material(imat)%ispec(ibasis)).eq.0) then
               write(*,*)
               write(*,1000) imat,material(imat)%ispec(ibasis)
               write(*,*)
     &              '   But only the following species are defined'
               write(*,*) '   in the EAM files:'
               do i=1,103
                  if(MapSpecies(i).ne.0) write(*,*) i
               enddo
               write(*,*)
               write(*,*) '*** Run Stopped due to an ERROR'
               stop
            endif
         enddo
      enddo
 1000 format('***ERROR: Material ',i2,' contains species ',i2)

      end subroutine ReadConstitutive

c**-----------------------------------------------------------------------
c** Nonlocal : find energy, force and stiffness contribution due to
c**            a nonlocal representative atom
c**
c**    Parameters :
c**          Look at the element routine
c**
c**
c**    Algorithm :
c**          Make Representative crystallite.
c**
c**
c**
c**    Notes :
c**           Specific to embedded atom.
c**
c**
c**
c**    WARNING : This will work only for NDA = 1
c**
c**
c--
      subroutine NonLocal(id,x,ix,f,b,dr,fls,flw,irep)

      use mod_material
      use mod_global
      use mod_grain
      use mod_dynamo
      implicit none

c--Variables transferred
      integer irep

      double precision  b(ndf,*),x(nxdm,*),f(ndf,*), dr(ndf,*)

      integer  id(ndf,*), ix(nen1,*)

      logical fls,flw,wn

c--Local Variables

      double precision  x1(2),x2(2),x3(2),s0(3)
      double precision, allocatable :: p1(:),p2(:)

      integer ineigh,idf,idx,i0,isp1,isp0,i,j,k,node,l
     $     ,m,ndnn

      double precision rhosum, phisum, r, ef,vol,
     $     df, d2f, ph, dp, d2p,em, dem, d2em, druai,
     $     drubj,sed

      double precision force(3)
      common/debugger/debug
      logical debug

c     hack: assume all atoms are the same species for now, assume only
c     one material defined     
      if(irep.ne.numnpp1) then 
         isp0=atomSpecie(irep)
      else
         isp0=-1
      endif

c     Ensure there is sufficient storage allocated for out-of-balance
c     force vector and stiffness matrix and initialize them.
      ndnn=ndf*(nneips+1)
      if(fls) then
         allocate(p1(ndnn),p2(ndnn))
         p1=0.d0
         p2=0.d0
      endif
c
c     Loop over all atoms in local crystal summing energy
c     and force contributions
      rhosum=0.d0
      phisum=0.d0
      do ineigh = 1, nneips
         if(jneigh(ineigh).ne.numnpp1) then
            isp1=atomSpecie(jneigh(ineigh))
         else
            isp1=-1
         endif
         r=sqrt(rneigh(ineigh))
c         call edens(r,isp1,ef,df,d2f,.true.,fls,.false.)
c         call pair(r,isp0,isp1,ph,dp,d2p,.true.,fls,.false.)
	 Debugflag=0
         call edens(r,isp1,ef,df,d2f,.true.,fls,.false.)
         call pair(r,isp0,isp1,ph,dp,d2p,.true.,fls,.false.)

         rhosum=rhosum+ef
         phisum=phisum+0.5*ph

         if (fls) then
            do idf=1,ndf
               i0=nneips*ndf+idf
               idx=(ineigh-1)*ndf+idf
               druai=dneigh(idf,ineigh)/r
               p1(i0)=p1(i0)+df*druai
               p2(i0)=p2(i0)+0.5d0*dp*druai
               p1(idx)=p1(idx)-df*druai
               p2(idx)=p2(idx)-0.5d0*dp*druai
            enddo
         endif
      enddo


      !Compute embedding energy
      call embed(rhosum,em,dem,d2em,isp0,fls,.false.)

      !Compute  energy    !CHANGED for nalpha (May 10,96)
                          !Changed for atomistic (Aug 25, 96)
      energy(irep) = em+phisum

c     Compute contribution to out-of-balance force
      if (fls) then
         do ineigh=1,nneips+1
            if(ineigh.gt.nneips) then
               node = irep
            else
               node = jneigh(ineigh)
            endif
            do idf=1,ndf
               idx = (ineigh-1)*ndf+idf
               force(idf) = (dem*p1(idx)+p2(idx))
               dr(idf,node)=dr(idf,node)-force(idf)
            end do
         end do

C     Compute Virial Stress FCC material Only

CCJSTest         vol=(material(1)%a0**3)/4.d0
CCJSTest: for Hex Al only
         vol=20.5052

         do ineigh=1,nneips
            node = jneigh(ineigh)
            do idf=1,ndf
               idx = (ineigh-1)*ndf+idf
               force(idf) = (dem*p1(idx)+p2(idx))
               virst(idf,1,node)=virst(idf,1,node)-force(idf)*
     $                        dneigh(1,ineigh)/vol
               virst(idf,2,node)=virst(idf,2,node)-force(idf)*
     $                        dneigh(2,ineigh)/vol
               virst(idf,3,node)=virst(idf,3,node)-force(idf)*
     $                        dneigh(3,ineigh)/vol
           end do
         end do
C
      endif

      if (fls) deallocate(p1,p2)
      return
      end subroutine nonlocal

c**---------------------------------------------------------------------
C     Calculate Electron Density as a function of Distance between Atoms
      subroutine edens(r,ispec,f,df,d2f,lf,ldf,ld2f)

      implicit none

      !** Transferred Variables **!
      double precision, intent(in) :: r
      integer, intent(in) :: ispec
      double precision, intent(out) :: f,df,d2f
      logical, intent(in) :: lf,ldf,ld2f

      !local variables
      double precision CutOffRadius
      f=0.d0
      df=0.d0
      d2f=0.d0
      call rhfl(r,ispec,f,df,d2f,lf,ldf,ld2f)

      return
      end subroutine edens

c**---------------------------------------------------------------------
C     Pair Potential at Distance R
      subroutine pair(r,isp1,isp2,f,df,d2f,lf,ldf,ld2f)

      implicit none

      !** Transferred Variables **!
      double precision, intent(in) :: r
      integer, intent(in) :: isp1,isp2
      double precision, intent(out) :: f,df,d2f
      logical, intent(in) :: lf,ldf,ld2f

      !local variables

      f=0.d0
      df=0.d0
      d2f=0.d0
      call v2fl(r,isp1,isp2,f,df,d2f,lf,ldf,ld2f)

      return
      end subroutine pair

c**---------------------------------------------------------------------
C     Embedding function and derivatives wrt rho
      subroutine femb(rho,ispec,f,df,d2f,lf,ldf,ld2f)

      implicit none

      double precision,intent(in) :: rho
      integer,intent(in) :: ispec
      double precision,intent(out) :: f,df,d2f
      logical,intent(in) :: lf,ldf,ld2f

      call uufl(rho,ispec,f,df,d2f,lf,ldf,ld2f) !VBS 23 Jan 96

      return
      end subroutine femb

c**---------------------------------------------------------------------
C     Calculate Energy to Embed and Atom in Energy Density RHO
      subroutine embed(rho,ff,scon,dcon,ispec,fls,flt)
      implicit none
      double precision,intent(in) :: rho
      double precision,intent(out) :: ff,scon,dcon
      integer,intent(in) :: ispec
      logical,intent(in) :: fls,flt
      call femb(rho,ispec,ff,scon,dcon,.true.,(fls.or.flt),flt)
      return
      end subroutine embed

c**---------------------------------------------------------------------
c** rhfl - gives electron density
c**
c**     Parameters :
c**              r  (in) : distance r
c**            isp  (in) : species of the atom
c**              f (out) : electron density
c**             df (out) : first derivative of electron density
c**            d2f (out) : second derivative of electron density
c**             lf  (in) : true to compute f
c**            ldf  (in) : true to compute df
c**           ld2f  (in) : true to compute ld2f
c**
      subroutine rhfl(r,isp,f,df,d2f,lf,ldf,ld2f)
      implicit double precision (a-h,o-z)
      data one/1.0/
      logical lf,ldf,ld2f

      if(isp.lt.0) then
         if(lf) f=0.d0
         if(ldf) df=0.d0
         if(ld2f) d2f=0.d0
         return
      endif

      p = r/drar + 1.0
      k = p
      k = min0(k,nrar-1)
      p = p - k
      p = min(p,one)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (lf) f =  ((rhorar3(k,isp)*p
     &             + rhorar2(k,isp))*p
     &             + rhorar1(k,isp))*p
     &             + rhorar(k,isp)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (ldf) df = (rhorar6(k,isp)*p
     &            + rhorar5(k,isp))*p
     &            + rhorar4(k,isp)

      if (ld2f) d2f =  rhorar7(k,isp)
     &              + (rhorar7(k+1,isp)
     &              - rhorar7(k,isp))*p

      return
      end subroutine rhfl

c**---------------------------------------------------------------------
c** v2fl - gives pair potential
c**
c**     Parameters :
c**              r  (in) : distance r
c**            ity  (in) : species of the atom1
c**            jty  (in) : species of the atom2
c**              f (out) : pair potential
c**             df (out) : first derivative of pair potential
c**            d2f (out) : second derivative of pair potential
c**             lf  (in) : true to compute f
c**            ldf  (in) : true to compute df
c**           ld2f  (in) : true to compute ld2f
c**
      subroutine v2fl(r,ity,jty,f,df,d2f,lf,ldf,ld2f)
      use mod_global
      implicit none
      double precision one,f,df,d2f,r,p
      integer ity,jty,k

      data one/1.d0/
      logical lf,ldf,ld2f
      common/debugger/debug
      logical debug

      if(ity.lt.0.and.jty.lt.0) then
         if(lf) f=0.d0
         if(ldf) df=0.d0
         if(ld2f) d2f=0.d0
         return
      else if(ity.lt.0.or.jty.lt.0) then
         if(r.gt.INDRAD) stop 'oops'
         if(lf) f=CON1*(r-INDRAD)**2
         if(ldf) df=2*CON1*(r-INDRAD)
         if(ld2f) d2f=2*CON1
         return
      endif

      p = r/drar + 1.0
      k = p
      k = min0(k,nrar-1)
      p = p - k
      p = min(p,one)
      
      if(r .le. 1.0d-6) then
        print *,'**WARNING -- To small r in pair potential',r
CC--Jun Song Test
        Debugflag=1
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (lf.or.ldf.or.ld2f) then
         f =  ((z2rar3(k,ity,jty)*p
     &        + z2rar2(k,ity,jty))*p
     &        +z2rar1(k,ity,jty))*p
     &        + z2rar(k,ity,jty)
         f = f/r
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (ldf.or.ld2f) then
         df = (z2rar6(k,ity,jty)*p
     &      + z2rar5(k,ity,jty))*p +
     &        z2rar4(k,ity,jty)
         df = (df - f)/r
      endif

      if (ld2f) then
         d2f = z2rar7(k,ity,jty)
     &       + ( z2rar7(k+1,ity,jty)
     &       - z2rar7(k,ity,jty) )*p
         d2f = (d2f - 2.*df)/r
      endif

      return
      end subroutine v2fl


c**---------------------------------------------------------------------
c** uufl - gives embedding energy
c**
c**     Parameters :
c**             rh  (in) : electron density
c**            isp  (in) : species of the atom
c**              f (out) : electron density
c**             df (out) : first derivative of electron density
c**            d2f (out) : second derivative of electron density
c**             lf  (in) : true to compute f
c**            ldf  (in) : true to compute df
c**           ld2f  (in) : true to compute ld2f
c**
      subroutine uufl(rh,isp,f,df,d2f,lf,ldf,ld2f)
      implicit double precision (a-h,o-z)
      data one/1.0/
      logical lf,ldf,ld2f

      if(isp.lt.0) then
         if(lf) f=0.d0
         if(ldf) df=0.d0
         if(ld2f) d2f=0.d0
         return
      endif

      p = rh/drhoar + 1.0
      k = p
      k = min0(k,nrar-1)
      p = p - k
      p = min(p,one)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (lf)
     &   f =  ((frhoar3(k,isp)*p
     &     + frhoar2(k,isp))*p
     &     + frhoar1(k,isp))*p
     &        + frhoar(k,isp)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (ldf)
     &   df = (frhoar6(k,isp)*p
     &      + frhoar5(k,isp))*p
     &      + frhoar4(k,isp)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (ld2f)
     &   d2f = frhoar7(k,isp)
     &       + (frhoar7(k+1,isp)
     &       - frhoar7(k,isp))*p

      return
      end subroutine uufl

c**---------------------------------------------------------------------
c** NumAtomSpec() : get the number of atomic species
c**
      integer function NumAtomSpec()

      NumAtomSpec = ntypes

      return
      end function numatomspec

c**--------------------------------------------------------------------
c** ReadPoten - Read the inter atomic potentials from specified file
c**
c**     Parameters :-
c**
c**     Algorithm :-
c**                 Whatever is there in the dynamo code!!
c**
c**
      subroutine ReadPoten
      use mod_file
c      use mod_global
      implicit double precision (a-h,o-z)
      implicit integer (i-n)
      dimension frhoin(ngrid,nelmax),drhoin(nelmax),nrhoin(nelmax)
      dimension rhorin(ngrid,nelmax),zrin(ngrid,nelmax),drin(nelmax),
     1 nrin(nelmax)
      dimension zrtemp(ngrid,nelmax)
      dimension fheader(10,nelmax),rcut(nelmax),blat(nelmax)
      character*8 llat(nelmax),latty
      character*80 sheader(3)
      data conmas/1.0365e-4/
      data two/2.0/
      logical  setflag
      character*80 funcfl(100),setfl
      data funcfl/100*'none'/,setfl/'none'/



c--Read in the function file name
      ipinter = 0
      read(5,*)setflag

      if( .not. setflag) then
         read(5,*)numelmts
         if(numelmts .gt. nelmax) then
             print*,'**ERROR : Too many func files'
             stop
         endif
         do i = 1, numelmts
            read(5,78)funcfl(i)
         end do
 78      format(a)
        i = 0
 10     continue
        i = i + 1
        if (i.gt.nelmax) goto 50
        if (funcfl(i)(1:4).eq.'NONE'.or.
     $      funcfl(i)(1:4).eq.'none') goto 50
        if (.not.FileExists(funcfl(i),.true.)) stop
        call iofile(funcfl(i),'formatted  ',iunit,.true.)
        read(iunit,20)(fheader(j,i),j=1,10)
 20     format(10a8)
        read(iunit,30) ielement(i), amass(i), blat(i), llat(i)
 30     format(i5,2g15.5,a8)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CC--Jun Song: use unformated input for more flexibility
        read(iunit,*) nrhoin(i), drhoin(i), nrin(i), drin(i), rcut(i)
cc9901    format(i5,e24.16,i5,2e24.16)
c
c  assume that z(r) and rho(r) grids coincide
c
        read(iunit,*) (frhoin(j,i),j=1,nrhoin(i))
        read(iunit,*) (zrin(j,i),j=1,nrin(i))
        read(iunit,*) (rhorin(j,i),j=1,nrin(i))
CC9902    format(5e24.16)
c
c close the interaction file
c
        close(iunit)
        go to 10
 50     continue
        ntypes = i - 1
      if (ntypes.gt.nelmax) then
         write(6,*)'error: number of types greater than nelmax'
         stop
       endif
c
c determine common grid spacings and number
c
c
c take largest grid spacing
c take smallest maximum
c
        drad = drin(1)
        drho = drhoin(1)
        rmax = (nrin(1)-1)*drin(1)
        rhomax = (nrhoin(1)-1)*drhoin(1)
        do 80 i1=2,ntypes
        drad = max(drad,drin(i1))
        drho = max(drho,drhoin(i1))
        rmaxi = (nrin(i1)-1)*drin(i1)
        rhomaxi = (nrhoin(i1)-1)*drhoin(i1)
        rmax = max(rmax,rmaxi)
        rhomax = max(rhomax,rhomaxi)
80      continue
        nr = nint(rmax/drad)
        nrho = nint(rhomax/drho)
c
c set up the z(r) and rho(r) grids
c
        do 90 i1=1,ntypes
        do 85 j=1,nr
        r = (j-1)*drad
c
c  do four-point lagrange interpolation
c
        p = r/drin(i1) + 1.0
        k = p
        k = min0(k,nrin(i1)-2)
        k = max0(k,2)
        p = p - k
c       make sure that p is less than 2.0
c       then if r is out of range, p = 2.0 and rhor = last value of rhorin
        p = min(p,two)
        cof1 = -0.166666667*p*(p-1.)*(p-2.)
        cof2 = 0.5*(p**2-1.)*(p-2.)
        cof3 = -0.5*p*(p+1.)*(p-2.)
        cof4 = 0.166666667*p*(p**2-1.)
        rhor(j,i1) = cof1*rhorin(k-1,i1)
     1      + cof2*rhorin(k,i1)
     2      + cof3*rhorin(k+1,i1)
     3      + cof4*rhorin(k+2,i1)
        zrtemp(j,i1) = cof1*zrin(k-1,i1)
     1      + cof2*zrin(k,i1)
     2      + cof3*zrin(k+1,i1)
     3      + cof4*zrin(k+2,i1)
85      continue
90      continue
c
c set up the f(rho) grid
c
        do 100 i1=1,ntypes
        do 95 j=1,nrho
        r = (j-1)*drho
c
c  do four-point lagrange interpolation
c
        p = r/drhoin(i1) + 1.0
        k = p
        k = min0(k,nrhoin(i1)-2)
        k = max0(k,2)
        p = p - k
c       make sure that p is less than 2.0
c       then if r is out of range, p = 2.0 and rhor = last value of rhorin
        p = min(p,two)
        cof1 = -0.166666667*p*(p-1.)*(p-2.)
        cof2 = 0.5*(p**2-1.)*(p-2.)
        cof3 = -0.5*p*(p+1.)*(p-2.)
        cof4 = 0.166666667*p*(p**2-1.)
        frho(j,i1) = cof1*frhoin(k-1,i1)
     1      + cof2*frhoin(k,i1)
     2      + cof3*frhoin(k+1,i1)
     3      + cof4*frhoin(k+2,i1)
95      continue
100     continue
c
c set up the z2 grid (zi*zj)
c
        do 110 i1=1,ntypes
        do 110 i2=1,ntypes
        do 105 j=1,nr
105     z2r(j,i1,i2) = 27.2*0.529*zrtemp(j,i1)*zrtemp(j,i2)
110     continue
        rcutsq = 0.0
        do 120 i=1,ntypes
120     rcutsq = max(rcutsq,rcut(i))
        if(rcutsq.eq.0.0)rcutsq = 5.0
c       print out types
        write(6,9000)
 9000   format(/,/' ******   interactions defined  ')
        write(6,9001)ntypes
 9001   format(1x,i5,' particle types')
        write(6,9002)
 9002   format('   type  element       amass              alat        ',
     1' lattype',/,
     2       '   ----  --------  --------------  ----------------    ',
     3'----------')
        write(6,9003)(i,ielement(i),amass(i),blat(i),llat(i),i=1,ntypes)
 9003   format(1x,i4,i9,g15.5,5x,g15.5,10x,a8)
      do 130 i=1,ntypes
c       print out header
      write(6,9004)
 9004 format('   type  file name  header',/,
     1       '   ----  ---------  ---------------')
      write(6,9005)i,funcfl(i),(fheader(j,i),j=1,10)
 9005 format(1x,i4,5x,a80,3x,10a8)
      write(6,9006)rcut(i)
 9006 format('   cut-off distance =',g15.5)
      if(ipinter.gt.0)then
         write(6,9007)
 9007 format('    r          z        rho              rho        f',/,
     1       ' --------  --------  ----------       --------  --------')
         do 140 j=1,nr
         jr = min0(j,nr)
         jrho = min0(j,nrho)
         r = (jr-1)*drad
         rhotmp = (jrho-1)*drho
         write(6,9008)r,zrtemp(jr,i),rhor(jr,i),rhotmp,frho(jrho,i)
 9008 format(1x,g15.7,1x,g15.7,1x,g15.7,1x,g15.7,1x,g15.7)
 140   continue
      endif
 130  continue


      else
        read(5,78)setfl
        if (.not.FileExists(setfl,.true.)) stop
        call iofile(setfl,'formatted  ',iunit,.true.)
        read(iunit,'(a80)')sheader(1)
        read(iunit,'(a80)')sheader(2)
        read(iunit,'(a80)')sheader(3)
        read(iunit,150) ntypes
 150    format(i5)
      if (ntypes.gt.nelmax) then
         write(6,*)'error: number of types greater than nelmax'
         stop
       endif
        read(iunit,*) nrho, drho, nr, drad, rcutall
        rcutsq = rcutall
        do 160 i=1,ntypes
        read(iunit,30) ielement(i),amass(i),blat(i),llat(i)
        read(iunit,*) (frho(j,i),j=1,nrho)
        read(iunit,*) (rhor(j,i),j=1,nr)
 160    continue
        do 170 i1=1,ntypes
        do 170 i2=1,i1
170     read(iunit,*) (z2r(j,i1,i2),j=1,nr)
        do 176 i1=1,ntypes
        do 175 i2=i1+1,ntypes
        do 174 j=1,nr
        z2r(j,i1,i2) = z2r(j,i2,i1)
174     continue
175     continue
176     continue
        close(iunit)
        write(6,9000)
        write(6,9001)ntypes
        write(6,9002)
        write(6,9003)(i,ielement(i),amass(i),blat(i),llat(i),i=1,ntypes)
        write(6,9111)setfl
9111    format(1x,'    file name  ',a80)
c  print out header
        write(6,9112)
9112    format('   header',/,
     1         '   ________________')
        write(6,'(4x,a80)')(sheader(i),i=1,3)
9115    format(4x,10a8)
        write(6,9006)rcutsq
        if(ipinter.gt.0)then
           do 190 i=1,ntypes
           write(6,9116)i
9116       format(' type ',i5,/,
     1           ' _____________')
           write(6,9117)
9117       format('      r          rho              rho        f',/,
     1 '   --------  ----------       --------  --------')
           do 200 j=1,max0(nr,nrho)
           jr = min0(j,nr)
           jrho = min0(j,nrho)
           r = (jr-1)*drad
           rhotmp = (jrho-1)*drho
           write(6,9118)r,rhor(jr,i),rhotmp,frho(jrho,i)
9118       format(1x,g15.7,1x,g15.7,1x,g15.7,1x,g15.7)
200        continue
190        continue
           do 210 i1=1,ntypes
           do 210 i2=1,i1
           write(6,9126)i1,i2
9126        format(' types ',i5,i5,/,
     1            ' ______________________')
           write(6,9127)
9127       format('      r            z**2',/,
     1           '   _______      _____________')
           do 220 j=1,nr
           r = (j-1)*drad
           write(6,9128)r,z2r(j,i1,i2)
9128       format(1x,g15.7,1x,g15.7)
220        continue
210        continue
        endif
      endif




c
c       this is the only place where the mass is in amu
c       here we convert to eV-psec**2/angstrom**2
c       this is the unit used throughout the program
c       restart assumes this mass unit
c
      do 75 i=1,ntypes
75    amass(i) = conmas*amass(i)
      write(6,9009)rcutsq
 9009 format('  use this cut-off distance: ',g15.5)
      sqrtrcutsq=rcutsq
      rcutsq = rcutsq**2
c
c       set the lattice constant to that for type 1 by default
      alat = blat(1)
      latty = llat(1)
c
c
c  now set up the dense grids
c
      nrhoar = nrho
      drhoar = drho
      nrar = nr
      drar = drad
c
c  compute the value and slope at the computational grid points
c  the methodology is as follows.
c  At each point, the slope is estimated from a 5-point Lagrange interpolation
c  of the data.
c  Between each pair of points, a cubic polynomial is used which is fit to the
c  value and slope at the two points.
c  This yields an interpolation which gives continuous value and first derivatives
c  without introducing the long-range effects of glitches in the data that
c  results from using splines.  (i.e. the procedure is local)
c
      do 500 i = 1,ntypes
         do 510 j = 1,nrhoar
            frhoar(j,i) = frho(j,i)
  510       continue
         do 515 j = 1,nrar
            rhorar(j,i) = rhor(j,i)
  515       continue
         frhoar1(1,i) = frhoar(2,i)-frhoar(1,i)
         frhoar1(2,i) = 0.5*(frhoar(3,i)-frhoar(1,i))
         frhoar1(nrhoar-1,i)=0.5*(frhoar(nrhoar,i)-frhoar(nrhoar-2,i))
         frhoar1(nrhoar,i) = frhoar(nrhoar,i)-frhoar(nrhoar-1,i)
         frhoar7(1,i) = ( frhoar(3,i)+frhoar(1,i)
     $                        -2.*frhoar(2,i))/(drhoar*drhoar)
         frhoar7(2,i) = ( frhoar(3,i)+frhoar(1,i)
     $                        -2.*frhoar(2,i))/(drhoar*drhoar)
         frhoar7(nrhoar-1,i) = ( frhoar(nrhoar-2,i)+frhoar(nrhoar,i)
     $                        -2.*frhoar(nrhoar-1,i))/(drhoar*drhoar)
         frhoar7(nrhoar,i) = ( frhoar(nrhoar-2,i)+frhoar(nrhoar,i)
     $                        -2.*frhoar(nrhoar-1,i))/(drhoar*drhoar)

         do 520 j = 3,nrhoar-2
            frhoar1(j,i) = ((frhoar(j-2,i)-frhoar(j+2,i)) +
     $                  8.*(frhoar(j+1,i)-frhoar(j-1,i)))/12.
            frhoar7(j,i) = ( 16.*(frhoar(j-1,i)+frhoar(j+1,i)) -
     $                      (frhoar(j+2,i)+frhoar(j-2,i)) -
     $                       30.*frhoar(j,i) )/(12.*drhoar*drhoar)
  520       continue
         do 525 j = 1,nrhoar-1
            frhoar2(j,i) = 3.*(frhoar(j+1,i)-frhoar(j,i))
     $                   - 2.*frhoar1(j,i) - frhoar1(j+1,i)
            frhoar3(j,i) = frhoar1(j,i) + frhoar1(j+1,i)
     $                   - 2.*(frhoar(j+1,i)-frhoar(j,i))
  525       continue
         frhoar2(nrhoar,i) = 0.
         frhoar3(nrhoar,i) = 0.
         do 528 j = 1,nrhoar
            frhoar4(j,i) = frhoar1(j,i)/drhoar
            frhoar5(j,i) = 2.*frhoar2(j,i)/drhoar
            frhoar6(j,i) = 3.*frhoar3(j,i)/drhoar
  528       continue
         rhorar1(1,i) = rhorar(2,i)-rhorar(1,i)
         rhorar1(2,i) = 0.5*(rhorar(3,i)-rhorar(1,i))
         rhorar1(nrar-1,i) = 0.5*(rhorar(nrar,i)-rhorar(nrar-2,i))
         rhorar1(nrar,i) = 0.
         rhorar7(1,i) = ( rhorar(3,i)+rhorar(1,i)
     $                        -2.*rhorar(2,i))/(drar*drar)
         rhorar7(2,i) = ( rhorar(3,i)+rhorar(1,i)
     $                        -2.*rhorar(2,i))/(drar*drar)
         rhorar7(nrar-1,i) = ( rhorar(nrar-2,i)+rhorar(nrar,i)
     $                        -2.*rhorar(nrar-1,i))/(drar*drar)
         rhorar7(nrar,i) = ( rhorar(nrar-2,i)+rhorar(nrar,i)
     $                        -2.*rhorar(nrar-1,i))/(drar*drar)
         do 530 j = 3,nrar-2
            rhorar1(j,i) = ((rhorar(j-2,i)-rhorar(j+2,i)) +
     $                  8.*(rhorar(j+1,i)-rhorar(j-1,i)))/12.
            rhorar7(j,i) = ( 16.*(rhorar(j-1,i)+rhorar(j+1,i)) -
     $                      (rhorar(j+2,i)+rhorar(j-2,i)) -
     $                       30.*rhorar(j,i) )/(12.*drar*drar)
  530       continue

         do 535 j = 1,nrar-1
            rhorar2(j,i) = 3.*(rhorar(j+1,i)-rhorar(j,i))
     $                   - 2.*rhorar1(j,i) - rhorar1(j+1,i)
            rhorar3(j,i) = rhorar1(j,i) + rhorar1(j+1,i)
     $                   - 2.*(rhorar(j+1,i)-rhorar(j,i))
  535       continue
         rhorar2(nrar,i) = 0.
         rhorar3(nrar,i) = 0.
         do 538 j = 1,nrar
            rhorar4(j,i) = rhorar1(j,i)/drar
            rhorar5(j,i) = 2.*rhorar2(j,i)/drar
            rhorar6(j,i) = 3.*rhorar3(j,i)/drar
  538       continue
      i1 = i
      do 540 i2 = 1,ntypes
         do 550 j = 1,nrar
            z2rar(j,i1,i2) = z2r(j,i1,i2)
  550       continue
         z2rar1(1,i1,i2) = z2rar(2,i1,i2)-z2rar(1,i1,i2)
         z2rar1(2,i1,i2) = 0.5*(z2rar(3,i1,i2)-z2rar(1,i1,i2))
         z2rar1(nrar-1,i1,i2) =
     $        0.5*(z2rar(nrar,i1,i2)-z2rar(nrar-2,i1,i2))
         z2rar1(nrar,i1,i2) = 0.

         z2rar7(1,i1,i2) = ( z2rar(3,i1,i2)+z2rar(1,i1,i2)
     $                        -2.*z2rar(2,i1,i2))/(drar*drar)
         z2rar7(2,i1,i2) = ( z2rar(3,i1,i2)+z2rar(1,i1,i2)
     $                        -2.*z2rar(2,i1,i2))/(drar*drar)
         z2rar7(nrar-1,i1,i2) = (z2rar(nrar-2,i1,i2)+z2rar(nrar,i1,i2)
     $                        -2.*z2rar(nrar-1,i1,i2) )/(drar*drar)
         z2rar7(nrar,i1,i2) = (z2rar(nrar-2,i1,i2)+z2rar(nrar,i1,i2)
     $                        -2.*z2rar(nrar-1,i1,i2) )/(drar*drar)

         do 560 j = 3,nrar-2
            z2rar1(j,i1,i2) = ((z2rar(j-2,i1,i2)-z2rar(j+2,i1,i2)) +
     $          8.*(z2rar(j+1,i1,i2)-z2rar(j-1,i1,i2)))/12.
          z2rar7(j,i1,i2) = ( 16.*(z2rar(j-1,i1,i2)+z2rar(j+1,i1,i2))-
     $                      (z2rar(j+2,i1,i2)+z2rar(j-2,i1,i2)) -
     $                       30.*z2rar(j,i1,i2) )/(12.*drar*drar)
  560       continue

         do 565 j = 1,nrar-1
            z2rar2(j,i1,i2) = 3.*(z2rar(j+1,i1,i2)-z2rar(j,i1,i2))
     $                   - 2.*z2rar1(j,i1,i2) - z2rar1(j+1,i1,i2)
            z2rar3(j,i1,i2) = z2rar1(j,i1,i2) + z2rar1(j+1,i1,i2)
     $                   - 2.*(z2rar(j+1,i1,i2)-z2rar(j,i1,i2))
  565       continue
         z2rar2(nrar,i1,i2) = 0.
         z2rar3(nrar,i1,i2) = 0.
         do 568 j = 1,nrar
            z2rar4(j,i1,i2) = z2rar1(j,i1,i2)/drar
            z2rar5(j,i1,i2) = 2.*z2rar2(j,i1,i2)/drar
            z2rar6(j,i1,i2) = 3.*z2rar3(j,i1,i2)/drar
  568       continue
  540    continue
  500    continue

c
c   output dense grids
c
      if(ipinter.gt.1)then
           do 790 i=1,ntypes
           write(6,9116)i
           write(6,9317)
9317       format('      r          rho         rhop',
     1            '        rho        f         fp',/,
     2            '   --------  ----------    --------',
     3            '      --------  --------  ----------')
           do 800 j=1,max0(nrar,nrhoar)
           jr = min0(j,nrar)
           jrho = min0(j,nrhoar)
           r = (jr-1)*drar
           rhotmp = (jrho-1)*drhoar
           write(6,9318)r,rhorar(jr,i),rhorar1(jr,i),
     1        rhotmp,frhoar(jrho,i),frhoar1(jrho,i)
9318    format(1x,g15.7,1x,g15.7,1x,g15.7,1x,g15.7,1x,g15.7,
     1  1x,g15.7)
800        continue
790        continue
           do 810 i1=1,ntypes
           do 810 i2=1,i1
           write(6,9126)i1,i2
           write(6,9327)
9327       format('      r            z**2          z**2p',/,
     1            '   _______      _____________  ____________')
           do 820 j=1,nrar
           r = (j-1)*drar
           write(6,9328)r,z2rar(j,i1,i2),z2rar1(j,i1,i2)
9328       format(1x,g15.7,1x,g15.7,1x,g15.7)
820        continue
810        continue
      endif
      return
      end subroutine readpoten
      end module mod_poten

c**---------------------------------------------------------------------
c** CutoffRadius() : Get the cut off radius
c**
      double precision function CutoffRadius(i)
      use mod_poten
      use mod_global
      implicit none
      integer i

      if(i.gt.0) then
         CutoffRadius = sqrtrcutsq
      else
         CutoffRadius=INDRAD
      endif

      return
      end function CutoffRadius

c**---------------------------------------------------------------------
c** CutoffR2() : Get the cut off radius
c**
      double precision function CutoffR2(i)
      use mod_poten
      use mod_global
      implicit none
      integer i

      if(i.gt.0) then
         CutoffR2 = rcutsq
      else
         CutoffR2 = INDRADSQ
      endif
      return
      end function CutoffR2

