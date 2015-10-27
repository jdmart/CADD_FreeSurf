!> --- This module contains the entire DD methodology
!! --- Generate slip planes
!! --- Generate sources to a given density 
!! --- Generate obstacles to a given spacing lobs
!! --- Calculate the P-K force (direct or multipole) -- currently
!!     only direct
!! --- Calculate the displacement at any location   
!! --- Calculate the stress at any given point
!! --- Handling of pinned dislocations 

      module mod_dd_slip
      integer  mxnslp, mxnsys, nsys, nslp 
      parameter(mxnslp=50000)
      parameter(mxnsys=3)
!<    mxnslp --> max. number of slip planes
!<    mxnsys --> max. number of slip systems (default = 3)
!<    nsys   --> Current number of slip systems (default = 3)
!<    nslp   --> Current number of slip planes
      double precision :: dtime_dd, tnuc, tincr, tincr_sav
      parameter(tnuc = 1.d-9)
      parameter(tincr_sav = 1.d-11)
!<    tincr ---> DD timer (default = 1.d-11)
!<    tnuc --->  Nucleation timer (default = 1.d-9). 
!<    dtime_dd ---> running DD time step
!!     Currently tnuc and tincr needs to be changed in the code

      integer ::  locphi(mxnslp)
      integer :: nntot,tot_source,n_active_slip !! 1/100 of all slip planes
      double precision :: process_xmin, process_xmax,
     $     process_ymin, process_ymax, process_area
      double precision :: source_den, avg_source_str,
     $     sd_source_str
      double precision :: atom_xmin, atom_xmax,
     $     atom_ymin, atom_ymax
      double precision :: e_dd, mu_dd, nu_dd
      double precision :: phisys, phislp(mxnsys), sinphi(mxnsys), 
     &                    cosphi(mxnsys)
      double precision,allocatable :: xslp(:), yslp(:),xendslp(:),
     &                                yendslp(:)
      integer,allocatable :: nshift_slip(:)
      integer, allocatable :: elem_slip(:)
!<     For elements in the DB this contains the slip plane number
      integer :: slipspc, inslip(3)  !< in burgers vectors
      parameter (slipspc = 100)
      double precision :: dx_slip, dy_slip
      double precision :: sdis_out_tol, disl_range_tol
!      common/slpsys/phisys,phislp,nsys
!      common/sliphi/locphi,sinphi,cosphi
!      common/slipxy/xslp,yslp,nslp
!      common/slipend/xendslp, yendslp

!    Nculeation common block       
      integer ::  mxnnuc
      integer, allocatable :: nnuc(:)
      parameter(mxnnuc = 100)
      double precision::  Bdrag, xLe, d, blen, tobs, dstar, vcutoff,
     & 	      	          tref, htinc, epsx, epss

      parameter(vcutoff = 1e13) ! A/s
      parameter(Bdrag = 5.d-15) ! ev/A^3 s
!      parameter(Bdrag = 0.0d0) ! ev/A^3 s
      double precision,allocatable :: t_FR(:,:), xLnuc(:,:),
     &                    tnlaps(:,:), taui(:,:),snuc(:,:), rnuc(:,:,:)
!      double precision, allocatable :: nuc_pk_stress(:,:), nuc_pk_force(:,:)
      integer, allocatable :: elem_source(:,:)
  
!      Obstacles common block
      integer :: tot_obs, mxnobs
      parameter(mxnobs=100)
      integer, allocatable :: nobs(:)
      double precision,allocatable :: tau_obs(:,:), sobs(:,:)
      double precision :: lobs, lobs_max, lobs_min 
!     (obstacle spacing Angstroms)
    
!      common/dprop/Bdrag,xLe,d,blen,tau_obs,dstar,vcutoff
!      common/nuclea/nnuc,snuc
!      common/nucles/t_FR,xLnuc
!      common/nuclock/tnlaps
!      common/timpar/tref,tincr,htinc,epsx,epss,tnuc
!      common/frankr/taui

!c     Dislocations Common block 
      
      integer  :: mxndis, npile, nd, nloop, tot_disl
      integer, allocatable ::  ndis(:)
      parameter (mxndis = 100)
      double precision, allocatable ::  sdis(:,:), b_dd(:,:),
     &	     	          vprev(:,:), bout(:,:),sdis_out(:,:),
     $     sdis_out1(:,:)
      double precision, allocatable :: rdis(:,:,:), rold(:,:,:)
      double precision :: rpile, dlim
      integer,allocatable :: iobpin(:,:), jcnpin(:,:), idple(:,:),
     $     elem_dis(:,:)

 !     common/disl/ndis,sdis,idple
 !     common/disv/vprev
 !     common/burger/b_dd


  !    common/pindis/iobpin,jcnpin
  !    common/dpile/rpile,npile,nd,nloop
  !    common/nuclim/dlim
      integer :: tot_size       !> Total size of dislocations + nucleation sites

      double precision :: range_ !> Range for multipole calculations 
      type DD
!<     Discrete dislocation object 
!!     Contains all the variables relating to dislocations 
      double complex :: xy      !< (x,y) coordinate of dislocation 
      double precision :: slipangle !< Angle of slip plane 
      double precision, dimension(3) :: burgers_dd !< Burgers vector of dislocation 
      double precision, dimension(3) :: stress !< Stress on dislocation
!!     stress(1) --> s11, stress(2) --> s22, stress(3) --> s12
      double complex :: force, forced, disp, dispd, velocity,
     $     forcedg, forceg
!<     force --> PK force on dislocation 
!!     forced --> PK force on dislocation due to direct calculation
!!     disp ---> displacement due to dislocation 
!!     forceg ---> gradient of PK force along the slip plane
		double complex  :: alpha, alpha1, alpha2, alpha3, alpha4
     $               , alpha5, alpha6, alpha7
		! --- Gradient Direct terms
		double complex  :: alpha1g, alpha2g, bialphag,
     $               beta1g, beta2g, bibetag
		double complex 	:: alphad1, betad1, alphad2,
     $               betad2
		double complex  :: beta, beta1, beta2, beta3,
     $               beta4, beta5, beta6, beta7
		double complex 	:: bialpha, bibeta, phi1,
     $               phi2, phi11, phid1
!>    Multipole terms 
		! --- Gradient terms
		double complex  :: phi1g, phi2g, phi11g
		double complex 	:: dphi1, dphi2, dphi11, dphid1
		integer		:: stype, slipplane_no
        end type DD




      contains
      subroutine gen_slip_planes(xmin, xmax,ymin,ymax,
     &     atomxmin, atomxmax, atomymin,atomymax, slip_angle,
     &     burg_len, xslip_start,yslip_start,dxslip, dyslip, pad_width)
      implicit none

      double precision :: xmin, xmax, ymin, ymax, burg_len
      double precision :: atomxmin, atomymin, atomxmax, atomymax
      double precision :: xslip_start, yslip_start, dxslip, dyslip
      double precision :: slip_angle(3)
!     Local Variables 
      integer :: i,j,k,l,islp,ii,nshift, jslp, li,lj
      double precision :: xstart, ystart, xend, yend, tmp
      double precision :: xi,xj
      double precision :: pad_width


      
c      disl_range_tol = pad_width
c      sdis_out_tol = pad_width/2.0d0

      disl_range_tol = 6.0d0
      sdis_out_tol = disl_range_tol

      dtime_dd = 0.0d0
      nsys = 3
      dx_slip = dxslip
      dy_slip = dyslip
      print *, '---------------------------------------'
      print *, 'Generating slip planes'
      print *, xmin, xmax, atomxmin, atomxmax
      print *, ymin, ymax, atomymin, atomymax
      print *, slip_angle(1), slip_angle(2), slip_angle(3)
      print *, burg_len
      print *, dxslip, dyslip
      print *, xslip_start, yslip_start
      print *, '---------------------------------------'
!     Generate bottom and top halves separately ? 
!     Currently it is set up so that both are generated separately
!     to account for the crack face 

      process_xmin = 0.75*xmin
      process_xmax = 0.75*xmax
      process_ymin = 0.75*ymin
      process_ymax = 0.75*ymax
      
      xmin = 0.75*xmin
      xmax = 0.75*xmax
      ymin = 0.75*ymin
      ymax = 0.75*ymax

      process_area = (process_xmax - process_xmin)
     $     *(process_ymax-process_ymin) ! Angstrom^2

      atom_xmin = atomxmin
      atom_xmax = atomxmax
      atom_ymin = atomymin
      atom_ymax = atomymax

      blen = burg_len
      xLe = 6.d0*blen
      nslp = 0
      do i = 1, 3
         phislp(i) = slip_angle(i)
         cosphi(i) = cos(slip_angle(i))
         sinphi(i) = sin(slip_angle(i))
         print *, 'Slip angle', phislp(i)
         if (i .ne. 3) then 
            tmp = ymax/abs(tan(slip_angle(i)))
            inslip(i) = 2*int((xmax-xmin)/dxslip)
         else
            inslip(i) = 2*int((ymax)/(dyslip))
         endif
         nslp = nslp + inslip(i)
         print *, 'No. of slip systems =', inslip(i)
      end do
      print *, 'Total no. of slip systems', nslp

      allocate(xslp(nslp),yslp(nslp),xendslp(nslp),yendslp(nslp))
      allocate(nshift_slip(nslp))
      allocate(elem_slip(nslp))
      xslp = 0.0d0
      yslp = 0.0d0
      xendslp = 0.0d0
      yendslp = 0.0d0
      

!     Now actually generate slip_planes through the whole sample
!     Count the actual no. of slip planes and compare with nslp
      do i = 1, 3
         if (i .ne. 3) then 
            if (i == 1) then 
               xstart = xmin + dxslip
!     ystart is hard coded to be 0.0 for crack
               ystart = 0.01d-3
               islp = 0
               do while (xstart < xmax) 
                  xend = xstart + (ymax-ystart)/tan(slip_angle(i)) 
                  if (xend > xmax) then 
                     xend = xmax 
                     yend = ystart + (xend-xstart)*tan(slip_angle(i))
                  else
                     yend = ymax
                  endif
                  islp = islp + 1
                  locphi(islp) = 1
                  xslp(islp) = xstart
                  xendslp(islp) = xend 
                  yslp(islp) = ystart 
                  yendslp(islp) = yend
                  xstart = xstart + dxslip
!                  xend = xstart + (ymax-ystart)/tan(slip_angle(i))
c$$$                  if (mod(islp,25) .eq. 0) then 
c$$$                     write(*,fmt='(I7,1X,4(1X,E15.7))')
c$$$     $                    islp, xslp(islp), xendslp(islp), yslp(islp),
c$$$     $                    yendslp(islp)
c$$$                  endif
               end do
               print *, 'End of first slip =', islp
               inslip(i) = islp
            endif
            if (i == 2) then 
               xstart = xmin + dxslip
!     ystart is hard coded to be 0.0 for crack
               ystart = 0.01d-3
!               xend =xstart+(ymax-ystart)/tan(slip_angle(i))
!               print *, xstart, xend
               !islp = 0
               j = 0
               do while (xstart < xmax) 
                  j = j + 1
                  xend = xstart + (ymax-ystart)/tan(slip_angle(i))
                  if (xend < xmin) then 
                     xend = xmin
                     yend = ystart + (xend-xstart)*tan(slip_angle(i))
                  else
                     yend = ymax
                  endif
                  islp = islp + 1
                  locphi(islp) = 2
                  xslp(islp) = xstart
                  xendslp(islp) = xend
                  yslp(islp) = ystart
                  yendslp(islp) = yend
                  xstart = xstart + dxslip
!                  xend = xstart + (ymax-ystart)/tan(slip_angle(i))
c$$$                  if (mod(islp,25) .eq. 0) then 
c$$$                     write(*,fmt='(I7,1X,4(1X,E15.7))')
c$$$     $                    islp, xslp(islp), xendslp(islp), yslp(islp),
c$$$     $                    yendslp(islp)
c$$$                  endif

                  !print *, islp, xslp(islp), xendslp(islp)
               end do
               print *, 'End of second slip =', islp
               inslip(i) = islp
            endif
         else
            xstart = xmin
            xend = xmax
            ystart = dyslip
c$$$            nqslp = islp
            j = 0
            do while(ystart < ymax) 
               islp = islp + 1
               j = j+1
               locphi(islp) = 3
               xslp(islp) = xstart
               xendslp(islp) = xend
               yslp(islp) = ystart
               yendslp(islp) = ystart
               ystart = ystart + dyslip
c$$$                  if (mod(islp,25) .eq. 0) then 
c$$$                     write(*,fmt='(I7,1X,4(1X,E15.7))')
c$$$     $                    islp, xslp(islp), xendslp(islp), yslp(islp),
c$$$     $                    yendslp(islp)
c$$$                  endif

            end do

             print *, 'End of third slip =', islp,j
             inslip(i) = islp
         endif
      end do
      print *, inslip(1), inslip(2), inslip(3)
      nshift = inslip(1)+inslip(2)+inslip(3)
!     Now mirror slip planes
!----------------------------------------------------------------------

      do i = 1, 3
         if (i .ne. 3) then 
            if (i == 1) then 
               xstart = xmin + dxslip
!     ystart is hard coded to be 0.0 for crack
               ystart = -0.01d-3
               do while (xstart < xmax) 
                  xend = xstart + (ymin-ystart)/tan(slip_angle(i)) 
                  if (xend < xmin) then 
                     xend = xmin 
                     yend = ystart + (xend-xstart)*tan(slip_angle(i))
                  else
                     yend = ymin
                  endif
                  islp = islp + 1
                  locphi(islp) = 1
                  xslp(islp) = xstart
                  xendslp(islp) = xend 
                  yslp(islp) = ystart 
                  yendslp(islp) = yend
                  xstart = xstart + dxslip
!                  xend = xstart + (ymax-ystart)/tan(slip_angle(i))
c$$$                  if (mod(islp,25) .eq. 0) then 
c$$$                     write(*,fmt='(I7,1X,4(1X,E15.7))')
c$$$     $                    islp, xslp(islp), xendslp(islp), yslp(islp),
c$$$     $                    yendslp(islp)
c$$$                  endif
               end do
               print *, 'End of first slip =', islp
               inslip(i) = islp
            endif
            if (i == 2) then 
               xstart = xmin + dxslip
!     ystart is hard coded to be 0.0 for crack
               ystart = -0.01d-3
!               xend =xstart+(ymax-ystart)/tan(slip_angle(i))
!               print *, xstart, xend
               !islp = 0
               j = 0
               do while (xstart < xmax) 
                  j = j + 1
                  xend = xstart + (ymin-ystart)/tan(slip_angle(i))
                  if (xend > xmax) then 
                     xend = xmax
                     yend = ystart - (xstart-xend)*tan(slip_angle(i))
                  else
                     yend = ymin
                  endif
                  islp = islp + 1
                  locphi(islp) = 2
                  xslp(islp) = xstart
                  xendslp(islp) = xend
                  yslp(islp) = ystart
                  yendslp(islp) = yend
                  xstart = xstart + dxslip
!                  xend = xstart + (ymax-ystart)/tan(slip_angle(i))
c$$$                  if (mod(islp,25) .eq. 0) then 
c$$$                     write(*,fmt='(I7,1X,4(1X,E15.7))')
c$$$     $                    islp, xslp(islp), xendslp(islp), yslp(islp),
c$$$     $                    yendslp(islp)
c$$$                  endif

                  !print *, islp, xslp(islp), xendslp(islp)
               end do
               print *, 'End of second slip =', islp
               inslip(i) = islp
            endif
         else
            xstart = xmin
            xend = xmax
            ystart = dyslip
c$$$            nqslp = islp
            j = 0
            do while(ystart < ymax) 
               islp = islp + 1
               j = j+1
               locphi(islp) = 3
               xslp(islp) = xstart
               xendslp(islp) = xend
               yslp(islp) = -ystart
               yendslp(islp) = -ystart
               ystart = ystart + dyslip
c$$$                  if (mod(islp,25) .eq. 0) then 
c$$$                     write(*,fmt='(I7,1X,4(1X,E15.7))')
c$$$     $                    islp, xslp(islp), xendslp(islp), yslp(islp),
c$$$     $                    yendslp(islp)
c$$$                  endif


            end do

             print *, 'End of third slip =', islp,j
             inslip(i) = islp
         endif
      end do
      


      nslp = islp
      n_active_slip = nslp/slipspc
      print *, 'No. of active slip systems = ',n_active_slip, nslp
c$$$     Calculate shifting parameters for slip planes
c$$$      Shifting is used to control dislocations leaving upper and
c$$$     $     entering lowere
      nshift_slip = 0
      do islp = 1, nslp
         xi = xslp(islp)
         li = locphi(islp)
         if (li < 3) then 
            do jslp = 1, nslp
               xj = xslp(jslp)
               lj = locphi(jslp)
               if (lj < 3) then 
                  if (islp .ne. jslp) then 
                     if (li == lj) then 
                        if (xi == xj) then 
                           nshift_slip(islp) = jslp
                        endif
                     endif
                  endif
               endif
            enddo
         endif
c$$$         print *, 'nshift', islp, nshift_slip(islp)
      end do
      !stop
      end subroutine gen_slip_planes


      subroutine remove_disl_global(idisl_slip)
      implicit none
      include 'disl_parameters.par'
      integer :: i,j,k,l,islp, idisl_slip
      k = 0
      do islp = 1, nslp
         if (ndis(islp) > 0) then 
            do i = 1, ndis(islp)
               k = k + 1
               if (k == idisl_slip) then 
                  call rmdisl_g(i,islp)
!                  call rmdisl_continuum(i,islp)
!                  call put_in_bucket(i,islp)
! Dislocation is removed from the continuum on slip plane level and 
! global level.. 
               endif
            end do
         endif
      end do
      call assign_disloc_global
      end subroutine remove_disl_global

      subroutine put_in_bucket(i,islp)
      implicit none
      integer :: i,j,k,l,islp,jslp
      double precision :: hold
      bout(1,islp) = bout(1,islp) + b_dd(i,islp)
      sdis(i,islp) = sdis_out(1,islp)
      rdis(1,i,islp) = xslp(islp) + sdis_out(1,islp) *
     $     cosphi(locphi(islp))
      rdis(2,i,islp) = yslp(islp) + sdis_out(1,islp) *
     $     sinphi(locphi(islp))
      elem_dis(i,islp) = 0
      print *, 'Dislocation' , i, 'on', islp, 'put in bucket'
      end subroutine put_in_bucket
      

      subroutine calc_bout
      implicit none
      integer :: i,j,k,l,islp, bplus, bminus
      do islp = 1, nslp
         bplus = 0
         bminus = 0
         if (ndis(islp) > 0) then 
            do i = 1,ndis(islp)
               if (b_dd(i,islp) > 0) then 
                  bplus = bplus + 1
               else
                  bminus = bminus + 1
               endif
            end do
            bout(1,islp) = -(bplus-bminus)*blen
            bout(2,islp) = 0
         print *, 'Bucket total =', bout(1,islp)
         endif
      end do
      end subroutine calc_bout


      subroutine rmdisl(i,islp, ar1, ar2, ar3, ar4)
      implicit none
      integer :: i,j,k,l,islp,jslp
      double precision :: hold
      double precision :: ar1(mxndis,nslp), ar2(mxndis,nslp),
     $     ar3(mxndis,nslp), ar4(mxndis,nslp)
      hold = sdis(i,islp)
      do j=i+1,ndis(islp)
         sdis(j-1,islp) =sdis(j,islp)
         b_dd(j-1,islp)=b_dd(j,islp)
         rdis(1:3,j-1,islp) = rdis(1:3,j,islp)
         elem_dis(j-1,islp) = elem_dis(j,islp)
         vprev(j-1,islp) = vprev(j,islp)
         iobpin(j-1,islp) = iobpin(j,islp)
         ar1(j-1,islp) = ar1(j,islp)
         ar2(j-1,islp) = ar2(j,islp)
         ar3(j-1,islp) = ar3(j,islp)
         ar4(j-1,islp) = ar4(j,islp)
      end do
      print *, 'Removing dislocation', i, ' on ', islp
      ndis(islp) = ndis(islp)-1
      tot_disl = tot_disl -1
      call assign_disloc_global
      end subroutine rmdisl

      subroutine rmdisl_g(i,islp)
      implicit none
      integer :: i,j,k,l,islp,jslp
      double precision :: hold, roldq(3)
      double precision :: v(mxndis,nslp), c(mxndis,nslp),
     $     tau(mxndis,nslp), sold(mxndis,nslp), bb
      print *, 'Removing dislocation', i, ' on ', islp
      roldq(1:3) = rdis(1:3,i,islp)
      hold = sdis(i,islp)
      do j=i+1,ndis(islp)
         sdis(j-1,islp) =sdis(j,islp)
         b_dd(j-1,islp)=b_dd(j,islp)
         rdis(1:3,j-1,islp) = rdis(1:3,j,islp)
         elem_dis(j-1,islp) = elem_dis(j,islp)
         vprev(j-1,islp) = vprev(j,islp)
         iobpin(j-1,islp) = iobpin(j,islp)
      end do
      ndis(islp) = ndis(islp)-1
      tot_disl = tot_disl -1
      call assign_disloc_global
      
      end subroutine rmdisl_g

      subroutine assign_new_global(xtip_move)
      implicit none
      include 'disl_parameters.par'
      integer :: i,j,islp,ii,k,li, btot,fe_locate, elem_old
      double precision :: pi, xend, yend,cphi,sphi
      double precision :: xstart, ystart, bb, bsign, ss, ss1, sd
      double precision :: xtip_move
      integer :: ndisl1
      pi = 2.d0*asin(1.0d0)
      tot_disl = 0

      do islp = 1, nslp
         do i = 1, nnuc(islp)
            rnuc(1,i,islp) = rnuc(1,i,islp) - xtip_move
            elem_old = elem_source(i,islp)
            elem_source(i,islp) = fe_locate(rnuc(:,i,islp), elem_old)
         end do
      end do

      do islp = 1,nslp
         if (ndis(islp) > 0) then 
            tot_disl = tot_disl + ndis(islp)
         end if
      end do
      
      call calc_bout
      ! --- Set all quantities to zero
      do islp = 1, nslp
         li = locphi(islp)
         cphi = cosphi(li)
         sphi = sinphi(li)

         do i = 1, ndis(islp)
            sd = sdis(i,islp)
            if (li < 3) then 
               if (yendslp(islp) > 0) then 
                  rdis(1,i,islp) = xslp(islp) + sd*cphi
                  rdis(2,i,islp) = yslp(islp) + sd*sphi
               else
                  if (li == 1) then 
                     rdis(1,i,islp) = xslp(islp) - sd*cphi
                  else
                     rdis(1,i,islp) = xslp(islp) + sd*cphi
                  end if
                  rdis(2,i,islp) = -(yslp(islp) + sd*sphi)
               end if
            else
               rdis(1,i,islp) = xslp(islp) + sd*cphi
               rdis(2,i,islp) = yslp(islp) + sd*sphi   
            end if
            rold(1,i,islp) = rold(1,i,islp) - xtip_move
            rold(2:3,i,islp) = rold(2:3,i,islp)
            elem_old = elem_dis(i,islp)
            elem_dis(i,islp) = fe_locate(rdis(:,i,islp), elem_old)
         end do
      end do
      call assign_disloc_global
      
      end subroutine assign_new_global

      subroutine assign_disloc_global
      implicit none
      include 'disl_parameters.par'
      integer :: i,j,islp,ii,k,li, btot,fe_locate
      double precision :: pi, xend, yend,cphi,sphi
      double precision :: xstart, ystart, bb, bsign, ss, ss1
      integer :: ndisl1
      pi = 2.d0*asin(1.0d0)
      tot_disl = 0
      do islp = 1,nslp
         if (ndis(islp) > 0) then 
            tot_disl = tot_disl + ndis(islp)
         end if
      end do
      call calc_bout
      ! --- Set all quantities to zero
      burgers = 0.0d0
      burg_length = blen
      theta_e = 0.0d0
      theta_s = 0.0d0
      r_disl = 0.0d0
      r_old = 0.0d0 
      pk_force = 0.0d0
      pk_f = 0.0d0
      elem_disl = 0
      ndisl = 0
      ndisl_dd = 0
      ndisl1 = 0
      k = 0
      do islp = 1, nslp
         li = locphi(islp)
         sphi = sinphi(li)
         cphi = cosphi(li)
         xend = xendslp(islp)
         yend = yendslp(islp)
         xstart = xslp(islp)
         ystart = yslp(islp)
         ss = sdis_out1(1,islp) + disl_range_tol
c$$$         if (sdis_out(1,islp) .ne. sdis_out(1,islp)) then
c$$$            ss1 = sdis_out1(1,islp) -6.d0*disl_range_tol
c$$$         else
            ss1 =sdis_out(1,islp)
c$$$         end if
         if (abs(bout(1,islp)) > 0.0d0) then 
            bb = blen
            btot = nint(abs(bout(1,islp))/blen)
            if (bout(1,islp) < 0) then 
               bb = -blen
            endif
!            print *, 'Bucket dislocations'
            do i = 1, btot
               k = k + 1
               burg_length(k) = bb
               burgers(1,k) = bb*cphi
               burgers(2,k) = bb*sphi
               burgers(3,k) = 0.0d0
               if (bb > 0) then 
                  theta_e(k) = -pi
               else
                  theta_e(k) = 0.0d0
               end if
c$$$               if (bb > 0) then 
c$$$                  if (yend > 0.0d0) then 
c$$$                     theta_e(k) = pi
c$$$                  else
c$$$                     theta_e(k) = -pi
c$$$                  endif
c$$$               else
c$$$                  if (yend > 0.0d0) then 
c$$$                     theta_e(k) = 0.0d0
c$$$                  else
c$$$                     theta_e(k) = 0.0d0
c$$$                  endif
c$$$               endif
c     Try making the position of the bucket dislocation to be at the
c     edge of the atomistic region 
c     This means replacing the next 3 lines here with sdis_out1
               r_disl(1,k) = xstart + (ss1)*cphi
               r_disl(2,k) = ystart + (ss1)*sphi
               r_disl(3,k) = 0.0d0


               r_old(1:3,k) = r_disl(1:3,k)
               elem_disl(k) = 0
               disl_index(k) = 2
               disl_range(1,k) = atom_xmax + disl_range_tol
               disl_range(2,k) = atom_ymax + disl_range_tol
               print *, 'Buck', i, islp, r_disl(1,k), r_disl(2,k)
            end do
         endif

         if (ndis(islp) > 0) then 
            do i = 1, ndis(islp)
               k = k + 1
               ndisl1 = ndisl1 + 1
               ndisl_dd(k) = ndisl1
               burg_length(k) = b_dd(i,islp)
               burgers(1,k) = b_dd(i,islp)*cosphi(locphi(islp))
               burgers(2,k) = b_dd(i,islp)*sinphi(locphi(islp))
               burgers(3,k) = 0.0d0
               bb = b_dd(i,islp)
               !print *, 'qqq',k, ndisl1, ndisl_dd(k)
c     Do not switch cut plane directions in switching from bottom to top
c     This is only to ensure that the currect disl_u subroutine is
c     satisfied 
c     Also since there is always going to be a bucket dislocation 
c     The direction of the cut plane is not relevant (???)
               if (bb > 0) then 
                  theta_e(k) = -pi
               else
                  theta_e(k) = 0.0d0
               end if
c$$$               if (bb > 0) then 
c$$$                  if (yend > 0.0d0) then 
c$$$                     theta_e(k) = -pi
c$$$                  else
c$$$                     theta_e(k) = 0.0d0
c$$$                  endif
c$$$               else
c$$$                  if (yend > 0.0d0) then 
c$$$                     theta_e(k) = 0.0d0
c$$$                  else
c$$$                     theta_e(k) = -pi
c$$$                  endif
c$$$               endif

               r_disl(1:3,k) = rdis(1:3,i,islp)
               r_old(1:3,k) = rold(1:3,i,islp)
               elem_disl(k) = elem_dis(i,islp)
               disl_index(k) = 2
               disl_range(1,k) = xstart + abs(ss*cphi)
               disl_range(2,k) = abs(ystart) + abs(ss*sphi)
            end do
         endif
      end do
      ndisl = k
!      ndisl = k-1
      print *, 'Global Disloc total = ', ndisl
      print *, 'Gloabal Dislocations are'
      do i = 1,ndisl
        write(*,fmt='(I7,6(1X,E15.8),1X,2I7,3(1X,E15.8))')
     $        i, r_old(1:2,i), r_disl(1:2,i), disl_range(1:2,i),
     $        elem_disl(i), ndisl_dd(i), burgers(1:2,i), theta_e(i)
      end do
      end subroutine assign_disloc_global


      subroutine move_disloc(rhs)
      implicit none
      integer :: i,j,k,l,ii,islp, li, elem_old, fe_locate
      double precision :: si,bi,sini,cosi, tauij,taug,
     $     dsi,sout(3), Fg, Fgg,snew, roldq(3), phi2, rhs(*)
      double precision :: sinit, sfinal, tauii, tau, xi, yi
      double precision :: dist, eps,dlim, flim, dstar, s1, xstart,ystart
      double precision, allocatable :: tau1(:,:), v1(:,:), c(:,:),
     $     sold(:,:)
c$$$      double precision :: tau1(mxndis,nslp), v1(mxndis,nslp),
c$$$     $     c(mxndis,nslp), sold(mxndis,nslp)
      double precision :: cnew, ddt, dt, dx, dv, ddtm, vi, vpp,vj, xqone
      double precision :: tlce, rpile, ds, temp_s(mxndis,nslp)
      logical :: move_bottom
      include 'disl_parameters.par'

! Include obstacle parameters here
      type(DD) :: dd1(tot_disl)
      dt = tincr
      ddt = dt
      eps = 1.0d-3*blen
      flim = eps
      dlim = flim/5.0d0
      dstar = 2.d0*blen
      rpile = dstar
      
      allocate(tau1(mxndis,nslp), v1(mxndis,nslp), c(mxndis,nslp))
      allocate(sold(mxndis,nslp))
      tau1 = 0.0d0
      v1 = 0.0d0
      c = 0.0d0
      sold = 0.0d0
      if (tot_disl > 0) then 
         call assign_disloc_only(dd1)
         call calcForce(dd1,1,tot_disl)
      end if
      k = 1
      do islp = 1, nslp
         li = locphi(islp)
         cosi = cosphi(li)
         sini = sinphi(li)
         sinit = sdis_out1(1,islp)
         sfinal = sdis_out(2,islp)
         if (ndis(islp) > 0) then 
            do i = 1, ndis(islp)
               phi2=2.d0*phislp(locphi(islp))
               bi = b_dd(i,islp)
               si = sdis(i,islp)
               sold(i,islp) = si
               elem_old = elem_dis(i,islp)
               roldq(1:3) = rdis(1:3,i,islp)
               Fg = real(dd1(k)%forced)
               Fgg = real(dd1(k)%forcedg)
               if (elem_dis(i,islp) > 0) then 
                  call fe_stress(elem_dis(i,islp), rhs,
     $                 sout)
               endif
               phi2=2.d0*phislp(locphi(islp))
               taug = 0.5d0*(sout(2)-sout(1))*sin(phi2)
     $              + sout(3)*cos(phi2)
               tauij = Fg/bi
               tauii = tauij + taug
               vi = tauii*bi/Bdrag

               vpp = vi/(1-Fgg*dt/Bdrag)
               if (vi/vpp > 0) then 
                  vi = vpp
               else
                  if (abs(vi) > abs(vpp)) then 
                     vi = vi + vpp
                  endif
               endif
               if (abs(vi) > vcutoff) then 
                  vi = sign(vcutoff,vi)
               endif
               if (yendslp(islp) < 0.0d0) then 
                  vi = -vi
               end if
               tau1(i,islp) = tauii
               dsi = vi*dt
               snew = si + dsi
               if (iobpin(i,islp) == 0) then 
               write(*,fmt='(A5,1X,2I7,1X,5(1X,E15.8))') 'Disl.',
     $              i, islp, vi, si, dsi, snew, dt
            end if
! Check against obstacles
               if (iobpin(i,islp) > 0) then 
                  snew = si
               else
                  do j =1,nobs(islp)
                     dist = sobs(j,islp)-si
                     if(abs(dist) .ge. eps) then 
                        if (dsi/dist > 1.0d0) then 
                           snew = sobs(j,islp)
                           iobpin(i,islp) = -j
                           print *, 'Disloc ', i, ' on ', islp,
     $                          'about to be pinned at obstacle', j
                        endif
                     endif
                  end do
               endif
               v1(i,islp) = (snew-si)/dt
               k = k + 1
            end do
         endif
      end do
!     Determine underrelaxation parameters --- not really necessary but
!     still there 
!      call undrex(v1,c)
      print *, 'Underrelaxation factors'
      c = 1.0d0
      !v1 = vi
      do islp = 1, nslp
         if (ndis(islp) > 0) then 
 128        continue
            do i = 2, ndis(islp)
               print *, 'undrex', i,islp, c(i,islp), v1(i,islp)
               dx = sdis(i,islp)-sdis(i-1,islp)
               if (dx < 0.0d0) then 
                  write(*,*) 'Wrong order of dislocations', dx, i,islp
                  stop
               endif
               vi = c(i,islp)*v1(i,islp)
               vj = c(i-1,islp)*v1(i-1,islp)
               dv = vi - vj
               if (dv .ne. 0.0d0) then 
                  ddtm = -dx/dv
                  if (ddtm > 0.0d0) then 
                     if (b_dd(i,islp)*b_dd(i-1,islp) > 0.0d0) then 
                        xqone = 1.0d0 - 1.0d-10
                        if (dx/rpile < xqone) then 
                           print *, 'too short in pileup', dx/rpile, i
     $                          ,islp
                        endif
                        dx = dx - rpile
                        if (dx < 0.0d0 .and. tlce(abs(dx)/eps) < 1.0d0)
     $                       then 
                           dx = 0.0d0
                        endif
                        ddtm = -dx/dv
                     endif
                     ddtm = abs(ddtm)
                     if (tlce(ddtm/ddt) <= 1.0d0) then 
                        ddt = ddt/2.0d0
                        cnew = 0.8d0*ddtm/dt
                        c(i,islp)=c(i,islp)*cnew
                        c(i-1,islp) = c(i-1,islp)*cnew
                        goto 128
                     endif
                  endif
               endif
            end do
         endif
      end do
        
      do islp = 1, nslp
         li = locphi(islp)
         cosi = cosphi(li)
         sini = sinphi(li)
         if (sdis_out1(1,islp) .ne. sdis_out(1,islp)) then 
            sinit = sdis_out1(1,islp)
         else
            sinit = sdis_out(1,islp)
         endif
         sfinal = sdis_out(2,islp)
         if (ndis(islp) > 0) then 
            print *, 'Sinit = ', sinit, ' Sfinal = ',sfinal 
            do i = 1, ndis(islp)
               s1 = sdis(i,islp) + c(i,islp)*v1(i,islp)*ddt
               bi = b_dd(i,islp)
               dsi = c(i,islp)*v1(i,islp)*ddt
               sdis(i,islp) = sdis(i,islp)+dsi
               print *, 'Disloc', i, 'Old = ', sdis(i,islp),
     $         ' new = ', s1
               rdis(3,i,islp) = 0.0d0        
               snew = sdis(i,islp)
               temp_s(i,islp) = snew
!     Pin disloc at the  interface
!               if (elem_dis(i,islp) > 0) then
                     if (snew < sinit .or. snew > sfinal) then
                        if (snew < sinit) then 
                           print *,
     $                          'Dislocation close to interface 
     $                          or atomisitcs',snew
                           snew = sinit
                        else
                           snew = sfinal
                        endif
                     endif
!               endif

               sdis(i,islp) = snew
               rold(1:3,i,islp) = rdis(1:3,i,islp)
               if (li < 3) then 
                  if (yendslp(islp) > 0) then 
                     rdis(1,i,islp)=xslp(islp) + snew*cosi
                     rdis(2,i,islp)=yslp(islp) + snew*sini
                  else
                     if (li == 1) then 
                        rdis(1,i,islp)=xslp(islp)-snew
     $                    *cosi
                     else 
                        rdis(1,i,islp)=xslp(islp)+abs(snew
     $                    *cosi)
                     endif
                     rdis(2,i,islp)=-(yslp(islp)+snew
     $                    *sini)
                  endif
               else
                  rdis(1,i,islp) = xslp(islp) + snew*cosi
                  rdis(2,i,islp) = yslp(islp) + snew*sini
               endif

               elem_dis(i,islp) = fe_locate(rdis(:,i,islp),elem_old)
               print *, 'Dislocation', i, 'on', islp, 'with',
     $              bi, 'moved from',
     $              rold(1,i,islp),rold(2,i,islp),'to',
     $              rdis(1,i,islp),rdis(2,i,islp)
     $              , si, dsi, elem_dis(i,islp)


c$$$               print *, 'pk-force =', i, islp, tauii*bi,
c$$$     $              real(dd1(k)%forced)
c$$$  print *, 'Dislocation ', i, ' on ', islp, 'moved from',
c$$$  $                elem_old, ' to ', elem_dis(i,islp)
c$$$               rold(1:3,i,islp) = rdis(i:3,i,islp)
               k = k + 1
            end do
         endif
      end do

      
      do islp =1,nslp
         if (ndis(islp) > 1) then 
            i = 2
            do while (i <= ndis(islp)) 
               if (b_dd(i,islp)*b_dd(i-1,islp)<0.0d0) then 
                  if (tlce(sdis(i-1,islp)/sdis(i,islp)) > 1.0d0) then 
                     si = sdis(i,islp)
                     write(*,*) 'Wrong order', i,islp,
     $                    sdis(i,islp),sdis(i-1,islp)
                  else
                     ds = sdis(i,islp) - sdis(i-1,islp)
                     if (ds <= xLe) then 
                        call rmdisl(i,islp,v1, tau1,c,sold)
                        call rmdisl(i-1,islp,v1,tau1,c,sold)
                        i = i -1
                     endif
                  endif
               endif
               i = i + 1
            end do
         endif
      end do

      do islp = 1, nslp
         do i=1,ndis(islp)
            do j = 1, nobs(islp)
               if (iobpin(i,islp)+j .eq. 0) then 
                  dist = dabs(sobs(j,islp)-sdis(i,islp))
                  if (dist .lt. eps) then 
                     iobpin(i,islp) = j
                     write(*,*) 'Dislocation ', i, ' on ', islp,
     $                    ' pinned at obstacle ', j
                  else
                     iobpin(i,islp) = 0
                  endif
               endif
            end do
         end do
      end do
      

      call rlsdis(tau1)

!     --- Pass dislocations from top to bottom or simply put in bucket
!     if xstart < atom_xmin + tol
      move_bottom = .true.
      do islp = 1,nslp
         li = locphi(islp)
         cosi = cosphi(li)
         sini = sinphi(li)
         xstart = xslp(islp)
         ystart = yslp(islp)
         if (sdis_out1(1,islp) .ne. sdis_out(1,islp)) then 
            move_bottom = .false. 
            sinit = sdis_out(1,islp)
         else
            sinit = sdis_out(1,islp)
            move_bottom = .true.
         endif
         sfinal = sdis_out(2,islp)
         i = 1
         k = 1
         do while(i<ndis(islp)) 
            if (ndis(islp) > 0) then 
               if (temp_s(i,islp) <  sinit) then 
                  if (move_bottom) then 
                     l = islp + nshift_slip(islp)
                     print *, 'Removing from top and putting in bottom',
     $                    temp_s(i,islp), sinit
                     if (xslp(islp) < atom_xmin) then 
                        call rmdisl_g(i,islp)                        
                     else
                        call crdisl(temp_s(i,islp), sinit, l, b_dd(i
     $                       ,islp),1)
                        call rmdisl_g(i,islp)
                     endif
                  else
                     i = i + 1
                  endif
               else
                  i = i + 1
               endif
               k = k + 1
            endif
         end do
      end do


      do islp = 1, nslp
         if (ndis(islp) > 0) then 
            do i = 1, ndis(islp)
               vprev(i,islp) = v1(i,islp)
            end do
         endif
      end do

      call assign_disloc_global
c$$$      deallocate(tau1,v1,c)
c$$$      deallocate(sold)
      return
      end subroutine move_disloc


      

      subroutine nucleate_disl(rhs)
      implicit none
      integer islp,jslp,i,j,k,l
      double precision s, xi, xLh, bsign, sout(3), rhs(*), phi2,tau
      double precision ev_convert, fact, tau_source, tauij_s, taug_s
      double precision qsign
      type(DD) :: dd1(tot_size)
      logical create
      bsign = 1.0d0
      create = .false.
            ev_convert = 1.602176462
            fact = 1.d0/ev_convert/1.d-5
      call assign_disloc_only(dd1)
      call assign_source_only(dd1)
      call calcForce(dd1,tot_disl+1,tot_size)
      k = tot_disl + 1
      do islp = 1, nslp
         if (nnuc(islp) > 0) then 
            do i = 1, nnuc(islp)
               if (elem_source(i,islp) > 0) then 
                    call fe_stress(elem_source(i,islp), rhs,
     $                   sout)
                 endif
                 phi2=2.d0*phislp(locphi(islp))
                 taug_s = 0.5d0*(sout(2)-sout(1))*sin(phi2)
     $                + sout(3)*cos(phi2)
                 tauij_s = real(dd1(k)%forced)/blen
                 tau_source = taug_s + tauij_s
c$$$                 taui(i,islp) = taui(i,islp) + real(dd1(k)%forced)/blen
c$$$     $                + tau 

                 sout = sout*fact
               write(*,fmt='(3I7,1X,3(1X,E15.7))') i,islp,k,tau_source,
     $                real(dd1(k)%forced)/blen*fact, taui(i,islp)*fact


               if (tnlaps(i,islp) < 0.5*tincr) then 
                  tnlaps(i,islp) = 0.0d0
               end if
               qsign = 1.0d0
            
               if (tau_source*taui(i,islp) < 0.0d0) then 
                  qsign = -1.0d0
               end if
            
               if (abs(tau_source) >= t_FR(i,islp)) then 
                  tnlaps(i,islp) = tnlaps(i,islp) + qsign*tincr
                  print *, 'Source ', i, ' on ', islp, 'tnlaps',
     $                 tnlaps(i,islp)
               else
                  if (tnlaps(i,islp) > 0) then 
                     tnlaps(i,islp) = tnlaps(i,islp) - tincr
                  end if
               end if
               if (tnlaps(i,islp) < 0.5*tincr) then 
                  tnlaps(i,islp) = 0.0d0
               end if
               k = k + 1
               taui(i,islp) = tau_source
            end do
         endif
      end do
      do islp = 1, nslp
         if (nnuc(islp) > 0) then 
            do i = 1, nnuc(islp)
               if (tnlaps(i,islp) > tnuc) then
                  xLh = xLnuc(i,islp)/2.d0
                  xi = snuc(i,islp)
                  if (taui(i,islp) < 0.0d0) bsign = -1.0d0
                  if (yendslp(islp) < 0.0d0) then 
                     bsign = -bsign
                  end if
                  !print *, 'Generating dislocations on slip plane',islp
                  call crdisl(xi-xLh, xi,islp,-blen*bsign,i)
                  call crdisl(xi+xLh, xi,islp,blen*bsign,i)
                  create = .true.
                  tnlaps(i,islp) = 0.0d0
               endif
            end do
         endif
      enddo
      !if (create) stop
      call assign_disloc_global
      end subroutine nucleate_disl

      subroutine nucleate_atomistic_disl (islp, s_disl, th_e, th_s, bv)
c     Nucleate a dislocation passed in from atomistics
c     Automatically handles dislocations to be put into the bucket 

c     islp -- Slip plane to nucleate dislocation
c     s_disl --- Position on slip plane to nucleate dislocation
c     bv -- Burgers vector of dislocation
c     th_e, th_s --- Angles of dislocation to be inserted
      implicit none
      integer :: islp, i, j, k, l, li
      double precision :: s_disl, th_e, th_s, bv(3)
      double precision :: bsign, bmag
      double precision, parameter :: dtol=1.d-3
c     Angle is already checked when finding slip plane
c     Now we only re-check with different criteria

c     Assign correct sign for dislocation 

      li = locphi(islp)
      if (abs(bv(1)-blen*cosphi(li))< dtol .and.
     $     abs(bv(2) - blen*sinphi(li))< dtol) then 
         bsign = 1.0d0
      else
         bsign = -1.0d0
      end if


      call crdisl(s_disl, s_disl-6.d0*blen, islp, bsign*blen)
      call assign_disloc_global 
      end subroutine nucleate_atomistic_disl

      subroutine crdisl (si,slim,islp,bi,nuc_num)
      implicit none
      integer i,j,k,l,m,n,islp, fe_locate, el_source
      integer, optional :: nuc_num
      double precision si,slim,bi,shold,xcl
      double precision xint, dstr,dstar,sold, stest
      dstar = 2*blen
      n = ndis(islp)
      shold = si
      xcl = si-slim
      xint = 3.d0*2.d0*blen
      dstr = 2.d0*dstar
      rpile = 2.d0*blen
      i = 0
      if (present(nuc_num)) then 
         el_source = elem_source(nuc_num,islp)
      end if
      ! --- Make sure dislocations are in order
cTODO ! -- Make sure dislocations are atleast 6*b away from each other 
      i = 1
      do while (shold > sdis(i,islp) .and. i <= ndis(islp))
         i = i + 1
      end do
      k = i
      !print *, 'Insertion location', k
      ndis(islp) = ndis(islp) + 1
      tot_disl = tot_disl + 1
      tot_size = tot_size + 1
      do i=ndis(islp),k+1,-1
         sdis(i,islp) = sdis(i-1,islp)
         b_dd(i,islp)=b_dd(i-1,islp)
         rdis(1:3,i,islp) = rdis(1:3,i-1,islp)
         elem_dis(i,islp)=elem_dis(i-1,islp)
         rold(1:3,i,islp) = rold(1:3,i-1,islp)
         vprev(i,islp)=vprev(i-1,islp)
         iobpin(i,islp) = iobpin(i-1,islp)
      end do
      
      sdis(k,islp) = shold
      b_dd(k,islp) = bi
      iobpin(k,islp) = 0
      if (yendslp(islp) > 0) then 
         rdis(1,k,islp)=xslp(islp)+si*cosphi(locphi(islp))
         rdis(2,k,islp)=yslp(islp)+si*sinphi(locphi(islp))
      else
         if (locphi(islp) == 1) then 
            rdis(1,k,islp)=xslp(islp)-si*cosphi(locphi(islp))
         else
            rdis(1,k,islp)=xslp(islp)+abs(si*cosphi(locphi(islp)))
         end if
         rdis(2,k,islp)=-(yslp(islp)+si*sinphi(locphi(islp)))
      endif
      rdis(3,k,islp) = 0.0d0
      elem_dis(k,islp) = fe_locate(rdis(1:3,k,islp),1)
      write(*, fmt='(A30,1X,I7,3(1X,E15.7),1X, A15,I7,1X,E15.8)')
     $     'Dislocation Generated on'
     $     , islp, sdis(k,islp),rdis(1,k,islp), rdis(2,k,islp),
     $     'in element', elem_dis(k,islp), b_dd(k,islp)
      end subroutine crdisl

      subroutine gen_slip_ends
      implicit none
      integer :: i,j,k,l,islp, li
      double precision :: xstart, xend, ystart, yend,xend1,yend1
      double precision :: sphi, cphi, tphi
      double precision :: dlmax, dlmax1
      double precision :: xmax1, xmin1, ymin1, ymax1

c     Initialize the dislocation end variables 

      sdis_out = 0.0d0
      sdis_out1 = 0.0d0

c$$$      xmax1 = atom_xmax + sdis_out_tol
c$$$      xmin1 = (atom_xmin - sdis_out_tol)/2.0d0
c$$$      ymax1 = atom_ymax + sdis_out_tol
c$$$      ymin1 = atom_ymin - sdis_out_tol

      xmax1 = atom_xmax 
      xmin1 = atom_xmin
      ymax1 = atom_ymax
      ymin1 = atom_ymin

      do islp = 1,nslp
         li = locphi(islp)
         cphi = cosphi(li)
         sphi = sinphi(li)
         tphi = sphi/cphi
         xstart = xslp(islp)
         ystart = yslp(islp)
         xend = xendslp(islp)
         yend = yendslp(islp)
         if (li .ne. 3) then 
            dlmax = abs(yend-ystart)/sphi
         else
            dlmax = xend-xstart
         endif
         sdis_out(1,islp) = 0.0d0 + 1.d-3 
         sdis_out(2,islp) = dlmax - 1.d-3
         if (xstart > xmin1 .and. xstart < xmax1) then 
            if (li .ne. 3) then 
               if (yend > 0) then 
                  xend1 = xstart + (ymax1-ystart)/tphi
                  yend1 = atom_ymax
                  if (xend1 > xmax1) then 
                     xend1 = xmax1
                     yend1 = abs(xend1-xstart)*abs(tphi)
                  endif
                  if (xend1 < xmin1) then 
                     xend1 = xmin1
                     yend1 = abs(xend1-xmin1)*abs(tphi)
                  endif
               endif
               dlmax1 = sqrt((xstart-xend1)**2 + (ystart-yend1)**2) +
     $              sdis_out_tol
               sdis_out1(1,islp) = dlmax1 - 2.0d0
               write(*,fmt='(2I7,5(1X,E15.8))'),
     $              li, islp, xslp(islp), xend1, yend1,sdis_out1(1,islp)
            endif
         endif
      end do
      end subroutine gen_slip_ends

      subroutine gen_sources
      implicit none
      double precision :: xe, xnu, xlambda, xmu, pi
      double precision :: cc(6,6)
      integer :: i_elas
      common /elastic/ xe, xnu, xlambda, xmu, cc, i_elas
      integer :: i,j,k,l, islp, iseed1, iseed2, fe_locate
      logical :: accept
      double precision :: ev_convert, e, mu, nu, lambda, lnuc
      double precision, allocatable :: source_dist(:)
      ! --- Obtain Properties of the material from main input
      !---- Calculate total no. of sources based on density and area
      ! --- Generate an array of source strengths
      ! --- Assign each source the strength and compute L_nuc
      ! --- Place the source according to L_nuc

      double precision :: dlmax, dlengthslp, xstart, xend, ystart,yend
      double precision :: qsrc, rand(2), rand1(2), sphi, cphi
      integer :: nsrc, ii, li
      double precision :: atomistic_exclusion, lnuc_exclusion, factor,
     $     tol



!     currently hard coded -- need to include in input file      
      pi = 3.1415926535898
      ev_convert = 1.602176462
      source_den = 66 ! per /um^2
      source_den = 66.0d-8  ! per Angstrom^2
      tot_source = int(process_area*source_den) + 10
      qsrc = tot_source*1.d0/(n_active_slip*1.d0) ! -- no. of sources per slip plane
      print *, 'QSRC =', qsrc, tot_source



!      tot_source = 1000
      avg_source_str = 50 ! MPa
      sd_source_str = 5 !MPa
      avg_source_str = avg_source_str*ev_convert*1.d-5
      sd_source_str = sd_source_str*ev_convert*1.d-5
      print *, 'Source strength =', avg_source_str, xmu, xnu
      print *, 'Source Strength MPa', avg_source_str/ev_convert/1.d-5,
     $     xmu/ev_convert/1.d-2, xe/ev_convert/1.d-2
      factor = xmu*blen/(2.d0*pi*(1.d0-xnu))
      lnuc = factor/avg_source_str
      tol = lnuc/4.d0
      print *, 'Lnuc =', lnuc, blen, tot_source
      iseed1 = 1234
      iseed2 = 567
      
      call rmarin(iseed1,iseed2)
      allocate(source_dist(1000))
      call rgauss(source_dist, 1000,avg_source_str,sd_source_str)
c$$$      do i = 1, tot_source
c$$$         write(*,fmt='(I7,1X,2E15.8)'), i, source_dist(i),
c$$$     $        avg_source_str
c$$$      end do
      ! Generate sources
      ! For now assume that qsrc always less than 1
c$$$      if (qsrc > 1.0d0) then 
c$$$         nsrc = int(qsrc*10.d0)/10
c$$$         qsrc = qsrc-nsrc
c$$$      endif
      
!      atom_exclusion = (atom_xmax-atom_xmin)/tan(phislp(1))
      allocate(nnuc(nslp))
      allocate(snuc(mxnnuc,nslp))
      allocate(rnuc(3,mxnnuc,nslp))
      allocate(t_FR(mxnnuc,nslp))
      allocate(xLnuc(mxnnuc,nslp))
      allocate(tnlaps(mxnnuc,nslp))
      allocate(taui(mxnnuc,nslp))
      allocate(elem_source(mxnnuc,nslp))
!      allocate(nuc_pk_stress(3,1,nslp))
      
      allocate(ndis(nslp))
      allocate(sdis(mxndis,nslp))
      allocate(b_dd(mxndis,nslp))
      allocate(vprev(mxndis,nslp))
      allocate(iobpin(mxndis,nslp))
      allocate(rdis(3,mxndis,nslp))
      allocate(rold(3,mxndis,nslp))
      allocate(elem_dis(mxndis,nslp))
      allocate(bout(2,nslp))
      allocate(sdis_out(2,nslp))
      allocate(sdis_out1(2,nslp))

      allocate(nobs(nslp))
      allocate(sobs(mxnobs,nslp))
      allocate(tau_obs(mxnobs,nslp))
      


      ndis = 0
      sdis = 0.0d0
      b_dd = 0.0d0
      vprev = 0.0d0
      iobpin = 0
      rdis = 0.0d0
      rold = 0.0d0
      elem_dis = 0

      nnuc = 0
      rnuc = 0.0d0
      snuc = 0.0d0
      t_FR = 0.0d0
      xLnuc = 0.0d0
      tnlaps = 0.0d0
      taui = 0.0d0
      bout = 0.0d0
      sdis_out = 0.0d0
      sdis_out1 = 0.0d0

      nobs = 0
      sobs = 0.0d0
      tau_obs = 0.0d0
      ii = 0
      print *, 'QSRC = ', qsrc

      do islp = 1,nslp
         li=locphi(islp)
         cphi = cosphi(locphi(islp))
         sphi = sinphi(locphi(islp))
         xstart = xslp(islp)
         ystart = yslp(islp)
         xend = xendslp(islp)
         yend = yendslp(islp)
         if (li < 3) then 
            dlmax = abs(yend-ystart)/abs(sphi)
         else
            dlmax = xend-xstart
         endif
c$$$         write(*,fmt='(A20,1X,2I7,5(1X,E15.8))'), 'Slip plane',
c$$$     $        islp, li, xstart, xend, ystart, yend, dlmax
      end do
      print *, 'Source generation'

      islp = 2
      print *, 'Total slip planes = ', nslp
c$$$      do while (islp <= nslp)
      do islp = 1, nslp
         li = locphi(islp)
         cphi = cosphi(li)
         sphi = sinphi(li)
         xstart = xslp(islp)
         ystart = yslp(islp)
         xend = xendslp(islp)
         yend = yendslp(islp)
         if (li < 3) then 
            dlmax = abs(yend-ystart)
     $           /abs(sphi)
         else
            dlmax = xend -xstart
         endif
         call ranmar(rand,2)

  !        Hard code one source on +60 slip plane a little distacnce   away
c$$$         if (locphi(islp) .ne. 3) then 
c$$$            if (xslp(islp) > 5.d0*blen-1.d0 .and.
c$$$     $           xslp(islp) < 5.d0*blen+1.d0) then 
c$$$                  if (yend > 0.0d0) then 
c$$$                     if (yendslp(islp) > 0.0d0) then 
c$$$                        if (locphi(islp) == 1) then 
c$$$                           dlengthslp = lnuc + 100.0d0 + rand(1)*lnuc
c$$$                        else
c$$$                           dlengthslp = 2.d0*lnuc + 100.0d0 + rand(1)
c$$$     $                          *lnuc
c$$$                        end if
c$$$                     else
c$$$                        dlengthslp = 4.d0*lnuc + 100.0d0 + rand(1)*lnuc
c$$$                     end if
c$$$                     if (locphi(islp) == 1) then 
c$$$                        if (yendslp(islp) < 0.0d0) then 
c$$$                           cycle
c$$$                        end if
c$$$                     end if
c$$$                     k = 1
c$$$                     nnuc(islp) = k
c$$$                     t_FR(k,islp) = avg_source_str
c$$$                     xLnuc(k,islp) = lnuc
c$$$                     ii = ii + 1
c$$$                     snuc(k,islp) = dlengthslp 
c$$$                     if (yend > 0) then 
c$$$                        rnuc(1,k,islp) = xstart + dlengthslp*cphi
c$$$                        rnuc(2,k,islp) = ystart + dlengthslp*sphi
c$$$                        rnuc(3,k,islp) = 0.0
c$$$                     else
c$$$                        if (locphi(islp) ==1) then 
c$$$                           rnuc(1,k,islp) = xstart - dlengthslp*cphi
c$$$                        else
c$$$                           rnuc(1,k,islp) = xstart+abs(dlengthslp*cphi)
c$$$                        endif
c$$$                        rnuc(2,k,islp) = -(ystart + dlengthslp*sphi)
c$$$                        rnuc(3,k,islp) = 0.0
c$$$                     endif
c$$$                     elem_source(k,islp) = fe_locate(rnuc(:,k,islp),1)
c$$$                     write( *,fmt='(A10,1X,3(1x,E15.7))') 'Coords are',
c$$$     $                    snuc(k,islp),rnuc(1,k,islp), rnuc(2,k,islp)
c$$$                     print *, 'Located in Element ', elem_source(k,islp)
c$$$                     print *, 'Islp = ', islp, locphi(islp), cphi, sphi
c$$$c$$$                     exit
c$$$                  endif 
c$$$            else
c$$$               goto 10
c$$$            endif
c$$$         else 
c$$$            goto 10
c$$$         endif


c$$$!     Generate source randomly from gaussian distribution
c$$$  
c$$$         if (rand(1) <= qsrc) then 
c$$$            accept = .false.
c$$$            k = 1
c$$$            ii = ii + 1
c$$$            if (yend > 0.0d0) then 
c$$$               t_FR(k,islp) = source_dist(ii)
c$$$            else
c$$$               t_FR(k,islp) = source_dist(ii)*10
c$$$            endif
c$$$            xLnuc(k,islp) = factor/t_FR(k,islp)
c$$$            if (dlmax < 2.d0*xLnuc(k,islp)) then 
c$$$c$$$  print *, 'Rejecting slip plane', ii, islp, li,
c$$$c$$$  $              dlmax,(2.d0*xLnuc(k,islp)), t_FR(k,islp)
c$$$               t_FR(k,islp) = 0.0d0
c$$$               xLnuc(k,islp) = 0.0d0
c$$$               ii = ii -1
c$$$               islp = islp + 1
c$$$               goto 10
c$$$            endif
c$$$            nnuc(islp) = k
c$$$
c$$$            do while (accept == .false.)
c$$$               call ranmar(rand1,2)
c$$$               dlengthslp = rand1(1)*(dlmax)
c$$$               if (locphi(islp) == 3) then 
c$$$                  if (abs(ystart) > atom_ymax + xLnuc(k,islp)*2.0d0)
c$$$     $                 then 
c$$$                     if (xstart + dlengthslp < atom_xmin .or.
c$$$     $                    xstart + dlengthslp > atom_xmax) then
c$$$                        accept = .true.
c$$$                     endif
c$$$                  endif
c$$$               endif
c$$$               if (dlengthslp > xLnuc(k,islp)/2.d0 + tol/2.d0) then
c$$$                  if (dlengthslp < dlmax - xLnuc(k,islp)/2.d0-tol/2.d0)
c$$$     $                 then 
c$$$                     accept = .true.
c$$$                  endif
c$$$               endif
c$$$            end do
c$$$            snuc(k,islp) = dlengthslp
c$$$            print *, 'Source no. generated = ', ii, ' on ', islp,
c$$$     $           li, t_FR(k,islp)/ev_convert/1.d-5, snuc(k,islp)
c$$$
c$$$
c$$$            if ( locphi(islp) .ne. 3) then 
c$$$               if (yend > 0) then 
c$$$                  rnuc(1,k,islp) = xstart + dlengthslp*cphi
c$$$                  rnuc(2,k,islp) = ystart + dlengthslp*sphi
c$$$                  rnuc(3,k,islp) = 0.0
c$$$               else
c$$$                  if (locphi(islp) == 1) then 
c$$$                     rnuc(1,k,islp) = xstart - dlengthslp*cphi
c$$$                  else
c$$$                     rnuc(1,k,islp) = xstart + abs(dlengthslp*cphi)
c$$$                  endif
c$$$                  rnuc(2,k,islp) = -(ystart + dlengthslp*sphi)
c$$$                  rnuc(3,k,islp) = 0.0
c$$$               endif
c$$$            else
c$$$               rnuc(1,k,islp) = xstart + dlengthslp*cphi
c$$$               rnuc(2,k,islp) = ystart
c$$$            endif
c$$$            
c$$$            elem_source(k,islp) = fe_locate(rnuc(:,k,islp),1)
c$$$            write( *,fmt='(A10,1X,3(1x,E15.7))') 'Coords are',
c$$$     $           snuc(k,islp),rnuc(1,k,islp), rnuc(2,k,islp)
c$$$            print *, 'Located in Element ', elem_source(k,islp)
c$$$            islp = islp + 100
c$$$         else
c$$$            islp = islp + 100
c$$$         endif
 10      continue
      end do

      do islp = 1, nslp
         if (nnuc(islp) > 0) then
            do i = 1, nnuc(islp)
               print *, 'Source = ',i, islp, elem_source(i,islp)
            end do
         endif
      end do
      print *, 'Total no. of sources = ', ii

      nntot = ii
      tot_source = ii
      tot_size = tot_source

       call gen_slip_ends

c$$$      stop
      end subroutine gen_sources
      
      subroutine plot_slip
      implicit none
      integer :: i,j,k,l,li,islp,ifile
      integer :: kk(3)
      double precision :: x, y, xstart,xend, ystart, yend, sn,x1,y1
      double precision :: cphi,sphi,tphi
      character*80 filename
      integer :: tot_obs, tot_source, tot_pts
      
c$$$      filename='slip_lines.plt'
c$$$      ifile = 999
c$$$      open(unit=ifile,file=filename,status='unknown')
c$$$      write(999,*) 'Variables = X Y'
c$$$      do islp = 1,nslp
c$$$         if (nnuc(islp) > 0) then
c$$$            li =locphi(islp)
c$$$            cphi = cosphi(li)
c$$$            sphi = sinphi(li)
c$$$            tphi = sphi/cphi
c$$$            xstart = xslp(islp)
c$$$            xend = xendslp(islp)
c$$$            ystart = yslp(islp)
c$$$            yend = yendslp(islp)
c$$$            write(999,*) "Zone T = line"
c$$$            write(999,fmt='(2(1X,E18.11))') xstart, ystart
c$$$            write(999,fmt='(2(1X,E18.11))') xend, yend
c$$$         endif
c$$$      end do
c$$$      close(ifile)
      filename = ''
      filename = 'source_obs.vtk'
      tot_obs = 0
      tot_source = 0
      tot_pts = 0
      do islp = 1, nslp
         tot_source = tot_source + nnuc(islp)
         tot_obs = tot_obs + nobs(islp)
         tot_pts = tot_pts + nnuc(islp) + nobs(islp)
      end do
      
      open(unit=ifile,file=filename,status='unknown')
      write(ifile,fmt='(A)') '# vtk DataFile Version 2.0'
      write(ifile,fmt='(A)') 'Strains from CADD'
      write(ifile,fmt='(A)') 'ASCII'
      write(ifile,fmt='(A)') 'DATASET UNSTRUCTURED_GRID'
      write(ifile,fmt='(A6,1x,I7,1x,A5)') 'POINTS',tot_pts,'float'
      do islp = 1, nslp
         li = locphi(islp)
         cphi = cosphi(li)
         sphi = sinphi(li)
         yend = yendslp(islp)
!     Sources are written first
         do i = 1, nnuc(islp)
            sn = snuc(i,islp)
            if (yend > 0.0d0) then 
               x1 = xslp(islp) + sn*cphi
               y1 = yslp(islp) + sn*sphi
            else
               if (li == 1) then 
                  x1 = xslp(islp) - sn*cphi
               else
                  x1 = xslp(islp) + abs(sn*cphi)
               end if
               y1 = -(yslp(islp) + sn*sphi)
            end if
            write(ifile,'(3(1X,E14.6))') x1, y1, 0.0d0
         end do
!     Obstacles are written next
         do i = 1, nobs(islp)
            sn = sobs(i,islp)
            if (yend > 0.0d0) then 
               x1 = xslp(islp) + sn*cphi
               y1 = yslp(islp) + sn*sphi
            else
               if (li == 1) then 
                  x1 = xslp(islp) - sn*cphi
               else
                  x1 = xslp(islp) + abs(sn*cphi)
               end if
               y1 = -(yslp(islp) + sn*sphi)
            end if
            write(ifile,'(3(1X,E14.6))') x1, y1, 0.0d0          
         end do
      end do
      write(ifile,*)
      write(ifile,'(A5,1X,I7,1X,I7)') 'CELLS', tot_pts, 2*tot_pts
      do i =1, tot_pts
         write(ifile,'(2(1X,I7))') 1, i-1
      end do
      write(ifile,*)
      write(ifile,'(A10,1X,I7)') 'CELL_TYPES', tot_pts
      do i = 1, tot_pts
         write(ifile,fmt='(5(1x,I7))') 1
      end do
      write(ifile, *) 'POINT_DATA', tot_pts
      write(ifile, *) 'SCALARS SO integer 1'
      write(ifile, *) 'LOOKUP_TABLE default'
      do islp = 1, nslp
         do i = 1, nnuc(islp)
            write(ifile,'(I7)') 1
         end do
         do i = 1, nobs(islp)
            write(ifile,'(I7)') 2
         end do
      end do

      
c$$$      write(999,*) 'Variables = X Y'
c$$$
c$$$      kk = 0
c$$$      write(ifile, *) 'Zone T = nuc1'
c$$$      do islp = 1,nslp
c$$$         if (nnuc(islp) > 0) then 
c$$$            li =locphi(islp)
c$$$            if (li .eq. 1) then 
c$$$               cphi = cosphi(li)
c$$$               sphi = sinphi(li)
c$$$               tphi = sphi/cphi
c$$$               xstart = xslp(islp)
c$$$               xend = xendslp(islp)
c$$$               ystart = yslp(islp)
c$$$               yend = yendslp(islp)
c$$$               do i = 1,nnuc(islp)
c$$$                  si = snuc(i,islp)
c$$$                  if ( yend > 0) then 
c$$$                     x = xstart + cphi*si
c$$$                     y = ystart + sphi*si
c$$$                  else
c$$$                     x = xstart - cphi*si
c$$$                     y = ystart - sphi*si
c$$$                  endif
c$$$                  write(ifile,fmt='(2(1X,E18.11))') x,y
c$$$               end do
c$$$            endif
c$$$         endif
c$$$      end do
c$$$      write(ifile, *) 'Zone T = nuc2'
c$$$
c$$$      do islp = 1,nslp
c$$$         if (nnuc(islp) > 0) then 
c$$$            li =locphi(islp)
c$$$            if (li .eq. 2) then 
c$$$               cphi = cosphi(li)
c$$$               sphi = sinphi(li)
c$$$               tphi = sphi/cphi
c$$$               xstart = xslp(islp)
c$$$               xend = xendslp(islp)
c$$$               ystart = yslp(islp)
c$$$               yend = yendslp(islp)
c$$$               do i = 1,nnuc(islp)
c$$$                  si = snuc(i,islp)
c$$$                  if ( yend > 0) then 
c$$$                     x = xstart + cphi*si
c$$$                     y = ystart + sphi*si
c$$$                  else
c$$$                     x = xstart - cphi*si
c$$$                     y = ystart - sphi*si
c$$$                  endif
c$$$                  write(ifile,fmt='(2(1X,E18.11))') x,y
c$$$               end do
c$$$            endif
c$$$         endif
c$$$      end do
c$$$      write(ifile, *) 'Zone T = nuc3'
c$$$
c$$$      do islp = 1,nslp
c$$$         if (nnuc(islp) > 0) then 
c$$$            li =locphi(islp)
c$$$            if (li .eq. 3) then 
c$$$               cphi = cosphi(li)
c$$$               sphi = sinphi(li)
c$$$               tphi = sphi/cphi
c$$$               xstart = xslp(islp)
c$$$               xend = xendslp(islp)
c$$$               ystart = yslp(islp)
c$$$               yend = yendslp(islp)
c$$$               do i = 1,nnuc(islp)
c$$$                  si = snuc(i,islp)
c$$$                     x = xstart + cphi*si
c$$$                     y = ystart + sphi*si
c$$$                  write(ifile,fmt='(2(1X,E18.11))') x,y
c$$$               end do
c$$$            endif
c$$$         endif
c$$$      end do
c$$$
c$$$
c$$$      write(ifile, *) 'Zone T = obs1'
c$$$      do islp = 1,nslp
c$$$         if (nobs(islp) > 0) then 
c$$$            li =locphi(islp)
c$$$            if (li .eq. 1) then 
c$$$               cphi = cosphi(li)
c$$$               sphi = sinphi(li)
c$$$               tphi = sphi/cphi
c$$$               xstart = xslp(islp)
c$$$               xend = xendslp(islp)
c$$$               ystart = yslp(islp)
c$$$               yend = yendslp(islp)
c$$$               do i = 1,nobs(islp)
c$$$                  si = sobs(i,islp)
c$$$                  if ( yend > 0) then 
c$$$                     x = xstart + cphi*si
c$$$                     y = ystart + sphi*si
c$$$                  else
c$$$                     x = xstart - cphi*si
c$$$                     y = ystart - sphi*si
c$$$                  endif
c$$$                  write(ifile,fmt='(2(1X,E18.11))') x,y
c$$$               end do
c$$$            endif
c$$$         endif
c$$$      end do
c$$$      write(ifile, *) 'Zone T = obs2'
c$$$
c$$$      do islp = 1,nslp
c$$$         if (nobs(islp) > 0) then 
c$$$            li =locphi(islp)
c$$$            if (li .eq. 2) then 
c$$$               cphi = cosphi(li)
c$$$               sphi = sinphi(li)
c$$$               tphi = sphi/cphi
c$$$               xstart = xslp(islp)
c$$$               xend = xendslp(islp)
c$$$               ystart = yslp(islp)
c$$$               yend = yendslp(islp)
c$$$               do i = 1,nobs(islp)
c$$$                  si = sobs(i,islp)
c$$$                  if ( yend > 0) then 
c$$$                     x = xstart + cphi*si
c$$$                     y = ystart + sphi*si
c$$$                  else
c$$$                     x = xstart - cphi*si
c$$$                     y = -(ystart + sphi*si)
c$$$                  endif
c$$$                  write(ifile,fmt='(2(1X,E18.11))') x,y
c$$$               end do
c$$$            endif
c$$$         endif
c$$$      end do
c$$$      write(ifile, *) 'Zone T = obs3'
c$$$
c$$$      do islp = 1,nslp
c$$$         if (nobs(islp) > 0) then 
c$$$            li =locphi(islp)
c$$$            if (li .eq. 3) then 
c$$$               cphi = cosphi(li)
c$$$               sphi = sinphi(li)
c$$$               tphi = sphi/cphi
c$$$               xstart = xslp(islp)
c$$$               xend = xendslp(islp)
c$$$               ystart = yslp(islp)
c$$$               yend = yendslp(islp)
c$$$               do i = 1,nobs(islp)
c$$$                  si = sobs(i,islp)
c$$$                     x = xstart + cphi*si
c$$$                     y = ystart + sphi*si
c$$$                  write(ifile,fmt='(2(1X,E18.11))') x,y
c$$$               end do
c$$$            endif
c$$$         endif
c$$$      end do








      close(ifile)
      end subroutine plot_slip

      subroutine gen_obstacles
      implicit none
      double precision :: xe, xnu, xlambda, xmu, pi
      double precision :: cc(6,6)
      integer :: i_elas, inuc
      common /elastic/ xe, xnu, xlambda, xmu, cc, i_elas
      integer :: i,j,k,l,li, ii, islp, iseed1, iseed2, fe_locate
      logical :: accept
      double precision :: ev_convert, e, mu, nu, lambda, lnuc
      double precision :: tobs, dlmax, dlengthslp
      double precision :: x, x1, rand(2)
      double precision :: sphi, cphi, xstart, ystart,xend,yend
      double precision :: x_lobs_nuc, xhold, xtemp
      integer :: nobs_slip

      pi = 3.1415926535898d0
      ev_convert = 1.602176462d0
      tobs = 1500.0d0 ! MPa
      tobs = tobs*ev_convert*1.d-5
      lobs = 2000.0d0 ! Angstrom
!     take care the lobs > 2*Lnuc
      lobs_max = 1.5d0*lobs
      lobs_min = 0.5d0*lobs
      ii = 0
      iseed1 = 111
      iseed2 = 222
      call rmarin(iseed1,iseed2)

      do islp = 1,nslp
         i = 0
         li = locphi(islp)
         cphi = cosphi(li)
         sphi = sinphi(li)
         xstart = xslp(islp)
         ystart = yslp(islp)
         xend = xendslp(islp)
         yend = yendslp(islp)
         if (locphi(islp) < 3) then 
            dlmax = abs(yend-ystart)/abs(sphi)
         else
            dlmax = xend -xstart
         endif
! ------ Put special conditions for slip-plane 3  
         if (li .ne. 3) then 
            x = (atom_ymax+15.0d0)/sphi
         else
            x = 0.0d0
         endif
         if (nnuc(islp) > 0) then 
            if (dlmax > 2.d0*lobs) then 
               do while (x < dlmax)
                  !call ranmar(rand,2)
                  accept = .false. 
                  k = 0
                  do while (accept == .false.)
                      call ranmar(rand,2)
                     if (i > 0) then 
                        x1 = lobs + lobs*(rand(1)-0.5)
                     else
                        x1 = lobs_min*rand(1)
                     endif
                     x = x + x1
                     accept = .true.
                     call xmin1(x,snuc(1:nnuc(islp),islp),nnuc(islp),
     $                    xhold,inuc)
                     x_lobs_nuc = max(lobs/3.0d0, xLnuc(inuc,islp))
c$$$                     write( *,fmt='(A20,2I7,5(1X,E15.9))')
c$$$     $                    'obstacle trial ', islp,k, x,
c$$$     $                    xhold,snuc(inuc,islp), x_lobs_nuc
                     if (xhold < x_lobs_nuc) then 
                        accept = .false.
                        x = x-x1
                        k = k + 1
                        if (k > 20) then 
                           x = x + lobs_min
                        endif
                     endif
                     if (li .eq. 3) then 
                        if (abs(ystart) < atom_ymax + 15.0d0) then 
                           xtemp = xstart + x
                           if (xtemp >= atom_xmin .and. xtemp <=
     $                          atom_xmax)then
                              accept = .false.
                           endif
                        endif
                     endif
                  end do
                  if (x <= dlmax) then 
                     i = i + 1
                     nobs(islp) = nobs(islp) + 1
                     sobs(i,islp) = x
                     tau_obs(i,islp) = tobs
                     print *, islp,i,x, tau_obs(i,islp)
                  else
                     exit
                  endif
               end do
            endif
            print *, 'Total obstatcles on ', islp, ' = ', nobs(islp)
         endif
      end do
      
      end subroutine gen_obstacles

      subroutine xmin1(x,xlist,n,xhold,inuc)
      implicit none
      double precision :: xlist(n), x, xhold, d
      integer :: inuc, n, i
      xhold = abs(x-xlist(1))
      inuc = 1
      do i = 2, n
         d = abs(x-xlist(i))
         if (d < xhold) then 
            xhold = d
            inuc = i
         endif
      end do
 
      return
      end subroutine xmin1


      subroutine rgauss(x, len, xmean, sd)
!     !$  --------------------------------------------------------------
!     -------------
!     !$  Generate a random array of variables x with a Gaussian
!     distribution
!     !$  with a mean value xmean and a standard deviation sd. The
!     routine uses a trick 
!     !$  found in 
!     !$  Roes, P.B.M Van Oorschot, H.J.L : Kansrekening and statistek.
!     DUM 1978
!     !$  The routine makes use of ranmar to generate random numbers and
!     assumes
!     !$  ramrin has been called before to initialize ranmar
!     !$  --------------------------------------------------------------
!     -------------
      implicit double precision (a-h, o-z)
      double precision :: pi 
      integer :: len
      double precision :: rand(2), xmean, sd, x(len)
      double precision :: Pi2
      pi = 3.1415926535898
      Pi2 = 2.d0*pi
      n = len/2
      if ((2*n) .ne. len) n = n+1
      ix = 0
      do i = 1, n
           call ranmar(rand,2)
           v1 = rand(1)
           v2 = rand(2)
           scale = sd*sqrt(-2.0d0*log(v1))
           ix = ix + 1
           x(ix) = xmean + scale*cos(Pi2*v2)
           ix = ix + 1
         if (ix > len) then 
            return
         endif
         x(ix) = xmean + scale*sin(Pi2*v2)
      end do
      end subroutine rgauss

c$$$      function xmin1(x,xlist,n)
c$$$      implicit double precision (a-h, o-z)
c$$$      double precision :: xlist(n), x
c$$$!     !$    print *, 'xmin1'
c$$$!     !$    do i = 1,n
c$$$!     !$       print *, i, xlist(i)
c$$$!     !$    end do
c$$$      xmin1 = abs(x-xlist(1))
c$$$      do i = 2, n
c$$$         d = abs(x-xlist(i))
c$$$         if (d < dmin) then 
c$$$            xmin1 = d
c$$$         end if
c$$$      end do
c$$$      return
c$$$      end function xmin1



! This random number generator originally appeared in "Toward a Universa
! Random Number Generator" by George Marsaglia and Arif Zaman.          
! Florida State University Report: FSU-SCRI-87-50 (1987)                
!                                                                       
! It was later modified by F. James and published in "A Review of Pseudo
! random Number Generators"                                             
!                                                                       
! THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.             
!       (However, a newly discovered technique can yield                
!         a period of 10^600. But that is still in the development stage
!                                                                       
! It passes ALL of the tests for random number generators and has a peri
!   of 2^144, is completely portable (gives bit identical results on all
!   machines with at least 24-bit mantissas in the floating point       
!   representation).                                                    
!                                                                       
! The algorithm is a combination of a Fibonacci sequence (with lags of 9
!   and 33, and operation "subtraction plus one, modulo one") and an    
!   "arithmetic sequence" (using subtraction).                          
!                                                                       
! On a Vax 11/780, this random number generator can produce a number in 
!    13 microseconds.                                                   
!=======================================================================
                                                                        
      subroutine RMARIN(IJ,KL) 
! This is the initialization routine for the random number generator RAN
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328 
!                                                      0 <= KL <= 30081 
! The random number sequences created by these two seeds are of sufficie
! length to complete an entire calculation with. For example, if sveral 
! different groups are working on different parts of the same calculatio
! each group could be assigned its own IJ seed. This would leave each gr
! with 30000 choices for the second seed. That is to say, this random   
! number generator can create 900 million different subsequences -- with
! each subsequence having a length of approximately 10^30.              
                                                                        
                                                                        
      double precision :: U(97), C, CD, CM 
      integer I97, J97 
      logical TEST 
      common /raset1/ U, C, CD, CM, I97, J97, TEST 
                                                                        
      TEST = .FALSE. 
                                                                        
      if( IJ .lt. 0  .or.  IJ .gt. 31328  .or.                          &
     &    KL .lt. 0  .or.  KL .gt. 30081 ) then                         
          print '(A)', ' The first random number seed must have a value &
     &between 0 and 31328'                                              
          print '(A)',' The second seed must have a value between 0 and &
     &30081'                                                            
            stop 
      endif 
                                                                        
      i = mod(IJ/177, 177) + 2 
      j = mod(IJ    , 177) + 2 
      k = mod(KL/169, 178) + 1 
      l = mod(KL,     169) 
                                                                        
      do 2 ii = 1, 97 
         s = 0.0 
         t = 0.5 
         do 3 jj = 1, 24 
            m = mod(mod(i*j, 179)*k, 179) 
            i = j 
            j = k 
            k = m 
            l = mod(53*l+1, 169) 
            if (mod(l*m, 64) .ge. 32) then 
               s = s + t 
            endif 
            t = 0.5 * t 
    3    continue 
         U(ii) = s 
    2 continue 
                                                                        
      C = 362436.0 / 16777216.0 
      CD = 7654321.0 / 16777216.0 
      CM = 16777213.0 /16777216.0 
                                                                        
      I97 = 97 
      J97 = 33 
                                                                        
      TEST = .TRUE. 
      return 
      END subroutine rmarin
                                                                        
                                                                        
      subroutine ranmar(RVEC, LEN) 
! This is the random number generator proposed by George Marsaglia in   
! Florida State University Report: FSU-SCRI-87-50                       
! It was slightly modified by F. James to produce an array of pseudorand
! numbers.                                                              
      double precision ::  RVEC(*) 
      double precision :: U(97), C, CD, CM 
      integer I97, J97 
      logical TEST 
      common /raset1/ U, C, CD, CM, I97, J97, TEST 
                                                                        
      integer ivec 
                                                                        
      if( .NOT. TEST ) then 
         print '(A)',' Call the init routine (RMARIN) before calling RAN&
     &MAR'                                                              
         stop 
      endif 
                                                                        
      do 100 ivec = 1, LEN 
         uni = U(I97) - U(J97) 
         if( uni .lt. 0.0 ) uni = uni + 1.0 
         U(I97) = uni 
         I97 = I97 - 1 
         if(I97 .eq. 0) I97 = 97 
         J97 = J97 - 1 
         if(J97 .eq. 0) J97 = 97 
         C = C - CD 
         if( C .lt. 0.0 ) C = C + CM 
         uni = uni - C 
         if( uni .lt. 0.0 ) uni = uni + 1.0 
         RVEC(ivec) = uni 
  100 continue 
      return 
      END subroutine ranmar                          

      subroutine assign_disloc_only(dd1)
      implicit none
      type(dd),dimension(tot_disl):: dd1

!-----------------------------------------------------
! Local variables
!------------------------
      integer :: islp, i,j, k,l,Ntot, Ntotmax,iib
      double precision :: x,y,s,si
      double precision :: maxx, maxy, minx,miny, pii
!-------------------------------------------------------
      maxx = 0.d0; maxy = 0.d0; minx = 0.d0; miny = 0.d0
      range_ = 2.d0*process_xmax + 1000.0d0
      pii = 2.d0*asin(1.d0)
      Ntotmax = mxnslp*(mxndis+mxnnuc)
      print *, 'Size of dd1 entering assign is ',size(dd1), tot_disl
      k = 1;
      print *,'Dislocations are '
! ----- Read in the dislocations-------------------------------------
      do islp = 1,nslp
         if (ndis(islp) > 0) then 
            do i = 1,ndis(islp)
!               if (elem_dis(i,islp) .ne. 0) then 
               iib = locphi(islp)
               si = sdis(i,islp)
c$$$               x = xslp (islp) + si * cosphi (iib) 
c$$$               y = yslp (islp) + si * sinphi (iib)  
               x = rdis(1,i,islp)
               y = rdis(2,i,islp)
               dd1(k)%xy=dcmplx(x,y)
               dd1(k)%slipangle = phislp(iib)
               dd1(k)%burgers_dd(1) = b_dd(i,islp)
               dd1(k)%burgers_dd(2) = b_dd(i,islp)*cosphi(iib)
               dd1(k)%burgers_dd(3) = b_dd(i,islp)*sinphi(iib) 
               dd1(k)%slipplane_no = islp
               dd1(k)%stype = 1
               dd1(k)%velocity = vprev(i,islp)*
     $              cmplx(cosphi(iib),sinphi(iib))
               print *, i,islp, sdis(i,islp), iobpin(i,islp)
               k = k + 1
!               if (k >= tot_disl ) print *,k
!               end do
            end do
         endif
      end do
      print *, 'End Dislocations'
      print *, '-------------------------------------------------------'
      !tot_disloc = k-1
      end subroutine assign_disloc_only
    
      subroutine assign_source_only(dd1)
      implicit none
      type(dd),dimension(tot_size) :: dd1
!----------------------------------------------------
! Local variables
!-----------------------------------------------------
      integer :: islp, i,j, k,l,Ntot, Ntotmax,iib
      double precision :: x,y,s,si
!----------------------------------------------------
      k = tot_disl + 1
      do islp = 1,nslp
         if (nnuc(islp) > 0) then 
            do i = 1,nnuc(islp)
               iib = locphi(islp)
               si = snuc(i,islp)
c$$$               x = xslp (islp) + si * cosphi (iib) 
c$$$               y = yslp (islp) + si * sinphi (iib)  
               x = rnuc(1,i,islp)
               y = rnuc(2,i,islp)
c$$$               if (yendslp(islp) > 0.0d0) then 
c$$$                  x = xslp(islp) + si*cosphi(iib)
c$$$                  y = yslp(islp) + si*sinphi(iib)
c$$$               else
c$$$                  if (iib == 1) then 
c$$$                     x = xslp(islp) - si*cosphi(iib)
c$$$                  else
c$$$                     x = xslp(islp) + abs(si*cosphi(iib))
c$$$                  end if
c$$$                  y = -(yslp(islp) + si*sinphi(iib))
c$$$               end if
               if (y == 0 ) y = 1.e-9 
               if (x == 0) x = 1.e-9
               dd1(k)%xy=cmplx(x,y)
               dd1(k)%slipangle = phislp(iib)
               dd1(k)%burgers_dd(1) = blen
               dd1(k)%burgers_dd(2) = blen*cosphi(iib)
               dd1(k)%burgers_dd(3) = blen*sinphi(iib) 
               dd1(k)%slipplane_no = islp
               dd1(k)%stype = 2
               dd1(k)%velocity = 0.d0
               k = k + 1					
            end do
         endif
      end do
      print *, 'Total size =', tot_size, k-1
      !tot_size = k -1;
      end subroutine assign_source_only

      subroutine calcForce(dd1, start, finish)
      implicit none
      double precision :: xe, xnu, xlambda, xmu
      double precision :: cc(6,6)
      integer :: i_elas
      common /elastic/ xe, xnu, xlambda, xmu, cc, i_elas
      type(dd), dimension(:), intent(inout) :: dd1
      double complex :: b,  bc,t,ct, t1, ct1, tmp1,im,tc,zph, zmh, czmh,
     $     czph, zstar, z2, z3, zx
      double complex :: phi1, phi2, phi11, z, z0, zc, z0c, phi2c,zi,z0i,
     $     zic, z0ic, th, tmp2, dzds
      double precision :: sig11, sig12, sig22, bi, sin2i, cos2i, fg, fc
     $     ,fgg, s11, s12, s22
      double complex :: phi1z, phi2z, phi11z, tmp1z, tmp2z, gz
      double precision :: sig11z, sig12z, sig22z
      integer :: i,j,k, start, start1,finish, finish1
      double precision :: e, mu, nu, factor, h, rcore,x, y, x0,pii
      double precision, dimension(2) :: zstar1(2),zstar2(2), z1(2)
      
      rcore = 2.d0*blen
      pii = 2.d0*asin(1.d0)
c$$$      start = 1
c$$$      finish = size(dd1)
      start1 = 1
      finish1 = tot_disl
      im = cmplx(0.,1.)
      mu = xmu
      nu = xnu
      factor = mu/(4.*pii*(1-nu))
!     print *, 'Factor = ', factor
      do i = start,finish
         dd1(i)%forced = 0.d0
         dd1(i)%forcedg = 0.d0
         dd1(i)%phi1 = 0.d0; dd1(i)%phi2 = 0.d0; dd1(i)%phi11 = 0.d0;
         dd1(i)%dphi1 = 0.d0; dd1(i)%dphi2 = 0.d0; dd1(i)%dphi11 = 0.d0;
         phi1 = 0.d0; phi2 = 0.d0; phi11 = 0.d0; sig11 = 0.0d0; sig12 =
     $        0.0d0; sig22 = 0.d0;
         phi1z = 0.d0; phi2z = 0.d0; phi11z = 0.d0; sig11z = 0.d0;
     $        sig12z = 0.d0; sig22z = 0.d0;
         
         z = dd1(i)%xy
         zc = conjg(z)
!     write (*,fmt='(I2,1x,2E12.5,2x,2e12.5)')i, z,
!     cmplx(dd1(i)%burgers_dd(2),dd1(i)%burgers_dd(3))
         sig11 = 0.; sig12 = 0.; sig22 = 0.
         do j = start1, finish1
            tmp1z = 0.d0; tmp2z = 0.d0;
            z0 = dd1(j)%xy
            z0c = conjg(z0)
            b = cmplx(dd1(j)%burgers_dd(2),dd1(j)%burgers_dd(3))
            bc = conjg(b)
  
            if (abs(z-z0) > 2*blen) then
               phi1 = -im*factor*b/(z-z0);
               phi11 =  im*factor*b/((z-z0)*(z-z0));
               phi2 = im*factor*bc/(z-z0)
               tmp1 = 2.d0*(phi1 + conjg(phi1));
               tmp2 = -2.d0*((z-z0)*conjg(phi11) + conjg(phi2))
               
               
               phi1z = im*factor*b/((z-z0)*(z-z0));
               phi11z =  -2.d0*im*factor*b/((z-z0)*(z-z0)*(z-z0));
               phi2z = -im*factor*bc/((z-z0)*(z-z0))
               tmp1z = 2.d0*(phi1z + conjg(phi1z));
               tmp2z = -2.d0*((z-z0)*conjg(phi11z) + conjg(phi2z))
               
               sig11z = sig11z + real(0.5d0*(tmp1z + tmp2z))
               sig22z = sig22z + real(0.5d0*(tmp1z - tmp2z))
               sig12z = sig12z + aimag(0.5d0*(tmp2z))
               sig11 = sig11 + real(0.5d0*(tmp1 + tmp2))
               sig22 = sig22 + real(0.5d0*(tmp1 - tmp2))
               sig12 = sig12 + aimag(0.5d0*(tmp2))

            else
               t = z-z0
               if (i .ne. j) then 
                  if (real(t) == 0) then
                     th = dcmplx(0,asin(sin(dd1(j)%slipangle)))
                  else
                     th = dcmplx(0,aimag(log(t)))
                  endif 
                  zstar = rcore*exp(th)
                  phi1 = -im*factor*b/zstar
                  phi11 = im*factor*b/(zstar*zstar)
                  phi2 = im*factor*bc/(zstar)
                  tmp1 = 2.d0*(phi1 + conjg(phi1));
                  tmp2 = -2.d0*((z-z0)*conjg(phi11) + conjg(phi2))
                  sig11 = sig11 + real(0.5d0*(tmp1 + tmp2))
                  sig22 = sig22 + real(0.5d0*(tmp1 - tmp2))
                  sig12 = sig12 + aimag(0.5d0*(tmp2))
               endif
            endif
         enddo
         sin2i = sin(2*dd1(i)%slipangle)
         cos2i = cos(2*dd1(i)%slipangle)
         bi = dd1(i)%burgers_dd(1)
         fg = dd1(i)%burgers_dd(1)*((sig22-sig11)*0.5*sin2i +
     $        sig12*cos2i)
         dd1(i)%forced = dd1(i)%forced + fg
         fgg = dd1(i)%burgers_dd(1)*((sig22z-sig11z)*0.5*sin2i + sig12z
     $        *cos2i)
         dd1(i)%forcedg= dd1(i)%forcedg + fgg
         !print *, 'xxxxx', i, dd1(i)%burgers_dd,dd1(i)%forced
      enddo			
      end subroutine calcForce


		
      subroutine dislp (r,u)
!-----------------------------------------------------------------------
!     Determine displacement vector u at node n caused by dislocations
!     
!     bout(1,islp) = total length of Burgers vectors of dislocations having
!     moved out of the bottom side (1) of the cell
!-----------------------------------------------------------------------
!     Edition log:
!     06/10/94 The field of escaped dislocations is corrected to b/4
!-----------------------------------------------------------------------
      implicit none
      double precision u(3), r(3)
      double precision dxx, dyy, r2, rcore, dxi, dyi, cosi, sini, r1
      double precision dx2, dy2, a, pi, t, u1, u2, dx0, dy0, x,y,dxy2
      double precision xnu, xmu, cc(6,6), xlambda, xe
      integer i_elas
      common /elastic/ xe, xnu, xlambda, xmu, cc, i_elas
      integer i, j, islp, li
      double precision sign1(2)
      x = r(1)
      y = r(2)
      pi = 2.d0*asin(1.0d0)
      print *, 'dislp',nslp,xnu

!     
!     
!.....For each node, loop over all active dislocations inside the cell
      u(1)=0.0d0
      u(2)=0.0d0
      do 2 islp=1,nslp
         li=locphi(islp)
         sini=sinphi(li)
         cosi=cosphi(li)
         dx0=x-xslp(islp)
         dy0=y-yslp(islp)
         dxi= dx0*cosi+dy0*sini
         dyi=-dx0*sini+dy0*cosi
         u1=0.0d0
         u2=0.0d0
         do 4 i=1,ndis(islp)
            rcore=2.0d0*dabs(b_dd(i,islp))
            dxx=dxi-sdis(i,islp)
            dyy=dyi
            dx2=dxx*dxx
            dy2=dyy*dyy
            r2=dx2+dy2
!     Cut off when closer than core radius
            if(r2.lt.rcore**2) then
               r1=sqrt(r2)
               dxx=dxx/r1*rcore
               dyy=dyy/r1*rcore
               dx2=dxx*dxx
               dy2=dyy*dyy
               r2=rcore**2
            endif
            dxy2=dxx**2+dyy**2
            a=b_dd(i,islp)/2.0d0/pi/(1.0d0-xnu)
            if(dyy.eq.0.0d0) then
               t=dsign(pi/2.0d0,dxx)
            else
               t=datan(dxx/dyy)
            endif
!.....first the local displacements
            u1=u1+a*(dxx*dyy/2.0d0/dxy2-(1.0d0-xnu)*t)
            u2=u2+a*(dyy**2/2.0d0/dxy2-(1.0d0-2.0d0*xnu)/4.0d0*
     $           log(dxy2/b_dd(i,islp)**2))
 4    continue
!.....transform to global displacements
      u(1)=u(1)+cosi*u1-sini*u2
      u(2)=u(2)+sini*u1+cosi*u2
      u(3) = 0.0d0
 2    continue

!.....Dislocations that have escaped from the cell previously:
!.....The displacement field generated by a + dislocation moved out at the
!.....right-hand side (2) is
!.....          b/4 for y>yd
!.....   u(1)=
!.....         -b/4 for y<yd
!.....Dislocations having moved out at the bottom side (1) generate
!.....displacements with the opposite sign.
      sign1(1)=-1.0d0
      sign1(2)= 1.0d0
      do 3 islp=1,nslp
      li=locphi(islp)
      sini=sinphi(li)
      cosi=cosphi(li)
      dx0=x-xslp(islp)
      dy0=y-yslp(islp)
      dxi= dx0*cosi+dy0*sini
      dyi=-dx0*sini+dy0*cosi
      u1=0.0d0
      u2=0.0d0
 !     sout(1,islp) = -1.d10
 !     sout(2,islp) = 1.d10
 !      soutq = 1.d10
      do 11 i=1,2
        if(bout(i,islp).eq.0) goto 11
        !soutq = sign1(i)*soutq
        rcore=2.d0*blen
        dxx=dxi-sdis_out(i,islp)
        dyy=dyi
        dx2=dxx*dxx
        dy2=dyy*dyy
        r2=dx2+dy2
!       Cut off when closer than core radius
        if(r2.lt.rcore**2) then
          r1=dsqrt(r2)
          dxx=dxx/r1*rcore
          dyy=dyy/r1*rcore
          dx2=dxx*dxx
          dy2=dyy*dyy
          r2=rcore**2
        endif
        dxy2=dxx**2+dyy**2
        a=bout(i,islp)/2.0d0/pi/(1.0d0-xnu)
        if(dyy.eq.0.0d0) then
          t=dsign(pi/2.0d0,dxx)
        else
          t=datan(dxx/dyy)
        endif
!.....first the local displacements
        if(dyy.gt.0.0d0) then
          u1=u1+sign1(i)*bout(i,islp)/4.0d0
        else
          u1=u1-sign1(i)*bout(i,islp)/4.0d0
        endif
        u2=u2+a*(dyy**2/2.0d0/dxy2-(1.0d0-2.0d0*xnu)/4.0d0*dlog(dxy2
     $       /bout(i,islp)**2))
 11     continue
!.....Now transform to global displacements
      u(1)=u(1)+cosi*u1-sini*u2
      u(2)=u(2)+sini*u1+cosi*u2

 3    continue

      return
      end subroutine dislp


      subroutine rlsdis(tau)
!-----------------------------------------------------------------------
! Release pinned dislocations when the resolved shear stress `tau'
! exceeds a critical value
! 08/04/01 dislocations pinned at junctions skipped at present
!-----------------------------------------------------------------------
      double precision :: tau(mxndis, nslp), tauii
      integer :: islp, i,j

      do islp=1,nslp
         do i=1,ndis(islp)
            if (iobpin(i,islp) > 0) then
               tauii=dabs(tau(i,islp))
               if (tauii.gt.tau_obs(iobpin(i,islp),islp)) then
                  write(6,625) i,islp,iobpin(i,islp)
 625              format(' Disloc. ',i3,' on slip plane ',i4
     $                 ,' released from obstacle ',i2)
                  iobpin(i,islp)=0
               endif
            endif
         end do
      end do
      return
      end subroutine rlsdis

      subroutine undrex(v,c)
!-----------------------------------------------------------------------
! Determine underrelaxation factor for updating dislocation positions.
! Currently, underrelaxation with a factor 0.75 is used when the
! velocity of two like-signed dislocations (in pile-ups) changes sign
! as compared to the previous increment.
! In addition, underrelaxation with a factor 0.5 is used for ANY
! sign change in velocity
!-----------------------------------------------------------------------
! Edition log
! 02/09/96 Underrelaxation implemented also for any change of sig
!-----------------------------------------------------------------------
      implicit none
      double precision :: v(mxndis,nslp), c(mxndis,nslp)
      integer :: i,j,k,l, islp
      double precision :: bi
      do  islp=1,nslp
         do i=1,ndis(islp)
            c(i,islp)=1.0d0
            if((v(i,islp)*vprev(i,islp))< 0.0d0) then 
               c(i,islp)=0.5d0*c(i,islp)
               bi=b_dd(i,islp)
               j=i-1
               if(j.ge.1) then
                  if((bi*b_dd(j,islp)).gt.0.0d0) then 
                     c(i,islp)=0.75d0*c(i,islp)
                  endif
               endif
               j=i+1
               if(j.le.ndis(islp)) then
                  if((bi*b_dd(j,islp)).gt.0.0d0)  then 
                     c(i,islp)=0.75d0*c(i,islp)
                  endif
               endif
            endif
         end do
      end do
      return
      end subroutine undrex


      end module mod_dd_slip
      function tlce(xx)
      implicit double precision(a-h,o-z)
      tlce=dint(xx*1.0d6)/1.0d6
      return
      end function tlce
