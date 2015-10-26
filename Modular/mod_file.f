!****************************************************************************
!**   
!**   MODULE mod_file : contains routines relevant to file handling
!**   
!**   
!**   Variable Definitions:
!**   ---------------------
!**   
!**   none
!**   
!**   
!**   Contains Routines:
!**   ------------------
!**   
!**   FileExists   - check for the existence of a file
!**   GetFreeUnit  - get a unit number not being currently used
!**   iofile       - opens files and returns a unit number (uses getfreeunit)
!**   
!****************************************************************************

      module mod_file

      contains

!------------------------------------------------------------------
! FileExists -- returns .true. if file 'filename' exists in the current
!               directory
!
!      Passed Parameters :
!          filename (in) : filename
!             print (in) : logical variable, set to .true. if function
!                          should print a message if file is not found
!
!      Module Parameters :
!                      none
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
      logical function FileExists(filename,print)

      implicit none

!** Transferred Variables **|
      character(len=*), intent(in) :: filename
      logical, intent(in) :: print

      inquire(file=filename,exist=FileExists)
      if (.not.FileExists.and.print) then
         print *,'***ERROR: File ',trim(filename),' not found.'
      end if
      return
      end function FileExists


!------------------------------------------------------------------
! GetFreeUnit -- return the number of a unit number currently not
!                being used.
!
!      Passed Parameters :
!                      none
!
!      Module Parameters :
!                      none
!
!      Algorithm :
!            Loop over all possible unit numbers until an unallocated
!            one is found, otherwise report an error and stop.
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
      integer function GetFreeUnit()

      implicit none

!** Local Variables **|
      logical InUse
!     
!     reserved file numbers: 5 = std in
!     6 = std out
!     99 = abort.dat
!     
      do GetFreeUnit=7,98
         inquire(unit=GetFreeUnit,opened=InUse)
         if (.not.InUse) return
      enddo
      print *,'***ERROR: Could not obtain a free unit handle.'
      stop
      end function GetFreeUnit
!------------------------------------------------------------------
! iofile -- if filename is open, returns the unit number, otherwise
!           gets a free unit number and opens filename.
!
!      Passed Parameters :
!          filename (in) : filename
!            format (in) : character string for file format
!            logic (out) : unit number
!
!      Module Parameters :
!                      uses GetFreeUnit
!
!      Algorithm :
!            self explanatory
!
!      Notes :
!          none
!
!      Author :
!            R.E.Miller and Feap (09/25/95)
!
!      Revisions :
!              none
!
!--
      subroutine iofile(filename,format,logic,verbose)
      implicit none
!     
!---- io file management: open file and return unit number
!     
      character*80 filename
      character*11 format
      integer logic
      logical ofl,verbose
!     :
      inquire(file=filename,opened=ofl)
      if (.not.ofl) then
         logic=GetFreeUnit()
         open(logic,file=filename,status='unknown',form=format)
         if(verbose) then
            write(*,*) '  ** Opening ',format,' file:'
            write(*,*) '        ',filename
         endif
      else
         inquire(file=filename,number=logic)
         if(verbose) then
            write(*,*) '  ** Accessing previously opened file:'
            write(*,*) '        ',filename
         endif
      end if
      return
      end subroutine iofile

      End module mod_file


