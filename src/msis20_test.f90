module test20

  contains
!##############################################################################
! MSIS� (NRL-SOF-014-1) SOFTWARE
!
! MSIS� is a registered trademark of the Government of the United States of 
! America, as represented by the Secretary of the Navy. Unauthorized use of 
! the trademark is prohibited. 
!
! The MSIS� Software (hereinafter Software) is property of the United States 
! Government, as represented by the Secretary of the Navy. Methods performed
! by this software are covered by U.S. Patent Number 10,641,925. The Government
! of the United States of America, as represented by the Secretary of the Navy, 
! herein grants a non-exclusive, non-transferable license to the Software for 
! academic, non-commercial, purposes only. A user of the Software shall not: 
! (i) use the Software for any non-academic, commercial purposes, (ii) make 
! any modification or improvement to the Software, (iii) disseminate the 
! Software or any supporting data to any other person or entity who will use 
! the Software for any non-academic, commercial purposes, or (iv) copy the 
! Software or any documentation related thereto except for (a) distribution 
! among the user�s personal computer systems, archival, or emergency repair 
! purposes, or (b) distribution for non-commercial, academic purposes, without 
! first obtaining the written consent of IP Counsel for the Naval Research 
! Laboratory. 
!
! As the owner of MSIS�, the United States, the United States Department of 
! Defense, and their employees: (1) Disclaim any warranties, express, or 
! implied, including but not limited to any implied warranties of 
! merchantability, fitness for a particular purpose, title or non-infringement, 
! (2) Do not assume any legal liability or responsibility for the accuracy, 
! completeness, or usefulness of the software, (3) Do not represent that use of 
! the software would not infringe privately owned rights, (4) Do not warrant 
! that the software will function uninterrupted, that is error-free or that any 
! errors will be corrected.
!
! BY USING THIS SOFTWARE YOU ARE AGREEING TO THE ABOVE TERMS AND CONDITIONS.  
!##############################################################################

!!! ===========================================================================
!!! NRLMSIS 2.0:
!!! Neutral atmosphere empirical model from the surface to lower exosphere
!!! John Emmert (john.emmert@nrl.navy.mil)
!!! Doug Drob (douglas.drob@nrl.navy.mil)
!!! ===========================================================================

!==================================================================================================
! MSISTEST: Test program for NRLMSIS 2.0
!==================================================================================================
subroutine msistest()

  use msis_init, only          : msisinit
  use msis_gtd8d

  implicit none

  integer, parameter          :: nrec = 200

  integer                     :: iyd, mass
  real(4)                     :: sec, alt, glat, glong, stl, f107a, f107, ap(7), apd
  real(4)                     :: d(9),t(2)
  
  integer                     :: i
  character(128)              :: dummy

  !Initialize model
  call msisinit(parmpath='data', parmfile='/msis20.parm')

  !Open input and output files, loop through records, and call model
  open(77,file='msis20_test_in.txt',status='old')
  open(78,file='msis20_test_out.txt',status='replace')
  read(77,*) dummy
  write(78,'(9a7,9a13,a8)') &
    'iyd','sec','alt','glat','glong','stl','f107a','f107','Ap','He','O','N2','O2','Ar','rho','H','N','O*','T'
  do i = 1,200
	  read(77,*) iyd,sec,alt,glat,glong,stl,f107a,f107,apd
    ap(1) = apd
    call gtd8d(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)
    write(78,'(2i7,3f7.1,f7.2,3f7.1,9e13.4,f8.2)')  &
      iyd,int(sec),alt,glat,glong,stl,f107a,f107,ap(1),d(1:9),t(2)
  enddo
  close(77)
  close(78)

  stop

end subroutine msistest

end module