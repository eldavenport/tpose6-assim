#include "COST_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

      subroutine cost_final_restore ( mythid, lastdiva )

c     ==================================================================
c     SUBROUTINE cost_final_restore
c     ==================================================================
c
c     Restore cost function value
c     This is required in conjunction with DIVA
c     heimbach@mit.edu 5-Mar-2003
c
c     ==================================================================
c     SUBROUTINE cost_final_restore
c     ==================================================================

      implicit none

c     == global variables ==

#include "EEPARAMS.h"
#include "SIZE.h"

#include "cost.h"

c     == routine arguments ==

      integer mythid
      logical lastdiva

CMM # ifdef ALLOW_DIVIDED_ADJOINT

c     == local variables ==
      logical exst

      integer idivbeg, idivend
c     == end of interface ==

c--   Each process has calculated the global part for itself.
cph      IF (myProcId .eq. 0) THEN
         _BEGIN_MASTER( mythid )
c
           inquire(file='divided.ctrl',exist=exst)
           if (exst) then
              open(unit=76,file='divided.ctrl',form='formatted')
              read(unit=76,fmt=*) idivbeg,idivend
              close(unit=76)
           else
              idivbeg = -2
           endif
c
CMM           if ( idivbeg .EQ. 0 ) then
              lastdiva = .TRUE.
              open(unit=76,file='costfinal')
              read(76,*) fc
              close(76)
CMM           else
CMM              lastdiva = .FALSE.
CMM           endif
c
         _END_MASTER( mythid )
cph      ENDIF
      _BARRIER

CMM #endif

      end
