      SUBROUTINE REVINM(EVIVAR,FILNAM,FULLPATH)
c
c_begin_doc
c  RCS ID string
c  $Id: revinm.F,v 1.1.1.1 1997/03/28 20:58:33 marki Exp $
c
c  Documentation for subroutine REVINM
c
c  Purpose:  Given the enviroment variable EVIVAR and the filename FILNAM
c  --------  returns the full path FULLPATH to the file.
c
c  Calling Sequence:
c  ----------------
c
c  Input Parameters:  (EVIVAR - C - Name of enviroment variable)
c  ----------------    FILNAM - C - Name of file
c
c  Output Parameters:  (FULLPATH - C - EVIVAR/FILNAM )
c  -----------------
c
c  Called from:
c  ------------
c
c  Other routines:
c  ---------------
c
c  Notes:
c  ------
c
c  Author:   Arne Freyberger      Created:  Wed Mar 15 13:39:27 EST 1995
c  -------
c
c  Major revisions:
c  ----------------
c  4/4/95   Removed messaging so that code can be used with the
c  apf      bulk of the RECUTL/RECSIS package.
c_end_doc
c
      IMPLICIT NONE
      SAVE
c
c_begin_inc
c  include files :
c  ---------------------
c BOS common block  uncomment the next line for BOS include file
c#include "bcs.inc"
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
c
c  Local pre-defined variables:
c  ---------------------------
c  RCS information:
      CHARACTER*132  CFILE, CREVIS, CSTATE, CDATE, CAUTHO, CRCSID
      PARAMETER (CFILE=  '$RCSfile: revinm.F,v $')
      PARAMETER (CREVIS= '$Revision: 1.1.1.1 $')
      PARAMETER (CSTATE= '$State: Exp $')
      PARAMETER (CDATE=  '$Date: 1997/03/28 20:58:33 $')
      PARAMETER (CAUTHO= '$Author: marki $')
      DATA CRCSID/
     1'$Id: revinm.F,v 1.1.1.1 1997/03/28 20:58:33 marki Exp $'
     3/
c  Module information:
      CHARACTER*(*)  CRNAME, CRAUTH
      CHARACTER*100  CRMESS
      PARAMETER (CRNAME='REVINM')
      PARAMETER (CRAUTH='Arne Freyberger')
c
c  Local User defined variables:
c  -----------------------------
      INTEGER LENOCC
      EXTERNAL LENOCC
      CHARACTER*(*) EVIVAR, FILNAM, FULLPATH
      INTEGER ICEND, ICBEG, IADD, ICEND2
c_end_var
c
c  executable code for routine REVINM:
c  -------------------------------------
c
c
      CALL GETENV(EVIVAR,FULLPATH)
      ICEND = LENOCC(FULLPATH)
      IF (FULLPATH(1:1) .EQ. ' ') THEN
        FULLPATH = FILNAM
        ICEND2 = LENOCC(EVIVAR)
        WRITE(CRMESS,10)EVIVAR(1:ICEND2)
c        CALL RECLOG(CRNAME,'W',CRMESS)
        RETURN
      ENDIF
      ICEND = ICEND + 1
c
      FULLPATH(ICEND:ICEND) = '/'
c
      IADD  = LENOCC(FILNAM)
      ICBEG = ICEND + 1
      ICEND = ICBEG + IADD
      FULLPATH(ICBEG:ICEND) = FILNAM(1:IADD)
c
10    FORMAT('Enviroment variable:  ',A,'  not defined')
      RETURN
      END
c
c------------------------------------------------------------------------------
