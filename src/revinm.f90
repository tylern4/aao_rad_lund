SUBROUTINE REVINM(EVIVAR, FILNAM, FULLPATH)
    !
    !_begin_doc
    !  RCS ID string
    !  $Id: revinm.F,v 1.1.1.1 1997/03/28 20:58:33 marki Exp $
    !
    !  Documentation for subroutine REVINM
    !
    !  Purpose:  Given the enviroment variable EVIVAR and the filename FILNAM
    !  --------  returns the full path FULLPATH to the file.
    !
    !  Calling Sequence:
    !  ----------------
    !
    !  Input Parameters:  (EVIVAR - C - Name of enviroment variable)
    !  ----------------    FILNAM - C - Name of file
    !
    !  Output Parameters:  (FULLPATH - C - EVIVAR/FILNAM )
    !  -----------------
    !
    !  Called from:
    !  ------------
    !
    !  Other routines:
    !  ---------------
    !
    !  Notes:
    !  ------
    !
    !  Author:   Arne Freyberger      Created:  Wed Mar 15 13:39:27 EST 1995
    !  -------
    !
    !  Major revisions:
    !  ----------------
    !  4/4/95   Removed messaging so that code can be used with the
    !  apf      bulk of the RECUTL/RECSIS package.
    !_end_doc
    !
    IMPLICIT NONE
    SAVE
    !
    !_begin_inc
    !  include files :
    !  ---------------------
    ! BOS common block  uncomment the next line for BOS include file
    !#include "bcs.inc"
    !_end_inc
    !
    !_begin_var
    !  input/output variables:
    !  -----------------------
    !
    !  Local pre-defined variables:
    !  ---------------------------
    !  RCS information:
    CHARACTER*132  CFILE, CREVIS, CSTATE, CDATE, CAUTHO, CRCSID
    PARAMETER (CFILE = '$RCSfile: revinm.F,v $')
    PARAMETER (CREVIS = '$Revision: 1.1.1.1 $')
    PARAMETER (CSTATE = '$State: Exp $')
    PARAMETER (CDATE = '$Date: 1997/03/28 20:58:33 $')
    PARAMETER (CAUTHO = '$Author: marki $')
    DATA CRCSID/&
            '$Id: revinm.F,v 1.1.1.1 1997/03/28 20:58:33 marki Exp $'&
            /
    !  Module information:
    CHARACTER*(*)  CRNAME, CRAUTH
    CHARACTER*100  CRMESS
    PARAMETER (CRNAME = 'REVINM')
    PARAMETER (CRAUTH = 'Arne Freyberger')
    !
    !  Local User defined variables:
    !  -----------------------------
    INTEGER LENOCC
    EXTERNAL LENOCC
    CHARACTER*(*) EVIVAR, FILNAM, FULLPATH
    INTEGER ICEND, ICBEG, IADD, ICEND2
    !_end_var
    !
    !  executable code for routine REVINM:
    !  -------------------------------------
    !
    !
    CALL GETENV(EVIVAR, FULLPATH)
    ICEND = LENOCC(FULLPATH)
    IF (FULLPATH(1:1) .EQ. ' ') THEN
        FULLPATH = FILNAM
        ICEND2 = LENOCC(EVIVAR)
        WRITE(CRMESS, 10)EVIVAR(1:ICEND2)
        !        CALL RECLOG(CRNAME,'W',CRMESS)
        RETURN
    ENDIF
    ICEND = ICEND + 1
    !
    FULLPATH(ICEND:ICEND) = '/'
    !
    IADD = LENOCC(FILNAM)
    ICBEG = ICEND + 1
    ICEND = ICBEG + IADD
    FULLPATH(ICBEG:ICEND) = FILNAM(1:IADD)
    !
    10    FORMAT('Enviroment variable:  ', A, '  not defined')
    RETURN
END
!
!------------------------------------------------------------------------------