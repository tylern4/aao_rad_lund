!      include file BCS
!
!_begin_doc
!  RCS ID string
!  $Id: bcs.inc,v 1.1 2000/05/18 23:07:48 lcsmith Exp $
!
!  Documentation for include file BCS
!
!  Purpose: Standard common for BOS banks
!  --------
!
!  Used in routines: all
!  ----------------
!
!  Notes:
!  ------
!
!  Author:   Dieter Cords      Created:  Mon Jun 12 17:50:51 EDT 1995
!  -------
!
!  Major revisions:
!  JAM - changed the parameter to 700,000
!  ----------------
!
!_end_doc
!
!
!  Remember NO NESTED INCLUDES
!
!  COMMON BLOCK NAME should be the same as include file name!
!
!_begin_var
!  variables for include file BCS:
!  -----------------------------------
!  RCS information: 
CHARACTER*(132)  CIBCS
DATA CIBCS/&
        '$Id: bcs.inc,v 1.1 2000/05/18 23:07:48 lcsmith Exp $'&
        /
!----6----------------------------------------------------------------72
!  
!  Variables contained in BCS common block
!  --------------------------------------------
!
!_end_var
!
INTEGER Nbcs
PARAMETER (Nbcs = 700000)
!
INTEGER IW(Nbcs), JUNK(5)
REAL RW(Nbcs)
INTEGER*2 IW16(2 * Nbcs)
BYTE IW08(4 * Nbcs)
EQUIVALENCE (IW(1), RW(1), IW08(1), IW16(1))
! 
! IWCOL and IWROW provide a quick way to find the number of columns
! and rows for a particular bank
!
!
COMMON/BCS/JUNK, IW
INTEGER IWROW(Nbcs), IWCOL(Nbcs), IWNAME(Nbcs), IWNR(Nbcs), &
        IWNXT(Nbcs)
EQUIVALENCE (IW(1), IWCOL(6))
EQUIVALENCE (IW(1), IWROW(5))
EQUIVALENCE (IW(1), IWNAME(4))
EQUIVALENCE (IW(1), IWNR(3))
EQUIVALENCE (IW(1), IWNXT(2))
!
!
save /BCS/
!----6----------------------------------------------------------------72



