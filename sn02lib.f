!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn02lib.f
!
!     snTitle   snInit   snSpec   snchkA   snJac
!     snMem     snMemA   snMemA0  snMemB
!     snLog     snLog2   snSTOP
!     snSet     snSeti   snSetr
!     snGet     snGetc   snGeti   snGetr
!     snRetH
!
!     09 Mar 2004: snSolF implemented.
!     17 Jun 2004: snSolF always flags infeasible jbInf1 as I.
!     13 Dec 2013: Moved snEXIT, snWRAP, snSolF to sn04wrap.f.
!     24 Oct 2014: Added snGetStats, which prints problem statistics.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snTitle( title )

      character
     &     title*30

!     ==================================================================
!     snTitle sets the title for snopt.
!     ==================================================================

      title  = 'S N O P T  7.5-1.1  (Jun 2015)'
!---------------123456789|123456789|123456789|--------------------------

      end ! subroutine snTitle

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snInit
     &   ( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iPrint, iSumm, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     snInit  is called by the user to do the following:
!     1. Open default files (Print, Summary).
!     2. Initialize title.
!     3. Set options to default values.
!
!     15 Nov 1991: First version.
!     14 Jul 1997: Thread-safe version.
!     02 Oct 1997: Character workspace added.
!     15 Oct 2003: snEXIT and snPRNT added.
!     18 Jun 2007: First 500 elements of cw, iw and rw are initialized
!     18 Jun 2008: Global call-status values added.
!     ==================================================================
      external
     &     s1outpt
      character
     &     Solver*6, str*80, str2*80, title*30
      integer
     &     idummy, inform, iSpecs, iStdo, lvlTim, maxcu, maxcw, maxiu,
     &     maxiw, maxru, maxrw, nnCon, nnJac, nnL, nnObjU, npStatus,
     &     qpStatus, s1outpt
!     ------------------------------------------------------------------
      parameter         (maxru     =   2) ! start of SNOPT part of rw
      parameter         (maxrw     =   3) ! end   of SNOPT part of rw
      parameter         (maxiu     =   4) ! start of SNOPT part of iw
      parameter         (maxiw     =   5) ! end   of SNOPT part of iw
      parameter         (maxcu     =   6) ! start of SNOPT part of cw
      parameter         (maxcw     =   7) ! end   of SNOPT part of cw
      parameter         (nnJac     =  21) ! # nonlinear Jac, variables
      parameter         (nnObjU    =  22) ! # user-defined nnObj vars
      parameter         (nnCon     =  23) ! # of nonlinear constraints
      parameter         (nnL       =  24) ! nonlinear vars
      parameter         (lvlTim    = 182) ! Timing level
      parameter         (qpStatus  = 235) ! QP user-routine call-status
      parameter         (npStatus  = 236) ! NP user-routine call-status
      parameter         (idummy    =-11111)
!     ------------------------------------------------------------------
      character          dashes*30
      data               dashes /'=============================='/
!     ------------------------------------------------------------------
      Solver = 'SNINIT'

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
         !--------------------------------------------------------------
         ! Not enough workspace to do ANYTHING!
         ! Print and exit without accessing the work arrays.
         !--------------------------------------------------------------
         inform = 81       ! Work arrays must have at least 500 elements
         call snWRAP( inform, Solver, str, str2, iw, leniw )
         go to 999
      end if

      !-----------------------------------------------------------------
      ! Initialize cw, iw, rw so that they may be copied safely.
      !
      ! This also sets the options to a specific "undefined" state.
      ! snopt  will check the options later and maybe print them.
      !-----------------------------------------------------------------
      call s3unsetAll
     &   ( cw, lencw, iw, leniw, rw, lenrw )

      cw(1)  = Solver//'  '

      !-----------------------------------------------------------------
      ! Initialize some default values.
      !-----------------------------------------------------------------
      iSpecs    = 0
      iStdo     = s1outpt( )
      iw( 10)   = iStdo   ! Standard Output
      iw( 11)   = iSpecs  ! Specs file (default)
      iw( 12)   = iPrint  ! Print file
      iw( 13)   = iSumm   ! Summary file

      iw(maxcu) = 500
      iw(maxiu) = 500
      iw(maxru) = 500
      iw(maxcw) = lencw
      iw(maxiw) = leniw
      iw(maxrw) = lenrw

      !-----------------------------------------------------------------
      ! These dimensions need to be initialized for an MPS run.
      !-----------------------------------------------------------------
      iw(nnCon ) = 0
      iw(nnJac ) = 0
      iw(nnObjU) = 0
      iw(nnL   ) = 0

      call snTitle( title )
      call s1init ( title, iw, leniw, rw, lenrw )

      call snPRNT (11, '         '//dashes, iw, leniw )
      call snPRNT ( 1, '         '//title , iw, leniw )
      call snPRNT ( 1, '         '//dashes, iw, leniw )

      call snPRNT (12, ' '//dashes, iw, leniw )
      call snPRNT ( 2, ' '//title , iw, leniw )
      call snPRNT ( 2, ' '//dashes, iw, leniw )

      !-----------------------------------------------------------------
      ! Initialize some global values.
      !-----------------------------------------------------------------
      iw(qpStatus) = idummy
      iw(npStatus) = idummy
      iw(lvlTim)   = 3

  999 return

      end ! subroutine snInit

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snSpec
     &   ( iSpecs, iExit, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iSpecs, iExit, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     snSpec  may be called by the user to read a Specs file.
!
!     07 Feb 1998: First version of snSpec.
!     01 Aug 2003: s3file now has a "title" parameter.  Use ' '.
!     13 Jul 2005: Included error count in value of iExit.
!     22 Apr 2007: Exit 107 added (some keywords unrecognized)
!     ==================================================================
      character
     &     Solver*6, str*80, str2*80
      integer
     &     Errors, iPrint, iSumm, Calls
      external
     &     s3opt
!     ------------------------------------------------------------------
      Solver = 'SNSPEC'

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
!        ---------------------------------------------------------------
!        Not enough workspace to do ANYTHING!
!        Print and exit without accessing the work arrays.
!        ---------------------------------------------------------------
         iExit = 81        ! Work arrays must have at least 500 elements
         call snWRAP( iExit, Solver, str, str2, iw, leniw )
         go to 999
      end if

      cw(1)  = Solver//'  '

      if (iSpecs .le. 0  .or.  iSpecs .gt. 99) then
         iExit = 131      ! iSPECS out of range
         go to 800
      end if

      iw( 11)   = iSpecs  ! Specs (options) file

      iPrint    = iw( 12) ! Print file
      iSumm     = iw( 13) ! Summary file

      Calls     = 1

!     ------------------------------------------------------------------
!     Read the Specs file.
!     snopt  will check the options later and maybe print them.
!     ------------------------------------------------------------------
      call s3file
     &   ( iExit, Calls, iSpecs, s3opt, ' ', iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

  800 if (iExit .eq. 0) then
         if (Errors .eq. 0) then
            iExit = 101         ! SPECS file read successfully
         else
            iExit = 107         ! some SPECS keywords not recognized
         end if
      end if

      call snWRAP( iExit, Solver, str, str2, iw, leniw )

  999 return

      end ! subroutine snSpec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snchkA
     &   ( iExit, nF, n, lvlChk, userfg,
     &     iGfun, jGvar, lenG, neG, x,
     &     mincw, miniw, minrw,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     userfg
      integer
     &     iExit, lenG, lencu, lencw, leniu, leniw, lenru, lenrw,
     &     lvlChk, mincw, miniw, minrw, neG, nF, n,
     &     iGfun(lenG), jGvar(lenG), iu(leniu), iw(leniw)
      double precision
     &     x(n), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     snchkA  is a stand-alone derivative checker for problems in
!     snOptA format.
!
!     snchkA is different in two ways from the built-in derivative
!     checker in snopta:

!       1. The derivatives are checked before the variables and
!          constraints are reordered to conform with snoptb format.
!
!       2. The derivatives are checked at the point defined by the
!          input argument x.  A feasible point is NOT computed before
!          the check is done.
!
!
!     lvlChk has the following meaning:
!
!       -1         do not perform any check.
!        0         do the cheap test only.
!       >0         do both cheap and full test on problem derivatives.
!
!     ------------------------------------------------------------------
!     NOTE: Before calling snchkA, there MUST be a call to the
!     initialization routine, i.e.,
!     call snInit( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )
!     This sets the default values of the optional parameters. You can
!     also alter the default values of iPrint and iSumm before snchkA
!     is used.  iPrint = 0, etc, is OK.
!     ------------------------------------------------------------------
!
!     EXIT snchkA 100 -- completed successfully
!     EXIT INFO   105 -- user-supplied derivatives appear to be correct
!     EXIT INFO   106 -- no derivatives were checked
!
!     EXIT snchkA  50 -- error in the user-supplied functions
!     EXIT INFO    55 -- incorrect derivatives
!
!     EXIT snchkA  70 -- user requested termination
!     EXIT INFO    71 -- terminated during function evaluation
!
!     EXIT snchkA  80 -- insufficient storage allocated
!     EXIT INFO    81 -- work arrays must have at least 500 elements
!     EXIT INFO    83 -- not enough integer storage
!
!     SNOPT package maintained by Philip E. Gill,
!     Dept of Mathematics, University of California, San Diego.
!
!     01 Jun 2006: First version of snchkA
!     22 Apr 2007: INFO 106 added
!     18 Feb 2015: Call-status reset on exit.
!     ==================================================================
      character
     &     Solver*6, str*80, str2*80
      integer
     &     idummy, iGmax, lElem, lF, lF1, lG, lG1, lSet,
     &     lw, lw1, lx1, ly, lz, maxcw, maxiw,
     &     maxrw, nextcw, nextiw, nextrw, npStatus
      double precision
     &     eGmax
!     ------------------------------------------------------------------
      parameter         (npStatus = 236) ! NP user-routine call-status
      parameter         (idummy   =-11111)
!     ------------------------------------------------------------------
      Solver = 'SNCHKA'
      iExit  = 0

!     ------------------------------------------------------------------
!     Check memory limits and fetch the workspace starting positions.
!     ------------------------------------------------------------------
      call s2Mem0
     &   ( iExit, Solver, lencw, leniw, lenrw, iw,
     &     mincw, miniw, minrw, maxcw, maxiw, maxrw,
     &     nextcw, nextiw, nextrw )
      if (iExit .ne. 0) go to 999

      cw(1)  = Solver//'  '

      lElem  = nextiw
      lSet   = lElem  + nF
      miniw  = lSet   + n  - 1

      lF     = nextrw
      lG     = lF     + nF
      lx1    = lG     + lenG
      lF1    = lx1    + n
      lG1    = lF1    + nF
      lw     = lG1    + lenG
      lw1    = lw     + nF
      ly     = lw1    + nF
      lz     = ly     + n
      minrw  = lz     + n  - 1

      if (miniw .gt. maxiw  .or.  minrw .gt. maxrw) then
!        ---------------------------------------------------------------
!        Not enough space to check the derivatives.
!        Exit with an (over) estimate of the additional space needed.
!        ---------------------------------------------------------------
         if (miniw .gt. maxiw) then
            write(str, 9010) miniw
            call snPRNT( 11, str, iw, leniw )
            iExit = 83
         end if

         if (minrw .gt. maxrw) then
            write(str, 9020) minrw
            call snPRNT( 11, str, iw, leniw )
            iExit = 84
         end if
         go to 800
      end if

!     ------------------------------------------------------------------
!     Go for it.
!     ------------------------------------------------------------------
      call s7checkA
     &   ( iExit, lvlChk, userfg,
     &     nF, n, iGmax, eGmax, iw(lElem), iw(lSet),
     &     iGfun, jGvar, lenG, neG,
     &     x, rw(lF), rw(lG),
     &     rw(lx1), rw(lF1), rw(lG1), rw(lw), rw(lw1), rw(ly),
     &     iw, leniw, rw, lenrw, cu, lencu,
     &     iu, leniu, ru, lenru )
      if (iExit .ne. 0) go to 800

!     Print the exit conditions.

  800 if (iExit .eq. 0) then
         iExit = 105            ! all derivatives appear to be correct
      end if

      iw(npStatus) = idummy     ! Reset call status for next solver

      call snWRAP( iExit, Solver, str, str2, iw, leniw )

  999 return

 9010 format(' Total integer   workspace  should be significantly',
     &       ' more than', i8)
 9020 format(' Total real      workspace  should be significantly',
     &       ' more than', i8)

      end ! subroutine snchkA

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snJac
     &   ( iExit, nF, n, userfg,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     x, xlow, xupp, mincw, miniw, minrw,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     userfg
      integer
     &     iExit, nF, n, neA, lenA, neG, lenG, mincw,
     &     miniw, minrw, lencu, lencw, leniu, leniw, lenru, lenrw,
     &     iAfun(lenA), jAvar(lenA), iGfun(lenG), jGvar(lenG),
     &     iu(leniu), iw(leniw)
      double precision
     &     A(lenA), x(n), xlow(n), xupp(n), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     snJac  computes the coordinates of the Jacobian.
!
!     All calculations are based on a point defined by moving the input
!     x  inside its upper and lower bounds.  snJac is terminated if
!     the problem functions are undefined at this point.
!
!     ------------------------------------------------------------------
!     NOTE: Before calling snJac, your calling program MUST call the
!     initialization routine using the call:
!     call snInit( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )
!     This sets the default values of the optional parameters. You can
!     also alter the default values of iPrint and iSumm before snJac
!     is used.  iPrint = 0, etc, is OK.
!     ------------------------------------------------------------------
!
!     EXIT snJac  100 -- completed successfully (for auxiliary routines)
!     EXIT INFO   102 -- Jacobian structure estimated
!
!     EXIT snJac  12x -- Errors while estimating Jacobian structure
!     EXIT INFO   121 -- cannot estimate Jacobian structure at given point
!     EXIT INFO    71 -- terminated during function evaluation
!     EXIT INFO    81 -- work arrays must have at least 500 elements
!     EXIT INFO    83 -- not enough integer storage
!
!     SNOPT package maintained by Philip E. Gill,
!     Dept of Mathematics, University of California, San Diego.
!
!     26 Oct 2002: First version of snJac.
!     27 Sep 2003: More thorough checks for feasibility
!     15 Jun 2008: Call-status implemented correctly.
!     18 Feb 2015: Call-status reset on exit.
!     ==================================================================
      character
     &     Solver*6, str*80, str2*80
      logical
     &     prtAll
      integer
     &     DerOpt, idummy, imaxJ, iPrt, iSum, lF, lFw, lFy, lFz, lG,
     &     lGcolw, lGcoly, lGcolz, lw, lx, lxlow, lxupp, ly, lz,
     &     maxcw, maxiw, maxrw, nextcw, nextiw, nextrw, npStatus,
     &     lcoltype, lrowtype
      double precision
     &     infBnd, emaxJ
!     ------------------------------------------------------------------
      parameter         (npStatus = 236) ! NP user-routine call-status
      parameter         (idummy   =-11111)

      double precision   zero
      parameter         (zero = 0.0d+0)
      double precision   ok
      parameter         (ok   = 0.1d+0)
!     ------------------------------------------------------------------
      Solver = 'SNJAC '
      iExit  = 0

!     ------------------------------------------------------------------
!     Check memory limits and fetch the workspace starting positions.
!     ------------------------------------------------------------------
      call s2Mem0
     &   ( iExit, Solver, lencw, leniw, lenrw, iw,
     &     mincw, miniw, minrw, maxcw, maxiw, maxrw,
     &     nextcw, nextiw, nextrw )
      if (iExit .ne. 0) go to 999

      cw(1)    = Solver//'  '

      lrowtype = nextiw
      lcoltype = lrowtype + nF
      miniw    = lcoltype + n  - 1

      lF       = nextrw
      lG       = lF     + nF
      lx       = lG     + lenG
      lxlow    = lx     + n
      lxupp    = lxlow  + n
      lw       = lxupp  + n
      ly       = lw     + n
      lz       = ly     + n
      lFw      = lz     + n
      lFy      = lFw    + nF
      lFz      = lFy    + nF
      lGcolw   = lFz    + nF
      lGcoly   = lGcolw + nF
      lGcolz   = lGcoly + nF
      minrw    = lGcolz + nF - 1

      if (miniw .gt. maxiw  .or.  minrw .gt. maxrw) then
!        ---------------------------------------------------------------
!        Not enough space to build the Jacobian.
!        Provide the user an (over) estimate of what is needed.
!        ---------------------------------------------------------------
         if (miniw .gt. maxiw) then
            write(str, 9010) miniw
            call snPRNT( 11, str, iw, leniw )
            iExit = 83
         end if

         if (minrw .gt. maxrw) then
            write(str, 9020) minrw
            call snPRNT( 11, str, iw, leniw )
            iExit = 84
         end if
         go to 800
      end if

!     ------------------------------------------------------------------
!     Go for it.
!     ------------------------------------------------------------------
      prtAll = .true.           ! ignore fixed variables for diffs

      infBnd = rw( 70)          ! definition of an infinite bound
      if (infBnd .lt. zero) then
         infBnd = 1.0d+20       ! User hasn't assigned it
      end if

!     ------------------------------------------------------------------
!     Make sure that snOptA does not let userf compute derivatives.
!     The parameters iPrt and iSum may refer to the Print and Summary
!     file respectively.  Setting them to 0 suppresses printing.
!     ------------------------------------------------------------------
      DerOpt = 0
      iPrt   = 0
      iSum   = 0
      call snSeti
     &   ( 'Derivative option', DerOpt, iPrt, iSum, iExit,
     &     cw, lencw, iw, leniw, rw, lenrw )

      call dcopy
     &   ( n, x, 1, rw(lx), 1 )
      call dcopy
     &   ( n, xlow, 1, rw(lxlow), 1 )
      call dcopy
     &   ( n, xupp, 1, rw(lxupp), 1 )

      call s7Jac
     &   ( iExit, userfg, prtAll,
     &     nF, n, infBnd, imaxJ, emaxJ,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     iw(lrowtype), iw(lcoltype), rw(lx), rw(lxlow), rw(lxupp),
     &     rw(lF), rw(lG), rw(lw), rw(ly), rw(lz),
     &     rw(lFw), rw(lFy), rw(lFz), rw(lGcolw), rw(lGcoly),rw(lGcolz),
     &     iw, leniw, cu, lencu, iu, leniu, ru, lenru )
      if (iExit .ne. 0) go to 800

      write(str, 1000) neA+neG
      call snPRNT( 13, str, iw, leniw )
      write(str, 1010) neG, neA
      call snPRNT(  3, str, iw, leniw )

      if (emaxJ .gt. ok) then
         iExit = 121            ! unable to estimate Jacobian structure
         write(str, 2000) emaxJ, imaxJ
         call snPRNT( 11, str, iw, leniw )
      end if

!     Print the exit conditions.

  800 if (iExit .eq. 0) then
         iExit = 102            ! Jacobian structure estimated
      end if

      iw(npStatus) = idummy     ! Reset call status for next solver

      call snWRAP( iExit, Solver, str, str2, iw, leniw )

  999 return

 1000 format(' Nonzero derivs  Jij  ',  i8)
 1010 format(' Non-constant    Jij''s', i8, 5x, 'Constant Jij''s',4x,i8)
 2000 format(' -->  largest error in the estimated Jacobian is',
     &                1p, e12.2, '  in row', i6)
 9010 format(' Total integer   workspace  should be significantly',
     &       ' more than', i8)
 9020 format(' Total real      workspace  should be significantly',
     &       ' more than', i8)

      end ! subroutine snJac

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snMem
     &   ( iExit, m, n, neJ, negCon,
     &     nnCon, nnJac, nnObjU,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, m, n, neJ, negCon, nnCon, nnJac, nnObjU,
     &     mincw, miniw, minrw, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     snMem   estimates the memory requirements for snOptB,
!     using the values:
!        m     , n    , neJ    negCon,
!        nnObjU, nnCon, nnJac
!
!     These values are used to compute the minimum required storage:
!     mincw, miniw, minrw.
!
!     Note:
!     1. The initialization routine snInit MUST be called before snMem.
!
!     2. All default parameters must be set before calling snMem,
!        since some values affect the amount of memory required.
!
!     3. The arrays rw and iw hold  constants and work-space addresses.
!        They must have dimension at least 500.
!
!     4. This version of snMem does not allow user accessible
!        partitions of cw, iw and rw.
!
!     Exit messages:
!
!     SNMEMB EXIT  80 -- insufficient storage allocated
!     SNMEMB INFO  81 -- work arrays must have at least 500 elements
!
!     SNMEMB EXIT 100 -- finished successfully
!     SNMEMB INFO 104 -- requirements estimated
!
!     29 Mar 1998: First version.
!     15 Oct 2003: iExit added as an argument.
!     15 Oct 2003: Current version of snMem.
!     ==================================================================

      call snMemB
     &   ( iExit, m, n, neJ, negCon,
     &     nnCon, nnJac, nnObjU,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snMem

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snMemA
     &   ( iExit, nF, n, nxname, nFname, neA, neG,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, lencw, leniw, lenrw, mincw, miniw, minrw, n, neA, neG,
     &     nF, nFname, nxname, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     snMemA   estimates the memory requirements for snOptA,
!     using the values:
!
!        nF, n, nxname, nFname, neA, neG
!
!     These values are used to compute the minimum required storage:
!     mincw, miniw, minrw.
!
!     Note:
!     1. The initialization routine snInit MUST be called before snMemA.
!
!     2. Some optional parameter settings affect the amount of memory
!        needed, so weird things may happen if some optional parameters
!        are set after the call to snMemA.
!
!     3. The arrays rw and iw hold constants and work-space addresses.
!        They must have dimension at least 500.
!
!     4. This version of snMemA does not allow user accessible
!        partitions of cw, iw and rw.
!
!     Exit messages:
!
!     SNMEMA EXIT  80 -- insufficient storage allocated
!     SNMEMA INFO  81 -- work arrays must have at least 500 elements
!
!     SNMEMA EXIT 100 -- finished successfully
!     SNMEMA INFO 104 -- requirements estimated
!
!     01 Aug 2002: First version based on snMem.
!     31 Jul 2003: snEXIT and snPRNT adopted.
!     15 Oct 2003: iExit added as an argument.
!     09 Nov 2004: Optional printing added.
!     11 Sep 2014: nnObjU and nnObj separated for FP mode.
!     ==================================================================
      integer
     &     iPrint, iSumm
!     ------------------------------------------------------------------
      iPrint = iw(12)
      iSumm  = iw(13)
      call snMemA0
     &   ( iExit, iPrint, iSumm, nF, n, nxname, nFname, neA, neG,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snMemA

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snMemA0
     &   ( iExit, lPrint, lSumm, nF, n, nxname, nFname, neA, neG,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, lPrint, lSumm, lencw, leniw, lenrw, mincw, miniw,
     &     minrw, n, neA, neG, nF, nFname, nxname, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     snMemA0   does the work for snMemA.
!
!     09 Nov 2004: First version of snMemA0
!     11 Sep 2014: nnObjU and nnObj separated for FP mode.
!     ==================================================================
      character
     &     Solver*6, str*80, str2*80, Usercw(130)*8
      logical
     &     PrintMem
      integer
     &     inform, iObjU, iPrint, iSumm, lenR, liwEst, lrwEst,
     &     llenrw, lleniw, llencw, lvlHes, maxcw, maxiw, maxrw,
     &     maxR, maxS, mQNmod, nextcw, nextiw, nextrw, nkx, nName,
     &     m, neJ, negCon, nnCon, nnJac, nnH, nnObj, nnObjU,
     &     Useriw(130)
      double precision
     &     Userrw(130)
!     ------------------------------------------------------------------
      iPrint = iw(12)
      iSumm  = iw(13)
      iw(12) = lPrint           ! Print (Log) file
      iw(13) = lSumm            ! Specs (options) file

      Solver = 'SNMEMA'
      iExit  = 0

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
!        ---------------------------------------------------------------
!        Not enough workspace to do ANYTHING!
!        Print and exit without accessing the work arrays.
!        ---------------------------------------------------------------
         iExit = 81        ! Work arrays must have at least 500 elements
         call snWRAP( iExit, Solver, str, str2, iw, leniw )
         go to 999
      end if

      cw(1)  = Solver//'  '

!     Save the user's option choices  (weird choices get overwritten).

      call chcopy( 130, cw(51), 1, Usercw, 1 )
      call icopy ( 130, iw(51), 1, Useriw, 1 )
      call dcopy ( 130, rw(51), 1, Userrw, 1 )

!     Assign fake values for lencw, leniw, lenrw.
!     This will force s2Mem to estimate the memory requirements.

      llenrw  = 500
      lleniw  = 500
      llencw  = 500

!     An obligatory call to snInit has `undefined' all options.
!     Check the user-defined values and assign undefined values.
!     s8Defaults needs various problem dimensions in iw.

!     Allocate temporary work arrays for s3size.

      nkx    = n + nF

!     Provide the user an (over) estimate of what is needed.

      neJ    = neA  + neG
      m      = nF

      if (nxname .eq. 1  .and.  nFname .eq. 1) then
         nName = 1
      else
         nName = n + m
      end if

      nnCon  = m
      nnJac  = n
      nnObjU = n
      nnObj  = n
      nnH    = n
      negCon = neJ

      iw( 15) = n      ! copy of the number of columns
      iw( 16) = m      ! copy of the nmber of rows
      iw( 17) = neJ    ! copy of the number of nonzeros in Jcol
      iw( 21) = nnJac  ! # nonlinear Jacobian variables
      iw( 22) = nnObjU ! # variables in user-defined gObj
      iw( 23) = nnCon  ! # of nonlinear constraints

      iObjU = 0        ! Temporary value
      call s8Defaults
     &   ( m, n, nnCon, nnJac, nnObjU, iObjU,
     &     cw, llencw, iw, lleniw, rw, llenrw )

      nextcw   = 501
      nextiw   = 501
      nextrw   = 501

      maxcw   = lencw
      maxiw   = leniw
      maxrw   = lenrw

      maxR    = iw( 52) ! max columns of R.
      maxS    = iw( 53) ! max # of superbasics
      mQNmod  = iw( 54) ! (ge 0) max # of BFGS updates
      lvlHes  = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian

      lenR    = maxR*(maxR + 1)/2  +  (maxS - maxR)
      nkx     = n + m

      call s8Map
     &   ( m, n, negCon, nkx, nnCon, nnJac, nnObjU, nnObj, nnH,
     &     lenR, maxR, maxS,  mQNmod, lvlHes,
     &     nextcw, nextiw, nextrw, iw, leniw )
      call s3mapA
     &   ( m, n, neJ, nF, neG, negCon, nkx, nnJac, nName,
     &     nextcw, nextiw, nextrw, iw, leniw )
      call s2Bmap
     &   ( m, n, neJ, maxS,
     &     nextiw, nextrw, maxiw, maxrw, liwEst, lrwEst, iw, leniw )
      PrintMem = .false.          ! Suppress messages from s2Mem
      call s2Mem
     &   ( inform, PrintMem, liwEst, lrwEst,
     &     nextcw, nextiw, nextrw,
     &     maxcw, maxiw, maxrw, llencw, lleniw, llenrw,
     &     mincw, miniw, minrw, iw )

!     mincw = mincw
      miniw = liwEst
      minrw = lrwEst

!     Restore the user's choices of options.

      call chcopy( 130, Usercw, 1, cw(51), 1 )
      call icopy ( 130, Useriw, 1, iw(51), 1 )
      call dcopy ( 130, Userrw, 1, rw(51), 1 )

!     Print the exit conditions.

      if (iExit .eq. 0) then
         iExit = 104            ! memory requirements estimated
      end if

      call snWRAP( iExit, Solver, str, str2, iw, leniw )

  999 iw(12) = iPrint
      iw(13) = iSumm

      end ! subroutine snMemA0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snMemB
     &   ( iExit, m, n, neJ, negCon,
     &     nnCon, nnJac, nnObjU,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, lencw, leniw, lenrw, m, mincw, miniw, minrw, n, neJ,
     &     negCon, nnCon, nnJac, nnObjU, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     snMemB   estimates the memory requirements for snoptB,
!     using the values:
!        m     , n    , neJ   negCon,
!        nnObjU, nnCon, nnJac
!
!     These values are used to compute the minimum required storage:
!     mincw, miniw, minrw.
!
!     Note:
!     1. The initialization routine snInit MUST be called before snMemB.
!
!     2. All default parameters must be set before calling snMemB,
!        since some values affect the amount of memory required.
!
!     3. The arrays rw and iw hold  constants and work-space addresses.
!        They must have dimension at least 500.
!
!     4. This version of snMemB does not allow user-accessible
!        partitions of cw, iw and rw.
!
!     Exit messages:
!
!     SNMEMB EXIT  80 -- insufficient storage allocated
!     SNMEMB INFO  81 -- work arrays must have at least 500 elements
!
!     SNMEMB EXIT 100 -- finished successfully
!     SNMEMB INFO 104 -- requirements estimated
!
!     29 Mar 1998: First version.
!     31 Jul 2003: snEXIT and snPRNT adopted.
!     15 Oct 2003: iExit added as an argument.
!     19 Feb 2004: Current version of snMemB.
!     ==================================================================
      character
     &     Solver*6, str*80, str2*80, Usercw(130)*8
      logical
     &     PrintMem
      double precision
     &     Userrw(130)
      integer
     &     inform, iObjU, lenR, liwEst, lrwEst,
     &     llenrw, lleniw, llencw, lvlHes, maxcw, maxiw, maxrw, maxR,
     &     maxS, mQNmod, nextcw, nextiw, nextrw, nkx, nnObj, nnH,
     &     Useriw(130)
!     ------------------------------------------------------------------
      Solver = 'SNMEMB'
      iExit  = 0

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
!        ---------------------------------------------------------------
!        Not enough workspace to do ANYTHING!
!        Print and exit without accessing the work arrays.
!        ---------------------------------------------------------------
         iExit = 81        ! Work arrays must have at least 500 elements
         call snWRAP( iExit, Solver, str, str2, iw, leniw )
         go to 999
      end if

      cw(1)  = Solver//'  '

!     Save the user's option choices  (weird choices get overwritten).

      call chcopy( 130, cw(51), 1, Usercw, 1 )
      call icopy ( 130, iw(51), 1, Useriw, 1 )
      call dcopy ( 130, rw(51), 1, Userrw, 1 )

!     Assign fake values for lencw, leniw, lenrw.
!     This will force s2Mem to estimate the memory requirements.

      llenrw  = 500
      lleniw  = 500
      llencw  = 500

!     An obligatory call to snInit has `undefined' all options.
!     Check the user-defined values and assign undefined values.
!     s8Defaults needs various problem dimensions in iw.

      iw( 15) = n      ! copy of the number of columns
      iw( 16) = m      ! copy of the number of rows
      iw( 17) = neJ    ! copy of the number of nonzeros in Jcol
      iw( 21) = nnJac  ! # nonlinear Jacobian variables
      iw( 22) = nnObjU ! # variables in gObj
      iw( 23) = nnCon  ! # of nonlinear constraints

      iObjU = 0
      call s8Defaults
     &   ( m, n, nnCon, nnJac, nnObjU, iObjU,
     &     cw, llencw, iw, lleniw, rw, llenrw )

      nextcw   = 501
      nextiw   = 501
      nextrw   = 501

      maxcw   = lencw
      maxiw   = leniw
      maxrw   = lenrw

      maxR    = iw( 52) ! max columns of R.
      maxS    = iw( 53) ! max # of superbasics
      mQNmod  = iw( 54) ! (ge 0) max # of BFGS updates
      lvlHes  = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian

      nnH     = max( nnJac, nnObjU )
      nnObj   = nnH             ! in case FP mode is selected
      lenR    = maxR*(maxR + 1)/2  +  (maxS - maxR)
      nkx     = n + m

      call s8Map
     &   ( m, n, negCon, nkx, nnCon, nnJac, nnObjU, nnObj, nnH,
     &     lenR, maxR, maxS,  mQNmod, lvlHes,
     &     nextcw, nextiw, nextrw, iw, leniw )
      call s2Bmap
     &   ( m, n, neJ, maxS,
     &     nextiw, nextrw, maxiw, maxrw, liwEst, lrwEst, iw, leniw )
      PrintMem = .false.          ! Suppress messages from s2Mem
      call s2Mem
     &   ( inform, PrintMem, liwEst, lrwEst,
     &     nextcw, nextiw, nextrw,
     &     maxcw, maxiw, maxrw, llencw, lleniw, llenrw,
     &     mincw, miniw, minrw, iw )

!     mincw = mincw
      miniw = liwEst
      minrw = lrwEst

!     Restore the user's choices of options.

      call chcopy( 130, Usercw, 1, cw(51), 1 )
      call icopy ( 130, Useriw, 1, iw(51), 1 )
      call dcopy ( 130, Userrw, 1, rw(51), 1 )

!     Print the exit conditions.

      if (iExit .eq. 0) then
         iExit = 104            ! memory requirements estimated
      end if

      call snWRAP( iExit, Solver, str, str2, iw, leniw )

  999 return

      end ! subroutine snMemB

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snLog
     &   ( iAbort,
     &     KTcond, mjrPrtlvl, minimize,
     &     n, nb, nnCon0, nnObj,
     &     nS, itn, nMajor, nMinor, nSwap,
     &     condZHZ, iObj, scaleObj, objAdd,
     &     fMerit, penParm, step,
     &     primalInf, dualInf, maxVi, maxViRel, hs,
     &     neJ, nlocJ, locJ, indJ, Jcol,
     &     scales, bl, bu, fCon, yCon, x,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     KTcond(2)
      integer
     &     iAbort, iObj, lencu, lencw, leniu, leniw,
     &     lenru, lenrw, mjrPrtlvl, minimize, n, neJ, nb, nlocJ, nnCon0,
     &     nnObj, nS, itn, nMajor, nMinor, nSwap, hs(nb), locJ(nlocJ),
     &     indJ(neJ), iu(leniu), iw(leniw)
      double precision
     &     condZHZ, scaleObj, objAdd, fMerit, maxViRel, maxVi,
     &     step, primalInf, dualInf, penParm(4), scales(nb),
     &     bl(nb), bu(nb), fCon(nnCon0), Jcol(neJ), yCon(nnCon0), x(nb),
     &     ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     snLog  prints the major iteration log.
!
!     The end-of-line summary is as follows:
!
!     Tag             | Tag value
!     ----------------|-------------------------------------------------
!     infoTags        | -1      0      1        2       3      4      5
!     ----------------+-------------------------------------------------
!     infoTags(QNInfo)| noBFGS BFGS   ssBFGS
!     infoTags(MdInfo)|               modA    modA+B
!     infoTags(LSInfo)|        UnLim  violLim UserLim UserRej
!     infoTags(FPInfo)|        Feas   inFeas
!     infoTags(QPInfo)|               sbLim   itMax   unbndd  weak SBlim
!     infoTags(FDInfo)|               Forwrd  Central
!     infoTags(HDInfo)| UnSet  Normal Diag    Unit
!
!
!     Tags            | printed values
!     ----------------|-------------------------------------------------
!     i   Tags(i)     | -1      0      1        2       3      4      5
!     ----------------+-------------------------------------------------
!     1  Tags(QNInfo) |  n             s
!     2  Tags(MdInfo) |                M        m
!     3  Tags(LSInfo) |                d        l       D
!     4  Tags(FPInfo) |                i
!     5  Tags(QPInfo) |                T        t       u       w     z
!     6  Tags(FDInfo) |                c
!     7  Tags(HDInfo) |                R        r
!
!     15 Nov 1991: First version.
!     19 Jul 1997: Thread-safe version.
!     02 Dec 2000: Reordered and sparsified.
!     28 Dec 2000: Row and column permutations added.
!     31 Jul 2003: snPRNT adopted.
!     11 Sep 2014: New argument nnObj.
!     ==================================================================
      external
     &     s2ColN, s2RowN
      character
     &     cflag*1, KTflag(4)*1,
     &     MjrMsg*8, buffS*85, buffP*115, str*80, form*40
      logical
     &     PrintLog, PrintC, PrintSum, NonlinearCon, NonlinearObj,
     &     prtHdg, FirstMajor, PrintHeader, prtx, prtl, prtf, prtj,
     &     scaled
      integer
     &     i, ir, Itns, j, k, k1, k2, lenL,
     &     lenU, LU, l, lvlScale, Mjrs, Mnrs, mjrHdP, mjrHdS, mline,
     &     nLine, nnCon, nnJac, nfCon(4), nfObj(4), s2ColN,
     &     s2RowN
      integer
     &     QNInfo, MdInfo, LSInfo, FPInfo, QPInfo, FDInfo, HDInfo
      double precision
     &     merit, infBnd, penNorm, signObj
!     ------------------------------------------------------------------
      integer            Diag,       Unit
      parameter         (Diag   = 1, Unit   = 2)

      integer            UnLim ,     violLim,     UserLim
      parameter         (UnLim  = 0, violLim = 1, UserLim = 2)
      integer                                     UserRej
      parameter         (                         UserRej = 3)

      integer            Scale,      UnScale
      parameter         (Scale  = 0, UnScale = 1)

      double precision   zero
      parameter         (zero   = 0.0d+0)

      parameter         (mLine  =  20)
      parameter         (mjrHdP = 224) ! >0 => Major heading for iPrint
      parameter         (mjrHdS = 226) ! >0 => Major heading for iSumm

      parameter         (QNInfo = 237) ! TagInfo(1): QN update type
      parameter         (MdInfo = 238) ! TagInfo(2):
      parameter         (LSInfo = 239) ! TagInfo(3): Line search outcome
      parameter         (FPInfo = 240) ! TagInfo(4): QP Feasibility
      parameter         (QPInfo = 241) ! TagInfo(5)  QP Optimality
      parameter         (FDInfo = 242) ! TagInfo(6)
      parameter         (HDInfo = 243) ! TagInfo(7): Approx Hessian type
!     ------------------------------------------------------------------
      character          flag(0:1)*1, key(-1:4)*2
      data               flag /' ', ' '/
      data               key  /'fr', 'lo', 'up', 'sb', 'bs', 'fx'/
!     ------------------------------------------------------------------
      nnCon     = iw( 23) ! # of nonlinear constraints
      nnJac     = iw( 21) ! # nonlinear Jacobian variables
      lvlScale  = iw( 75) ! scale option
      lenL      = iw(173) ! size of current  L
      lenU      = iw(174) ! size of current  U

      nfCon(1)  = iw(189) ! number of calls of fCon
      nfCon(2)  = iw(190) ! number of calls of fCon
      nfCon(3)  = iw(191) ! number of calls of fCon
      nfCon(4)  = iw(192) ! number of calls of fCon

      nfObj(1)  = iw(194) ! number of calls of fObj
      nfObj(2)  = iw(195) ! number of calls of fObj
      nfObj(3)  = iw(196) ! number of calls of fObj
      nfObj(4)  = iw(197) ! number of calls of fObj

      infBnd    = rw( 70) ! definition of an infinite bound

      iAbort    = 0
      LU        = lenL + lenU

      NonlinearCon = nnCon     .gt. 0
      NonlinearObj = nnObj     .gt. 0
      PrintC       = mjrPrtlvl .ge. 100
      PrintLog     = mjrPrtlvl .ge. 1
      PrintSum     = mjrPrtlvl .ge. 1
      FirstMajor   = nMajor    .eq. 0

      if (FirstMajor  .and.  PrintLog) then
         call snPRNT( 1, ' ', iw, leniw )
      end if

      penNorm = penParm(3)

      Itns   = mod( itn   , 1000000 )
      Mnrs   = mod( nMinor, 1000000 )
      Mjrs   = mod( nMajor, 1000000 )
      MjrMsg = '_'              ! Used to detect major iteration summary
      buffP  = ' '
      buffS  = ' '

      ! Add the alphabet soup.

      if (     iw(QNInfo) .eq.-1) then
         mjrMsg(2:2) = 'n'      ! No update could be made
      else if (iw(QNInfo) .eq. 1) then
         mjrMsg(2:2) = 's'      ! Scaled BFGS
      end if

      if (     iw(MdInfo) .eq. 1) then
         mjrMsg(3:3) = 'M'      ! BFGS + qN Hessian mod. 1
      else if (iw(MdInfo) .eq. 2) then
         mjrMsg(3:3) = 'm'      ! BFGS + qN Hessian mods. 1 + 2
      end if

      if (     iw(HDInfo) .eq. Diag) then
         mjrMsg(4:4) = 'R'      ! H set to a diagonal
      else if (iw(HDInfo) .eq. Unit) then
         mjrMsg(4:4) = 'r'      ! H set to the identity
      end if

      if (     iw(LSInfo) .eq. violLim) then
         mjrMsg(5:5) = 'd'      ! Violation limited via step
      else if (iw(LSInfo) .eq. UserRej) then
         mjrMsg(5:5) = 'D'      ! User rejected steps
      else if (iw(LSInfo) .eq. UserLim) then
         mjrMsg(5:5) = 'l'      ! Vars limited via step
      end if

      if (     iw(FPInfo) .eq. 1) then
         mjrMsg(6:6) = 'i'      ! QP infeasible
      end if

      if (     iw(QPInfo) .eq. 1) then
         mjrMsg(7:7) = 'T'      ! terminated by SBlimit
      else if (iw(QPInfo) .eq. 2) then
         mjrMsg(7:7) = 't'      ! terminated by itMax
      else if (iw(QPInfo) .eq. 3) then
         mjrMsg(7:7) = 'u'      ! QP unbounded
      else if (iw(QPInfo) .eq. 4) then
         mjrMsg(7:7) = 'w'      ! Weak QP solutions
      else if (iw(QPInfo) .eq. 5) then
         mjrMsg(7:7) = 'z'      ! superbasic limit reached
      end if

      if (     iw(FDInfo) .eq. 2) then
         mjrMsg(8:8) = 'c'      ! Central differences
      end if

      ! Put ( ) around small primal and dual infeasibilities.

      do k = 1, 4
         KTflag(k) = ' '
      end do

      k    = 1
      do j = 1, 2
         if (KTcond(j)) then
            KTflag(k  ) = '('
            KTflag(k+1) = ')'
         end if
         k = k + 2
      end do

      signObj = minimize
      merit   = signObj*(objAdd + fMerit)

      if (PrintLog) then
!        ------------------------------------------
!        Terse line for the Print file.
!        ------------------------------------------
         prtHdg      = iw(mjrHdP) .gt. 0
         nLine       = mod( Mjrs, mLine )
         PrintHeader = nLine  .eq. 0  .or.  prtHdg

         if (PrintHeader) then
            cflag      = flag(min( nLine, 1 ))
            iw(mjrHdP) = 0
         end if

         if (NonlinearCon) then
            if (PrintHeader) then
               call snPRNT( 11,
     &              '   Itns Major Minors    Step   nCon'
     &           // ' Feasible  Optimal  MeritFunction'
     &           // '     L+U BSwap     nS condZHZ Penalty'
     &           // cflag, iw, leniw )
            end if
            write(buffP, 3000) Itns, Mjrs, Mnrs, step, nfCon(2),
     &                         KTflag(1), primalInf, KTflag(2),
     &                         KTflag(3), dualInf, KTflag(4),
     &                         merit, LU, nSwap, nS, condZHZ, penNorm,
     &                         flag(0), MjrMsg
            if (condZHZ .eq. zero) buffP( 90: 97) = ' '
            if (penNorm .eq. zero) buffP( 98:105) = ' '

         else if (NonlinearObj) then
            if (PrintHeader) then
               call snPRNT( 11,
     &              '   Itns Major Minors    Step   nObj'
     &           // ' Feasible  Optimal      Objective'
     &           // '     L+U BSwap     nS condZHZ'
     &           // cflag, iw, leniw )
            end if
            write(buffP, 3100) Itns, Mjrs, Mnrs, step, nfObj(2),
     &                         KTflag(3), dualInf, KTflag(4),
     &                         merit, LU, nSwap, nS, condZHZ,
     &                         flag(0), MjrMsg
            if (condZHZ .eq. zero) buffP( 90: 97) = ' '
         else
            if (PrintHeader) then
               if (nS .gt. 0) then
                  call snPRNT( 11,
     &                 '   Itns Major Minors    Step       '
     &              // ' Feasible  Optimal    LPobjective'
     &              // '     L+U BSwap     nS'
     &              // cflag, iw, leniw )
               else
                  call snPRNT( 11,
     &                 '   Itns Major Minors    Step       '
     &              // ' Feasible  Optimal    LPobjective'
     &              // '     L+U             '
     &              // cflag, iw, leniw )
               end if
            end if
            write(buffP, 3200) Itns, Mjrs, Mnrs, step, nfObj(2),
     &                         KTflag(3), dualInf, KTflag(4),
     &                         merit, LU, nSwap, nS,
     &                         flag(0), MjrMsg
            buffP( 29: 35) = ' '
         end if

         if (step   .eq. zero) buffP( 21: 28) = ' '
         if (nSwap  .eq.    0) buffP( 77: 82) = ' '
         if (nS     .eq.    0) buffP( 83: 89) = ' '

         call snPRNT( 1, buffP, iw, leniw )
      end if

      if (PrintSum) then
!        --------------------------------------------
!        Terse line for the Summary file.
!        --------------------------------------------
         MjrMsg(1:1) = ' '
         prtHdg      = iw(mjrHdS)    .gt. 0

         PrintHeader = mod(Mjrs, 10) .eq. 0  .or.  prtHdg
     &                                       .or.  FirstMajor

         if (PrintHeader) then
            iw(mjrHdS) = 0
         end if

         if (NonlinearCon) then
            if (PrintHeader) then
               call snPRNT( 12,
     &              ' Major Minors     Step   nCon'
     &           // ' Feasible  Optimal  MeritFunction     nS'
     &           // ' Penalty', iw, leniw )
            end if
            write(buffS, 5000) Mjrs, Mnrs, step, nfCon(2),
     &                         KTflag(1), primalInf, KTflag(2),
     &                         KTflag(3), dualInf, KTflag(4),
     &                         merit, nS, penNorm, MjrMsg
            if (penNorm .eq. zero) buffS(70:77) = ' '

         else if (NonlinearObj) then
            if (PrintHeader) then
               call snPRNT( 12,
     &              ' Major Minors     Step   nObj'
     &           // ' Feasible  Optimal      Objective     nS',
     &              iw, leniw )
            end if
            write(buffS, 5100) Mjrs, Mnrs, step, nfObj(2),
     &                         KTflag(3), dualInf, KTflag(4),
     &                         merit, nS,         MjrMsg
         else
            if (PrintHeader) then
               if (nS .gt. 0) then
                  call snPRNT( 12,
     &                 ' Major Minors     Step       '
     &              // ' Feasible  Optimal    LPobjective     nS',
     &                 iw, leniw )
               else    ! Once zero, nS remains zero.
                  call snPRNT( 12,
     &                 ' Major Minors     Step       '
     &              // ' Feasible  Optimal    LPobjective',
     &                 iw, leniw )
               end if
            end if
            write(buffS, 5100) Mjrs, Mnrs, step, nfObj(2),
     &                         KTflag(3), dualInf, KTflag(4),
     &                         merit, nS,          MjrMsg
            buffS( 29: 35) = ' '   ! Zap nObj for LPs.
         end if

         if (step   .eq. zero) buffS(14:22) = ' '
         if (nS     .eq.    0) buffS(63:69) = ' '

         call snPRNT( 2, buffS, iw, leniw )
      end if

      if (PrintC  .and.  nnCon .gt. 0) then
!        ---------------------------------------------------------------
!        Output heading for detailed log.
!        ---------------------------------------------------------------
         call s1page( 0, iw, leniw )
         if (FirstMajor) call snPRNT( 1, ' ', iw, leniw )

!        Unscale everything if necessary.

         scaled = lvlScale .ge. 2
         if (scaled) then
            call s2applyScales
     &         ( UnScale, nnCon, n, nb, iObj, infBnd, scaleObj,
     &           neJ, nlocJ, locJ, indJ, Jcol,
     &           scales, bl, bu, yCon, x )
            call ddscl
     &         ( nnCon, scales(n+1), 1, fCon, 1 )
         end if

         i      = 0             ! Keeps ftnchek happy
         l      = mjrPrtlvl/100
         prtx   = mod( l,10 ) .gt. 0
         l      = l/10
         prtl   = mod( l,10 ) .gt. 0  .and.  Mjrs .gt. 0
         l      = l/10
         prtf   = mod( l,10 ) .gt. 0
         l      = l/10
         prtj   = mod( l,10 ) .gt. 0
         form   = '(1p, i9, e13.5)'

         if (prtx) then
            call snPRNT( 11, ' Jacobian variables', iw, leniw )
            call snPRNT(  1, ' ------------------', iw, leniw )
            do j = 1, nnJac
               write(str, form) s2ColN( j, leniw, iw ), x(j)
               call snPRNT( 1, str, iw, leniw )
            end do
         end if

         if ( prtl ) then
            call snPRNT( 11, ' Multiplier estimates', iw, leniw )
            call snPRNT(  1, ' --------------------', iw, leniw )
            do i = 1, nnCon
               write(str, form) s2RowN( i, leniw, iw ), yCon(i)
               call snPRNT( 1, str, iw, leniw )
            end do
         end if

         if ( prtf ) then
            call snPRNT( 11, ' Constraint functions', iw, leniw )
            call snPRNT(  1, ' --------------------', iw, leniw )
            do i = 1, nnCon
               write(str, form) s2RowN( i, leniw, iw ), fCon(i)
               call snPRNT( 1, str, iw, leniw )
            end do
            write(str, 7600) maxVi, maxViRel
            call snPRNT( 11, str, iw, leniw )
         end if

         if ( prtj ) then
            call snPRNT( 11, ' x  and  Jacobian', iw, leniw )
            call snPRNT(  1, ' ----------------', iw, leniw )
            do 160 j  = 1, nnJac
               l  = hs(j)
               write(str, 7410) s2ColN( j, leniw, iw ), x(j), key(l)
               call snPRNT( 11, str, iw, leniw )

               k1 = locJ(j)
               k2 = locJ(j+1) - 1
               do k  = k1, k2
                  ir = indJ(k)
                  if (ir .gt. nnCon) go to 160
                  write(str, 7420) s2RowN( indJ(k), leniw, iw ), Jcol(k)
                  call snPRNT( 1, str, iw, leniw )
               end do
  160       continue
         end if

!        Scale again if necessary.

         if (scaled) then
            call s2applyScales
     &         ( Scale, nnCon, n, nb, iObj, infBnd, scaleObj,
     &           neJ, nlocJ, locJ, indJ, Jcol,
     &           scales, bl, bu, yCon, x )
            call dddiv ( nnCon, scales(n+1), 1, fCon, 1 )
         end if
      end if
      return

!     Major log,  Print file.

 3000 format(i7, i6, i7, 1p, e8.1, i7, 1x,  2(a,e7.1,a), e14.7,
     &       i8, i6, i7, e8.1, e8.1, a1, a)
 3100 format(i7, i6, i7, 1p, e8.1, i7, 10x, 1(a,e7.1,a), e14.7,
     &       i8, i6,       i7, e8.1, a1, a)
 3200 format(i7, i6, i7, 1p, e8.1, i7, 10x, 1(a,e7.1,a), e14.7,
     &       i8, i6,       i7,       a1, a)

!     Major log,  Summary file.

 5000 format(i6, i7, 1p, e9.1, i7,  1x, 2(a,e7.1,a), e14.7, i7, e8.1, a)
 5100 format(i6, i7, 1p, e9.1, i7, 10x, 1(a,e7.1,a), e14.7, i7, a)
 7410 format(' x(', i6, ')', 1p, e13.5, 1x, a2)
 7420 format('   ', i6, ' ', 1p, e13.5)
 7600 format(' Maximum constraint violation    =', 1p, e12.4,
     &       4x, ' ( =', e11.4, ' normalized)' )

      end ! subroutine snLog

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snLog2
     &   ( probType, probTag,
     &     Elastic, GotR, FirstFeas, Feasible, JustPhase1,
     &     m, mBS, nnH, nS, jSq, jBr, jSr,
     &     linesP, linesS, itn, itQP, kPrc, lvlObjE,
     &     pivot, step, nInf, sInf, nInfE, sInfE, wtInf,
     &     nonOpt, objPrt, condZHZ, djqPrt, rgNorm, kBS, xBS,
     &     iw, leniw )

      implicit
     &     none
      character
     &     probTag*20
      logical
     &     Elastic, FirstFeas, Feasible, JustPhase1, GotR
      integer
     &     probType, m, mBS, nnH, nS, jSq, jBr, jSr, itn, itQP, kPrc,
     &     leniw, linesP, linesS, lvlObjE, nonOpt, nInf, nInfE,
     &     kBS(mBS), iw(leniw)
      double precision
     &     condZHZ, djqPrt, objPrt, pivot, rgNorm, step, sInf, sInfE,
     &     wtInf, xBS(mBS)

!     ==================================================================
!     snLog2  prints the minor iteration log.
!
!     mBS = m + maxS  .ge.  m + nS.
!
!     mnrHdP is  1  if a new heading is required for some reason other
!     than frequency (e.g., after a basis factorization).
!
!     The output consists of a number of ``sections'' of one line
!     summaries, with each section preceded by a header message.
!     linesP and linesS count the number of lines remaining to be
!     printed in each section of the print and summary files
!     respectively.   They too may force a new heading.
!
!     01 Dec 1991: First version based on Minos routine m5log.
!     17 Nov 2000: First version for SQP minor itns.
!     02 Dec 2000: Reordered and sparsified.
!     28 Dec 2000: Row and column permutations added.
!     01 Aug 2003: cgItn added to Print log.
!     04 Jul 2005: Current version of snLog2.
!     06 Dec 2014: nonOpt added as an argument.
!     29 Dec 2014: Log and summary output reorganized to print nonOpt.
!     14 May 2015: Added JustPhase1 to argument list.
!     14 May 2015: Fixed printing bug when feasible for phase 1.
!     ==================================================================
      external
     &     s1intmx, s2VarN
      character
     &     buffP*138, buffS*77, str*80
      logical
     &     NewSet,  Phase1,
     &     PrintLog, PrintSum, PrintHead, PrintHeadP, PrintHeadS
      integer
     &     cgItn, Itns, jSqN, jBrN, jSrN, k, lenL, lenU, lprDbg, lvlSys,
     &     mjrHdP, mjrHdS, mnrHdP, mnrHdS, mnrPrtlvl, minors, numInf,
     &     ncp, printP, printS, QPmode, s1intmx, s2VarN, width
      double precision
     &     mxwdth, prgNorm, rmxint, sumInf
!     ------------------------------------------------------------------
      integer            CG
      parameter         (CG     = 1)
      integer            FP,         FPE,        FPS
      parameter         (FP     = 0, FPE    = 3, FPS    = 4)
      integer            YES
      parameter         (YES    = 1)
      integer            mLineP,      mLineS
      parameter         (mLineP = 40, mLineS = 10)
      double precision   zero
      parameter         (zero   = 0.0d+0)
      parameter         (mnrHdP = 223) ! >0 => Minor heading for Log
      parameter         (mjrHdP = 224) ! >0 => Major heading for Log
      parameter         (mnrHdS = 225) ! >0 => Minor heading for Summary
      parameter         (mjrHdS = 226) ! >0 => Major heading for Summary
!     ------------------------------------------------------------------
      lvlSys    = iw( 71) ! > 0   => print system info
      lprDbg    = iw( 85) ! > 0   => private debug print
      mnrPrtlvl = iw( 93) ! Minor print level
      ncp       = iw(176) ! no. of LU compressions
      lenL      = iw(173) ! size of current  L
      lenU      = iw(174) ! size of current  U
      QPmode    = iw(208) ! Current QP solver   (based on QPslvr)
      printP    = iw(218) ! (on/off) current log     file output status
      printS    = iw(219) ! (on/off) current summary file output status
      cgItn     = iw(387) ! symmlq itns for the last QP minor itn

      PrintLog  = printP     .eq. YES
      PrintSum  = printS     .eq. YES

      mxwdth = 1.0d+7           ! Integers printed i7
      rmxint = s1intmx()        ! Largest integer without overflow
      if (mxwdth .gt. rmxint) mxwdth = rmxint
      width  = mxwdth

      Itns   = mod( itn , width )
      minors = mod( itQP, width )

      buffP  = ' '
      buffS  = ' '

      if (rgNorm .lt. 1.0d-99) then

!        Print zero if rgNorm is underflowing gradually

         prgNorm = zero
      else
         prgNorm = rgNorm
      end if

      Phase1 = .not. Feasible .or. JustPhase1

      if (Phase1) then
         numInf = nInf
         sumInf = sInf
      else
         numInf = nInfE
         sumInf = sInfE
      end if
                                ! If  newly feasible, print something.
      if (FirstFeas  .and.  mnrPrtlvl .ge. 10) then

         if (.not. Elastic) then
                                ! Constraints feasible in Normal mode.
                                ! Print a message.
                                ! probTag is one of the following:
                                ! probTag = 'QP problem'
                                ! probTag = 'LP problem'
                                ! probTag = 'QP subproblem'
                                ! probTag = 'Linear constraints'
            if (probType .ne. FPS  .and. probType .ne. FP
     &                             .and. probType .ne. FPE) then
               write(str, 8010) itn, probTag
               call snPRNT( 23, str, iw, leniw )
            end if
         else
                                ! Elastic mode
                                ! Elastic Phase 1 has completed.
            write(str, 8020) itn
            call snPRNT( 23, str, iw, leniw )
            if (lvlObjE .eq. 2) then
                                ! Infinite weight on sumInf.
                                ! Minimize the infeasible elastics.
               write(str, 8030) itn
               call snPRNT( 23, str, iw, leniw )

            else if (lvlObjE .eq. 1) then
                                ! Finite nonzero weight on sumInf
                                ! Minimize a weighted objective.
               write(str, 8040) itn
               call snPRNT( 23, str, iw, leniw )
            end if
         end if

         if (lvlSys .gt. 0) then
            iw(mnrHdP) = 1      ! Print the header to the print   file
            iw(mnrHdS) = 1      ! Print the header to the summary file
         end if
      end if

      PrintHeadP = iw(mnrHdP) .gt. 0
      PrintHeadS = iw(mnrHdS) .gt. 0

      if (PrintLog) then
!        --------------------------------------
!        Terse line for the Print file.
!        --------------------------------------
         NewSet    = linesP .eq. 0
         PrintHead = PrintHeadP  .or.  NewSet

         if (PrintHead) then
            iw(mnrHdP) = 0
            linesP     = mlineP
         end if

         iw(mjrHdP) = 1
         linesP     = linesP - 1

         jSqN  = s2VarN( jSq , leniw, iw )
         jSrN  = s2VarN( jSr , leniw, iw )
         jBrN  = s2VarN( jBr , leniw, iw )

         if (nnH .gt. 0) then
            if (PrintHead) then
               buffP = '    Itn     QP mult  QP step'
     &              // '   rgNorm          NonOpt   QP Objective'
     &              // '   +SBS   -SBS    -BS    Pivot'
     &              // '     L+U ncp    nS condZHZ'
               if (Phase1) then
                  buffP( 13: 14) = 'FP'
                  buffP( 22: 23) = 'FP'
                  buffP( 48: 53) = 'NumInf'
                  buffP( 57: 62) = ' '
                  buffP( 63: 68) = 'SumInf'
               end if

               if (Elastic) then
                  buffP( 40: 46) = 'SumInfE'
                  if (Feasible) then
                  buffP( 55: 68) = 'Elastic QP obj'
                  end if
               end if

               if (QPmode .eq. CG) buffP(126:131) = 'cgItns'
               call snPRNT( 11, buffP, iw, leniw )

            end if

            if (Phase1) then
               write(buffP, 3000) Itns, djqPrt, step,
     &                            prgNorm, sInfE, numInf, sumInf,
     &                            jSqN, jSrN, jBrN, pivot,
     &                            lenL+lenU, ncp, nS, condZHZ, cgItn
               if (numinf .le. 0) buffP(52:53) = ' '
            else
               write(buffP, 3000) Itns, djqPrt, step,
     &                            prgNorm, sInfE, nonOpt, objPrt,
     &                            jSqN, jSrN, jBrN, pivot,
     &                            lenL+lenU, ncp, nS, condZHZ, cgItn
               if (nonOpt .le. 0) buffP(52:53) = ' '
            end if

         else  ! nnH == 0
            if (PrintHead) then
               buffP = '    Itn     LP mult  LP step'
     &              // '                   NonOpt   LP Objective'
     &              // '   +SBS   -SBS    -BS    Pivot'
     &              // '     L+U ncp'

               if (Phase1) then
                  buffP( 13: 14) = 'FP'
                  buffP( 22: 23) = 'FP'
                  buffP( 48: 53) = 'NumInf'
                  buffP( 57: 62) = ' '
                  buffP( 63: 68) = 'SumInf'
               end if

               if (Elastic) then
                  if (Feasible) then
                  buffP( 55: 68) = 'Elastic LP obj'
                  end if
                  buffP( 40: 46) = 'SumInfE'
               end if

               if (nS     .gt.  0) buffP(32 : 37) = 'rgNorm'
               if (nS     .gt.  0) buffP(115:116) = 'nS'
               call snPRNT( 11, buffP, iw, leniw )
            end if

            if (Phase1) then
               write(buffP, 3000) Itns, djqPrt, step,
     &                            prgNorm, sInfE, numInf, sumInf,
     &                            jSqN, jSrN, jBrN, pivot,
     &                            lenL+lenU, ncp, nS
               if (numInf .le. 0)  buffP(52:53) = ' '
            else
               write(buffP, 3000) Itns, djqPrt, step,
     &                            prgNorm, sInfE, nonOpt, objPrt,
     &                            jSqN, jSrN, jBrN, pivot,
     &                            lenL+lenU, ncp, nS
               if (nonOpt .le. 0) buffP(52:53) = ' '
            end if
         end if

         if (djqPrt  .eq. zero) buffP( 11: 19) = ' '
         if (step    .eq. zero) buffP( 20: 28) = ' '
         if (nS      .eq.    0) then
                                buffP( 29: 37) = ' ' ! rgNorn
                                buffP(111:116) = ' ' ! nS
         end if
         if (nInfE   .eq.    0) buffP( 38: 46) = ' ' ! sInfE
         if (jSq     .eq.    0) buffP( 69: 75) = ' '
         if (jSr     .eq.    0) buffP( 76: 82) = ' '
         if (jBr     .eq.    0) buffP( 83: 89) = ' '
         if (pivot   .eq. zero) buffP( 90: 98) = ' '
         if (ncp     .eq.    0) buffP(107:110) = ' '
         if (condZHZ .eq. zero) buffP(117:124) = ' ' ! condZHZ
         if (cgItn   .eq.    0) buffP(125:131) = ' '
         call snPRNT( 1, buffP, iw, leniw )
      end if

      if (PrintSum) then
!        --------------------------------
!        Terse line for the Summary file.
!        --------------------------------
         NewSet    = linesS .eq. 0
         PrintHead = PrintHeadS  .or.  NewSet

         if (PrintHead) then
            iw(mnrHdS) = 0
            linesS = mlineS
         end if

         iw(mjrHdS) = 1
         linesS     = linesS - 1

         if (nnH .gt. 0) then
            if (PrintHead) then
               buffS = '        Minor NonOpt  QP mult  QP step'
     &              // '   rgNorm   QP objective     nS'
               if (Phase1) then
                  buffS( 15: 20) = 'NumInf'
                  buffS( 23: 24) = 'FP'
                  buffS( 32: 33) = 'FP'
                  buffS( 51: 56) = ' '
                  buffS( 57: 62) = 'SumInf'
               else
                  if (Elastic) then
                  if (lvlObjE .eq. 1) buffS(49:62) = 'Elastic QP obj'
                  if (lvlObjE .eq. 2) buffS(49:62) = '       SumInfE'
                  end if
                  if (QPmode  .eq.CG) buffS(72:77) = 'cgItns'
               end if
               call snPRNT( 12, buffS, iw, leniw )
            end if

            if (Phase1) then
               write(buffS, 5000) minors, numInf, djqPrt, step, prgNorm,
     &                            sumInf, nS, cgItn
               if (numInf .le. 0)  buffS(14:20) = ' '
            else
               write(buffS, 5000) minors, nonOpt, djqPrt, step, prgNorm,
     &                            objPrt, nS, cgItn
               if (nonOpt .le. 0) buffS(14:20) = ' '
            end if

         else  ! nnH == 0
            if (PrintHead) then
               buffS = '        Minor NonOpt  LP mult  LP Step'
     &              // '            LP objective'
               if (Phase1) then
                  buffS( 15: 20) = 'NumInf'
                  buffS( 23: 24) = 'FP'
                  buffS( 32: 33) = 'FP'
                  buffS( 51: 56) = ' '
                  buffS( 57: 62) = 'SumInf'
               else
                  if (Elastic) then
                  if (lvlObjE .eq. 1) buffS(49:62) = 'Elastic LP obj'
                  if (lvlObjE .eq. 2) buffS(49:62) = '       SumInfE'
                  end if
               end if
               call snPRNT( 12, buffS, iw, leniw )
            end if

            if (Phase1) then
               write(buffS, 5010) minors, numInf, djqPrt, step, prgNorm,
     &                            sumInf
               if (numInf .le. 0) buffS(14:20) = ' '
            else
               write(buffS, 5010) minors, nonOpt, djqPrt, step, prgNorm,
     &                            objPrt
               if (nonOpt .le. 0) buffS(14:20) = ' '
            end if
         end if

         if (djqPrt  .eq. zero)   buffS(21:29) = ' '
         if (step    .eq. zero)   buffS(30:38) = ' '
         if (nS      .eq.    0) then
                                  buffS(39:47) = ' '
                                  buffS(63:69) = ' '
         end if
         if (cgItn   .eq.    0)   buffS(70:77) = ' '

         call snPRNT( 2, buffS, iw, leniw )
      end if

!     ------------------------------------------------------------------
!     Debug output.
!     ------------------------------------------------------------------
      if (lprDbg .eq. 100) then
         call snPRNT( 11, ' BS values...', iw, leniw )
         do k = 1, m
            write(buffP, 6000) s2VarN( kBS(k), leniw,iw ), xBS(k)
            call snPRNT( 1, buffP, iw, leniw )
         end do

         call snPRNT( 11, ' SB values...', iw, leniw )
         do k = m+1, m+nS
            write(buffP, 6000) s2VarN( kBS(k), leniw,iw ), xBS(k)
            call snPRNT( 1, buffP, iw, leniw )
         end do
      end if

      return

                                ! Minor log,  Print   file.
 3000 format(1p, i7, 3x, 4e9.1, i7, e15.7,
     &          3i7, e9.1, i8, i4, i6, e8.1, i7 )
                                ! Minor log,  Summary file.
 5000 format(1p, i13, i7, 3e9.1, e15.7, i7, 1x, i7)
 5010 format(1p, i13, i7, 3e9.1, e15.7)
 6000 format(i7, g17.8)
 8010 format(  ' Itn', i7, ': Feasible ', a)
 8020 format(  ' Itn', i7, ': Feasible non-elastics')
 8030 format(  ' Itn', i7, ': Elastic Phase 2 -- minimizing',
     &                     ' elastic variables')
 8040 format(  ' Itn', i7, ': Elastic Phase 2 -- minimizing',
     &                     ' obj + weighted elastics')

      end ! subroutine snLog2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snSTOP
     &   ( iAbort,
     &     KTcond, mjrPrtlvl, minimize,
     &     m, maxS, n, nb, nnCon0, nnCon, nnObj0, nnObj, nS,
     &     itn, nMajor, nMinor, nSwap,
     &     condZHZ, iObj, scaleObj, objAdd,
     &     fMerit, penParm, step,
     &     primalInf, dualInf, maxVi, maxViRel, hs,
     &     neJ, nlocJ, locJ, indJ, Jcol, negCon,
     &     scales, bl, bu, fCon, gCon, gObj,
     &     yCon, pi, rc, rg, x,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     KTcond(2)
      integer
     &     iAbort, iObj, itn,
     &     lencu, lencw, leniu, leniw, lenru, lenrw,
     &     mjrPrtlvl, minimize, m, maxS, n, nb, neJ, negCon, nlocJ,
     &     nnCon0, nnCon, nnObj0, nnObj, nMajor, nMinor, nS, nSwap,
     &     hs(nb), locJ(nlocJ), indJ(neJ), iu(leniu), iw(leniw)
      double precision
     &     condZHZ, scaleObj, objAdd, fMerit, penParm(4),
     &     maxViRel, maxVi, step, primalInf, dualInf,
     &     scales(nb), bl(nb), bu(nb),
     &     fCon(nnCon0), gCon(negCon), gObj(nnObj0), Jcol(neJ), pi(m),
     &     rc(nb), rg(maxS), yCon(nnCon0), x(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     snSTOP is called every major iteration.
!     If iAbort > 0 on exit, the run is terminated.
!     By specifying a custom version of snSTOP, the user can arrange for
!     snopt to be terminated at any given major iteration.
!
!     14 Oct 2004: First version of   snSTOP.
!     29 Aug 2007: Parameter list extended.
!     ==================================================================

      iAbort    = 0
!     Relax

      end ! subroutine snSTOP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snSet
     &   ( buffer, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     iPrint, iSumm, Errors, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     snSet  decodes the option contained in  buffer.
!
!     The buffer is output to file iPrint, minus trailing blanks.
!     Error messages are output to files iPrint and iSumm.
!     Buffer is echoed to iPrint but normally not to iSumm.
!     It is echoed to iSumm before any error msg.
!
!     On entry,
!     iPrint is the print   file.  no output occurs if iPrint .le 0.
!     iSumm  is the Summary file.  no output occurs if iSumm  .le 0.
!     Errors is the number of errors so far.
!
!     On exit,
!     Errors is the number of errors so far.
!
!     27 Nov 1991: first version of snSet.
!     03 Nov 2000: current version.
!     ==================================================================
      integer
     &     ivalue
      double precision
     &     rvalue
      character
     &     cvalue*8, key*16
!     ------------------------------------------------------------------
      call s3opt
     &   ( .true., buffer, key, cvalue, ivalue, rvalue,
     &     iPrint, iSumm, Errors, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snSet

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snSeti
     &   ( buffer, ivalue, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     ivalue, iPrint, iSumm, Errors, lencw, leniw, lenrw,
     &     iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     snSeti decodes the option contained in  buffer // ivalue.
!     The parameters other than ivalue are as in snSet.
!
!     27 Nov 1991: first version of snSeti.
!     03 Nov 2000: current version.
!     ==================================================================
      integer
     &     ivalxx, lenbuf
      double precision
     &     rvalue
      character
     &     cvalue*8, key*16, buff72*72
!     ------------------------------------------------------------------
      write(key, '(i16)') ivalue
      lenbuf = len(buffer)
      buff72 = buffer
      buff72(lenbuf+1:lenbuf+16) = key
      ivalxx = ivalue
      call s3opt
     &   ( .true., buff72, key, cvalue, ivalxx, rvalue,
     &     iPrint, iSumm, Errors, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snSeti

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snSetr
     &   ( buffer, rvalue, iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     iPrint, iSumm, Errors, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rvalue, rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     snSetr decodes the option contained in  buffer // rvalue.
!     The parameters other than rvalue are as in snSet.
!
!     27 Nov 1991: first version of snSetr.
!     03 Nov 2000: current version.
!     ==================================================================
      integer
     &     ivalue, lenbuf
      character
     &     cvalue*8, key*16, buff72*72
      double precision
     &     rvalxx
!     ------------------------------------------------------------------
      write(key, '(1p, e16.8)') rvalue
      lenbuf = len(buffer)
      buff72 = buffer
      buff72(lenbuf+1:lenbuf+16) = key
      rvalxx = rvalue
      call s3opt
     &   ( .true., buff72, key, cvalue, ivalue, rvalxx,
     &     iPrint, iSumm, Errors, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snSetr

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer function snGet
     &   ( buffer, Errors, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     Errors, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     snGet  decodes the option contained in  buffer
!     and returns 1 if the option has previously been set, else 0.
!     For example,
!     i = snGet ( 'Maximize', Errors, cw, lencw, iw, leniw, rw, lenrw )
!
!     01 Aug 2003: First version of snGet.  Needed because
!                  snGetc, snGeti, snGetr were not well defined
!                  for strings that had no numerical value.
!     01 Aug 2003: Current version of snGet.
!     ==================================================================
      integer
     &     ivalue
      double precision
     &     rvalue
      character
     &     cvalue*8, key*16
!     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cvalue, ivalue, rvalue,
     &     0, 0, Errors, cw, lencw, iw, leniw, rw, lenrw )

      snGet  = ivalue

      end ! integer function snGet

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snGetc
     &   ( buffer, cvalue, Errors, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     Errors, lencw, leniw, lenrw, iw(leniw)
      character
     &     cvalue*8, cw(lencw)*8
      double precision
     &     rw(lenrw)

!     ==================================================================
!     snGetc gets the value of the option contained in  buffer.
!     The parameters other than cvalue are as in snSet.
!
!     17 May 1998: first version of snGetc.
!     03 Nov 2000: current version.
!     ==================================================================
      integer
     &     ivalue
      double precision
     &     rvalue
      character
     &     key*16
!     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cvalue, ivalue, rvalue,
     &     0, 0, Errors, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snGetc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snGeti
     &   ( buffer, ivalue, Errors, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     ivalue, Errors, lencw, leniw, lenrw, iw(leniw)
      character
     &     cw(lencw)*8
      double precision
     &     rw(lenrw)

!     ==================================================================
!     snGeti gets the value of the option contained in  buffer.
!     The parameters other than ivalue are as in snSet.
!
!     17 May 1998: first version of snGeti.
!     03 Nov 2000: current version.
!     ==================================================================
      double precision
     &     rvalue
      character
     &     key*16, cvalue*8
!     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cvalue, ivalue, rvalue,
     &     0, 0, Errors, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snGeti

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snGetr
     &   ( buffer, rvalue, Errors, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     Errors, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rvalue, rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     snGetr gets the value of the option contained in  buffer.
!     The parameters other than rvalue are as in snSet.
!
!     17 May 1998: first version of snGetr.
!     03 Nov 2000: current version.
!     ==================================================================
      integer
     &     ivalue
      character
     &     key*16, cvalue*8
!     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cvalue, ivalue, rvalue,
     &     0, 0, Errors, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snGetr

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snRetH
     &   ( Errors, lenH, H, nnH,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Errors, lencw, lenH, leniw, lenrw, nnH,
     &     iw(leniw)
      double precision
     &     H(lenH), rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     snRetH  retrieves the SNOPT approximate Hessian from memory.
!
!     snRetH must be called immediately after a call to npOpt, snoptB
!     or snoptC  with the option "Hessian full memory" selected.
!
!     On entry:
!        lenH      (ge 1) is the length the array H, i.e., H(lenH).
!
!     On exit
!        Errors    contains the number of errors.
!                  if Errors gt 0, then error messages are printed on
!                  the standard output
!        H(lenH)   contains the upper triangular part of the approximate
!                  Hessian, stored by columns.
!        nnH       contains the number of columns of H.
!
!
!     03 Sep 2006: First version of snRetH
!     ==================================================================
      character
     &     Solver*6, str*80, str2*80
      integer
     &     iExit, lenU,
     &     lvlHes, lU, ly, ly1
!     ------------------------------------------------------------------
      integer            FM
      parameter         (FM     = 1)
!     ------------------------------------------------------------------
      Solver = 'SNRETH'
      Errors = 0

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
!        ---------------------------------------------------------------
!        Not enough workspace to do ANYTHING!
!        Print and exit without accessing the work arrays.
!        ---------------------------------------------------------------
         iExit  = 81       ! Work arrays must have at least 500 elements
         Errors = Errors + 1
         call snWRAP( iExit, Solver, str, str2, iw, leniw )
         go to 999
      end if

      cw(1)     = Solver//'  '

!     ------------------------------------------------------------------
!     Retrieve the approximate Hessian.
!     ------------------------------------------------------------------
      lvlHes    = iw( 72)       ! LM, FM or Exact Hessian
      nnH       = iw( 24)       ! max( nnObj, nnJac )
      ly        = iw(311)       ! y (nb)      =  real work vector
      ly1       = iw(312)       ! y1(nb)      =  real work vector

      lU        = iw(391)       ! U(lenU), BFGS Hessian H = U'U
      lenU      = iw(392)       !

!     Check pointers, lengths, etc,  retrieved from memory.

      call s4chkP ( Errors, 'lvlHes', lvlHes, iw, leniw )
      call s4chkP ( Errors, '    ly',     ly, iw, leniw )
      call s4chkP ( Errors, '   ly1',    ly1, iw, leniw )
      call s4chkP ( Errors, '    lU',     lU, iw, leniw )

      if (lenH .lt. lenU) then
         Errors = Errors + 1
         write(str, 9990) lenU, lenH
         call snPRNT( 5, str, iw, leniw )
      end if

      if (lvlHes .ne. FM) then
         Errors = Errors + 1
         write(str, 9991)
         call snPRNT( 5, str, iw, leniw )
      end if

      if (Errors .eq. 0) then
         call s8getH
     &      ( nnH, lenH, rw(lU), H, rw(ly), rw(ly1) )
      end if

  999 return

 9990 format(' XXX  lenH too small: needs to be at least', i6 )
 9991 format(' XXX  Full-memory Hessian not requested' )

      end ! subroutine snRetH

