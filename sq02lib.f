!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sq02lib.f
!
!     sqTitle   sqInit   sqSpec   sqHx     sqMem    sqlog
!
!     sqSet     sqSeti       sqSetr
!     sqGet     sqGetc       sqGeti   sqGetr
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqTitle( title )

      character
     &     title*30

!     ==================================================================
!     sqtitl sets the title.
!     ==================================================================

      title  = 'S Q O P T  7.5-1.1  (Jun 2015)'
!---------------123456789|123456789|123456789|--------------------------

      end ! subroutine sqTitle

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqInit
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
!     sqInit  is called by the user to do the following:
!     1. Open default files (Print, Summary).
!     2. Initialize title.
!     3. Set options to default values.
!
!     15 Nov 1991: First version.
!     14 Jul 1997: Thread-safe version.
!     21 Mar 1997: First version based on snopt routine snInit
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
     &     maxiw, maxru, maxrw, nnCon, nnJac, nnL, nnObj, npStatus,
     &     qpStatus, s1outpt
!     ------------------------------------------------------------------
      parameter         (maxru     =   2) ! start of SNOPT part of rw
      parameter         (maxrw     =   3) ! end   of SNOPT part of rw
      parameter         (maxiu     =   4) ! start of SNOPT part of iw
      parameter         (maxiw     =   5) ! end   of SNOPT part of iw
      parameter         (maxcu     =   6) ! start of SNOPT part of cw
      parameter         (maxcw     =   7) ! end   of SNOPT part of cw
      parameter         (nnJac     =  21) ! # nonlinear Jac, variables
      parameter         (nnObj     =  22) ! # variables in gObj
      parameter         (nnCon     =  23) ! nonlinear constraints
      parameter         (nnL       =  24) ! nonlinear vars
      parameter         (lvlTim    = 182) ! Timing level
      parameter         (qpStatus  = 235) ! QP user-routine call-status
      parameter         (npStatus  = 236) ! NP user-routine call-status
      parameter         (idummy    =-11111)
!     ------------------------------------------------------------------
      character          dashes*30
      data               dashes /'=============================='/
!     ------------------------------------------------------------------
      Solver = 'SQINIT'

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

      cw(1)     = Solver//'  '

      !-----------------------------------------------------------------
      ! Initialize some default values.
      !-----------------------------------------------------------------
      iSpecs    = 0
      iStdo     = s1outpt( )
      iw( 10)   = iStdo    ! Standard Output
      iw( 11)   = iSpecs
      iw( 12)   = iPrint   ! Print file
      iw( 13)   = iSumm    ! Summary file

      iw(maxcu) = 500
      iw(maxiu) = 500
      iw(maxru) = 500
      iw(maxcw) = lencw
      iw(maxiw) = leniw
      iw(maxrw) = lenrw

      !-----------------------------------------------------------------
      ! These dimensions need to be initialized for an MPS run.
      !-----------------------------------------------------------------
      iw(nnCon) = 0
      iw(nnJac) = 0
      iw(nnObj) = 0
      iw(nnL  ) = 0

      call sqTitle( title )
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

      end ! subroutine sqInit

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqSpec
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
!     sqSpec  is called by the user to read the Specs file.
!
!     07 Feb 1998: First version.
!     01 Aug 2003: s3file now has a "title" parameter.  Use ' '.
!     27 Oct 2003: Current version of sqSpec.
!     ==================================================================
      character
     &     Solver*6, str*80, str2*80
      external
     &     s3opt
      integer
     &     Errors, Calls, iPrint, iSumm
!     ------------------------------------------------------------------
      Solver = 'SQSPEC'

      if (lencw .lt. 500 .or. leniw .lt. 500 .or. lenrw .lt. 500) then
!        ---------------------------------------------------------------
!        Not enough workspace to do ANYTHING!
!        Print and exit without accessing the work arrays.
!        ---------------------------------------------------------------
         iExit = 81        ! Work arrays must have at least 500 elements
         call snWRAP( iExit, Solver, str, str2, iw, leniw )
         go to 999
      end if

      cw(1)     = Solver//'  '

      if (iSpecs .le. 0) then
         iExit  = 131
         go to 800
      end if

      iw( 11)   = iSpecs

      iPrint    = iw( 12)
      iSumm     = iw( 13)

      iExit     = 0
      Calls     = 1

!     ------------------------------------------------------------------
!     Read the Specs file.
!     sqopt  will check the options later and maybe print them.
!     ------------------------------------------------------------------
      call s3file
     &   ( iExit, Calls, iSpecs, s3opt, ' ', iPrint, iSumm, Errors,
     &     cw, lencw, iw, leniw, rw, lenrw )

  800 if (iExit .eq. 0) then
         iExit = 101            ! SPECS file read successfully
      end if

  999 return

      end ! subroutine sqSpec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqHx
     &   ( usrHx, nnH,
     &     neH, nlocH, locH, indH, Hcol,
     &     x, Hx, sqStat,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     usrHx
      integer
     &     lencu, leniu, lenru, lencw, leniw, lenrw, neH, nlocH, nnH,
     &     sqStat, indH(neH), locH(nlocH), iu(leniu), iw(leniw)
      double precision
     &     Hcol(neH), Hx(nnH), x(nnH), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     sqHx  computes the user-defined product  Hx  and scales it.
!
!     15 Mar 1999: First   version of sqHx
!     16 Jun 2008: Call-status implemented correctly.
!     05 Nov 2014: neH, indH, locH, Hcol added as arguments.
!     ==================================================================
      logical
     &     Scaled
      integer
     &     lvlScale, lscales, lxScaled, Status
!     ------------------------------------------------------------------
      lvlScale  = iw( 75) ! scale option
      lscales   = iw(296) ! scales(nb)  = row and column scales
      lxScaled  = iw(302) ! xScaled(n)  = copy of scaled x(nnL)

      Scaled    = lvlScale .gt. 0

      ! Determine the user-function call-status.

      call s5CallStatus( Status, iw, leniw )

      if (Scaled) then
         call dcopy ( nnH,           x, 1, rw(lxScaled), 1 )
         call ddscl ( nnH, rw(lscales), 1,            x, 1 )
      end if

      call usrHx
     &   ( nnH, x, Hx, Status,
     &     cu, lencu, iu, leniu, ru, lenru )

      if (Scaled) then
         call dcopy ( nnH, rw(lxScaled), 1, x , 1 )
         call ddscl ( nnH, rw(lscales) , 1, Hx, 1 )
      end if

      end ! subroutine sqHx

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqMem
     &   ( iExit, m, n, neA,
     &     lencObj, ncolH,
     &     mincw, miniw, minrw,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, m, n, neA, lencObj, ncolH, mincw, miniw, minrw,
     &     lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     sqMem   estimates the memory requirements for sqopt,
!     using the values:
!        m       , n    , ne
!        lencObj , ncolH
!
!     These values are used to compute the minimum required storage:
!     mincw, miniw, minrw.
!
!     Note:
!     1. All default parameters must be set before calling sqMem,
!        since some values affect the amount of memory required.
!
!     2. The arrays rw and iw hold  constants and work-space addresses.
!        They must have dimension at least 500.
!
!     3. This version of sqMem does not allow user accessible
!        partitions of cw, iw and rw.
!
!     01 May 1998: First version.
!     31 Jul 2003: snPRNT adopted.
!     25 Oct 2003: Current version of sqMem.
!     ==================================================================
      character
     &     Solver*6, str*80, str2*80, Usercw(130)*8
      logical
     &     PrintMem
       integer
     &     inform, lenR, liwEst, lrwEst, llenrw, lleniw, llencw,
     &     maxcw, maxiw, maxrw, maxR, maxS, nextcw, nextiw, nextrw,
     &     nnCon, nnJac, nnObj, nkx, Useriw(130)
      double precision
     &     Userrw(130)
!     ------------------------------------------------------------------
      Solver = 'SQMEM '
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
!     This will force s5Mem to estimate the memory requirements.

      llenrw  = 500
      lleniw  = 500
      llencw  = 500

!     An obligatory call to sqInit has `undefined' all options.
!     Check the user-defined values and assign undefined values.
!     s5Defaults needs various problem dimensions in iw.

      nnCon   = 0       ! Not used in sqopt
      nnObj   = 0       ! ditto
      nnJac   = 0       ! ditto

      iw( 15) = n       ! copy of the number of columns
      iw( 16) = m       ! copy of the number of rows
      iw( 17) = neA     ! copy of the number of nonzeros in Acol

      iw( 21) = nnJac   ! # nonlinear Jacobian variables
      iw( 22) = nnObj   ! # variables in gObj
      iw( 23) = nnCon   ! # of nonlinear constraints

      iw( 26) = lencObj ! length of QP constant vector
      iw( 27) = ncolH   ! # QP Hessian columns

      call s5Defaults
     &   ( m, n, lencObj, ncolH,
     &     cw, llencw, iw, lleniw, rw, llenrw )

      nextcw   = 501
      nextiw   = 501
      nextrw   = 501

      maxcw   = lencw
      maxiw   = leniw
      maxrw   = lenrw

      maxR    = iw( 52) ! max columns of R.
      maxS    = iw( 53) ! max # of superbasics
      lenR    = maxR*(maxR + 1)/2  +  (maxS - maxR)
      nkx     = n + m

      call s5Map
     &   ( m, n, nkx, lencObj, ncolH,
     &     lenR, maxR, maxS,
     &     nextcw, nextiw, nextrw, iw, leniw )
      call s2Bmap
     &   ( m, n, neA, maxS,
     &     nextiw, nextrw, maxiw, maxrw, liwEst, lrwEst, iw, leniw )
      PrintMem = .false.        ! Print all messages in s2Mem
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

      end ! subroutine sqMem

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqLog
     &   ( probType, probTag,
     &     Elastic, GotR, FirstFeas, Feasible, JustPhase1,
     &     m, mBS, nnH, nS, jSq, jBr, jSr,
     &     linesP, linesS, itn, itQP, kPrc, lvlObjE,
     &     pivot, step, nInf, sInf, nInfE, sInfE, wtInf,
     &     nonOpt, ObjPrt, condZHZ, djqPrt, rgNorm, kBS, xBS,
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
     &     condZHZ, djqPrt, ObjPrt, pivot, rgNorm, step, sInf, sInfE,
     &     wtInf, xBS(mBS)

!     ==================================================================
!     sqLog  prints the LP/QP iteration log.
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
!     28 Dec 2000: Row and column permutations added.
!     01 Aug 2003: cgItn added to Print log.
!     26 Mar 2005: Reordered and sparsified.
!     20 Dec 2005: LP/QP Headings made consistent.
!     06 Dec 2014: nonOpt added as an argument.
!     14 May 2015: Added JustPhase1 to argument list.
!     14 May 2015: Fixed printing bug when feasible for phase 1.
!     ==================================================================
      external
     &     s1intmx, s2VarN
      character
     &     buffP*138, buffS*80, str*80
      logical
     &     NewSet, Phase1,
     &     PrintLog, PrintSum, PrintHeadP, PrintHeadS, PrintHead
      integer
     &     cgItn, Itns, jSqN, jBrN, jSrN, k, lenL, lenU, lprDbg,
     &     mnrHdP, mnrHdS, minors, ncp, numInf,
     &     PrintP, printS, QPmode, s1intmx, s2VarN, width
      double precision
     &     mxwdth, rmxint, sumInf
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
      parameter         (mnrHdS = 225) ! >0 => Minor heading for Summary
!     ------------------------------------------------------------------
      lprDbg   = iw( 85) ! > 0    => private debug print
      ncp      = iw(176) ! no. of LU compressions
      lenL     = iw(173) ! size of current  L
      lenU     = iw(174) ! size of current  U
      QPmode   = iw(208) ! Current QP solver   (based on QPslvr)
      printP   = iw(218) ! (on/off) current log     file output status
      printS   = iw(219) ! (on/off) current summary file output status
      cgItn    = iw(387) ! symmlq itns for the last QP minor itn

      PrintLog = printP .eq. YES
      PrintSum = printS .eq. YES

      mxwdth   = 1.0d+7           ! Integers printed i7
      rmxint   = s1intmx()        ! Largest integer without overflow
      if (mxwdth .gt. rmxint) mxwdth = rmxint
      width    = mxwdth

      Itns     = mod( itn , width )
      minors   = mod( itQP, width )
      buffP    = ' '
      buffS    = ' '

      Phase1   = .not. Feasible .or. JustPhase1

      if (Phase1) then
         numInf = nInf
         sumInf = sInf
      else
         numInf = nInfE
         sumInf = sInfE
      end if

                                ! If  newly feasible, print something.
      if (FirstFeas) then
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
               call snPRNT( 13, str, iw, leniw )
            end if
         else
                                ! Elastic mode
                                ! Elastic Phase 1 has completed.
            if (lvlObjE .eq. 2) then
                                ! Infinite weight on sumInf.
                                ! Minimize the infeasible elastics.
               write(str, 8030) itn
               call snPRNT( 13, str, iw, leniw )

            else if (lvlObjE .eq. 1) then
                                ! Finite nonzero weight on sumInf
                                ! Minimize a weighted objective.
               write(str, 8040) itn
               call snPRNT( 13, str, iw, leniw )
            end if
         end if

         iw(mnrHdP) = 1         ! Print the header to the print   file
         iw(mnrHdS) = 1         ! Print the header to the summary file
      end if

      PrintHeadP = iw(mnrHdP) .gt. 0
      PrintHeadS = iw(mnrHdS) .gt. 0

!     ==================================================================
!     Print (log) file
!     ==================================================================
      if (PrintLog) then          ! Print one line to the print file
         NewSet    = linesP .eq. 0
         PrintHead = PrintHeadP  .or.  NewSet

         if (PrintHead) then
            iw(mnrHdP) = 0
            linesP     = mlineP
         end if

         linesP = linesP - 1

         jSqN   = s2VarN( jSq , leniw, iw )
         jSrN   = s2VarN( jSr , leniw, iw )
         jBrN   = s2VarN( jBr , leniw, iw )

         if (nnH .gt. 0) then
            if (PrintHead) then
               buffP =     '    Itn pp  QP mult   +SBS   -SBS'
     &                  // '    -BS     Step    Pivot NonOpt'
     &                  // '    QP objective     L+U ncp'
     &                  // '  rgNorm    nS condZHZ'
               if (Phase1) then
                  buffP( 13: 14) = 'FP'
                  buffP( 60: 65) = 'NumInf'
                  buffP( 69: 75) = ' '
                  buffP( 76: 81) = 'SumInf'
               else if (Feasible .and. Elastic) then
                  if (lvlObjE .eq. 1) buffP( 69: 81) = 'Elastic QPobj'
                  if (lvlObjE .eq. 2) buffP( 69: 81) = '      SumInfE'
                  if (QPmode  .eq.CG) buffP(117:122) = 'cgItns'
               end if
               call snPRNT( 11, buffP, iw, leniw )
            end if

            if (Phase1) then
               write(buffP, 3000) Itns, kPrc, djqPrt, jSqN, jSrN, jBrN,
     &                            step, pivot, numInf, sumInf,
     &                            lenL+lenU, ncp,rgnorm,nS,condZHZ,cgItn
            else
               write(buffP, 3000) Itns, kPrc, djqPrt, jSqN, jSrN, jBrN,
     &                            step, pivot, nonOpt, ObjPrt,
     &                            lenL+lenU, ncp,rgnorm,nS,condZHZ,cgItn
            end if

         else  ! nnH == 0
            if (PrintHead) then
               buffP =     '    Itn pp  LP mult   +SBS   -SBS'
     &                  // '    -BS     Step    Pivot NonOpt'
     &                  // '    LP objective     L+U ncp'

               if (Phase1) then
                  buffP( 13: 14) = 'FP'
                  buffP( 60: 65) = 'NumInf'
                  buffP( 69: 75) = ' '
                  buffP( 76: 81) = 'SumInf'
               else
                  if (lvlObjE .eq. 1) buffP(69:81) = 'Elastic LPobj'
                  if (lvlObjE .eq. 2) buffP(69:81) = '      SumInfE'
               end if
               if (nS .gt. 0)         buffP( 94:107) = '  rgNorm    nS'
               call snPRNT( 11, buffP, iw, leniw )
            end if

            if (Phase1) then
               write(buffP, 3000) Itns, kPrc, djqPrt, jSqN, jSrN, jBrN,
     &                            step, pivot, numInf, sumInf,
     &                            lenL+lenU, ncp, rgnorm, nS
            else
               write(buffP, 3000) Itns, kPrc, djqPrt, jSqN, jSrN, jBrN,
     &                            step, pivot, nonOpt, ObjPrt,
     &                            lenL+lenU, ncp, rgnorm, nS
            end if
         end if

         if (kPrc   .eq. zero) buffP(  8: 10) = ' '
         if (djqPrt .eq. zero) buffP( 11: 19) = ' '
         if (jSq    .eq.    0) buffP( 20: 26) = ' '
         if (jSr    .eq.    0) buffP( 27: 33) = ' '
         if (jBr    .eq.    0) buffP( 34: 40) = ' '
         if (step   .eq. zero) buffP( 41: 49) = ' '
         if (pivot  .eq. zero) buffP( 50: 58) = ' '

         if (ncp    .eq.    0) buffP( 90: 93) = ' '
         if (rgnorm .eq. zero) buffP( 95:101) = ' '
         if (nS     .eq.    0) buffP(102:107) = ' '
         if (condZHZ.eq. zero) buffP(109:115) = ' ' ! condZHZ
         if (cgItn  .eq.    0) buffP(117:122) = ' '
         call snPRNT( 1, buffP, iw, leniw )
      end if

!     ==================================================================
!     Summary file
!     ==================================================================
      if (PrintSum) then        ! Print one line to the summary file

         NewSet    = linesS .eq. 0
         PrintHead = PrintHeadS  .or.  NewSet

         if (PrintHead) then
            iw(mnrHdS) = 0
            linesS     = mlineS
         end if

         linesS = linesS - 1

         if (nnH .gt. 0) then
            if (PrintHead) then
               buffS = '    Itn  QP mult  NonOpt'
     &              // '   QP objective    nS    rgNorm'

               if (Phase1) then
                  buffS( 10: 11) = 'FP'
                  buffS( 19: 24) = 'NumInf'
                  buffS( 27: 33) = ' '
                  buffS( 34: 39) = 'SumInf'
               else if (Feasible .and. Elastic) then
                  if (lvlObjE .eq. 1) buffS(27:39) = 'Elastic QPobj'
                  if (lvlObjE .eq. 2) buffS(27:39) = '      SumInfE'
                  if (QPmode  .eq.CG) buffS(58:63) = 'cgItns'
               end if
               call snPRNT( 12, buffS, iw, leniw )
            end if

            if (Phase1) then
               write(buffS, 5000) minors, djqPrt, numInf, sumInf,
     &                            nS, rgNorm, cgItn
            else
               write(buffS, 5000) minors, djqPrt, nonOpt, ObjPrt,
     &                            nS, rgNorm, cgItn
            end if

         else  ! nnH == 0
            if (PrintHead) then
               buffS = '    Itn  LP mult  NonOpt'
     &              // '   LP objective'
               if (Phase1) then
                  buffS( 10: 11) = 'FP'
                  buffS( 19: 24) = 'NumInf'
                  buffS( 27: 33) = ' '
                  buffS( 34: 39) = 'SumInf'
               else if (Feasible .and. Elastic) then
                  if (lvlObjE .eq. 1) buffS(27:39) = 'Elastic LPobj'
                  if (lvlObjE .eq. 2) buffS(27:39) = '      SumInfE'
               end if
               call snPRNT( 12, buffS, iw, leniw )
            end if

            if (Phase1) then
               write(buffS, 5010) minors, djqPrt, numInf, sumInf
            else
               write(buffS, 5010) minors, djqPrt, nonOpt, ObjPrt
            end if
         end if

         if (djqPrt .eq. zero)     buffS( 9:16) = ' '
         if (nS     .eq.    0)     buffS(40:45) = ' '
         if (rgNorm .eq. zero)     buffS(46:55) = ' '
         if (cgItn  .eq.    0)     buffS(56:63) = ' '

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

!     Minor log,  Print file.

 3000 format(1p, i7, i3, e9.1, 3i7, 2e9.1, i7, e16.8,
     &           i8, i4, e8.1, i6, e8.1, i7 )

!     Minor log,  Summary file.

!                itn  mul    ninf    obj  ns  rg     cg
 5000 format(1p, i7, e9.1, 1x, i7, e15.7, i6, e10.1, i7)
 5010 format(1p, i7, e9.1, 1x, i7, e15.7)
 6000 format(i7, g17.8)
 8010 format(  ' Itn', i7, ': Feasible ', a)
 8030 format(  ' Itn', i7, ': Elastic Phase 2 -- minimizing',
     &                     ' elastic variables')
 8040 format(  ' Itn', i7, ': Elastic Phase 2 -- minimizing',
     &                     ' obj + weighted elastics')

      end ! subroutine sqLog

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqSet
     &   ( buffer, iPrint, iSumm, iExit,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     iPrint, iSumm, iExit, lencw, leniw, lenrw,
     &     iw(leniw)
      character
     &     cw(lencw)*8
      double precision
     &     rw(lenrw)

!     ==================================================================
!     sqSet  decodes the option contained in  buffer.
!
!     The buffer is output to file iPrint, minus trailing blanks.
!     Error messages are output to files iPrint and iSumm.
!     Buffer is echoed to iPrint but normally not to iSumm.
!     It is echoed to iSumm before any error msg.
!
!     On entry,
!     iPrint is the print   file.  no output occurs if iPrint .le 0.
!     iSumm  is the Summary file.  no output occurs if iSumm  .le 0.
!     iExit  is the number of errors so far.
!
!     On exit,
!     iExit  is the number of errors so far.
!
!     27 Nov 1991: first version of sqSet.
!     20 Sep 1998: current version.
!     ==================================================================
      integer
     &     ival
      double precision
     &     rval
      character
     &     cval*8, key*16
!     ------------------------------------------------------------------
      call s3opt
     &   ( .true., buffer, key, cval, ival, rval,
     &     iPrint, iSumm, iExit,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine sqSet

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqSeti
     &   ( buffer, ivalue, iPrint, iSumm, iExit,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     ivalue, iPrint, iSumm, iExit, lencw, leniw, lenrw,
     &     iw(leniw)
      character
     &     cw(lencw)*8
      double precision
     &     rw(lenrw)

!     ==================================================================
!     sqSeti decodes the option contained in  buffer // ivalue.
!     The parameters other than ivalue are as in sqSet.
!
!     27 Nov 1991: first version of sqSeti.
!     20 Sep 1998: current version.
!     ==================================================================
      integer
     &     ival, lenbuf
      double precision
     &     rval
      character
     &     cval*8, key*16
      character
     &     buff72*72
!     ------------------------------------------------------------------
      write(key, '(i16)') ivalue
      lenbuf = len(buffer)
      buff72 = buffer
      buff72(lenbuf+1:lenbuf+16) = key
      call s3opt
     &   ( .true., buff72, key, cval, ival, rval,
     &     iPrint, iSumm, iExit,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine sqSeti

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqSetr
     &   ( buffer, rvalue, iPrint, iSumm, iExit,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     iPrint, iSumm, iExit, lencw, leniw, lenrw,
     &     iw(leniw)
      double precision
     &     rvalue, rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     sqSetr decodes the option contained in  buffer // rvalue.
!     The parameters other than rvalue are as in sqSet.
!
!     27 Nov 1991: first version of sqSetr.
!     20 Sep 1998: current version.
!     ==================================================================
      integer
     &     ival, lenbuf
      double precision
     &     rval
      character
     &     cval*8, key*16, buff72*72
!     ------------------------------------------------------------------
      write(key, '(1p, e16.8)') rvalue
      lenbuf = len(buffer)
      buff72 = buffer
      buff72(lenbuf+1:lenbuf+16) = key
      call s3opt
     &   ( .true., buff72, key, cval, ival, rval,
     &     iPrint, iSumm, iExit,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine sqSetr

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer function sqGet
     &   ( buffer, iExit, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     iExit, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     sqGet  decodes the option contained in  buffer
!     and returns 1 if the option has previously been set, else 0.
!     For example,
!     i = sqGet ( 'Maximize', iExit, cw, lencw, iw, leniw, rw, lenrw )
!
!     01 Aug 2003: First version of sqGet.  Needed because
!                  sqGetc, sqGeti, sqGetr were not well defined
!                  for strings that had no numerical value.
!     01 Aug 2003: Current version of sqGet.
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
     &     0, 0, iExit, cw, lencw, iw, leniw, rw, lenrw )

      sqGet  = ivalue

      end ! integer function sqGet

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqGetc
     &   ( buffer, cvalue, iExit, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     iExit, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cvalue*8, cw(lencw)*8

!     ==================================================================
!     sqGetc gets the value of the option contained in  buffer.
!     The parameters other than cvalue are as in sqSet.
!
!     17 May 1998: first version of sqGetc.
!     20 Sep 1998: current version.
!     ==================================================================
      integer
     &     ival
      double precision
     &     rval
      character
     &     key*16
!     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cvalue, ival, rval,
     &     0, 0, iExit,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine sqGetc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqGeti
     &   ( buffer, ivalue, iExit, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     ivalue, iExit, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     sqGeti gets the value of the option contained in  buffer.
!     The parameters other than ivalue are as in sqSet.
!
!     17 May 1998: first version of sqGeti.
!     20 Sep 1998: current version.
!     ==================================================================
      double precision
     &     rval
      character
     &     cval*8, key*16
!     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cval, ivalue, rval,
     &     0, 0, iExit,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine sqGeti

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sqGetr
     &   ( buffer, rvalue, iExit, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     buffer
      integer
     &     iExit, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rvalue, rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     sqGetr gets the value of the option contained in  buffer.
!     The parameters other than rvalue are as in sqSet.
!
!     17 May 1998: first version of sqGetr.
!     20 Sep 1998: current version.
!     ==================================================================
      integer
     &     ival
      character
     &     cval*8, key*16
!     ------------------------------------------------------------------
      call s3opt
     &   ( .false., buffer, key, cval, ival,
     &     rvalue, 0, 0, iExit,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine sqGetr

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine nullHx
     &   ( nnH, x, Hx, Status,
     &     cu, lencu, iu, leniu, ru, lenru )

      implicit
     &     none
      integer
     &     lencu, leniu, lenru, nnH, Status, iu(leniu)
      double precision
     &     x(nnH), Hx(nnH), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     This is the dummy (empty) version of the routine qpHx.
!     It should never be called by SQOPT.
!
!     Warn the user (on the standard output) that it has been called.
!     ==================================================================
      integer
     &     nOut
!     ------------------------------------------------------------------
      nOut = 6
      if (Status .eq. 1) then
         if (nOut .gt. 0) write(nOut, 1000)
      end if

      return

 1000 format(//
     &     ' XXX  The default (dummy) version of subroutine Hx',
     &     '     has been called. '
     &    /' XXX  A user-defined version is required when solving a QP')

      end ! subroutine nullHx

