!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     file  sn56qncg.f
!
!     s5QN
!     s5QNgetp  s5QNitn   s5Sswap   s5ZHZv   s5ZHZv1   s5Zswap
!     SYMMLQ    s5Msolv
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5QN
     &   ( iExit,
     &     probType, probTag, subOptimize,
     &     qpLog, Hprod, Hprod1, HvCalls,
     &     Elastic, GotR, NeedLU, typeLU, Needx,
     &     lenR, m, maxS, mBS, n, nb, nDegen,
     &     ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &     itQP, itQPmax, itQPtarget, itn,
     &     eMode, lvlObjE, printLevel,
     &     minimize, iObj, scaleObj, objAdd, objQP,
     &     condZHZ, condZmax, tolOptFP, tolOptQP, tolx,
     &     nInf, sInf, elastics, nInfE, sInfE, wtInf, piNorm,
     &     neA, nlocA, locA, indA, Acol,
     &     neH, nlocH, locH, indH, Hcol,
     &     eType, eState, feasType, hs, kBS,
     &     bl, bu, blQP, buQP, blBS, buBS,
     &     gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg, rg2,
     &     nrhs0, nrhs, rhs, scales,
     &     lenx0, nx0, x0, x, xBS, xFrozen,
     &     iy, iy1, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     Elastic, GotR, NeedLU, Needx
      character
     &     probTag*20
      external
     &     qpLog, Hprod, Hprod1
      integer
     &     HvCalls, iExit, iObj, itQP, itQPmax, itQPtarget, itn,
     &     lencu, leniu, lenru, lencw, leniw, lenrw, lenR, lenx0, eMode,
     &     lvlObjE, m, maxS, mBS, minimize, n, nb, nDegen, neA, neH,
     &     elastics, ngQP0, ngQP, ngObj0, ngObj, nInf, nInfE,
     &     nlocA, nlocH, nnH0, nnH, nS, nrhs0, nrhs, nx0, printLevel,
     &     probType, subOptimize, typeLU, locA(nlocA), locH(nlocH),
     &     indA(neA), indH(neH), eType(nb), eState(nb), hs(nb),
     &     kBS(mBS), feasType(mBS), iy(nb), iy1(nb),
     &     iu(leniu), iw(leniw)
      double precision
     &     condZHZ, condZmax, objAdd, objQP, piNorm, scaleObj,
     &     sInf, sInfE, tolOptFP, tolOptQP, tolx, wtInf, Acol(neA),
     &     bl(nb), bu(nb), blQP(nb), buQP(nb), blBS(mBS), buBS(mBS),
     &     gBS(mBS), gObj(*), gQP(ngQP0), Hcol(neH), Hdx(nnH0),
     &     pBS(mBS), pi(m), rc(nb), rhs(nrhs0), R(lenR), rg(maxS),
     &     rg2(maxS), scales(nb), x0(lenx0), x(nb), xBS(mBS),
     &     xFrozen(nb), y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

      !=================================================================
      !     s5QN   solves a linear or quadratic program.  If the objective
      !     is quadratic, a BFGS quasi-Newton method is used.
      !
      ! The optimization can pass through the following phases:
      !
      !   Phase 1               find a feasible point for all variables
      !
      !   Elastic Phase 1       make the non-elastic variables feasible
      !                         while allowing infeasible elastics
      !
      !   Phase 2               minimize the objective
      !
      !   Elastic Phase 2       minimize a composite objective while
      !                         keeping the non-elastics feasible
      !
      !                         In this phase, lvlObjE means:
      !
      !            lvlObjE = 0  zero     weight on the infeasibilities
      !                                  (infeasibillities are ignored)
      !                      1  finite   weight on the infeasibilities
      !                      2  infinite weight on the infeasibilities
      !                                  (the objective is ignored)
      !
      ! On entry:
      ! ---------
      ! elastics        is the number of elastic variables allowed to
      !                 move outside their bounds.
      ! nInfE           is the number of free elastics.
      !
      ! subOptimize =-1 means no suboptimization.
      ! subOptimize = 0 allows optimization if necessary.
      !                 During suboptimization, a subset of the
      !                 variables are frozen at their initial values.
      !                 Suboptimization is triggered if:
      !                 (1) the number of new superbasics exceeds a
      !                     preassigned limit.
      !                 (2) the number of iterations exceeds itQPTarget.
      !
      ! problemType     the type of QP being solved.
      !                 Val
      !                 ---
      !                  0  FP   feasible point only
      !                  1  LP   LP problem
      !                  2  QP   QP problem
      !                  3  FPE  feasible point for equalities only
      !                  4  FPS  feasible point for QP subProblem
      !                  5  QPS  QP subproblem
      !                  6  QPP  FP subproblem with proximal point obj
      !
      ! ngQP = max( nnH, ngObj )
      !
      ! The array kBS is a permutation on the column indices.
      ! kBS(1  :m )    holds the col. indices of the basic variables.
      ! kBS(m+1:m+nS)  holds the col. indices of the superbasic variables.
      !                These nS columns indices must have hs(j) = 2.
      !
      ! On exit:
      ! ---------
      ! subOptimize > 0 means the QP was suboptimized.
      !               1 new superbasic limit exceeded.
      !               2 itQP > itQPtarget
      !
      !  iExit       Result
      !  -----       ------
      !   >0         Fatal LU error
      !    0         QP solution found
      !   -1         QP is infeasible
      !   -2         Void
      !   -3         Too many iterations
      !   -4         Void
      !   -5         Too many superbasics
      !   -6         Void
      !   -7         Void
      !   -8         Ill-conditioned null-space basis
      !  -10         Too many subspace iterations
      !
      !
      ! 22 Jun 2001: First s5QN  based on Qpsol routine qpcore.
      ! 29 Jul 2003: Make sure R is never too ill-conditioned.
      ! 30 Jul 2003: Superbasic slacks allowed.
      ! 02 Aug 2003: snEXIT and snPRNT adopted.
      ! 07 May 2006: s4ksave handles negative values of hs.
      ! 03 Mar 2013: piNorm initialized to one instead of zero.
      ! 26 May 2013: infBnd used to identify infinite bounds.
      ! 03 Nov 2014: neH, indH, locH, Hcol added as arguments.
      ! 31 Dec 2014: rg tolerance fixed for FP phase.
      ! 02 May 2015: Basis no longer refactorized on return to phase 1.
      ! 09 May 2015: Max-time test moved inside subspace itns loop.
      !=================================================================
      character
     &     str*120
      external
     &     dnormi
      logical
     &     BadCG, Checkx, CheckFeas, Checkpi, conv(3), Converged, Done,
     &     Feasible, GotgQP, GotH, HitCon, Increase, FirstFeas,
     &     JustPhase1, LUok, ModLU, Needf, Needv, NeedR, NewB,
     &     NewLU, NewSB, Newx, Optimal, PartHess, Prtlvl10, PrintLog,
     &     PrintSum, QPdone, Stationary, UsegQP
      integer
     &     cgItn, inform, itnfix, itnlim,
     &     jq, jBq, jBr, jSq, jSr, kchk, kDegen, kfac, klog, kObj,
     &     ksav, kSumm, kp, kPrc, kPrPrt, linesP, linesS, lenL0, lenL,
     &     lenU0, lenU, LUitn, lvlTol, LUsiz0, LUmax, LUrequest, maxR,
     &     mnrHdP, mnrHdS, mNewSB, mUncon, nBS, nFac,
     &     nfix(2), nfmove, frozen, LUmod, nonOpt,
     &     nParPr, nParPrLP, nParPrQP, nSmax, nSwap,
     &     nUncon, PreCon, printP, printS, qpMode, runTime, sqStat,
     &     toldj1, toldj2, toldj3
      double precision
     &     Bgrwth, Bold, c6, deltaf, deltax, djq0, djq,
     &     djqmod, djqPrt, dnormi, dRmax, dRmin, eps0, eps2, etarg,
     &     featol, fObj0, fObj, ftoli, ftol(2), gSNorm, Hcondbnd,
     &     objPrint, infBnd, maxTime, normA, normg,
     &     pivot, pSNorm, pSNrm1, rgNorm, rgTest, rgTol(2),
     &     rowError, signObj, step, tolinc, tolrg, tolx0, weight,
     &     objSlack, xSNorm, xSNrm1, xtoli, xtol(2)
!     ------------------------------------------------------------------
      double precision   zero,             one
      parameter         (zero    = 0.0d+0, one = 1.0d+0)

      integer            LOAD
      parameter         (LOAD    = 1)
      integer            loose,       tight
      parameter         (loose   = 1, tight  = 2)
      integer            CG,          QN
      parameter         (CG      = 1, QN     = 2)
      integer            Intern,      Extern
      parameter         (Intern  = 0, Extern = 1)
      integer            Check
      parameter         (Check   = 1)
      integer            WithB
      parameter         (WithB   = 1)
      integer            FP,          FPE,        FPS
      parameter         (FP      = 0, FPE    = 3, FPS   = 4)
      integer            BS,          BT
      parameter         (BS      = 2, BT     = 3)
      integer            Init,        Optml,      Cycle
      parameter         (Init    = 0, Optml  = 1, Cycle = 2)

      parameter         (nParPr  = 94 ) ! partial pricing in use
      parameter         (nParPrLP= 99 ) ! partial pricing for LPs
      parameter         (nParPrQP= 100) ! partial pricing for QPs
      parameter         (lenL0   = 171) ! size of L0
      parameter         (lenU0   = 172) ! size of initial  U
      parameter         (lenL    = 173) ! size of current  L
      parameter         (lenU    = 174) ! size of current  U
      parameter         (toldj1  = 184) ! phase 1 dj tol for p.p.
      parameter         (toldj2  = 185) ! phase 2 dj tol for p.p.
      parameter         (toldj3  = 186) ! current optimality tol
      parameter         (kObj    = 205) ! xBS(kObj) is the obj. slack
      parameter         (LUitn   = 215) ! itns since last factorize
      parameter         (LUmod   = 216) ! number of LU mods
      parameter         (printP  = 218) ! (on/off) log     status
      parameter         (printS  = 219) ! (on/off) summary status
      parameter         (linesP  = 220) ! # lines in log     file
      parameter         (linesS  = 221) ! # lines in summary file
      parameter         (mnrHdP  = 223) ! >0 => Minor head for iPrint
      parameter         (mnrHdS  = 225) ! >0 => Minor head for iSumm
      parameter         (cgItn   = 387) ! symmlq itns for last minor
      parameter         (runTime = 462) ! Solve time
!     ------------------------------------------------------------------
      eps0      = rw(  2) ! eps**(4/5)       IEEE DP  3.00e-13
      eps2      = rw(  4) ! eps**(1/2)       IEEE DP  1.49e-08
      infBnd    = rw( 70) ! definition of an infinite bound
      maxTime   = rw( 79) ! max time allowed
      etarg     = rw( 83) ! rgNorm tolerance
      Hcondbnd  = rw( 85) ! bound on the condition of Hz

      maxR      = iw( 52) ! max columns of R
      kchk      = iw( 58) ! check (row) frequency
      kfac      = iw( 59) ! factorization frequency
      ksav      = iw( 60) ! save basis map
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      kDegen    = iw( 63) ! max. expansions of featol
      itnlim    = iw( 89) ! limit on total iterations
      mNewSB    = iw( 95) ! max # of new superbasics
      nFac      = iw(210) ! # of LU factorizations
      qpMode    = iw(208) ! Current QP solver
      PreCon    = iw(209) ! Current precon mode (based on QPslvr)

      iExit     = 0
      sqStat    = 0
      mUncon    = 20

      c6        = max( 1.0d-6, eps2 )

      if (nFac .gt. 0) then
         LUsiz0    = iw(lenL0) + iw(lenU0)
         LUmax     = 2*LUsiz0
      end if

      Prtlvl10 =        printLevel .ge. 10
      PrintLog =        printLevel .ge.  1  .and.
     &        ((mod( itQP,  klog ) .eq.  0  .and.  itQP   .ne. 0)  .or.
     &                      klog   .eq.  1                             )
      PrintSum =        printLevel .ge.  1  .and.
     &        ((mod( itQP,  kSumm) .eq.  0  .and.  itQP   .ne. 0)  .or.
     &                      kSumm  .eq.  1                             )
      iw(printP) = 0
      iw(printS) = 0
      if (PrintLog) iw(printP) = 1
      if (PrintSum) iw(printS) = 1

!     ------------------------------------------------------------------
!     s5QN operates in either ``Normal'' or ``Elastic'' mode.
!     Everything is normal unless a weighted sum is being minimized or
!     the constraints are infeasible.
!     The logical Feasible refers to the non-elastic variables.
!     Note that Feasible can be false while in elastic mode.
!     wtInf  is the optional parameter Infeasibility Weight.
!     ------------------------------------------------------------------
      Feasible   = .false.

      GotH       =  nnH .gt. 0
      GotgQP     = ngQP .gt. 0

!     JustPhase1 = stop at the end of phase 1 (in normal or elastic mode)

      JustPhase1 = probType .eq. FP  .or.
     &             probType .eq. FPE .or.
     &             probType .eq. FPS

!     The phase 2 objective is F1 + wtInf*F2.

      if (Elastic) then
         Needf = lvlObjE .ne. 2  ! F1 required in phase 2
         Needv = lvlObjE .ne. 0  ! F2 required in phase 2
      else
         Needf = .true.
         Needv = .false.
      end if

!     NeedR    = true if R needs to be set or reset.

      NeedR    =  .not.  JustPhase1
     &           .and.  (qpMode .eq. QN    .or.
     &                  (qpMode .eq. CG    .and.  PreCon .eq. 1))
      PartHess =             nS .gt. maxR  .and.  GotR

!     call dcopy ( 2, rw( ftol1), 1,  ftol, 1 )
!     call dcopy ( 2, rw( xtol1), 1,  xtol, 1 )
!     call dcopy ( 2, rw(rgTol1), 1, rgTol, 1 )

      xtol(1)  =  0.1d+0
      xtol(2)  =  1.0d-6
      ftol(1)  =  xtol(1)*0.1d+0
      ftol(2)  =  xtol(2)**2

      rgTol(loose) = 1.0d-3             ! relaxed   rgTol
      rgTol(tight) = min( tolOptQP, c6) ! stringent rgTol
!     rgTol(tight) = 1.0d-7             ! relaxed   rgTol
      lvlTol       = loose              ! working   rgTol

      condZHZ  =  zero
      objQP    =  zero
      objSlack =  zero
      pivot    =  zero
      step     =  zero

      tolrg    =  zero

      sInfE    =  zero
      nInfE    =  0
      nonOpt   = -1

      normA    =  one

      jq       =  0
      djq      =  zero
      djq0     =  zero
      djqPrt   =  zero

      jBq      =  0             ! x(jBq) is the incoming   BS
      jBr      =  0             ! x(jBr) is the outgoing   BS
      jSq      =  0             ! x(jSq) is the incoming SBS
      jSr      =  0             ! x(jSr) is the outgoing SBS
      kPrPrt   =  0
      signObj  =  minimize

      iw(cgItn ) = 0
      rw(toldj1) = 100.0d+0
      rw(toldj2) = 100.0d+0

      if ( justPhase1 ) then
         iw(nParPr) = iw(nParPrLP)
      else
         iw(nParPr) = iw(nParPrQP)
      end if

      kPrc      = 0             ! last sec scanned in part. prc
      LUrequest = 0

      CheckFeas = .true.        ! Check that x is feasible.
      Checkpi   = .true.
      Converged = .false.
      NewLU     = .true.
      Newx      = .false.
      QPdone    = .false.

!     nUncon  counts the number of unconstrained (i.e., Newton) steps.
!             If the test for a minimizer were scale-independent,
!             Uncon would never be larger than 1.
!     nfmove  counts the number of times that the QP obj is decreased,

      nfmove = 0
      nUncon = 0

!     subOptimize ne 0 implies that optimization occurs with a subset of
!     the variables frozen at their initial values.
!     During suboptimization, frozen = the number of frozen variables.

      frozen = 0
      nSmax  = nS + mNewSB
      call s5hs  ( Intern, nb, blQP, buQP, hs, x )
      call s5degen
     &   ( inform, Init, printLevel, nb, nInf, itn,
     &     featol, tolx, tolinc, hs, blQP, buQP, x,
     &     itnfix, nfix, tolx0, iw, leniw, rw, lenrw )

!!    ======================Start of main loop==========================
!+    do while (.not. QPdone  .and.  iExit .eq. 0)
  100 if       (.not. QPdone  .and.  iExit .eq. 0) then
!        ===============================================================
!        Check the initial  x  and move it onto  ( A  -I )*x = b.
!        If NeedLU is true, this will require a basis factorization.
!        ===============================================================
!        If necessary,  factorize the basis  ( B = LU ) and set x.
!        If NeedLU is false on entry to s5QN, the first call to s2Bfac
!        will try to use existing factors.
!        If NeedLU is true on entry to s5QN, an LU factorization of
!        type typeLU is computed.
!
!        The reason for the requested LU is as follows.
!
!        LUrequest =  0  First LU for a given subproblem
!        LUrequest =  1  Frequency
!        LUrequest =  2  LU nonzeros increased
!        LUrequest =  3
!        LUrequest =  4
!        LUrequest =  5  Singular after LU mod
!        LUrequest =  6  Unstable LU mod (growth in new column of U)
!        LUrequest =  7  Not enough memory
!        LUrequest =  8
!        LUrequest =  9
!        LUrequest = 10  Row error in setx
!        LUrequest = 11  Big  dx   in setx
!
!        LUrequest = 20
!        LUrequest = 21  Iterative refinement failed in QP
!        LUrequest = 22  Unbounded QP
!        LUrequest = 23  Infeasibility after refactorization
!        LUrequest = 24  Small directional derivative in QP
!        LUrequest = 25  Ill-conditioned Z in QP
!        LUrequest = 26  Indefinite Z'HZ in QP
!        LUrequest = 27  R singular after bound swap in QP
!        ---------------------------------------------------------------
         FirstFeas = .false.
         if (LUrequest .gt. 0) NeedLU = .true.

         if (Needx  .or.  NeedLU) then
            call s2Bfac
     &         ( iExit, typeLU, NeedLU, NewLU, NewB,
     &           iObj, itn, printLevel, LUrequest,
     &           m, mBS, n, nb, nnH, nS, nSwap,
     &           neA, nlocA, locA, indA, Acol,
     &           kBS, hs, blQP, buQP, blBS, buBS,
     &           nrhs0, nrhs, rhs, x, xBS,
     &           iy, iy1, y, y1,
     &           iw, leniw, rw, lenrw )
            if (NewLU) then
               LUsiz0 = iw(lenL0) + iw(lenU0)
               LUmax  = 2*LUsiz0
               if (nSwap .gt. 0) GotR = .false. ! Reset R.
               if (Prtlvl10) iw(mnrHdP) = 1
            end if

            if (iExit .ne. 0) go to 100

            PartHess  = nS .gt. maxR  .and.  GotR

            Converged = .false.
            NeedLU    = .false.
            Needx     = .false.
            Newx      = .true.
            CheckFeas = .true.
            Checkpi   = .true.  ! Check for NaNs when getting piNorm

            pivot     = zero
            nUncon    = 0
         end if

         NewSB   = .false.
         Optimal = .false.

         nBS     = m + nS

         nInf    = 0
         sInf    = zero

         if (CheckFeas) then

!           In Phase 1 or just after a factorize, check the feasibility
!           of the basic and superbasic non-elastics.

            if (Elastic) then

!              Check that the elastic variables satisfy blQP and buQP.

               call s5eReset
     &            ( nBS, nb, elastics, featol, infBnd,
     &              eType, eState, kBS,
     &              bl, bu, blQP, buQP, blBS, buBS, x )
            end if

!           FirstFeas  indicates that we have just become feasible.
!           FirstFeas is turned off once a step is taken.

            call dload ( nBS, zero, gBS, 1 )
            normg  = one
            call s5Inf
     &         ( nBS, featol, infBnd,
     &           nInf, sInf, feasType, blBS, buBS, gBS, xBS )

            if (nInf .gt. 0) then

!              Non-elastics are infeasible.
!              If necessary, switch back to the feasibility phase, after
!              refactorization.
!              Print something if the basis has just been refactorized.

               if (Prtlvl10  .and.  iw(LUitn) .eq. 0) then
                  write(str, 1030) itn, nInf, sInf
                  call snPRNT( 21, str, iw, leniw )
               end if
               Feasible = .false.
            end if

!           Feasible => the nonelastics are feasible.
!           Feasible => normal or elastic Phase 2

            if (.not. Feasible) then
               FirstFeas = nInf .eq. 0
            end if
            Feasible  = nInf .eq. 0
            CheckFeas = nInf .gt. 0
         end if ! if CheckFeas

         if (Elastic) then
            !-----------------------------------------------------------
            ! Compute the sum of infeasibilities of the elastics.
            ! If  nInfE .ne. elastics, then there is at least one
            ! nonbasic elastic fixed at its current value.
            !-----------------------------------------------------------
            call s5eInf
     &         ( nb, nBS, eState, kBS, featol, nInfE, sInfE,
     &           bl, bu, x )
         end if

         if (Feasible  .and.  JustPhase1) then
            ! The non-elastics are feasible, prepare to exit.
            condZHZ = zero
            djqPrt  = zero
            rgNorm  = zero
            call dload ( m, zero, pi, 1 )
            piNorm  = one        ! piNorm = max(norm(pi), 1.0)
         else

            if (Feasible) then
!              ---------------------------------------------------------
!              The nonelastics are feasible.
!              (Elastic = false means no elastics.)
!              ---------------------------------------------------------
!              If just feasible, compute the QP obj, grad and R.

               if (FirstFeas  .or.  Newx) then
                  if (Needf) then
!                    ===================================================
!                    Initialize the QP objective and gradient.
!                    objQP is the explicit linear plus quadratic obj.
!                    objQP is not scaled by signObj.
!                    objQP is updated after each QP step.
!                    ===================================================
                     if (GotgQP) then
                        if (HvCalls .eq. 0) sqStat = 1
                        call s5QPfg
     &                     ( Hprod, Hprod1, HvCalls,
     &                       ngQP, ngObj0, ngObj, nnH,
     &                       neH, nlocH, locH, indH, Hcol,
     &                       sqStat, objQP,
     &                       gObj, gQP, lenx0, nx0, x0, x, y,
     &                       cu, lencu, iu, leniu, ru, lenru,
     &                       cw, lencw, iw, leniw, rw, lenrw )
                        sqStat = 0
                     end if

                     if (GotH  .and.  NeedR  .and.  .not. GotR) then
!                       ------------------------------------------------
!                       Load the reduced Hessian.
!                       This happens after every LU factorize.
!                       ------------------------------------------------
                        call s6Rset
     &                     ( LOAD, maxR, nS, lenR, R, y, condZHZ )
                        GotR     = .true.
                        PartHess = nS .gt. maxR
                     end if ! GotH and not GotR
                  end if ! Needf

!                 ------------------------------------------------------
!                 Gather the QP gradient in BS order.
!                 Assign the nonzero components of gBS.
!                 ------------------------------------------------------
                  if (GotgQP) then
                     call s2gather
     &                  ( ngQP, nBS, kBS, signObj, gQP, gBS )
                  end if
                  if (iObj .gt. 0) gBS(iw(kObj)) = signObj*scaleObj

                  if (Elastic  .and.  nInfE .gt. 0  .and.  Needv) then
                     call s5eGrad
     &                  ( nb, nBS, wtInf, eState, kBS, gBS )
                  end if
               end if ! FirstFeas .or. Newx

!              ---------------------------------------------------------
!              See if its time to suboptimize.
!              No suboptimization if all steps have been degenerate.
!              ---------------------------------------------------------
               if (subOptimize .ne. 0  .or.  nfmove .eq. 0) then
!                 Relax
               else
                  if (nS  .ge. nSmax) then
                     subOptimize = 1
                     if (Prtlvl10) then
                        write(str, 1610) itn, mNewSB
                        call snPRNT( 21, str, iw, leniw )
                     end if
                  else if (itQP .ge. itQPtarget) then
                     subOptimize = 2
                     if (Prtlvl10) then
                        write(str, 1620) itn, itQPtarget
                        call snPRNT( 21, str, iw, leniw )
                     end if
                  end if
               end if
            end if ! Feasible

!           ============================================================
!           Check for an approximate stationary point.
!           If the gradient has changed, compute the reduced gradient
!           ============================================================
            if (.not. Feasible  .or.  FirstFeas  .or.  Newx) then
               normg = max(dnormi( nBS, gBS, 1 ), one)
               call dcopy ( m, gBS, 1, y, 1 )
               call s5setpi
     &            ( inform,
     &              m, Checkpi, condZmax, normg, piNorm, y, pi,
     &              iw, leniw, rw, lenrw )
               Checkpi = .false.

               if (inform .ne. 0) then
                  if (inform .gt. 0) then
                     iExit     =  inform
                  else          ! pi is infinite or contains a NaN/Inf.
                     LUrequest = -inform
                     call s2tryLU
     &                  ( itn, LUrequest, nS, LUrequest, LUok, typeLU,
     &                    iw, leniw, rw, lenrw )
                     if (.not. LUok) iExit = 44 ! Large cond(Z) estimate
                  end if
                  go to 100
               end if

               rgNorm = zero
               if (nS .gt. 0) then
                  call s5rg
     &               ( m, nBS, n, nS, eps0,
     &                 neA, nlocA, locA, indA, Acol,
     &                 gBS, pi, rg, rgNorm, kBS )
               end if
            end if ! .not. Feasible  .or.  FirstFeas  .or.  Newx

            if (tolrg .eq. zero) tolrg = etarg*rgNorm

            if (Feasible) then
               rw(toldj3) = tolOptQP
            else
               rw(toldj3) = tolOptFP
            end if

            if (PartHess) then
               gSNorm = dnormi( maxR, rg, 1 )
            else
               gSNorm = rgNorm
            end if

            rgTest = max( piNorm, normg )

            if (Feasible) then
               Stationary =         gSNorm .le.     0.1d+0  *tolrg
     &                      .or.    gSNorm .le. rgTol(tight)*rgTest
     &                      .or. Converged
            else
               Stationary = gSNorm .le. rgTol(tight)*rgTest
            end if

            if (PartHess  .and.  Stationary) then

!              Swap the largest reduced gradient in Z2 into the front of
!              Z2 and see if it is significantly large.
!              Reset Stationary if necessary.

               call s5Zswap
     &            ( Stationary, m, maxR, maxS, nBS, lenR, nS,
     &              tolrg, kBS, blBS, buBS, gBS, R, rg, xBS, rw, lenrw )
            end if

            if (GotR) then
               call s6Rcnd
     &            ( maxR, nS, lenR, R, dRmax, dRmin, condZHZ )
               if (condZHZ .gt. Hcondbnd) then
                  call s6Rset
     &               ( LOAD, maxR, nS, lenR, R, y, condZHZ )
               end if
            end if

            kPrPrt = kPrc
            jq     = 0

            if (Stationary) then
!              ---------------------------------------------------------
!              Compute Lagrange multipliers.
!              ---------------------------------------------------------
               djq0   = djq     ! save djq in case of bad Stationary
               djq    = zero
               nUncon = 0
               UsegQP = Feasible  .and.  GotgQP
               weight = zero
               if (Elastic  .and.  Feasible) then
                  weight = wtInf
               end if

               call s5price
     &            ( Elastic, Feasible, Increase, UsegQP, subOptimize,
     &              itn, m, n, nb, ngQP0, ngQP,
     &              frozen, nonOpt, weight, signObj, piNorm,
     &              jq, djq, kPrc, rw(toldj1),
     &              neA, nlocA, locA, indA, Acol,
     &              eType, hs, gQP, pi, rc, x, xFrozen,
     &              iw, leniw, rw, lenrw )

               Optimal = nonOpt .eq. 0
               NewSB   = nonOpt .gt. 0

               if ( lvlTol .eq. loose .and.
     &             (nS     .ge. maxS  .or.  Optimal)) then
                  lvlTol = tight
                  tolrg  = rw(toldj3)*piNorm
                  if ( Prtlvl10 ) then
                     write(str, 1700) itn, tolrg
                     call snPRNT( 21, str, iw, leniw )
                  end if
                  if (rgNorm .gt. tolrg) then
                     Optimal = .false.
                     NewSB   = .false.
                  end if
               end if
            end if ! Stationary
         end if ! JustPhase1

         QPdone = Optimal .or. (Feasible .and. JustPhase1)

         if (QPdone) then
!           ------------------------------------------------------------
!           Apparently we are finished.
!           See if any nonbasics have to be set back on their bounds.
!           ------------------------------------------------------------
            call s5degen
     &         ( inform, Optml, printLevel, nb, nInf, itn,
     &           featol, tolx, tolinc, hs, blQP, buQP, x,
     &           itnfix, nfix, tolx0,
     &           iw, leniw, rw, lenrw )

            QPdone = inform .eq. 0

            if (QPdone) then
!              ---------------------------------------------------------
!              So far so good.  Now check the row residuals.
!              ---------------------------------------------------------
               if (iw(LUitn) .gt. 0) then
                  call s5setx
     &               ( inform, Check, itn,
     &                 m, n, nb, nBS, rowError,
     &                 neA, nlocA, locA, indA, Acol,
     &                 kBS, xBS, nrhs0, nrhs, rhs, x, y, y2,
     &                 iw, leniw, rw, lenrw )

                  QPdone    = inform .eq. 0
                  LUrequest = inform
                  if (LUrequest .gt. 0) typeLU = BS
               end if
            end if

            if (.not. QPdone) then
               Needx     = .true.
               Converged = .false.
               go to 100
            end if

            if (FirstFeas  .and.  JustPhase1) then
!              Relax, we are about to exit without printing anything.
            else
!              =========================================================
!              Print the details of the final iteration.
!              =========================================================
               objPrint = zero
               if (Feasible .and. .not. JustPhase1) then
                  if (Needf) then
                     if (iObj .ne. 0) then
                        objSlack = xBS(iw(kObj))*scaleObj
                     end if
                     objPrint = objAdd + objSlack + objQP
                  end if
                  if (Needv) then
                     objPrint = objPrint + signObj*wtInf*sInfE
                  end if
               end if

               call qpLog
     &            ( probType, probTag,
     &              Elastic, GotR, FirstFeas, Feasible, JustPhase1,
     &              m, mBS, nnH, nS, jSq, jBr, jSr,
     &              iw(linesP), iw(linesS), itn, itQP, kPrPrt, lvlObjE,
     &              pivot, step, nInf, sInf, nInfE, sInfE, wtInf,
     &              nonOpt, objPrint, condZHZ, djqPrt, rgNorm, kBS, xBS,
     &              iw, leniw )
            end if

            jBq       = 0
            jBr       = 0
            jSq       = 0
            jSr       = 0
            kPrPrt    = 0
            iw(cgItn) = 0
            djqPrt    = zero
            condZHZ   = zero

!           ------------------------------------------------------------
!           Convergence.
!           ------------------------------------------------------------
            if (nInf .gt. 0) then

!              No feasible point.
!              Stop or continue in elastic mode, depending on the
!              specified level of infeasibility.

               if (eMode .eq. 0  .or.  Elastic) then
!                 ------------------------------------------------------
!                 The nonelastic constraints cannot be satisfied.
!                 ------------------------------------------------------
                  iExit = -1    ! Infeasible nonelastics

               else
!                 ------------------------------------------------------
!                 Infeasible nonelastics in Normal mode.
!                 Print a message and start elastic Phase 1.
!                 ------------------------------------------------------
                  if (Prtlvl10) then
                     write(str, 8050) itn, probTag
                     call snPRNT( 23, str, iw, leniw )
                     write(str, 8060) itn
                     call snPRNT( 23, str, iw, leniw )
                     iw(mnrHdP) = 1
                     iw(mnrHdS) = 1
                  end if

                  Elastic   = .true.
                  CheckFeas = .true. ! call s5eReset
                  QPdone    = .false.

                  Needf   = lvlObjE .ne. 2 ! F1 used in phase 2
                  Needv   = lvlObjE .ne. 0 ! F2 used in phase 2
                  djq     = zero
                  step    = zero
               end if
               go to 100
            end if

            if (Prtlvl10 .and. .not. JustPhase1) then
               if (jq .ne. 0) then
                  djqprt = signObj*djq
                  if (Prtlvl10  .and.  klog .eq. 1) then
                     write(str, 1010) djq, jq, rgnorm, piNorm
                     call snPRNT( 31, str, iw, leniw )
                  end if
               else
                  if (Prtlvl10  .and.  klog .eq. 1) then
                     write(str, 1020)          rgnorm, piNorm
                     call snPRNT( 31, str, iw, leniw )
                  end if
               end if
            end if
         else ! not QPdone
!           ============================================================
!           Take a series of reduced gradient steps until the new point
!           either approximately minimizes the objective on the working
!           set or moves to the boundary of a new constraint.
!+          ============================================================
!           Repeat            (until  HitCon or Converged or Stationary)

  500          Converged   = .false.
               objPrint    = zero
               if (Feasible) then
                  if (Needf) then
                     if (iObj .ne. 0) then
                        objSlack = xBS(iw(kObj))*scaleObj
                     end if
                     objPrint = objAdd + objSlack + objQP
                  end if
                  if (Needv) then
                     objPrint = objPrint + signObj*wtInf*sInfE
                  end if
               end if

               call qpLog
     &            ( probType, probTag,
     &              Elastic, GotR, FirstFeas, Feasible, JustPhase1,
     &              m, mBS, nnH, nS, jSq, jBr, jSr,
     &              iw(linesP), iw(linesS), itn, itQP, kPrPrt, lvlObjE,
     &              pivot, step, nInf, sInf, nInfE, sInfE, wtInf,
     &              nonOpt, objPrint, condZHZ, djqPrt, rgNorm, kBS, xBS,
     &              iw, leniw )
               jBq       = 0
               jBr       = 0
               jSq       = 0
               jSr       = 0
               kPrPrt    = 0
               iw(cgItn) = 0
               djqPrt    = zero
               condZHZ   = zero

               if (NewSB) then
!                 ------------------------------------------------------
!                 A nonbasic has been selected to become superbasic.
!                 Compute the vector y such that B y = column jq.
!                 ------------------------------------------------------
!                 Set the level to which rgNorm must be reduced in the
!                 current subspace before we consider moving off another
!                 constraint.

                  djqmod = abs( djq )

                  rgNorm = max( rgNorm, djqmod )
                  tolrg  = etarg*djqmod

                  if (nS+1 .gt. maxS) then
                     iExit = -5
                     go to 100
                  end if

                  djqPrt = djq

!                 ------------------------------------------------------
!                 Compute the vector pBS such that B pB = column jq.
!                 pBS is a multiple of part of the new column of  Z  and
!                 is used to define the QN search direction.
!                 ------------------------------------------------------
!                 Unpack column jq into  y1  and solve  B*pB = y1.
!                 The altered  y1  satisfies  L*y1 = ajq.
!                 It is used later in s5QNitn to modify L and U.

                  call s2unpack
     &               ( jq, m, n, neA, normA, nlocA, locA, indA, Acol,y1)
                  call s2Bsol
     &               ( iExit, WithB, m, y1, pBS, iw, leniw, rw, lenrw )
                  if (iExit .ne. 0) return
               end if ! NewSB

               if (itn  .ge. itnlim  .or.  itQP .ge. itQPmax) then
                  iExit = -3
                  go to 100
               end if

               itQP   = itQP   + 1
               itn    = itn    + 1

!              Decide if we want to print something this iteration.

               PrintLog = printLevel .ge. 1 .and. mod(itQP,klog ) .eq. 0
               PrintSum = printLevel .ge. 1 .and. mod(itQP,kSumm) .eq. 0

               iw(printP) = 0
               iw(printS) = 0
               if (PrintLog) iw(printP) = 1
               if (PrintSum) iw(printS) = 1

               if (Feasible) then
                  fObj0 = zero
                  if (Needf) then
                     fObj0 = signObj*(objSlack + objQP)
                  end if
                  if (Needv) then
                     fObj0 = fObj0 + wtInf*sInfE
                  end if
               end if

!              ---------------------------------------------------------
!              Take a ``reduced gradient'' step.
!              ---------------------------------------------------------
               call s5QNitn
     &            ( inform,
     &              Hprod, Hprod1, HvCalls,
     &              Elastic, Feasible, GotgQP, GotH, GotR,
     &              HitCon, Increase, ModLU, Needf, Needv, NewSB,
     &              itn, lenR, m, mBS, maxR, maxS,
     &              n, nb, nnH0, nnH, nS, ngQP0, ngQP, nDegen,
     &              elastics, LUrequest,
     &              kp, jBq, jSq, jBr, jSr,
     &              jq, lvlTol, nfmove, nUncon,
     &              djq0, djq, minimize, iObj, scaleObj, objQP,
     &              condZHZ, condZmax, featol, pivot,
     &              piNorm, rgNorm, step, tolinc,
     &              nInfE, sInfE, wtInf, pSNorm, pSNrm1, xSNorm, xSNrm1,
     &              neA, nlocA, locA, indA, Acol,
     &              neH, nlocH, locH, indH, Hcol,
     &              eType, eState, feasType, hs, kBS,
     &              bl, bu, blQP, buQP, blBS, buBS, gBS,
     &              gQP, Hdx, pBS, pi, rg, rg2, R,
     &              x, xBS, y, y1, y2,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )

!              inform  Meaning                             Action
!              ------  -------                             ------
!                -5    Refinement step hit deleted bound.  (LUrequest)
!                -4    Bad directional derivative.         (LUrequest)
!                -3    Unable to update the LU factors.    (LUrequest)
!                -2    Z'HZ not positive semidefinite      (exit)
!                -1    Unbounded                           (exit)
!                 0    Normal exit
!                >0    Fatal LU error                      (exit)

               if (inform .ne. 0) then
                  if      (inform .gt.  0) then
                     iExit = inform ! Fatal LU error
                  else if (inform .eq. -1) then

!                    Unbounded. As Z'HZ is positive definite, this can
!                    happen only because the norm of p must be very big.
!                    Refactor B, possibly with a BS factorize.

                     if (iw(LUitn) .eq. 0) then
                        iExit = -2 ! unbounded
                     else
                        call s6Rset
     &                     ( LOAD, maxR, nS, lenR, R, y, condZHZ )
                        GotR = .true.

                        call s2tryLU
     &                     ( itn, 25, nS, LUrequest, LUok,typeLU,
     &                       iw, leniw, rw, lenrw )
                        if (.not. LUok) then
                           iExit = 43
                        end if
                        inform = 0
                     end if
                  end if
                  go to 100
               end if

               nBS  = m + nS

!              ---------------------------------------------------------
!              Stop if there are errors.
!              ---------------------------------------------------------
               Done  = LUrequest .ne. 0
               BadCG = nUncon .ge. mUncon  .and.  lvlTol .eq. tight

               if (.not. Done) then
                  if (.not. HitCon) then
!                    ---------------------------------------------------
!                    The QP step was unconstrained.
!                    Test for convergence on the current working set.
!                    ---------------------------------------------------
                     PartHess  = nS .gt. maxR  .and.  GotR
                     gSNorm    = rgNorm

                     if (PartHess) then
                        pSNorm = pSNrm1
                        xSNorm = xSNrm1
                        gSNorm = dnormi( maxR, rg, 1 )
                     end if

                     if (Feasible) then
                        fObj = zero
                        if (Needf) then
                           if (iObj .gt. 0) then
                              objSlack = xBS(iw(kObj))*scaleObj
                           end if
                           fObj = signObj*(objSlack + objQP)
                        end if

                        if (Needv) then
                           fObj = fObj + wtInf*sInfE
                        end if

                        ftoli   = ftol(lvlTol)
                        xtoli   = xtol(lvlTol)

                        deltax  = step * pSNorm
                        deltaf  = fObj0 - fObj

                        conv(1) = deltax .le. xtoli*(one + xSnorm )
                        conv(2) = deltaf .le. ftoli*(one + abs(fObj0))
                        conv(3) = gSNorm .le. tolrg

                        Converged = conv(1) .and. conv(2) .and. conv(3)
                     end if

                     Stationary =  gSNorm .le.       0.1d+0*tolrg  .or.
     &                             gSNorm .le. rgTol(tight)*piNorm
                  end if ! .not. HitCon
               end if

!              Check time limit

               if (maxTime .gt. zero .and. mod(itn, 20) .eq. 0) then
                  call s1time ( -2, 0, iw, leniw, rw, lenrw )
                  call s1time (  2, 0, iw, leniw, rw, lenrw )
                  if (rw(runTime) .gt. maxTime) then
                     iExit = 34 ! time limit reached
                     go to 100
                  end if
               end if

!+          until    (HitCon .or. Stationary .or. Converged
!                            .or.       Done .or. BadCG    )
            if (.not.(HitCon .or. Stationary .or. Converged
     &                       .or.       Done .or. BadCG)   )
     &      go to 500
!           ------------------------------------------------------------

            if (BadCG) then
               iExit = -10       ! too many subspace iterations
               go to 100
            end if

            if (LUrequest .gt. 0) then
               call s2tryLU
     &            ( itn, LUrequest, nS, LUrequest, LUok, typeLU,
     &              iw, leniw, rw, lenrw )
               if (.not. LUok) then
                  iExit = 43
                  go to 100
               end if
            end if

            iw(LUitn) = iw(LUitn)  + 1 ! itns since the last LU
            NewLU     = .false.
            Newx      = .false.
            Checkpi   = .false.

!           ============================================================
!           Test for error condition and/or frequency interrupts.
!           ============================================================
            if (mod(itn,ksav) .eq. 0) then
               call s4ksave
     &            ( minimize, m, n, nb, nS, mBS,
     &              itn, nInf, sInf, objQP, kBS, hs,
     &              scales, blQP, buQP, x, xBS, cw, lencw, iw, leniw )
            end if

!           Increment featol every iteration.

            featol = featol + tolinc

!           Every kdegen iterations, reset featol and move nonbasic
!           variables onto their bounds if they are very close.

            if (mod( itn, kdegen ) .eq. 0) then
               call s5degen
     &            ( inform, Cycle, printLevel, nb, nInf, itn,
     &              featol, tolx, tolinc, hs, blQP, buQP, x,
     &              itnfix, nfix, tolx0,
     &              iw, leniw, rw, lenrw )
               Needx  = inform .gt. 0
            end if

!           Refactorize the basis if it has been modified too much.

            if (LUrequest .eq. 0) then
               if (     iw(LUmod) .ge. kfac-1) then
                  LUrequest = 1
               else if (iw(LUmod) .ge. 20  .and.
     &                 iw(lenL)+iw(lenU) .gt. LUmax) then
                  Bgrwth = iw(lenL) + iw(lenU)
                  Bold   = LUsiz0
                  Bgrwth = Bgrwth/Bold
                  if (Prtlvl10) then
                     write(str, 1000) Bgrwth
                     call snPRNT( 21, str, iw, leniw )
                  end if
                  LUrequest = 2
               end if
               if (LUrequest .gt. 0) typeLU = BT
            end if

!           Update the LU factors.

            if (LUrequest .eq. 0  .and.  ModLU) then
               iw(LUmod) = iw(LUmod) + 1
               call s2Bmod2
     &            ( inform, kp, m, y1, iw, leniw, rw, lenrw )
               if (inform .eq. -1) LUrequest = 5 ! Singular after LU mod
               if (inform .eq.  2) LUrequest = 6 ! Unstable LU mod
               if (inform .eq.  7) LUrequest = 7 ! Insufficient free memory
               if (LUrequest  .gt.  0) then
                  call s2tryLU
     &               ( itn, LUrequest, nS, LUrequest, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )
                  if (.not. LUok) then
                     iExit = 43
                     go to 100
                  end if
               end if
            end if

            if (LUrequest .eq. 0) then
               Checkx = mod(iw(LUitn),kchk) .eq. 0
               if (Checkx  .and.  .not. Needx) then
                  call s5setx
     &               ( inform, Check, itn,
     &                 m, n, nb, nBS, rowError,
     &                 neA, nlocA, locA, indA, Acol,
     &                 kBS, xBS, nrhs0, nrhs, rhs, x, y, y2,
     &                 iw, leniw, rw, lenrw )
                  LUrequest = inform
                  Converged = .false.
               end if
               if (LUrequest .gt. 0) typeLU = BT
            end if
         end if ! not QPdone

         go to 100
!+    end while
      end if
!     ======================end of main loop============================
!
      call s5hs  ( Extern, nb, blQP, buQP, hs, x )

      if (subOptimize .gt. 0) then
         if (frozen .gt. 0) then
!           Relax
         else
            subOptimize = 0
         end if
      end if

      return

 1000 format(' ==> LU file has increased by a factor of', f6.1)
 1010 format(' Biggest dj =', 1p, e11.3, ' (variable', i7, ')',
     &       '    norm rg =',     e11.3, '   norm pi =', e11.3)
 1020 format(   ' Norm rg =', 1p, e11.3, '   norm pi =', e11.3)
 1030 Format(' Itn', i7, ': Infeasible nonelastics.  Num =', i5, 1p,
     &                   '  Sum of Infeasibilities =', e8.1)
 1610 format(' Itn', i7, ': Suboptimize: ', i7, ' new superbasics')
 1620 format(' Itn', i7, ': Suboptimize: ', i7, ' minor iterations')
 1700 format(' Itn', i7, ': Subspace tolerance = ', 1p, e11.3)
 8050 format(' Itn', i7, ': Infeasible ', a)
 8060 format(' Itn', i7, ': Elastic Phase 1 -- making ',
     &                   'non-elastic variables feasible')

      end ! subroutine s5QN

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5QNgetp
     &   ( Hprod, Hprod1,
     &     Feasible, Exactp, itn, QPterm, condZHZ, rgNorm,
     &     maxR, nS, lenR, R, rg, p, w, gp,
     &     neA, nlocA, locA, indA, Acol,
     &     neH, nlocH, locH, indH, Hcol,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      logical
     &     Exactp, Feasible, QPterm
      integer
     &     itn, lencu, lencw, leniu, leniw, lenR, lenru, lenrw, maxR,
     &     neA, neH, nlocA, nlocH, nS, locA(nlocA), locH(nlocH),
     &     indA(neA), indH(neH), iu(leniu), iw(leniw)
      double precision
     &     condZHZ, rgNorm, gp, Acol(neA), Hcol(neH), R(lenR), rg(nS),
     &     p(nS), ru(lenru), rw(lenrw), w(nS)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5QNgetp  computes a search direction  p  for the superbasic
!     variables, using the current reduced gradient  g.
!
!     16 Nov 2001: First version of s5QNgetp.
!     26 Jul 2003: cgItn added.
!     14 Mar 2004: Argument QPterm added for lvlObjE = 2.
!     04 Dec 2004: Current version.
!     ==================================================================
      character
     &     str*120
      external
     &     ddot, s5ZHZv, s5Msolv
      logical
     &     CheckA, Goodb, SYMpre
      integer
     &     cgItn, cgItns, cgItmx, iStop, lprDbg,
     &     lr1, lr2, ls1, ls2, ls3, nout, PreCon, qpMode
      double precision
     &     ddot, Hznorm, rNorm, rtol, shift, tolCG, ynorm
!     ------------------------------------------------------------------
      integer            CG,            QN
      parameter         (CG        = 1, QN     = 2)
      integer            WithR,         WithRt
      parameter         (WithR     = 0, WithRt = 1)
      parameter         (cgItns    = 386) ! Total symmlq iterations
      parameter         (cgItn     = 387) ! symmlq itns for last minor
      double precision   zero,               one
      parameter         (zero      = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      lprDbg    = iw( 85) ! >0    => private debug print
      cgItmx    = iw( 97) ! CG iteration limit
      nout      = iw(151) ! unit # for printed messages
      qpMode    = iw(208) ! Current QP solver   (based on QPslvr)
      PreCon    = iw(209) ! Current precon mode (based on QPslvr)

      lr1       = iw(353) ! r1(maxS) SYMMLQ real work vector
      lr2       = iw(354) ! r1(maxS) SYMMLQ real work vector
      ls1       = iw(355) ! s1(maxS) SYMMLQ real work vector
      ls2       = iw(356) ! s2(maxS) SYMMLQ real work vector
      ls3       = iw(357) ! s3(maxS) SYMMLQ real work vector

      tolCG     = rw( 54) ! cg tolerance

      if (lprDbg .eq. 0)
     &nout      =  0

      CheckA    = .false.
      Exactp    = .true.
      Goodb     = .false.
      SYMpre    = .false.          ! No preconditioning done by SYMMLQ
      iw(cgItn) =  0

      shift     =  -1.0d-6         ! symmlq solves with  (A - shift*I) !!!
      rtol      = tolCG*min( one, rgNorm )

!     Steepest-descent direction.

      call dcopy ( nS, rg, 1, p, 1 )

      if (Feasible  .and.  QPterm) then

!        Preconditioned CG.  Save  w  such that  R'*w = rg.

         if ((qpMode .eq. CG  .and.  PreCon .eq. 1)  .or.
     &        qpMode .eq. QN                            ) then
            call s6Rsol( WithRt, maxR, nS, lenR, R, p )
         end if

         call dcopy ( nS, p, 1, w, 1 )

         if ((qpMode .eq. QN  .and.  nS .gt. maxR)  .or.
     &        qpMode .eq. CG                            ) then
            Exactp = .false.
            call SYMMLQ
     &         ( nS, w, rw(lr1), rw(lr2), rw(ls1), rw(ls2),p,rw(ls3),
     &           s5ZHZv, s5Msolv, CheckA, Goodb, SYMpre, shift,
     &           nout , cgItmx, rtol,
     &           iStop, iw(cgItn), Hznorm, condZHZ, rNorm, ynorm,
     &           Hprod, Hprod1,                        ! Added for SNOPT
     &           neA, nlocA, locA, indA, Acol,         ! Added for SNOPT
     &           neH, nlocH, locH, indH, Hcol,         ! Added for SNOPT
     &           cu, lencu, iu, leniu, ru, lenru,      ! Added for SNOPT
     &           cw, lencw, iw, leniw, rw, lenrw )     ! Added for SNOPT
            if (iStop .eq. -1) then
               call dcopy ( nS, w, 1, p, 1 )
            end if
         end if

         call dcopy ( nS, p, 1, w, 1 )

         if ((qpMode .eq. CG  .and.  PreCon .eq. 1)  .or.
     &        qpMode .eq. QN                            ) then
            call s6Rsol( WithR , maxR, nS, lenR, R, p )
         end if
      end if ! Feasible

!     Change the sign of p.

      call dscal ( nS, (-one), p, 1 )
      gp   = ddot  ( nS, rg, 1, p, 1 )

      if (gp .ge. zero) then
         write(str, 1000) itn, gp, rtol
         call snPRNT( 23, str, iw, leniw )
      end if

      iw(cgItns) = iw(cgItns) + iw(cgItn )

      return

 1000 format(' Itn', i7, ': CG gives  gp = ', 1p, e8.1,
     &                   ' rtol = ', e8.1 )

      end ! subroutine s5QNgetp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5QNitn
     &   ( iExit,
     &     Hprod, Hprod1, HvCalls,
     &     Elastic, Feasible, GotgQP, GotH, GotR,
     &     HitCon, Increase, ModLU, Needf, Needv, NewSB,
     &     itn, lenR, m, mBS, maxR, maxS,
     &     n, nb, nnH0, nnH, nS, ngQP0, ngQP, nDegen,
     &     elastics, LUrequest,
     &     kp, jBq, jSq, jBr, jSr,
     &     jq, lvlTol, nfmove, nUncon,
     &     djq0, djq, minimize, iObj, scaleObj, objQP,
     &     condZHZ, condZmax, featol, pivot,
     &     piNorm, rgNorm, step, tolinc,
     &     nInfE, sInfE, wtInf, pSNorm, pSNrm1, xSNorm, xSNrm1,
     &     neA, nlocA, locA, indA, Acol,
     &     neH, nlocH, locH, indH, Hcol,
     &     eType, eState, feasType, hs, kBS,
     &     bl, bu, blQP, buQP, blBS, buBS, gBS,
     &     gQP, Hdx, pBS, pi, rg, rg2, R,
     &     x, xBS, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      logical
     &     Elastic, Feasible, GotgQP, GotH, GotR, HitCon, Increase,
     &     ModLU, Needf, Needv, NewSB
      integer
     &     elastics, HvCalls, iExit, iObj, itn, lenR, LUrequest, lvlTol,
     &     m, maxR, maxS, mBS, minimize, n, nb, nDegen, neA, neH, nInfE,
     &     nlocA, nlocH, nnH0, nnH, ngQP0, ngQP, nS, kp,
     &     jBq, jBr, jq, jSq, jSr, nfmove, nUncon,
     &     lencu, lencw, leniu, leniw, lenru, lenrw,
     &     eType(nb), eState(nb), feasType(mBS), hs(nb), kBS(mBS),
     &     locA(nlocA), locH(nlocH), indA(neA), indH(neH),
     &     iu(leniu), iw(leniw)
      double precision
     &     condZHZ, condZmax, djq0, djq, featol, objQP, pivot,
     &     pSNorm, pSNrm1, rgNorm, scaleObj, sInfE, step,
     &     tolinc, wtInf, xSNorm, xSNrm1,
     &     Acol(neA), bl(nb), bu(nb), blQP(nb), buQP(nb), blBS(mBS),
     &     buBS(mBS), gBS(mBS), gQP(ngQP0), Hcol(neH), Hdx(nnH0),
     &     pBS(mBS), pi(m), R(lenR), rg(maxS), rg2(maxS), x(nb),
     &     xBS(mBS), y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5QNitn performs a QP step.
!
!     On entry,
!        NewSB = true implies that variable jq just went superbasic.
!                In this case:
!                pBS(1:m)  satisfies B pBS = a(jq).
!                 y1(1:m)  satisfies L  y1 = a(jq).
!
!      iExit       Result
!      -----       ------
!       -1         unbounded
!        0         normal exit
!       >0         Fatal LU error

!     On exit,
!        pBS contains the most recent QP search direction.
!
!     Increase says if the new variable should increase or not.
!            It is used only in linear mode when s5Price is moving free
!            nonbasics toward their closest bound (and djq = zero).
!
!     GotR   says if a useful quasi-Newton preconditioner  R  exists.
!            It will be false if no preconditioning is being used.
!            If  preconditioned CG (QN) is being used, then GotR may
!            be true even if the current itn is infeasible.
!            This allows  R  to be updated in during a temporary loss
!            of feasibility.
!
!     PartHess is true if  R  is not big enough to store a full
!            preconditioner for the reduced Hessian.
!            The null space is then  Z = ( Z1  Z2 ),  and  R  estimates
!            the Hessian only in the subspace  Z1.  A diagonal estimate
!            is used (also in  R) for the subspace  Z2.
!
!     12 Jun 2001: First   version of s5QNitn based on s5Qpit
!                  and MINOS routine m7rgit.
!     14 Jul 2001: Different gp needed for PartHess.
!     02 Aug 2003: snPRNT adopted.
!     19 Jun 2008: eState accessed in elastic mode only.
!     19 Jun 2006: s5Zp added to compute Z*p.
!     11 Dec 2014: Last diagonal element of R is argument of s6Radd.
!     06 May 2015: x(jr) set exactly on its bound.
!     06 May 2015: y1 not recomputed if jq .eq. jSr.
!     06 May 2015: ModLU added to argument list.
!     30 May 2015: Removed bogus check on unbounded QP objective.
!     ==================================================================
      character
     &     str*120
      external
     &     ddot, dnormi, dnrm2
      logical
     &     Checkpi, HitLow, Exactp, Move, OnBound, PartHess, Unbounded,
     &     Uncon
      integer
     &     inform, infpiv, jqeState, jqState, jr, jreState, jrState,
     &     kBSq, kObj, kSq, m1, modR1, mtry, nBS,
     &     ntry, nZ1, nZ2, sqStat
      double precision
     &     normA, bigdx, bound, eps, eps0, eps2, exact,
     &     gp, gpQP, infBnd, normg, pBS1, pHp, pHpQP,
     &     pNorm, piNorm, pSNrm2, rgdel, rgNorm2, sclPiv, signObj,
     &     StepB, stepmx, stepP, t, tolP0, tolP, tolpiv, xSNrm2
      double precision
     &     ddot, dnormi, dnrm2
!     ------------------------------------------------------------------
      parameter         (kObj   = 205) ! xBS(kObj) is the obj. slack
      integer            tight
      parameter         (tight  = 2)
      integer            xBStox
      parameter         (xBStox = 1)
      parameter         (mtry   = 6)
      integer            LOAD
      parameter         (LOAD   = 1)
      integer            WithL,           WithBt
      parameter         (WithL  = 0,      WithBt = 2)
      double precision   zero,            half,          one
      parameter         (zero   = 0.0d+0, half = 0.5d+0, one = 1.0d+0)
      double precision   ten
      parameter         (ten    =10.0d+0)
!     ------------------------------------------------------------------
      eps       = rw(  1) ! unit round-off.  IEEE DP  2.22e-16
      eps0      = rw(  2) ! eps**(4/5)       IEEE DP  3.00e-13
      eps2      = rw(  4) ! eps**(1/2)       IEEE DP  1.49e-08
      tolpiv    = rw( 60) ! excludes small elements of pBS.
      infBnd    = rw( 70) ! definition of an infinite bound
      bigdx     = rw( 72) ! unbounded step.

      iExit     = 0
      inform    = 0
      sqStat    = 0

      Checkpi   = .false.
      Unbounded = .false.
      signObj   = minimize

      m1        = m + 1
      nBS       = m + nS

      if (NewSB) then
!        ---------------------------------------------------------------
!        New superbasic.
!        ---------------------------------------------------------------
         nS        = nS   + 1
         nBS       = nBS  + 1

         kBS (nBS) =      jq
         xBS (nBS) =   x (jq)
         blBS(nBS) = blQP(jq)
         buBS(nBS) = buQP(jq)
         jqState   =   hs(jq)

         feasType(nBS) = 0

         if (Feasible .and. Needf .and. GotgQP .and. jq .le. ngQP) then
            gBS(nBS) = signObj*gQP(jq)
         else
            gBS(nBS) = zero
         end if

!        ===============================================================
!        Set eState(jq) and the elastic parts of blBS and buBS.
!        ===============================================================
         if (Elastic) then

            jqeState = eState(jq) ! Saved value

!           Check if the new superbasic is an elastic variable that
!           wants to move infeasible. If so, set its elastic state and
!           modify the working bounds.

            if (eType(jq) .gt. 0  .and.  eState(jq) .eq. 0) then
               if (Increase) then
                  if (jqState .eq. 1  .or.  jqState .eq. 4) then
                     elastics    =  elastics + 1
                     eState(jq)  =  2
                     blQP(jq)    =  buQP(jq)
                     buQP(jq)    = +infBnd
                     blBS(nBS)   =  blQP(jq)
                     buBS(nBS)   =  buQP(jq)
                     if (Feasible  .and.  Needv) then
                        gBS(nBS) = gBS(nBS) + wtInf
                     end if
                  end if
               else ! the variable is decreasing
                  if (jqState .eq. 0  .or.  jqState .eq. 4) then
                     elastics    =  elastics + 1
                     eState(jq)  =  1
                     buQP(jq)    =  blQP(jq)
                     blQP(jq)    = -infBnd
                     blBS(nBS)   =  blQP(jq)
                     buBS(nBS)   =  buQP(jq)
                     if (Feasible  .and.  Needv) then
                        gBS(nBS) = gBS(nBS) - wtInf
                     end if
                  end if
               end if
            end if
         end if ! Elastic

!        ---------------------------------------------------------------
!        In phase 1, or in phase 2 for an LP, price can select nonbasics
!        floating free between their bounds with zero reduced cost.
!        We have to check that dqj is not zero.
!        ---------------------------------------------------------------
         rg(nS) = djq
         if (.not. Feasible  .or.  (Needf .and.  .not. GotH)) then
            if (hs(jq) .eq. -1) then
               if (Increase) then
                  rg(nS) = -one
               else
                  rg(nS) =  one
               end if
            end if
         end if
         jSq    = jq
         hs(jq) = 2

!        Add a unit column to R at position nS.

         if (GotR) then
            rgNorm2 = dnrm2( nS, rg, 1 )
            rgNorm2 = max( one, sqrt(rgNorm2) )
            call s6Radd( maxR, lenR, nS, rgNorm2, R )
         end if

      end if ! NewSB

!     ------------------------------------------------------------------
!     Compute the search direction for the superbasics.
!     ------------------------------------------------------------------
      PartHess = nS .gt. maxR  .and.  GotR

      if (PartHess) then
         nZ1 = maxR
      else
         nZ1 = nS
      end if

      nZ2 = nS - maxR

!     ------------------------------------------------------------------
!     Store the free components of the search direction in pBS(1:nBS).
!     First, find the search direction pS for the superbasics.
!     The vector v such that R*v = rg is held in y2 for the BFGS update.
!     ------------------------------------------------------------------
  100 call s5QNgetp
     &   ( Hprod, Hprod1,
     &     Feasible, Exactp, itn, Needf, condZHZ, rgNorm,
     &     maxR, nS, lenR, R, rg, pBS(m1), y2, gp,
     &     neA, nlocA, locA, indA, Acol,
     &     neH, nlocH, locH, indH, Hcol,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (gp .ge. zero) then
         GotR      = .false.
         LUrequest = 24
         go to 900
      end if

      pBS1 = pBS(m+nS)

      if (NewSB  .and.  Feasible) then

!        Check that pBS(nBS) is feasible wrt blBS(nBS) and buBS(nBS).
!        A large rgTol may give a pBS(nBS) with the wrong sign.

         if (djq*pBS1 .gt. zero) then

            write(str, 1000) itn
            call snPRNT( 23, str, iw, leniw )

            nS     = nS  - 1
            nBS    = nBS - 1
            hs(jq) = jqState

            if (Elastic) then
!              If a variable went elastic, reset it as not elastic.

               if (eState(jq) .ne. jqeState) then
                  elastics   = elastics - 1
                  blQP(jq)   = bl(jq)
                  buQP(jq)   = bu(jq)
                  estate(jq) = jqeState
               end if
            end if

            jq     =  0
            djq    =  djq0
            jSq    = -jSq
            NewSB  = .false.
            if (lvlTol .eq. tight  .and.  GotR) then
               call s6Rset
     &            ( LOAD, maxR, nS, lenR, R, y, condZHZ )
            else
               lvlTol = tight
            end if
            go to 100
         end if
      end if

!     ------------------------------------------------------------------
!     Find the norm of  pS.  Put the search direction for the basic
!     variables in pBS(1)  ,...,pBS(m).
!     ------------------------------------------------------------------
      xSNrm1 = dnormi( nZ1, xBS(m1), 1 )
      pSNrm1 = dnormi( nZ1, pBS(m1), 1 )

      xSNrm2 = zero
      pSNrm2 = zero

      if (PartHess) then
         xSNrm2 = dnormi( nZ2, xBS(m1+maxR), 1 )
         pSNrm2 = dnormi( nZ2, pBS(m1+maxR), 1 )
      end if

      pSNorm = max ( pSNrm1, pSNrm2 )
      xSNorm = max ( xSNrm1, xSNrm2 )

      call s5Zp
     &   ( iExit, m, mBS, n, nb, nS, eps0, pNorm,
     &     neA, nlocA, locA, indA, Acol,
     &     kBS, pBS, y, iw, leniw, rw, lenrw )
      if (iExit .ne. 0) go to 900

      pHp = zero

      if (Feasible) then
!        ---------------------------------------------------------------
!        Compute y = pBS(scattered) and Hdx(scattered).
!        The vector Hdx is used to update the objective and gradient of
!        the QP.  Form  gpQP  and  pHpQP  for the quadratic term.
!        gp = gpQP - pBS(kObj) + terms from the elastic gradient.
!        ---------------------------------------------------------------
         if (Needf  .and.  (GotgQP  .or.  GotH)) then
            call s2scatter
     &         ( ngQP, nBS, kBS, one, pBS, y )

            if (GotgQP) then
               gpQP = ddot ( ngQP, gQP, 1, y, 1 )
            end if

            if (GotH  ) then
               call Hprod
     &            ( Hprod1, nnH,
     &              neH, nlocH, locH, indH, Hcol,
     &              y, Hdx, sqStat,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )
               HvCalls = HvCalls + 1
               pHpQP   = ddot  ( nnH, y, 1, Hdx, 1 )
               pHp     = signObj*pHpQP
            end if
         end if
      end if ! Feasible

      if (gp .ge. zero) then
         write(str, 1100) itn, gp
         call snPRNT( 23, str, iw, leniw )
         LUrequest = 24
         go to 900
      end if

      NewSB  = .false.

!     ------------------------------------------------------------------
!     Find   the nearest constraint in direction  x + step*pBS (step > 0).
!     Exact  is the step that takes xBS(kp) exactly onto bound.
!            It may be positive or slightly negative.
!            (Not defined if Unbounded.)
!
!     If OnBound  is true, step is a step that reaches a bound exactly.
!     xBS(kp) reaches the value bound.  If we take a constrained step,
!     bound is used to put the new nonbasic variable x(jr) exactly on
!     its bound.
!
!     If Unbounded is true, step = stepmx.
!     ------------------------------------------------------------------
      stepmx = bigdx /pNorm
      sclPiv = one
      tolP0  = tolpiv
      tolP   = tolpiv*pNorm
      ntry   = 0

!+    Repeat
  200    tolP  = tolP /sclPiv
         tolP0 = tolP0/sclPiv
         call s5step
     &      ( nBS, nDegen,
     &        featol, infBnd, stepmx, tolinc, tolP,
     &        feasType, blBS, buBS, xBS, pBS,
     &        HitLow, Move, OnBound, Unbounded,
     &        infpiv, kp, bound, exact, stepB, stepP )

!        Find if the step is constrained or unconstrained.
!        If R has been flagged as singular, we double check by trying
!        to compute the QP minimizer along pBS.  If the minimizer
!        exists,  the singularity tolerance must be too large.

         if (Feasible) then
            Uncon     = stepP*pHp .gt.  (- gp)
            Unbounded =    stepmx .le.  one    .or.
     &                 (Unbounded .and. .not. Uncon)
         else
            Uncon = .false.
         end if

         sclPiv = ten
         ntry   = ntry + 1

!+    until     (   infpiv .eq. 0 .and. (.not. Unbounded .or. Feasible)
!+             .or. ntry   .ge. mtry)
      if (.not.((  infpiv  .eq. 0 .and. (.not. Unbounded .or. Feasible))
     &         .or. ntry   .ge. mtry)) go to 200

      if (Unbounded) then
         iExit = -1
         go to 900
      end if

      HitCon = .not. Uncon

      if (HitCon) then
         nUncon =  0
         step   =  stepB
         pivot  = -pBS(kp)
      else
         nUncon =  nUncon + 1
         step   = (-gp)/pHp
         pivot  =  zero
      end if

!     ------------------------------------------------------------------
!     Note: pHpQP and pHp are defined before and after scaling by
!     signObj.
!     Update the value and gradient of the quadratic obj term.
!     ------------------------------------------------------------------
      if (Feasible  .and.  Needf  .and.  step .gt. zero) then
         if (GotgQP) then
            objQP = objQP + step*gpQP
         end if
         if (GotH  ) then
            objQP = objQP + half*pHpQP*step**2
            call daxpy ( nnH, step, Hdx, 1, gQP, 1 )
         end if
      end if

      if (Feasible  .and.  Move) nfmove = nfmove + 1

!     ------------------------------------------------------------------
!     Update the basic variables xBS.
!     ------------------------------------------------------------------
      call daxpy ( nBS, step, pBS, 1, xBS, 1 )
      call s5BSx ( xBStox, nBS, nb, kBS, x, xBS )

      if (Feasible  .and.  step .gt. zero) then
!        ===============================================================
!        gQP has been changed.
!        Compute the new gBS and reduced gradient rg2.
!        Update the reduced Hessian.
!        ===============================================================
         call dload ( nBS, zero, gBS, 1 )

         if (Needf) then
            if (GotgQP) then
               call s2gather
     &            ( ngQP, nBS, kBS, signObj, gQP, gBS )
            end if
            if (iObj .gt. 0) gBS(iw(kObj)) = signObj*scaleObj
         end if

         if (Elastic  .and.  Needv) then
            call s5eInf
     &         ( nb, nBS, eState, kBS, featol, nInfE, sInfE,
     &           blQP, buQP, x )
            call s5eGrad
     &         ( nb, nBS, wtInf, eState, kBS, gBS )
         end if

         normg = max(dnormi( nBS, gBS, 1 ), one)

         call dcopy ( m, gBS, 1, y, 1 )
         call s5setpi
     &      ( iExit,
     &        m, Checkpi, condZmax, normg, piNorm, y, pi,
     &        iw, leniw, rw, lenrw )
         if (iExit .ne. 0) go to 900

         call s5rg
     &      ( m, nBS, n, nS, eps0,
     &        neA, nlocA, locA, indA, Acol,
     &        gBS, pi, rg2, rgNorm, kBS )

         if (GotR) then
            call s6Rqn
     &         ( inform, Exactp, maxR, nS,
     &           lenR, R, gp, rg, rg2, pBS(m1), y2, step, eps2, eps0 )

!           inform  = 0  if no update was performed,
!                   = 1  if the update was successful,
!                   = 2  if it was nearly singular.
         end if
      else
         call dcopy ( nS, rg, 1, rg2, 1 )
      end if

      if (Uncon) then
!        ===============================================================
!        The step is unconstrained.
!        ===============================================================
         call dcopy ( nS, rg2, 1, rg, 1 )
         ModLU = .false.        ! No LU update

      else                      ! hit constraint
!        ===============================================================
!        There is a blocking variable.
!        It could be a fixed variable, whose new state must be 4.
!        ===============================================================
         jr    =   kBS(kp)

         if (OnBound) then
            x(jr) = bound
         else if (HitLow) then
            x(jr) = min( x(jr), blQP(jr) )
         else
            x(jr) = max( x(jr), buQP(jr) )
         end if

         if (Elastic) then
            jreState = eState(jr)
         else
            jreState = 0
         end if

         if (jreState .eq. 0) then
            if (blBS(kp) .eq. buBS(kp)) then
               jrstate = 4
            else if (HitLow) then
               jrstate = 0
            else
               jrstate = 1
            end if
         else

!           Elastic x(jr) hits its bound.
!           Reset the true upper and lower bounds.

            eState(jr) = 0
            elastics   = elastics - 1
            blQP(jr)   = bl(jr)
            buQP(jr)   = bu(jr)

            if (jreState .eq. 1) then
               if (blQP(jr)   .eq. buQP(jr)) then
                  jrstate =  4
               else if (OnBound) then
                  jrstate =  0
               else if (x(jr) .lt. buQP(jr)) then
                  jrstate = -1
               else
                  jrstate =  1
               end if

            else !   js .eq. 2
               if (blQP(jr) .eq. buQP(jr)) then
                  jrstate =  4
               else if (OnBound) then
                  jrstate =  1
               else if (x(jr) .gt. blQP(jr)) then
                  jrstate = -1
               else
                  jrstate =  0
               end if
            end if
         end if

         if (kp .le. m) then
!           ============================================================
!           A variable in B hit a bound.
!           Find column kSq = kBSq-m  of S to replace column kp of B.
!           If nS = 1, it must be the entering SB column.
!           ============================================================
!           if (nS .eq. 1) then
!              kBSq  = nBS
!              pivot = pivot/pBS1
!           else
               call dload ( m, zero, y2, 1 )
               y2(kp) = one     ! Set      y2 = ep
                                ! Solve  B'yB = ep
               call s2Bsol
     &            ( inform, WithBt, m, y2, y, iw, leniw, rw, lenrw )
               call s5chzq
     &            ( m, mBS, n, nb, nS, kBSq, pivot, tolP0,
     &              neA, nlocA, locA, indA, Acol,
     &              kBS, blQP, buQP, xBS, y, iw, leniw, rw, lenrw )
               if (kBSq .le. 0) then
                  write(str, 9999) itn
                  call snPRNT( 23, str, iw, leniw )
                  kBSq   = nBS
               end if
!           end if

            kSq        = kBSq - m

            hs(jr)     = jrState
            jBr        = jr                     ! Outgoing basic
            jSr        =  kBS(kBSq)             ! Outgoing superbasic
            kBS (kBSq) = jBr
            jBq        = jSr                    ! Incoming basic
            kBS (kp)   = jSr
            blBS(kp)   = blBS(kBSq)
            buBS(kp)   = buBS(kBSq)
            xBS (kp)   = xBS (kBSq)
            gBS (kp)   = gBS (KBSq)
            hs(jBq)    = 3

!           Finish computing yS = (y(m+1), ..., y(m+nS)).

            y(kBSq) = - (one + pivot)
            call dscal ( nS, (one/pivot), y(m1), 1 )

            if (GotR  .and.  kSq .le. maxR) then
               call s6Rswap( maxR, nZ1, lenR, R, y2, y(m1), kSq, eps0 )
            end if

            if (Feasible) then

!              Modify  pi  using  y  where  B' y = e(kp).

               t      = rg2(kSq) / pivot
               call daxpy ( m, t, y, 1, pi, 1 )
               piNorm = max( dnormi( m, pi, 1 ), one )
            end if

!           ------------------------------------------------------------
!           Get a new  y1, used to modify L and U.  If the outgoing
!           superbasic just came in, we already have it.
!           ------------------------------------------------------------
            if (jSr .ne. jq) then
               call s2unpack
     &            ( jBq, m, n, neA, normA, nlocA, locA, indA, Acol, y1 )
               call s2Bsol
     &            ( iExit, WithL, m, y1, y, iw, leniw, rw, lenrw )
               if (iExit .ne. 0) return
            end if
            ModLU  = .true.      ! Update the LU on exit

         else
!           ============================================================
!           A variable in S hit a bound.
!           ============================================================
            hs(jr) = jrState
            jSr    = jr
            kBSq   = kp
            kSq    = kBSq - m
            ModLU  = .false.     ! No LU update
         end if

         if (Feasible) then
            call s5rg
     &         ( m, nBS, n, nS, eps0,
     &           neA, nlocA, locA, indA, Acol,
     &           gBS, pi, rg, rgNorm, kBS )
         end if

!        ===============================================================
!        If necessary, swap the largest reduced-gradient in  Z2  to the
!        front of  Z2,  so it will end up at the end of  Z1.
!        ===============================================================
         rgdel  = abs( rg(kSq) )
         if (PartHess  .and.  kSq .le. maxR) then
            call s5Sswap
     &         ( m, maxR, lenR, nS, nBS,
     &           kBS, blBS, buBS, gBS, R, rg, xBS )
         end if

!        ---------------------------------------------------------------
!        Delete the kSq-th superbasic, updating R if it exists and
!        adjusting all arrays in BS order.
!        ---------------------------------------------------------------
         if (GotR  .and.  kSq .lt. nS) then
            call s6Rdel
     &         ( kSq, maxR, nS, lenR, R, eps )
         end if

         call s5Sdel
     &      ( kSq, m, nS, nBS, kBS, blBS, buBS, gBS, rg, xBS )

         nS     = nS  - 1
         nBS    = m   + nS

         if (rgNorm .le. rgdel) then
            rgNorm = zero
            if (nS .gt. 0) then
               rgNorm = dnormi( nS, rg, 1 )
            end if
         end if
      end if ! if unCon

  900 return

 1000 format(' Itn', i7, ': Bad direction after adding a superbasic.')
 1100 format(' Itn', i7, ': CG gives  gp = ', 1p, e8.1)
 9999 format(' Itn', i7, ': Chzq failed in s5QNitn!!')

      end ! subroutine s5QNitn

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Sswap
     &   ( m, maxR, lenR, nS, nBS,
     &     kBS, blBS, buBS, gBS, R, rg, xBS )

      implicit
     &     none
      integer
     &     m, maxR, lenR, nS, nBS, kBS(nBS)
      double precision
     &     blBS(nBS), buBS(nBS), gBS(nBS), R(lenR), rg(nS), xBS(nBS)

!     ==================================================================
!     s5Sswap  (superbasic swap)  finds the largest reduced gradient in
!     the range  rg(maxR+1), ..., rg(nS)  and swaps it into position
!     maxR + 1  (so we know where it is).
!
!     16 Jun 2001: First version based on MINOS routine m6swap.
!     16 Jun 2001: Current version of s5Sswap.
!     ==================================================================
      external
     &     idamax
      integer
     &     idamax, j, j1, k, k1, k2, lastR, ldiag1, ldiag2, nz2
      double precision
     &     bl1, bu1, gBS1, rdiag1, rg1, xBS1
!     ------------------------------------------------------------------
      k1  = maxR + 1

      if (nS .gt. k1) then
         nz2    = nS - maxR
         k2     = maxR + idamax( nz2, rg(k1), 1 )
         if (k2 .gt. k1) then
            j      = m + k1
            k      = m + k2
            lastR  = maxR*k1/2
            ldiag1 = lastR + 1
            ldiag2 = lastR + (k2 - maxR)

            rdiag1    = R(ldiag1)
            rg1       = rg(k1)
            j1        = kBS(j)
            bl1       = blBS(j)
            bu1       = buBS(j)
            gBS1      = gBS(j)
            xBS1      = xBS(j)

            R(ldiag1) = R(ldiag2)
            rg(k1)    = rg(k2)
            kBS(j)    = kBS(k)
            blBS(j)   = blBS(k)
            buBS(j)   = buBS(k)
            gBS(j)    = gBS(k)
            xBS(j)    = xBS(k)

            R(ldiag2) = rdiag1
            rg(k2)    = rg1
            kBS(k)    = j1
            blBS(k)   = bl1
            buBS(k)   = bu1
            gBS(k)    = gBS1
            xBS(k)    = xBS1
         end if
      end if

      end ! subroutine s5Sswap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5ZHZv
     &   ( nS, x, ZHZx,
     &     Hprod, Hprod1,
     &     neA, nlocA, locA, indA, Acol,
     &     neH, nlocH, locH, indH, Hcol,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      integer
     &     lencu, lencw, leniu, leniw, lenru, lenrw, neA, neH,
     &     nlocA, nlocH,nS, locA(nlocA), locH(nlocH), indA(neA),
     &     indH(neH), iu(leniu), iw(leniw)
      double precision
     &     Acol(neA), Hcol(neH), x(nS), ZHZx(nS), ru(lenru), rw(lenrw)
      character
     &    cu(lencu)*8,  cw(lencw)*8

!     ==================================================================
!     s5ZHZv  is called by SYMMLQ to compute the vector
!       ZHZx = Z'HZ x  or  R^(-T) Z'HZ R^(-1) x.
!     (Preconditioning with R'R is done explicitly,
!     not via SYMMLQ's Msolve.)
!
!     16 Nov 2001: First version of s5ZHZv.
!     17 Aug 2004: Current version.
!     ==================================================================
      integer
     &     HvCalls, lenR, lkBS, lR, ly, ly3, mBS, m, maxR, maxS,
     &     n, nb, nnH
!     ------------------------------------------------------------------
      parameter         (HvCalls    = 188) ! number of Hx products
!     ------------------------------------------------------------------
      lenR      = iw( 28) ! R(lenR) is the reduced Hessian factor
      n         = iw( 15) ! copy of the number of columns
      m         = iw( 16) ! copy of the number of rows
      nnH       = iw( 24) !    max( nnObj, nnJac )
      maxR      = iw( 52) ! max columns of R
      maxS      = iw( 53) ! max # of superbasics
      lkBS      = iw(292) ! kBS(mBS)    = ( B  S ) list

      lR        = iw(295) ! R(lenR)     = factor of Z'HZ
      ly        = iw(311) ! y (nb)      =  real work vector
      ly3       = iw(314) ! y3(nb)      =  real work vector

      nb    = n + m
      mBS   = m + maxS

      call s5ZHZv1
     &   ( Hprod, Hprod1, iw(HvCalls),
     &     nnH, m, mBS, n, nb, nS, maxR,
     &     neA, nlocA, locA, indA, Acol,
     &     neH, nlocH, locH, indH, Hcol,
     &     iw(lkBS), x, ZHZx, rw(ly), rw(ly3), lenR, rw(lR),
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine s5ZHZv

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5ZHZv1
     &   ( Hprod, Hprod1, HvCalls,
     &     nnH, m, mBS, n, nb, nS, maxR,
     &     neA, nlocA, locA, indA, Acol,
     &     neH, nlocH, locH, indH, Hcol,
     &     kBS, x, ZHZx, y, y1, lenR, R,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      integer
     &     HvCalls, lencu, lencw, leniu, leniw, lenru, lenrw, lenR,
     &     maxR, mBS, n, nb, nnH, nS, neA, neH, nlocA, nlocH,
     &     indA(neA), indH(neH), locA(nlocA), locH(nlocH), kBS(mBS),
     &     iu(leniu), iw(leniw)
      double precision
     &     Acol(neA), Hcol(neH), R(lenR), y(nb), y1(nb), x(nS),
     &     ZHZx(nS), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5ZHZv1  does the work for s5ZHZv.
!     Computes the vector ZHZx = Z'HZ * x  or  R^(-T) Z'HZ R^(-1) x.
!
!     16 Nov 2001: First version of s5ZHZv1.
!     17 Aug 2004: Current version.
!     ==================================================================
      integer
     &     inform, m, m1, minimize, nBS, PreCon, sqStat
      double precision
     &     eps0, signObj
!     ------------------------------------------------------------------
      integer            Normal,     Transp
      parameter         (Normal = 0, Transp = 1)
      integer            WithR,      WithRt
      parameter         (WithR  = 0, WithRt = 1)
      integer            WithB,      WithBt
      parameter         (WithB  = 1, WithBt = 2)

      double precision   zero,           one
      parameter         (zero = 0.0d+0,  one = 1.0d+0)
!     ------------------------------------------------------------------
      eps0      = rw(  2) ! eps**(4/5)       IEEE DP  3.00e-13
      minimize  = iw(199) ! (-1)(+1)    => (max)(min)
      PreCon    = iw(209) ! Current precon mode (based on QPslvr)

      nBS       = m + nS
      m1        = m + 1

      signObj   = minimize
      sqStat    = 0

!     ------------------------------------------------------------------
!     Set y = Z*x.
!     ------------------------------------------------------------------
      call dcopy ( nS, x, 1, y(m1), 1 )

      if (PreCon .eq. 1) then
         call s6Rsol( WithR, maxR, nS, lenR, R, y(m1) )
      end if

!     Compute  y1 = - S*yS.

      call s2Bprod
     &   ( Normal, eps0, n, nS, kBS(m1),
     &     neA, nlocA, locA, indA, Acol,
     &     (-one), y(m1), nS, zero, y1, m )

!     Solve  B*y = y1.

      call s2Bsol
     &   ( inform, WithB, m, y1, y, iw, leniw, rw, lenrw  )

!     ------------------------------------------------------------------
!     Set  y = H*y = H*Z*x.
!     ------------------------------------------------------------------
      call s2scatter
     &   ( nnH, nBS, kBS, one, y, y1 )
!     call s8Hx
!    &   ( HvCalls, nnH, y1, y, cw, lencw, iw, leniw, rw, lenrw )

      call Hprod
     &   ( Hprod1, nnH,
     &     neH, nlocH, locH, indH, Hcol,
     &     y1, y, sqStat,       ! y1 = dx,  y = Hdx
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )
      HvCalls = HvCalls + 1

!     ------------------------------------------------------------------
!     Compute  ZHZx = Z'y.
!     ------------------------------------------------------------------
!     Gather y1 = yBS and solve  B'y = y1B.

      call s2gather
     &   ( nnH, nBS, kBS, signObj, y, y1 )
      call s2Bsol
     &   ( inform, WithBt, m, y1, y, iw, leniw, rw, lenrw )

!     Set ZHZx = y1S - S'y.

      call dcopy
     &   ( nS, y1(m1), 1, ZHZx, 1 )
      call s2Bprod
     &   ( Transp, eps0, n, nS, kBS(m1),
     &     neA, nlocA, locA, indA, Acol,
     &     (-one), y, m, one, ZHZx, nS )

      if (PreCon .eq. 1) then
         call s6Rsol
     &      ( WithRt, maxR, nS, lenR, R, ZHZx )
      end if

      end ! subroutine s5Hz1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Zswap
     &   ( Stationary, m, maxR, maxS, nBS, lenR, nS,
     &     tolrg, kBS, blBS, buBS, gBS, R, rg, xBS, rw, lenrw )

      implicit
     &     none
      logical
     &     Stationary
      integer
     &     lenrw, lenR, m, maxR, maxS, nBS, nS,
     &     kBS(nBS)
      double precision
     &     tolrg, blBS(nBS), buBS(nBS),
     &     gBS(nBS), R(lenR), rg(maxS), xBS(nBS), rw(lenrw)

!     ==================================================================
!     s5Zswap  (subspace convergence)  decides whether or not
!     optimization should continue on the current subspace.
!     On exit,  Stationary = .false.  means it should,
!               Stationary = .true    means it should not.
!
!     R  is a partial Hessian for the
!     first  maxR  superbasic variables.  The superbasics are then
!     in two sets  Z1  (containing   nZ1 = maxR       variables)
!             and  Z2  (containing   nZ2 = nS - maxR  variables).
!
!     The null-space matrix is similarly partitioned as  Z = ( Z1  Z2 ).
!
!     The normal convergence test is first applied to  Z1.  If it looks
!     like we are near enough to an optimum in that restricted
!     subspace, we find the largest reduced gradient in  Z2  (index k2).
!     If this is less than  tolrg  we exit with  nxtphs = 3  (to ask for
!     price).  Otherwise we make room in  Z1  for the corresponding
!     variable by the moving the superbasic with the smallest reduced
!     gradient in  Z1  (index  k1)  to the end of  Z2.
!
!     16 Jun 2001: First version based on MINOS routine m7sscv.
!     19 Jul 2001: Current version of s5Sswap.
!     ==================================================================
      integer
     &     jZ1, k, k1, k2, lastR, ldiag1, lR
      double precision
     &     blBS1, buBS1, eps, gBS1, Rdiag1, rg1,
     &     rgmin1, rgNrm2, xBS1
!     ------------------------------------------------------------------
      eps       = rw(  1) ! unit round-off.  IEEE DP  2.22e-16

!     Swap the largest reduced gradient in  Z2  to the front of  Z2
!     and see if it is significantly large.

      call s5Sswap
     &   ( m, maxR, lenR, nS, nBS,
     &     kBS, blBS, buBS, gBS, R, rg, xBS )

      k2     = maxR + 1
      rgNrm2 = abs( rg(k2) )

      if (rgNrm2 .gt. tolrg) then

!        Find the smallest component of  Z1'g.

         rgmin1 = abs( rg(1) )
         k1     = 1
         do  k  = 1, maxR
            if (rgmin1 .ge. abs( rg(k) )) then
               rgmin1 = abs( rg(k) )
               k1     = k
            end if
         end do

         if (rgmin1 .lt. rgNrm2) then

!           Save the relevant values.

            Stationary = .false.
            lastR  = maxR*(maxR + 1)/2
            ldiag1 = (k1 - 1)*maxR + (3 - k1)*k1/2  ! Magic formula!
            Rdiag1 = R(ldiag1)
            rg1    = rg(k1)
            k      = m + k1
            jZ1    = kBS(k)
            blBS1  = blBS(k)
            buBS1  = buBS(k)
            gBS1   = gBS(k)
            xBS1   = xBS(k)

!           Delete the k1-th variable from  Z1,  and shift the remaining
!           superbasics in  Z1  and  Z2  one place to the left.

            call s6rdel
     &         ( k1, maxR, nS, lenR, R, eps )
            call s5Sdel
     &         ( k1, m, nS, nBS, kBS, blBS, buBS, gBS, rg, xBS )

!           Put the old k1-th superbasic in at the very end.

            lR      = lastR + (nS - maxR)
            R(lR)   = Rdiag1
            rg(nS)  = rg1
            kBS(nBS)  = jZ1
            blBS(nBS) = blBS1
            buBS(nBS) = buBS1
            gBS(nBS)  = gBS1
            xBS(nBS)  = xBS1
         end if
      end if

      end ! subroutine s5Zswap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine SYMMLQ
     &   ( n, b, r1, r2, v, w, x, y,
     &     Aprod, Msolve, CheckA, Goodb, Precon, shift,
     &     nout , itnlim, rtol,
     &     iStop, itn, normA, aCond, rnorm, ynorm,
     &     Hprod, Hprod1,                              ! Added for SNOPT
     &     neA, nlocA, locA, indA, Acol,               ! Added for SNOPT
     &     neH, nlocH, locH, indH, Hcol,               ! Added for SNOPT
     &     cu, lencu, iu, leniu, ru, lenru,            ! Added for SNOPT
     &     cw, lencw, iw, leniw, rw, lenrw )           ! Added for SNOPT

      implicit           none
      external           Aprod, Msolve
      integer            n, nout, itnlim, iStop, itn
      logical            CheckA, Goodb, Precon
      double precision   shift, rtol, normA, aCond, rnorm, ynorm,
     &                   b(n), r1(n), r2(n), v(n), w(n), x(n), y(n)

      external                                         ! Added for SNOPT
     &     Hprod, Hprod1                               ! Added for SNOPT
      integer                                          ! Added for SNOPT
     &     neA, neH, nlocA, nlocH,                     ! Added for SNOPT
     &     lencu, lencw, leniu, leniw, lenru, lenrw,   ! Added for SNOPT
     &     locA(nlocA), locH(nlocH),                   ! Added for SNOPT
     &     indA(neA), indH(neH),                       ! Added for SNOPT
     &     iu(leniu), iw(leniw)                        ! Added for SNOPT
      double precision                                 ! Added for SNOPT
     &     Acol(neA), Hcol(neH), ru(lenru), rw(lenrw)  ! Added for SNOPT
      character                                        ! Added for SNOPT
     &     cu(lencu)*8, cw(lencw)*8                    ! Added for SNOPT
!     ------------------------------------------------------------------
!
!     SYMMLQ  is designed to solve the system of linear equations
!
!                Ax = b
!
!     where A is an n by n symmetric matrix and b is a given vector.
!     The matrix A is not required to be positive definite.
!     (If A is known to be definite, the method of conjugate gradients
!     might be preferred, since it will require about the same number of
!     iterations as SYMMLQ but slightly less work per iteration.)
!
!
!     The matrix A is intended to be large and sparse.  It is accessed
!     by means of a subroutine call of the form
!
!                call Aprod ( n, x, y )
!
!     which must return the product y = Ax for any given vector x.
!
!
!     More generally, SYMMLQ is designed to solve the system
!
!                (A - shift*I) x = b
!
!     where  shift  is a specified scalar value.  If  shift  and  b
!     are suitably chosen, the computed vector x may approximate an
!     (unnormalized) eigenvector of A, as in the methods of
!     inverse iteration and/or Rayleigh-quotient iteration.
!     Again, the matrix (A - shift*I) need not be positive definite.
!     The work per iteration is very slightly less if  shift = 0.
!
!
!     A further option is that of preconditioning, which may reduce
!     the number of iterations required.  If M = C C' is a positive
!     definite matrix that is known to approximate  (A - shift*I)
!     in some sense, and if systems of the form  My = x  can be
!     solved efficiently, the parameters precon and Msolve may be
!     used (see below).  When  precon = .true., SYMMLQ will
!     implicitly solve the system of equations
!
!             P (A - shift*I) P' xbar  =  P b,
!
!     i.e.                  Abar xbar  =  bbar
!     where                         P  =  C**(-1),
!                                Abar  =  P (A - shift*I) P',
!                                bbar  =  P b,
!
!     and return the solution       x  =  P' xbar.
!     The associated residual is rbar  =  bbar - Abar xbar
!                                      =  P (b - (A - shift*I)x)
!                                      =  P r.
!
!     In the discussion below, eps refers to the machine precision.
!     eps is computed by SYMMLQ.  A typical value is eps = 2.22e-16
!     for IBM mainframes and IEEE double-precision arithmetic.
!
!     Parameters
!     ----------
!
!     n       input      The dimension of the matrix A.
!
!     b(n)    input      The rhs vector b.
!
!     r1(n)   workspace
!     r2(n)   workspace
!     v(n)    workspace
!     w(n)    workspace
!
!     x(n)    output     Returns the computed solution  x.
!
!     y(n)    workspace
!
!     Aprod   external   A subroutine defining the matrix A.
!                        For a given vector x, the statement
!
!                              call Aprod ( n, x, y )
!
!                        must return the product y = Ax
!                        without altering the vector x.
!
!     Msolve  external   An optional subroutine defining a
!                        preconditioning matrix M, which should
!                        approximate (A - shift*I) in some sense.
!                        M must be positive definite.
!                        For a given vector x, the statement
!
!                              call Msolve( n, x, y )
!
!                        must solve the linear system My = x
!                        without altering the vector x.
!
!                        In general, M should be chosen so that Abar has
!                        clustered eigenvalues.  For example,
!                        if A is positive definite, Abar would ideally
!                        be close to a multiple of I.
!                        If A or A - shift*I is indefinite, Abar might
!                        be close to a multiple of diag( I  -I ).
!
!                        NOTE.  The program calling SYMMLQ must declare
!                        Aprod and Msolve to be external.
!
!     CheckA  input      If CheckA = .true., an extra call of Aprod will
!                        be used to check if A is symmetric.  Also,
!                        if precon = .true., an extra call of Msolve
!                        will be used to check if M is symmetric.
!
!     Goodb   input      Usually, Goodb should be .false.
!                        If x is expected to contain a large multiple of
!                        b (as in Rayleigh-quotient iteration),
!                        better precision may result if Goodb = .true.
!                        See Lewis (1977) below.
!                        When Goodb = .true., an extra call to Msolve
!                        is required.
!
!     Precon  input      If Precon = .true., preconditioning will
!                        be invoked.  Otherwise, subroutine Msolve
!                        will not be referenced; in this case the
!                        actual parameter corresponding to Msolve may
!                        be the same as that corresponding to Aprod.
!
!     shift   input      Should be zero if the system Ax = b is to be
!                        solved.  Otherwise, it could be an
!                        approximation to an eigenvalue of A, such as
!                        the Rayleigh quotient b'Ab / (b'b)
!                        corresponding to the vector b.
!                        If b is sufficiently like an eigenvector
!                        corresponding to an eigenvalue near shift,
!                        then the computed x may have very large
!                        components.  When normalized, x may be
!                        closer to an eigenvector than b.
!
!     nout    input      A file number.
!                        If nout .gt. 0, a summary of the iterations
!                        will be printed on unit nout.
!
!     itnlim  input      An upper limit on the number of iterations.
!
!     rtol    input      A user-specified tolerance.  SYMMLQ terminates
!                        if it appears that norm(rbar) is smaller than
!                              rtol * norm(Abar) * norm(xbar),
!                        where rbar is the transformed residual vector,
!                              rbar = bbar - Abar xbar.
!
!                        If shift = 0 and Precon = .false., SYMMLQ
!                        terminates if norm(b - A*x) is smaller than
!                              rtol * norm(A) * norm(x).
!
!     iStop   output     An integer giving the reason for termination...
!
!              -1        beta2 = 0 in the Lanczos iteration; i.e. the
!                        second Lanczos vector is zero.  This means the
!                        rhs is very special.
!                        If there is no preconditioner, b is an
!                        eigenvector of A.
!                        Otherwise (if Precon is true), let My = b.
!                        If shift is zero, y is a solution of the
!                        generalized eigenvalue problem Ay = lambda My,
!                        with lambda = alpha1 from the Lanczos vectors.
!
!                        In general, (A - shift*I)x = b
!                        has the solution         x = (1/alpha1) y
!                        where My = b.
!
!               0        b = 0, so the exact solution is x = 0.
!                        No iterations were performed.
!
!               1        Norm(rbar) appears to be less than
!                        the value  rtol * norm(Abar) * norm(xbar).
!                        The solution in  x  should be acceptable.
!
!               2        Norm(rbar) appears to be less than
!                        the value  eps * norm(Abar) * norm(xbar).
!                        This means that the residual is as small as
!                        seems reasonable on this machine.
!
!               3        Norm(Abar) * norm(xbar) exceeds norm(b)/eps,
!                        which should indicate that x has essentially
!                        converged to an eigenvector of A
!                        corresponding to the eigenvalue shift.
!
!               4        aCond (see below) has exceeded 0.1/eps, so
!                        the matrix Abar must be very ill-conditioned.
!                        x may not contain an acceptable solution.
!
!               5        The iteration limit was reached before any of
!                        the previous criteria were satisfied.
!
!               6        The matrix defined by Aprod does not appear
!                        to be symmetric.
!                        For certain vectors y = Av and r = Ay, the
!                        products y'y and r'v differ significantly.
!
!               7        The matrix defined by Msolve does not appear
!                        to be symmetric.
!                        For vectors satisfying My = v and Mr = y, the
!                        products y'y and r'v differ significantly.
!
!               8        An inner product of the form  x' M**(-1) x
!                        was not positive, so the preconditioning matrix
!                        M does not appear to be positive definite.
!
!                        If iStop .ge. 5, the final x may not be an
!                        acceptable solution.
!
!     itn     output     The number of iterations performed.
!
!     normA   output     An estimate of the norm of the matrix operator
!                        Abar = P (A - shift*I) P',   where P = C**(-1).
!
!     aCond   output     An estimate of the condition of Abar above.
!                        This will usually be a substantial
!                        under-estimate of the true condition.
!
!     rnorm   output     An estimate of the norm of the final
!                        transformed residual vector,
!                           P (b  -  (A - shift*I) x).
!
!     ynorm   output     An estimate of the norm of xbar.
!                        This is sqrt( x'Mx ).  If Precon is false,
!                        ynorm is an estimate of norm(x).
!
!
!
!     To change precision
!     -------------------
!
!     Alter the words
!            double precision,
!            daxpy, dcopy, ddot, dnrm2
!     to their single or double equivalents.
!     ------------------------------------------------------------------
!
!
!     This routine is an implementation of the algorithm described in
!     the following references:
!
!     C.C. Paige and M.A. Saunders,  Solution of Sparse Indefinite
!          Systems of Linear Equations,
!          SIAM J. Numer. Anal. 12, 4, September 1975, pp. 617-629.
!
!     J.G. Lewis,  Algorithms for Sparse Matrix Eigenvalue Problems,
!          Report STAN-CS-77-595, Computer Science Department,
!          Stanford University, Stanford, California, March 1977.
!
!     Applications of SYMMLQ and the theory of preconditioning
!     are described in the following references:
!
!     D.B. Szyld and O.B. Widlund,  Applications of Conjugate Gradient
!          Type Methods to Eigenvalue Calculations,
!          in R. Vichnevetsky and R.S. Steplman (editors),
!          Advances in Computer Methods for Partial Differential
!          Equations -- III, IMACS, 1979, 167-173.
!
!     D.B. Szyld,  A Two-level Iterative Method for Large Sparse
!          Generalized Eigenvalue Calculations,
!          Ph. D. dissertation, Department of Mathematics,
!          New York University, New York, October 1983.
!
!     P.E. Gill, W. Murray, D.B. Ponceleon and M.A. Saunders,
!          Preconditioners for indefinite systems arising in
!          optimization, SIMAX 13, 1, 292--311, January 1992.
!          (SIAM J. on Matrix Analysis and Applications)
!     ------------------------------------------------------------------
!
!
!     SYMMLQ development:
!            1972: First version.
!            1975: John Lewis recommended modifications to help with
!                  inverse iteration:
!                  1. Reorthogonalize v1 and v2.
!                  2. Regard the solution as x = x1  +  bstep * b,
!                     with x1 and bstep accumulated separately
!                     and bstep * b added at the end.
!                     (In inverse iteration, b might be close to the
!                     required x already, so x1 may be a lot smaller
!                     than the multiple of b.)
!            1978: Daniel Szyld and Olof Widlund implemented the first
!                  form of preconditioning.
!                  This required both a solve and a multiply with M.
!            1979: Implemented present method for preconditioning.
!                  This requires only a solve with M.
!            1984: Sven Hammarling noted corrections to tnorm and x1lq.
!                  SYMMLQ added to NAG Fortran Library.
!     15 Sep 1985: Final F66 version.  SYMMLQ sent to "misc" in netlib.
!     16 Feb 1989: First F77 version.
!
!     22 Feb 1989: Hans Mittelmann observed beta2 = 0 (hence failure)
!                  if Abar = const*I.  iStop = -1 added for this case.
!
!     01 Mar 1989: Hans Mittelmann observed premature termination on
!                  ( 1  1  1 )     (   )                   ( 1  1    )
!                  ( 1  1    ) x = ( 1 ),  for which  T3 = ( 1  1  1 ).
!                  ( 1     1 )     (   )                   (    1  1 )
!                  T2 is exactly singular, so estimating cond(A) from
!                  the diagonals of Lbar is unsafe.  We now use
!                  L       or  Lbar         depending on whether
!                  lqnorm  or  cgnorm       is least.
!
!     03 Mar 1989: eps computed internally instead of coming in as a
!                  parameter.
!     07 Jun 1989: ncheck added as a parameter to say if A and M
!                  should be checked for symmetry.
!                  Later changed to CheckA (see below).
!     20 Nov 1990: Goodb added as a parameter to make Lewis's changes
!                  an option.  Usually b is NOT much like x.  Setting
!                  Goodb = .false. saves a call to Msolve at the end.
!     20 Nov 1990: Residual not computed exactly at end, to save time
!                  when only one or two iterations are required
!                  (e.g. if the preconditioner is very good).
!                  Beware, if precon is true, rnorm estimates the
!                  residual of the preconditioned system, not Ax = b.
!     04 Sep 1991: Parameter list changed and reordered.
!                  integer ncheck is now logical CheckA.
!     22 Jul 1992: Example from Lothar Reichel and Daniela Calvetti
!                  showed that beta2 = 0 (iStop = -1) means that
!                  b is an eigenvector when M = I.
!                  More complicated if there is a preconditioner;
!                  not clear yet how to describe it.
!     20 Oct 1999: Bug.  alfa1 = 0 caused normA = 0, divide by zero.
!                  Need to estimate normA from column of Tk.
!     18 Nov 2001: bnorm printed in first line, since rnorm at itn = 0
!                  isn't bnorm as we might expect -- it's cgnorm after
!                  the proverbial step to the CG point!
!     02 Aug 2003: eps grabbed from rw (special to SNOPT).
!     23 Dec 2003: Normal backward error exit (iStop = 1) terminates
!                  too early when rtol = 0.01 say (rather large).
!                  Keep it with   rtol = eps  (iStop = 2)
!                  but replace normal test with ||r|| < rtol*||b||
!                  like most other CG solvers.  This fits better with
!                  the inexact Newton context.  Note that it is correct
!                  only with Precon = .false.
!
!     08 Apr 2004: Always move to the CG point.
!
!     Michael A. Saunders                   na.msaunders@na-net.ornl.gov
!     Systems Optimization Laboratory       saunders@stanford.edu
!     Dept of Management Science and Engineering
!     Stanford University
!     Stanford, CA 94305-4026                             (650) 723-1875
!     ------------------------------------------------------------------
!
!
!     Subroutines and functions
!
!     USER       Aprod, Msolve
!     BLAS       daxpy, dcopy, ddot , dnrm2
!
!
!     Intrinsics and local variables

      intrinsic          abs, max, min, mod, sqrt
      double precision   ddot, dnrm2
      double precision   alfa, b1, beta, beta1, bnorm, bstep, cs,
     &                   cgnorm, dbar, delta, denom, diag,
     &                   eps, epsa, epsln, epsr, epsx,
     &                   gamma, gbar, gmax, gmin,
     &                   lqnorm, oldb, qrnorm, rhs1, rhs2,
     &                   s, sn, snprod, t, tnorm,
     &                   x1cg, ynorm2, zbar, z
      integer            i

      double precision   zero         ,  one
      parameter        ( zero = 0.0d+0,  one = 1.0d+0)

      character          enter*16, exit*16, msg(-1:9)*52

      data               enter /' Enter SYMMLQ.  '/,
     &                   exit  /' Exit  SYMMLQ.  '/

      data               msg
     & / 'beta2 = 0.  If M = I, b and x are eigenvectors of A',
     &   'beta1 = 0.  The exact solution is  x = 0',
     &   'Requested accuracy achieved, as determined by rtol',
     &   'Reasonable accuracy achieved, given eps',
     &   'x has converged to an eigenvector',
     &   'aCond has exceeded 0.1/eps',
     &   'The iteration limit was reached',
     &   'Aprod  does not define a symmetric matrix',
     &   'Msolve does not define a symmetric matrix',
     &   'Msolve does not define a pos-def preconditioner',
     &   'x is not a descent direction' /
!     ------------------------------------------------------------------

      eps       = rw(  1) ! unit round-off.  IEEE DP  2.22e-16

!     Print heading and initialize.

      bnorm  = dnrm2 ( n, b, 1 )
      beta1  = bnorm
      if (nout .gt. 0) then
         write(nout, 1000) enter, beta1,
     &                     n, CheckA, Goodb, Precon,
     &                     itnlim, rtol, shift
      end if
      iStop  = 0
      itn    = 0
      normA  = zero
      aCond  = zero
      rnorm  = zero
      ynorm  = zero

      call dload ( n, zero, x, 1 )

!     Set up y for the first Lanczos vector v1.
!     y is really beta1 * P * v1  where  P = C**(-1).
!     y and beta1 will be zero if b = 0.

      call dcopy ( n, b, 1, y , 1 )
      call dcopy ( n, b, 1, r1, 1 )
      if ( Precon ) call Msolve( n, r1, y )
      if ( Goodb  ) then
         b1  = y(1)
      else
         b1  = zero
      end if
      beta1  = ddot  ( n, r1, 1, y, 1 )

!     See if Msolve is symmetric.

      if (CheckA  .and.  Precon) then
         call Msolve( n, y, r2 )
         s      = ddot  ( n, y, 1, y, 1 )
         t      = ddot  ( n,r1, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333D+0
         if (z .gt. epsa) then
            iStop = 7
            go to 900
         end if
      end if

!     Test for an indefinite preconditioner.

      if (beta1 .lt. zero) then
         iStop = 8
         go to 900
      end if

!     If b = 0 exactly, stop with x = 0.

      if (beta1 .eq. zero) then
         go to 900
      end if

!     Here and later, v is really P * (the Lanczos v).

      beta1  = sqrt( beta1 )
      s      = one / beta1

      call dcopy ( n,    y, 1, v, 1 )
      call dscal ( n, s, v, 1 )

      call Aprod
     &   ( n, v, y,
     &     Hprod, Hprod1,                              ! Added for SNOPT
     &     neA, nlocA, locA, indA, Acol,               ! Added for SNOPT
     &     neH, nlocH, locH, indH, Hcol,               ! Added for SNOPT
     &     cu, lencu, iu, leniu, ru, lenru,            ! Added for SNOPT
     &     cw, lencw, iw, leniw, rw, lenrw )           ! Added for SNOPT

!     See if Aprod  is symmetric.

      if (CheckA) then
         call Aprod
     &      ( n, y, r2,
     &        Hprod, Hprod1,                           ! Added for SNOPT
     &        neA, nlocA, locA, indA, Acol,            ! Added for SNOPT
     &        neH, nlocH, locH, indH, Hcol,            ! Added for SNOPT
     &        cu, lencu, iu, leniu, ru, lenru,         ! Added for SNOPT
     &        cw, lencw, iw, leniw, rw, lenrw )        ! Added for SNOPT
         s      = ddot  ( n, y, 1, y, 1 )
         t      = ddot  ( n, v, 1,r2, 1 )
         z      = abs( s - t )
         epsa   = (s + eps) * eps**0.33333D+0
         if (z .gt. epsa) then
            iStop = 6
            go to 900
         end if
      end if

!     Set up y for the second Lanczos vector.
!     Again, y is beta * P * v2  where  P = C**(-1).
!     y and beta will be zero or very small if b is an eigenvector.

      call daxpy ( n, (- shift), v, 1, y, 1 )
      alfa   = ddot  ( n, v, 1, y, 1 )
      call daxpy ( n, (- alfa / beta1), r1, 1, y, 1 )

!     Make sure  r2  will be orthogonal to the first  v.

      z      = ddot  ( n, v, 1, y, 1 )
      s      = ddot  ( n, v, 1, v, 1 )
      call daxpy ( n, (- z / s), v, 1, y, 1 )

      call dcopy ( n, y, 1, r2, 1 )
      if ( Precon ) call Msolve( n, r2, y )
      oldb   = beta1
      beta   = ddot  ( n, r2, 1, y, 1 )

      if (beta .lt. zero) then
         iStop = 8
         go to 900
      end if

!     Cause termination (later) if beta is essentially zero.

      beta   = sqrt( beta )
      if (beta .le. eps) then
         iStop = -1
      end if

!     See if the local reorthogonalization achieved anything.

      denom  = sqrt( s ) * dnrm2( n, r2, 1 )  +  eps
      s      = z / denom
      t      = ddot  ( n, v, 1, r2, 1 ) / denom
      if (nout .gt. 0  .and.  Goodb) then
         write(nout, 1100) beta1, alfa, s, t
      end if

!     Initialize other quantities.

      cgnorm = beta1
      gbar   = alfa
      dbar   = beta
      rhs1   = beta1
      rhs2   = zero
      bstep  = zero
      snprod = one
      tnorm  = alfa**2 + beta**2
      ynorm2 = zero
      gmax   = abs( alfa ) + eps
      gmin   = gmax

      if ( Goodb ) then
         call dload ( n, zero, w, 1 )
      else
         call dcopy ( n, v, 1, w, 1 )
      end if

!     ------------------------------------------------------------------
!     Main iteration loop.
!     ------------------------------------------------------------------
!+    Repeat                                          Until iStop .ne. 0

!        Estimate various norms and test for convergence.

  100    normA  = sqrt( tnorm  )
         ynorm  = sqrt( ynorm2 )
         epsa   = normA * eps
         epsx   = normA * ynorm * eps
         epsr   = normA * ynorm * rtol
         diag   = gbar
         if (diag .eq. zero) diag = epsa

         lqnorm = sqrt( rhs1**2 + rhs2**2 )
         qrnorm = snprod * beta1
         cgnorm = qrnorm * beta / abs( diag )

!        Estimate  cond(A).
!        In this version we look at the diagonals of  L  in the
!        factorization of the tridiagonal matrix,  T = L*Q.
!        Sometimes, T(k) can be misleadingly ill-conditioned when
!        T(k+1) is not, so we must be careful not to overestimate aCond.

         if (lqnorm .le. cgnorm) then
            aCond  = gmax / gmin
         else
            denom  = min( gmin, abs( diag ) )
            aCond  = gmax / denom
         end if

!        See if any of the stopping criteria are satisfied.
!        In rare cases, iStop is already -1 from above (Abar = const * I).

         if (iStop .eq. 0) then
            if (itn    .ge. itnlim    ) iStop = 5
            if (aCond  .ge. 0.1d+0/eps) iStop = 4
            if (epsx   .ge. beta1     ) iStop = 3
            if (cgnorm .le. epsx      ) iStop = 2
            if (cgnorm .le. bnorm*rtol) iStop = 1 ! Replaces next line
!           if (cgnorm .le. epsr      ) iStop = 1 ! for inexact Newton.
         end if

!        ===============================================================
!        See if it is time to print something.

         if ((n      .le. 40             .or.
     &        itn    .le. 10             .or.
     &        itn    .ge. itnlim - 10    .or.
     &        mod(itn,10) .eq.   0       .or.
     &        cgnorm .le. 10.0d+0*epsx   .or.
     &        cgnorm .le. 10.0d+0*epsr   .or.
     &        aCond  .ge. 0.01d+0/eps    .or.
     &        iStop  .ne. 0                  ) .and.
     &        nout   .gt. 0                  ) then

!           Print a line for this iteration.

            zbar   = rhs1 / diag
            z      = (snprod * zbar  +  bstep) / beta1
!           x1lq   = x(1)  +  b1 * bstep / beta1
            x1cg   = x(1)  +  w(1) * zbar  +  b1 * z

            if (    itn     .eq. 0) write(nout, 1200)
            write(nout, 1300) itn, x1cg, cgnorm, bstep/beta1,normA,aCond
            if (mod(itn,10) .eq. 0) write(nout, 1500)
         end if
!        ===============================================================

!        Obtain the current Lanczos vector  v = (1 / beta)*y
!        and set up  y  for the next iteration.

         if (iStop .eq. 0) then
            s      = one / beta

            call dcopy ( n,    y, 1, v, 1 )
            call dscal ( n, s, v, 1 )

            call Aprod
     &         ( n, v, y,
     &           Hprod, Hprod1,                        ! Added for SNOPT
     &           neA, nlocA, locA, indA, Acol,         ! Added for SNOPT
     &           neH, nlocH, locH, indH, Hcol,         ! Added for SNOPT
     &           cu, lencu, iu, leniu, ru, lenru,      ! Added for SNOPT
     &           cw, lencw, iw, leniw, rw, lenrw )     ! Added for SNOPT

            call daxpy ( n, (- shift), v, 1, y, 1 )
            call daxpy ( n, (- beta / oldb), r1, 1, y, 1 )
            alfa   = ddot( n, v, 1, y, 1 )
            call daxpy ( n, (- alfa / beta), r2, 1, y, 1 )
            call dcopy ( n, r2, 1, r1, 1 )
            call dcopy ( n, y, 1, r2, 1 )

            if ( Precon ) call Msolve( n, r2, y )
            oldb   = beta
            beta   = ddot  ( n, r2, 1, y, 1 )
            if (beta .lt. zero) then
               iStop = 6
               go to 800
            end if
            beta   = sqrt( beta )
            tnorm  = tnorm  +  alfa**2  +  oldb**2  +  beta**2

!           Compute the next plane rotation for  Q.

            gamma  = sqrt( gbar**2 + oldb**2 )
            cs     = gbar / gamma
            sn     = oldb / gamma
            delta  = cs * dbar  +  sn * alfa
            gbar   = sn * dbar  -  cs * alfa
            epsln  = sn * beta
            dbar   =            -  cs * beta

!           Update  x.

            z      = rhs1 / gamma
            s      = z * cs
            t      = z * sn

            do i = 1, n
               x(i) = (w(i) * s   +   v(i) * t)  +  x(i)
               w(i) =  w(i) * sn  -   v(i) * cs
            end do

!           Accumulate the step along the direction  b,
!           and go round again.

            bstep  = snprod * cs * z  +  bstep
            snprod = snprod * sn
            gmax   = max( gmax, gamma )
            gmin   = min( gmin, gamma )
            ynorm2 = z**2  +  ynorm2
            rhs1   = rhs2  -  delta * z
            rhs2   =       -  epsln * z
            itn    = itn   +  1
         end if

!+    until    (iStop .ne. 0)
      if (.not.(iStop .ne. 0)) go to 100

!     ------------------------------------------------------------------
!     End of main iteration loop.
!     ------------------------------------------------------------------

!     Move to the CG point if it seems better.
!     In this version of SYMMLQ, the convergence tests involve
!     only cgnorm, so we're unlikely to stop at an LQ point,
!     EXCEPT if the iteration limit interferes.
!
!     April 8, 2004. Always move to the CG point for SNOPT application

  800 zbar   = rhs1 / diag
      bstep  = snprod * zbar  +  bstep
      ynorm  = sqrt( ynorm2  +  zbar**2 )
      rnorm  = cgnorm
      call daxpy ( n, zbar, w, 1, x, 1 )

      if ( Goodb ) then

!        Add the step along  b.

         bstep  = bstep / beta1
         call dcopy ( n, b, 1, y, 1 )
         if ( Precon ) call Msolve( n, b, y )
         call daxpy ( n, bstep, y, 1, x, 1 )
      end if

!     ==================================================================
!     Display final status.
!     ==================================================================
  900 if (nout  .gt. 0) then
         write(nout, 2000) exit, iStop, itn,
     &                     exit, normA, aCond,
     &                     exit, rnorm, ynorm
         write(nout, 3000) exit, msg(iStop)
      end if

      return

!     ------------------------------------------------------------------
 1000 format(// 1p,    a, 5x, 'Solution of symmetric   Ax = b', 7x,
     &          'bnorm  =', e10.2
     &       / ' n      =', i7, 5x, 'CheckA =', l4, 12x,
     &          'Goodb  =', l4, 7x, 'Precon =', l4
     &       / ' itnlim =', i7, 5x, 'rtol   =', e11.2, 5x,
     &          'shift  =', e23.14)
 1100 format(/ 1p, ' beta1  =', e10.2, 3x, 'alpha1 =', e10.2
     &       / ' (v1,v2) before and after ', e14.2
     &       / ' local reorthogonalization', e14.2)
 1200 format(// 5x, 'itn', 7x, 'x1(cg)', 10x,
     &         'norm(r)', 5x, 'bstep', 7x, 'norm(A)', 3X, 'cond(A)')
 1300 format(1p, i8, e19.10, e11.2, e14.5, 2e10.2)
 1500 format(1x)
 2000 format(/ 1p, a, 6x, 'istop =', i3,   15x, 'itn   =', i8
     &       /     a, 6x, 'normA =', e12.4, 6x, 'aCond =', e12.4
     &       /     a, 6x, 'rnorm =', e12.4, 6x, 'ynorm =', e12.4)
 3000 format(      a, 6x, a )
!     ------------------------------------------------------------------
!     end of SYMMLQ
      end

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Msolv
     &   ( n, x, y )

      integer
     &     n
      double precision
     &     x(n), y(n)

!     ==================================================================
!     s5Msolv  dummy Msolve for SYMMLQ.
!
!     04 Dec 2004: First version of s5Msolv.
!     04 Dec 2004: Current version.
!     ==================================================================

!     Relax

      end ! subroutine s5Msolv
