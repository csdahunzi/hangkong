!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     file  sn55qp.f
!
!     s5QP
!     s5checkp  s5chzq   s5getp   s5QPfg          s5QPitn   s5Rcheck
!     s5Rcol    s5rg     s5Sdel   s5setCondZmax   s5ZHZ     s5ZHZeig
!     s5ZHZfac  s5Zp
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5QP
     &   ( iExit,
     &     probType, probTag, subOptimize,
     &     qpLog, Hprod, Hprod1, HvCalls, eigH,
     &     Elastic, GotR, NeedLU, typeLU, Needx,
     &     lenR, m, maxS, mBS, n, nb, nDegen,
     &     ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &     itQP, itQPmax, itQPtarget, itn,
     &     eMode, lvlObjE, printLevel,
     &     minimize, iObj, scaleObj, objAdd, objQP,
     &     condZmax, condZHZmax, tolOptFP, tolOptQP, tolx,
     &     nInf, sInf, elastics, nInfE, sInfE, wtInf,
     &     piNorm, rgNorm,
     &     neA, nlocA, locA, indA, Acol,
     &     neH, nlocH, locH, indH, Hcol,
     &     eType, eState, feasType, hs, kBS,
     &     bl, bu, blQP, buQP, blBS, buBS,
     &     gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg,
     &     nrhs0, nrhs, rhs, scales,
     &     lenx0, nx0, x0, x, xBS, xFrozen,
     &     iy, iy1, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     Elastic, GotR, NeedLU, Needx
      external
     &     Hprod, Hprod1, qpLog
      integer
     &     elastics, eigH, HvCalls, iExit, itQP, itQPmax, itQPtarget,
     &     itn, iObj, lencu, leniu, lenru, lencw, leniw, lenrw, lenR,
     &     lenx0, eMode, lvlObjE, m, maxS, mBS, minimize, ngObj0,
     &     ngObj, n, nb, neA, neH, ngQP0, ngQP, nlocA, nlocH, nInf,
     &     nInfE, nnH0, nnH, nS, nDegen, nrhs0, nrhs, nx0, probType,
     &     printLevel, subOptimize, typeLU, locA(nlocA), locH(nlocH),
     &     indA(neA), indH(neH), eType(nb), eState(nb), hs(nb),
     &     kBS(mBS), feasType(mBS), iy(nb), iy1(nb),
     &     iu(leniu), iw(leniw)
      double precision
     &     condZmax, objAdd, objQP, piNorm, rgNorm, scaleObj,
     &     sInf, sInfE, condZHZmax, tolOptFP, tolOptQP, tolx, wtInf,
     &     Acol(neA), bl(nb), bu(nb), blQP(nb), buQP(nb), rc(nb),
     &     blBS(mBS), buBS(mBS), gBS(mBS), gObj(*), gQP(ngQP0),
     &     Hcol(neH), Hdx(nnH0), pBS(mBS), pi(m), rhs(nrhs0),
     &     R(lenR), rg(maxS), scales(nb), x0(lenx0), x(nb), xBS(mBS),
     &     xFrozen(nb), y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     probTag*20, cu(lencu)*8, cw(lencw)*8

      !=================================================================
      ! s5QP   solves a linear or quadratic program.
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
      !   -2         QP is unbounded
      !   -3         Too many iterations
      !   -4         Weak QP minimizer
      !   -5         Too many superbasics
      !   -6         QP Hessian not positive semidefinite after pricing
      !   -7         Z'g could not be made sufficiently small
      !
      ! 30 Sep 1991: First version of s5QP  based on Qpsol routine qpcore.
      ! 29 Oct 1993: QP objective computed separately.
      ! 19 May 1995: Bordered Hessian updated.
      ! 30 Jul 1995: Border updates removed.
      ! 04 Jan 1996: Positive semi-definite H treated correctly.
      ! 20 Jul 1996: Slacks changed to be the row value.
      ! 09 Aug 1996: First Min Sum version.
      ! 15 Jul 1997: Thread-safe version.
      ! 02 Feb 1998: Piecewise linear line search added.
      ! 07 Nov 1998: Explicit Hessian option added.
      ! 24 Dec 1999: Sub-optimization option added.
      ! 25 Jul 2001: Exit on nS > maxR activated.
      ! 30 Jul 2003: Superbasic slacks allowed.
      ! 02 Aug 2003: snEXIT and snPRNT adopted.
      ! 24 Dec 2003: pi checked for NaN and Inf entries.
      ! 07 May 2006: s4ksave handles negative values of hs.
      ! 16 May 2006: Explicit target itQP added
      ! 26 May 2013: infBnd used to identify infinite bounds.
      ! 02 Nov 2014: Refactor B and ZHZ before iExit = -6.
      ! 03 Nov 2014: neH, indH, locH, Hcol added as arguments.
      ! 22 Nov 2014: iExit set correctly when eMode = 0.
      ! 02 May 2015: Redundant variable BndSwap removed.
      ! 02 May 2015: Basis no longer refactorized on return to phase 1.
      !=================================================================
      character
     &     str*132
      external
     &     dnormi
      logical
     &     Checkx, CheckFeas, Checkpi, Deadpoint,
     &     Feasible, GotgQP, GotH, Increase, FirstFeas, JustPhase1,
     &     LUok, Optimal, Stationary, MaxRef, Needf, Needv, NeedLM,
     &     Needpi, NewB, NewLU, NewSB, Newx, Prtlvl10,
     &     PrintLog, PrintSum, QPdone, Rcheck, UsegQP
      integer
     &     eigZHZ, inform, itnfix, itnlim, jq, jBq, jBr, jSq, jSr,
     &     jqSave, kchk, kDegen, kfac, klog, kObj, ksav, kSumm, kp,
     &     kPrc, kPrPrt, linesP, linesS, lenL0, lenL, lenU0, lenU,
     &     LUitn, LUmod, LUsiz0, LUmax, LUrequest, lvlTol, maxR,
     &     mnrHdP, mnrHdS, mNewSB, nBS, nFac, nfmove, frozen,
     &     nonOpt, nParPr, nParPrLP, nParPrQP,
     &     nSmax, nSwap, nUncon, printP, printS, sqStat,
     &     rankZHZ, runTime, toldj1, toldj2, toldj3, nfix(2)
      double precision
     &     Bgrowth, Bold, c6, condZHZ,
     &     djq0, djq, djqPrt, dnormi, dRmax, dRmin,
     &     eps0, eps2, featol, ZHZmin, ZHZmax, infBnd,
     &     maxTime, normAj, normg, objPrint,
     &     objSlack, pivot, rgTest, rgTol(2), rowError,
     &     signObj, step, tolx0, tolinc, weight
!     ------------------------------------------------------------------
      integer            CHOL
      parameter         (CHOL   = 0)
      integer            loose,      tight
      parameter         (loose  = 1, tight  = 2)
      integer            Intern,     Extern
      parameter         (Intern = 0, Extern = 1)
      integer            mUncon,     Check
      parameter         (mUncon = 3, Check  = 1)
      integer            WithB
      parameter         (WithB  = 1)
      integer            FP,         FPE,        FPS
      parameter         (FP     = 0, FPE    = 3, FPS    = 4)
      integer            BS        , BT
      parameter         (BS     = 2, BT     = 3)
      integer            Init,       Optml,      Cycle
      parameter         (Init   = 0, Optml  = 1, Cycle  = 2)
      integer            SEMDEF,     POSDEF
      parameter         (SEMDEF = 0, POSDEF = 1)

      parameter         (nParPr  =  94) ! partial pricing in use
      parameter         (nParPrLP=  99) ! partial pricing for LPs
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
      parameter         (mnrHdP  = 223) ! >0 => Minor head in log file
      parameter         (mnrHdS  = 225) ! >0 => Minor head for iSumm
      parameter         (runTime = 462) ! Solve time

      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      eps0      = rw(  2) ! eps**(4/5)       IEEE DP  3.00e-13
      eps2      = rw(  4) ! eps**(1/2)       IEEE DP  1.49e-08
      infBnd    = rw( 70) ! definition of an infinite bound
      maxTime   = rw( 79) ! max time allowed

      maxR      = iw( 52) ! max columns of R.
      kchk      = iw( 58) ! check (row) frequency
      kfac      = iw( 59) ! factorization frequency
      ksav      = iw( 60) ! save basis map
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      kDegen    = iw( 63) ! max. expansions of featol
      itnlim    = iw( 89) ! limit on total iterations
      mNewSB    = iw( 95) ! max # of new superbasics
      nFac      = iw(210) ! # of LU factorizations

      iExit     = 0
      sqStat    = 0

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
!     s5QP operates in either ``Normal'' or ``Elastic'' mode.
!     Everything is normal unless a weighted sum is being minimized or
!     the constraints are infeasible.
!     The logical Feasible refers to the feasibility of the nonelastics.
!     wtInf  is the optional parameter Infeasibility Weight.
!     ------------------------------------------------------------------
      Feasible   = .false.

      GotH       = nnH  .gt. 0
      GotgQP     = ngQP .gt. 0

!     JustPhase1 = stop at the end of phase1 (either regular or elastic)

      JustPhase1 = probType .eq. FP   .or.
     &             probType .eq. FPE  .or.
     &             probType .eq. FPS

!     The phase 2 objective is F1 + wtInf*F2.

      if (Elastic) then
         Needf = lvlObjE .ne. 2 ! F1 required in phase 2
         Needv = lvlObjE .ne. 0 ! F2 required in phase 2
      else
         Needf = .true.
         Needv = .false.
      end if

      objQP    =  zero
      objSlack =  zero
      pivot    =  zero
      step     =  zero

      sInfE    =  zero
      nInfE    =  0
      nonOpt   = -1

      jq       =  0
      djq      =  zero
      djq0     =  zero
      djqPrt   =  zero
      jBq      =  0             ! x(jBq) is the incoming  BS
      jBr      =  0             ! x(jBr) is the outgoing  BS
      jSq      =  0             ! x(jSq) is the incoming SBS
      jSr      =  0             ! x(jSr) is the outgoing SBS
      jqSave   =  0
      kPrPrt   =  0
      signObj  =  minimize

      rgTol(loose) = min( tolOptQP, c6) ! relaxed   rgTol
      rgTol(tight) = eps0               ! stringent rgTol
      lvlTol = tight                    ! working   rgTol

      rw(toldj1) = 100.0d+0
      rw(toldj2) = 100.0d+0

      if (JustPhase1) then
         iw(nParPr) = iw(nParPrLP)
      else
         iw(nParPr) = iw(nParPrQP)
      end if

      kPrc      = 0             ! last sec scanned in part. prc
      LUrequest = 0

      eigZHZ    =  SEMDEF

      CheckFeas = .true.        ! Check that x is feasible.
      Checkpi   = .true.        ! Check norm pi in next call of s5setpi
      Deadpoint = .false.
      Needpi    = .true.
      NewLU     = .true.
      Newx      = .false.
      QPdone    = .false.

!     nUncon  counts the number of unconstrained (i.e., Newton) steps.
!             If the test for a minimizer were scale-independent,
!             Uncon would never be larger than 1.
!     nfmove  counts the number of times that the QP obj is decreased,

      nfmove = 0
      nUncon = 0

!     subOptimize ne 0 forces optimization with a subset of the
!     variables frozen at their initial values.

      frozen = 0
      nSmax  = nS + mNewSB
      call s5hs  ( Intern, nb, blQP, buQP, hs, x )
      call s5degen
     &   ( inform, Init, printLevel, nb, nInf, itn,
     &     featol, tolx, tolinc, hs, blQP, buQP, x,
     &     itnfix, nfix, tolx0, iw, leniw, rw, lenrw )

!     ======================Start of main loop==========================
!+    do while (.not. QPdone  .and.  iExit .eq. 0)
  100 if       (.not. QPdone  .and.  iExit .eq. 0) then

         !==============================================================
         ! Check the initial  x  and move it onto  ( A  -I )*x = b.
         ! If NeedLU is true, this will require a basis factorization.
         !==============================================================
         ! If necessary,  factorize the basis  ( B = LU ) and set x.
         ! If NeedLU is false on entry to s5QP, the first call to s2Bfac
         ! will try to use existing factors.
         ! If NeedLU is true on entry to s5QP, an LU factorization of
         ! type typeLU is computed.
         !
         ! The reason for the requested LU is as follows.
         !
         ! LUrequest =  0  First LU for a given subproblem
         ! LUrequest =  1  Frequency
         ! LUrequest =  2  LU nonzeros increased
         ! LUrequest =  3
         ! LUrequest =  4
         ! LUrequest =  5  Singular after LU mod
         ! LUrequest =  6  Unstable LU mod (growth in new column of U)
         ! LUrequest =  7  Not enough memory
         ! LUrequest =  8
         ! LUrequest =  9
         ! LUrequest = 10  Row error      in setx
         ! LUrequest = 11  Big  dx or pi  in setx or setpi
         !
         ! LUrequest = 20
         ! LUrequest = 21  Iterative refinement failed in QP
         ! LUrequest = 22  Unbounded QP
         ! LUrequest = 23  Infeasibility after refactorization
         ! LUrequest = 24  Small directional derivative in QP
         ! LUrequest = 25  Ill-conditioned Z in QP
         ! LUrequest = 26  Z'HZ not positive definite in QP
         ! LUrequest = 27  R singular at a QP stationary point.
         !--------------------------------------------------------------
         FirstFeas = .false.
         if (LUrequest .gt. 0) NeedLU = .true.

         if (Needx  .or.  NeedLU) then
            call s2Bfac
     &         ( iExit,
     &           typeLU, NeedLU, NewLU, NewB,
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
               GotR   = .false. ! Reset R.
               if (Prtlvl10) iw(mnrHdP) = 1
            end if

            if (iExit .ne. 0) go to 100

            Needpi    = .true.  ! Recalculate the pi's.
            NeedLU    = .false.
            Needx     = .false.
            Newx      = .true.
            CheckFeas = .true.
            Checkpi   = .true.  ! Check for NaNs when getting piNorm

            pivot     = zero
            jqSave    = 0
            nUncon    = 0
         end if

         NewSB   = .false.
         Optimal = .false.

         nBS     = m + nS

         nInf    = 0
         sInf    = zero

         call dload ( nBS, zero, gBS, 1 )
         normg  = one

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
!           FirstFeas  is turned off once a step is taken.

            call s5Inf
     &         ( nBS, featol, infBnd,
     &           nInf, sInf, feasType, blBS, buBS, gBS, xBS )

            if (nInf .gt. 0) then

!              Non-elastics are infeasible.
!              If necessary, switch back to the feasibility phase, after
!              refactorization (possibly with tighter tols).
!              Print something if the basis has just been refactorized.

               if (Prtlvl10  .and.  iw(LUitn) .eq. 0) then
                  write(str, 1030) itn, nInf, sInf
                  call snPRNT( 21, str, iw, leniw )
               end if
               Feasible = .false.
               GotR     = .false. ! Is this needed?
            end if

!           Feasible => the nonelastics are feasible.
!           Feasible => normal or elastic Phase 2

            if (.not. Feasible) then
               FirstFeas = nInf .eq. 0 ! Feasible for the first time
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
     &         ( nb, nBS, eState, kBS, featol, nInfE, sInfE, bl, bu, x )
          end if

         if (Feasible  .and.  JustPhase1) then
            ! The non-elastics are feasible, prepare to exit.
            condZHZ = zero
            djqPrt  = zero
            rgNorm  = zero
            call dload ( m, zero, pi, 1 )
            piNorm  = one       ! piNorm = max(norm(pi), 1.0)
         else

            if (Feasible) then
!              ---------------------------------------------------------
!              Feasible for the nonelastics.
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

                     if (GotH  .and. .not. GotR) then
!                       ------------------------------------------------
!                       Load and factor the reduced Hessian.
!                       This happens after every LU factorize.
!                       If the reduced Hessian is not positive definite,
!                       reduce the LU factor tolerances to get a better
!                       conditioned Z.
!                       ------------------------------------------------
                        if (nS .gt. 0) then
                           call s5ZHZ
     &                        ( inform,
     &                          Hprod, Hprod1, HvCalls,
     &                          maxR, lenR, minimize, m, mBS,
     &                          n, nb, nnH, nS,
     &                          neA, nlocA, locA, indA, Acol,
     &                          neH, nlocH, locH, indH, Hcol,
     &                          condZmax, ZHZmin, ZHZmax,
     &                          kBS, R, y, y1, y2,
     &                          cu, lencu, iu, leniu, ru, lenru,
     &                          cw, lencw, iw, leniw, rw, lenrw )

!                          inform  Meaning
!                          ------  -------
!                            >0    Fatal LU error
!                             0    Normal exit. Z'HZ computed
!                            -1    Large estimate of cond(Z)

                           if (inform .ne. 0) then
                              if (inform .gt. 0) then
                                 iExit = inform
                              else
                                 call s5setCondZmax
     &                              ( condZmax, rw, lenrw )
                                 call s2tryLU
     &                              ( itn, 11, nS,LUrequest,LUok,typeLU,
     &                                iw, leniw, rw, lenrw )
                                 if (.not. LUok) then
                                    iExit = 44 ! Large cond(Z) estimate.
                                 end if
                              end if
                              go to 100
                           end if

                           call s5ZHZfac
     &                        ( inform,
     &                          eigH, CHOL, itn, lenR, m,
     &                          maxR, mBS, nb, nS,
     &                          condZHZmax, ZHZmin, ZHZmax, rankZHZ,
     &                          hs, kBS, iy,
     &                          blQP, buQP, blBS, buBS, x, xBS, R, y1,
     &                          iw, leniw )

!                          inform  Meaning
!                          ------  -------
!                            >0    Fatal LU error
!                             0    Normal exit. R computed
!                            -1    Z'HZ appears to be singular
!                            -2    Z'HZ appears to be indefinite

                           if (inform .ne. 0) then
                              if (inform .gt. 0) then
                                 iExit = inform
                              else
                                 call s2tryLU
     &                              ( itn, 26, nS,LUrequest,LUok,typeLU,
     &                                iw, leniw, rw, lenrw )
                                 if (.not. LUok) then
                                    iExit = -6 ! Z'HZ indefinite
                                 end if
                              end if
                              go to 100
                           end if
                        end if ! nS > 0
                        GotR = .true.
                     end if ! GotH and not GotR
                  end if ! Needf

                  Newx     = .false. ! objQP and gQP are now updated
                  nBS      = m + nS
                  if (GotR  .or.  nS .eq. 0) then
                     eigZHZ = POSDEF
                  else
                     eigZHZ = SEMDEF
                  end if
               end if ! FirstFeas .or. Newx

               Rcheck = .false. !  Rcheck = .true. recomputes R.
               if (Rcheck  .and.  GotR) then
                  call s5Rcheck
     &               ( iExit,
     &                 Hprod, Hprod1, HvCalls, eigZHZ,
     &                 itn, minimize,
     &                 maxR, lenR, m, mBS, n, nb, nnH, nS,
     &                 neA, nlocA, locA, indA, Acol,
     &                 neH, nlocH, locH, indH, Hcol,
     &                 kBS, R, y, y1, pBS,
     &                 cu, lencu, iu, leniu, ru, lenru,
     &                 cw, lencw, iw, leniw, rw, lenrw )
                  if (iExit .ne. 0) go to 100
               end if

!              ---------------------------------------------------------
!              Gather the QP gradient in BS order.
!              Assign the nonzero components of gBS.
!              ---------------------------------------------------------
               if (Needf) then
                  if (GotgQP) then
                     call s2gather
     &                  ( ngQP, nBS, kBS, signObj, gQP, gBS )
                  end if

                  if (iObj .gt. 0) then
                     gBS(iw(kObj)) = signObj*scaleObj
                  end if
               end if

               if (Elastic  .and.  nInfE .gt. 0  .and.  Needv) then
                  call s5eGrad
     &               ( nb, nBS, wtInf, eState, kBS, gBS )
                end if

               normg = max(dnormi( nBS, gBS, 1 ), one)

!              ---------------------------------------------------------
!              See if it's time to suboptimize.
!              NOTE: We must not suboptimize if all steps have been
!              degenerate.
!              ---------------------------------------------------------
               if (subOptimize .ne. 0  .or.  nfmove .eq. 0) then
!                 Relax
               else
                  if (nS .ge. nSmax) then
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

            if (Needpi) then
               call dcopy ( m, gBS, 1, y, 1 )
               call s5setpi
     &            ( inform,
     &              m, Checkpi, condZmax, normg, piNorm, y, pi,
     &              iw, leniw, rw, lenrw )
               Checkpi = .false.

               if (inform .ne. 0) then
                  if (inform .gt. 0) then
                     iExit     =  inform
                  else          ! pi is large or contains a NaN/Inf.
                     LUrequest = -inform
                     call s2tryLU
     &                  ( itn, LUrequest, nS, LUrequest, LUok, typeLU,
     &                    iw, leniw, rw, lenrw )
                     if (.not. LUok) iExit = 44 ! Large cond(Z) estimate
                  end if
                  go to 100
               end if
               Needpi = .false.
            end if

            rgNorm = zero
            if (nS .gt. 0) then
               call s5rg
     &            ( m, nBS, n, nS, eps0,
     &              neA, nlocA, locA, indA, Acol,
     &              gBS, pi, rg, rgNorm, kBS )
            end if

!           ============================================================
!           Determine if the reduced Hessian is positive definite.
!           ============================================================
            condZHZ = zero
            if (GotR) then
               call s6Rcnd
     &            ( maxR, nS, lenR, R, dRmax, dRmin, condZHZ )
            end if

!           ============================================================
!           Check for a stationary point.  Use a stringent rgTol after
!           a constrained step to help avoid false stationary points.
!
!           If x is a subspace minimizer,  reduced costs are calculated.
!           ============================================================
            if (Feasible) then
               rw(toldj3) = tolOptQP
            else
               rw(toldj3) = tolOptFP
            end if

            rgTest = max( piNorm, normg )

            if (.not. Feasible) then
               Stationary = rgNorm .le. rgTol(loose) *rgTest
            else if (nUncon .ge. 1) then
               Stationary = rgNorm .le. rgTol(lvlTol)*rgTest
            else
               Stationary = rgNorm .le. rgTol(tight) *rgTest
            end if

            if (Feasible) then

               MaxRef = nUncon .gt. mUncon

               if (MaxRef  .and.  .not. Stationary) then

!                 If this point should be stationary but isn't.
!                 If possible, relax the reduced-gradient tolerance.

                  if (lvlTol .eq. tight) then
                     lvlTol     = loose
                     Stationary = rgNorm .le. rgTol(lvlTol)*rgTest
                  end if

                  if (MaxRef  .and. .not. Stationary) then
                     call s2tryLU
     &                  ( itn, 21, nS, LUrequest, LUok, typeLU,
     &                    iw, leniw, rw, lenrw )
                     if (.not. LUok) then
                        iExit = -7 ! Large Z'g
                     end if
                     go to 100
                  end if
               end if
               Deadpoint = Stationary .and. Needf
     &                                .and. eigZHZ .eq. SEMDEF
            end if

            if (Stationary) then
               jqSave  = 0
               nUncon  = 0

               if (GotR  .and.  eigZHZ .eq. SEMDEF) then
                  write(str, 1600) itn
                  call snPRNT( 23, str, iw, leniw )
                  call s2tryLU
     &               ( itn, 27, nS, LUrequest, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )
                  if (.not. LUok) iExit = 44 ! Large cond(Z) estimate.
                  go to 100
               end if
            end if

            NeedLM = Stationary

            kPrPrt = kPrc
            jq     = 0

            if (NeedLM) then
!              ---------------------------------------------------------
!              Compute Lagrange multipliers.
!              ---------------------------------------------------------
               djq0   = djq     ! save djq in case of bad statpt
               djq    = zero
               nUncon = 0
               UsegQP = Feasible  .and.  Needf  .and.  GotgQP
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
            end if ! NeedLM
         end if ! Feasible and JustPhase1

         QPdone = Optimal .or. Deadpoint .or.(Feasible .and. JustPhase1)

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
               !--------------------------------------------------------
               ! So far so good.  Now check the row residuals.
               !--------------------------------------------------------
               if (iw(LUitn) .gt. 0) then
                  call s5setx
     &               ( inform, Check, itn,
     &                 m, n, nb, nBS, rowError,
     &                 neA, nlocA, locA, indA, Acol,
     &                 kBS, xBS, nrhs0, nrhs, rhs, x, y, y1,
     &                 iw, leniw, rw, lenrw )
                  QPdone    = inform .eq. 0
                  LUrequest = inform
                  if (LUrequest .gt. 0) typeLU = BS
               end if
            end if

            if (QPdone) then
               if (Deadpoint) iExit = -4
            else
               Needx  = .true.
               Needpi = .true.
               go to 100
            end if
         end if ! done

         if (FirstFeas  .and.  JustPhase1) then
!           Relax, we are about to exit without printing anything.
         else
!           ============================================================
!           Print the details of this iteration.
!           ============================================================
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
     &         ( probType, probTag,
     &           Elastic, GotR, FirstFeas, Feasible, JustPhase1,
     &           m, mBS, nnH, nS, jSq, jBr, jSr,
     &           iw(linesP), iw(linesS), itn, itQP, kPrPrt, lvlObjE,
     &           pivot, step, nInf, sInf, nInfE, sInfE, wtInf,
     &           nonOpt, objPrint, condZHZ, djqPrt, rgNorm, kBS, xBS,
     &           iw, leniw )
         end if
         jBq    = 0
         jBr    = 0
         jSq    = 0
         jSr    = 0
         kPrPrt = 0
         djqPrt = zero

         if (QPdone) then
!           ------------------------------------------------------------
!           Optimal .or. Deadpoint  .or. (Feasible .and. JustPhase1)
!           ------------------------------------------------------------
            if (nInf .gt. 0) then

               ! No feasible point for the nonelastic variables
               ! Stop or continue in elastic mode, depending on the
               ! specified level of infeasibility.

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
                  CheckFeas = .true.         ! call s5eReset
                  QPdone    = .false.

                  Needf     = lvlObjE .ne. 2 ! Need F1 in e-phase 2
                  Needv     = lvlObjE .ne. 0 ! Need F2 in e-phase 2
                  Needpi    = .true.
                  djq       = zero
                  step      = zero
               end if
               go to 100
            end if

            if (Prtlvl10 .and. .not. JustPhase1) then
               if (jq .ne. 0) then
                  djqprt = signObj*djq
                  if (Prtlvl10  .and.  klog .eq. 1) then
                     write(str, 1010) djq, jq, rgnorm, piNorm
                     call snPRNT( 11, str, iw, leniw )
                  end if
               else
                  if (Prtlvl10  .and.  klog .eq. 1) then
                     write(str, 1020)          rgnorm, piNorm
                     call snPRNT( 11, str, iw, leniw )
                  end if
               end if
            end if
         else
            ! ----------------------------------------------------------
            ! Do another QP iteration.
            ! A nonbasic has been selected to become superbasic.
            ! Compute the vector y such that B y = column jq.
            ! ----------------------------------------------------------
            if (NewSB) then
               !--------------------------------------------------------
               ! The price has selected a nonbasic to become superbasic.
               !--------------------------------------------------------
               if (nS+1 .gt. maxR) then
                  iExit = -5
                  go to 100
               end if

               lvlTol = tight
               djqPrt = djq

!              ---------------------------------------------------------
!              Compute the vector pBS such that B pB = column jq.
!              pBS is a multiple of part of the new column of  Z  and
!              is used to define the QP search direction and update R.
!              ---------------------------------------------------------
!              Unpack column jq into  y1  and solve  B*y = y1.
!              The solve computes  y1  such that  L*y1 = ajq.
!              It is used below to modify L and U in s5QPitn.

               call s2unpack
     &            ( jq, m, n, neA, normAj, nlocA, locA, indA, Acol, y1 )
               call s2Bsol
     &            ( iExit, WithB, m, y1, pBS, iw, leniw, rw, lenrw )
               if (iExit .ne. 0) return
            end if

!           ============================================================
!           Take a step.
!           ============================================================
            if (itn  .ge. itnlim  .or.  itQP .ge. itQPmax) then
               iExit = -3
               go to 100
            end if

            itQP   = itQP   + 1
            itn    = itn    + 1

!           Decide if we want to print something this iteration.

            PrintLog = printLevel .ge. 1 .and. mod(itQP, klog ) .eq. 0
            PrintSum = printLevel .ge. 1 .and. mod(itQP, kSumm) .eq. 0

            iw(printP) = 0
            iw(printS) = 0
            if (PrintLog) iw(printP) = 1
            if (PrintSum) iw(printS) = 1

!           ------------------------------------------------------------
!           Take a reduced gradient step.
!           The new  x  will either minimize the objective on the
!           working set or lie on the boundary of a new constraint.
!           ------------------------------------------------------------
            call s5QPitn
     &         ( inform,
     &           Hprod, Hprod1, HvCalls, eigH, eigZHZ,
     &           Elastic, Feasible,
     &           GotgQP, GotH, GotR, Increase, Needf, Needv,
     &           Needpi, NewSB, itn, lenR,
     &           m, mBS, maxR, maxS, n, nb,
     &           nnH0, nnH, nS, ngQP0, ngQP, nDegen, elastics,
     &           LUrequest, kp, jBq, jSq, jBr, jSr,
     &           jq, jqSave, nfmove, nUncon,
     &           djq0, djq, minimize, objQP,
     &           featol, pivot, step, tolinc, wtInf,
     &           neA, nlocA, locA, indA, Acol,
     &           neH, nlocH, locH, indH, Hcol,
     &           eType, eState, feasType, hs, kBS,
     &           bl, bu, blQP, buQP, blBS, buBS, gBS,
     &           gQP, Hdx, pBS, rg, R,
     &           x, xBS, y, y1, y2,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

!           inform  Meaning                             Action
!           ------  -------                             ------
!             -5    Refinement step hit deleted bound.  (LUrequest)
!             -4    Bad directional derivative.         (exit)
!             -3    Unable to update the LU factors.    (LUrequest)
!             -2    Z'HZ not positive semidefinite      (exit)
!             -1    Unbounded                           (exit)
!              0    Normal exit
!             >0    Fatal LU error                      (exit)
!
            if (inform .ne. 0) then
               iExit = 0
               if      (inform .gt.  0) then
                  iExit = inform ! Fatal LU error
               else if (inform .eq. -1) then
                  iExit = -2     ! unbounded
               else if (inform .eq. -2) then
                  iExit = -6     ! Z'HZ appears to be indefinite
               else if (inform .eq. -4) then
                  iExit = -7     ! Bad QP  directional derivative.
               end if
               if (iExit .ne. 0) go to 100
            end if

            if (LUrequest  .gt.  0) then

!                  5   Singular after LU mod
!                  6   Unstable LU mod
!                  7   Insufficient free memory
!                 24   small directional deriv

               call s2tryLU
     &            ( itn, LUrequest, nS, LUrequest, LUok, typeLU,
     &              iw, leniw, rw, lenrw )
               if (.not. LUok) then
                  iExit = 43    ! Cannot satisfy Ax - s = b
                  go to 100
               end if
            end if

            iw(LUitn) = iw(LUitn)  + 1

!           Increment featol every iteration.

            featol = featol + tolinc

!           ============================================================
!           Test for error condition and/or frequency interrupts.
!           ============================================================
!           (1) Save a basis map (frequency controlled).
!           (2) Every kdegen iterations, reset featol and move nonbasic
!               variables onto their bounds if they are very close.
!           (3) Refactorize the basis if it has been modified too many
!               times.
!           (4) Update the LU factors of the basis if requested.
!           (5) Check row error (frequency controlled).

            if (mod(itn,ksav) .eq. 0) then
               call s4ksave
     &            ( minimize, m, n, nb, nS, mBS,
     &              itn, nInf, sInf, objQP, kBS, hs,
     &              scales, blQP, buQP, x, xBS, cw, lencw, iw, leniw )
            end if

            if (mod( itn, kdegen ) .eq. 0) then
               call s5degen
     &            ( inform, Cycle, printLevel, nb, nInf, itn,
     &              featol, tolx, tolinc, hs, blQP, buQP, x,
     &              itnfix, nfix, tolx0,
     &              iw, leniw, rw, lenrw )
               Needx  = inform .gt. 0
            end if

            if (LUrequest .eq. 0) then
               if (     iw(LUmod) .ge. kfac-1) then
                  LUrequest = 1
               else if (iw(LUmod) .ge. 20  .and.
     &                                iw(lenL)+iw(lenU) .gt. LUmax) then
                  Bgrowth = iw(lenL) + iw(lenU)
                  Bold    = LUsiz0
                  Bgrowth = Bgrowth/Bold
                  if (Prtlvl10) then
                     write(str, 1000) Bgrowth
                     call snPRNT( 21, str, iw, leniw )
                  end if
                  LUrequest = 2
               else
                  Checkx = mod(iw(LUitn),kchk) .eq. 0
                  if (Checkx  .and.  .not. Needx) then
                     call s5setx
     &                  ( inform, Check, itn,
     &                    m, n, nb, nBS, rowError,
     &                    neA, nlocA, locA, indA, Acol,
     &                    kBS, xBS, nrhs0, nrhs, rhs, x, y, y1,
     &                    iw, leniw, rw, lenrw )
                     LUrequest = inform ! 0, 10, 11, 142 (s2Bsol)
                  end if
               end if
               if (LUrequest .gt. 0) typeLU = BT
            end if

            ! Check time limit

            if (maxTime .gt. zero .and. mod(itn, 20) .eq. 0) then
               call s1time ( -2, 0, iw, leniw, rw, lenrw )
               call s1time (  2, 0, iw, leniw, rw, lenrw )
               if (rw(runTime) .gt. maxTime) then
                  iExit = 34    ! time limit reached
               end if
            end if

         end if ! not done

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
 1600 format(' Itn', i7, ': Singularity after a ',
     &                   'bound swap.  Basis refactorized')
 1610 format(' Itn', i7, ': Suboptimize: ', i7, ' new superbasics')
 1620 format(' Itn', i7, ': Suboptimize: ', i7, ' minor iterations')
 8050 format(' Itn', i7, ': Infeasible ', a)
 8060 format(' Itn', i7, ': Elastic Phase 1 -- making ',
     &                   'non-elastic variables feasible')

      end ! subroutine s5QP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5checkp
     &   ( iExit,
     &     itn, nBS, jqSave, kBS, gtp, pBS, iw, leniw )

      implicit
     &     none
      integer
     &     iExit, itn, nBS, jqSave, leniw, kBS(nBS), iw(leniw)
      double precision
     &     gtp, pBS(nBS)

!     ==================================================================
!     s5checkp  makes  pBS  a feasible direction.
!
!     If iExit = 1, the directional derivative is positive.
!
!     16 Jun 1995: First version of s5checkp.
!     02 Aug 2003: snPRNT adopted.
!     02 Aug 2003: Current version of s5checkp.
!     ==================================================================
      character
     &     str*80
      integer
     &     kSave, j, jq, k
      double precision
     &     pSave
!     ------------------------------------------------------------------
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      iExit = 0

!     ------------------------------------------------------------------
!     Find the element of  pBS  corresponding to the most recently freed
!     variable. Usually, it will be pBS(nBS).
!     ------------------------------------------------------------------
      jq    = abs(jqSave)

      kSave = 0
      do k =  nBS, 1, -1
         j = kBS(k)
         if (j .eq. jq) then
            kSave = k
            go to 100
         end if
      end do

!     ------------------------------------------------------------------
!     Choose the sign of  pBS  so that the most recently freed
!     variable continues to increase or decrease.
!     ------------------------------------------------------------------
  100 if (kSave .gt. 0) then
         pSave = pBS(kSave)

         if ((jqSave .lt. 0  .and.  pSave .gt. zero)  .or.
     &       (jqSave .gt. 0  .and.  pSave .lt. zero)      ) then
            call dscal ( nBS, (-one), pBS, 1 )
            gtp  = - gtp
         end if

         if (gtp .gt. zero) then
!           ------------------------------------------------------------
!           Looks as though the sign of gtp cannot be relied upon.
!           In later versions we'll fix this variable.
!           For now, we just print a warning and stop.
!           ------------------------------------------------------------
            write(str, 1000) itn, gtp
            call snPRNT( 23, str, iw, leniw )
            iExit = 1           ! Bad directional derivative
         end if
      else
!        ---------------------------------------------------------------
!        Couldn't find the index of the most recently freed variable.
!        This should never happen!
!        ---------------------------------------------------------------
         write(str, 9000) jqSave
         call snPRNT( 23, str, iw, leniw )
      end if

      return

 1000 format(' Itn', i7, ': Bad directional derivative ', 1p, e9.1 )
 9000 format(' XXX  s5checkp.  kSave not found. jqSave = ', i5 )

      end ! subroutine s5checkp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5chzq
     &   ( m, mBS, n, nb, nS, kBSq, pivot, tolpiv,
     &     neA, nlocA, locA, indA, Acol,
     &     kBS, bl, bu, xBS, y, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     m, mBS, n, nb, neA, nlocA, nS, kBSq, leniw, lenrw,
     &     locA(nlocA), indA(neA), kBS(mBS), iw(leniw)
      double precision
     &     pivot, tolpiv, Acol(neA), bl(nb), bu(nb), xBS(mBS), y(mBS),
     &     rw(lenrw)

!     ==================================================================
!     s5chzq  selects a superbasic to replace the kp-th basic variable.
!     On entry,  y  contains the kp-th row of B(inverse).
!     On exit, pivot and  y(m+1), ..., y(m+nS) define the S-part of
!     the modifying vector w.
!
!     01 Dec 1991: First version based on Minos routine m7chzq.
!     02 Aug 2003: snPRNT adopted.
!     30 Jun 2005: Current version of s5chzq.
!     ==================================================================
      character
     &     str*80
      integer
     &     j, k, m1, idamax
      double precision
     &     d1, d2, dpiv, eps0, tol, xj
!     ------------------------------------------------------------------
      double precision   zero,          point1,          one
      parameter         (zero = 0.0d+0, point1 = 0.1d+0, one = 1.0d+0)
      integer            Transp
      parameter         (Transp = 1)
!     ------------------------------------------------------------------
      eps0      = rw(  2) ! eps**(4/5)

!     Set yS = 0 -  S'*y.

      m1        = m  + 1
      call s2Bprod
     &   ( Transp, eps0, n, nS, kBS(m1),
     &     neA, nlocA, locA, indA, Acol,
     &     (-one), y, m, zero, y(m1), nS )

      kBSq   = m  +  idamax( nS, y(m1), 1 )
      pivot  = abs( y(kBSq) )

!     Exit if the pivot is too small.

      if (pivot .lt. tolpiv) then
         write(str, 1000)  pivot
         call snPRNT( 31, str, iw, leniw )
         kBSq   = - (m + nS)
      else

!        Choose one away from its bounds if possible.

         tol    =   point1*pivot
         dpiv   = - one

         do k = m1, m+nS
            if (abs( y(k) ) .ge. tol) then
               j     = kBS(k)
               xj    = xBS(k)
               d1    = xj    - bl(j)
               d2    = bu(j) - xj
               d1    = min( abs( d1 ), abs( d2 ) )
               if (dpiv .le. d1) then
                  dpiv  = d1
                  kBSq  = k
               end if
            end if
         end do

         pivot   = - y(kBSq)

      end if ! pivot .ge. tolpiv

      return

 1000 format(' XXX  s5chzq.  Max pivot is too small:', 1p, e11.1 )

      end ! subroutine s5chzq

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5getp
     &   ( Feasible, GotR, eigZHZ,
     &     maxR, lenR, nS, R, rg, p, gp, pHp )

      implicit
     &     none
      logical
     &     GotR, Feasible
      integer
     &     eigZHZ, maxR, lenR, nS
      double precision
     &     gp, pHp, R(lenR), rg(nS), p(nS)

!     ==================================================================
!     s5getp  computes a search direction  p  for the superbasic
!     variables, using the current reduced gradient  g.
!
!     29 Mar 2001: R stored by rows.
!     20 May 2001: Current version.
!     ==================================================================
      external
     &     ddot
      integer
     &     lRlast
      double precision
     &     delta, Rlast
      double precision
     &     ddot
!     ------------------------------------------------------------------
      integer            POSDEF
      parameter         (POSDEF = 1)
      integer            PULL
      parameter         (PULL   = 1)
      integer            WithR,      WithRt
      parameter         (WithR  = 0, WithRt = 1)

      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------

      call dcopy ( nS, rg, 1,  p, 1 )
      call dscal ( nS, (-one), p, 1 )

      if (Feasible  .and.  GotR) then
         if (eigZHZ .eq. POSDEF) then
            !--------------------------------------------------------------
            ! The Newton direction.
            !--------------------------------------------------------------
            call s6Rsol( WithRt, maxR, nS, lenR, R, p )
            gp   = -ddot( nS, p, 1, p, 1 )
            pHp  = -gp
            call s6Rsol( WithR , maxR, nS, lenR, R, p )
         else
            !--------------------------------------------------------------
            ! A direction of zero or negative curvature.
            !--------------------------------------------------------------
            lRlast    = (nS-1)*maxR + (3-nS)*nS/2 ! Magic formula!
            Rlast     = R(lRlast)                 ! the last diag of R.
            R(lRlast) = one

            delta     = Rlast*Rlast
            if (Rlast .ge. zero) then
               delta  =  Rlast*Rlast
            else
               delta  = -Rlast*Rlast
            end if

            call s6Rsol( WithRt, maxR, nS, lenR, R, p )
            if (nS .gt. 1)
     &      call dscal ( nS-1, delta, p, 1 )
            call s6Rsol( WithR , maxR, nS, lenR, R, p )

            gp        = -delta*ddot( nS, p, 1, p, 1 ) - p(nS)**2
            pHp       = -delta*gp
            R(lRlast) =  Rlast
         end if
      else
         !--------------------------------------------------------------
         ! A direction of steepest-descent.
         !--------------------------------------------------------------
         pHp = zero
         gp  = ddot( nS, rg, 1, p, 1 )
      end if ! Feasible and GotR

      if (gp .gt. zero) then
         call dscal ( nS, (-one), p, 1 )
         gp  = -gp
         pHp = -pHp
      end if

      end ! subroutine s5getp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5QPfg
     &   ( Hprod, Hprod1, HvCalls,
     &     ngQP, ngObj0, ngObj, nnH,
     &     neH, nlocH, locH, indH, Hcol,
     &     sqStat, fQP,
     &     gObj, gQP, lenx0, nx0, x0, x, dx,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      integer
     &     HvCalls, lenx0, lencu, lencw, leniu, leniw, lenru, lenrw,
     &     neH, ngQP, ngObj0, ngObj, nlocH, nnH, nx0, sqStat,
     &     indH(neH), locH(nlocH), iu(leniu), iw(leniw)
      double precision
     &     fQP, gObj(ngObj0), gQP(ngQP), Hcol(neH),
     &     x0(lenx0), x(ngQP), dx(ngQP), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5QPfg  computes various quantities associated with the LP/QP.
!
!       1.  fQP =  gObj'*(x-x0)  + half*(x - x0)'*H*(x - x0)
!       2.  gQP =  gradient of fQP
!
!     On entry,
!     ngQP         is max( ngObj, nnH )
!     x(ngQP)      are the nonlinear variables
!     x0(ngQP)     is the base point x0 (scaled if necessary)
!     gObj(ngObj)  defines the explicit QP linear term
!
!     On exit,
!     fQP          is the QP quadratic term (1) above
!     gQP(ngQP)    is the gradient of fQP
!     dx(ngQP)     is  x-x0
!
!     02 May 1992: First version of s5QPfg.
!     23 Oct 1993: Hx added as an argument.
!     29 Oct 1993: Modified to compute only the QP objective.
!     07 Oct 1994: gQP added as an argument.
!     09 Dec 2004: f77+ version.
!     02 Mar 2013: HvCalls added.
!     03 Nov 2014: neH, indH, locH, Hcol added as arguments.
!     ==================================================================
      external
     &     ddot
      integer
     &     nzero
      double precision
     &     ddot
!     ------------------------------------------------------------------
      double precision   zero,          half,          one
      parameter         (zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      if (ngQP .le. 0) return

      call dcopy ( ngQP,         x , 1, dx, 1 )
      if (nx0 .gt. 0)
     &call daxpy ( ngQP, (-one), x0, 1, dx, 1 )

      fQP  = zero

      if (nnH .gt. 0) then
         call Hprod
     &      ( Hprod1, nnH,
     &        neH, nlocH, locH, indH, Hcol,
     &        dx, gQP, sqStat,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         HvCalls = HvCalls + 1
         fQP     = half*ddot( nnH, dx, 1, gQP, 1 )
      end if

      nzero = ngQP - nnH
      if (nzero .gt. 0) call dload ( nzero, zero, gQP(nnH+1), 1 )

      if (ngObj .gt. 0) then
         fQP = fQP + ddot( ngObj, gObj, 1,  dx, 1 )
         call daxpy ( ngObj, one, gObj, 1, gQP, 1 )
      end if

      end ! subroutine s5QPfg

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5QPitn
     &   ( iExit,
     &     Hprod, Hprod1, HvCalls, eigH, eigZHZ,
     &     Elastic, Feasible,
     &     GotgQP, GotH, GotR, Increase, Needf, Needv,
     &     Needpi, NewSB, itn, lenR,
     &     m, mBS, maxR, maxS, n, nb,
     &     nnH0, nnH, nS, ngQP0, ngQP, nDegen, elastics,
     &     LUrequest, kp, jBq, jSq, jBr, jSr,
     &     jq, jqSave, nfmove, nUncon,
     &     djq0, djq, minimize, objQP,
     &     featol, pivot, step, tolinc, wtInf,
     &     neA, nlocA, locA, indA, Acol,
     &     neH, nlocH, locH, indH, Hcol,
     &     eType, eState, feasType, hs, kBS,
     &     bl, bu, blQP, buQP, blBS, buBS, gBS,
     &     gQP, Hdx, pBS, rg, R,
     &     x, xBS, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      logical
     &     Elastic, Feasible, GotgQP, GotH, GotR, Increase,
     &     Needf, Needv, Needpi, NewSB
      integer
     &     eigH, eigZHZ, elastics, HvCalls, iExit, itn, jBq, jBr,
     &     jq, jqSave, kp, lenR, lencu, lencw, leniu, leniw, lenru,
     &     lenrw, LUrequest, m, maxR, maxS, mBS, minimize, n, nb,
     &     nDegen, neA, neH, nfmove, nlocA, nlocH, nnH0, nnH,
     &     ngQP0, ngQP, nS, nUncon, locA(nlocA), locH(nlocH),
     &     indA(neA), indH(neH), eType(nb), eState(nb), hs(nb),
     &     feasType(mBS), kBS(mBS), iu(leniu), iw(leniw)
      double precision
     &     djq0, djq, objQP, featol, pivot, step, tolinc, wtInf,
     &     Acol(neA), bl(nb), bu(nb), blQP(nb), buQP(nb), blBS(mBS),
     &     buBS(mBS), gBS(mBS), gQP(ngQP0), Hcol(neH), Hdx(nnH0),
     &     pBS(mBS), R(lenR), rg(maxS), x(nb), xBS(mBS), y(nb),
     &     y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5QPitn performs a QP step.
!
!     On entry,
!        NewSB = true implies that variable jq just went superbasic.
!                In this case:
!                pBS  satisfies B pBS = a(jq).
!                y1   satisfies L  y1 = a(jq).
!
!     On exit,
!        pBS contains the most recent QP search direction.
!
!      iExit       Result
!      -----       ------
!       -5         Refinement step hit deleted bound.
!       -4         Bad directional derivative.
!       -3         Unable to update the LU factors.
!       -2         Z'HZ not positive semidefinite
!       -1         unbounded
!        0         normal exit
!       >0         Fatal LU error
!
!     25 Nov 1991: First version of s5QPitn.
!     05 Jan 1996: Positive semidefinite R treated correctly.
!     29 Aug 1996: First min sum version added.
!     27 Jul 1997: Thread-safe version.
!     02 Feb 1998: Piecewise linear line search added.
!     23 Mar 2000: gQP  and  H  scaled.
!     16 Oct 2000: Reverted to non-bordered version of s5QPitn.
!     04 Dec 2000: R converted to row-wise storage.
!     02 Aug 2003: snPRNT adopted.
!     07 May 2006: s5Zp added to compute Z*p.
!     08 Apr 2008: eState accessed in elastic mode only.
!     04 Jul 2008: Modify both bl and bu in elastic phase 1.
!     02 Mar 2013: HvCalls added.
!     02 Nov 2014: Set LUrequest = 25 when Z'HZ is indefinite.
!     03 Nov 2014: neH, indH, locH, Hcol added as arguments.
!     10 Dec 2014: Removed unused argument obj.
!     06 May 2015: x(jr) set exactly on its bound.
!     23 May 2015: Z'HZ not positive definite triggers refactorization.
!     ==================================================================
      character
     &     str*80
      external
     &     ddot
      logical
     &     HitCon, HitLow, Move, OnBound, TestRd, Unbounded, Uncon,
     &     Rcheck
      integer
     &     inform, infpiv, jqeState, jqState, jr, jreState,
     &     jrState, jsq, jsr, kBSq, ksq, lRlast, lRs, LUitn, LUmod,
     &     mtry, nBS, nBS1, nS1, ntry, sqStat
      double precision
     &     normA, bigdx, bound, ddot, eps, eps0, exact,
     &     gp, gpQP, infBnd, pBS1, pHp, pHpQP, pNorm, Rlast, sclPiv,
     &     signObj, stepB, stepMax, stepP, tolpiv, tolP0, tolP
!     ------------------------------------------------------------------
      integer            INDEF,      SEMDEF,     POSDEF
      parameter         (INDEF  =-1, SEMDEF = 0, POSDEF = 1)
       integer            xBStox
      parameter         (xBStox = 1)
      parameter         (mtry   = 6)
      integer            WithL,      WithBt
      parameter         (WithL  = 0, WithBt = 2)

      parameter         (LUitn  = 215) ! itns since last factorize
      parameter         (LUmod  = 216) ! number of LU mods

      double precision   zero,            half,          one
      parameter         (zero   = 0.0d+0, half = 0.5d+0, one = 1.0d+0)
      double precision   ten
      parameter         (ten    =10.0d+0)
!     ------------------------------------------------------------------
      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      eps0      = rw(  2) ! eps**(4/5)          IEEE DP  3.00e-13
      tolpiv    = rw( 60) ! excludes small elements of pBS.
      infBnd    = rw( 70) ! definition of an infinite bound
      bigdx     = rw( 72) ! unbounded step.

      iExit     = 0
      LUrequest = 0
      sqStat    = 0

      Unbounded = .false.
      signObj   = minimize

      nBS       = m + nS

      if (NewSB) then
!        ---------------------------------------------------------------
!        New superbasic.
!        Definite must be true if there is a new superbasic.
!        ---------------------------------------------------------------
         nS1        = nS   + 1
         nBS1       = nBS  + 1

         kBS (nBS1) =      jq
         xBS (nBS1) =    x(jq)
         blBS(nBS1) = blQP(jq)
         buBS(nBS1) = buQP(jq)
         jqState    =   hs(jq)

         eigZHZ     = SEMDEF

         if (GotR) then
!           ------------------------------------------------------------
!           Add the new column to R at position nS+1.
!           Check for a singular or indefinite reduced Hessian.
!           ------------------------------------------------------------
            kBS(nBS1) = jq
            call s5Rcol
     &         ( iExit,
     &           Hprod, Hprod1, HvCalls,
     &           minimize, jq, nS1, Rlast, lRlast,
     &           maxR, lenR, m, mBS, n, nb, nnH, nS1,
     &           neA, nlocA, locA, indA, Acol,
     &           neH, nlocH, locH, indH, Hcol,
     &           kBS, R, y, y2, pBS,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) go to 900

            ! Test sign(Rlast) in s5ZHZeig

            TestRd = .true.
            call s5ZHZeig
     &         ( TestRd, eigH, eigZHZ, itn,
     &           maxR, lenR, nS1, Rlast, R, iw, leniw, rw, lenrw )

            if (Feasible) then
               if (eigZHZ .ne. POSDEF) then
                  if (eigZHZ .eq. SEMDEF  .and.  iw(LUitn) .eq. 0  .or.
     &                rlast  .eq. zero  ) then
!                    Relax
                  else
                     iExit     = -2
                     LUrequest = 25
                     go to 900
                  end if
               end if
            else

!              ZHZ needs to be pos def in phase 1
!              If not, stop updating it.

               GotR   = eigZHZ .eq. POSDEF
               eigZHZ = SEMDEF ! ZHZ isn't used to define pBS
            end if
         end if ! GotR

         !-->  R may be checked here by setting  Rcheck = .true.

         Rcheck = .false.

         if (Rcheck  .and.  GotR) then
            call s5Rcheck
     &         ( iExit,
     &           Hprod, Hprod1, HvCalls, eigH,
     &           itn, minimize,
     &           maxR, lenR, m, mBS, n, nb, nnH, nS1,
     &           neA, nlocA, locA, indA, Acol,
     &           neH, nlocH, locH, indH, Hcol,
     &           kBS, R, y, y2, pBS,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) go to 900
         end if

         if (Increase) then
            jqSave =  jq
         else
            jqSave = -jq
         end if

         nS     =  nS1
         nBS    = nBS1

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
               else
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
!        In phase 1, or phase 2 for an LP, price can select nonbasics
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
      end if ! newSB

!     ------------------------------------------------------------------
!     Compute the search direction for the superbasics.
!     Store the free components of the search direction in pBS(1:nBS).
!     First, find the search direction pS for the superbasics, store
!     it in  pBS(m+1:nBS), and find its norm.  Store the search
!     direction for the basic variables in pBS(1)  ,...,pBS(m).
!     ------------------------------------------------------------------
  100 call s5getp
     &   ( Feasible, GotR, eigZHZ,
     &     maxR, lenR, nS, R, rg, pBS(m+1), gp, pHp )

      pBS1   = pBS(m+nS)

      call s5Zp
     &   ( iExit,
     &     m, mBS, n, nb, nS, eps0, pNorm,
     &     neA, nlocA, locA, indA, Acol,
     &     kBS, pBS, y2, iw, leniw, rw, lenrw )
      if (iExit .ne. 0) go to 900

      if (Feasible) then
         !--------------------------------------------------------------
         ! If R is singular, ensure that pBS is a feasible direction.
         ! A nonzero exit value of inform implies that the directional
         ! derivative is too small to be relied upon.
         !--------------------------------------------------------------
         if (GotR  .and.  eigZHZ .eq. SEMDEF) then
            call s5checkp
     &         ( inform, itn, nBS, jqSave, kBS, gp, pBS, iw, leniw )
            if (inform .gt. 0) then
               iExit     = -4
               LUrequest = 24
               go to 900
            end if
         end if

         if (NewSB  .and.  GotR) then
            ! Check for a feasible direction.
            ! A large  rgTol  may give a pBS(nBS) with the wrong sign.
            ! If so, continue minimizing with the old superbasic set.

            if (djq*pBS1 .gt. zero) then

               write(str, 1000) itn
               call snPRNT( 23, str, iw, leniw )

               nS     = nS  - 1
               nBS    = nBS - 1
               hs(jq) = jqState

               if (Elastic) then
!                 If a variable went elastic, reset it as not elastic.

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
               jqSave =  0
               eigZHZ =  POSDEF
               NewSB  = .false.
               go to 100
            end if
         end if

!        ---------------------------------------------------------------
!        Compute y = pBS(scattered) and Hdx(scattered).
!        The vector Hdx is used to update the objective and gradient of
!        the QP.  Form  gpQP  and  pHpQP  for the quadratic.
!        gp = gpQP - pBS(kObj) + terms from the elastic gradient.
!        ---------------------------------------------------------------
         if (Needf  .and.  (GotgQP  .or.  GotH)) then
            call s2scatter
     &         ( ngQP, nBS, kBS, one, pBS, y )

            if (GotgQP) then
               gpQP  = ddot ( ngQP, gQP, 1, y, 1 )
            end if

            if (GotH) then
               pHpQP = zero

               call Hprod
     &            ( Hprod1, nnH,
     &              neH, nlocH, locH, indH, Hcol,
     &              y, Hdx, sqStat,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )
               HvCalls = HvCalls + 1
               pHpQP   = pHpQP + ddot( nnH, y, 1, Hdx, 1 )
            end if
         end if
      end if ! Feasible

!     ------------------------------------------------------------------
!     Find the nearest constraint in direction  x + step*pBS (step > 0).
!     Exact  is the step that takes xBS(kp) exactly onto bound. It may
!     be positive or slightly negative. (Not defined if Unbounded.)
!
!     If OnBound  is true, step is a step that reaches a bound exactly.
!     xBS(kp) reaches the value bound.  If we take a constrained step,
!     bound is used to put the new nonbasic variable x(jr) exactly on
!     its bound.
!
!     If Unbounded is true, step = stepMax.
!     ------------------------------------------------------------------
      stepMax = bigdx /pNorm
      sclPiv  = one
      tolP0   = tolpiv
      tolP    = tolpiv*pNorm
      ntry    = 0

!+    Repeat
  200    tolP  = tolP /sclPiv
         tolP0 = tolP0/sclPiv
         call s5step
     &      ( nBS, nDegen,
     &        featol, infBnd, stepMax, tolinc, tolP,
     &        feasType, blBS, buBS, xBS, pBS,
     &        HitLow, Move, OnBound, Unbounded,
     &        infpiv, kp, bound, exact, stepB, stepP )

!        Find if the step is constrained or unconstrained.
!        If R has been flagged as singular, we double check by trying
!        to compute the QP minimizer along pBS.  If the minimizer
!        exists,  the singularity tolerance must be too large.

         if (Feasible) then
            if (eigZHZ .eq. POSDEF) then
               Uncon = stepP .gt. one

            else if (eigZHZ .eq. SEMDEF) then
               if (pHp .le. zero) then
                  Uncon = .false.
               else
                  Uncon = pHp .gt. zero  .and.  stepP*pHp .gt. (- gp)
               end if
            end if

            Unbounded = (Unbounded .and. .not. Uncon) .or.
     &                   stepMax .le. one
         else ! infeasible
            Uncon = .false.
         end if

         sclPiv = ten
         ntry   = ntry + 1

!+    until    ( infpiv .eq. 0 .and. (.not.Unbounded .or. Feasible) .or.
!+                 ntry .ge. mtry)
      if (.not.((infpiv .eq. 0 .and. (.not.Unbounded .or. Feasible)).or.
     &             ntry .ge. mtry)) go to 200

      if (Unbounded) then
         iExit = -1
         go to 900
      end if

      HitCon = .not. Uncon
      Needpi = .true.

      if (HitCon) then
         if (abs(jqSave) .eq. kBS(kp)  .and.  nUncon .gt. 0) then
            LUrequest = 28
            iExit     = -5
            go to 900
         end if
         nUncon = 0
         step   = stepB
      else
         pivot  = zero
         if (eigZHZ .eq. POSDEF) then
            nUncon = nUncon + 1
            step   = one
         else if (eigZHZ .eq. SEMDEF) then
            nUncon = nUncon + 1
            step   = (- gp)/pHp
            eigZHZ = POSDEF
         else
            step   = one
         end if
      end if

      !-----------------------------------------------------------------
      ! Compute the new objective function.
      ! Note: pHp = signObj*pHpQP
      !-----------------------------------------------------------------
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

      if (HitCon) then
!        ===============================================================
!        There is a blocking variable.
!        It could be a fixed variable, whose new state must be 4.
!        ===============================================================
         pivot   = -pBS(kp)
         jr      =  kBS(kp)

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
               if (blQP(jr) .eq. buQP(jr)) then
                  jrstate =  4
               else if (OnBound) then
                  jrstate =  0
               else if (x(jr) .lt. buQP(jr)) then
                  jrstate = -1
               else
                  jrstate =  1
               end if
            else ! jreState =  2
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
!           If nS = 1 there is no choice.
!           ============================================================
            if (nS .eq. 1) then
               kBSq  = nBS
               pivot = pivot/pBS1
            else
               call dload ( m, zero, y2, 1 )
               y2(kp) = one     ! Set      y2 = ep
                                ! Solve  B'yB = ep
               call s2Bsol
     &            ( iExit,
     &              WithBt, m, y2, y, iw, leniw, rw, lenrw )
               if (iExit .ne. 0) return
               call s5chzq
     &            ( m, mBS, n, nb, nS, kBSq, pivot, tolP0,
     &              neA, nlocA, locA, indA, Acol,
     &              kBS, blQP, buQP, xBS, y, iw, leniw, rw, lenrw )
               if (kBSq .le. 0) then
                  write(str, 9999) itn
                  call snPRNT( 23, str, iw, leniw )
                  kBSq   = nBS
               end if
            end if

            kSq        = kBSq - m

            hs(jr)     = jrState
            jBr        = jr                     ! Outgoing basic
            jSr        = kBS(kBSq)              ! Outgoing superbasic
            kBS (kBSq) = jBr
            jBq        = jSr                    ! Incoming basic
            kBS (kp)   = jSr
            blBS(kp)   = blBS(kBSq)
            buBS(kp)   = buBS(kBSq)
            xBS (kp)   = xBS (kBSq)
            hs(jBq)    = 3

            if (nS .gt. 1  .and.  GotR) then

!              Finish computing y(m+1), ..., y(m+nS).

               y(kBSq) = - (one + pivot)
               call dscal ( nS, (one/pivot), y(m+1), 1 )
               call s6Rswap( maxR, nS, lenR, R, y2, y(m+1), kSq, eps0 )
            end if

!           ------------------------------------------------------------
!           Get a new  y1, used to modify L and U.  If the outgoing
!           superbasic just came in, we already have it.
!           ------------------------------------------------------------
            if (jSr .ne. jq) then
               call s2unpack
     &            ( jBq, m, n, neA, normA, nlocA, locA, indA, Acol, y1 )
               call s2Bsol
     &            ( iExit,
     &              WithL, m, y1, y, iw, leniw, rw, lenrw )
               if (iExit .ne. 0) return
            end if

!           Update the LU factors.

            iw(LUmod)  = iw(LUmod) + 1

            call s2Bmod2
     &         ( inform, kp, m, y1, iw, leniw, rw, lenrw )

!           Check if LU factors got updated. If not, refactor on exit.

            if (inform .ne. 0) then
               if (inform .eq. -1) then
                  LUrequest = 5 ! Singular after LU mod
               else if (inform .eq.  2) then
                  LUrequest = 6 ! Unstable LU mod
               else if (inform .eq.  7) then
                  LUrequest = 7 ! Insufficient free memory
               end if
               iExit = -3
            end if
         else
!           ============================================================
!           A variable in S hit a bound.
!           ============================================================
            hs(jr) = jrState
            jSr    = jr
            kBSq   = kp
            kSq    = kBSq - m
         end if

!        Delete the kSq-th superbasic and adjust all arrays in BS order.

         call s5Sdel
     &      ( kSq, m, nS, nBS, kBS, blBS, buBS, gBS, rg, xBS )

         nS  = nS  - 1
         nBS = nBS - 1

         if (iExit .eq. 0) then
            if (GotR) then
!              ---------------------------------------------------------
!              Cyclically demote column kSq of R to position nS+1.
!              ---------------------------------------------------------
               if (kSq .lt. nS+1) then
                  call s6Rdel( kSq, maxR, nS+1, lenR, R, eps )
               end if
            end if ! Feasible and GotH

            if (eigZHZ .eq. SEMDEF) then

!              Recheck the last diagonal of R
!              It should not increase in magnitude

               if (Feasible  .and.  GotR) then
                  if (nS .gt. 0) then
                     lRs   = (nS-1)*maxR + (3-nS)*nS/2 ! Magic formula!
                     Rlast = R(lRs)
                  end if
                  TestRd = .false. ! ignore the sign of Rlast
                  call s5ZHZeig
     &               ( TestRd, eigH, eigZHZ, itn,
     &                 maxR, lenR, nS, Rlast, R, iw, leniw, rw, lenrw )
               else if (nS .eq. 0) then
                  eigZHZ = POSDEF
               else
                  eigZHZ = SEMDEF
               end if
            end if

*-->        R can be checked here.

         end if
      end if ! HitCon

  900 return

 1000 format(' Itn', i7, ': Bad direction after adding a superbasic.')
 9999 format(' Itn', i7, ': Choose q failed in s5QPitn!!')

      end ! subroutine s5QPitn

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Rcheck
     &   ( iExit,
     &     Hprod, Hprod1, HvCalls, eigH,
     &     itn, minimize,
     &     maxR, lenR, m, mBS, n, nb, nnH, nS,
     &     neA, nlocA, locA, indA, Acol,
     &     neH, nlocH, locH, indH, Hcol,
     &     kBS, R, v, w, y,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      integer
     &     eigH, HvCalls, iExit, itn, lencu, leniu, lenru, lencw, leniw,
     &     lenrw, lenR, maxR, m, minimize, mBS, n, nb, neA, neH,
     &     nlocA, nlocH, nnH, nS, locA(nlocA), locH(nlocH), indA(neA),
     &     indH(neH), kBS(mBS), iu(leniu), iw(leniw)
      double precision
     &     Acol(neA), Hcol(neH), R(lenR), v(nb), w(nb), y(mBS),
     &     ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s5Rcheck  computes the Cholesky factor R such that
*     R'R = Z'HZ.  The update corresponds to the addition of a new
*     column to Z.
*
*     On entry,
*        R     holds the columns of the factor associated with the
*              first jRadd-1 columns of Q.
*
*        nS    is the number of columns in R.
*
*     14 Mar 2001: First version (s6Rchk) based on SNOPT routine s5Rcol.
*     20 Jun 2008: Renamed s5Rcheck.
!     29 Jun 2008: Revamped.
!     02 Mar 2013: HvCalls added.
!     03 Nov 2014: neH, indH, locH, Hcol added as arguments.
*     ==================================================================
      external
     &     ddot
      logical
     &     TestRd
      integer
     &     eigZHZ, jq, jS, lRlast, lencol, nBS, status
      double precision
     &     normAj, dRsq, eps0, Rlast, Rnormsq, signObj, wHw, ddot
*     ------------------------------------------------------------------
      integer            INDEF
      parameter         (INDEF  =-1)
      integer            PUSH
      parameter         (PUSH  =  0)
      integer            Transp
      parameter         (Transp = 1)
      integer            WithRt
      parameter         (WithRt = 1)
      integer            WithB,      WithBt
      parameter         (WithB  = 1, WithBt = 2)
      double precision   zero,          one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      eps0    = rw(  2)

      nBS     = m   + nS
      signObj = minimize
      iExit   = 0

*     ------------------------------------------------------------------
*     Main loop to find a column of Z'HZ.
*     ------------------------------------------------------------------
      do jS = 1, nS
         lencol    = min( jS-1, nnH )
*        ---------------------------------------------------------------
*        Get the nonlinear elements of the column of Z.
*        Find y such that B y = column jq.
*        Scatter the nonlinear part of y into w.
*        ---------------------------------------------------------------
         jq  = kBS(m+jS)
         call s2unpack
     &      ( jq, m, n, neA, normAj, nlocA, locA, indA, Acol, w )
         call s2Bsol
     &      ( iExit,
     &        WithB, m, w, y, iw, leniw, rw, lenrw )
         if (iExit .ne. 0) return
         call s2scatter
     &      ( nnH, m, kBS, (-one), y, w )
         if (jq .le. nnH) w(jq) = one

*        ---------------------------------------------------------------
*        Compute  H*w  and  w'*H*w.
*        ---------------------------------------------------------------
         wHw = zero

         if (nnH .gt. 0) then
            status = 0
            call Hprod
     &         ( Hprod1, nnH,
     &           neH, nlocH, locH, indH, Hcol,
     &           w, v, status,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            HvCalls = HvCalls + 1
            wHw     = wHw + ddot ( nnH, w, 1, v, 1 )

            if (minimize .lt. 0) then
               call dscal ( nnH, signObj, v, 1 )
               wHw = signObj*wHw
            end if
         end if

         Rnormsq = zero

         if (jS .gt. 1) then
*           ------------------------------------------------------------
*           Gather the nonlinear elements of v in w (= vBS).
*           Compute Z'w  (solve  B'vB = wB and form  wS = wS - S'vB).
*           ------------------------------------------------------------
            call s2gather
     &         ( nnH, nBS, kBS, one, v, w )
            call s2Bsol
     &         ( iExit, WithBt, m, w, v, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) return

            if (nS .gt. 0) then
               call s2Bprod
     &            ( Transp, eps0, n, nS, kBS(m+1),
     &              neA, nlocA, locA, indA, Acol,
     &              (-one), v, m, one, w(m+1), nS )
            end if

*           ------------------------------------------------------------
*           Solve  R'v = Z(j)'Hw.  Store v in w(m+1:m+jS).
*           ------------------------------------------------------------
            call s6Rsol
     &         ( WithRt, maxR, lencol, lenR, R, w(m+1) )
            Rnormsq = ddot  ( lencol, w(m+1), 1, w(m+1), 1 )
         end if

         if (jS .le. nnH) then
            dRsq = wHw - Rnormsq ! square of new diagonal of R.
            if (dRsq .ge. zero) then
               Rlast =  sqrt(dRsq)
            else
               Rlast = -sqrt(abs(dRsq))
            end if
         else
            dRsq  = zero
            Rlast = zero
         end if
         w(m+jS)  = Rlast

*        Insert w(m+1:m+jS) as column jS of R.

         call s6Rcol
     &      ( PUSH, jS, maxR, jS, lenR, R, w(m+1), lRlast )

         TestRd = .true.        ! the sign of Rlast is relevant
         call s5ZHZeig
     &      ( TestRd, eigH, eigZHZ, itn,
     &        maxR, lenR, nS, Rlast, R, iw, leniw, rw, lenrw )

         if (eigZHZ .eq. INDEF) then
            iExit = 6
            return
         end if
      end do

      end ! subroutine s5Rcheck

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Rcol
     &   ( iExit,
     &     Hprod, Hprod1, HvCalls,
     &     minimize, jq, jRadd, Rlast, lRlast,
     &     maxR, lenR, m, mBS, n, nb, nnH, nS,
     &     neA, nlocA, locA, indA, Acol,
     &     neH, nlocH, locH, indH, Hcol,
     &     kBS, R, v, w, y,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      integer
     &     HvCalls, iExit, jq, jRadd, lRlast, lenR, lencu, leniu, lenru,
     &     lencw, leniw, lenrw, maxR, m, minimize, mBS, n, nb, neA, neH,
     &     nlocA, nlocH, nnH, nS, locA(nlocA), locH(nlocH), indA(neA),
     &     indH(neH), kBS(mBS), iu(leniu), iw(leniw)
      double precision
     &     Rlast, Acol(neA), Hcol(neH), R(lenR), v(nb), w(nb), y(mBS),
     &     ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5Rcol  computes column jRadd of the Cholesky factor R such that
!     R'R = Z'HZ.  The update corresponds to the addition of a new
!     column to Z.
!
!     On entry,
!        R     holds the columns of the factor associated with the
!              first jRadd-1 columns of Q.
!
!        y     is the vector such that B y = a(jq).
!
!        nS    is the number of columns in R.
!
!     11 Dec 1991: First version based on Qpsol routine Qpcolr.
!     24 Apr 1994: Columns of Nx no longer in Q.
!     27 Oct 2000: Previous version of s5Rcol.
!     04 Dec 2000: R converted to row-wise storage.
!     09 Dec 2004: Current version of s5Rcol.
!     02 Mar 2013: HvCalls added.
!     03 Nov 2014: neH, indH, locH, Hcol added as arguments.
!     ==================================================================
      external
     &     ddot
      integer
     &     lencol, nBS, sqStat
      double precision
     &     dRsq, eps0, Rnormsq, signObj, wHw, ddot
!     ------------------------------------------------------------------
      integer            PUSH
      parameter         (PUSH  =  0)
      integer            Transp
      parameter         (Transp = 1)
      integer            WithRt
      parameter         (WithRt = 1)
      integer            WithBt
      parameter         (WithBt = 2)
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      eps0      = rw(  2)

      iExit     = 0
      nBS       = m   + nS
      signObj   = minimize
      lencol    = min( jRadd-1, nnH )

!     ------------------------------------------------------------------
!     Get w, the vector of nonlinear components of the new column of Z.
!     ------------------------------------------------------------------
!     The input vector y satisfies B y = column jq.
!     Scatter the nonlinear components of y into w.

      call s2scatter
     &   ( nnH, m, kBS, (-one), y, w )
      if (jq .le. nnH) w(jq) = one

!     ------------------------------------------------------------------
!     Compute  H*w  and  w'*H*w.
!     ------------------------------------------------------------------
      wHw = zero

      if (nnH .gt. 0) then
         sqStat = 0
         call Hprod
     &      ( Hprod1, nnH,
     &        neH, nlocH, locH, indH, Hcol,
     &        w, v, sqStat,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         HvCalls = HvCalls + 1
         wHw     = wHw + ddot ( nnH, w, 1, v, 1 )

         if (minimize .lt. 0) then
            call dscal ( nnH, signObj, v, 1 )
            wHw = signObj*wHw
         end if
      end if

      Rnormsq = zero

      if (jRadd .gt. 1) then
!        --------------------------------------------------------------
!        Gather the nonlinear elements of v in w (= vBS).
!        Compute Z'w  (solve  B'vB = wB and form  wS = wS - S'vB).
!        --------------------------------------------------------------
         call s2gather
     &      ( nnH, nBS, kBS, one, v, w )
         call s2Bsol
     &      ( iExit,
     &        WithBt, m, w, v, iw, leniw, rw, lenrw  )
         if (iExit .ne. 0) return

         if (nS .gt. 0) then
            call s2Bprod
     &         ( Transp, eps0, n, nS, kBS(m+1),
     &           neA, nlocA, locA, indA, Acol,
     &           (-one), v, m, one, w(m+1), nS )
         end if

!        --------------------------------------------------------------
!        Solve  R'v = Z(j)'Hw.  Store v in w(m+1:m+jRadd).
!        --------------------------------------------------------------
         call s6Rsol
     &      ( WithRt, maxR, lencol, lenR, R, w(m+1) )
         Rnormsq = ddot  ( lencol, w(m+1), 1, w(m+1), 1 )
      end if

      if (jRadd .le. nnH) then
         dRsq = wHw - Rnormsq   ! square of the new diagonal of R.
         if (dRsq .ge. zero) then
            Rlast =  sqrt(dRsq)
         else
            Rlast = -sqrt(abs(dRsq))
         end if
      else
         dRsq  = zero
         Rlast = zero
      end if
      w(m+jRadd) = Rlast

!     Insert w(m+1:m+jRadd) as column jRadd of R.

      call s6Rcol
     &   ( PUSH, jRadd, maxR, jRadd, lenR, R, w(m+1), lRlast )

      end ! subroutine s5Rcol

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5rg
     &   ( m, nBS, n, nS, tolz,
     &     neA, nlocA, locA, indA, Acol,
     &     gBS, pi, rg, rgNorm, kBS )

      implicit
     &     none
      integer
     &     m, nBS, n, neA, nlocA, nS, locA(nlocA), indA(neA), kBS(nBS)
      double precision
     &     tolz, rgNorm, Acol(neA), gBS(nBS), pi(m), rg(nS)

!     ==================================================================
!     s5rg    calculates the reduced gradient  rg = gS - S'*pi.
!
!     23 Nov 1991: First version based on Minos routine m7rg.
!     16 Nov 2001: Current version.
!     ==================================================================
      external
     &     dnormi
      double precision
     &     dnormi
!     ------------------------------------------------------------------
      integer            Transp
      parameter         (Transp = 1)
      double precision   one
      parameter         (one = 1.0d+0)
!     ------------------------------------------------------------------
      call dcopy
     &   ( nS, gBS(m+1), 1, rg, 1 )

      call s2Bprod
     &   ( Transp, tolz, n, nS, kBS(m+1),
     &     neA, nlocA, locA, indA, Acol,
     &     (-one), pi, m, one, rg, nS )

      rgNorm = dnormi( nS, rg, 1 )

      end ! subroutine s5rg

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Sdel
     &   ( kSq, m, nS, nBS, kBS, blBS, buBS, gBS, rg, xBS )

      implicit
     &     none
      integer
     &     kSq, m, nS, nBS, kBS(nBS)
      double precision
     &     blBS(nBS), buBS(nBS), gBS(nBS), rg(nS), xBS(nBS)

!     ==================================================================
!     s5Sdel  deletes the kSqth superbasic variable from the arrays
!     kBS, blBS, blBS, gBS, rg and xBS.
!
!     16 Jun 2001: First version of s5Bswp.
!     16 Jun 2001: Current version.
!     ==================================================================
      integer
     &     j, k
!     ------------------------------------------------------------------
!     Shift all the arrays one place to the left.

      do j  = kSq, nS-1
         k  = m + j
         kBS (k) = kBS (k+1)
         blBS(k) = blBS(k+1)
         buBS(k) = buBS(k+1)
         gBS (k) = gBS (k+1)
         xBS (k) = xBS (k+1)
         rg (j)  = rg  (j+1)
      end do

      end ! subroutine s5Sdel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5setCondZmax
     &   ( condZmax, rw, lenrw )

      implicit
     &     none
      integer
     &     lenrw
      double precision
     &     condZmax, rw(lenrw)

!     ==================================================================
!     s5setCondZmax increases condZmax towards condZbnd.
!     In particular, condZmax = half*sqrt(condZmax*condZbnd)
!
!     12 May 2015: First version of s5setCondZmax.
!     ==================================================================
      double precision
     &     rtcondZbnd
!     ------------------------------------------------------------------
      double precision   half
      parameter         (half       = 0.5d+0)
!     ------------------------------------------------------------------
      rtcondZbnd = rw(192) ! square root of the max condition est of Z

      condZmax   = half*sqrt(condZmax)*rtcondZbnd

      end ! subroutine s5setCondZmax

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5ZHZ
     &   ( iExit,
     &     Hprod, Hprod1, HvCalls,
     &     maxR, lenR, minimize, m, mBS,
     &     n, nb, nnH, nS,
     &     neA, nlocA, locA, indA, Acol,
     &     neH, nlocH, locH, indH, Hcol,
     &     condZmax, ZHZmin, ZHZmax,
     &     kBS, R, v, w, y,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1
      integer
     &     HvCalls, iExit, lenR, lencu, leniu, lenru, lencw, leniw,
     &     lenrw, m, maxR, mBS, minimize, n, nb, neA, neH, nlocA, nlocH,
     &     nnH, nS, locA(nlocA), locH(nlocH), indA(neA), indH(neH),
     &     kBS(mBS), iu(leniu), iw(leniw)
      double precision
     &     condZmax, ZHZmax, ZHZmin, Acol(neA), Hcol(neH), R(lenR),
     &     y(nb), v(nb), w(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5ZHZ computes the reduced Hessian and loads it by columns into
!     the upper triangle R.
!
!      iExit       Status
!      -----       ------
!        0         reduced Hessian computed successfully
!       >0         Fatal error in LU solve
!
!     13 Oct 1992: First version based on QPSOL routine Qpcrsh.
!     15 Oct 1994: Dependent columns fixed at their current value.
!     04 Dec 2000: R converted to row-wise storage.
!     02 Mar 2013: HvCalls added.
!     01 Nov 2014: ZHZmin added as an argument.
!     01 Nov 2014: Removed exit based on large condZ.
!     03 Nov 2014: neH, indH, locH, Hcol added as arguments.
!     12 May 2015: condZmax added as argument.
!     ==================================================================
      external
     &     dnormi
      integer
     &     jq, jS, lRlast, nBS, sqStat
      double precision
     &     condZ, diag, dnormi, eps0, flmax, normAj, signObj
!     ------------------------------------------------------------------
      integer            Transp
      parameter         (Transp = 1)
      integer            WithB,           WithBt
      parameter         (WithB  = 1,      WithBt = 2)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      eps0    = rw(  2) ! eps**(4/5)       IEEE DP  3.00e-13
      flmax   = rw(  8) ! est. of the largest pos. real

      iExit   = 0

      ZHZmin  = flmax
      ZHZmax  = zero

      if (nS .eq. 0) return

      signObj = minimize
      nBS     = m + nS
      sqStat  = 0

!     ------------------------------------------------------------------
!     Main loop to find a column of Z'HZ.
!     ------------------------------------------------------------------
      do jS = 1, nS
         !--------------------------------------------------------------
         ! Get the nonlinear elements of the column of Z.
         ! Find y such that B y = column jq.
         ! Scatter the nonlinear part of y into w.
         !--------------------------------------------------------------
         jq  = kBS(m+jS)
         call s2unpack
     &      ( jq, m, n, neA, normAj, nlocA, locA, indA, Acol, w )
         call s2Bsol
     &      ( iExit, WithB, m, w, y, iw, leniw, rw, lenrw )
         if (iExit .gt. 0) return

         condZ = max ( dnormi( m, y, 1 )/normAj, one )
         if (condZ .ge. condZmax) then
            iExit = -1
            return
         end if

         call s2scatter
     &      ( nnH, m, kBS, (-one), y, w )
         if (jq .le. nnH) w(jq) = one

         ! Set v = H w.

         if (nnH .gt. 0) then
            call Hprod
     &         ( Hprod1, nnH,
     &           neH, nlocH, locH, indH, Hcol,
     &           w, v, sqStat,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            HvCalls = HvCalls + 1
            if (minimize .lt. 0) then
               call dscal ( nnH, signObj, v, 1 )
            end if
         end if

         !--------------------------------------------------------------
         ! Gather w = vBS and compute v = Z' w.
         ! Solve  B' vB = wB  and  form  wS = wS - S' vB.
         !--------------------------------------------------------------
         call s2gather
     &      ( nnH, nBS, kBS, one, v, w )
         call s2Bsol
     &      ( iExit,
     &        WithBt, m, w, v, iw, leniw, rw, lenrw )
         if (iExit .gt. 0) return

         call s2Bprod
     &      ( Transp, eps0, n, nS, kBS(m+1),
     &        neA, nlocA, locA, indA, Acol,
     &        (-one), v, m, one, w(m+1), nS )

         !--------------------------------------------------------------
         ! Store w(1:nS) in the jS-th row of R.
         ! R is NO LONGER SYMMETRIZED.
         !--------------------------------------------------------------
         call s6Rrow
     &      ( jS, maxR, nS, lenR, R, w(m+1), lRlast )
         diag   = R(lRlast)
         ZHZmin = min( ZHZmin, abs(diag) )
         ZHZmax = max( ZHZmax, abs(diag) )
      end do

      end ! subroutine s5ZHZ

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5ZHZeig
     &   ( TestRd, eigH, eigZHZ, itn,
     &     maxR, lenR, nS, Rlast, R, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     TestRd
      integer
     &     eigH, eigZHZ, itn, maxR, lenR, nS, leniw, lenrw, iw(leniw)
      double precision
     &     Rlast, R(lenR), rw(lenrw)

      !=================================================================
      ! s5ZHZeig  estimates the inertia of the current reduced Hessian.
      !
      ! On entry,
      !   eigH   encodes the inertia of the Hessian H
      !              eigH   Hessian
      !             -----   ---------------------
      !              -1                indefinite
      !               0     positive semidefinite
      !               1     positive semidefinite
      !
      !   Rlast   is the last diagonal element of R
      !   TestRd  indicates if the sign of Rlast is to be tested.
      !           If TestRd = true, then
      !
      !                  (>)                       (Pos Def )
      !            Rlast (=) 0 implies that ZHZ is (singular)
      !                  (<)                       (indef   )
      !
      ! On exit,
      !   eigZHZ  encodes the inertia of the reduced Hessian.
      !
      !
      ! 15 Jul 1995: First version of s5ZHZeig.
      ! 02 Aug 2003: snPRNT adopted.
      ! 16 Apr 2005: Signed square root of SC stored in Rlast.
      ! 26 Dec 2014: Fixed bug in the test for a singular R.
      !=================================================================
      character
     &     str*110
      double precision
     &     condH, dRsq, dRsqmin, dRmax, dRmin, condZHZbnd
      ! ------------------------------------------------------------------
      integer            INDEF,     SEMDEF,     POSDEF
      parameter         (INDEF =-1, SEMDEF = 0, POSDEF = 1)
      double precision   zero
      parameter         (zero  = 0.0d+0)
      ! ------------------------------------------------------------------
      condZHZbnd = rw( 85) ! bound on the condition of ZHZ

      if (nS .eq. 0) then
         ! ZHZ is positive definite at a vertex.

         eigZHZ = POSDEF

      else
         ! Compute dRsqmin, the square of the
         ! smallest allowable diagonal of a
         ! positive-definite ZHZ.

         call s6Rcnd
     &      ( maxR, nS-1, lenR, R, dRmax, dRmin, condH )
         dRsqmin = dRmax*(dRmax/condZHZbnd)

         dRsq    = Rlast**2
         if (TestRd  .and.  Rlast .lt. zero) dRsq = -dRsq

         if      (dRsq .ge.  dRsqmin) then
            eigZHZ = POSDEF
         else if (dRsq .ge. -dRsqmin) then
            eigZHZ = SEMDEF
         else
            eigZHZ = INDEF
         end if

         if (eigH .eq. SEMDEF  .or.  eigH .eq. POSDEF) then
            ! should be a positive semidefinite ZHZ

            if (eigZHZ .eq. INDEF) then
               write(str, 1000) itn, dRsq, dRsqmin
               call snPRNT( 21, str, iw, leniw )
            end if

            if (eigH .eq. POSDEF  .and.  eigZHZ .eq. SEMDEF) then
               ! should be a positive-definite ZHZ

               write(str, 2000) itn, dRsq, dRsqmin
               call snPRNT( 21, str, iw, leniw )
            end if
         end if
      end if

      return

 1000 format(' Itn', i7, ': Reduced Hessian is indefinite.',
     &                   ' Square of diag, min diag = ', 1p, 2e9.1 )
 2000 format(' Itn', i7, ': Reduced Hessian is semidefinite.',
     &                   ' Square of diag, min diag = ', 1p, 2e9.1 )

      end ! subroutine s5ZHZeig

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5ZHZfac
     &   ( iExit,
     &     eigH, typeFac, itn, lenR, m,
     &     maxR, mBS, nb, nS,
     &     condZHZmax, ZHZmin, ZHZmax, rankZHZ,
     &     hs, kBS, perm,
     &     bl, bu, blBS, buBS, x, xBS, R, E,
     &     iw, leniw )

      implicit
     &     none
      integer
     &     eigH, iExit, itn, leniw, lenR, m, maxR, mBS, nb, nS,
     &     rankZHZ, typeFac, hs(nb), kBS(mBS), perm(maxR), iw(leniw)
      double precision
     &     condZHZmax, ZHZmin, ZHZmax, blBS(mBS), buBS(mBS),
     &     bl(nb), bu(nb), E(maxR), xBS(mBS), x(nb), R(lenR)

!     ==================================================================
!     s5ZHZfac  factors the reduced Hessian Z'HZ.
!
!       On entry,
!         eigH     defines the class of QP Hessian:  positive definite,
!                  semidefinite or indefinite.
!         R(lenR)  holds the upper-triangular part of Z'HZ.
!
!       On exit,
!          R(lenR) holds the factor of the largest positive-definite
!                  subset of the rows and columns of ZHZ.  Excluded
!                  superbasics are made nonbasic at their current value.
!
!         iExit    Result
!         -----    ------
!          -2      H  singular when it should be positive definite
!          -1      H  indefinite when it should be positive definite
!           0      positive-definite factor computed successfully
!
!     13 Oct 1992: First version based on Qpsol routine Qpcrsh.
!     15 Oct 1994: Dependent columns fixed at their current value.
!     02 Aug 2003: snPRNT adopted.
!     01 Nov 2014: ZHZmin and normZ added as arguments.
!     03 Nov 2014: Modified Cholesky added as an option.
!     03 Nov 2014: rankZHZ added to the argument list.
!     10 May 2015: normZ removed from the argument list.
!     19 May 2015: Always factor with interchanges if rank(ZHZ) < nS.
!     ==================================================================
      character
     &     str*90
      logical
     &     NeedChol
      integer
     &     inform, j, jmax, jS, k, kmax, ksave, nmodH, nSsave, pivot
      double precision
     &     dpiv, eps, Hdmin, s
!     ------------------------------------------------------------------
      integer            INDEF,      SEMDEF,     POSDEF
      parameter         (INDEF = -1, SEMDEF = 0, POSDEF  = 1)
      integer            CHOL,       MDCHOL
      parameter         (CHOL  =  0, MDCHOL = 1)
      integer            NoPiv,     Piv
      parameter         (NoPiv = 0, Piv     = 1)
      double precision   one
      parameter         (one   = 1.0d+0)
!     ------------------------------------------------------------------
      iExit   = 0

      eps     = max ( ZHZmax, one )/condZHZmax
      Hdmin   = max ( ZHZmax/condZHZmax, eps )

      if (eigH .eq. POSDEF) then
         Pivot = NoPiv
      else if (eigH .eq. INDEF  .or.  eigH .eq. SEMDEF) then
         Pivot =   Piv
      end if

      NeedChol = .true.

!     -----------------------------------------------------------------
!+    while NeedChol  do
  100 if   (NeedChol) then

         if (     typeFac .eq.   CHOL) then
            call s6chol
     &         ( inform,
     &           Pivot, maxR, nS, lenR, R,
     &           Hdmin, dpiv, rankZHZ, perm )
         else if (typeFac .eq. MDCHOL) then
            call s6mchl
     &         ( inform,
     &           Pivot, maxR, nS, lenR, R,
     &           Hdmin, eps, dpiv, rankZHZ, nmodH, perm, E )
         end if

!        inform > 0 implies rankZHZ < nS.
!        if rankZHZ < nS, then refactor if no interchanges were used.

         if (inform .eq. 0  .or.  Pivot .eq. Piv) then
            NeedChol = .false.
         else
            Pivot    = Piv
         end if

         go to 100
      end if
!+    ------------------------------------------------------------------

      if (Pivot .eq. Piv) then
!        -----------------------
!        Apply any interchanges.
!        -----------------------
         do j = 1, min(rankZHZ,nS)
            jmax = perm(j)
            if (jmax .gt. j) then
               kmax       = m + jmax
               k          = m + j

               ksave      = kBS(kmax)
               kBS(kmax)  = kBS(k)
               kBS(k)     = ksave

               s          = xBS(kmax)
               xBS(kmax)  = xBS(k)
               xBS(k)     = s

               s          = blBS(kmax)
               blBS(kmax) = blBS(k)
               blBS(k)    = s

               s          = buBS(kmax)
               buBS(kmax) = buBS(k)
               buBS(k)    = s
            end if
         end do
      end if

      if (dpiv .lt. Hdmin) then
!        ---------------------------------------
!        Apparently, Z'HZ is not positive definite.
!        rankZHZ < nS
!        ---------------------------------------
         write(str, 9000) itn, dpiv, Hdmin
         call snPRNT( 21, str, iw, leniw )

         if (dpiv .ge. -Hdmin) then ! Z'HZ appears to be singular.
            iExit = -1
            write(str, 9100) itn, nS-rankZHZ, Hdmin
            call snPRNT( 21, str, iw, leniw )
         else                       ! Z'HZ appears to be indefinite.
            iExit = -2
         end if

!        Set hs to match the positive-definite part of Z'HZ.

         nSsave = nS
         do jS = rankZHZ+1, nSsave
            k  = m + jS
            j  = kBS(k)

!           Make variable  j  nonbasic (it is already feasible).
!           hs(j) = -1 means x(j) is strictly between its bounds.

            if      (x(j) .le. bl(j)) then
               x(j)  =  bl(j)
               hs(j) =  0
            else if (x(j) .ge. bu(j)) then
               x(j)  =  bu(j)
               hs(j) =  1
            else
               hs(j) = -1
            end if
            if (bl(j) .eq. bu(j)) hs(j) = 4

            nS = nS - 1
         end do
         nS = min( nS, rankZHZ )
      end if

      return

 9000 format(' Itn', i7, ': Reduced Hessian is not positive definite.',
     &         ' dpiv, Hdmin = ', 1p, e9.2, ',', e9.2 )
 9100 format(' Itn', i7, ': Reduced Hessian appears to have ',
     &         i6, ' small eigenvalues.  PD tol = ', 1p, e9.2 )

      end ! subroutine s5ZHZfac

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Zp
     &   ( iExit,
     &     m, mBS, n, nb, nS, eps0, pNorm,
     &     neA, nlocA, locA, indA, Acol,
     &     kBS, pBS, y, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iExit, leniw, lenrw, m, mBS, n, nb, neA, nlocA, nS,
     &     locA(nlocA), indA(neA), kBS(mBS), iw(leniw)
      double precision
     &     eps0, pNorm, Acol(neA), pBS(mBS), y(nb), rw(lenrw)

!     ==================================================================
!     s5Zp computes the free components of the search direction
!     p = Z pS, where pS is the search direction for the superbasics,
!     stored in  pBS(m+1:nBS)
!
!     On exit, the  free components of the search direction are stored
!     in pBS(1:nBS). The search direction for the basic variables is
!     stored in pBS(1),...,pBS(m).
!
!     20 Dec 2005: First version of s5Zp.
!     01 Dec 2012: Added pNorm <= 0.0 test.
!     ==================================================================
      external
     &     dnormi
      double precision
     &     dnormi
      integer
     &     nBS
!     ------------------------------------------------------------------
      integer            Normal
      parameter         (Normal = 0)
      integer            WithB
      parameter         (WithB  = 1)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      iExit  = 0
      nBS    = m + nS
      pNorm  = dnormi( nS, pBS(m+1), 1 )

      if (pNorm .le. zero) pNorm = one

      ! First, compute  y = - S*pS and prepare to solve  B*pB = y
      ! for pB, the search direction for the basic variables.
      ! We first normalize y so the LU solver won't ignore
      ! too many "small" elements while computing pB.

      call dscal
     &   ( nS, (one/pNorm), pBS(m+1), 1 )
      call s2Bprod
     &   ( Normal, eps0, n, nS, kBS(m+1),
     &     neA, nlocA, locA, indA, Acol,
     &     (-one), pBS(m+1), nS, zero, y, m )

      ! Solve  B*pBS = y  and unnormalize all of pBS.

      call s2Bsol
     &   ( iExit,
     &     WithB, m, y, pBS, iw, leniw, rw, lenrw  )
      if (iExit .ne. 0) return

      call dscal
     &   ( nBS, pNorm, pBS, 1 )
      pNorm  = dnormi( nBS, pBS, 1 )

      end ! subroutine s5Zp
