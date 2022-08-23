!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn50lp.f
!
!     s5LP     s5BSx    s5degen   s5eGrad  s5eInf  s5eReset
!     s5erc    s5FixS   s5FixX    s5getB   s5hs    s5Inf
!     s5LG     s5LPit   s5price   s5rc     s5setpi s5setx   s5step
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5LP
     &   ( iExit,
     &     probType, probTag, Elastic,
     &     subOptimize, lpLog, NeedLU, Needx,
     &     m, n, nb, nDegen, itLP, itLPmax, itn,
     &     eMode, lvlObjE, printLevel,
     &     minimize, iObj, scaleObj, objAdd,
     &     condZmax, tolOptFP, tolOptLP, tolx,
     &     nInf, sInf, elastics, nInfE, sInfE, wtInf,
     &     piNorm, rgNorm,
     &     neA, nlocA, locA, indA, Acol,
     &     eType, eState, feasType, hs, kBS,
     &     bl, bu, blQP, buQP, blBS, buBS,
     &     gBS, pi, rc, nrhs0, nrhs, rhs, scales,
     &     x, xBS, xFrozen,
     &     iy, iy1, y, y1,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     lpLog
      logical
     &     Elastic, NeedLU, Needx
      integer
     &     probType, iExit, iObj, itLP, itLPmax, itn, lencw, leniw,
     &     lenrw, m, minimize, nrhs0, nrhs, n, nb, neA, nDegen,
     &     nlocA, nInf, elastics, nInfE, eMode, lvlObjE, printLevel,
     &     subOptimize, hs(nb), locA(nlocA), indA(neA), kBS(m+1),
     &     feasType(m+1), eType(nb), eState(nb), iy(nb), iy1(nb),
     &     iw(leniw)
      double precision
     &     condZmax, scaleObj, objAdd, tolOptFP, tolOptLP, tolx,
     &     sInf, sInfE, wtInf, piNorm, rgNorm,
     &     Acol(neA), bl(nb), bu(nb), blQP(nb), buQP(nb),
     &     blBS(m+1), buBS(m+1), gBS(m+1), pi(m), rc(nb), rhs(nrhs0),
     &     scales(nb), x(nb), xBS(m+1), xFrozen(nb), y(nb), y1(nb),
     &     rw(lenrw)
      character
     &     probTag*20, cw(lencw)*8

      !=================================================================
      ! s5LP   solves a linear program.
      !
      ! The optimization can pass through the following phases:
      !
      ! Phase 1               find a feasible point for all variables
      !
      ! Elastic Phase 1       make the nonelastic variables feasible
      !                       while allowing infeasible elastics
      !
      ! Phase 2               minimize the objective
      !
      ! Elastic Phase 2       minimize a composite objective while
      !                       keeping the nonelastics feasible
      !
      !                       In this phase, lvlObjE means:
      !
      !          lvlObjE = 0  zero     weight on the infeasibilities
      !                                (infeasibillities are ignored)
      !                    1  finite   weight on the infeasibilities
      !                    2  infinite weight on the infeasibilities
      !                                (the objective is ignored)
      !
      ! On entry:
      ! ---------
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
      !                  3  FPE  feasible point for equalities only
      !                  4  FPS  feasible point for QP subProblem
      !
      !
      ! On exit:
      ! --------
      !   iExit         Status
      !   -----         ------
      !     -3          Too many iterations
      !     -2          LP is unbounded
      !     -1          Nonelastic variables are infeasible
      !      0          LP solution found
      !     >0          Fatal error
      !
      !  The array kBS is a permutation on the column indices.
      !  kBS(1  :m )  holds the column indices of the basic variables.
      !  Superbasics have been temporarily fixed at their current value.
      !
      ! 30 Sep 1991: First version of s5LP based on Qpsol's lpcore.
      ! 20 Jul 1996: Slacks changed to be the row value.
      ! 06 Aug 1996: First Min Sum version.
      ! 14 Jul 1997: Thread-safe version.
      ! 24 Dec 1999: Suboptimization option added.
      ! 01 Aug 2003: snEXIT and snPRNT adopted.
      ! 24 Dec 2003: pi checked for NaN and Inf entries.
      ! 07 May 2006: s4ksave handles negative values of hs.
      ! 17 Jun 2008: Real workspace reorganized.
      ! 23 Oct 2010: piNorm initialized to one instead of zero.
      ! 26 May 2013: infBnd used to identify infinite bounds.
      ! 20 Sep 2014: bl and bu added for elastic mode.
      ! 17 Nov 2014: Time limit added.
      ! 02 May 2015: Basis no longer refactorized on return to phase 1.
      !=================================================================
      character
     &     str*115
      logical
     &     Checkx, CheckFeas, Checkpi, Feasible, Gotg, GotR, Increase,
     &     FirstFeas, JustPhase1, LPdone, LUok, Needf, Needv, Needpi,
     &     NewB, NewLU, Optimal, Prtlvl10, PrintLog, PrintSum
      integer
     &     inform, itnfix, itnlim, jq, jBq, jBr, jSq, jSr,
     &     kchk, kDegen, kfac, klog, kObj, kp, kPrc, kPrPrt, ksav,
     &     kSumm, leng, lenL0, lenL, lenU0, lenU, LUitn, LUmod, LUsiz0,
     &     LUmax, LUrequest, linesL, linesS, mBS, mnrHdP, mnrHdS,
     &     nBS, nFac, frozen, nonOpt, nParPr, nParPrLP,
     &     nS, nSwap, nnH, neg,
     &     nfix(2), printL, printS, runTime, toldj1, toldj2, toldj3,
     &     typeLU
      double precision
     &     Anorm, Bgrowth, Bold, featol, condZHZ, djq, djqPrt, infBnd,
     &     maxTime, normg, objLP, objPry, pivot, rowError, signObj,
     &     step, tolx0, tolinc, weight, dummy(1)
      double precision
     &     dnormi
!     ------------------------------------------------------------------
      integer            Intern,     Extern
      parameter         (Intern = 0, Extern = 1)
      integer            Check
      parameter         (Check  = 1)
      integer            FP,         FPE,        FPS
      parameter         (FP     = 0, FPE    = 3, FPS    = 4)
      integer            B
      parameter         (B      = 0)
      integer            WithB
      parameter         (WithB  = 1)
      integer            Init,       Optml,      Cycle
      parameter         (Init   = 0, Optml  = 1, Cycle = 2)

      parameter         (nParPr  =  94) ! partial pricing in use
      parameter         (nParPrLP=  99) ! partial pricing for LPs
      parameter         (toldj1  = 184) ! phase 1 dj tol for p.p.
      parameter         (toldj2  = 185) ! phase 2 dj tol for p.p.
      parameter         (toldj3  = 186) ! current optimality tol
      parameter         (lenL0   = 171) ! size of L0
      parameter         (lenU0   = 172) ! size of initial  U
      parameter         (lenL    = 173) ! size of current  L
      parameter         (lenU    = 174) ! size of current  U
      parameter         (kObj    = 205) ! xBS(kObj) is the obj. slack
      parameter         (LUitn   = 215) ! itns since last factorize
      parameter         (LUmod   = 216) ! number of LU mods
      parameter         (printL  = 218) ! (on/off) log     status
      parameter         (printS  = 219) ! (on/off) summary status
      parameter         (linesL  = 220) ! # lines in log     file
      parameter         (linesS  = 221) ! # lines in summary file
      parameter         (mnrHdP  = 223) ! >0 => Minor heading for iPrint
      parameter         (mnrHdS  = 225) ! >0 => Minor heading for iSumm
      parameter         (runTime = 462) ! Solve time

      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one   = 1.0d+0)
!     ------------------------------------------------------------------
      infBnd    = rw( 70) ! definition of an infinite bound
      maxTime   = rw( 79) ! max time allowed

      kchk      = iw( 58) ! check (row) frequency
      kfac      = iw( 59) ! factorization frequency
      ksav      = iw( 60) ! save basis map
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      kDegen    = iw( 63) ! max. expansions of featol
      itnlim    = iw( 89) ! limit on total iterations
      nFac      = iw(210) ! # of LU factorizations

      if (nFac .gt. 0) then
         LUsiz0    = iw(lenL0) + iw(lenU0)
         LUmax     = 2*LUsiz0
      end if

      nnH    = 0                ! local value of nnH
      neg    = 0
      leng   = 1
      nS     = 0                ! local value of nS
      nBS    = m
      mBS    = m + 1            ! simplex steps

      Prtlvl10 =        printLevel .ge. 10
      PrintLog =        printLevel .ge.  1  .and.
     &        ((mod( itLP,  klog ) .eq.  0  .and.  itLP   .ne. 0)  .or.
     &                      klog   .eq.  1                      )
      PrintSum =        printLevel .ge.  1  .and.
     &        ((mod( itLP,  kSumm) .eq.  0  .and.  itLP   .ne. 0)  .or.
     &                      kSumm  .eq.  1                      )
      iw(printL) = 0
      iw(printS) = 0
      if (PrintLog) iw(printL) = 1
      if (PrintSum) iw(printS) = 1

      kPrc       = 0            ! last section scanned in part. pricing
      iExit      = 0
      LUrequest  = 0

      CheckFeas  = .true.       ! Check that x is feasible.
      Checkpi    = .true.       ! Compute norm pi in s5setpi
      Feasible   = .false.
      GotR       = .false.
      FirstFeas  = .false.

      !-----------------------------------------------------------------
      ! s5LP operates in either ``Normal'' or ``Elastic'' mode.
      ! Everything is normal unless a weighted sum is being minimized or
      ! the constraints are infeasible.
      ! The logical Feasible refers to the nonelastic variables.
      ! wtInf  is the optional parameter Infeasibility Weight.
      !-----------------------------------------------------------------
      ! JustPhase1 = stop after phase1 (either normal or elastic)

      JustPhase1 = probType .eq. FP   .or.
     &             probType .eq. FPE  .or.
     &             probType .eq. FPS

      ! The phase 2 objective is F1 + wtInf*F2.

      if (Elastic) then
         Needf = lvlObjE .ne. 2 ! F1 required in phase 2
         Needv = lvlObjE .ne. 0 ! F2 required in phase 2
      else
         Needf = .true.
         Needv = .false.
      end if

      Needpi  = .true.
      NewLU   = .true.
      Optimal = .false.
      LPdone  = .false.

      condZHZ =  zero
      objLP   =  zero
      pivot   =  zero
      rgnorm  =  zero
      step    =  zero

      sInfE   =  zero
      nInfE   =  0
      nonOpt  = -1

      jq      =  0
      djq     =  zero
      jBq     =  0              ! x(jBq) is the incoming  BS
      jBr     =  0              ! x(jBr) is the outgoing  BS
      jSq     =  0              ! x(jSq) is the incoming SBS
      jSr     =  0              ! x(jSr) is the outgoing SBS
      kPrPrt  =  0
      signObj =  minimize
      typeLU  =  B

      frozen  =  0

      dummy(1)   = zero
      rw(toldj1) = 100.0d+0     ! Used only for LP partial pricing
      rw(toldj2) = 100.0d+0     !

      iw(nParPr) = iw(nParPrLP)


      call s5hs
     &   ( Intern, nb, blQP, buQP, hs, x )
      call s5degen
     &   ( inform, Init, printLevel, nb, nInf, itn,
     &     featol, tolx, tolinc, hs, blQP, buQP, x,
     &     itnfix, nfix, tolx0,
     &     iw, leniw, rw, lenrw )

!!    ======================Start of main loop==========================
!+    do while (.not. LPdone  .and.  iExit .eq. 0)
  100 if       (.not. LPdone  .and.  iExit .eq. 0) then
         !==============================================================
         ! Check the initial  x  and move it onto  ( A -I )*x = b.
         ! If NeedLU is true, this will require a basis factorization.
         !==============================================================
         ! If necessary,  factorize the basis  ( B = LU ) and set x.
         ! If NeedLU is false on entry to s5LP, the first call to s2Bfac
         ! will try to use existing factors.
         ! If NeedLU is true on entry to s5LP, an LU factorization of
         ! type typeLU is computed.
         !
         ! LUrequest =  1  Frequency
         ! LUrequest =  2  LU nonzeros increased
         ! LUrequest =  3
         ! LUrequest =  4
         ! LUrequest =  5  Singular after LU mod
         ! LUrequest =  6  Unstable LU mod (growth in new column of U)
         ! LUrequest =  7  Not enough memory
         ! LUrequest =  8
         ! LUrequest =  9
         ! LUrequest = 10  Row error in setx
         ! LUrequest = 11  Big  dx   in setx or setpi
         ! LUrequest = 23  Infeasibility after refactorization
         !
         ! In elastic mode, reset eState, blQP, buQP, blBS and buBS for
         ! the elastics. elastics gives the number of elastic variables.
         !--------------------------------------------------------------
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
               if (Prtlvl10) iw(mnrHdP) = 1
            end if

            if (iExit .gt. 0) go to 100

            Needpi    = .true.     ! Recalculate the pi's.
            Needx     = .false.
            Checkpi   = .true.
            CheckFeas = .true.
         end if

         Optimal = .false.
         nBS     = m + nS

         nInf    = 0
         sInf    = zero

         call dload ( nBS, zero, gBS, 1 )
         normg   = one

         if (CheckFeas) then

            ! In Phase 1 or just after a factorize, check that the basic
            ! and superbasic nonelastics are feasible.
            ! Any elastic variables are checked first

            if (Elastic) then

!              Check that the elastic variables satisfy blQP and buQP.
!              Change eState if necessary.

               call s5eReset
     &            ( nBS, nb, elastics, featol, infBnd,
     &              eType, eState, kBS,
     &              bl, bu, blQP, buQP, blBS, buBS, x )
            end if

            ! FirstFeas = true means that we have just become feasible.
            ! FirstFeas is turned off once a step is taken.

            call s5Inf
     &         ( nBS, featol, infBnd,
     &           nInf, sInf, feasType, blBS, buBS, gBS, xBS )

            if (nInf .gt. 0) then

               ! The nonelastics are infeasible.  If necessary, switch
               ! back to the feasibility phase,  after refactorization.
               ! Print something if the basis has been refactorized.

               if (Prtlvl10  .and.  iw(LUitn) .eq. 0) then
                  write(str, 1030) itn, nInf, sInf
                  call snPRNT( 21, str, iw, leniw )
               end if
               Feasible = .false.
            end if

!           Feasible => the nonelastics are feasible.
!           Feasible => normal or elastic Phase 2

            if (.not. Feasible) then
               FirstFeas = nInf .eq. 0 ! Feasible for the first time
            end if
            Feasible  = nInf .eq. 0
            CheckFeas = nInf .gt. 0
         end if ! CheckFeas

         if (Elastic) then
            !-----------------------------------------------------------
            ! Find the sum of infeasibilities of the elastic variables.
            ! If  nInfE .ne. elastics, then there is at least one
            ! nonbasic elastic fixed at its current value.
            !-----------------------------------------------------------
            call s5eInf
     &         ( nb, nBS, eState, kBS, featol, nInfE, sInfE, bl, bu, x )
         end if

         objLP  = zero
         if (iObj .gt. 0) then
            objLP  = xBS(iw(kObj))*scaleObj
         end if

         if (Feasible  .and.  JustPhase1) then
            ! The non-elastics are feasible.  Prepare to exit.
            djqPrt = zero
            call dload ( m, zero, pi, 1 )
            piNorm  = one       ! piNorm = max(norm(pi), 1.0)
         else

            if (Feasible) then
               !--------------------------------------------------------
               ! Feasible for the nonelastics.
               ! (Elastic = false means no elastics.)
               !--------------------------------------------------------
               if (Needf) then
                  if (iObj .ne. 0) then
                     gBS(iw(kObj)) = signObj*scaleObj
                  end if
               end if

               if (Elastic  .and.  elastics .gt. 0  .and.  Needv) then
                  call s5eGrad
     &               ( nb, nBS, wtInf, eState, kBS, gBS )
               end if
               normg = max( dnormi( nBS, gBS, 1 ), one )
            end if

            if (Needpi) then
               !--------------------------------------------------------
               ! Compute pi, the multipliers for Ax - s = b.
               !--------------------------------------------------------
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
                     if (.not. LUok) iExit = 43
                  end if
                  go to 100
               end if
               Needpi = .false.
            end if

            !===========================================================
            ! Check for optimality.
            ! Find the reduced costs.
            !===========================================================
            if (Feasible) then
               rw(toldj3) = tolOptLP
            else
               rw(toldj3) = tolOptFP
            end if

            kPrPrt = kPrc
            jq     = 0
            djqPrt = djq
            djq    = zero
            Gotg   = .false.
            weight = zero
            if (Elastic  .and.  Feasible) then
               weight = wtInf
            end if

            call s5price
     &         ( Elastic, Feasible, Increase, Gotg, subOptimize,
     &           itn, m, n, nb, leng, neg,
     &           frozen, nonOpt, weight, signObj, piNorm,
     &           jq, djq, kPrc, rw(toldj1),
     &           neA, nlocA, locA, indA, Acol,
     &           eType, hs, dummy, pi, rc, x, xFrozen,
     &           iw, leniw, rw, lenrw )
            Optimal = nonOpt .eq. 0
         end if ! Feasible and JustPhase1

         LPdone = Optimal .or. (Feasible .and. JustPhase1)

         if (LPdone) then
            !-----------------------------------------------------------
            ! Apparently we are optimal.
            ! See if any nonbasics have to be set back on their bounds.
            !-----------------------------------------------------------
            call s5degen
     &         ( inform, Optml, printLevel,
     &           nb, nInf, itn,
     &           featol, tolx, tolinc, hs, blQP, buQP, x,
     &           itnfix, nfix, tolx0,
     &           iw, leniw, rw, lenrw )

            LPdone = inform .eq. 0

            if (LPdone) then
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
                  LPdone    = inform .eq. 0
                  LUrequest = inform
               end if
            end if

            ! If x is not optimal, set  x  so that ( A  -I )*x = b
            ! and check feasibility.

            if (.not. LPdone) then
               Needx  = .true.
               go to 100
            end if
         end if ! LPdone

!        ============================================================
!        Print the details of this iteration.
!        ============================================================
         if (FirstFeas .and. JustPhase1) then
            ! Relax, we are about to exit without printing
         else
            objPry = zero
            if (Feasible) then
               if (Needf) then
                  objPry = objAdd + objLP
               end if
               if (Needv) then
                  objPry = objPry + signObj*wtInf*sInfE
               end if
            end if

            call lpLog
     &         ( probType, probTag,
     &           Elastic, GotR, FirstFeas, Feasible, JustPhase1,
     &           m, mBS, nnH, nS, jSq, jBr, jSr,
     &           iw(linesL), iw(linesS), itn, itLP, kPrPrt, lvlObjE,
     &           pivot, step, nInf, sInf, nInfE, sInfE, wtInf,
     &           nonOpt, objPry, condZHZ, djqPrt, rgNorm, kBS, xBS,
     &           iw, leniw )
         end if

         jBq    = 0
         jBr    = 0
         jSq    = 0
         jSr    = 0
         kPrPrt = 0

         if (LPdone) then
            !-----------------------------------------------------------
            ! Optimal  .or. (Feasible .and. JustPhase1)
            !-----------------------------------------------------------
            if (nInf .gt. 0) then

               ! No feasible point.
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
                  LPdone    = .false.

                  Needf     = lvlObjE .ne. 2 ! Use F1 in e-phase 2
                  Needv     = lvlObjE .ne. 0 ! Use F2 in e-phase 2
                  Needpi    = .true.
                  djq       = zero
                  step      = zero
               end if
               go to 100
            end if

            if (Prtlvl10 .and. .not. JustPhase1) then
               If (jq .ne. 0) then
                  djq = signObj*djq
                  if (klog .eq. 1) then
                     write(str, 1010) djq, jq, rgnorm, piNorm
                     call snPRNT( 31, str, iw, leniw )
                  end if
               else
                  if (klog .eq. 1) then
                     write(str, 1020)          rgnorm, piNorm
                     call snPRNT( 31, str, iw, leniw )
                  end if
               end if
            end if
         else
            !-----------------------------------------------------------
            ! Do another LP iteration.
            ! A nonbasic has been selected to become superbasic.
            ! Compute the vector y such that B y = column jq.
            !-----------------------------------------------------------
            ! Unpack column jq into  y1  and solve  B*y = y1.
            ! The solve computes  y1  such that  L*y1 = ajq.
            ! It is used below to modify L and U in s5LPit.

            call s2unpack
     &         ( jq, m, n, neA, Anorm, nlocA, locA, indA, Acol, y1 )
            call s2Bsol
     &         ( iExit, WithB, m, y1, y, iw, leniw, rw, lenrw )
            if (iExit .gt. 0) return

            !===========================================================
            ! Take a simplex step.  A variable will become nonbasic
            ! at the new x.
            !===========================================================
            if (itn  .ge. itnlim  .or.  itLP .ge. itLPmax) then
               iExit = -3       ! Excess iterations
               go to 100
            end if

            itLP      = itLP       + 1
            itn       = itn        + 1
            iw(LUitn) = iw(LUitn)  + 1
            NewLU     = .false.

            ! Decide if we want to print something this iteration.

            PrintLog = printLevel .ge. 1 .and. mod( itLP, klog  ) .eq. 0
            PrintSum = printLevel .ge. 1 .and. mod( itLP, kSumm ) .eq. 0

            iw(printL) = 0
            iw(printS) = 0
            if (PrintLog) iw(printL) = 1
            if (PrintSum) iw(printS) = 1

            !-----------------------------------------------------------
            ! Take a simplex step.
            ! The new x will still be at a vertex (possibly temporary)
            ! Check for unboundedness (inform = 1).
            !-----------------------------------------------------------
            call s5LPitn
     &         ( inform,
     &           Feasible, Increase, Needpi, Elastic,
     &           m+1, m, nb, nDegen, elastics, LUrequest,
     &           kp, jBq, jSq, jBr, jSr, jq,
     &           featol, pivot, step, tolinc,
     &           eType, eState, feasType, hs, kBS,
     &           bl, bu, blQP, buQP, blBS, buBS,
     &           x, xBS, y, y1,
     &           iw, leniw, rw, lenrw )
            if (inform .eq. 1) then
               iExit = -2       ! Unbounded direction
               go to 100
            end if

            ! Increment featol every iteration.

            featol = featol + tolinc

            !===========================================================
            ! Test for error condition and/or frequency interrupts.
            !===========================================================
            !(1) Save a basis map (frequency controlled).
            !(2) Every kdegen iterations, reset featol and move nonbasic
            !    variables onto their bounds if they are very close.
            !(3) Refactorize the basis if it has been modified
            !    too many times.
            !(4) Update the LU factors of the basis if requested.
            !(5) Check row error (frequency controlled).

            if (mod(itn,ksav) .eq. 0) then
               call s4ksave
     &            ( minimize, m, n, nb, nS, mBS,
     &              itn, nInf, sInf, objLP, kBS, hs,
     &              scales, blQP, buQP, x, xBS, cw, lencw, iw, leniw )
            end if

            if (mod( itn, kdegen ) .eq. 0) then
               call s5degen
     &            ( inform, Cycle, printLevel,
     &              nb, nInf, itn,
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
                  Checkx    = mod(iw(LUitn),kchk) .eq. 0
                  if (Checkx  .and.  .not. Needx) then
                     call s5setx
     &                  ( inform, Check, itn,
     &                    m, n, nb, nBS, rowError,
     &                    neA, nlocA, locA, indA, Acol,
     &                    kBS, xBS, nrhs0, nrhs, rhs, x, y, y1,
     &                    iw, leniw, rw, lenrw )
                     LUrequest = inform
                  end if
               end if
               if (LUrequest .gt. 0) typeLU = B
            end if

            ! Check time limit

            if (maxTime .gt. zero .and. mod(itn, 50) .eq. 0) then
               call s1time ( -2, 0, iw, leniw, rw, lenrw )
               call s1time (  2, 0, iw, leniw, rw, lenrw )
               if (rw(runTime) .gt. maxTime) then
                  iExit = 34    ! time limit reached
               end if
            end if
         end if ! not Optimal

         go to 100
!+    end while
      end if
!     ======================end of main loop============================
!
      call s5hs
     &   ( Extern, nb, blQP, buQP, hs, x )

      return

 1000 format(' ==> LU file has increased by a factor of', f6.1)
 1010 format(' Biggest dj =', 1p, e11.3, ' (variable', i7, ')',
     &       '    norm rg =',     e11.3, '   norm pi =', e11.3)
 1020 format(   ' Norm rg =', 1p, e11.3, '   norm pi =', e11.3)
 1030 Format(' Itn', i7, ': Infeasible nonelastics.  Num =', i7, 1p,
     &       '   Sum of Infeasibilities =', e8.1 )

 8050 format(' Itn', i7, ': Infeasible ', a)
 8060 format(' Itn', i7, ': Elastic Phase 1 -- making',
     &                   ' nonelastic variables feasible')

      end ! subroutine s5LP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5BSx
     &   ( Task, nBS, nb, kBS, x, xBS )

      implicit
     &     none
      integer
     &     Task, nBS, nb, kBS(nBS)
      double precision
     &     xBS(nBS), x(nb)

!     =================================================================
!     s5BSx   copies free variables from  xBS  into  x  or vice versa,
!             depending on whether  Task is 'xBS to x' or 'x to xBS'.
!
!     07 Nov 1991: First version based on Minos routine m5bsx.
!     21 Aug 1999: Current version of s5BSx.
!     =================================================================
      integer
     &     j, k
!     ------------------------------------------------------------------
      integer            xtoxBS,     xBStox
      parameter         (xtoxBS = 0, xBStox = 1)
!     ------------------------------------------------------------------

      if (Task .eq. xBStox) then
         do k = 1, nBS
            j     = kBS(k)
            x(j)  = xBS(k)
         end do

      else if (Task .eq. xtoxBS) then
         do k = 1, nBS
            j      = kBS(k)
            xBS(k) = x(j)
         end do
      end if

      end ! subroutine s5BSx

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5degen
     &   ( inform, Task, printLevel, nb, nInf, itn,
     &     featol, tolx, tolinc, hs, bl, bu, x,
     &     itnfix, nfix, tolx0,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Task, inform, itnfix, nb, nInf, itn, leniw, lenrw,
     &     printLevel, hs(nb), nfix(2), iw(leniw)
      double precision
     &     featol, tolx, tolinc, tolx0, bl(nb), bu(nb), x(nb), rw(lenrw)

!     ==================================================================
!     s5degen performs most of the tasks associated with degeneracy.
!     The degeneracy-resolving strategy operates in the following way.
!
!     Over a cycle of iterations, the feasibility tolerance featol
!     increases slightly (from tolx0 to tolx1 in steps of tolinc).
!     This ensures that all steps taken will be positive.
!
!     After kdegen consecutive iterations, nonbasic variables within
!     featol of their bounds are set exactly on their bounds and the
!     basic variables are recomputed to satisfy ( A  -I )*x = b.
!     featol is then reduced to tolx0 for the next cycle of iterations.
!
!
!     If Task = Init, s5degen initializes the parameters:
!
!     featol  is the current feasibility tolerance.
!     tolx0   is the minimum feasibility tolerance.
!     tolx1   is the maximum feasibility tolerance.
!     tolinc  is the increment to featol.
!     kdegen  is the expand frequency (specified by the user).
!             it is the frequency of resetting featol to tolx0.
!     nDegen  counts the number of degenerate steps (not used here, but
!             incremented by s5step).
!     itnfix  is the last iteration at which an 'Optimal' or 'Cycle'
!             entry set nonbasics onto their bound.
!     nfix(j) counts the number of times an 'Optimal' entry has
!             set nonbasics onto their bound,
!             where j=1 if infeasible, j=2 if feasible.
!
!     tolx0 and tolx1 are both close to the feasibility tolerance tolx
!     specified by the user.  (They must both be less than tolx.)
!
!
!     If Task = Cycle,  s5degen has been called after a cycle of
!     kdegen iterations.  Nonbasic x(j)s are examined to see if any are
!     off their bounds by an amount approaching featol.  inform returns
!     how many.  Deviations as small as tolz (e.g. 1.0d-11) are not
!     counted. If inform is positive, the basic variables are
!     recomputed.  It is assumed that s5LP or s5QP will then continue
!     iterations.
!
!     itnfix, nfix, tolx0 could be treated as SAVED variables.
!     They are included as arguments to prevent threading problems in a
!     multiprocessing environment.
!
!     If Task = Optml,  s5degen is being called after a subproblem
!     has been judged optimal, infeasible or unbounded.
!     Nonbasic x(j)s are examined as above.
!
!     07 Nov 1991: First version based on Minos routine m5dgen.
!     12 Jul 1997: Thread-safe version.
!     11 May 2001: Removed summary file printing.
!     01 Aug 2003: snPRNT adopted.
!     01 Aug 2003: Current version of s5degen.
!     ==================================================================
      character
     &     str*80
      integer
     &     j, kDegen, maxfix
      double precision
     &     b1, b2, d1, d2, eps1, tolx1, tolz
!     ------------------------------------------------------------------
      integer            Init,     Optml,     Cycle
      parameter         (Init = 0, Optml = 1, Cycle = 2)
      double precision   zero
      parameter         (zero = 0.0d+0)
!     ------------------------------------------------------------------
      eps1      = rw(  3) ! eps**(2/3)
      kDegen    = iw( 63) ! max. expansions of featol

      inform    = 0
      if (Task .eq. Init) then

!        Task = Initialize.
!        Initialize at the start of each major iteration.
!        kdegen is the expand frequency      and
!        tolx   is the feasibility tolerance
!        (specified by the user).  They are not changed.
!        nDegen counts the total number of degenerate steps, and is
!        initialized by s5solve or s8solve.

         itnfix  = 0
         nfix(1) = 0
         nfix(2) = 0
         tolx0   = 0.5d+0 *tolx
         tolx1   = 0.99d+0*tolx
         if (kdegen .lt. 99999999) then
            tolinc = (tolx1 - tolx0) / kdegen
         else
            tolinc = zero
         end if
         featol  = tolx0
      else if (Task .eq. Cycle  .or.  Task .eq. Optml) then
!        ---------------------------------------------------------------
!        Task = 'E'nd of cycle or 'O'ptimal.
!        initialize local variables maxfix and tolz.
!        ---------------------------------------------------------------
         maxfix = 2
         tolz   = eps1
         if (Task .eq. Optml) then

!           Task = Optimal.
!           Return with inform = 0 if the last call was at the
!           same itn, or if there have already been maxfix calls
!           with the same state of feasibility.

            if (itnfix .eq. itn   ) return
            if (nInf   .gt.   0   ) then
               j = 1
            else
               j = 2
            end if
            if (nfix(j) .ge. maxfix) return
            nfix(j) = nfix(j) + 1
         end if

!        Set nonbasics on their nearest bound if they are within
!        the current featol of that bound.

         itnfix = itn

         do j = 1, nb
            if (hs(j) .le. 1  .or.  hs(j) .eq. 4) then
               b1    = bl(j)
               b2    = bu(j)
               d1    = abs( x(j) - b1 )
               d2    = abs( x(j) - b2 )
               if (d1 .gt. d2) then
                  b1   = b2
                  d1   = d2
               end if
               if (d1 .le. featol) then
                  if (d1 .gt. tolz) inform = inform + 1
                  x(j) = b1
               end if
            end if
         end do

!        Reset featol to its minimum value.

         featol = tolx0
         if (inform .gt. 0) then

!           The basic variables will be reset.

            if (printLevel .ge. 10) then
               write(str, 1000) itn, inform
               call snPRNT( 21, str, iw, leniw )
            end if
         end if
      end if

      return

 1000 format(' Itn', i7, ': Basics recomputed after ', i7,
     &       '  nonbasics set on bound')

      end ! subroutine s5degen

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5eGrad
     &   ( nb, nBS, wtInf, eState, kBS, gBS )

      implicit
     &     none
      integer
     &     nb, nBS, eState(nb), kBS(nBS)
      double precision
     &     wtInf, gBS(nBS)

!     ==================================================================
!     s5eGrad  is called when elastic variables are allowed to violate
!     their true bounds bl and bu.  It updates the gradient gBS
!     to include the gradient of the elastic penalty term.
!
!     On exit,
!       gBS(nBS)    is the rhs for the equations for pi.
!
!     08 Oct 1996: First version of s5Egrd.
!     21 Apr 1999: Current version.
!     ==================================================================
      integer
     &     j, eStatej, k
!     ------------------------------------------------------------------
      do k       = 1, nBS
         j       = kBS(k)
         eStatej = eState(j)
         if      (eStatej .eq. 0) then
            ! Relax
         else if (eStatej .eq. 1) then
            gBS(k) = gBS(k) - wtInf
         else !  (eStatej .eq. 2)
            gBS(k) = gBS(k) + wtInf
         end if
      end do

      end ! subroutine s5eGrad

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5eInf
     &   ( nb, nBS, eState, kBS, featol, nInfE, sInfE,
     &     bl, bu, x )

      implicit
     &     none
      integer
     &     nb, nBS, nInfE, eState(nb), kBS(nBS)
      double precision
     &     featol, sInfE, bl(nb), bu(nb), x(nb)

!     ==================================================================
!     s5eInf  computes the sum and number of elastic infeasibilities.
!
!     On exit,
!     nInfE is the number of elastics allowed to go infeasible.
!
!     20 Aug 1996: First version of s5eInf.
!     15 Oct 2014: Current version.
!     ==================================================================
      integer
     &     j, eStatej, k, numInf
      double precision
     &     res, sumInf
!     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
      !-----------------------------------------------------------------
      ! Find the number and sum of the elastic infeasibilities.
      !-----------------------------------------------------------------
      numInf = 0
      sumInf = zero

      do k       = 1, nBS
         j       = kBS(k)
         eStatej = eState(j)

         if      (eStatej .eq. 0) then

            ! Relax, this variable is not elastic

         else if (eStatej .eq. 1) then
            numInf = numInf + 1
            res    = bl(j) - x(j)

            if (res .gt. featol) then
               sumInf = sumInf + res
            end if
         else if (eStatej .eq. 2) then
            numInf = numInf + 1
            res    = x(j) - bu(j)

            if (res .gt. featol) then
               sumInf = sumInf + res
            end if
         end if
      end do

      sInfE = sumInf
      nInfE = numInf

      end ! subroutine s5eInf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5eReset
     &   ( nBS, nb, elastics, featol, infBnd,
     &     eType, eState, kBS,
     &     bl, bu, blQP, buQP, blBS, buBS, x )

      implicit
     &     none
      integer          ! In
     &     nBS, nb, eType(nb), kBS(nBS)
      integer          ! Out
     &     elastics
      integer          ! InOut
     &     eState(nb)
      double precision ! In
     &     infBnd, featol, bl(nb), bu(nb), x(nb)
      double precision ! InOut
     &     blQP(nb), buQP(nb), blBS(nBS), buBS(nBS)

!     ==================================================================
!     s5eReset is called after a factorize in elastic mode.
!
!     The upper and lower bounds on the elastic variables are reset for
!     Elastic Phase 1 and 2.
!
!     x = x(feasible) + w - v
!
!     s5eReset is called in elastic mode after a factorize.
!     The bounds blBS and buBS are redefined for the basic elastics.
!
!     In the case of a singular basis, variables associated with
!     dependent basic columns are fixed inside their bounds. If one of
!     these basics is elastic, it will be fixed inside its working
!     bounds but outside bl or bu.  These variables are not counted in
!     elastic infeasibilities because they will be made basic in
!     s5price.
!
!     23 Aug 1996: First version of e5Reset (s5Eset).
!     12 Jun 2000: Elastic mode cleaned up.
!     15 Oct 2014: blQP and buQP added to handle fixed elastics.
!     ==================================================================
      integer
     &     j, eTypej, eStatej, k
      double precision
     &     blj, buj, violL, violU, xj
!     ------------------------------------------------------------------
      elastics   = 0

      do k       = 1, nBS
         j       = kBS(k)
         eTypej  = eType(j)
         eStatej = eState(j)

         blj     = bl(j)
         buj     = bu(j)
         xj      =  x(j)

         if (eStatej .gt. 0) then
!           ------------------------------------------------------------
!           x(j) is elastic.
!
!           EXPAND may give a slightly feasible elastic variable, i.e.,
!               x(j) le bl(j) + featol  or x(j) ge bu(j) - featol.
!           x(j) is kept elastic regardless.
!           ------------------------------------------------------------
            if (eStatej .eq. 1) then

!              eStatej predicts that xj violates its true lower bound,
!              i.e., eTypej must be either 1 or 3.
!
!              In this case  x(j) = bl(j) - v(j), with  v(j) ge -featol

               violL = blj - xj

               if (violL .gt. -featol) then

!                 x(j) violates its true lower bound as predicted.
!                 Keep eState and the QP bounds as they are.

                  elastics = elastics + 1
               else

!                 x(j) satisfies its true lower bound.
!                 Reset eState and the QP bounds to their true values.

                  eStatej = 0
                  blQP(j) = blj
                  buQP(j) = buj

!                 x(j) is made nonelastic unless it violates its
!                 (elastic) bu(j) by more than +featol.

                  if (eTypej .eq. 3) then
!
!                    x(j) = bu(j) + w(j),  with  w(j) ge +featol

                     violU = xj - buj
                     if (violU .gt. +featol) then
                        elastics =  elastics + 1
                        eStatej  =  2
                        blQP(j)  =  buj
                        buQP(j)  =  infBnd
                     end if
                  end if
                  blBS(k) = blQP(j)
                  buBS(k) = buQP(j)
               end if

            else if (eStatej .eq. 2) then

!              x(j) is predicted to violate its true upper bound.
!              eTypej must be either 2 or 3.
!              x(j) = bu(j) + w(j),  with  w(j) ge -featol

               violU = xj - buj

               if (violU .gt. -featol) then

!                 x(j) violates its true upper bound as predicted.
!                 Keep eState and the QP bounds as they are.

                  elastics = elastics + 1
               else

!                 x(j) satisfies its true upper bound.
!                 Reset the QP bounds to their true values.

                  eStatej = 0
                  blQP(j) = blj
                  buQP(j) = buj

!                 x(j) is made nonelastic unless it violates its
!                 (elastic) bl(j) by more than +featol.

                  if (eTypej .eq. 3) then
                     violL = blj - xj

                     if (violL .gt. +featol) then
                        elastics =  elastics + 1
                        eStatej  =  1
                        blQP(j)  = -infBnd
                        buQP(j)  =  blj
                     end if
                  end if
                  blBS(k) = blQP(j)
                  buBS(k) = buQP(j)
               end if
            end if

         else ! eStatej = 0
!           -----------------------------------------------------------
!           Check if any basic or superbasic elastic variables that are
!           predicted to be feasible violate their elastic bounds.
!           -----------------------------------------------------------
            if (eTypej .eq. 1) then

!              If xj violates its lower bound by more than featol,
!              make it elastic.

               violL = blj - xj
               if (violL .gt. +featol) then
                  elastics =  elastics + 1
                  eStatej  =  1
                  blQP(j)  = -infBnd
                  buQP(j)  =  blj
                  blBS(k)  =  blQP(j)
                  buBS(k)  =  buQP(j)
               end if

            else if (eTypej .eq. 2) then
               violU = xj - buj
               if (violU .gt. featol) then
                  elastics =  elastics + 1
                  eStatej  =  2
                  blQP(j)  =  buj
                  buQP(j)  =  infBnd
                  blBS(k)  =  blQP(j)
                  buBS(k)  =  buQP(j)
               end if
            else if (eTypej .eq. 3) then
               violL  = blj  - xj
               violU  = xj   - buj

               if (violL .gt. featol) then
                  elastics =  elastics + 1
                  eStatej  =  1
                  blQP(j)  = -infBnd
                  buQP(j)  =  blj
                  blBS(k)  =  blQP(j)
                  buBS(k)  =  buQP(j)
               else if (violU .gt. featol) then
                  elastics =  elastics + 1
                  eStatej  =  2
                  blQP(j)  =  buj
                  buQP(j)  =  infBnd
                  blBS(k)  =  blQP(j)
                  buBS(k)  =  buQP(j)
               end if
            end if
         end if

!        Mark the basic and superbasic elastics.

         eState(j) = eStatej

      end do

      end ! subroutine s5eReset

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5erc
     &   ( j1, j2, Gotg, m, n, leng, neg, signObj,
     &     neA, nlocA, locA, indA, Acol,
     &     eType, hs, g, pi, rc )

      implicit
     &     none
      logical
     &     Gotg
      integer
     &     j1, j2, m, n, nlocA, leng, neg, neA,
     &     eType(n), hs(n), locA(nlocA), indA(neA)
      double precision
     &     signObj, Acol(neA), g(leng), pi(m), rc(n)

!     ==================================================================
!     s5erc  computes reduced costs rc(j) in the range j = j1 to j2
!     for fixed nonbasic columns that have one or more elastic bounds.
!     It is called by s5price.
!
!     07 Feb 1998: First version based on s5rc.
!     20 Sep 2014: Compute reduced costs for  hs(j) = -1.
!     ==================================================================
      integer
     &     i, j, l
      double precision
     &     dj
!     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
!     ------------------------------------------------------------------
      do j = j1, j2
         if ((hs(j) .eq. 4 .or.      hs(j) .eq. -1)
     &                     .and.  eType(j) .gt.  0) then
            dj    = zero
            do l  = locA(j), locA(j+1)-1
               i  = indA(l)
               dj = dj  +  pi(i) * Acol(l)
            end do
            rc(j) = - dj
         end if
      end do

      ! Include the nonlinear gradient term if present.

      if (Gotg) then
         if (j1 .le. neg) then
            do j = j1, min( j2, neg )
               if ((hs(j) .eq. 4 .or.      hs(j) .eq. -1)
     &                           .and.  eType(j) .gt.  0) then
                  rc(j) = rc(j) + signObj*g(j)
               end if
            end do
         end if
      end if

      end ! subroutine s5erc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5fixS
     &   ( Task, m, maxS, mBS, n, nb, nS, hs, kBS,
     &     bl, bu, blBS, buBS, x, xBS )

      implicit
     &     none
      integer
     &     Task, m, maxS, mBS, n, nb, nS, hs(nb), kBS(mBS)
      double precision
     &     bl(nb), bu(nb), x(nb), blBS(mBS), buBS(mBS), xBS(mBS)

!     ==================================================================
!     s5fixS   concerns temporary bounds on superbasic variables.
!     If Task = Fix,  s5fixS sets hs(j) = -1, 0, 1 or 4 for certain
!     superbasic variables.
!
!     If Task = Free, s5fixS changes -1 values to hs(j) = 2.
!
!     30 May 1995: First version of s5fixS.
!     12 Jul 2001: Current version.
!     ==================================================================
      integer
     &     j, k
!     ------------------------------------------------------------------
      integer            Fix,     Free
      parameter         (Fix = 0, Free = 1)
!     ------------------------------------------------------------------

      if (Task .eq. Fix) then
         if (nS .gt. 0) then
!           ------------------------------------------------------------
!           Change superbasic hs(j) to be temporarily fixed.
!           ------------------------------------------------------------
            nS = 0
            do j = 1, nb
               if (hs(j) .eq. 2) then
                  if (bl(j) .eq. bu(j)) then
                     hs(j) =  4
                  else if (x(j) .le. bl(j)) then
                     hs(j) =  0
                  else if (x(j) .ge. bu(j)) then
                     hs(j) =  1
                  else
                     hs(j) = -1
                  end if
               end if
            end do
         end if

      else if (Task .eq. Free) then
!        ---------------------------------------------------------------
!        Free the temporarily fixed structurals.
!        Load the superbasic variables/bounds into xBS, blBS, buBS.
!        ---------------------------------------------------------------
         j = 1
!+       while (j .le. n  .and.  nS .lt. maxS) do
  100    if    (j .le. n  .and.  nS .lt. maxS) then
            if (hs(j) .eq. -1) then
               nS      = nS + 1
               k       = m  + nS
               hs(j)   = 2
               xBS (k) = x(j)
               blBS(k) = bl(j)
               buBS(k) = bu(j)
               kBS (k) = j
            end if
            j  = j + 1
            go to 100
!+       end while
         end if
      end if

      end ! subroutine s5fixS

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5FixX
     &   ( Task, b1, b2, tolx, hs, bl, bu, x )

      implicit
     &     none
      integer
     &     Task, b1, b2, hs(b2)
      double precision
     &     tolx, bl(b2), bu(b2), x(b2)

!     ==================================================================
!     s5FixX  ensures that variables satisfy their simple bounds.
!
!     If Task = xBound, variables x(b1) through x(b2) are made to
!                       satisfy their bounds.
!     If Task = xMove , variables x(b1) through x(b2) are made to
!                       satisfy their bounds. In addition, any nonbasic
!                       variable close to its bound is moved onto it.
!
!     29 Apr 1999: First version of s5FixX.
!     29 Apr 1999: Current version.
!     ==================================================================
      integer
     &     j
      double precision
     &     xj
!     ------------------------------------------------------------------
      integer            xBound,     xMove
      parameter         (xBound = 0, xMove = 1)
!     ------------------------------------------------------------------
      if (Task .eq. xBound) then
         do j = b1, b2
            x(j) = max( x(j), bl(j) )
            x(j) = min( x(j), bu(j) )
         end do
      else if (Task .eq. xMove) then
         do j  = b1, b2
            xj = max( x(j), bl(j) )
            xj = min( xj  , bu(j) )
            if (hs(j) .le. 1) then
               if (xj .le. bl(j) + tolx) xj = bl(j)
               if (xj .ge. bu(j) - tolx) xj = bu(j)
            end if
            x(j)   = xj
         end do
      end if

      end ! subroutine s5FixX

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5getB
     &   ( iExit, Start, lpLog, NeedB, m, maxS, mBS,
     &     n, nb, nnCon, nnJac, nnObj, nNames, nS,
     &     itQP, itQPmax, itn,
     &     nDegen, numLC, numLIQ, tolOptFP, tolOptQP, tolx,
     &     nInf, sInf, wtInf,
     &     iObj, scaleObj, piNorm, rgNorm,
     &     neA, nlocA, locA, indA, Acol,
     &     eType, eState, RowTypes, feasType,
     &     hs, kBS, Names,
     &     bl, bu, blQP, buQP, blBS, buBS, blSave, buSave,
     &     gBS, pi, rc, nrhs0, nrhs, rhs, scales,
     &     lenx0, nx0, x0, x, xBS,
     &     iy, iy1, y, y1, y2,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     lpLog
      logical
     &     NeedB
      integer
     &     Start, iExit, iObj, itn, itQP, itQPmax, lenx0, nrhs0, nrhs,
     &     nx0, m, maxS, mBS, n, nb, neA, nlocA, nNames, nnCon, nnJac,
     &     nnObj, nDegen, nInf, nS, numLC, numLIQ,
     &     lencw, leniw, lenrw,
     &     eType(nb), eState(nb), feasType(mBS), RowTypes(nb), hs(nb),
     &     kBS(mBS), locA(nlocA), indA(neA), iy(nb), iy1(nb), iw(leniw)
      double precision
     &     tolx, tolOptFP, tolOptQP, sInf, scaleObj, wtInf, piNorm,
     &     rgNorm, Acol(neA), scales(nb), bl(nb), bu(nb),
     &     blQP(nb), buQP(nb), blSave(nb), buSave(nb),
     &     rc(nb), x0(lenx0), x(nb), blBS(mBS), buBS(mBS),
     &     gBS(mBS), xBS(mBS), pi(m), rhs(nrhs0), y(nb),
     &     y1(nb), y2(nb), rw(lenrw)
      character
     &     Names(nNames)*8, cw(lencw)*8

!     ==================================================================
!     s5getB   finds an initial basis kBS(1:m) for the linear
!     constraints and bounds.  First, we attempt to find a  feasible
!     point for the bounds and general linear equalities. This may fail,
!     so there is the chance that the initial basis is not feasible.
!     This difficulty is taken care of in subsequent calls to s5LP from
!     s5solveLP, s5solveQP, s5solveQN or s8feasLC.
!
!     1. The linear constraints are (optionally) scaled.
!
!     2. Elements x(n+1:n+m) of the initial x are assigned.
!        A scaled copy of the initial point is stored in x0.
!
!     3. An LP is used to find a feasible x for the bounds and
!        linear equality constraints.
!
!
!     On entry,
!     ---------
!     Start        is 0,1,2,3: an integer version of the solver's Start
!                  (character for sqopt, snoptb, snoptc, npopt
!                  integer    for snopta).
!
!     RowTypes(nb) is the vector of row types defined in s2Amat.
!
!     bl, bu     contain the user-defined lower and upper bounds
!
!     x(nb)        contains the initial x. Only the first n elements
!                  have been defined (by the user).  Here, elements
!                  x(n+1:n+m) are assigned, (optionally) scaled and
!                  saved in x0.
!
!
!     blSave, buSave are used to save the scaled lower and upper bounds.
!
!     iExit  Status
!     -----  ------
!        0   basis found.
!            nInf = 0 => linear equalities are     satisfied.
!            nInf > 0 +> linear equalities are not satisfied
!                        (unbounded problem with emode > 0)
!
!       12   no feasible point and emode = 0.
!       31   Iterations limit.
!       >0   fatal error
!
!     31 Jul 1996: First version of s5getB.
!     12 Jul 1997: Thread-safe version.
!     01 Aug 2003: snEXIT and snPRNT adopted.
!     08 Mar 2004: Implemented gotHes, gotScl for Hot starts.
!     17 Jun 2008: Real workspace reorganized.
!     20 Nov 2014: Normal exit if infeasible and eMode > 0.
!     ==================================================================
      logical
     &     Elastic, NeedLU, Needx
      character
     &     probTag*20, str*120
      integer
     &     eMode, eModeLEQ, iInsrt, iLoadB, inform, iOldB, iStart, j,
     &     lc1, lCrash, lsSave, lvlObjE, minimize, mjrPrtlvl, mnrPrtlvl,
     &     elastics, nInfE, nrhsLP0, nrhsLP, nnL, numLEQ, lvlScale,
     &     iCrash, iObjFP, subOptimize
      double precision
     &     condZmax0, condZmax, objAdd, infBnd, sInfE, tCrash
!     ------------------------------------------------------------------
      integer            COLD,         BASIS,      WARM,       HOT
      parameter         (COLD     = 0, BASIS  = 1, WARM   = 2, HOT  = 3)
      integer            RowTyp
      parameter         (RowTyp   = 0)
      integer            NoSubOpt
      parameter         (NoSubOpt =-1)
      integer            Scale
      parameter         (Scale    = 0)
      integer            FPE
      parameter         (FPE      = 3)
      integer            Fix
      parameter         (Fix      = 0)
      integer            Intern,       Extern
      parameter         (Intern   = 0, Extern = 1)
      integer            xBound
      parameter         (xBound   = 0)

      double precision   zero,              one
      parameter         (zero     = 0.0d+0, one = 1.0d+0)

      parameter         (lvlScale =  75) ! scale option
      parameter         (iCrash   =  88) ! Crash option
      integer            gotHes,       gotScl
      parameter         (gotHes   = 231, gotScl = 232)
!     ------------------------------------------------------------------
      character          line*4
      data               line /'----'/
!     ------------------------------------------------------------------
      iLoadB    = iw(122) ! load file
      iInsrt    = iw(125) ! insert file
      iOldB     = iw(126) ! old basis file

      eMode     = iw( 56) ! >0   => use elastic mode
      mjrPrtlvl = iw( 92) ! Major print level
      mnrPrtlvl = iw( 93) ! Minor print level

      infBnd    = rw( 70) ! definition of an infinite bound
      tCrash    = rw( 62) ! crash tolerance.
      condZmax0 = rw( 86) ! Initial bound on the condition est of Z

      ! Initialize a few things.

      iExit     = 0
      inform    = 0
      minimize  = 1
      nnL       = max( nnObj, nnJac )
      nInf      = 0
      nInfE     = 0             ! Local to s5getB
      sInfE     = zero          ! Local to s5getB
      objAdd    = zero          ! Local to s5getB
      iObjFP    = 0             ! Used for FP calculation
      scaleObj  = one
      numLEQ    = numLC - numLIQ

      condZmax  = condZmax0

      ! Initialize the row and column scales to "unit scaling"

      if (iw(lvlScale) .gt. 0) then
         call dload ( nb, one, scales, 1 )
      end if

      !=================================================================
      ! Decode Start.
      !=================================================================
      NeedB  = .true.

      if (Start .eq. COLD  .or.  Start .eq. BASIS) then
         !-------------------------------
         ! Cold start  or  Basis file.
         !-------------------------------
         iStart = 0
         NeedB  = max( ioldB, iInsrt, iloadB ) .le. 0
         nS     = 0

      else if (Start .eq. WARM) then
         !-------------------------------
         ! Warm start.
         !-------------------------------
         iStart = 1
         NeedB  = .false.

      else if (Start .eq. HOT ) then
         !-------------------------------
         ! Hot start.
         !-------------------------------
         iStart = 1
         NeedB  = .false.
      end if

      call s1page( 1, iw, leniw )

      if (iStart .eq. 0) then
         !--------------------------------------------------------------
         ! Cold start, or Basis file provided.
         ! Input a basis file if one exists, thereby defining hs and x.
         ! (Otherwise, s2crash will be called later to define hs.)
         !--------------------------------------------------------------
         ! Initialize x(n+1:nb) and pi(1:m) before scaling the problem.
         ! The basis files initialize all of x.
         ! One day they may load pi for nonlinear problems.

         call dload ( m, zero, x(n+1), 1 )
         call dload ( m, zero, pi    , 1 )

         call snPRNT( 21, ' Initial basis', iw, leniw )
         call snPRNT( 21, ' -------------', iw, leniw )

         if (NeedB) then
            call snPRNT( 31, ' No basis file supplied', iw, leniw )

            if (iw(iCrash) .eq. 0) then
               NeedB  = .false.
               lCrash = 0
               call s2crash
     &            ( lCrash, mnrPrtlvl, m, n, nb, nnCon,
     &              iw(iCrash), tCrash,
     &              neA, nlocA, locA, indA, Acol,
     &              kBS, hs, RowTypes, bl, bu, x,
     &              iw, leniw, rw, lenrw )
            end if
         else
            call s4getB
     &         ( iExit, m, n, nb, nNames, nS, iObj,
     &           hs, bl, bu, x, Names, iw, leniw, rw, lenrw )

            if (iExit .gt. 0) then
               ! Set blSave and buSave to avoid uninitialized copying.

               call dcopy ( nb, bl, 1, blSave, 1 )
               call dcopy ( nb, bu, 1, buSave, 1 )
               go to 900
            end if
         end if
      end if ! iStart = 0

      !-----------------------------------------------------------------
      ! Move x inside its bounds.
      !-----------------------------------------------------------------
      call s5FixX
     &   ( xBound, 1, n, tolx, hs, bl, bu, x )

      !-----------------------------------------------------------------
      ! Scale the linear part of the constraints.
      ! (Any nonlinear elements in A contain fake nonzeros.)
      !-----------------------------------------------------------------
      if (iw(lvlScale) .gt. 0  .and.  numLC .gt. 0) then
         if (iw(gotScl) .eq. 0) then
            iw(gotScl)   = 1
            lsSave       = iw(lvlScale)
            iw(lvlScale) = 2
            call s2GetScales
     &         ( mjrPrtlvl, m, n, nb, nnL, nnCon, nnJac, RowTypes,
     &           neA, nlocA, locA, indA, Acol,
     &           scales, bl, bu, y, y1,
     &           iw, leniw, rw, lenrw )
            iw(lvlScale) = lsSave
         end if

         call s2ApplyScales
     &      ( Scale, m, n, nb, iObj, infBnd, scaleObj,
     &        neA, nlocA, locA, indA, Acol,
     &        scales, bl, bu, pi, x )
      end if

      ! Save the (scaled) bounds.  Save the scaled initial point in x0.

      call dcopy ( nb , bl, 1, blSave, 1 )
      call dcopy ( nb , bu, 1, buSave, 1 )
      if (nx0 .gt. 0)
     &call dcopy ( nx0, x  , 1, x0   , 1 )

      !-----------------------------------------------------------------
      ! Prepare to get feasible for the linear constraints.
      ! Relax any nonlinear rows.
      !-----------------------------------------------------------------
      if (nnCon .gt. 0  .and.  numLC .gt. 0) then
         call dload ( nnCon, (-infBnd), bl(n+1), 1 )
         call dload ( nnCon,   infBnd , bu(n+1), 1 )
      end if

      !-----------------------------------------------------------------
      ! Compute a starting basis for the linear constraints.
      ! This means attempting to get feasible for the linear equalities.
      ! If the E rows are infeasible, s5LP  will take care of it.
      !-----------------------------------------------------------------
      if (mnrPrtlvl .gt. 10) then
         write(str, 1332) (line, j=1,28)
         call snPRNT( 31, str, iw, leniw )
         write(str, 1316) (line, j=1,16)
         call snPRNT( 32, str, iw, leniw )
      end if

      !-----------------------------------------------------------------
      ! Define the working bounds.
      !-----------------------------------------------------------------
      call dcopy ( nb, bl, 1, blQP, 1 )
      call dcopy ( nb, bu, 1, buQP, 1 )

      if (NeedB) then
         !--------------------------------------------------------------
         ! Crash is needed to find a basis.
         !--------------------------------------------------------------
         ! Treat Crash 1 the same as Crash 2.

         iw(iCrash)  = max( iw(iCrash), 2 )
         iw(iCrash)  = min( iw(iCrash), 3 )

         if (numLC .gt. 0) then
            !===========================================================
            ! Find a feasible point for the linear constraints.
            !===========================================================
            ! Crash 2 finds a basis for ALL the LINEAR constraints
            !         (equalities and inequalities).
            ! Crash 3 treats linear EQUALITIES separately.
            !-----------------------------------------------------------

            if (iw(iCrash) .eq. 2) then
               ! Relax.

            else if (iw(iCrash) .eq. 3) then
               if (numLEQ .gt. 0) then

                  ! Find a basis for the linear EQUALITIES.
                  ! Linear inequality rows are made to appear to be free.

                  if (mnrPrtlvl .ge. 10) then
                     write(str, 2100) itn
                     call snPRNT( 23, str, iw, leniw )
                  end if
               end if ! numLEQ > 0
            end if ! iCrash = 2 or 3
         end if ! numLC > 0

         ! Call Crash even if there are no linear constraints.
         ! kBS(1:m) is used as workspace.
         ! We haven't done any solve yet, so we should NOT need
         ! call s5hs  ( Extern, nb, bl, bu, hs, x ) .

         lCrash = iw(iCrash)
         call s2crash
     &      ( lCrash, mnrPrtlvl, m, n, nb, nnCon,
     &        iw(iCrash), tCrash,
     &        neA, nlocA, locA, indA, Acol,
     &        kBS, hs, RowTypes, bl, bu, x,
     &        iw, leniw, rw, lenrw )
      end if ! NeedB

      !=================================================================
      ! 1. Set nS to match hs(*).
      ! 2. Set kBS(m+1:m+nS) to define the initial superbasics.
      ! 3. Set all nonbasic x to be within bounds, which may change some
      !    hs values from 0 to 1.
      ! 4. Set nonbasic x to be exactly on nearly satisfied bounds.
      !    (Some nonbasics may still be between bounds.)
      !=================================================================
      call s4checkhs
     &   ( m, maxS, mBS, n, nb, NeedB, iw(gotHes),
     &     nS, iObj, hs, kBS, bl, bu, x,
     &     iw, leniw, rw, lenrw )

      if (NeedB  .and.  numLC .gt. 0) then

         ! Fix SBs. This forces simplex steps.

         call s5fixS
     &      ( Fix, m, maxS, mBS, n, nb, nS, hs, kBS,
     &        bl, bu, blBS, buBS, x, xBS )

         if (numLEQ .gt. 0) then
            if (numLIQ .eq. 0) then
               nrhsLP = 0
            else

               lc1 = n + nnCon + 1

               ! Relax the bounds bl and bu associated with the
               ! linear INEQUALITIES.
               ! The rhs stord in y2  makes the relaxed rows satisfied
               ! at the current x.

               call s5LG
     &            ( m, n, nb, nnCon, nrhsLP,
     &              neA, nlocA, locA, indA, Acol,
     &              bl, bu, nrhs0, nrhs, rhs, y2, x, y,
     &              rw, lenrw )

               ! Initialize the working upper and lower bounds.
               ! bl and bu are unaltered in this call because eMode = 0.

               call dcopy ( numLC, bl(lc1), 1, blQP(lc1), 1 )
               call dcopy ( numLC, bu(lc1), 1, buQP(lc1), 1 )
            end if

            probTag     = 'linear equalities'
            nrhsLP0     =  max(nrhsLP, 1)
            Elastic     = .false.
            lvlObjE     =  0
            eModeLEQ    =  0    ! No Elastic mode for E-row phase 1
            NeedLU      = .true.
            Needx       =  NeedLU
            subOptimize =  NoSubOpt

            ! Set hs = 4, -1 for fixed variables and nonbasics
            ! between their bounds.

            call s5hs
     &         ( Intern, nb, bl, bu, hs, x )
            call s5LP
     &         ( inform, FPE, probTag, Elastic,
     &           subOptimize, lpLog, NeedLU, Needx,
     &           m, n, nb, nDegen, itQP, itQPmax, itn,
     &           eModeLEQ, lvlObjE, mnrPrtlvl,
     &           minimize, iObjFP, scaleObj, objAdd,
     &           condZmax, tolOptFP, tolOptQP, tolx,
     &           nInf, sInf, elastics, nInfE, sInfE, wtInf,
     &           piNorm, rgNorm,
     &           neA, nlocA, locA, indA, Acol,
     &           eType, eState, feasType, hs, kBS,
     &           bl, bu, blQP, buQP, blBS, buBS,
     &           gBS, pi, rc, nrhsLP0, nrhsLP, y2, scales,
     &           x, xBS, x,
     &           iy, iy1, y, y1,
     &           cw, lencw, iw, leniw, rw, lenrw )

            ! Reset the original bounds on the linear inequality rows.

            if (numLIQ .gt. 0) then
               lc1    = n + nnCon + 1
               call dcopy ( numLC, blSave(lc1), 1, bl(lc1), 1 )
               call dcopy ( numLC, buSave(lc1), 1, bu(lc1), 1 )
               call dcopy ( numLC, blSave(lc1), 1, blQP(lc1), 1 )
               call dcopy ( numLC, buSave(lc1), 1, buQP(lc1), 1 )
               call s5hs  ( Intern, nb, bl, bu, hs, x )
            end if ! numLIQ > 0

            if (mnrPrtlvl .ge. 10) then
               if (nInf .eq. 0) then
                  ! The linear E rows are now satisfied.

                  write(str, 2200) itn
                  call snPRNT( 21, ' ', iw, leniw )
                  call snPRNT( 23, str, iw, leniw )
               else
                  ! The linear E rows are infeasible.

                  write(str, 2300) itn
                  call snPRNT( 23, str, iw, leniw )
               end if
            end if
         end if ! numLEQ > 0

         ! Potential inform values from s5LP are:
         !   -3    Too many iterations
         !   -2    variable is unbounded (this should not happen)
         !   -1    linear constraints are infeasible
         !    0    feasible point found
         !   >0    Fatal error

         if (inform .gt. 0) then

            iExit = inform      !  fatal error. Exit
         else if ( inform .eq. -3) then

            iExit = 31          ! Too many iterations
         else if ((inform .eq. -1  .or. inform .eq. -2) .and.
     &                                   eMode .eq.  0) then
            ! Infeasible linear equalities and no elastic mode.

            iExit = 12          ! infeasible linear equalities
         else
            ! A feasible point has been found.
            !
            ! Define a basis that includes the linear INEQUALITIES.
            ! If the linear equalities are infeasible, the feasibility
            ! phase for the inequalities (with an appropriate value of
            ! eMode) will determine what happens next.

            iExit = 0

            if (numLIQ .gt. 0) then
               if (mnrPrtlvl .ge. 10) then
                  if (numLEQ .gt. 0) then
                     write(str, 2400) itn
                     call snPRNT( 23, str, iw, leniw )
                  else
                     write(str, 2410) itn
                     call snPRNT( 23, str, iw, leniw )
                  end if
               end if

               call s5hs
     &            ( Extern, nb, bl, bu, hs, x )
               call s2Amat
     &            ( RowTyp, mnrPrtlvl, m, n, nb,
     &              nnCon, nnJac, nnObj, iObj, numLC, numLIQ,
     &              neA, nlocA, locA, indA, Acol,
     &              bl, bu, RowTypes,
     &              iw, leniw, rw, lenrw )
               lCrash = 4
               call s2crash
     &            ( lCrash, mnrPrtlvl, m, n, nb, nnCon,
     &              iw(iCrash), tCrash,
     &              neA, nlocA, locA, indA, Acol,
     &              kBS, hs, RowTypes, bl, bu, x,
     &              iw, leniw, rw, lenrw )
            end if ! numLIQ > 0
         end if
      end if ! NeedB and numLC > 0

  900 return

 1316 format(1x, 16a4)
 1332 format(1x, 28a4)
 2100 format(' Itn', i7, ': Making linear equality rows feasible')
 2200 format(' Itn', i7, ': Feasible linear equality rows')
 2300 format(' Itn', i7, ': Infeasible linear equality rows')
 2400 format(' Itn', i7, ': Making all linear rows feasible')
 2410 format(' Itn', i7, ': Making the linear rows feasible')

      end ! subroutine s5getB

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5hs
     &   ( Mode, nb, bl, bu, hs, x )

      implicit
     &     none
      integer
     &     Mode, nb, hs(nb)
      double precision
     &     bl(nb), bu(nb), x(nb)

!     ==================================================================
!     s5hs   sets the state vector hs.
!     if mode = 'Internal', s5hs sets hs(j) = -1 or 4 for certain
!        nonbasic variables.  This allows s5price to operate more
!        efficiently.  The internal values of hs are now as follows:
!
!        hs(j) = -1  Nonbasic between bounds (bl     <  x <  bu    )
!        hs(j) =  0  Nonbasic on lower bound (bl-tol <  x <= bl    )
!        hs(j) =  1  Nonbasic on upper bound (bu     <= x <  bu+tol)
!        hs(j) =  2  Superbasic
!        hs(j) =  3  Basic
!        hs(j) =  4  Nonbasic and fixed      (bl     = x  =  bu    )
!
!        where 0 <= tol < the feasibility tolerance.
!
!     if mode = 'External', s5hs changes -1 or 4 values to hs(j) = 0,
!        ready for basis saving and the outside world.
!
!     08 Apr 1992: First version of s5hs.
!     21 Aug 1999: Current version.
!     ==================================================================
      integer
     &     j
!     ------------------------------------------------------------------
      integer            Intern,     Extern
      parameter         (Intern = 0, Extern = 1)
!     ------------------------------------------------------------------
      if (Mode .eq. Intern) then
!        ---------------------------------------------------------------
!        Change nonbasic hs(j) to internal values (including 4 and -1).
!        This may change existing internal values if bl and bu have been
!        changed -- e.g. at the start of each major iteration.
!        ---------------------------------------------------------------
         do j = 1, nb
            if (hs(j) .le. 1) then
               if (bl(j) .eq. bu(j)) then
                  hs(j) =  4
               else if (x(j) .le. bl(j)) then
                  hs(j) =  0
               else if (x(j) .ge. bu(j)) then
                  hs(j) =  1
               else
                  hs(j) = -1
               end if
            end if
         end do

      else if (Mode .eq. Extern) then
!        ---------------------------------------------------------------
!        Change hs to external values.
!        Some nonbasic hs(j) may be 4 or -1.  Change them to 0.
!        ---------------------------------------------------------------
         do j = 1, nb
            if (hs(j) .eq. 4  .or.  hs(j) .eq. -1) hs(j) = 0
         end do
      end if

      end ! subroutine s5hs

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Inf
     &   ( nBS, featol, infBnd,
     &     nInf, sInf, feasType, blBS, buBS, gBS, xBS )

      implicit
     &     none
      integer
     &     nBS, nInf, feasType(nBS)
      double precision
     &     featol, infBnd, sInf,
     &     blBS(nBS), buBS(nBS), gBS(nBS), xBS(nBS)

!     ==================================================================
!     s5Inf  computes the sum and number of infeasibilities.
!
!     On entry:
!     --------
!     featol        is the current feasibility tolerance.
!     infBnd        is the magnitude of an infinite bound.
!     gBB(nBS)      is the zero vector.
!     xBS(nBS)      is the vector of basics and superbasics
!
!     On exit:
!     --------
!     nInf          is the number violated nonelastic bounds.
!     gBS(nBS)      is the rhs for the equations for pi.
!     feasType(nBS) defines the type of infeasibility.
!
!     feasType     x(j)                          Meaning
!     --------     -----                         -------
!      -2      infeasible                           x(j) .le. bl(j)-tol
!       0        feasible            bl(j)-tol .le. x(j) .le. bu(j)+tol
!      +2      infeasible                           x(j) .ge. bu(j)+tol
!
!     feasType is used in s5step.
!
!     29 Oct 1993: First version of s5Inf.
!     26 May 2013: infBnd used to identify infinite bounds.
!     ==================================================================
      integer
     &     infType, k
      double precision
     &     infLow, infUpp, res, xk
!     ------------------------------------------------------------------
      double precision   zero,          one,          point9
      parameter         (zero = 0.0d+0, one = 1.0d+0, point9 = 9.9d-1)
!     ------------------------------------------------------------------
      infUpp =  point9*infBnd
      infLow = -point9*infBnd

      nInf   =  0
      sInf   =  zero

      do k = 1, nBS
         infType = 0
         xk      = xBS(k)

         ! Check if the lower bound (if any) is violated.

         if (blBS(k) .gt. infLow) then
            res = blBS(k) - xk
            if (res .gt. featol) then
               gBS(k)  = -one
               inftype = -2
            end if
         end if

         ! Check if the upper bound (if any) is violated.

         if (infType .eq. 0  .and.  buBS(k) .lt. infUpp) then
            res = xk - buBS(k)
            if (res .gt. featol) then
               gBS(k)  =  one
               infType =  2
            end if
         end if

         if (infType .ne. 0) then
            nInf    = nInf + 1
            sInf    = sInf + res
         end if

         feasType(k) = infType

      end do

      end ! subroutine s5Inf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5LG
     &   ( m, n, nb, nnCon, nrhsLP,
     &     neA, nlocA, locA, indA, Acol,
     &     bl, bu, nrhs0, nrhs, rhs, rhsLP, x, y, rw, lenrw )

      implicit
     &     none
      integer
     &     lenrw, m, n, nb, neA, nlocA, nnCon, nrhs0, nrhs, nrhsLP,
     &     locA(nlocA), indA(neA)
      double precision
     &     Acol(neA), bl(nb), bu(nb), x(nb),
     &     rhs(nrhs0), rhsLP(m), y(m), rw(lenrw)

!     ==================================================================
!     s5LG  relaxes the linear inequality constraints for an LP that
!     gets feasible for the linear equality constraints. A right-hand
!     side is computed that makes the relaxed rows satisfied at the
!     current x.   Then x will not be disturbed more than
!     necessary during a Warm start.
!
!     02 Apr 2005: First version of s5LG.
!     ==================================================================
      integer
     &     i, j
      double precision
     &     eps0, infBnd
!     ------------------------------------------------------------------
      integer            Normal
      parameter         (Normal = 0)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      eps0   = rw(  2) ! eps**(4/5)
      infBnd = rw( 70) ! definition of an infinite bound

      nrhsLP = m
      call s2Aprod
     &   ( Normal, eps0, neA, nlocA, locA, indA, Acol,
     &     one, x, n, zero, y, m )
      call daxpy( nrhsLP, (-one), x(n+1), 1, y, 1 )
      if (nrhs .gt. 0)
     &call daxpy( nrhs  , (-one), rhs   , 1, y, 1 )
      call dload( nrhsLP,  zero , rhsLP , 1 )

      do i = nnCon+1, m
         j = n + i
         if (bl(j) .lt. bu(j)) then
            bl(j)    = - infBnd
            bu(j)    = + infBnd
            rhsLP(i) =   y(i)
         end if
      end do

      end ! subroutine s5LG

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5LPitn
     &   ( iExit,
     &     Feasible, Increase, Needpi, Elastic,
     &     m1, m, nb, nDegen, elastics, LUrequest,
     &     kp, jBq, jSq, jBr, jSr, jq,
     &     featol, pivot, step, tolinc,
     &     eType, eState, feasType, hs, kBS,
     &     bl, bu, blQP, buQP, blBS, buBS,
     &     x, xBS, pBS, y1,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     Feasible, Increase, Needpi, Elastic
      integer
     &     iExit, m1, m, nb, nDegen, elastics, LUrequest, kp,
     &     jBq, jSq, jBr, jSr, jq, leniw, lenrw,
     &     feasType(m1), kBS(m1), eType(nb), eState(nb), hs(nb),
     &     iw(leniw)
      double precision
     &     featol, pivot, step, tolinc, bl(nb), bu(nb), blQP(nb),
     &     buQP(nb), blBS(m1), buBS(m1), x(nb), xBS(m1),
     &     pBS(m1), y1(m1), rw(lenrw)

!     ==================================================================
!     s5LPitn computes one step of the primal simplex method.
!     jq is the variable entering the basis and djq is its reduced cost.
!
!      iExit       Status
!      -----       ------
!        0         Normal exit
!        1         LP is unbounded
!
!     10 Sep 1991: First version based on Minos routine m5lpit.
!     02 Aug 1996: First min sum version added by PEG.
!     12 Jul 1997: Thread-safe version.
!     01 Aug 2003: snPRNT adopted.
!     08 Apr 2008: eState accessed in elastic mode only.
!     04 Jul 2008: Modify both bl and bu in elastic phase 1.
!     06 May 2015: x(jr) set exactly on its bound.
!     =================================================================
      character
     &     str*50
      logical
     &     HitLow, Move, OnBound, Unbounded
      integer
     &     inform, infpiv, jr, jeState, jrstate, jqState, LUmod,
     &     mtry, ntry
      double precision
     &     bound, exact, bigdx, infBnd, pNorm, sclPiv, stepP, stepMax,
     &     tolpiv, tolP0, tolP
      external
     &     dnormi
      double precision
     &     dnormi
!     ------------------------------------------------------------------
      integer            xBStox
      parameter         (xBStox = 1)
      parameter         (mtry   = 6)

      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one  = 1.0d+0)
      double precision   ten
      parameter         (ten    =10.0d+0)

      parameter         (LUmod  = 216) ! number of LU mods
!     ------------------------------------------------------------------
      tolpiv    = rw( 60) ! excludes small elements of pBS.
      infBnd    = rw( 70) ! definition of an infinite bound
      bigdx     = rw( 72) ! unbounded step.

      iExit     = 0

!     s5step assumes that the first m1 components of xBS can Move.

      jSq          =       jq   ! Entering superbasic
      jSr          =       jq   ! leaving  superbasic
      feasType(m1) =       0
      kBS  (m1)    =       jq
      xBS  (m1)    =     x(jq)
      blBS (m1)    =  blQP(jq)
      buBS (m1)    =  buQP(jq)

      !=================================================================
      ! Set eState(jq) and the elastic parts of blBS and buBS.
      !=================================================================
      if (Elastic) then

         ! If the new superbasic is an elastic variable
         ! and it wants to move infeasible, set its elastic state.

         if (eType(jq) .gt. 0  .and.  eState(jq) .eq. 0) then
            jqState = hs(jq)
            if (Increase) then
               if (jqState .eq. 1  .or.  jqState .eq. 4) then
                  elastics   =  elastics + 1
                  eState(jq) =  2
                  blQP(jq)   =  buQP(jq)
                  buQP(jq)   = +infBnd
                  blBS(m1)   =  blQP(jq)
                  buBS(m1)   =  buQP(jq)
               end if
            else
               if (jqState .eq. 0  .or.  jqState .eq. 4) then
                  elastics   =  elastics + 1
                  eState(jq) =  1
                  buQP(jq)   =  blQP(jq)
                  blQP(jq)   = -infBnd
                  blBS(m1)   =  blQP(jq)
                  buBS(m1)   =  buQP(jq)
               end if
            end if
         end if
      end if ! Elastic mode

!     ==================================================================
!     Select a variable to be dropped from B.
!     s5step  uses the (m+1)th element of  blBS, buBS, xBS and pBS.
!     ==================================================================
      if (Increase) then
         call dscal ( m, (-one), pBS, 1 )
         pBS(m1) = +one
      else
         pBS(m1) = -one
      end if

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
      pNorm   = dnormi( m1, pBS, 1 )
      stepMax = bigdx /pNorm
      sclPiv  = one
      tolP0   = tolpiv
      tolP    = tolpiv*pNorm
      ntry    = 0

!+    Repeat
  200    tolP  = tolP /sclPiv
         tolP0 = tolP0/sclPiv
         call s5step
     &      ( m1, nDegen,
     &        featol, infBnd, stepMax, tolinc, tolP,
     &        feasType, blBS, buBS, xBS, pBS,
     &        HitLow, Move, OnBound, Unbounded,
     &        infpiv, kp, bound, exact, step, stepP )

         sclPiv = ten
         ntry   = ntry + 1

!+    until    ( infpiv .eq. 0 .and. (.not.Unbounded .or. Feasible) .or.
!+                 ntry .ge. mtry)
      if (.not.((infpiv .eq. 0 .and. (.not.Unbounded .or. Feasible)).or.
     &             ntry .ge. mtry)) go to 200

      if (.not. Unbounded) then
         !--------------------------------------------------------------
         ! Update the basic variables xBS and copy them into x.
         !--------------------------------------------------------------
         jr      = kBS(kp)

         call daxpy
     &      ( m1, step, pBS, 1, xBS, 1 )
         call s5BSx
     &      ( xBStox, m1, nb, kBS, x, xBS )

         if (OnBound) then
            x(jr) = bound
         else if (HitLow) then
            x(jr) = min( x(jr), blQP(jr) )
         else
            x(jr) = max( x(jr), buQP(jr) )
         end if

         if (kp .eq. m1) then
            !-----------------------------------------------------------
            ! Variable jq reaches its opposite bound.
            !-----------------------------------------------------------
            if (Increase) then
               hs(jq) = 1
            else
               hs(jq) = 0
            end if
            feasType(kp) = 0
            pivot        = zero
            if (.not. Feasible  .or.  Elastic) Needpi = .true.

         else
            !-----------------------------------------------------------
            ! Variable jq replaces the kp-th variable of  B.
            ! It could be a fixed variable, whose new state must be 4.
            !-----------------------------------------------------------
            Needpi =  .true.
            jBq    =  jq
            jBr    =  jr
            hs(jq) =  3
            pivot  = -pBS(kp)

            if (Elastic) then
               jeState = eState(jr)
            else
               jeState = 0
            end if

            if (jeState .eq. 0) then

               ! Normal variable hits its bound.

               if (blBS(kp) .eq. buBS(kp)) then
                  jrstate = 4
               else if (HitLow) then
                  jrstate = 0
               else
                  jrstate = 1
               end if
            else

               ! Elastic x(jr) hits its bound.
               ! Reset the true upper and lower bounds.

               eState(jr) = 0
               elastics   = elastics - 1
               blQP(jr)   = bl(jr)
               buQP(jr)   = bu(jr)

               if (jeState .eq. 1) then

                  if (blQP(jr) .eq. buQP(jr)) then
                     jrstate =  4
                  else if (OnBound) then
                     jrstate =  0
                  else if (x(jr) .lt. buQP(jr)) then
                     jrstate = -1
                  else
                     jrstate =  1
                  end if

               else !   jeState .eq. 2

                  if (blQP(jr)   .eq. buQP(jr)) then
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

            hs(jr)   = jrstate
            kBS(kp)  = jq
            xBS(kp)  = xBS(m1)
            blBS(kp) = blBS(m1)
            buBS(kp) = buBS(m1)

            ! Update the LU factors.

            iw(LUmod)  = iw(LUmod) + 1
            call s2Bmod
     &         ( inform, kp, m, y1, iw, leniw, rw, lenrw )
            if (inform .eq. -1) LUrequest = 5 ! Singular after LU mod
            if (inform .eq.  2) LUrequest = 6 ! Unstable LU mod
            if (inform .eq.  7) LUrequest = 7 ! Insufficient free memory
         end if ! kp ne m1

      else

         ! Apparently the solution is unbounded.

         if (Increase) then
            write(str, 1000) jq
         else
            write(str, 1100) jq
         end if
         call snPRNT( 21, str, iw, leniw )
         iExit = 1              ! Unbounded direction
      end if

      return

 1000 format(' Variable', i6, '  can increase indefinitely')
 1100 format(' Variable', i6, '  can decrease indefinitely')

      end ! subroutine s5LPitn

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5price
     &   ( Elastic, Feasible, Increase, Gotg, subOptimize,
     &     itn, m, n, nb, leng, neg,
     &     frozen, nonOpt, weight, signObj, piNorm,
     &     jq, djq, kPrc, toldj,
     &     neA, nlocA, locA, indA, Acol,
     &     eType, hs, g, pi, rc, x, xFrozen,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     Elastic, Feasible, Gotg, Increase
      integer
     &     m, n, nb, neA, nlocA, frozen, leng, neg,
     &     nonOpt, itn, jq, kPrc, leniw, lenrw, subOptimize,
     &     eType(nb), hs(nb), locA(nlocA), indA(neA), iw(leniw)
      double precision
     &     weight, piNorm, signObj, djq, toldj(3), Acol(neA), g(leng),
     &     pi(m), rc(nb), x(nb), xFrozen(nb), rw(lenrw)

!     ==================================================================
!     s5price  selects a nonbasic variable to enter the basis,
!     using the reduced gradients  dj = g(j) - pi'*a(j).
!
!     This version does partial pricing on both structurals and slacks.
!     Dynamic tolerances are used if partial price is in effect.
!
!     Partial pricing here means sectional pricing, because the three
!     blocks of  (A  -I)  are sliced up into nParPr sections
!     of equal size.  (The last section of each may be a little bigger,
!     since nParPr is unlikely to divide evenly into  n  or  m.)
!
!     input    g      = gradient for nonlinear variables.
!              pi     = pricing vector.
!              kPrc   = the no. of the section where s5price last found
!                       a useful dj.
!                       (kPrc = 0 at the start of each major iteration.)
!              toldj(1:2) hold the current told if partial pricing, for
!                       phase 1 and 2 respectively.  told is used to
!                       determine if a dj is significant.
!              toldj(3) holds the specified optimality tolerance.
!              biggst   keeps track of the biggest dj found during the
!                       last scan of all sections of ( A  -I ).
!
!     output   kPrc   = the last section scanned.
!              nonOpt = the no. of useful djs found in that section.
!              jq     = best column found.
!              djq    = best dj.
!              toldj(1:2) save the current told if partial pricing.
!              Increase says if variable jq should increase or decrease.
!
!     In the code below,
!     the next section of  A  contains nPr1 structurals (j1+1 thru k1),
!     the next section of -I  contains nPr2 slacks      (j2+1 thru k2).
!     If  nParPr  is rather large, either nPr1 or nPr2 could be zero,
!     but not both.
!
!     If subOptimize > 0,  variables that haven't moved are not
!     priced, thereby limiting the number of superbasic variables.
!
!     ------------------------------------------------------------------
!     09 Aug 1992: First version of s5pric based on Minos 5.4 m5pric.
!     29 Jul 1996: Multiple pricing removed.
!     05 Aug 1996: First version with Elastic mode.
!     12 Jul 1997: Thread-safe version.
!     22 Dec 1999: subOptimize implemented.
!     31 Mar 2001: Free variables with hs(j)=-1 eligible for the basis.
!     19 Sep 2014: Name changed from s5pric to s5price.
!     19 Sep 2014: Added case for a temporarily fixed elastic.
!     21 Feb 2015: Deleted code that removes nonbasics floating free
!                  between their bounds. This code can cause a bogus
!                  unbounded indication in s5step.
!     ==================================================================
      character
     &     str*100
      integer
     &     j, eTypej, jj, js, jSlk, j1, j2, k1, k2,
     &     kPsav, lprDbg, lvldj, np, nParPr, nPrc, nParP, nPr1, nPr2
      double precision
     &     d, d1, d2, dj, dj1, dj2, djmax, infBnd, told, tolmin
!     ------------------------------------------------------------------
!-->  parameter         (zero = 0.0d+0, reduce = 0.25d+0)
      double precision   zero,          reduce
      parameter         (zero = 0.0d+0, reduce = 0.2d+0 )
!     ------------------------------------------------------------------
      infBnd    = rw( 70) ! definition of an infinite bound
      lprDbg    = iw( 85) ! > 0    => private debug print
      nParPr    = iw( 94) ! # of partial pricing sections

      djmax     = - infBnd
      djq       =   zero

      jq        = 0
      frozen    = 0
      nonOpt    = 0
      nPrc      = 0
      nParP     = nParPr
      nPr1      = n  / nParP
      nPr2      = m  / nParP
      if (max( nPr1, nPr2 ) .le. 0) nParP = 1

!     Set the tolerance for a significant dj.

      tolmin    = toldj(3) * piNorm
      if (Feasible) then
         lvldj  = 2
      else
         lvldj  = 1
      end if
      told      = toldj(lvldj)
      if (nParP .eq. 1) told = tolmin

!     Set pointers to the next section of  A  and  -I.
!     nPrc counts how many sections have been scanned in this call.
!     kPrc keeps track of which one to start with.

  100 nPrc = nPrc + 1
      kPrc = kPrc + 1
      if (kPrc .gt. nParP) kPrc = 1

      nPr1 = n  / nParP
      j1   =      (kPrc - 1)*nPr1
      k1   = j1 + nPr1
      if (kPrc .eq. nParP) k1 = n
      nPr1 = max( 0, k1-j1 )

      nPr2 = m  / nParP
      j2   = n  + (kPrc - 1)*nPr2
      k2   = j2 + nPr2
      if (kPrc .eq. nParP) k2 = nb
      nPr2 = max( 0, k2-j2 )

!     ------------------------------------------------------------------
!     Main loops for partial pricing (or full pricing).
!     Compute reduced costs rc(*)
!     for the kPrc-th section of structurals
!     and the kPrc-th section of slacks.
!     ------------------------------------------------------------------
      call s5rc
     &   ( j1+1, k1, Gotg, m, n, leng, neg, signObj,
     &     neA, nlocA, locA, indA, Acol,
     &     hs, g, pi, rc )

      do j = j2+1, k2
         rc(j)  = pi(j-n)
      end do

!     ------------------------------------------------------------------
!     Main loop for pricing structural and slack reduced costs.
!     dj is rc(j), the reduced cost.
!     d  is -dj or +dj, depending on which way x(j) can move.
!     We are looking for the largest d (which will be positive).
!     ------------------------------------------------------------------
      np   = nPr1 + nPr2
      j    = j1
      jSlk = nPr1 + 1

      do jj = 1, np
         if (jj .eq. jSlk) j = j2
         j       = j + 1
         js      = hs(j)

         if (js .le. 1) then
            dj     = rc(j)

            if      (js .eq. 0) then
!              xj  is allowed to increase.
               d      = - dj
            else if (js .eq. 1) then
!              xj  is allowed to decrease.
               d      =   dj
            else
!              js is -1.
!              xj  is free to move either way.
               d      = abs( dj )
            end if

            if (subOptimize .gt. 0) then
               if (x(j) .eq. xFrozen(j)) then
                  if (d  .gt. told) then
                     frozen = frozen + 1
                     d      = zero
                  end if
               end if
            end if

!           See if this dj is significant.
!           Also see if it is the biggest dj so far.

            if (d  .gt. told) nonOpt = nonOpt + 1
            if (djmax .lt. d) then
               djmax  = d
               djq    = dj
               jq     = j
               kPsav  = kPrc
            end if
         end if
      end do

      if (Elastic) then
!        ---------------------------------------------------------------
!        Scan this section again, looking for nonbasic elastics.
!        ---------------------------------------------------------------
!        Compute reduced costs rc(j) for fixed nonbasic columns.
!        (These columns are skipped in s5rc)

         call s5Erc
     &      ( j1+1, k1, Gotg, m, n, leng, neg, signObj,
     &        neA, nlocA, locA, indA, Acol,
     &        eType, hs, g, pi, rc )

         j    = j1
         do jj = 1, np
            if (jj .eq. jSlk) j = j2
            j      = j + 1
            eTypej = eType(j)

            if (eTypej .gt. 0) then
               js  = hs(j)
               dj  = rc(j)

               if      (js .eq. 0) then
!                 ------------------------------------------------------
!                 Nonbasic at its lower bound.
!                 An elastic xj can decrease through the bound.
!                 ------------------------------------------------------
                  if (eTypej .eq. 1  .or.  eTypej .eq. 3) then
                     dj  =   dj - weight
                     d   =   dj
                  end if

               else if (js .eq. 1) then
!                 ------------------------------------------------------
!                 Nonbasic at its upper bound.
!                 The default is to allow xj to decrease.
!                 However, an elastic xj can increase through the bound.
!                 ------------------------------------------------------
                  if (eTypej .eq. 2  .or.  eTypej .eq. 3) then
                     dj  =   dj + weight
                     d   = - dj
                  end if

               else if (js .eq. 4  .or. js .eq. -1) then
!                 ------------------------------------------------------
!                 Fixed elastic variable.
!                 xj is free to move either way.
!                 ------------------------------------------------------
                  if (eTypej .eq. 2) then
                     dj1 =   zero
                     d1  =   zero
                  else
                     dj1 =   dj - weight
                     d1  =   dj1
                  end if

                  if (eTypej .eq. 1) then
                     dj2 =   zero
                     d2  =   zero
                  else
                     dj2 =   dj + weight
                     d2  = - dj2
                  end if

                  if (d1 .ge. d2) then
!                    xj  is allowed to decrease.
                     dj =   dj1
                     d  =   d1
                  else
!                    xj  is allowed to increase.
                     dj =   dj2
                     d  =   d2
                  end if
               else
                  d  = zero
                  dj = zero
               end if

               if (subOptimize .gt. 0) then
                  if (x(j) .eq. xFrozen(j)) then
                     if (d  .gt. told) then
                        frozen = frozen + 1
                        d      = zero
                     end if
                  end if
               end if

!              See if this dj is significant.
!              Also see if it is the biggest dj so far.

               if (d  .gt. told) nonOpt = nonOpt + 1
               if (djmax .lt. d) then
                  djmax  = d
                  djq    = dj
                  jq     = j
                  kPsav  = kPrc
               end if
            end if
         end do
      end if

!     ------------------------------------------------------------------
!     End of loop looking for biggest dj in the kPrc-th section.
!     ------------------------------------------------------------------
      if (nonOpt .eq. 0) then
         if (nParP .gt. 1) then
!           ============================================================
!           No significant dj has been found.  (All are less than told.)
!           Price the next section, if any remain.
!           ============================================================
            if (nPrc .lt. nParP) go to 100

!           ============================================================
!           All sections have been scanned.  Reduce told
!           and grab the best dj if it is bigger than tolmin.
!           ============================================================
            if (djmax .gt. tolmin) then
               nonOpt = 1
               kPrc   = kPsav
               told   = max( reduce * djmax, tolmin  )
               toldj(lvldj) = told
               if (lprDbg .ge. 1) then
                 write(str, 1000) itn, told, piNorm, weight
                 call snPRNT( 23, str, iw, leniw )
               end if
            end if
         end if
      end if

      Increase = djq .lt. zero

 1000 format(' Itn', i7, ': toldj =', 1p, e8.1,
     &       '    Norm pi =', e8.1, '    weight = ', e8.1)

      end ! subroutine s5price

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5rc
     &   ( j1, j2, Gotg, m, n, leng, neg, signObj,
     &     neA, nlocA, locA, indA, Acol,
     &     hs, g, pi, rc )

      implicit
     &     none
      logical
     &     Gotg
      integer
     &     j1, j2, m, n, nlocA, leng, neg, neA,
     &     locA(nlocA), indA(neA), hs(n)
      double precision
     &     signObj, Acol(neA), g(leng), pi(m), rc(n)

!     ==================================================================
!     s5rc   computes reduced costs rc(j) for nonbasic columns of A
!     in the range j = j1 to j2.  It is called by s5price.
!
!     The loop for computing dj for each column could conceivably be
!     optimized on some machines.  However, there are seldom more than
!     5 or 10 entries in a column.
!
!     Note that we could skip fixed variables by passing in the bounds
!     and testing if bl(j) .eq. bu(j), but these are relatively rare.
!     But see comment for 08 Apr 1992 in m5pric.
!
!     31 Jan 1992: First version of s5rc.
!     08 Apr 1992: Internal values of hs(j) are now used, so fixed
!                  variables (hs(j) = 4) are skipped as we would like.
!     03 Apr 1999: Linear objective stored as row 0 of A.
!     11 Apr 1999: Current version.
!     ==================================================================
      integer
     &     i, j, l
      double precision
     &     dj
!     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
!     ------------------------------------------------------------------
      do j = j1, j2
         if (hs(j) .le. 1) then
            dj    = zero
            do l  = locA(j), locA(j+1)-1
               i  = indA(l)
               dj = dj  +  pi(i) * Acol(l)
            end do
            rc(j) = - dj
         end if
      end do

!     Include the nonlinear objective gradient if relevant.

      if (Gotg) then
         if (j1 .le. neg) then
            do j = j1, min( j2, neg )
               if (hs(j) .le. 1) rc(j) = rc(j) + signObj*g(j)
            end do
         end if
      end if

      end ! subroutine s5rc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5setpi
     &   ( iExit,
     &     m, Checkpi, condZmax, rhsNorm, piNorm, rhs, pi,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     Checkpi
      integer
     &     iExit, leniw, lenrw, m, iw(leniw)
      double precision
     &     condZmax, rhsNorm, piNorm, pi(m), rhs(m), rw(lenrw)

!     ==================================================================
!     s5setpi solves  B' pi = rhs.  Beware -- rhs is altered by s2Bsol.
!     If a new x has just been computed, the norm is computed by dnormj.
!
!     On Exit:
!      iExit  = -25 estimate of condZ > condZmax.
!      iExit  = -11 if pi contains a NaN or Inf.
!      iExit  =   0 if pi was computed successfully.
!      iExit  >   0 if there was an unexpected error in the solver.
!
!     08 Aug 1996: First version of s5setp.
!     16 Nov 2001: Infinity norm used for piNorm (no longer dnrm1s).
!     19 Sep 2014: Name changed to s5setpi.
!     10 May 2015: cond(Z) estimator added.
!     12 May 2015: New argument rhsNorm and new call to s5setCondZmax..
!     ==================================================================
      external
     &     dnormi, dnormj
      double precision
     &     condZ, dnormi, dnormj, flMax
!     ------------------------------------------------------------------
      integer            WithBt
      parameter         (WithBt = 2)

      double precision   one
      parameter         (one = 1.0d+0)
!     ------------------------------------------------------------------
      flMax = rw(  8) ! est. of the largest pos. real

      iExit = 0

      call s2Bsol
     &   ( iExit, WithBt, m, rhs, pi,  iw, leniw, rw, lenrw )
      if (iExit .gt. 0) return

      if (Checkpi) then
         piNorm  = dnormj( m, pi, 1 )
         if (piNorm .lt. flMax) then ! false if pi = inf, nan
            condZ   = piNorm/rhsNorm ! Crude est. of cond(Z)
            if (condZ .gt. condZmax) then
               iExit    = -25   ! defines LUrequest = 25
               call s5setCondZmax( condZmax, rw, lenrw )
            end if
            piNorm = max( piNorm, one ) ! Defined regardless of iExit
         else
            iExit  = -11        ! defines LUrequest = 11
         end if
      else
         piNorm = max(dnormi( m, pi, 1 ), one)
      end if

      end ! subroutine s5setpi

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5setx
     &   ( iExit, Task, itn,
     &     m, n, nb, nBS, rowError,
     &     neA, nlocA, locA, indA, Acol,
     &     kBS, xBS, nrhs0, nrhs, rhs, x, y, y1,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     Task, iExit, itn, nrhs0, m, n, nb, nBS, neA, nlocA, nrhs,
     &     leniw, lenrw, locA(nlocA), indA(neA), kBS(nBS), iw(leniw)
      double precision
     &     rowError, Acol(neA), rhs(nrhs0), xBS(nBS), x(nb), y(m),
     &     y1(m), rw(lenrw)

!     ==================================================================
!     s5setx performs the following functions:
!      Task            Function
!      ====           ========
!       0 (Resetx)    the basic components of x are computed to satisfy
!                     Ax - s = b; that is  (A -I)*x = b. Then a row
!                     check is performed to see how well  (A -I)*x = b
!                     is satisfied.  y is set to be the row residuals,
!                     y = b - Ax + s,  and the row error is norm(y).
!
!       1 (GetRes)    just get the row error.
!
!     The row error is a measure of how well x satisfies (A -I)*x = b.
!
!     18 Nov 1991: First version of s5setx based on Minos routine m5setx.
!     12 Jul 1997: Thread-safe version.
!     25 Jul 2003: Realized dx can sometimes be a NAN (or INF) but
!                  norm(NAN) > 1/eps is always false.
!     ==================================================================
      character
     &     str*110
      external
     &     jdamax, dnormi, dnormj
      logical
     &     BigRes, Goodx
      integer
     &     imax, lprDbg, jdamax
      double precision
     &     dnormi, dnormj, eps0, rmax, tolrow, xNorm, dxNorm
!     ------------------------------------------------------------------
      integer            Normal,     WithB
      parameter         (Normal = 0, WithB  = 1)
      integer            Resetx
      parameter         (Resetx = 0)
      integer            xtoxBS,     xBStox
      parameter         (xtoxBS = 0, xBStox = 1)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      eps0      = rw(  2) ! eps**(4/5)       IEEE DP  3.00e-13
      tolrow    = rw( 61) ! tolerance for the row error.
      lprDbg    = iw( 85) ! > 0    => private debug print

      iExit    = 0

      call s5BSx
     &   ( xtoxBS, nBS, nb, kBS, x, xBS )
      xNorm     = dnormi( nBS, xBS, 1 )
      dxNorm    = zero
      Goodx     = .true.

!     ------------------------------------------------------------------
!     Compute row residuals  y  =  rhs - (A -I)*x.
!     The slack columns are done separately.
!     ------------------------------------------------------------------
      if (nrhs .gt. 0)
     &   call dcopy ( nrhs, rhs, 1, y        , 1 )
      if (nrhs .lt. m)
     &   call dload ( m-nrhs, zero, y(nrhs+1), 1 )

      call s2Aprod
     &   ( Normal, eps0,
     &     neA, nlocA, locA, indA, Acol,
     &     (-one), x, n, one, y, m )
      call daxpy
     &   ( m, one, x(n+1), 1, y, 1 )

!     ------------------------------------------------------------------
!     Do a row check, perhaps after recomputing the basic x.
!     -----------------------------------------------------------------
      if (Task .eq. Resetx) then
!        ================================================
!        Extract xBS, the basics and superbasics, from x.
!        See if iterative refinement is worth doing.
!        ================================================
         rowError = dnormj( m, y, 1 )
         BigRes   = rowError .gt. eps0

         if (BigRes) then
!           ------------------------------------------------------------
!           Compute a correction to basic x from  B*y1 = y.
!           Extract the basic and superbasic variables from x.
!           Set basic x = x + y1.
!           Store the new basic variables in x.
!           ------------------------------------------------------------
            call s2Bsol
     &         ( iExit, WithB, m, y, y1, iw, leniw, rw, lenrw )
            if (iExit .ne. 0) return
            dxNorm = dnormj( m, y1, 1 )
            Goodx  = dxNorm*eps0 .le. one ! false if dxNorm = inf, nan

            if (Goodx) then
               call daxpy
     &            ( m, one, y1, 1, xBS, 1 )
               call s5BSx
     &            ( xBStox, m, nb, kBS, x, xBS )

!              Compute  y  =  rhs  -  (A -I)*x  again for the new x.

               if (nrhs .gt. 0)
     &            call dcopy ( nrhs, rhs, 1, y, 1 )
               if (nrhs .lt. m)
     &            call dload ( m-nrhs, zero, y(nrhs+1), 1 )

               call s2Aprod
     &            ( Normal, eps0,
     &              neA, nlocA, locA, indA, Acol,
     &              (-one), x, n, one, y, m )
               call daxpy
     &            ( m, one, x(n+1), 1, y, 1 )
            else
               iExit = 11      ! big dx (may be nan or inf)
            end if
         end if ! BigRes
      end if ! Task .eq. Reset

!     Find the norm of xBS, the basic and superbasic x.
!     Find the maximum row residual.

      imax   = jdamax( m, y, 1 )
      if (imax .gt. 0) then
         rmax     = abs( y(imax) )
         rowError = rmax / (one + xNorm )
      else
         imax     = -imax
         rmax     =  dnormj( m, y, 1 ) ! = flmax!
         rowError =  rmax
      end if

      BigRes = rowError .gt. tolRow
      if (BigRes) iExit = 10

      if (iExit .gt. 0  .or.  lprDbg .ge. 2) then
         write(str, 1000) itn, rmax, imax, xNorm, dxNorm
         call snPRNT( 21, str, iw, leniw )
         write(str, 1001) itn, rmax, imax
         call snPRNT( 22, str, iw, leniw )
      end if

      return

 1000 format(  ' Itn', i7, ': Row check',
     &         '.  Max residual =', 1p, e8.1, ' on row', i8,
     &         '.  Norm x =', e8.1, '.  Norm dx =', e8.1 )
 1001 format(  ' Itn', i7, ': Row check',
     &         '.  Max residual =', 1p, e8.1, ' on row', i8 )

      end ! subroutine s5setx

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5step
     &   ( nBS, nDegen,
     &     featol, infBnd, stepMax, tolinc, tolpiv,
     &     feasType, blBS, buBS, xBS, pBS,
     &     HitLow, Move, OnBound, Unbounded,
     &     infPiv, kp, bound, exact, stepB, stepP )

      implicit
     &     none
      logical
     &     HitLow, Move, OnBound, Unbounded
      integer
     &     infPiv, kp, nBS, nDegen, feasType(nBS)
      double precision
     &     featol, infBnd, stepMax, tolinc, tolpiv, bound, exact,
     &     stepB, stepP, blBS(nBS), buBS(nBS), xBS(nBS), pBS(nBS)

!     ==================================================================
!     s5step  finds a steplength  stepB  such that  xBS + stepB*pBS
!     reaches one of the bounds on  xBS.
!
!     In this  version of s5step,  when x  is infeasible, the  number of
!     infeasibilities  will never  increase.   If the  number stays  the
!     same, the  sum of  infeasibilities will  decrease.  If  the number
!     decreases by one or more,  the sum of infeasibilities will usually
!     decrease also, but  occasionally it will increase  after the step-
!     length stepB is taken.  (Convergence  is still assured because the
!     number has decreased.)
!
!     Two possible steps are computed as follows:
!
!     stepF = the maximum step that can be taken without violating
!             one of the bounds that are currently satisfied.
!
!     stepI = the maximum (nonzero) step that has the property of
!             reaching a  bound that  is currently violated,  subject to
!             the  pivot being  reasonably  close to  the maximum  pivot
!             among infeasible variables.  (stepI is not defined if x is
!             feasible.)
!
!     stepI  is needed  occasionally when  infeasible, to  prevent going
!     unnecessarily far when stepF is  quite large.  It will always come
!     into  effect when  x  is about  to become  feasible.   The sum  of
!     infeasibilities will  decrease initially  as stepB  increases from
!     zero, but may start increasing for larger steps.  choosing a large
!     stepI allows several elements of x  to become feasible at the same
!     time.
!
!     In the end, we take  stepB = stepF  if x is feasible, or if
!     stepI > stepP (where stepP is the perturbed step from pass 1).
!     Otherwise,  we take  stepB = stepI.
!
!     Input parameters
!     ----------------
!     nBS      is  m + 1  for s5lpit,  m + nBS  for s5QPit.
!     stepMax  defines what should be treated as an unbounded step.
!     infBnd   provides insurance for detecting unboundedness.
!              if stepB reaches a bound as large as infBnd, it is
!              classed as an unbounded step.
!     tolpiv   is a tolerance to exclude negligible elements of pBS.
!     featol   is the current feasibility tolerance used by s5QP.
!              Typically in the range 0.5*tolx to 0.99*tolx,
!              where tolx is the featol specified by the user.
!     tolinc   is used to determine stepMin (see below), the minimum
!              positive step.
!     feasType is set by  s5Inf  as follows:
!              feasType(j) = -2  if x(j) .lt. bl(j) - featol
!                          =  0  if x(j)  is feasible
!                          = +2  if x(j) .gt. bu(j) + featol
!     blBS     the lower bounds on the basic and superbasic variables.
!     buBS     the upper bounds on ditto.
!     xBS      the values of       ditto.
!     pBS      the search direction for the basics and superbasics.
!
!
!     Output parameters
!     -----------------
!     HitLow = true  if a lower bound restricted stepB.
!            = false otherwise.
!     Move   = true  if exact ge stepMin (defined at end of code).
!     OnBound = true  if stepB =  exact.
!                    this means that the step stepB moves xBS(kp)
!                    exactly onto one of its bounds, namely bound.
!            = false if the exact step would be too small
!                    ( exact lt stepMin ).
!              (with these definitions,  Move = OnBound).
!     Unbounded = true  if stepB = stepMax.  kp may possibly be zero.
!              the parameters HitLow, Move, OnBound, bound and exact
!              should not be used.
!     infPiv = the number of  indices such  that |pBS(i)| <  tolpiv and
!              xBS(i) + stepP*pBS(i) is infeasible by more than featol.
!     kp     = the index (if any) such that xBS(kp) reaches a bound.
!     bound  = the bound value blBS(kp) or buBS(kp) corresponding
!              to HitLow.
!     exact  = the step that would take xBS(kp) exactly onto bound.
!     stepB  = an allowable, positive steplength.
!              if Unbounded is true,  stepB = stepMax.
!              otherwise,          stepB = max(stepMin, exact).
!     stepP  = the perturbed steplength from pass 1.
!
!     07 Nov 1991: First version based on Minos routine m5chzr.
!     27 Dec 2003: infPiv added to monitor unwanted infeasibilities.
!     12 Dec 2012: Recoded to remove goto statements.
!     13 Dec 2012: Fixed definition of infPiv.
!     25 May 2013: Fixed missing else block giving undefined values of
!                  kp and HitLow.
!     26 May 2013: infBnd used to identify infinite bounds.
!     ==================================================================
      logical
     &     BlockF, BlockI
      integer
     &     j, jtype, jhitI, jhitF
      double precision
     &     delta, infLow, infUpp, pivot, pivabs, pivmaxI, pivmaxF,
     &     res, stepI, stepMin
!     ------------------------------------------------------------------
      double precision   zero,          gamma,          point9
      parameter         (zero = 0.0d+0, gamma = 1.0d-3, point9 = 9.9d-1)
!     ------------------------------------------------------------------
      infUpp =  point9*infBnd
      infLow = -point9*infBnd

      !-----------------------------------------------------------------
      ! First pass.
      ! For feasible variables, find the steplength stepP that reaches
      ! the nearest perturbed (expanded) bound.  stepP will be slightly
      ! larger than the step to the nearest true bound.
      ! For infeasible variables, find the maximum pivot pivmaxI.
      !-----------------------------------------------------------------
      delta   = featol
      stepP   = stepMax
      pivmaxI = zero
      jhitF   = 0

      do j   = 1, nBS
         pivot  = pBS(j)
         pivabs = abs( pivot )

         if (pivabs .gt. tolpiv) then
            jtype  = feasType(j)

            if (pivot .lt. zero) then  ! x is decreasing.
               if (jtype .ge. 0) then

                  ! The lower bound (if any) is satisfied.
                  ! Test for a smaller stepP.

                  if (blBS(j) .gt. infLow) then
                     res = xBS(j) - blBS(j) + delta

                     if (stepP*pivabs .gt. res) then
                        stepP = res / pivabs
                        jhitF = j
                     end if
                  end if

                  ! If the upper bound is violated,
                  ! test if this variable has a bigger pivot.

                  if (jtype .gt. 0) then
                     pivmaxI = max( pivmaxI, pivabs )
                  end if
               end if
            else                       ! x is increasing.
               if (jtype .le. 0) then

                  ! The upper bound (if any) is satisfied.
                  ! Test for a smaller stepP.

                  if (buBS(j) .lt. infUpp) then
                     res = buBS(j) - xBS(j) + delta

                     if (stepP*pivabs .gt. res) then
                        stepP = res / pivabs
                        jhitF = j
                     end if
                  end if

                  ! If the lower bound is violated,
                  ! test if this variable has a bigger pivot.

                  if (jtype .lt. 0) then
                     pivmaxI = max( pivmaxI, pivabs )
                  end if
               end if
            end if
         end if
      end do

      !-----------------------------------------------------------------
      ! Second pass.
      ! For feasible variables, compute the steps without perturbation.
      ! Choose the largest pivot element subject to the step being no
      ! greater than stepP.
      ! For infeasible variables, find the largest step subject to
      ! the pivot element being no smaller than gamma*pivmaxI.
      !-----------------------------------------------------------------
      stepI   = zero
      pivmaxF = zero
      pivmaxI = gamma*pivmaxI
      jhitI   = 0
      infPiv  = 0

      do j = 1, nBS
         pivot  = pBS(j)
         pivabs = abs( pivot )
         jtype  = feasType(j)

         if (pivabs .gt. tolpiv) then

            if (pivot .lt. zero) then  ! x is decreasing.
               if (jtype .ge. 0) then

                  ! The lower bound (if any) is satisfied.
                  ! Test for a bigger pivot.

                  if (pivabs .gt. pivmaxF .and. blBS(j).gt. infLow) then
                     res = xBS(j) - blBS(j)

                     if (pivabs*stepP .ge. res) then
                        pivmaxF = pivabs
                        jhitF   = j
                     end if
                  end if

                  if (jtype .gt. 0) then

                     ! An upper bound is present and violated.
                     ! Test for a bigger stepI.

                     if (pivabs .ge. pivmaxI) then
                        res = xBS(j) - buBS(j)

                        if (pivabs*stepI .lt. res) then
                           stepI = res / pivabs
                           jhitI = j
                        end if
                     end if
                  end if
               end if
            else                       ! x is increasing.
               if (jtype .le. 0) then

                  ! The upper bound (if any) is satisfied.
                  ! Test for a bigger pivot

                  if (pivabs .gt. pivmaxF .and. buBS(j).lt. infUpp) then
                     res = buBS(j) - xBS(j)

                     if (pivabs*stepP .ge. res) then
                        pivmaxF = pivabs
                        jhitF   = j
                     end if
                  end if

                  if (jtype .lt. 0) then

                     ! A lower bound is present and violated.
                     ! Test for a bigger stepI.

                     if (pivabs .ge. pivmaxI) then
                        res = blBS(j) - xBS(j)

                        if (pivabs*stepI .lt. res) then
                           stepI = res / pivabs
                           jhitI = j
                        end if
                     end if
                  end if
               end if
            end if

         else if (jtype .eq. 0  .and.  pivabs .gt. zero) then

            ! Feasible variable with a negligible pivot.
            ! Check the step makes us infeasible by no more than stepP.

            if (pivot  .lt. zero) then ! x is decreasing.
               if (blBS(j) .gt. infLow) then
                  res = xBS(j) - blBS(j) + delta
                  if (stepP*pivabs .gt. res) then
                     infPiv = infPiv + 1
                  end if
               end if
            else                       ! x is increasing.
               if (buBS(j) .lt. infUpp) then
                  res = buBS(j) - xBS(j) + delta
                  if (stepP*pivabs .gt. res) then
                     infPiv = infPiv + 1
                  end if
               end if
            end if
         end if
      end do

      !-----------------------------------------------------------------
      ! See if a feasible and/or infeasible variable blocks.
      !-----------------------------------------------------------------
      BlockF    = jhitF .gt. 0
      BlockI    = jhitI .gt. 0
      Unbounded = .not. (BlockF  .or.  BlockI)

      if (Unbounded) then
         stepB   = stepMax
         Move    = .true.
         OnBound = .false.
      else
         if (BlockF) then
            !-----------------------------------------------------------
            ! Variable hits a bound for which it is currently feasible.
            ! The step length stepF is not used, so no need to get it,
            ! but we know that stepF .le. stepP, the step from pass 1.
            !-----------------------------------------------------------
            kp     = jhitF
            pivot  = pBS(kp)
            HitLow = pivot .lt. zero
         else ! BlockI
            kp     = jhitI
            pivot  = pBS(kp)
            HitLow = pivot .gt. zero
         end if

         ! If there is a choice between stepF and stepI, it is probably
         ! best to take stepI (so that the infeasible variable jhitI can
         ! be kicked out of the basis).
         ! However, we can't if stepI is bigger than stepP.

         if (BlockI  .and.  stepI .le. stepP) then
            kp     = jhitI
            pivot  = pBS(kp)
            HitLow = pivot .gt. zero
         end if

         !--------------------------------------------------------------
         ! Try to step exactly onto a bound, but make sure the exact
         ! stepB is sufficiently positive (exact will be either stepF or
         ! stepI). As featol increases by tolinc each iteration, we know
         ! that a step as large as stepMin (below) will not cause  any
         ! feasible variables to become infeasible (where feasibility
         ! is measured by the current featol).
         !---------------------------------------------------------------
         if (HitLow) then
            bound = blBS(kp)
         else
            bound = buBS(kp)
         end if

         stepMin =            tolinc/abs(pivot)
         exact   = (bound - xBS(kp))/    pivot
         stepB   = max(stepMin, exact)
         OnBound = stepB .eq. exact
         Move    = exact .ge. stepMin

         if (.not. Move) then
            nDegen = nDegen + 1
         end if
      end if

      end ! subroutine s5step
