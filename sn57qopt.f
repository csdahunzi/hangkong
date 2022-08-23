!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn57qopt.f
!
!     s5solve
!     s5defaults  s5Map  s5CallStatus
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5solve
     &   ( iExit,
     &     Solver, startType,
     &     qpLog, Hprod, Hprod1, GotR,
     &     m, n, nb, nnH0, nnH,
     &     nNames, ngQP, ngObj0, ngObj,
     &     iObj, objAdd, objQP, objTrue,
     &     nInf, sInf, nInfE, sInfE,
     &     neA, nlocA, locA, indA, Acol,
     &     neH, nlocH, locH, indH, Hcol,
     &     bl, bu, gObj, Names,
     &     nrhs0, nrhs, rhs, lenx0, nx0, x0,
     &     eType, hs, x, pi, rc, nS,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1, qpLog
      logical
     &     GotR
      integer
     &     iExit, iObj, lencu, leniu, lenru, lencw, leniw,
     &     lenrw, m, n, nb, neA, neH, ngObj0, ngObj, ngQP, nInf, nInfE,
     &     nNames, nnH0, nnH, nlocA, nlocH, nrhs0, nrhs, lenx0, nx0,
     &     nS, startType, locA(nlocA), locH(nlocH), indA(neA),
     &     indH(neH), eType(nb), hs(nb), iu(leniu), iw(leniw)
      double precision
     &     objAdd, objQP, objTrue, sInf, sInfE, Acol(neA), Hcol(neH),
     &     bl(nb), bu(nb), gObj(ngObj0), rhs(nrhs0),
     &     x0(lenx0), x(nb), pi(m), rc(nb), ru(lenru), rw(lenrw)
      character
     &     Solver*6, Names(nNames)*8, cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s5solve solves the current problem.
!
!     On entry
!     ---------
!     the SPECS file has been read,
!     all data items have been loaded (including Acol, indA, locA, ...),
!     and workspace has been allocated within cw, iw and rw.
!     startType = lvlStart from s3argQ.
!
!     On exit,
!     --------
!     iExit  =  0 if an optimal solution was found,
!            =  1 if the problem was infeasible,
!            =  2 if the problem was unbounded,
!            =  3 if the Iteration limit was exceeded,
!           ge  4 if iterations were terminated by some other
!                 error condition (see the SQOPT user's guide).
!
!     01 Oct 1994: First version of s5solve.
!     06 Aug 1996: Min Sum option added.
!     14 Jul 1997: Thread-safe version.
!     02 Aug 2003: snEXIT and snPRNT adopted.
!     08 Mar 2004: Hot starts implemented.
!     16 May 2006: Explicit target itQP added.
!     18 Jun 2008: Added space for iy2, pBS and rg2.
!     07 Mar 2013: mnrPrint changed to mjrPrint in call to s2Amat.
!     20 Sep 2014: Upper and lower bounds copied for elastic mode.
!     03 Nov 2014: neH, indH, locH, Hcol added as arguments.
!     27 Dec 2014: Implemented switch to QN when H is indefinite.
!     ==================================================================
      character
     &     mProb*8, istate(3)*4,  probTag*20, str*132, str2*132
      logical
     &     BadZ, BadZg, BigItn, BigSB, Elastic, FPonly,
     &     IndefQP, Infeasible, NeedB, NeedLU, Needx, Optimal, SubItns,
     &     SwitchToQN, Terminated, Unbounded, UseQP, WeakMin
      integer
     &     cgItn, cgItns, eigH, elastics, gotFac, gotHes,
     &     HvCalls, inform, inewB,
     &     itn, itnlim, itQP, itQPtarget, j, k,
     &     lscales, lblQP, lbuQP,
     &     lblBS, lbuBS, lblSave, lbuSave, eMode, lenR, lgBS, lgQP,
     &     lHdx, leState, lfeasType, linesL, linesS, liy, liy1, liy2,
     &     lkBS, lkx, lpBS, lR, lrg, lrg2, lsSave, lvlObjE, lvlScale,
     &     lxBS, ly, ly1, ly2, maxR, maxS, mBS, minimize, minmax,
     &     mnrHdP, mnrHdS, mnrPrint, mjrPrint, nDegen, nFac, ngQP0,
     &     nkx, nnb, nnCon0, nnCon, nnObj, nnJac, numLC, numLIQ,
     &     preCon, printLevel, probType, QPmode, QPsolver, qpStat,
     &     sqStat, subOptimize, typeLU
      double precision
     &     condZHZ, condZmax0, condZmax, degen, dnormi, Hcondbnd,
     &     objInf, objLP, objMin, piNorm, pNorm1, pNorm2, rgNorm,
     &     scaleObj, signObj, targetH, tolOptFP, tolOptQP, tolx,
     &     vimax, wtInf0, wtInf, xNorm, Fx(1)
      external
     &     dnormi
!     ------------------------------------------------------------------
      integer            BS         , BT
      parameter         (BS      = 2, BT     = 3)
      integer            QPChol,      CG,         QN
      parameter         (QPChol  = 0, CG     = 1, QN   = 2)
      integer            FP,          LP,         QP
      parameter         (FP      = 0, LP     = 1, QP   = 2)
      integer            SEMDEF
      parameter         (SEMDEF  = 0)
      integer            Fix
      parameter         (Fix     = 0)
      integer            NO,          YES
      parameter         (NO      = 0, YES    = 1)
      integer            SaveB,       PrintS,     Wrap
      parameter         (SaveB   = 0, PrintS = 1, Wrap = 1)
      integer            Stats
      parameter        ( Stats   = 1 )

      parameter         (lvlScale =  75) ! scale option
      parameter         (HvCalls  = 188) ! Hessian- vector products
      parameter         (QPmode   = 208) ! Current QP solver
      parameter         (preCon   = 209) ! Current precon mode
      parameter         (nFac     = 210) ! # of LU factorizations
      parameter         (linesL   = 220) ! # lines in log     file
      parameter         (linesS   = 221) ! # lines in summary file
      parameter         (mnrHdP   = 223) ! >0 => Mnr heading for iPrint
      parameter         (mnrHdS   = 225) ! >0 => Minor heading for iSumm
      parameter         (gotFac   = 230) ! Save the LU factors
      parameter         (gotHes   = 231) ! Save the reduced Hessian
      parameter         (qpStat   = 235) ! QP user-routine call-status
      parameter         (cgItns   = 386) ! Number of symmlq iterations
      parameter         (cgItn    = 387) ! symmlq itns for last minor

      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      tolOptFP  = rw( 51) ! Minor Phase 1 Opt tol
      tolOptQP  = rw( 52) ! Minor Phase 2 Opt tol
      tolx      = rw( 56) ! Minor feasibility tolerance.
      Hcondbnd  = rw( 85) ! bound on the condition of Hz
      condZmax0 = rw( 86) ! Initial bound on the condition est of Z
      wtInf0    = rw( 88) ! infeasibility weight

      nnObj     = iw( 22) ! # of objective variables
      lenR      = iw( 28) ! R(lenR) is the reduced Hessian factor
      maxR      = iw( 52) ! max columns of R.
      maxS      = iw( 53) ! max # of superbasics
      QPsolver  = iw( 55) ! = 0:1:2   => QPChol:CG:QN QP solver
      eMode     = iw( 56) ! >0    => use elastic mode
      lvlObjE   = iw( 73) ! Elastic objective type
      minmax    = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
      itnlim    = iw( 89) ! limit on total iterations
      mnrPrint  = iw( 93) ! Minor print level
      mjrPrint  = iw( 92) ! Major print level
      iNewB     = iw(124) ! new basis file
      minimize  = iw(199) ! 1 (-1)    => minimize (maximize)
      nkx       = iw(247) ! dimension of kx and its inverse, kxN

      mProb     = cw( 51) ! Problem name

      ! Addresses

      lkx       = iw(251) ! j  = kx (jN) => col j of Jcol is variable jN
      lfeasType = iw(284) ! feasType(mBS) = feasibility types
      leState   = iw(285) ! eState(nb)    = status of elastics
      lkBS      = iw(292) ! kBS(mBS)      = ( B  S ) list
      lblQP     = iw(271) ! blQP(nb)      = working QP lower bounds
      lbuQP     = iw(272) ! buQP(nb)      = working QP upper bounds
      lblBS     = iw(273) ! blBS(mBS)     = lower bounds for xBS
      lbuBS     = iw(274) ! buBS(mBS)     = upper bounds for xBS
      lpBS      = iw(277) ! pBS(nb)       = search direction
      lxBS      = iw(301) ! xBS(mBS)      = basics, superbasics
      lgQP      = iw(290) ! gQP(ngQP)     = QP gradient
      lgBS      = iw(291) ! gBS(mBS)      = BS components of g
      lrg       = iw(293) ! rg (maxS)     = reduced gradient
      lrg2      = iw(294) ! rg2(maxS)     = reduced gradient copy
      lR        = iw(295) ! R(lenR)       = factor of Z'HZ
      lscales   = iw(296) ! scales(nb)    = row and column scales
      liy       = iw(308) ! iy (nb)       = integer work vector
      liy1      = iw(309) ! iy1(nb)       = integer work vector
      liy2      = iw(310) ! iy2(nb)       = integer work vector
      ly        = iw(311) !  y (nb)       = real work vector
      ly1       = iw(312) !  y1(nb)       = real work vector
      ly2       = iw(313) !  y2(nb)       = real work vector
      lHdx      = iw(288) ! Hdx(nnH)      = product of H with  x - x0
      lblSave   = iw(275) ! blSave        = bl for the basis-finding LP
      lbuSave   = iw(276) ! buSave        = bu for the basis-finding LP

      iExit     = 0

      mBS       = m + maxS

      ! Figure out what type of objective we have.

      if (minmax .eq. 0  .or.(eMode .eq. 2  .and.  lvlObjE .eq. 2)) then
         probType = FP
      else if (ngQP .eq. 0) then ! No explicit objective. Must be an LP.
         if (iObj .eq. 0) then
            probType = FP
         else
            probType = LP
         end if
      else !  Explicit objective. Check for quadratic term.
         if (nnH .gt. 0) then
            probType = QP
         else
            probType = LP
         end if
      end if

      FPonly     = probType .eq. FP
      probTag    = 'linear constraints'

      iw(mnrHdP) = 0            ! Print the header for the Print   file
      iw(mnrHdS) = 0            ! Print the header for the summary file
      iw(linesL) = 0            ! Line count for the print   file
      iw(linesS) = 0            ! Line count for the summary file

      printLevel = mnrPrint

!     Initialize counters based on gotHes and gotFac (set in s3prtQ)

      if (iw(gotFac) .le. 0) then
         iw(nFac) = 0
      end if

      if (iw(gotHes) .le. 0) then
         iw(HvCalls) = 0
      end if

      iw(cgItns)  = 0
      iw(cgItn )  = 0

      itn         = 0
      itQP        = 0
      nDegen      = 0
      nInf        = 0
      nInfE       = 0
      nnCon       = 0
      nnCon0      = 1
      nnJac       = 0
      numLC       = m
      ngQP0       = max( ngQP , 1 )

      iw(QPmode)  = QPsolver    ! Local value of QPslvr
      objQP       = zero
      sInf        = zero
      sInfE       = zero
      scaleObj    = one
      signObj     = minimize
      targetH     = Hcondbnd
      condZmax    = condZmax0   ! First bound on the condition est of Z
      wtInf       = wtInf0

      subOptimize = NO          ! No suboptimization
      itQPtarget  = itnlim      ! No suboptimization
      SwitchToQN  = .false.

!     Start recording the solve time.

      call s1time( 2, 0, iw, leniw, rw, lenrw )

!     Initialize quantities to avoid them being used before being set.

      call dload ( m    , zero, pi         , 1 )
      call dload ( ngQP0, zero, rw(lgQP)   , 1 )
      call iload ( nb   ,    0, iw(leState), 1 )

!     Make copies of the upper and lower bounds.

      call dcopy ( nb, bl, 1, rw(lblQP), 1 )
      call dcopy ( nb, bu, 1, rw(lbuQP), 1 )

      !-----------------------------------------------------------------
      ! Print the matrix statistics.
      ! Find the rowtypes for use in s5getB (they are held in iy2).
      !-----------------------------------------------------------------
      call s2Amat
     &   ( Stats, mjrPrint, m, n, nb,
     &     nnCon, nnJac, nnObj, iObj, numLC, numLIQ,
     &     neA, nlocA, locA, indA, Acol,
     &     rw(lblQP), rw(lbuQP), iw(liy2),
     &     iw, leniw, rw, lenrw )

      !=================================================================
      ! Find a basis kBS(1:m) for the linear constraints and bounds.
      !=================================================================
      ! s5getB does the following.
      !  1. The linear constraints are (optionally) scaled.
      !  2. Elements x(n+1:n+m) of the initial x are assigned.
      !  3. An LP is used to find a feasible x for the bounds and
      !     linear equality constraints.
      !  The base point x0 is not touched.

      call s5getB
     &   ( inform, startType, qpLog, NeedB, m, maxS, mBS,
     &     n, nb, nnCon, nnJac, nnObj, nNames, nS,
     &     itQP, itnlim, itn,
     &     nDegen, numLC, numLIQ, tolOptFP, tolOptQP, tolx,
     &     nInf, sInf, wtInf,
     &     iObj, scaleObj, piNorm, rgNorm,
     &     neA, nlocA, locA, indA, Acol,
     &     eType, iw(leState), iw(liy2), iw(lfeasType),
     &     hs, iw(lkBS), Names,
     &     bl, bu, rw(lblQP), rw(lbuQP),
     &     rw(lblBS), rw(lbuBS), rw(lblSave), rw(lbuSave),
     &     rw(lgBS), pi, rc, nrhs0, nrhs, rhs, rw(lscales),
     &     1, 0, x0, x, rw(lxBS),
     &     iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2),
     &     cw, lencw, iw, leniw, rw, lenrw )

!     Potential inform values are:
!        0   basis found.
!            nInf = 0 => linear equalities are     satisfied.
!            nInf > 0 +> linear equalities are not satisfied
!                        unbounded FP problem
!       >0   fatal error. No feasible point for the equalities

      if (inform .gt. 0) then
         iExit = inform         ! fatal error
         go to 900
      end if

      if (ngObj .gt. 0  .and.  iw(lvlScale) .gt. 0) then
         call ddscl ( ngObj, rw(lscales), 1, gObj, 1 )
      end if

      if (numLIQ .gt. 0) then
         NeedLU = .true.
      else
         NeedLU = .false.
      end if
      Needx     = NeedLU
      typeLU    = BT

!     ------------------------------------------------------------------
!     s5getB has already found a basis for the linear constraints that
!     may or may not be feasible.
!     ------------------------------------------------------------------
      Elastic  = eMode .eq. 2
      elastics = 0
      call iload ( nb, 0, iw(leState), 1 )

!     ==================================================================
!     Solve the problem.
!     ==================================================================
      UseQP = ( probType .eq. LP  .and.  ngObj .gt. 0)  .or.
     &          probType .eq. QP

      if (UseQP) then

         iw(mnrHdP) = 1         ! Refresh print   heading.
         iw(mnrHdS) = 1         ! Refresh summary heading

!        Unless CG was requested explicitly,  use maxS = maxR
!        with Cholesky and QN. Then switch to CG if necessary.

         if (   iw(QPmode) .eq. QPChol) then
         !-----------------------------------------------------------
         ! Solve the QP.
         !-----------------------------------------------------------
            eigH = SEMDEF
            gotR = .false.

            call s5QP
     &         ( inform,
     &           probType, probTag, subOptimize,
     &           qpLog, Hprod, Hprod1, iw(HvCalls), eigH,
     &           Elastic, GotR, NeedLU, typeLU, Needx,
     &           lenR, m, maxR, mBS, n, nb, nDegen,
     &           ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &           itQP, itnlim, itQPtarget, itn,
     &           eMode, lvlObjE, printLevel,
     &           minimize, iObj, scaleObj, objAdd, objQP,
     &           condZmax, targetH, tolOptFP, tolOptQP, tolx,
     &           nInf, sInf, elastics, nInfE, sInfE, wtInf,
     &           piNorm, rgNorm,
     &           neA, nlocA, locA, indA, Acol,
     &           neH, nlocH, locH, indH, Hcol,
     &           eType,  iw(leState), iw(lfeasType), hs, iw(lkBS),
     &           bl, bu, rw(lblQP), rw(lbuQP), rw(lblBS), rw(lbuBS),
     &           rw(lgBS), gObj, rw(lgQP), rw(lHdx), rw(lpBS),
     &           pi, rw(lR), rc, rw(lrg),
     &           nrhs0, nrhs, rhs,  rw(lscales),
     &           lenx0, nx0, x0, x, rw(lxBS), x, ! xFrozen = x, unused
     &           iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2),
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

         else if (iw(QPmode) .eq. QN) then

            GotR = iw(gotHes) .gt. 0

            call s5QN
     &         ( inform,
     &           probType, probTag, subOptimize,
     &           qpLog, Hprod, Hprod1, iw(HvCalls),
     &           Elastic, GotR, NeedLU, typeLU, Needx,
     &           lenR, m, maxR, mBS, n, nb, nDegen,
     &           ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &           itQP, itnlim, itQPtarget, itn,
     &           eMode, lvlObjE, printLevel,
     &           minimize, iObj, scaleObj, objAdd, objQP,
     &           condZHZ, condZmax, tolOptFP, tolOptQP, tolx,
     &           nInf, sInf, elastics, nInfE, sInfE, wtInf, piNorm,
     &           neA, nlocA, locA, indA, Acol,
     &           neH, nlocH, locH, indH, Hcol,
     &           eType, iw(leState), iw(lfeasType), hs, iw(lkBS),
     &           bl, bu, rw(lblQP), rw(lbuQP), rw(lblBS), rw(lbuBS),
     &           rw(lgBS), gObj, rw(lgQP), rw(lHdx), rw(lpBS), pi,
     &           rw(lR), rc, rw(lrg), rw(lrg2),
     &           nrhs0, nrhs, rhs, rw(lscales),
     &           lenx0, nx0, x0, x, rw(lxBS), x, ! xFrozen = x, unused
     &           iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2),
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
         end if

         if (inform .eq. -5  .and.  maxR .lt. maxS) then

!           Too many superbasics. Switch to a CG solver

            if (     iw(QPmode) .eq. QPChol) then
               iw(preCon) = NO  ! with no preconditioning
               GotR       = .false.
            else if (iw(QPmode) .eq. QN    ) then
               iw(preCon) = YES ! with QN preconditioning
            end if
            iw(QPmode)    = CG  ! Switch to CG

         else if  (inform .eq. -6) then

!           Hessian Indefinite or ill-conditioned Switch to QN solver.

            iw(QPmode)    = QN  ! Switch to QN
            SwitchToQN    = .true.
            GotR          = .true.

         end if

         if (iw(QPmode) .eq. CG  .or.  SwitchToQN) then
            if (QPsolver  .eq. CG) GotR = .false.
            call s5QN
     &         ( inform,
     &           probType, probTag, subOptimize,
     &           qpLog, Hprod, Hprod1, iw(HvCalls),
     &           Elastic, GotR, NeedLU, typeLU, Needx,
     &           lenR, m, maxS, mBS, n, nb, nDegen,
     &           ngQP0, ngQP, ngObj0, ngObj, nnH0, nnH, nS,
     &           itQP, itnlim, itQPtarget, itn,
     &           eMode, lvlObjE, printLevel,
     &           minimize, iObj, scaleObj, objAdd, objQP,
     &           condZHZ, condZmax, tolOptFP, tolOptQP, tolx,
     &           nInf, sInf, elastics, nInfE, sInfE, wtInf, piNorm,
     &           neA, nlocA, locA, indA, Acol,
     &           neH, nlocH, locH, indH, Hcol,
     &           eType, iw(leState), iw(lfeasType), hs, iw(lkBS),
     &           bl, bu, rw(lblQP), rw(lbuQP), rw(lblBS), rw(lbuBS),
     &           rw(lgBS), gObj, rw(lgQP), rw(lHdx), rw(lpBS), pi,
     &           rw(lR), rc, rw(lrg), rw(lrg2),
     &           nrhs0, nrhs, rhs, rw(lscales),
     &           lenx0, nx0, x0, x, rw(lxBS), x, ! xFrozen = x, unused
     &           iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2),
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

!           If the QP is unbounded and the reduced Hessian was flagged
!           as indefinite in s5QP, the problem is likely unbounded.

            if (inform .eq. -2  .and.  SwitchToQN) then
               inform = -6
            end if
         end if

      else
         !--------------------------------------------------------------
         ! LP with objective row in A.
         !--------------------------------------------------------------
         nS         = 0         ! Local value
         nnH        = 0         ! Local value
         call s5fixS
     &      ( Fix, m, maxS, mBS, n, nb, nS, hs, iw(lkBS),
     &        rw(lblQP), rw(lbuQP), rw(lblBS), rw(lbuBS), x, rw(lxBS) )

         iw(mnrHdP) = 0
         iw(mnrHdS) = 0
         call s5LP
     &      ( inform,
     &        probType, probTag, Elastic,
     &        subOptimize, qpLog, NeedLU, Needx,
     &        m, n, nb, nDegen, itQP, itnlim, itn,
     &        eMode, lvlObjE, printLevel,
     &        minimize, iObj, scaleObj, objAdd,
     &        condZmax, tolOptFP, tolOptQP, tolx,
     &        nInf, sInf, elastics, nInfE, sInfE, wtInf,
     &        piNorm, rgNorm,
     &        neA, nlocA, locA, indA, Acol,
     &        eType, iw(leState), iw(lfeasType), hs, iw(lkBS),
     &        bl, bu, rw(lblQP), rw(lbuQP), rw(lblBS), rw(lbuBS),
     &        rw(lgBS), pi, rc,
     &        nrhs0, nrhs, rhs, rw(lscales),
     &        x, rw(lxBS), x,   ! xFrozen = x, unused
     &        iw(liy), iw(liy1), rw(ly), rw(ly1),
     &        cw, lencw, iw, leniw, rw, lenrw )
      end if

      Terminated = inform .gt.   0 ! Fatal error (e.g., LU, time limit)
      Optimal    = inform .eq.   0 ! Optimal multipliers
      Infeasible = inform .eq.  -1 ! infeas nonelastics in elastic mode
      Unbounded  = inform .eq.  -2 ! LP is unbounded
      BigItn     = inform .eq.  -3 ! Too many iterations
      WeakMin    = inform .eq.  -4 ! Weak QP minimizer
      BigSB      = inform .eq.  -5 ! Too many superbasics
      IndefQP    = inform .eq.  -6 .or.
     &             inform .eq.  -9 ! QP Hessian not positive semidef
      BadZg      = inform .eq.  -7 ! Z'g could not be made small enough
      BadZ       = inform .eq.  -8 ! Ill-conditioned Z
      SubItns    = inform .eq. -10 ! Too many subiterations.

      ! Many possible outcomes!

      if (Terminated) then
         ! LU error or time limit

         iExit = inform

      else if (Optimal  .and.  nInf .gt. 0) then
         ! Optimal multipliers with infeasible constraints.

         if (Elastic) then
            iExit = 16          ! infeasible nonelastics
         else
            iExit = 11          ! infeasible linear constraints
         end if

      else if (Optimal  .and.  nInf .eq. 0) then
         ! Optimal multipliers with feasible nonelastic constraints.

         if (.not. Elastic) then
            ! Optimal multipliers in normal mode.

            iExit = 1           ! Optimal

         else if (nInfE .gt. 0) then
            ! Elastic mode with some infeasible elastic constraints.

            if      (lvlObjE .eq. 1) then
               iExit = 5        ! elastic objective minimized
            else if (lvlObjE .eq. 2) then
               if (FPonly) then
                  iExit =  6    ! elastic infeasibilities minimized
               else
                  iExit = 14    ! infeasibilities minimized
               end if
            end if
         else
            ! Elastic mode with feasible elastics and nonelastics.

            if (FPonly) then
               iExit = 2        ! Feasible point found.
            else
               iExit = 1        ! Optimal
            end if
         end if

      else if (Infeasible) then
         iExit = 16             ! infeasible nonelastics

      else if (Unbounded ) then
         iExit = 21             ! unbounded

      else if (BigItn    ) then
         iExit = 31             ! too many iterations

      else if (WeakMin   ) then
         iExit =  4             ! weak minimizer

      else if (BigSB     ) then
         iExit = 33             ! too many superbasics

      else if (IndefQP     ) then
         iExit = 53             ! QP Hessian is indefinite

      else
         iExit = 41             ! Current point cannot be improved
      end if

!     ==================================================================
!     Exit.
!     Set output variables and print a summary of the final solution.
!     objTrue is printed in s4newB
!     ==================================================================
  900 call snWRAP( iExit, Solver, str, str2, iw, leniw )

      call s1time(-2, 0, iw, leniw, rw, lenrw )

      degen = 100.0d+0 * nDegen / max( itn, 1 )

!     Print statistics.

      objTrue   = zero

      if (iObj .eq. 0) then
         objLP = objAdd
      else
         objLP = objAdd + x(n+iObj)*scaleObj
      end if

      Infeasible = nInf  .gt. 0
      xNorm      = dnormi( n , x, 1 )

      objMin     = zero
      objInf     = zero

      if (Elastic) then
         objTrue  = objLP
         if (ngQP .gt. 0) then
            objTrue = objTrue + objQP
         end if
         objMin   = signObj*objTrue + wtInf*sInfE
         objInf   =                         sInfE
      else ! Normal mode
         objTrue = objLP
         if (ngQP .gt. 0) then
            objTrue = objTrue + objQP
         end if
         objMin  = signObj*objTrue
      end if

      ! Count basic nonlinear variables (used only for printing).

      nnb   = 0
      do j  = 1, nnH
         if (hs(j) .eq. 3) nnb = nnb + 1
      end do

      if (inewB .gt. 0  .and.  iExit/10 .lt. 8) then
         k      = 1 + iExit/10
         call s4stat
     &      ( k, istate )
         call s4newB
     &      ( Wrap, iNewB, minimize, m, n, nb,
     &        nS, mBS, itn, nInf, sInf, objTrue, iw(lkBS), hs,
     &        rw(lscales), bl, bu, x, rw(lxBS), istate,
     &        cw, lencw, iw, leniw )
      end if

!     ------------------------------------------------------------------
!     Print statistics.
!     ------------------------------------------------------------------
      call snPRNT(13,
     &     ' Problem name                 '//mProb, iw, leniw )

      write(str, 1900) itn, objTrue
      call snPRNT( 3, str, iw, leniw )
      if (Infeasible) then
         write(str, 1910) nInf, sInf
         call snPRNT( 3, str, iw, leniw )
      end if
      if (probType .eq. QP) then
         write(str, 1920) iw(HvCalls), ObjLP
         call snPRNT( 3, str, iw, leniw )
         write(str, 1930) objQP
         call snPRNT( 3, str, iw, leniw )
      end if
      if (nS .gt. 0) then
         write(str, 1940) nS, nnb
         call snPRNT( 3, str, iw, leniw )
      end if
      if (Elastic) then
         if (nInfE .gt. 0) then
            write(str, 1950) wtInf, objMin
            call snPRNT( 3, str, iw, leniw )
            write(str, 1960) nInfE, objInf
            call snPRNT( 3, str, iw, leniw )
         end if
      else                      ! Normal mode
         if (iw(cgItns) .gt. 0) then
            write(str, 1970) iw(cgItns)
            call snPRNT( 3, str, iw, leniw )
         end if
      end if

      write(str, 1980) nDegen, degen
      call snPRNT( 3, str, iw, leniw )

!     ------------------------------------------------------------------
!     Unscale, save basis files and prepare to print the solution.
!     Clock 3 is "Output time".
!     ------------------------------------------------------------------
      call s1time( 3, 0, iw, leniw, rw, lenrw )
      call s4saveB
     &   ( inform, SaveB, minimize, m, n, nb, nkx,
     &     nnCon0, nnCon, ngQP0, ngQP, nNames, nS,
     &     itn, nInf, sInf, wtInf, vimax, iObj, scaleObj, objTrue,
     &     pNorm1, pNorm2, piNorm, xNorm,
     &     neA, nlocA, locA, indA, Acol, iw(lkx),
     &     iw(leState), hs, rw(lscales), bl, bu, Fx, rw(lgQP),
     &     Names, pi, rc, x,
     &     cw, lencw, iw, leniw, rw, lenrw )

      if (ngObj .gt. 0  .and.  iw(lvlScale) .gt. 0)
     &     call dddiv ( ngObj, rw(lscales), 1, gObj, 1 )

!     If task = PrintS, s4saveB prints the solution under the control
!     of lprSol (set by the  Solution  keyword in the SPECS file).
!     The printed solution may or may not be wanted, as follows:
!
!     lprSol = 0   means      No
!            = 1   means      If optimal, infeasible or unbounded
!            = 2   means      Yes
!            = 3   means      If error condition

      call s4saveB
     &   ( inform, PrintS, minimize, m, n, nb, nkx,
     &     nnCon0, nnCon, ngQP0, ngQP, nNames, nS,
     &     itn, nInf, sInf, wtInf, vimax,
     &     iObj, scaleObj, objTrue,
     &     pNorm1, pNorm2, piNorm, xNorm,
     &     neA, nlocA, locA, indA, Acol, iw(lkx),
     &     iw(leState), hs, rw(lscales), bl, bu, Fx, rw(lgQP),
     &     Names, pi, rc, x,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call s1time(-3, 0, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Set Obj for output.
!     Call  Hx  one last time with  nState .ge. 2.
!     Everything has been  unscaled, so we have to disable scaling.
!     ------------------------------------------------------------------
      lsSave       = iw(lvlScale)
      iw(lvlScale) = 0
      sqStat       = 2 + min( iExit/10,4 )
      iw(qpStat)   = sqStat

      ObjLP = 0
      if (probType .eq. FP) then
         objTrue = zero
      else if (probType .eq. LP  .or.  probType .eq. QP) then
         objTrue = objAdd

         if (iObj .gt. 0) then
            ObjLP   = x(n+iObj)*scaleObj
            objTrue = objTrue + ObjLP
         end if

         if (ngQP .gt. 0) then
            call s5QPfg
     &         ( Hprod, Hprod1, iw(HvCalls),
     &           ngQP, ngObj0, ngObj, nnH,
     &           neH, nlocH, locH, indH, Hcol,
     &           sqStat, objQP,
     &           gObj, rw(lgQP), lenx0, nx0, x0, x, rw(ly),
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            objTrue = objTrue + objQP
         end if
      end if

      iw(lvlScale) = lsSave

!     Save some things needed by solvers calling s5Solve

      rw(421) = objTrue         ! The true objective
      rw(422) = piNorm          ! Lagrange multiplier norm
      rw(423) = xNorm           ! Norm of the variables (for GAMS)
      rw(424) = wtInf           ! Infeasibility weight

      rw(433) = sInf + sInfE    ! Sum of infeasibilities
      rw(434) = objLP           ! Linear    objective term

      iw(421) = itn             ! Total iteration count
      iw(423) = maxS            ! max # of superbasics

!     No nonlinear constraints.

      rw(432) = zero            ! Inf norm of the constraint violation
      rw(435) = zero            ! Nonlinear objective term
      rw(436) = zero            ! Norm of penalty parameters

      iw(422) = 0               ! Major iterations

      return

 1900 format(
     &     ' No. of iterations', i20, 2x,
     &     ' Objective', 6x, 1p, e22.10)
 1910 format(
     &     ' No. of infeasibilities', i15, 2x,
     &     ' Sum of infeas', 1p, e24.10)
 1920 format(
     &     ' No. of Hessian products', i14, 2x,
     &     ' Linear    objective', 1p, e18.10)
 1930 format(
     &       40x,
     &     ' Quadratic objective', 1p, e18.10)
 1940 format(
     &     ' No. of superbasics', i19, 2x,
     &     ' No. of basic nonlinears', i14)
 1950 format(
     &     ' Elastic weight            ', 1p, e11.1, 2x,
     &     ' Elastic objective',      1p, e20.10)
 1960 format(
     &     ' No. of infeas elastics', i15, 2x,
     &     ' Elastic infeas   ', 1p, e20.10)
 1970 format(
     &     ' No. of CG iterations', i17)
 1980 format(
     &     ' No. of degenerate steps', i14, 2x,
     &     ' Percentage', f27.2)

      end ! subroutine s5solve

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5defaults
     &   ( m, n, lencObj, ncolH, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     lencObj, lencw, leniw, lenrw, m, n, ncolH, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     s5defaults checks and possibly prints the optional parameter values
!     for sqopt.
!
!     Optional parameters are checked and, if necessary,  changed to
!     reasonable values.
!
!     Note that parameters are checked before the amount of working
!     storage has been defined.
!
!     See  snworkspace.info  for full documentation of cw, iw and rw.
!
!     15 Nov 1991: first version.
!     02 Aug 2003: snPRNT adopted.
!     22 Jun 2004: Added default LU mod singularity tol
!     21 Dec 2004: Default LU tols fixed up.
!     02 May 2006: lvlTim removed.
!     01 Sep 2007: stickyOp added.
!     18 Jan 2010: PreCon initialized.
!     13 Jul 2013: Default lvldif set correctly.
!     25 Oct 2014: Bogus assignment of LUprnt removed (its set in s2BLU).
!     25 Oct 2014: nout set independently of lvlSys.
!     27 Dec 2014: tolOptFP set independently of tolOptQP.
!     ==================================================================
      logical
     &     linear, QP
      integer
     &     cgItmx, iCrash, iBack, iDump, iLoadB, iMPS, iNewB, iInsrt,
     &     iOldB, iPnch, iPrint, iReprt, iSoln, itnlim, kchk,
     &     kDegen, kFac, klog, kReset, ksav, kSumm, eMode, lprDbg,
     &     lprPrm, lprScl, lprSol, lvlObjE, lvlPre, lvlPiv,
     &     lvlScale, lvlSys, maxmn, maxCol, maxR, maxS, mflush,
     &     minimize, minmax, minPrc, mMinor, mjrPrint, mnrPrint,
     &     mSkip, mNewSB, never, nout, nParPrU, nParPrLP, nParPrQP,
     &     nPr1, nPr2, Precon, QPslvr, stickyOp, TPivot
      double precision
     &     bigdx, bigFx, c4, c6, chzbnd, condZmax0, Dens1, Dens2,
     &     Lmax1, Lmax2, eps, eps0, eps1, eps2, eps3, eps4, etarg,
     &     Hcondbnd, infBnd, maxTime, rmaxS, rtcondZbnd, scltol, small,
     &     tCrash, tolCG, tolCon, tolDcp, tolDdp, tolDpp, tolDrp,
     &     tolDup, toldj3, tolFac, tolOptFP, tolNLP, tolpiv, tolOptQP,
     &     tolRow, tolSwp, tolUpd, tolx, Uspace, Utol1, Utol1m, Utol2,
     &     Utol2m, wtInf0, xdlim
!     ------------------------------------------------------------------
      integer            QPChol,     CG
      parameter         (QPChol = 0, CG = 1)
      integer            idummy
      parameter         (idummy = -11111)
      double precision   zero,             one
      parameter         (zero   =  0.0d+0, one    = 1.0d+0)
      double precision   ten
      parameter         (ten    = 10.0d+0)
      double precision   tenp6,            hundrd
      parameter         (tenp6  = 1.0d+6,  hundrd = 100.0d+0)
!     ------------------------------------------------------------------
!     Set some local machine-dependent constants.

      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      eps0      = rw(  2) ! eps**(4/5)          IEEE DP  3.00e-13
      eps1      = rw(  3) ! eps**(2/3)          IEEE DP  3.67e-11
      eps2      = rw(  4) ! eps**(1/2)          IEEE DP  1.49e-08
      eps3      = rw(  5) ! eps**(1/3)          IEEE DP  6.05e-06
      eps4      = rw(  6) ! eps**(1/4)          IEEE DP  1.22e-04

!     ------------------------------------------------------------------
!     rw(51)--rw(150): optional parameters set via the specs file.
!     ------------------------------------------------------------------
      tolOptFP  = rw( 51) ! Minor Phase 1 Opt tol
      tolOptQP  = rw( 52) ! Minor Phase 2 Opt tol
      tolCG     = rw( 54) ! cg tolerance

      tolx      = rw( 56) ! Minor feasibility tolerance

      tolpiv    = rw( 60) ! excludes small elements of y
      tolrow    = rw( 61) ! tolerance for the row error
      tCrash    = rw( 62) ! crash tolerance
      Utol1m    = rw( 63) ! abs tol for small diag of U in LU mod
      Utol2m    = rw( 64) ! rel tol for small diag of U in LU mod
      tolswp    = rw( 65) ! LU swap tolerance
      tolFac    = rw( 66) ! LU factor tolerance
      tolUpd    = rw( 67) ! LU update tolerance
      infBnd    = rw( 70) ! definition of plus infinity
      bigFx     = rw( 71) ! unbounded objective
      bigdx     = rw( 72) ! unbounded step
      lvlPre    = iw( 77) ! >0    => QN preconditioned CG
      maxTime   = rw( 79) ! time limit
      xdlim     = rw( 80) ! Step limit
      etarg     = rw( 83) ! Quasi-Newton QP rg tolerance
      Hcondbnd  = rw( 85) ! bound on the condition of Hz
      condZmax0 = rw( 86) ! Initial bound on the condition est of Z
      wtInf0    = rw( 88) ! infeasibility weight

      scltol    = rw( 92) ! scale tolerance.
!     ------------------------------------------------------------------
!     rw(151)--rw(180) contain  parmLU  parameters for LUSOL.
!     ------------------------------------------------------------------
      Lmax1     = rw(151) ! max L-multiplier in factor
      Lmax2     = rw(152) ! max L-multiplier in update
      small     = rw(153) ! defn of small real
      Utol1     = rw(154) ! abs tol for small diag of U
      Utol2     = rw(155) ! rel tol for small diag of U
      Uspace    = rw(156) ! limit on waste space in U
      Dens1     = rw(157) ! switch to search maxcol columns and no rows
      Dens2     = rw(158) ! switch to dense LU
!     ------------------------------------------------------------------
!     rw(181)--rw(199) pass parameters into various routines.
!     ------------------------------------------------------------------
!     toldj3    = rw(186) ! current optimality tol
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     iw(1)--iw(50): I/O file numbers and dimensions.
!     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file
!     ------------------------------------------------------------------
!     iw(51)--iw(150): optional parameters set via the specs file.
!     ------------------------------------------------------------------
      maxR      = iw( 52) ! max columns of R.
      maxS      = iw( 53) ! max # of superbasics
      QPslvr    = iw( 55) ! 0(1) => QP(QN) QP solver
      eMode     = iw( 56) ! >0    => use elastic mode
      kchk      = iw( 58) ! check (row) frequency
      kFac      = iw( 59) ! factorization frequency
      ksav      = iw( 60) ! save basis map
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      kDegen    = iw( 63) ! max. expansions of featol
      kReset    = iw( 64) ! Hessian frequency
      mFlush    = iw( 66) ! Hessian flush
      mSkip     = iw( 67) ! # largest value of nSkip
!     lvlStart  = iw( 69) ! = 0:1:2:3 => cold:warm:basis:hot start
      lvlSys    = iw( 71) ! > 0   => print system info
      lvlObjE   = iw( 73) ! Elastic option
      lvlScale  = iw( 75) ! scale option
      lvlPiv    = iw( 80) ! 0(1) LU threshold partial(complete) pivoting
      lprPrm    = iw( 81) ! > 0    => parms are printed
      lprScl    = iw( 83) ! > 0    => print the scales
      lprSol    = iw( 84) ! > 0    => print the solution
      lprDbg    = iw( 85) ! > 0    => private debug print
      minmax    = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
      iCrash    = iw( 88) ! Crash option
      itnlim    = iw( 89) ! limit on total iterations
      mMinor    = iw( 91) ! limit on minor iterations
      mnrPrint  = iw( 93) ! Minor print level
      nParPrU   = iw(101) ! # of partial pricing sections
      mNewSB    = iw( 95) ! # of working set changes
      cgItmx    = iw( 97) ! CG iteration limit
      stickyOp  = iw(116) ! > 0 => optional parameters are sticky
      iBack     = iw(120) ! backup file
      iDump     = iw(121) ! dump file
      iLoadB    = iw(122) ! load file
      iMPS      = iw(123) ! MPS file
      iNewB     = iw(124) ! new basis file
      iInsrt    = iw(125) ! insert file
      iOldB     = iw(126) ! old basis file
      iPnch     = iw(127) ! punch file
      iReprt    = iw(130) ! Report file
      iSoln     = iw(131) ! Solution file
!     ------------------------------------------------------------------
!     iw(151)--iw(180) contain luparm parameters for LUSOL.
!     ------------------------------------------------------------------
      nout      = iw(151) ! unit # for printed messages
      maxcol    = iw(153) ! lu1fac: max. # columns
!     ------------------------------------------------------------------
      PreCon    = iw(209) ! Current precon mode (based on QPslvr)

      c4        = max( 1.0d-4, eps3 )
      c6        = max( 1.0d-6, eps2 )
      never     = 99999999
      QP        = ncolH .gt. 0
      linear    = .not. QP

!     ==================================================================
!     Check the optional parameters.
!     ==================================================================
      if (iBack  .eq. idummy ) iBack  =     0
      if (iDump  .eq. idummy ) iDump  =     0
      if (iLoadB .eq. idummy ) iLoadB =     0
      if (iNewB  .eq. idummy ) iNewB  =     0
      if (iInsrt .eq. idummy ) iInsrt =     0
      if (iOldB  .eq. idummy ) iOldB  =     0
      if (iPnch  .eq. idummy ) iPnch  =     0
      if (iReprt .eq. idummy ) iReprt =     0
      if (iSoln  .eq. idummy ) iSoln  =     0

!     Set unspecified frequencies or silly values to defaults.

      if (kchk   .eq. idummy ) kchk   =    60
      if (kfac   .le.    0   ) then
                               kfac   =   100
                     if ( QP ) kfac   =    50
      end if
      if (klog  .eq. idummy  ) klog   =   100
      if (kSumm .eq. idummy  ) kSumm  =   100
      if (ksav  .eq. idummy  ) ksav   =   100
      if (kDegen.eq. idummy  ) kDegen = 10000

!     Sometimes, frequency 0 means "almost never".

      if (kchk   .le. 0      ) kchk   = never
      if (klog   .le. 0      ) klog   = never
      if (ksav   .le. 0      ) ksav   = never
      if (kSumm  .le. 0      ) kSumm  = never
      if (kDegen .le. 0      ) kDegen = never
      if (iCrash .lt. 0      ) iCrash =  3
      if (minmax .eq. idummy ) minmax =  1
      if (minmax .eq. -1     ) then
                               minimize = -1
      else
                               minimize =  1
      end if

      if (mMinor    .lt. 0     ) mMinor   = max(1000, 5*max(n,m))
      if (mNewSB    .le. 0     ) mNewSB   = never
      if (lprDbg    .lt. 0     ) lprDbg   = 0
      if (lprPrm    .lt. 0     ) lprPrm   = 1
      if (lprScl    .lt. 0     ) lprScl   = 0
      if (lprSol    .lt. 0     ) lprSol   = 2
!     lvlStart is checked in s3argQ
!     if (lvlStart  .lt. 0     ) lvlStart = 0
      if (mnrPrint  .lt. 0     ) mnrPrint = 1
                                 mjrPrint = mnrPrint
      if (lvlObjE   .lt. 0  .or. lvlObjE .gt. 2
     &                         ) lvlObjE  = idummy
      if (lvlObjE   .eq. idummy) lvlObjE  = 2
      if (lvlSys    .lt. 0     ) lvlSys   = 0
      if (eMode     .lt. 0  .or.  eMode .gt. 2
     &                         ) eMode    = idummy
      if (eMode     .eq. idummy) eMode    = 1

      if (stickyOp  .lt. 0     ) stickyOp = 0

!     Check superbasics limit and reduced Hessian size.

      if ( QP ) then
         if (maxR .lt. 0     ) maxR   = min( 2000, ncolH+1 )
         if (maxS .lt. 0     ) maxS   =            ncolH+1
                               maxR   = max( min( maxR ,n ) , 0 )
                               maxS   = max( min( maxS ,n ) , 1 )
      else ! linear
         if (maxS   .le. 0   ) maxS   = 1
                               maxR   = 1
      end if

      if (maxS   .lt. maxR   ) maxS   = maxR

      if (QPslvr .lt. 0      ) QPslvr = QPChol
      if (maxR   .eq. 0      ) QPslvr = CG
      if (lvlPre .lt. 0 .or.
     &    lvlPre .gt. 1      ) then
                               lvlPre = 0
                               PreCon = 0
      else
                               PreCon = 1
      end if
      if (cgItmx .lt. 0      ) cgItmx = 100
      if (etarg  .lt. zero  .or.
     &    etarg  .gt. one    ) etarg  = 0.5d+0

!     Check other options.

      if (lvlScale.lt. 0     ) lvlScale = 0
                               lvlScale = min( lvlScale, 2 )

!     Partial pricing parameters
                               minPrc   = 10
                               maxmn    = max( m, n )

!     Temporarily go back to what we had before.
!     Set nParPrU to previous value.

      if (nParPrU .le. 0) then
         if (QP) then
                               nParPrU = 1
         else
                               nParPrU = 10
         end if
      end if
                               minPrc  = 10
                               nPr1    = n / nParPrU
                               nPr2    = m / nParPrU
      if (max( nPr1, nPr2 ) .lt. minPrc) then
                               maxmn   = max( m, n )
                               nParPrU = maxmn / min( maxmn, minPrc )
      end if

      nParPrLP = nParPrU
      nParPrQP = nParPrU

c$$$!     The new part follows
c$$$      if (nParPrU .le. 0     ) then
c$$$         if (nParPrLP .le. 0 ) nParPrLP = 10
c$$$                               nPr1     = n / nParPrLP
c$$$                               nPr2     = m / nParPrLP
c$$$         if (max( nPr1, nPr2 ) .lt. minPrc) then
c$$$                               nParPrLP = maxmn / min( maxmn, minPrc )
c$$$         end if
c$$$
c$$$         if (nParPrQP .le. 0 ) nParPrQP = 1
c$$$                               nPr1     = n / nParPrQP
c$$$                               nPr2     = m / nParPrQP
c$$$         if (max( nPr1, nPr2 ) .lt. minPrc) then
c$$$                               nParPrQP = maxmn / min( maxmn, minPrc )
c$$$         end if
c$$$      else
c$$$         nParPrLP = nParPrU
c$$$         nParPrQP = nParPrU
c$$$      end if

      if (maxTime  .lt. zero ) maxTime = zero

      rmaxS  = maxS
      cHzbnd = max ( one/(hundrd*eps*rmaxS), tenp6 )

      if (infBnd    .lt. zero ) infBnd     = 1.0d+20
      if (bigFx     .le. zero ) bigFx      = 1.0d+15
      if (bigdx     .le. zero ) bigdx      = infBnd
      if (Hcondbnd  .le. zero ) Hcondbnd   = cHzbnd
                                rtcondZbnd = 1.0d+8
      if (xdlim     .le. zero ) xdlim      = 2.0d+0
      if (condZmax0 .le. zero ) condZmax0  = 1.0d+6

      if (tCrash   .lt. zero  .or.
     &    tCrash   .ge. one   ) tCrash     = 0.1d+0

!     ---------------------------------------
!     Set up the parameters for lu1fac.
!     LUprnt > 0 gives LU output on unit nout.
!     LUprnt is set in s2BLU.
!     ----------------------------------------
      if (maxcol .lt.  0     ) maxcol =   5
      if (nout   .eq.  idummy) nout   = iPrint

      if (lvlPiv .le.  0     ) lvlPiv =  0
      if (lvlPiv .gt.  3     ) lvlPiv =  0
                               TPivot =  lvlPiv
      if (linear) then
                               tolDpp =  hundrd
                               tolDrp =  ten
                               tolDcp =  ten
                               tolDdp =  ten
                               tolDup =  ten
      else ! QP
                               tolDpp =  3.99d+0
                               tolDrp =  3.99d+0
                               tolDcp =  3.99d+0
                               tolDdp =  3.99d+0
                               tolDup =  3.99d+0
      end if
      if (tolFac .lt. one    ) then
         if (lvlPiv .eq.   0 ) tolFac =  tolDpp
         if (lvlPiv .eq.   1 ) tolFac =  tolDrp
         if (lvlPiv .eq.   2 ) tolFac =  tolDcp
         if (lvlPiv .eq.   3 ) tolFac =  tolDdp
      end if
      if (tolUpd    .lt. one ) tolUpd =  tolDup
                               Lmax1  =  tolFac
                               Lmax2  =  tolUpd
      if (Utol1     .le. zero) Utol1  =  eps1
      if (Utol2     .le. zero) Utol2  =  eps1
      if (Utol1m    .le. zero) Utol1m =  eps1
      if (Utol2m    .le. zero) Utol2m =  eps1
      if (Dens2     .lt. zero) Dens2  =  0.6d+0
      if (small     .le. zero) small  =  eps0
      if (Uspace    .le. zero) Uspace =  3.0d+0
      if (Dens1     .le. zero) Dens1  =  0.3d+0

!     Set some tolerances.
!     Set the optimality tolerance.
!     Solve the QP subproblems fairly accurately.

      if (tolCG     .le. zero) tolCG    =  1.0d-2
      if (tolOptQP  .le. zero) tolOptQP =  c6
      if (tolOptFP  .lt. zero) tolOptFP =  c6
      if (tolrow    .le. zero) tolrow   =  c4
      if (tolswp    .le. zero) tolswp   =  eps4
      if (tolx      .le. zero) tolx     =  c6
                               toldj3   =  tolOptQP
      if (scltol    .le. zero) scltol   =  0.90d+0
      if (scltol    .ge. one ) scltol   =  0.99d+0
      if (tolpiv    .le. zero) tolpiv   =  eps1

      if (wtInf0    .lt. zero) wtInf0   =  1.0d+0

      if (iBack     .eq.iNewB) iBack  = 0
      if (itnlim    .lt. 0   ) itnlim = max(10000, 10*max(n,m))

!     Load tolerances used to mark variables during printing in s4SavB.

      tolNLP  = tolOptQP
      tolCon  = tolx

!     ------------------------------------------------------------------
!     Re-assign the options to their respective work arrays.
!     ------------------------------------------------------------------
      rw( 51) = tolOptFP
      rw( 52) = tolOptQP
      rw( 53) = tolNLP
      rw( 54) = tolCG
      rw( 56) = tolx
      rw( 57) = tolCon
      rw( 60) = tolpiv
      rw( 61) = tolrow
      rw( 62) = tCrash
      rw( 65) = tolswp
      rw( 66) = tolFac
      rw( 67) = tolUpd
      rw( 70) = infBnd
      rw( 71) = bigFx
      rw( 72) = bigdx
      rw( 79) = maxTime
      rw( 80) = xdlim
      rw( 83) = etarg
      rw( 85) = Hcondbnd
      rw( 86) = condZmax0
      rw( 88) = wtInf0
      rw( 92) = scltol

      rw(151) = Lmax1  ! max L-multiplier in factor
      rw(152) = Lmax2  ! max L-multiplier in update
      rw(153) = small  ! defn of small real
      rw(154) = Utol1  ! abs tol for small diag of U
      rw(155) = Utol2  ! rel tol for small diag of U
      rw(156) = Uspace ! limit on waste space in U
      rw(157) = Dens1  ! switch to search maxcol columns and no rows
      rw(158) = Dens2  ! switch to dense LU
      rw(181) = tolDpp
      rw(182) = tolDcp
      rw(183) = tolDup
      rw(186) = toldj3
      rw(187) = tolDrp

      iw( 52) = maxR
      iw( 53) = maxS
      iw( 55) = QPslvr
      iw( 56) = eMode
      iw( 58) = kchk
      iw( 59) = kFac
      iw( 60) = ksav
      iw( 61) = klog
      iw( 62) = kSumm
      iw( 63) = kDegen
      iw( 64) = kReset
      iw( 66) = mFlush
      iw( 67) = mSkip
!     iw( 69) = lvlStart
      iw( 71) = lvlSys
      iw( 73) = lvlObjE
      iw( 75) = lvlScale
      iw( 77) = lvlPre
      iw( 80) = lvlPiv
      iw( 81) = lprPrm
      iw( 83) = lprScl
      iw( 84) = lprSol
      iw( 85) = lprDbg
      iw( 87) = minmax
      iw( 88) = iCrash
      iw( 89) = itnlim
      iw( 91) = mMinor
      iw( 92) = mjrPrint
      iw( 93) = mnrPrint
      iw( 95) = mNewSB
      iw( 97) = cgItmx
      iw(116) = stickyOp
      iw(120) = iBack
      iw(121) = iDump
      iw(122) = iLoadB
      iw(123) = iMPS
      iw(124) = iNewB
      iw(125) = iInsrt
      iw(126) = iOldB
      iw(127) = iPnch
      iw(130) = iReprt
      iw(131) = iSoln
      iw(151) = nout
      iw(153) = maxcol
      iw(156) = TPivot
      iw(199) = minimize

      iw(209) = PreCon   ! not optional parameters, but set here.
      iw( 99) = nParPrLP
      iw(100) = nParPrQP
      rw(186) = toldj3
      rw(192) = rtcondZbnd

      end ! subroutine s5defaults

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5Map
     &   ( m, n, nkx, ngObj, nnH,
     &     lenR, maxR, maxS,
     &     nextcw, nextiw, nextrw, iw, leniw )

      implicit
     &     none
      integer
     &     m, n, nkx, ngObj, nnH, lenR, maxR, maxS,
     &     nextcw, nextiw, nextrw, leniw, iw(leniw)

!     ==================================================================
!     s5Map   allocates all array storage for sqopt,
!     using the values:
!        m    , n    , neA
!        maxS                    Set in s5defaults.
!        ngObj, nnH              Set from the argument list.
!        lenR                    Set in the calling program.
!
!     15 Nov 1991: First version based on Minos 5.4 routine m2core.
!     12 Nov 1994: Converted to integer and real storage.
!     06 Aug 1996: First min sum version.
!     14 Jul 1997: Thread-safe version.
!     01 May 1998: First version called by sqMem. This simplified
!                  version may slightly overestimate needed memory.
!     02 Aug 2003: snPRNT adopted.
!     13 May 2005: Bug fix: ly3 assigned to iw correctly
!     18 Jun 2008: Added space for iy2, pBS and rg2.
!     20 Sep 2014: Added space for bl and bu.
!     ==================================================================
      integer
     &     mBS, nb, lblQP, lbuQP, lblBS, lbuBS, lblSave, lbuSave,
     &     ldx, lgBS, lgQP, lHdx, leState, lfeasType, liy, liy1, liy2,
     &     lkBS, lkx, lQPrhs, lpBS, lR, lr1, lr2, lrg, lrg2,
     &     ls1, ls2, ls3, lscales, lxBS, lxScaled, ly, ly1, ly2, ly3,
     &     ngQP
!     ------------------------------------------------------------------
      ngQP    = max( ngObj, nnH )
      mBS     = m     + maxS
      nb      = n     + m

!     sqopt can use all of cw, iw and rw
!     except the first user workspace partitions.

      lkx       = nextiw
      lfeasType = lkx       + nkx
      lkBS      = lfeasType + mBS
      leState   = lkBS      + mBS
      liy       = leState   + nb
      liy1      = liy       + nb
      liy2      = liy1      + nb
      nextiw    = liy2      + nb

!     Addresses for the double precision arrays.

      lscales   = nextrw
      ly        = lscales   + nb
      ly1       = ly        + nb
      ly2       = ly1       + nb
      if (maxR .lt. maxS) then  ! Define SYMMLQ workspace
         ly3    = ly2     + nb
         ls1    = ly3     + nb
         ls2    = ls1     + maxS
         ls3    = ls2     + maxS
         lr1    = ls3     + maxS
         lr2    = lr1     + maxS
         lblQP  = lr2     + maxS
      else
         ly3    = ly2     + nb
         ls1    = ly3
         ls2    = ls1
         ls3    = ls2
         lr1    = ls3
         lr2    = lr1
         lblQP    = lr2
      end if
      lbuQP    = lblQP    + nb
      lblBS    = lbuQP    + nb
      lbuBS    = lblBS    + mBS
      lxBS     = lbuBS    + mBS
      lxScaled = lxBS     + mBS
      lHdx     = lxScaled + nnH
      lpBS     = lHdx     + nnH
      lgQP     = lpBS     + nb
      lgBS     = lgQP     + ngQP
      lR       = lgBS     + mBS
      lrg      = lR       + lenR
      lrg2     = lrg      + maxS
      lblSave  = lrg2     + maxS
      lbuSave  = lblSave  + nb
      lQPrhs   = lbuSave  + nb
      ldx      = lQPrhs   + m
      nextrw   = ldx      + ngQP

!     ---------------------------
!     Store the addresses in iw.
!     ---------------------------
      iw(251) = lkx

      iw(271) = lblQP
      iw(272) = lbuQP
      iw(273) = lblBS
      iw(274) = lbuBS
      iw(275) = lblSave
      iw(276) = lbuSave
      iw(277) = lpBS
      iw(278) = lQPrhs

      iw(284) = lfeasType
      iw(285) = leState

      iw(287) = ldx
      iw(288) = lHdx
      iw(290) = lgQP
      iw(291) = lgBS
      iw(292) = lkBS
      iw(293) = lrg
      iw(294) = lrg2
      iw(295) = lR
      iw(296) = lscales

      iw(301) = lxBS
      iw(302) = lxScaled

      iw(308) = liy
      iw(309) = liy1
      iw(310) = liy2
      iw(311) = ly
      iw(312) = ly1
      iw(313) = ly2
      iw(314) = ly3

      iw(353) = lr1
      iw(354) = lr2
      iw(355) = ls1
      iw(356) = ls2
      iw(357) = ls3

      end ! subroutine s5Map

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s5CallStatus
     &   ( Status, iw, leniw)

      implicit
     &     none

      integer
     &     Status, leniw, iw(leniw)

      !=================================================================
      ! s5CallStatus fetches the call-status for the sqOpt user-defined
      ! matrix-vector product.
      !
      ! 16 Jun 2008: First version of s5CallStatus.
      !=================================================================
      character
     &     str*80
      integer
     &     qpStatus
      !-----------------------------------------------------------------
      parameter         (qpStatus = 235) ! QP user-routine call-status
      !-----------------------------------------------------------------

      if (     iw(qpStatus) .eq. 0) then
         ! Standard call

         Status       =  0
      else if (iw(qpStatus) .lt. 0) then
         ! First call

         Status       =  1
         iw(qpStatus) =  0
      else if (iw(qpStatus) .ge. 2) then
         ! Last orders please

         Status       = iw(qpStatus)
         iw(qpStatus) = -1
      else
         Status       = iw(qpStatus)
         write(str, 9999) Status
         call snPRNT( 3, str, iw, leniw )
      end if

      return

 9999 format(' XXX  user-function call-status not recognized.',
     &       ' Requested status =', i6)

      end ! subroutine s5CallStatus
