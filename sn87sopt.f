!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn87sopt.f
!
!     s8solve
!     s8defaults  s8firstCall
!     s8Map       s8SQP        s8callStatus
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8solve
     &   ( iExit, Solver, startType,
     &     funwrapper, funcon, funobj, userHv,
     &     mjrLog, mnrLog, snSTOP, GotR,
     &     m, n, nb, nnH, nnCon, nnJac, nnObj,
     &     nNames, iObj, objAdd, fObj, objTrue,
     &     nInf, sInf, nInfE, sInfE,
     &     neJ, nlocJ, locJ, indJ, Jcol,
     &     neH, nlocH, locH, indH, Hcol,
     &     bl, bu, Names,
     &     hs, x, pi, rc, majors, nS,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     funwrapper, funcon, funobj, userHv,
     &     mjrLog, mnrLog, snSTOP
      logical
     &     GotR
      integer
     &     iExit, iObj, lencu, lencw, leniu, leniw, lenru, lenrw,
     &     m, n, nb, neJ, neH, nlocJ, nlocH, nInf, nInfE, nNames,
     &     nnCon, nnJac, nnObj, nnH, nS, startType, locJ(nlocJ),
     &     locH(nlocH), indJ(neJ), indH(neH), hs(nb),
     &     iu(leniu), iw(leniw)
      double precision
     &     objAdd, objTrue, fObj, sInf, sInfE, Hcol(neH), Jcol(neJ),
     &     bl(nb), bu(nb), x(nb), pi(m), rc(nb), ru(lenru), rw(lenrw)
      character
     &     Solver*6, Names(nNames)*8, cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s8solve solves the current problem.
!
!     On entry,
!     the specs file has been read,
!     all data items have been loaded (including locJ, indJ, Jcol, ...),
!     and workspace has been allocated.
!
!     On exit,
!     iExit  =  0 if an optimal solution was found,
!            =  1 if the problem was infeasible,
!            =  2 if the problem was unbounded,
!            =  3 if the Iteration limit was exceeded,
!           ge  4 if iterations were terminated by some other
!                 error condition (see the SNOPT user's guide).
!
!     15 Nov 1991: First version based on Minos 5.4 routine misolv.
!     13 Feb 1994: Eliminated "Cycle" options.
!                  Simplified s4getB.
!     12 Nov 1994: Integer workspace added.
!     25 Jul 1996: Sign of the slacks changed.
!     28 Sep 1997: Character workspace added.
!     11 Nov 1997: Backtracking for undefined functions.
!     26 Dec 1997: Dummy Jacobian scaled in feasibility phase.
!     27 Aug 1998: Constant Jacobian elements handled correctly.
!     10 Oct 1998: Objective and constraint gradient checking merged.
!     11 Oct 1998: Facility to combine funobj and funcon added.
!     23 Dec 1999: Suboptimize option added.
!     30 Dec 2000: Housekeeping for first function call now in snwrapup.
!     03 Aug 2003: snEXIT and snPRNT adopted.
!     19 Mar 2006: Fx initialized for s8savB.
!     16 Jun 2008: Call-status implemented correctly.
!     18 Jun 2008: iy2, pBS and rg2 added to workspace.
!     23 Oct 2010: pinorm initialized at 1.0.
!     01 Dec 2012: eState initialized in s8solve.
!     05 Apr 2014: nnH added to argument list.
!     19 Oct 2014: Added user-defined Hessian.
!     09 Mar 2015: Removed (unused) eigH from workspace.
!     ==================================================================
      character
     &     istate(3)*4, mProb*8, str*133, str2*133
      external
     &     dnormi
      logical
     &     Elastic, FeasibleLC, FPonly, GotFuns, NeedB,
     &     Nonlinear, NonlinearCon, NonlinearObj
      integer
     &     cgItn, cgItns, gotG, iCrash, inewB, inform, itn,
     &     itnlim, j, k, lenR, lenx0, lCrash, lsSave, lvlDer, lvlDif,
     &     lvlHess, lvlScale, lvlSrch, lblQP, lbuQP,
     &     lblBS, lbuBS, lblSave, lbuSave, ldyCon, ldx, lFv, lFx,
     &     lfCon, lfCon1, lfCon2, lgCon, lgCon1, lgCon2, lgConU,
     &     lgObj1, lgObj2, lgBS, lgObj, lgQP,
     &     ldg, lHD, lHdx, leType, leState, lfeasType, linesL, linesS,
     &     liy, liy1, liy2, lkBS, lkx, lyCon, lyCon1, lyCon2, llocG,
     &     lpBS, lQPrhs, lR, lrg, lrg2, lscales, lUdx, lx1, lxBS,
     &     lxPen, lxQP, lxQP0, ly, ly1, ly2, maxS, imaxVi, mBS,
     &     minimize, minmax, mjrPrint, mnrPrint, modefg, nDegen, nFac,
     &     nfCon1, nfCon2, nfCon3, nfCon4, nfObj1, nfObj2, nfObj3,
     &     nfObj4, negCon, nkx, nlocG, majors, minors, nnb,
     &     nnH0, nnObj0, nnObj1, nnCon0, nnCon1, npStatus,
     &     nrhs0, nrhs, numLC, numLIQ, nx0
      double precision
     &     Degen, dnormi, duInf, eps0, fMerit, objInf, objLin, objMin,
     &     penParm(4), infBnd, piNorm, pNorm1, pNorm2, rgNorm, signObj,
     &     scaleObj, tolOptFP, tolOptQP, tolx, tCrash, viLim, maxVi,
     &     maxViRel, penNorm, supVi, wtInf, wtInf0, xNorm
!     ------------------------------------------------------------------
      integer            SEMDEF,       POSDEF
      parameter         (SEMDEF   = 0, POSDEF  = 1)
      integer            Scale,        UnScale
      parameter         (Scale    = 0, UnScale = 1)
      integer            RowTyp,       Stats
      parameter         (Rowtyp   = 0, Stats   = 1)
      integer            Wrap
      parameter         (Wrap     = 1)
      integer            SaveB,        PrintS
      parameter         (SaveB    = 0, PrintS  = 1)
      integer            LM,           FM
      parameter         (LM       = 0, FM      = 1)
      integer            xBound,       xMove
      parameter         (xBound   = 0, xMove   = 1)

      parameter         (lvlDer   =  70) ! = 0,1,2,3 or 4, deriv level
      parameter         (lvlHess  =  72) ! 0,1,2  => LM, FM, Newton
      parameter         (lvlScale =  75) ! scale option
      parameter         (minmax   =  87) ! 1, 0, -1  => MIN, FP, MAX
      parameter         (lvlDif   = 181) ! =1(2) for forwd(cntrl) diffs
      parameter         (gotG     = 184) ! > 0 => some exact derivs
      parameter         (nfCon1   = 189) ! number of calls of fCon
      parameter         (nfCon2   = 190) ! number of calls of fCon
      parameter         (nfCon3   = 191) ! number of calls of fCon
      parameter         (nfCon4   = 192) ! number of calls of fCon
      parameter         (nfObj1   = 194) ! number of calls of fObj
      parameter         (nfObj2   = 195) ! number of calls of fObj
      parameter         (nfObj3   = 196) ! number of calls of fObj
      parameter         (nfObj4   = 197) ! number of calls of fObj
      parameter         (nFac     = 210) ! # of LU factorizations
      parameter         (linesL   = 220) ! # lines in log     file
      parameter         (linesS   = 221) ! # lines in summary file
      parameter         (npStatus = 236) ! NP user-routine call-status
      parameter         (cgItns   = 386) ! Number of symmlq iterations
      parameter         (cgItn    = 387) ! symmlq itns for last minor

      double precision   zero,            one,          ten
      parameter         (zero   = 0.0d+0, one = 1.0d+0, ten = 10.0d+0)
!     ------------------------------------------------------------------
      iNewB     = iw(124) ! new basis file

      negCon    = iw( 20) ! # of nonzero elems in J
      lenR      = iw( 28) ! R(lenR) is the reduced Hessian factor
      maxS      = iw( 53) ! max # of superbasics
      lvlSrch   = iw( 76) ! >0     => use derivatives in the line search

      iCrash    = iw( 88) ! Crash option
      itnlim    = iw( 89) ! limit on total iterations

      mjrPrint  = iw( 92) ! Major print level
      mnrPrint  = iw( 93) ! Minor print level
      minimize  = iw(199) ! 1 (-1)    => minimize (maximize)
      nkx       = iw(247) ! dimension of kx and its inverse, kxN

      eps0      = rw(  2)
      tolOptFP  = rw( 51) ! Minor Phase 1 Opt tol
      tolOptQP  = rw( 52) ! Minor Phase 2 Opt tol
      tolx      = rw( 56) ! Minor feasibility tolerance.
      tCrash    = rw( 62) ! crash tolerance.
      infBnd    = rw( 70) ! definition of an infinite bound.
      vilim     = rw( 81) ! violation limit
      wtInf0    = rw( 88) ! infeasibility weight

      mProb     = cw( 51) ! Problem name

!     Addresses

      lkx       = iw(251) ! j  = kx (jN) => col j of Jcol is variable jN
      llocG     = iw(260) ! locG(nnJac+1) = column pointers for indG

      lblQP     = iw(271) ! blQP(nb)    = working lower bounds
      lbuQP     = iw(272) ! buQP(nb)    = working upper bounds
      lblBS     = iw(273) ! blBS(mBS)   = lower bounds for xBS
      lbuBS     = iw(274) ! buBS(mBS)   = upper bounds for xBS
      lblSave   = iw(275) ! blSave(nb)  = bl for LC feasibility
      lbuSave   = iw(276) ! buSave(nb)  = bu for LC feasibility
      lpBS      = iw(277) ! pBS(nb)     = search direction
      lQPrhs    = iw(278) ! QPrhs(nnCon)=  QP constraint rhs

      leType    = iw(283) ! eType(nb) list of elastic vars
      lfeasType = iw(284) ! feasType(mBS), feasibility types
      leState   = iw(285) ! eState(nb), status of elastics

      ldx       = iw(287) ! dx(nb)      = x1 - x
      lHdx      = iw(288) ! Hdx(nnH)    = product of H with  x1 - x
      ldg       = iw(289) ! dg(nnH)     = gradient difference
      lgQP      = iw(290) ! gQP(ngQP)   = QP gradient
      lgBS      = iw(291) ! gBS(mBS)    = BS components of g
      lkBS      = iw(292) ! kBS(mBS), ( B  S ) list
      lrg       = iw(293) ! rg (maxS)   = reduced gradient
      lrg2      = iw(294) ! rg2(maxS)   = reduced gradient
      lR        = iw(295) ! R(lenR)     = factor of Z'HZ
      lscales   = iw(296) ! scales(nb)  = row and column scales
      lgObj     = iw(297) ! gObj(nnObj) = Objective gradient
      lx1       = iw(300) ! x1(nb)      = new x, used to store x0
      lxBS      = iw(301) ! xBS(mBS)    = basics, superbasics
      lxPen     = iw(304) ! xPen(nnCon) = penalty params
      lxQP      = iw(305) ! xQP(nb)     = QP solution
      lxQP0     = iw(306) ! xQP0(nb)    = QP feasible pt.
      liy       = iw(308) ! iy(nb)      =  integer work vector
      liy1      = iw(309) ! iy1(nb)     =  integer work vector
      liy2      = iw(310) ! iy2(nb)     =  integer work vector
      ly        = iw(311) ! y(nb)       =  real work vector
      ly1       = iw(312) ! y1(nb)      =  real work vector
      ly2       = iw(313) ! y2(nb)      =  real work vector
      lfCon     = iw(316) ! fCon (nnCon) constraints at x
      lfCon1    = iw(317) ! fCon1(nnCon) constraints at x1
      lfCon2    = iw(318) ! fCon2(nnCon) work vector
      lgConU    = iw(319) ! record of unknown derivatives and constants
      lgCon     = iw(320) ! gCon (negCon)   constraint gradients at x
      lgCon1    = iw(321) ! gCon1(negCon)   constraint gradients at x1
      lgCon2    = iw(322) ! gCon2(negCon)   work vector
      lgObj1    = iw(324) ! gObj1(nnObj) objective gradients at x1
      lgObj2    = iw(325) ! gObj2(nnObj) work gObj

      lFx       = iw(336) ! Fx (nnCon)  = F(x) + A(linear)x
      lFv       = iw(337) ! Fv          = F(x) + A(linear)x - sN

      lUdx      = iw(345) ! Udx(nnH)      = product of U with dx
      lHD       = iw(347) ! Diagonal of BFGS Hessian
      lyCon     = iw(348) ! yCon (nnCon)  = multipliers for fCon
      lyCon1    = iw(349) ! yCon1(nnCon)  = yCon at x1
      lyCon2    = iw(350) ! yCon2(nnCon)  = work copy of yCon
      ldyCon    = iw(351) ! dyCon(nnCon)  = yCon1 - yCon

      iExit        = 0

      FeasibleLC   = .false.
      GotFuns      = .false.
      FPonly       = iw(minmax) .eq. 0
      signObj      = minimize

      NonlinearCon = nnCon  .gt. 0
      NonlinearObj = nnObj  .gt. 0
      Nonlinear    = nnH    .gt. 0

      nnObj0  = max( nnObj, 1 )
      nnCon0  = max( nnCon, 1 )
      nnH0    = max( nnH  , 1 )
      mBS     = m     + maxS
      nlocG   = nnJac + 1

      numLC   = m - nnCon

!     Initialize yCon from pi.
!     Zap the pi(i) to prevent them being printed without being set.

      if (NonlinearCon) then
         call dcopy ( nnCon, pi, 1, rw(lyCon)  , 1 )
      end if
      if (numLC .gt. 0) then
         call dload ( numLC, zero,  pi(nnCon+1), 1 )
      end if

!     Initialize a few things.

      iw(lvlDif) = 1
      iw(nFac)   = 0
      iw(cgItns) = 0
      iw(cgItn ) = 0
      nInf       = 0
      nInfE      = 0
      wtInf      = one
      piNorm     = one

      Elastic    = .false.

      duInf      = zero
      fMerit     = zero
      fObj       = zero
      sInf       = zero
      sInfE      = zero
      maxVi      = zero
      maxViRel   = zero

      call dload (  4, zero, penParm, 1 )

      iw(linesL) = 0            ! Line count for the print   file
      iw(linesS) = 0            ! Line count for the summary file

      call iload (  4, 0, iw(nfCon1), 1 )
      call iload (  4, 0, iw(nfObj1), 1 )

      itn        = 0
      nDegen     = 0
      majors     = 0
      minors     = 0

!     Start recording the solve time.

      call s1page( 1, iw, leniw )
      call s1time( 2, 0, iw, leniw, rw, lenrw )

!     Make copies of the upper and lower bounds.

      call dcopy ( nb, bl, 1, rw(lblQP), 1 )
      call dcopy ( nb, bu, 1, rw(lbuQP), 1 )

!     ------------------------------------------------------------------
!     Print the matrix statistics before the nonlinear part of Jcol is
!     loaded with random elements.
!     Find the rowtypes for use in s5getB (they are held in iy2).
!     ------------------------------------------------------------------
      call s2Amat
     &   ( Stats, mjrPrint, m, n, nb,
     &     nnCon, nnJac, nnObj, iObj, numLC, numLIQ,
     &     neJ, nlocJ, locJ, indJ, Jcol,
     &     rw(lblQP), rw(lbuQP), iw(liy2),
     &     iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Make a permanent copy in gConU of the constant Jacobian elements
!     stored in J.  Load the nonlinear part of J with random elements.
!     ------------------------------------------------------------------
      if (NonlinearCon) then
         call s8Gcopy
     &      ( nnCon, nnJac,
     &        neJ   , nlocJ,     locJ, indJ,
     &        neJ   , nlocJ,     locJ,       Jcol,
     &        negCon, nlocG, iw(llocG), rw(lgConU) )
         call s8rand
     &      ( negCon, negCon, rw(lgCon) )
         call s8Gcopy
     &      ( nnCon, nnJac, neJ, nlocJ, locJ, indJ,
     &        negCon, nlocG, iw(llocG), rw(lgCon),
     &        neJ   , nlocJ,     locJ ,     Jcol )
      end if

      !=================================================================
      ! Find a basis kBS(1:m) for the linear constraints and bounds.
      !=================================================================
      ! s5getB does the following.
      !  1. The linear constraints are (optionally) scaled.
      !  2. The bounds bl and bu on the nonlinear rows are relaxed.
      !  3. Elements x(n+1:n+m) of the initial x are assigned.
      !  4. An LP is used to find a feasible x for the bounds and
      !     linear equality constraints.
      !  5. x(nb) is (optionally) scaled and saved in x1.
      !
      !  The linear constraints are unscaled in s8getFeasLC after a
      !  feasible point for the linear constraints has been found.

      nrhs  = 0                 ! No QP rhs when finding the first basis
      nrhs0 = 1
      nx0   = nb                ! elements in x1(nb), the base point
      lenx0 = nb

      call iload ( nb, 0, iw(leType) , 1 ) ! placeholder for s5getB
      call iload ( nb, 0, iw(leState), 1 )

      call s5getB
     &   ( inform,
     &     startType, mnrLog, NeedB, m, maxS, mBS,
     &     n, nb, nnCon, nnJac, nnObj, nNames, nS,
     &     minors, itnlim, itn,
     &     nDegen, numLC, numLIQ, tolOptFP, tolOptQP, tolx,
     &     nInf, sInf, wtInf,
     &     iObj, scaleObj, piNorm, rgNorm,
     &     neJ, nlocJ, locJ, indJ, Jcol,
     &     iw(leType), iw(leState), iw(liy2), iw(lfeasType),
     &     hs, iw(lkBS), Names,
     &     bl, bu, rw(lblQP), rw(lbuQP), rw(lblBS), rw(lbuBS),
     &     rw(lblSave), rw(lbuSave),
     &     rw(lgBS), pi, rc, nrhs0, nrhs, rw(lQPrhs), rw(lscales),
     &     lenx0, nx0, rw(lx1), x, rw(lxBS),
     &     iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2),
     &     cw, lencw, iw, leniw, rw, lenrw )

!     Potential inform values are:
!        0   basis found.
!            nInf = 0 => linear equalities are     satisfied.
!            nInf > 0 +> linear equalities are not satisfied and the
!                        FP problem is unbounded.
!       >0   fatal error. No feasible point

      if (inform .gt. 0) then
         iExit = inform         ! fatal error
      end if

      !=================================================================
      ! An initial basis has been assigned to hs. The array kBS is
      ! defined if s5LP was used to find a feasible point for the linear
      ! equality constraints. Otherwise kBS is set after the first basis
      ! factorization.
      !
      ! Find a feasible point for all the linear constraints.
      ! The norm of x is minimized via a proximal-point QP.
      ! If there is no feasible point, the linear rows can be elastic.
      !=================================================================
      if (numLC .gt. 0) then
         if (iExit .eq. 0) then ! the E rows are feasible, check LG rows
            call s8getFeasLC
     &         ( iExit,
     &           startType, mnrLog, lenR, m, maxS, mBS,
     &           n, nb, nnCon0, nnCon, nnH0, nnH, nDegen, nS,
     &           numLC, numLIQ, itn, itnlim, minors, mnrPrint,
     &           scaleObj, tolOptQP, tolx,
     &           nInf, sInf, nInfE, sInfE, wtInf, piNorm, rgNorm,
     &           neJ, nlocJ, locJ, indJ, Jcol,
     &           neH, nlocH, locH, indH, Hcol,
     &           iw(leType), iw(leState), iw(lfeasType), hs, iw(lkBS),
     &           bl, bu, rw(lblQP), rw(lbuQP), rw(lblBS), rw(lbuBS),
     &           rw(lblSave), rw(lbuSave),
     &           rw(lgBS), rw(lgQP), rw(lHdx), rw(lpBS), pi,
     &           rw(lR), rc, rw(lrg), rw(lQPrhs),  rw(lscales),
     &           rw(lx1), x, rw(lxBS),
     &           iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2),
     &           cw, lencw, iw, leniw, rw, lenrw )
         end if

         !--------------------------------------------------------------
         ! Do some housekeeping in case we have to exit
         !--------------------------------------------------------------
         ! Reinstate the scaled bounds on the nonlinear rows.

         if (NonlinearCon) then
            call dcopy ( nnCon, rw(lblSave+n), 1, bl(n+1), 1 )
            call dcopy ( nnCon, rw(lbuSave+n), 1, bu(n+1), 1 )
         end if ! NonlinearCon

         ! Unscale any linear constraints that were scaled in s5getB.

         if (iw(lvlScale) .gt. 0) then
            call s2applyScales
     &         ( UnScale, m, n, nb, iObj, infBnd, scaleObj,
     &           neJ, nlocJ, locJ, indJ, Jcol,
     &           rw(lscales), bl, bu, pi, x )
         end if
      end if ! numLC > 0

!     Exit if the linear constraints are infeasible.

      FeasibleLC = iExit .eq. 0
      if (.not. FeasibleLC) go to 900

!     ------------------------------------------------------------------
!     Copy the constant Jacobian elements into gCon, gCon1 and gCon2.
!     Reset eType so that only nonlinear rows are elastic.
!     Make sure variables are not outside their bounds
!     (in particular, check the nonlinear slacks).
!     ------------------------------------------------------------------
      if (NonlinearCon) then
         call dcopy ( negCon, rw(lgConU), 1, rw(lgCon ), 1 )
         call dcopy ( negCon, rw(lgConU), 1, rw(lgCon1), 1 )
         call dcopy ( negCon, rw(lgConU), 1, rw(lgCon2), 1 )
         call iload ( nnCon, 3, iw(leType+n), 1 )
      end if ! NonlinearCon

      call s5FixX( xBound, 1, nb, tolx, hs, bl, bu, x )

!     ==================================================================
!     ==================================================================
!     The linear constraints have been satisfied!
!     Compute the problem functions at this all-important point.
!     No scaling yet.

!     Compute any missing derivatives.
!     ==================================================================
!     ==================================================================
      if (Nonlinear) then
         call s8firstCall
     &      ( iExit,
     &        FeasibleLC, GotFuns, NonlinearCon, NonlinearObj,
     &        n, nb, nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &        funwrapper, funcon, funobj, userHv,
     &        bl, bu, x, rw(lx1), rw(lyCon),
     &        neJ   , nlocJ, locJ, indJ,
     &        neH   , nlocH, locH, indH, Hcol,
     &        negCon, nlocG, iw(llocG), fObj,
     &        rw(lfCon) , rw(lgCon) , rw(lgObj) ,
     &        rw(lfCon2), rw(lgCon2), rw(lgObj2), rw(ly), rw(ly1),
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )

!        ---------------------------------------------------------------
!        Load the Jacobian gCon in  J.
!        ---------------------------------------------------------------
         if (GotFuns .and. NonlinearCon) then
            call s8Gcopy
     &         ( nnCon, nnJac,
     &           neJ   , nlocJ,     locJ, indJ,
     &           negCon, nlocG, iw(llocG), rw(lgCon),
     &           neJ, nlocJ, locJ, Jcol )
            call dcopy ( nnCon, rw(lfCon), 1, rw(lFx), 1 )
         end if
         if (iExit .ne. 0) go to 900

      end if

!     ==================================================================
!     Scale the problem.
!     ==================================================================
      if (iw(lvlScale) .gt. 0) then
!        ---------------------------------------------------------------
!        Recompute the vector of row types.
!        ---------------------------------------------------------------
         call s2Amat
     &      ( RowTyp, mjrPrint, m, n, nb,
     &        nnCon, nnJac, nnObj, iObj, numLC, numLIQ,
     &        neJ, nlocJ, locJ, indJ, Jcol,
     &        bl, bu, iw(liy2),
     &        iw, leniw, rw, lenrw )
         call s2getScales
     &      ( mjrPrint, m, n, nb, nnH, nnCon, nnJac, iw(liy2),
     &        neJ, nlocJ, locJ, indJ, Jcol,
     &        rw(lscales), bl, bu, rw(ly), rw(ly2),
     &        iw, leniw, rw, lenrw )
         call s2applyScales
     &      ( Scale, m, n, nb, iObj, infBnd, scaleObj,
     &        neJ, nlocJ, locJ, indJ, Jcol,
     &        rw(lscales), bl, bu, pi, x )

!        ---------------------------------------------------------------
!        The objective and constraint functions haven't been scaled yet.
!        Scale the constant elements in gCon1 and gCon2.
!        Don't forget the initial pi.
!        ---------------------------------------------------------------
         if (NonlinearCon) then
            call dddiv
     &         ( nnCon, rw(lscales+n), 1, rw(lfCon), 1 )
            if (iw(gotG) .gt. 0) then
               call s8scaleJ
     &            ( nnCon, nnJac, negCon, n, rw(lscales),
     &              neJ, nlocJ, locJ, indJ, rw(lgCon), rw, lenrw )
               call dcopy ( negCon, rw(lgCon), 1, rw(lgCon1), 1 )
               call dcopy ( negCon, rw(lgCon), 1, rw(lgCon2), 1 )
            end if
            call ddscl ( nnCon, rw(lscales+n), 1, rw(lyCon), 1 )
         end if

         if (NonlinearObj  .and.  iw(gotG) .gt. 0) then
            call s8scaleG
     &         ( nnObj, rw(lscales), rw(lgObj), rw, lenrw )
         end if
      end if ! iw(lvlScale) > 0

!     ==================================================================
!     s8Fx computes the nonlinear constraint values Fx.
!     Copy these into the slacks x(n+i) and make sure they are feasible.
!     Crash uses them to decide which slacks to grab for the basis
!     If any nonbasic nonlinear slacks are close to a bound,
!     move them exactly onto the bound to avoid very small steps.
!     ==================================================================
      if (NonlinearCon) then
         call s8Fx
     &      ( n, nnCon, nnJac, eps0,
     &        neJ, nlocJ, locJ, indJ, Jcol, rw(lfCon), x, rw(lFx) )
         call s2vmax
     &      ( n, nnCon, imaxVi, maxVi, bl, bu, rw(lFx) )
         supVi = vilim*max( ten, maxVi )
         call dcopy ( nnCon, rw(lFx), 1, x(n+1), 1 )

         call s5FixX
     &      ( xMove, n+1, n+nnCon, tolx, hs, bl, bu, x )

!        ===============================================================
!        Crash on the nonlinear rows.
!        hs(*) already defines a basis for the full problem,  but we
!        want to do better by not including all of the slacks.
!        ===============================================================
         if (NeedB) then

!           Load  iy2  with the row types.
!           s2crash uses kBS as workspace.  It may alter x(n+i) for
!           nonlinear slacks.

            call s2Amat
     &         ( RowTyp, mjrPrint, m, n, nb,
     &           nnCon, nnJac, nnObj, iObj, numLC, numLIQ,
     &           neJ, nlocJ, locJ, indJ, Jcol,
     &           bl, bu, iw(liy2),
     &           iw, leniw, rw, lenrw )
            lcrash = 5
            call s2crash
     &         ( lcrash, mjrPrint, m, n, nb, nnCon,
     &           iCrash, tCrash,
     &           neJ, nlocJ, locJ, indJ, Jcol,
     &           iw(lkBS), hs, iw(liy2), bl, bu, x,
     &           iw, leniw, rw, lenrw )
            NeedB = .false.
         end if ! NeedB
      end if ! NonlinearCon

!     ------------------------------------------------------------------
!     Solve the problem.
!     ------------------------------------------------------------------
      call s1page( 1, iw, leniw )

      call s8SQP
     &   ( iExit, funwrapper, funcon, funobj, userHv,
     &     mjrLog, mnrLog, snSTOP,
     &     Elastic, GotR, startType,
     &     itn, lenR, m, maxS, mBS, n, nb, nS,
     &     nnCon0, nnCon, nnObj0, nnObj, nnH0, nnH,
     &     majors, minors, nDegen, duInf,
     &     minimize, iObj, scaleObj, objAdd, fObj, fMerit,
     &     maxVi, maxViRel, supVi,
     &     nInf, sInf, nInfE, sInfE, wtInf0, wtInf,
     &     penParm, piNorm, xNorm,
     &     neJ   , nlocJ, locJ, indJ, Jcol,
     &     neH   , nlocH, locH, indH, Hcol,
     &     negCon, nlocG, iw(llocG),
     &     iw(leType), iw(leState), iw(lfeasType), hs, iw(lkBS),
     &     bl, bu, rw(lblQP), rw(lbuQP), rw(lblBS), rw(lbuBS),
     &     rw(lFv), rw(lFx), rw(lfCon), rw(lgCon), rw(lgObj),
     &     rw(lfCon1), rw(lgCon1), rw(lgObj1),
     &     rw(lfCon2), rw(lgCon2), rw(lgObj2),
     &     rw(lgBS), rw(lgQP), rw(ldyCon), rw(ldx), rw(ldg),
     &     rw(lUdx), rw(lHD), rw(lHdx), rw(lpBS),
     &     rw(lyCon), rw(lyCon1), rw(lyCon2), pi, rw(lQPrhs),
     &     rw(lR), rc, rw(lrg), rw(lrg2), rw(lscales),
     &     x, rw(lx1), rw(lxBS), rw(lxQP0), rw(lxQP), rw(lxPen),
     &     iw(liy), iw(liy1), rw(ly), rw(ly1), rw(ly2),
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     ==================================================================
!     Exit.
!     Set output variables and print a summary of the final solution.
!     ==================================================================
  900 call snWRAP( iExit, Solver, str, str2, iw, leniw )

      call s1time(-2, 0, iw, leniw, rw, lenrw )

!     Set the objective function  (minimized or maximized)
!     objTrue = objAdd + fObj + x(n+iObj)*scaleObj

      if (iObj .eq. 0) then
         objLin = objAdd
      else
         objLin = objAdd + x(n+iObj)*scaleObj
      end if

      objTrue = objLin
      objMin  = zero
      objInf  = zero

      if (GotFuns) then
         if (Elastic) then
            if (FPonly) then
               objMin = wtInf*sInfE
               objInf =       sInfE
            else
               if (nnObj .gt. 0) then
                  objTrue = objTrue + fObj
               end if
               objMin = signObj*objTrue + wtInf*sInfE
               objInf =                         sInfE
            end if
         else ! Normal mode
            if (FPonly) then
!              Relax
            else
               if (nnObj .gt. 0) then
                  objTrue = objTrue + fObj
               end if
               objMin = signObj*objTrue
            end if
         end if
      end if

!     ------------------------------------------------------------------
!     Print statistics.
!     ------------------------------------------------------------------
      degen   = 100.0d+0 * nDegen / max( itn, 1 )
      xNorm   = dnormi( n , x, 1 )
      penNorm = penParm(3)

!     Count basic nonlinear variables (used only for printing).

      nnb  = 0
      do j = 1, nnH
         if (hs(j) .eq. 3) nnb = nnb + 1
      end do

      if (inewB .gt. 0  .and.  iExit/10 .lt. 8) then
         k      = 1 + iExit/10
         call s4stat
     &      ( k, istate )
         call s4newB
     &      ( Wrap, iNewB, minimize, m, n, nb,
     &        nS, mBS, itn, nInf, sInf, fObj, iw(lkBS), hs,
     &        rw(lscales), bl, bu, x, rw(lxBS), istate,
     &        cw, lencw, iw, leniw )
      end if

      call snPRNT(13,
     &     ' Problem name                 '//mProb, iw, leniw )

      if (GotFuns) then
         if (Elastic) then
            write(str, 1900) itn, objTrue
            call snPRNT( 3, str, iw, leniw )

            write(str, 1920) majors, objLin
            call snPRNT( 3, str, iw, leniw )
            write(str, 1935) penNorm, fObj
            call snPRNT( 3, str, iw, leniw )

            if (nInfE .gt. 0) then ! objective is known
               write(str, 1915) wtInf, objMin
               call snPRNT( 3, str, iw, leniw )
               write(str, 2050) objInf
               call snPRNT( 3, str, iw, leniw )
            end if
         else ! Normal mode
            write(str, 1900) itn, objTrue
            call snPRNT( 3, str, iw, leniw )

            write(str, 1920) majors, objLin
            call snPRNT( 3, str, iw, leniw )

            if (nnCon .gt. 0) then ! Nonlinear constraints
               if (FPonly) then
                  write(str, 1930) penNorm, fObj
                  call snPRNT( 3, str, iw, leniw )
               else
                  write(str, 1935) penNorm, fObj
                  call snPRNT( 3, str, iw, leniw )
               end if
            else
               if (FPonly) then
!                 Relax. Nothing is printed
               else
                  write(str, 1940)          fObj
                  call snPRNT( 3, str, iw, leniw )
               end if
            end if
         end if

         if (Solver .eq. 'SNOPTA'  .or. Solver .eq. 'SNOPTC') then

!           SNOPTA and SNOPTC call one user-supplied function

            if (iw(lvlDer) .lt. 3  .or.
     &         (   lvlSrch .eq. 0 .and. Nonlinear)) then
               write(str, 1952) iw(nfObj1), iw(nfObj2)
               call snPRNT( 3, str, iw, leniw )
            else
               write(str, 1951) iw(nfObj1)
               call snPRNT( 3, str, iw, leniw )
            end if

            if (iw(lvlDer) .lt. 3) then
               write(str, 1953) iw(nfObj3), iw(nfObj4)
               call snPRNT( 3, str, iw, leniw )
            end if
         else

!           SNOPTB and NPOPT call separate obj and constraint functions.

            write(str, 1950) iw(nfObj1), iw(nfCon1)
            call snPRNT( 3, str, iw, leniw )

            if (iw(lvlDer) .lt. 3  .or.
     &         (   lvlSrch .eq. 0 .and. Nonlinear)) then
               write(str, 1955) iw(nfObj2), iw(nfCon2)
               call snPRNT( 3, str, iw, leniw )
            end if

            if (iw(lvlDer) .lt. 3) then
               if (iw(lvlDer) .eq. 0) then
                  write(str, 1960) iw(nfObj3), iw(nfCon3)
                  call snPRNT( 3, str, iw, leniw )

                  write(str, 1962) iw(nfObj4), iw(nfCon4)
                  call snPRNT( 3, str, iw, leniw )
               else if  (iw(lvlDer) .eq. 1) then
                  write(str, 1963) iw(nfCon3)
                  call snPRNT( 3, str, iw, leniw )

                  write(str, 1964) iw(nfCon4)
                  call snPRNT( 3, str, iw, leniw )
               else if  (iw(lvlDer) .eq. 2) then
                  write(str, 1965) iw(nfObj3)
                  call snPRNT( 3, str, iw, leniw )

                  write(str, 1966) iw(nfObj4)
                  call snPRNT( 3, str, iw, leniw )
               end if
            end if
         end if

         if (nS         .gt. 0) then
            write(str, 1970) nS, nnb
            call snPRNT( 3, str, iw, leniw )
         end if

         if (iw(cgItns) .gt. 0) then
            write(str, 1973) iw(cgItns)
            call snPRNT( 3, str, iw, leniw )
         end if

         write(str, 1975) nDegen, degen
         call snPRNT( 3, str, iw, leniw )

      else ! No functions computed
         if (nInf .gt. 0) then
            write(str, 1910) nInf, sInf
            call snPRNT( 3, str, iw, leniw )
         end if

         write(str, 1905) itn
         call snPRNT( 3, str, iw, leniw )

         if (Elastic) then
            if (nInfE .gt. 0) then
               write(str, 1980) nInfE, sInfE
               call snPRNT( 3, str, iw, leniw )
            end if
         end if
      end if

!     ------------------------------------------------------------------
!     Unscale, compute nonlinear constraint violations,
!     save basis files and prepare to print the solution.
!     Clock 3 is "Output time".
!     ------------------------------------------------------------------
      call s1time( 3, 0, iw, leniw, rw, lenrw )

!     Skip the functions if we don't have them.
!     Skip unscaling everything for infeasible linear constraints,
!     they have already been unscaled.

      lsSave  = iw(lvlScale)

      if (.not. GotFuns) then
         nnCon1 = 0
         nnObj1 = 0
         if (.not. FeasibleLC) then
            iw(lvlScale) = 0
         end if
      else
         nnCon1 = nnCon
         nnObj1 = nnObj
      end if

      call s4saveB
     &   ( iExit, SaveB, minimize, m, n, nb, nkx,
     &     nnCon0, nnCon1, nnH0, nnObj1, nNames, nS,
     &     itn, nInf, sInf, wtInf, maxVi, iObj, scaleObj, objTrue,
     &     pNorm1, pNorm2, piNorm, xNorm,
     &     neJ, nlocJ, locJ, indJ, Jcol, iw(lkx),
     &     iw(leState), hs, rw(lscales), bl, bu, rw(lFx), rw(lgObj),
     &     Names, pi, rc, x,
     &     cw, lencw, iw, leniw, rw, lenrw )

!     If task = 'Print', s4saveB prints the solution under the control
!     of lprSol (set by the  Solution  keyword in the SPECS file).
!     The printed solution may or may not be wanted, as follows:
!
!     lprSol = 0   means      No
!            = 2   means      Yes

      call s4saveB
     &   ( iExit, PrintS, minimize, m, n, nb, nkx,
     &     nnCon0, nnCon1, nnH0, nnObj, nNames, nS,
     &     itn, nInf, sInf, wtInf, maxVi, iObj, scaleObj, objTrue,
     &     pNorm1, pNorm2, piNorm, xNorm,
     &     neJ, nlocJ, locJ, indJ, Jcol, iw(lkx),
     &     iw(leState), hs, rw(lscales), bl, bu, rw(lFx), rw(lgObj),
     &     Names, pi, rc, x,
     &     cw, lencw, iw, leniw, rw, lenrw )
      iw(lvlScale) = lsSave

      call s1time(-3, 0, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     If the user hasn't already pulled the plug,
!     call the functions one last time with  Status .ge. 2.
!     Everything has been  unscaled, so we have to disable scaling.
!     modefg = 0  requests that no gradients are computed.
!     ------------------------------------------------------------------
      fObj = objLin

      if (GotFuns  .and.  iExit/10 .ne. 7  .and.  iExit/10 .ne. 6) then
         iw(npStatus) = 2 + min( iExit/10,4 )
         modefg       = 0

         lsSave       = iw(lvlScale)
         iw(lvlScale) = 0
         call funwrapper
     &      ( inform,
     &        modefg, NonlinearCon, NonlinearObj,
     &        n, negCon, nnCon0, nnCon,
     &        nnJac, nnH, nnObj0, nnObj,
     &        funcon, funobj, userHv,
     &        x, rw(lyCon),
     &        neJ, nlocJ, locJ, indJ,
     &        neH, nlocH, locH, indH, Hcol,
     &        rw(lfCon), fObj, rw(lgCon2), rw(lgObj2),
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         iw(lvlScale) = lsSave
      end if

      ! Save some things needed by solvers calling SNOPT

      rw(421) = objTrue      ! The objective (minimized or maximized)
      rw(422) = piNorm       ! Lagrange multiplier norm
      rw(423) = xNorm        ! Norm of the variables (for GAMS)
      rw(424) = wtInf        ! Infeasibility weight

      rw(433) = sInf + sInfE ! Sum of infeasibilities
      rw(434) = objLin       ! Linear    objective term

      iw(421) = itn          ! Total iteration count
      iw(423) = maxS         ! max # of superbasics

!     Nonlinear constraint information.

      rw(432) = maxVi        ! Inf norm of the constraint violation
      rw(435) = fObj         ! Nonlinear objective term
      rw(436) = penParm(3)   ! Norm of penalty parameters

      iw(422) = majors       ! Major iterations

      return

 1900 format(
     &     ' No. of iterations', i20, 2x,
     &     ' Objective', 6x, 1p, e22.10)
 1905 format(
     &     ' No. of iterations', i20 )
 1910 format(
     &     ' No. of infeasibilities', i15, 2x,
     &     ' Sum of infeas', 1p, e24.10)
 1915 format(
     &     ' Elastic weight', 12x, 1p, e11.1, 2x,
     &     ' Elastic   objective', 1p, e18.10)
 1920 format(
     &     ' No. of major iterations', i14, 2x,
     &     ' Linear    obj. term', 1p, e18.10)
 1930 format(
     &     ' Penalty parameter',   1p, e20.3, 2x,
     &     ' Norm (x - x0)**2   ', 1p, e18.10)
 1935 format(
     &     ' Penalty parameter', 1p, e20.3, 2x,
     &     ' Nonlinear obj. term', 1p, e18.10)
 1940 format(
     &     '                  ',           22x,
     &     ' Nonlinear obj. term', 1p, e18.10)
 1950 format(
     &     ' No. of calls to funobj', i15, 2x,
     &     ' No. of calls to funcon', i15)
 1951 format(
     &     ' User function calls (total)', i10)
 1952 format(
     &     ' User function calls (total)', i10, 2x,
     &     ' Calls with modes 1,2 (known g)', i7)
 1953 format(
     &     ' Calls for forward differencing', i7, 2x,
     &     ' Calls for central differencing', i7)
 1955 format(
     &     ' Calls with modes 1,2 (known g)', i7, 2x,
     &     ' Calls with modes 1,2 (known g)', i7)
 1960 format(
     &     ' Calls for forward differencing', i7, 2x,
     &     ' Calls for forward differencing', i7)
 1962 format(
     &     ' Calls for central differencing', i7, 2x,
     &     ' Calls for central differencing', i7)
 1963 format(
     &     '                  ',           22x,
     &     ' Calls for forward differencing', i7)
 1964 format(
     &     '                  ',           22x,
     &     ' Calls for central differencing', i7)
 1965 format(
     &     ' Calls for forward differencing', i7)
 1966 format(
     &     ' Calls for central differencing', i7)
 1970 format(
     &     ' No. of superbasics', i19, 2x,
     &     ' No. of basic nonlinears', i14)
 1973 format(
     &     ' No. of CG iterations', i17)
 1975 format(
     &     ' No. of degenerate steps', i14, 2x,
     &     ' Percentage', f27.2)
 1980 format(
     &     ' No. of infeas elastics', i15, 2x,
     &     ' Elastic infeas   ', 1p, e20.10)
 2050 format(
     &     '                  ',             22x,
     &     ' Elastic infeas   ',      1p, e20.10)

      end ! subroutine s8solve

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8defaults
     &   ( m, n, nnCon, nnJac, nnObj, iObj,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     iObj, nnCon, nnJac, nnObj, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     s8defaults checks the optional parameter values and possibly
!     changes them to reasonable values.
!
!     Note that checking occurs before the amount of working storage has
!     been defined.
!
!     See  snworkspace.info  for full documentation of cw, iw and rw.
!
!     15 Nov 1991: first version.
!     27 Apr 2001: wtMax introduced.
!     10 Dec 2002: Added defaults for LU Rook and Diagonal Pivoting.
!     31 Dec 2002: Added default MPS character names.
!     30 Jul 2003: Added default CG tolerance.
!     22 Jun 2004: Added default LU mod singularity tol
!     20 Dec 2004: Default LU tols reduced.
!     21 Dec 2004: Default LU tols fixed up.
!     09 Jun 2005: Default tolpiv same as in s5dflt.
!     02 Jul 2005: Default Utol's set back to eps1.
!     02 May 2006: lvlTim removed.
!     01 Sep 2007: stickyOp added.
!     25 Nov 2007: Hessian options added.
!     05 Apr 2014: Objective derivatives not verified in FP mode.
!     12 Sep 2014: FPonly set if nnObj .eq. 0 and iObj .eq. 0
!     25 Oct 2014: Bogus assignment of LUprnt removed (its set in s2BLU).
!     25 Oct 2014: nout set independently of lvlsys.
!     15 Nov 2014: No scaling is now the default.
!     21 Nov 2014: Default wtInf now 10^5.
!     18 Feb 2015: SNOPTA parameters set separately.
!     09 May 2015: Default bigFx now 10^14 (was 10^15).
!     10 May 2015: Zcndbd replaced by condZmax0.
!     12 May 2015: rtcondZbnd set here.
!     ==================================================================
      character
     &     blank*8, cdummy*8, mProb*8, mObj*8, mRhs*8, mRng*8, mBnd*8,
     &     Solver*8
      logical
     &     LinCon, Linear, NonlinearCon, Nonlinear
      integer
     &     cgItmx, DerOpt, eMode, iCrash, iBack, iDump, iLoadB,
     &     iInsrt, iNewB, iOldB, iPnch, iPrint, iReprt, iSoln, itnlim,
     &     jverf1, jverf2, jverf3, jverf4, jverf5, jverf6,
     &     kchk, kdegen, kFac, klog, kreset, ksav, kSumm,
     &     lprDbg, lprPrm, lprSch, lprScl, lprSol,
     &     lvlDer, lvlHess, lvlObjE, lvlPiv, lvlPre, lvlPPm, lvlSrch,
     &     lvlScale, lvlSys, lvlVer, m, maxmn, maxCol, maxR, maxS,
     &     mflush, minimize, minmax, minPrc, mjrPrint, mMajor,
     &     mMinor, mnrPrint, mQNmod, mskip, mNewSB, n, never, nnH,
     &     nout, nParPrU, nParPrLP, nParPrQP, nPr1, nPr2,
     &     ObjRow, QPsolver, stickyOp, TPivot
      double precision
     &     bigdx, bigFx, c4, c6, chzbnd, condUmax, condZmax0,
     &     Dens1, Dens2, eps, eps0, eps1, eps2, eps3, eps4, epsrf,
     &     etarg, fdint1, fdint2, Hcondbnd, Lmax1, Lmax2, maxTime,
     &     infBnd, proxWeight, rtcondZbnd, scltol, small, tCrash,
     &     tolCG, tolCon, tolDcp, tolDdp, tolDpp, tolDrp, tolDup,
     &     toldj3, tolFac, tolOptFP, tolOptNP, tolpiv, tolOptQP,
     &     tolRow, tolSwp, tolUpd, tolx, Uspace, Utol1, Utol1m,
     &     Utol2, Utol2m, viLim, wolfeG, wtInf0, wtMax, xdlim, xPen0
!     ------------------------------------------------------------------
      integer            QPChol,     CG,     QN
      parameter         (QPChol = 0, CG = 1, QN = 2)
      integer            LM        , FM
      parameter         (LM     = 0, FM = 1)
      parameter         (cdummy ='-1111111', blank ='        ')
      integer            idummy
      parameter         (idummy = -11111)
      double precision   zero,             one
      parameter         (zero   =  0.0d+0, one    =   1.0d+0)
      double precision   ten
      parameter         (ten    = 10.0d+0)
      double precision   tenp6,            hundrd
      parameter         (tenp6  = 1.0d+6,  hundrd = 100.0d+0)
!     ------------------------------------------------------------------
!     Set some local machine-dependent constants.

      eps        = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      eps0       = rw(  2) ! eps**(4/5)          IEEE DP  3.00e-13
      eps1       = rw(  3) ! eps**(2/3)          IEEE DP  3.67e-11
      eps2       = rw(  4) ! eps**(1/2)          IEEE DP  1.49e-08
      eps3       = rw(  5) ! eps**(1/3)          IEEE DP  6.05e-06
      eps4       = rw(  6) ! eps**(1/4)          IEEE DP  1.22e-04

      Solver     = cw(1)
!     ------------------------------------------------------------------
!     rw(51)--rw(150): optional parameters set via the specs file.
!     ------------------------------------------------------------------
      tolOptFP   = rw( 51) ! Minor Phase 1 Opt tol
      tolOptQP   = rw( 52) ! Minor Phase 2 Opt tol
      tolOptNP   = rw( 53) ! Major Optimality tolerance
      tolCG      = rw( 54) ! cg tolerance

      tolx       = rw( 56) ! Minor feasibility tolerance.
      tolCon     = rw( 57) ! Major feasibility tolerance.

      tolpiv     = rw( 60) ! excludes small elements of y
      tolrow     = rw( 61) ! tolerance for the row error
      tCrash     = rw( 62) ! crash tolerance
      Utol1m     = rw( 63) ! abs tol for small diag of U in LU mod
      Utol2m     = rw( 64) ! rel tol for small diag of U in LU mod
      tolswp     = rw( 65) ! LU swap tolerance
      tolFac     = rw( 66) ! User-defined LU factor tolerance
      tolUpd     = rw( 67) ! User-defined LU update tolerance
      infBnd     = rw( 70) ! definition of an infinite bound
      bigFx      = rw( 71) ! unbounded objective
      bigdx      = rw( 72) ! unbounded step
      epsrf      = rw( 73) ! relative function precision.
      fdint1     = rw( 76) ! (1) forwrd diff. interval
      fdint2     = rw( 77) ! (2) cntrl  diff. interval
      maxTime    = rw( 79) ! Time limit
      xdlim      = rw( 80) ! Step limit
      vilim      = rw( 81) ! violation limit
      etarg      = rw( 83) ! Quasi-Newton QP rg tolerance
      wolfeG     = rw( 84) ! line search tolerance.
      Hcondbnd   = rw( 85) ! bound on the condition of Hz
      condZmax0  = rw( 86) ! bound on the condition of Z
      condUmax   = rw( 87) ! max cond estimator for U with H = U'U
      wtInf0     = rw( 88) ! infeasibility weight
      xPen0      = rw( 89) ! initial penalty parameter.
      wtMax      = rw( 90) ! max     infeasibility weight
      proxWeight = rw( 91) ! Proximal-point weight

      scltol     = rw( 92) ! scale tolerance.
!     ------------------------------------------------------------------
!     rw(151)--rw(180) are parmLU parameters for LUSOL (some optional).
!     ------------------------------------------------------------------
      small      = rw(153) ! defn of small real.
      Utol1      = rw(154) ! abs tol for small diag of U.
      Utol2      = rw(155) ! rel tol for small diag of U.
      Uspace     = rw(156) ! limit on waste space in U.
      Dens1      = rw(157) ! switch to search maxcol columns and no rows.
      Dens2      = rw(158) ! switch to dense LU.
!     ------------------------------------------------------------------
!     rw(181)--rw(199) pass parameters into various routines.
!     ------------------------------------------------------------------
!     toldj3     = rw(186) ! current optimality tol
!     ------------------------------------------------------------------
!     iw(1)--iw(50): I/O file numbers and dimensions.
!     ------------------------------------------------------------------
      iPrint     = iw( 12) ! Print file
!     ------------------------------------------------------------------
!     iw(51)--iw(150): optional parameters set via the specs file.
!     ------------------------------------------------------------------
      maxR       = iw( 52) ! max columns of R.
      maxS       = iw( 53) ! max # of superbasics
      mQNmod     = iw( 54) ! (ge 0) max # of BFGS updates
      QPsolver   = iw( 55) ! 0(1) => QP(QN) QP solver
      eMode      = iw( 56) ! >0    => use elastic mode
      kchk       = iw( 58) ! check (row) frequency
      kFac       = iw( 59) ! factorization frequency
      ksav       = iw( 60) ! save basis map
      klog       = iw( 61) ! log/print frequency
      kSumm      = iw( 62) ! Summary print frequency
      kDegen     = iw( 63) ! max. expansions of featol
      kReset     = iw( 64) ! Hessian frequency
      mFlush     = iw( 66) ! Hessian flush
      mSkip      = iw( 67) ! # largest value of nSkip
!     lvlStart   = iw( 69) ! = 0:1:2:3 => cold:basis:warm:hot start
      lvlDer     = iw( 70) ! = 0, 1 or 2, the derivative level
      lvlSys     = iw( 71) ! > 0   => print system info
      lvlHess    = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian
      lvlObjE    = iw( 73) ! Elastic option
      lvlScale   = iw( 75) ! scale option
      lvlSrch    = iw( 76) ! >0     => use derivatives in the line search
      lvlPre     = iw( 77) ! >0    => QN preconditioned CG
      lvlVer     = iw( 78) ! Verify level
      lvlPPm     = iw( 79) ! 1(2)-norm proximal point method for x0
      lvlPiv     = iw( 80) ! 0(1 2 3) LU partial(rook complete diagonal) pivoting
      lprPrm     = iw( 81) ! > 0    => parms are printed
      lprSch     = iw( 82) ! line search debug starting itn
      lprScl     = iw( 83) ! > 0    => print the scales
      lprSol     = iw( 84) ! > 0    => print the solution
      lprDbg     = iw( 85) ! > 0    => private debug print
      minmax     = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
      iCrash     = iw( 88) ! Crash option
      itnlim     = iw( 89) ! limit on total iterations
      mMajor     = iw( 90) ! limit on major iterations
      mMinor     = iw( 91) ! limit on minor iterations
      mjrPrint   = iw( 92) ! Major print level
      mnrPrint   = iw( 93) ! Minor print level
      mNewSB     = iw( 95) ! maximum # of new superbasics per major
      cgItmx     = iw( 97) ! CG iteration limit
      nParPrU    = iw(101) ! # of partial pricing sections
      ObjRow     = iw(103) ! Objective row of user-defined F
      DerOpt     = iw(104) ! 0, 1, 2 => derivative option
      jverf1     = iw(110) ! col # to start derivative checking
      jverf2     = iw(111) ! col # to stop  derivative checking
      jverf3     = iw(112) ! col # to start derivative checking
      jverf4     = iw(113) ! col # to stop  derivative checking
      jverf5     = iw(114) ! start col for Hessian checking
      jverf6     = iw(115) ! stop  col for Hessian checking
      stickyOp   = iw(116) ! > 0 => optional parameters are sticky
      iBack      = iw(120) ! backup file
      iDump      = iw(121) ! dump file
      iLoadB     = iw(122) ! load file
      iNewB      = iw(124) ! new basis file
      iInsrt     = iw(125) ! insert file
      iOldB      = iw(126) ! old basis file
      iPnch      = iw(127) ! punch file
      iReprt     = iw(130) ! report file
      iSoln      = iw(131) ! solution file
!     ------------------------------------------------------------------
!     iw(151)--iw(180) are luparm parameters for LUSOL (some optional).
!     ------------------------------------------------------------------
      nout       = iw(151) ! unit # for printed messages
      maxcol     = iw(153) ! lu1fac: max. # columns
!     ------------------------------------------------------------------
!     Character  workspace.
!     cw(51)--cw(150): optional parameters
!     ------------------------------------------------------------------
      mProb      = cw( 51) ! Problem name
      mObj       = cw( 52) ! Objective name
      mRhs       = cw( 53) ! rhs name
      mRng       = cw( 54) ! range name
      mBnd       = cw( 55) ! bounds name
!     ------------------------------------------------------------------
      c4         = max( 1.0d-4, eps3 )
      c6         = max( 1.0d-6, eps2 )
      cHzbnd     = max ( one/(hundrd*eps), tenp6 )

      never      = 99999999

!     ===============================================================
!     Check the optional parameters.
!     ===============================================================
      if (nnCon .eq. 0) nnJac = 0
      if (nnJac .eq. 0) nnCon = 0
      nnH = max( nnJac, nnObj )

      LinCon       = nnCon   .eq. 0
      NonlinearCon = nnCon   .gt. 0
      Linear       = nnH     .eq. 0
      Nonlinear    = nnH     .gt. 0

      if (iBack  .eq. idummy ) iBack    =     0
      if (iDump  .eq. idummy ) iDump    =     0
      if (iLoadB .eq. idummy ) iLoadB   =     0
      if (iNewB  .eq. idummy ) iNewB    =     0
      if (iInsrt .eq. idummy ) iInsrt   =     0
      if (iOldB  .eq. idummy ) iOldB    =     0
      if (iPnch  .eq. idummy ) iPnch    =     0
      if (iReprt .eq. idummy ) iReprt   =     0
      if (iSoln  .eq. idummy ) iSoln    =     0

!     Set unspecified frequencies or silly values to defaults.

      if (kchk   .eq. idummy ) kchk     =    60
      if (kFac   .le.    0   ) then
                               kFac     =   100
              if (Nonlinear  ) kFac     =    50
      end if
      if (klog  .eq. idummy  ) klog     =   100
      if (kSumm .eq. idummy  ) kSumm    =   100
      if (ksav  .eq. idummy  ) ksav     =   100
      if (kDegen.eq. idummy  ) kDegen   = 10000
      if (mFlush.eq. idummy  ) mFlush   =     0

!     Sometimes, frequency 0 means "almost never".

      if (kchk   .le. 0      ) kchk     = never
      if (mFlush .le. 0      ) mFlush   = never
      if (klog   .le. 0      ) klog     = never
      if (ksav   .le. 0      ) ksav     = never
      if (kSumm  .le. 0      ) kSumm    = never
      if (kDegen .le. 0      ) kDegen   = never
      if (kReset .le. 0      ) kReset   = never

      if (iCrash .lt. 0      ) iCrash   =  3

      if (minmax .eq. idummy ) minmax   =  1
      if (nnObj  .eq.  0 .and. iObj .eq. 0 .or. ObjRow .eq. 0
     &                       ) minmax   =  0

      if (minmax .eq. -1) then
                               minimize = -1
                         else
                               minimize =  1
                         end if
      if (ObjRow .gt.  0 .and. minmax .eq. 0
     &                       ) ObjRow   =  0

      if (mjrPrint.eq. idummy) mjrPrint =  1
      if (mnrPrint.eq. idummy) mnrPrint =  1

!     if (mMinor   .lt. 0    ) mMinor   = max( 1000,5*max( n,m ) )
      if (mMinor   .lt. 0    ) mMinor   = 500
      if (mMajor   .lt. 0    ) mMajor   = max( 1000,3*max( n,m ) )
      if (mSkip    .lt. 0  .and.  LinCon
     &                       ) mSkip    = never
      if (mSkip    .lt. 0  .and.  NonlinearCon
     &                       ) mSkip    =  2
      if (mNewSB   .le. 0    ) mNewSB   = 99

      if (lprDbg   .lt. 0    ) lprDbg   =  0
      if (lprPrm   .lt. 0    ) lprPrm   =  1
      if (lprSch   .lt. 0    ) lprSch   = never
      if (lprScl   .lt. 0    ) lprScl   =  0
      if (lprSol   .lt. 0    ) lprSol   =  2

!     lvlStart is checked in s3argA or s3argB
!     if (lvlStart.lt. 0     ) lvlStart =  0

      if (Solver(1:6) .eq. 'SNOPTA') then
         if      (DerOpt .ne. idummy ) then
            if (DerOpt .lt. 0  .or.  DerOpt .gt. 1)
     &                            DerOpt   =  1
            if (DerOpt .eq. 0   ) lvlDer   =  0
            if (DerOpt .eq. 1   ) lvlDer   =  3
         else
                                  DerOpt   =  1
                                  lvlDer   =  3
         end if

         if (lvlSrch.lt. 0      ) then
               if (DerOpt .eq. 0) lvlSrch  =  0
               if (DerOpt .eq. 1) lvlSrch  =  1
         end if
      else
         if (lvlDer .ne. idummy ) then
            if (lvlDer .lt. 0  .or.  lvlDer .gt. 3)
     &                            lvlDer   =  3
         else
                                  lvlDer   =  3
         end if
         if (lvlSrch.lt. 0      ) then
               if (lvlDer .ne. 3) lvlSrch  =  0
               if (lvlSrch.lt. 0) lvlSrch  =  1
         end if
      end if

      if (lvlVer .eq. idummy    ) lvlVer   =  0
      if (lvlVer .lt. 0         ) lvlVer   = -1
      if (lvlVer .gt. 3         ) lvlVer   =  0
                                  lvlObjE  =  2

      if (minmax .eq. 0) then
         if (lvlDer .eq. 2      ) lvlDer   =  3
         if (lvlDer .eq. 0      ) lvlDer   =  1
         if (lvlVer .eq. 1      ) lvlVer   =  0
         if (lvlVer .eq. 3      ) lvlVer   =  2
      end if

!     Check  START and STOP  column numbers for derivative checking.

      if (Solver(1:6) .eq. 'SNOPTA') then
                                  jverf1 = 1
                                  jverf2 = n
         if (lvlVer .eq. 2  .or.
     &       lvlVer .eq. 0      ) jverf2 = 0
                                  jverf3 = 1
                                  jverf4 = n
         if (lvlVer .eq. 1  .or.
     &       lvlVer .eq. 0      ) jverf4 = 0
                                  jverf5 = 1
                                  jverf6 = n
         if (lvlVer .eq. 1  .or.
     &       lvlVer .eq. 0      ) jverf6 = 0
      else
         if (jverf1 .le. 0      ) jverf1 = 1
         if (jverf2 .lt. 0      ) jverf2 = n
         if (lvlVer .eq. 2  .or.
     &       lvlVer .eq. 0      ) jverf2 = 0

         if (jverf3 .le. 0      ) jverf3 = 1
         if (jverf4 .lt. 0      ) jverf4 = n
         if (lvlVer .eq. 1  .or.
     &       lvlVer .eq. 0      ) jverf4 = 0

         if (jverf5 .le. 0      ) jverf5 = 1
         if (jverf6 .lt. 0      ) jverf6 = n
         if (lvlVer .eq. 1  .or.
     &       lvlVer .eq. 0      ) jverf6 = 0
      end if

      if (lvlSys .lt. 0      ) lvlSys   =  0
      if (lvlPPm .lt. 0      ) lvlPPm   =  1
                               eMode    =  1

      if (stickyOp .lt. 0    ) stickyOp =  0

!     Check superbasics limit maxS and Hessian dimension maxR.

      if (Nonlinear) then
         if (maxR .lt. 0     ) maxR      = min( 2000, nnH+1 )
         if (maxS .lt. 0     ) maxS      =            nnH+1
                               maxR      = max( min( maxR ,n ) , 0 )
                               maxS      = max( min( maxS ,n ) , 1 )
      else ! Linear
         if (maxS   .le. 0   ) maxS      = 1
         if (maxR   .le. 0   ) maxR      = 1
      end if

      if (maxS     .lt. maxR   ) maxR     = maxS

      if (QPsolver .lt. 0      ) QPsolver = QPChol
      if (maxR     .eq. 0      ) QPsolver = CG

      if (QPsolver .eq. QN     ) lvlPre   = 0
      if (QPsolver .eq. QPChol ) lvlPre   = 0
      if (lvlPre   .gt. 1      ) lvlPre   = 1
      if (lvlPre   .lt. 0  .and.  QPsolver .eq. CG)
     &                           lvlPre   = 0

      if (cgItmx   .lt. 0      ) cgItmx   = 100

      if (QPsolver .eq. CG  .or.  maxR   .lt. maxS) then
           if (lvlHess.lt. 0   ) lvlHess  = LM
           if (mQNmod .lt. 0   ) mQNmod   = 10
      else
         if (lvlHess .lt. 0 .and.  nnH  .gt. 75  )
     &                           lvlHess  = LM
         if (lvlHess .lt. 0 .and.  nnH  .le. 75  )
     &                           lvlHess  = FM
         if (lvlHess .eq. FM   ) mQNmod   = kReset
         if (mQNmod  .lt. 0    ) mQNmod   = 10
      end if

!     ---------------------------------
!     CG QP optional parameters
!     ---------------------------------
      if (etarg    .lt. zero  .or.
     &    etarg    .gt. one    ) etarg   = 0.1d+0

!     Check other options.

!     if (lvlScale .lt. 0  .and.
!    &    nnCon    .eq. 0      ) lvlScale   =  2
!     if (lvlScale .lt. 0      ) lvlScale   =  1
!                                lvlScale   =  min( lvlScale, 2 )
!     if (lvlScale .eq. 1  .and.  nnJac .ge. n)
!    &                           lvlScale   = 0

      if (lvlScale .lt. 0      ) lvlScale   =  0
                                 lvlScale   =  min( lvlScale, 2 )
      if (lvlScale .eq. 1  .and.  nnJac .ge. n)
     &                           lvlScale   = 0

!     Partial pricing parameters
                               minPrc   = 10
                               maxmn    = max( m, n )

!     Temporarily go back to what we had before.
!     Set nParPrU to previous value.
!     Then nParPrLP/nParPrQP will take that value.

      if (nParPrU   .le. 0     ) then
                                 nParPrU = 10
         if (Nonlinear)          nParPrU =  1
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

      if (proxWeight.lt. zero ) proxWeight = 1.0d-4
      if (maxTime   .lt. zero ) maxTime    = zero
      if (infBnd    .lt. zero ) infBnd     = 1.0d+20
      if (bigFx     .le. zero ) bigFx      = 1.0d+10
      if (bigdx     .le. zero ) bigdx      = infBnd
      if (Hcondbnd  .le. zero ) Hcondbnd   = cHzbnd
                                rtcondZbnd = 1.0d+15
      if (condUmax  .le. zero ) condUmax   = sqrt(cHzbnd)
      if (tCrash    .lt. zero  .or.
     &    tCrash    .ge. one  ) tCrash     = 0.1d+0
      if (vilim     .le. zero ) vilim      = 1.0d+1
      if (wolfeG    .lt. zero  .or.
     &    wolfeG    .gt. one  ) wolfeG     = 0.9d+0
      if (wtMax     .lt. zero ) wtMax      = 1.0d+10
      if (xdlim     .le. zero ) xdlim      = 2.0d+0
      if (xPen0     .lt. zero ) xPen0      = zero
      if (condZmax0 .le. zero ) condZmax0  = 1.0d+6

!     ----------------------------------------
!     Set up the parameters for lu1fac.
!     LUprnt > 0 gives LU output on unit nout.
!     LUprnt is set in s2BLU.
!     ----------------------------------------
      if (maxcol .lt.  0     ) maxcol =  5
      if (nout   .eq.  idummy) nout   = iPrint

      if (lvlPiv .le.  0     ) lvlPiv =  0
      if (lvlPiv .gt.  3     ) lvlPiv =  0
                               TPivot =  lvlPiv
      if (Linear) then
                               tolDpp =  hundrd
                               tolDrp =  ten
                               tolDcp =  ten
                               tolDdp =  ten
                               tolDup =  ten
      else ! nonlinear
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

!     Set some SQP tolerances.
!     Set the minor and major optimality tolerances.
!     Solve the QP subproblems fairly accurately even if the
!     NLP Optimality Tolerance is big.

      if (tolOptNP .le. zero) then
                               tolOptNP = 2.0d+0*c6
         if (epsrf .gt. zero ) tolOptNP = max(tolOptNP, sqrt(ten*epsrf))
      end if
      if (tolOptQP .le. zero ) tolOptQP = min(c6      , tolOptNP/2.0d+0)
      if (tolOptFP .lt. zero ) tolOptFP =  c6
      if (tolCG    .le. zero ) tolCG    =  1.0d-2
      if (tolrow   .le. zero ) tolrow   =  c4
      if (tolswp   .le. zero ) tolswp   =  eps4
      if (tolx     .le. zero ) tolx     =  c6
      if (tolCon   .le. eps  ) tolCon   =  c6
                               toldj3   =  tolOptQP
      if (scltol   .le. zero ) scltol   =  0.90d+0
      if (scltol   .ge. one  ) scltol   =  0.99d+0
      if (tolpiv   .le. zero ) tolpiv   =  eps1

      if (LinCon) then
         if (wtInf0.lt. zero ) wtInf0 =  1.0d+0
      else
         if (wtInf0.lt. zero ) wtInf0 =  1.0d+5
      end if

      if (epsrf    .le. zero ) epsrf  = eps0
      if (fdint1.le. zero    ) fdint1 = sqrt(epsrf)
      if (fdint2.le. zero    ) fdint2 = epsrf**0.33333d+0

      if (iBack  .eq. iNewB  ) iBack  = 0
      if (itnlim .lt. 0      ) itnlim = max(10000, 10*max(n,m))

!     Set default names (they may be printed by the basis routines).

      if (mProb  .eq. cdummy ) mProb  = blank
      if (mObj   .eq. cdummy ) mObj   = blank
      if (mRhs   .eq. cdummy ) mRhs   = blank
      if (mRng   .eq. cdummy ) mRng   = blank
      if (mBnd   .eq. cdummy ) mBnd   = blank

!     ------------------------------------------------------------------
!     Done.
!     Re-assign the options to their respective work arrays.
!     ------------------------------------------------------------------
      rw( 51)  =  tolOptFP
      rw( 52)  =  tolOptQP
      rw( 53)  =  tolOptNP
      rw( 54)  =  tolCG
      rw( 56)  =  tolx
      rw( 57)  =  tolCon
      rw( 60)  =  tolpiv
      rw( 61)  =  tolrow
      rw( 62)  =  tCrash
      rw( 63)  =  Utol1m
      rw( 64)  =  Utol2m
      rw( 65)  =  tolswp
      rw( 66)  =  tolFac
      rw( 67)  =  tolUpd
      rw( 70)  =  infBnd
      rw( 71)  =  bigFx
      rw( 72)  =  bigdx
      rw( 73)  =  epsrf
      rw( 76)  =  fdint1
      rw( 77)  =  fdint2
      rw( 79)  =  maxTime
      rw( 80)  =  xdlim
      rw( 81)  =  vilim
      rw( 83)  =  etarg
      rw( 84)  =  wolfeG
      rw( 85)  =  Hcondbnd
      rw( 86)  =  condZmax0
      rw( 87)  =  condUmax
      rw( 88)  =  wtInf0
      rw( 89)  =  xPen0
      rw( 90)  =  wtMax
      rw( 91)  =  proxWeight
      rw( 92)  =  scltol
      rw(151)  =  Lmax1
      rw(152)  =  Lmax2
      rw(153)  =  small
      rw(154)  =  Utol1
      rw(155)  =  Utol2
      rw(156)  =  Uspace
      rw(157)  =  Dens1
      rw(158)  =  Dens2
!     Dependent parameters set in s8defaults.
      rw(181)  =  tolDpp
      rw(182)  =  tolDcp
      rw(183)  =  tolDup
      rw(186)  =  toldj3
      rw(187)  =  tolDrp

!     Addresses for integer quantities.

      iw( 52)  =  maxR
      iw( 53)  =  maxS
      iw( 54)  =  mQNmod
      iw( 55)  =  QPsolver
      iw( 56)  =  eMode
      iw( 58)  =  kchk
      iw( 59)  =  kFac
      iw( 60)  =  ksav
      iw( 61)  =  klog
      iw( 62)  =  kSumm
      iw( 63)  =  kDegen
      iw( 64)  =  kReset
      iw( 66)  =  mFlush
      iw( 67)  =  mSkip
!     iw( 69)  =  lvlStart
      iw( 70)  =  lvlDer
      iw( 71)  =  lvlSys
      iw( 72)  =  lvlHess
      iw( 73)  =  lvlObjE
      iw( 75)  =  lvlScale
      iw( 76)  =  lvlSrch
      iw( 77)  =  lvlPre
      iw( 78)  =  lvlVer
      iw( 79)  =  lvlPPm
      iw( 80)  =  lvlPiv
      iw( 81)  =  lprPrm
      iw( 82)  =  lprSch
      iw( 83)  =  lprScl
      iw( 84)  =  lprSol
      iw( 85)  =  lprDbg
      iw( 87)  =  minmax
      iw( 88)  =  iCrash
      iw( 89)  =  itnlim
      iw( 90)  =  mMajor
      iw( 91)  =  mMinor
      iw( 92)  =  mjrPrint
      iw( 93)  =  mnrPrint
      iw( 95)  =  mNewSB
      iw( 97)  =  cgItmx
      iw(103)  =  ObjRow
      iw(104)  =  DerOpt
      iw(110)  =  jverf1
      iw(111)  =  jverf2
      iw(112)  =  jverf3
      iw(113)  =  jverf4
      iw(114)  =  jverf5
      iw(115)  =  jverf6
      iw(116)  =  stickyOp
      iw(120)  =  iBack
      iw(121)  =  iDump
      iw(122)  =  iLoadB
      iw(124)  =  iNewB
      iw(125)  =  iInsrt
      iw(126)  =  iOldB
      iw(127)  =  iPnch
      iw(130)  =  iReprt
      iw(131)  =  iSoln
      iw(151)  =  nout
      iw(153)  =  maxcol
      iw(156)  =  TPivot
!     Dependent parameters set in s8defaults.
      iw( 99)  =  nParPrLP
      iw(100)  =  nParPrQP
      rw(192)  =  rtcondZbnd
      iw(199)  =  minimize

!     Addresses for character quantities.

      cw( 51)  =  mProb
      cw( 52)  =  mObj
      cw( 53)  =  mRhs
      cw( 54)  =  mRng
      cw( 55)  =  mBnd

      end ! subroutine s8defaults

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Map
     &   ( m, n, negCon, nkx, nnCon, nnJac, nnObjU, nnObj, nnH,
     &     lenR, maxR, maxS, mQNmod, lvlHess,
     &     nextcw, nextiw, nextrw, iw, leniw )

      implicit
     &     none
      integer
     &     leniw, lenR, lvlHess, m, maxR, maxS, mQNmod, n, negCon,
     &     nextcw, nextiw, nextrw, nkx, nnCon, nnJac, nnObjU, nnObj,
     &     nnH, iw(leniw)

!     ==================================================================
!     s8Map   allocates all array storage for snopt,
!     using the values:
!        m    , n    , neJ
!        maxS                          Set in s8defaults.
!        nnObj, nnObjU, nnCon, nnJac   Set in specs file or arguments.
!        lenR , negCon                 Set in calling program
!
!     On exit,
!        nextcw, nextiw, nextrw are pointers to the next elements of
!                               free space in cw, iw, and rw.
!
!     29 Dec 2000: First version of s8Map.
!     24 Jan 2003: Added workspace for SYMMLQ
!     18 Jun 2008: Added space for iy2, pBS and rg2.
!     11 Sep 2014: Added nnObjU and nnObj for FP mode.
!     ==================================================================
      integer
     &     lscales, lblQP, lbuQP, lblBS, lbuBS, lblSave, lbuSave,
     &     ldyCon, ldg, ldx, lenfR, leType, leState, lfCon, lfCon1,
     &     lfCon2, lfR, lFv, lFx, lgCon, lgCon1, lgCon2, lgConU,
     &     lgObj, lgObj1, lgObj2, lgObjU, lgSave,
     &     lgBS, lgQP, lHd, lHdx, lfeasType, liy, liy1, liy2,
     &     lkBS, lkx, lyCon, lyCon1, lyCon2, llocG, lpBS, lQPrhs, lR,
     &     lr1, lr2, lrg, lrg2, lS, ls1, ls2, ls3, lU0, lUx, lV, lx0,
     &     lx1, lxBS, lxScaled, lxPen, lxQP, lxQP0, ly, ly1, ly2, ly3,
     &     mBS, nb, ngQP, nlocG
!     ------------------------------------------------------------------
      integer            LM        , FM
      parameter         (LM     = 0, FM     = 1)
!     ------------------------------------------------------------------

!     All dimensions are computed from
!        m     , n    , ne
!        lenR  , maxS , mQMmod
!        nnObjU, nnCon, nnJac
!        negCon

      ngQP     = nnH
      mBS      = m        + maxS
      nb       = n        + m

!     Nonlinear constraints.

      nlocG    = nnJac    + 1

!     Addresses for the integer arrays.

      lkx       = nextiw
      lfeasType = lkx       + nkx
      lkBS      = lfeasType + mBS
      leState   = lkBS      + mBS
      leType    = leState   + nb
      liy       = leType    + nb
      liy1      = liy       + nb
      liy2      = liy1      + nb
      nextiw    = liy2      + nb

!     Addresses for the double precision arrays.

      lscales   = nextrw
      ly        = lscales   + nb
      ly1       = ly        + nb
      ly2       = ly1       + nb
      ly3       = ly2       + nb      ! SYMMLQ workspace
      ls1       = ly3       + nb      ! SYMMLQ workspace
      ls2       = ls1       + maxS    ! SYMMLQ workspace
      ls3       = ls2       + maxS    ! SYMMLQ workspace
      lr1       = ls3       + maxS    ! SYMMLQ workspace
      lr2       = lr1       + maxS    ! SYMMLQ workspace
      lblQP     = lr2       + maxS
      lbuQP     = lblQP     + nb
      lblBS     = lbuQP     + nb
      lbuBS     = lblBS     + mBS
      lxBS      = lbuBS     + mBS
      lxScaled  = lxBS      + mBS
      lgBS      = lxScaled  + nnH
      lgQP      = lgBS      + mBS
      lUx       = lgQP      + ngQP
      lHdx      = lUx       + nnH
      lpBS      = lHdx      + nnH
      lU0       = lpBS      + nb
      ldg       = lU0       + nnH
      lR        = ldg       + nnH
      lrg       = lR        + lenR
      lrg2      = lrg       + maxS
      lblSave   = lrg2      + maxS
      lbuSave   = lblSave   + nb
      nextrw    = lbuSave   + nb

!     Nonlinear Objective.

      lgObj     = nextrw
      lgObj1    = lgObj     + nnObj
      lgObj2    = lgObj1    + nnObj
      lgObjU    = lgObj2    + nnObj
      lgSave    = lgObjU    + nnObj
      nextrw    = lgSave    + nnObjU

!     Nonlinear constraints.

      llocG     = nextiw
      nextiw    = llocG     + nlocG

      lfCon     = nextrw
      lfCon1    = lfCon     + nnCon
      lfCon2    = lfCon1    + nnCon
      lFx       = lfCon2    + nnCon
      lFv       = lFx       + nnCon
      lyCon     = lFv       + nnCon
      lyCon1    = lyCon     + nnCon
      lyCon2    = lyCon1    + nnCon
      ldyCon    = lyCon2    + nnCon
      lxPen     = ldyCon    + nnCon
      lgCon     = lxPen     + nnCon
      lgCon1    = lgCon     + negCon
      lgCon2    = lgCon1    + negCon
      lgConU    = lgCon2    + negCon
      lQPrhs    = lgConU    + negCon
      ldx       = lQPrhs    + m
      lxQP      = ldx       + nb
      lxQP0     = lxQP      + nb
      lx0       = lxQP0     + nb
      lx1       = lx0       + nb
      nextrw    = lx1       + nb

!     Store the addresses in iw.

      iw(251) = lkx
      iw(260) = llocG

      iw(271) = lblQP
      iw(272) = lbuQP
      iw(273) = lblBS
      iw(274) = lbuBS
      iw(275) = lblSave
      iw(276) = lbuSave

      iw(277) = lpBS
      iw(278) = lQPrhs

      iw(283) = leType
      iw(284) = lfeasType
      iw(285) = leState

      iw(287) = ldx
      iw(288) = lHdx
      iw(289) = ldg
      iw(290) = lgQP
      iw(291) = lgBS
      iw(292) = lkBS
      iw(293) = lrg
      iw(294) = lrg2
      iw(295) = lR
      iw(296) = lscales
      iw(297) = lgObj
      iw(298) = lx0

      iw(300) = lx1
      iw(301) = lxBS
      iw(302) = lxScaled

      iw(304) = lxPen
      iw(305) = lxQP
      iw(306) = lxQP0

      iw(308) = liy
      iw(309) = liy1
      iw(310) = liy2
      iw(311) = ly
      iw(312) = ly1
      iw(313) = ly2
      iw(314) = ly3

      iw(316) = lfCon
      iw(317) = lfCon1
      iw(318) = lfCon2
      iw(319) = lgConU
      iw(320) = lgCon
      iw(321) = lgCon1
      iw(322) = lgCon2
      iw(323) = lgObjU
      iw(324) = lgObj1
      iw(325) = lgObj2

      iw(336) = lFx
      iw(337) = lFv
      iw(339) = lgSave
      iw(345) = lUx
      iw(346) = lU0

      iw(348) = lyCon
      iw(349) = lyCon1
      iw(350) = lyCon2
      iw(351) = ldyCon

      iw(353) = lr1
      iw(354) = lr2
      iw(355) = ls1
      iw(356) = ls2
      iw(357) = ls3

!     Allocate space for an approximate Hessian.
!     The amount will depend on the method selected.

      if (lvlHess .eq. LM) then
!        ---------------------------------------------------------------
!        Compute the addresses of the limited-memory arrays.
!        These are saved and used for subsequent entries.
!        ---------------------------------------------------------------
         lHd     = nextrw
         lS      = lHd    + nnH
         lV      = lS     + nnH*mQNmod
         nextrw  = lV     + nnH*mQNmod

         iw(347) = lHd
         iw(401) = lS
         iw(402) = lV

      else if (lvlHess .eq. FM) then
         lenfR   = nnH*(nnH + 1)/2

         lHd     = nextrw
         lfR     = lHd    + nnH
         nextrw  = lfR    + lenfR

         iw(347) = lHd
         iw(391) = lfR
         iw(392) = lenfR
      end if

      end ! subroutine s8Map

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8firstCall
     &   ( iExit,
     &     FeasibleLC, GotFuns, NonlinearCon, NonlinearObj,
     &     n, nb, nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &     funwrapper, funcon, funobj, userHv,
     &     bl, bu, x, x1, yCon,
     &     neJ   , nlocJ, locJ, indJ,
     &     neH   , nlocH, locH, indH, Hcol,
     &     negCon, nlocG, locG, fObj,
     &     fCon , gCon , gObj ,
     &     fCon2, gCon2, gObj2, y, y1,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     funwrapper, funcon, funobj, userHv
      logical
     &     FeasibleLC, GotFuns, NonlinearCon, NonlinearObj
      integer
     &     iExit, lencu, leniu, lenru, lencw, leniw, lenrw,
     &     n, nb, negCon, neH, neJ, nlocG, nlocH, nlocJ,
     &     nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &     indH(neH), indJ(neJ), locG(nlocG), locH(nlocH), locJ(nlocJ),
     &     iu(leniu), iw(leniw)
      double precision
     &     fObj, bl(nb), bu(nb), fCon(nnCon0),  fCon2(nnCon0),
     &     gObj(nnObj0), gObj2(nnObj0),
     &     gCon(negCon), gCon2(negCon), Hcol(neH), x(n), x1(n),
     &     yCon(nnCon0), y(nb), y1(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s8firstCall  computes the first set of problem functions
!
!     Optionally, the derivatives are verified using finite differences.
!
!     17 Oct 2014: First version based on dnopt routine dnFirstCall.
*     19 Oct 2014: Added user-defined Hessian.
!     ==================================================================
      integer
     &     inform, lvlScale, modefg,  savedlvlScl
!     ------------------------------------------------------------------
      parameter         (lvlScale = 75) ! scale option
!     ------------------------------------------------------------------

!     Scaling is turned off for now.

      savedlvlScl  = iw(lvlScale)
      iw(lvlScale) = 0

      modefg = 2
      call funwrapper
     &   ( inform,
     &     modefg, NonlinearCon, NonlinearObj,
     &     n, negCon, nnCon0, nnCon,
     &     nnJac, nnH, nnObj0, nnObj,
     &     funcon, funobj, userHv,
     &     x, yCon,
     &     neJ, nlocJ, locJ, indJ,
     &     neH, nlocH, locH, indH, Hcol,
     &     fCon, fObj, gCon, gObj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      GotFuns = inform .eq. 0

      if (.not. GotFuns) then
         if (inform .lt. 0) then
            if (FeasibleLC) then
               iExit = 61       ! Undefined fun at first feasible point
            else
               iExit = 62       ! Undefined fun at initial point
            end if
         else
            iExit = inform      ! User wants to stop
         end if
         go to 999
      end if

!     ------------------------------------------------------------------
!     Check derivatives.
!     (One day, we will do this on the SCALED problem.)
!
!     Individual objective derivatives are not checked in FP mode.
!     ------------------------------------------------------------------
      call s7checkG
     &   ( inform,
     &     n, nnCon0, nnCon, nnJac,
     &     nnH, nnObj0, nnObj,
     &     funwrapper, funcon, funobj, userHv,
     &     x, x1, bl, bu, fObj, gObj, yCon,
     &     neJ   , nlocJ, locJ, indJ,
     &     neH   , nlocH, locH, indH, Hcol,
     &     negCon, nlocG, locG,
     &     fCon, gCon, gObj2, fCon2, gCon2,
     &     y, y1, cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (inform .ne. 0) then
         iExit = inform
         go to 999
      end if

!     ------------------------------------------------------------------
!     Compute any missing derivatives.
!     ------------------------------------------------------------------
      call s6getMissing
     &   ( inform,
     &     n, negCon, nnCon0, nnCon, nnJac,
     &     nnH, nnObj0, nnObj,
     &     funwrapper, funcon, funobj, userHv,
     &     bl, bu, x, yCon,
     &     neJ, nlocJ, locJ, indJ,
     &     neH, nlocH, locH, indH, Hcol,
     &     fCon, fObj, gCon, gObj, y1,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )
      if (inform .ne. 0) then
         iExit = inform
         go to 999
      end if

  999 iw(lvlScale) = savedlvlScl
      return

      end ! subroutine s8firstCall

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8SQP
     &   ( iExit,
     &     funwrapper, funcon, funobj, userHv,
     &     mjrLog, mnrLog, snSTOP,
     &     Elastic, GotR, startType,
     &     itn, lenR, m, maxS, mBS, n, nb, nS,
     &     nnCon0, nnCon, nnObj0, nnObj, nnH0, nnH,
     &     majors, minors, nDegen, dualInf,
     &     minimize, iObj, scaleObj, objAdd, fObj, fMerit,
     &     maxVi, maxViRel, supVi,
     &     nInf, sInf, nInfE, sInfE, wtInf0, wtInf,
     &     penParm, piNorm, xNorm,
     &     neJ   , nlocJ, locJ, indJ, Jcol,
     &     neH   , nlocH, locH, indH, Hcol,
     &     negCon, nlocG, locG,
     &     eType, eState, feasType, hs, kBS,
     &     bl, bu, blQP, buQP, blBS, buBS, Fv, Fx,
     &     fCon , gCon , gObj ,
     &     fCon1, gCon1, gObj1,
     &     fCon2, gCon2, gObj2,
     &     gBS, gQP, dyCon, dx, dg,
     &     Udx, HD, Hdx, pBS,
     &     yCon, yCon1, yCon2, pi, QPrhs,
     &     R, rc, rg, rg2, scales,
     &     x, x1, xBS, xQP0, xQP, xPen,
     &     iy, iy1, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     funwrapper, funcon, funobj, userHv,
     &     mjrLog, mnrLog, snSTOP
      logical
     &     Elastic, GotR
      integer
     &     iExit, iObj, itn, lencu, lencw, leniu, leniw, lenru, lenrw,
     &     lenR, maxS, mBS, m, minimize, n, nb, nDegen, negCon,
     &     neH, neJ, nInf, nInfE, nlocG, nlocH, nlocJ, majors, minors,
     &     nnCon0, nnCon, nnH0, nnH, nnObj0, nnObj, nS, startType,
     &     locH(nlocH), locJ(nlocJ), indH(neH), indJ(neJ),
     &     eType(nb), hs(nb), eState(nb), feasType(mBS), locG(nlocG),
     &     kBS(mBS), iy(nb), iy1(nb), iu(leniu), iw(leniw)
      double precision
     &     dualInf, objAdd, fMerit, fObj, maxVi, maxViRel, supVi,
     &     sInf, sInfE, wtInf0, wtInf, piNorm, scaleObj,
     &     bl(nb), bu(nb), blQP(nb), buQP(nb), blBS(mBS), buBS(mBS),
     &     dg(nnH0), dx(nb), dyCon(nnCon0), Fv(nnCon0), Fx(nnCon0),
     &     gBS(mBS), gQP(nnH0), Hcol(neH), HD(nnH0), Hdx(nnH0),
     &     pBS(nb), Jcol(neJ),
     &     fCon(nnCon0) , gCon(negCon) , gObj(nnObj0) ,
     &     fCon1(nnCon0), gCon1(negCon), gObj1(nnObj0),
     &     fCon2(nnCon0), gCon2(negCon), gObj2(nnObj0),
     &     yCon(nnCon0), yCon1(nnCon0), yCon2(nnCon0),
     &     penParm(4), rc(nb), rg(maxS), rg2(maxS), scales(nb),
     &     x(nb), x1(nb), xBS(mBS), xQP(nb), xQP0(nb),
     &     xPen(nnCon0), pi(m), QPrhs(m), R(lenR), Udx(nnH0),
     &     y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s8SQP  solves a nonlinear programming problem.
!     A basis is assumed to be specified by nS, hs, x and the
!     superbasic parts of kBS.
!     In particular, there must be nS values hs(j) = 2, and the
!     corresponding j's must be listed in kBS(m+1) thru kBS(m+ns).
!     The ordering in kBS(m+1:m+nS) matches the reduced Hessian R.
!
!     On entry, if there are nonlinear constraints, Fx contains
!     the true nonlinear slacks (i.e., constraint values)
!     Fx  =  fCon + (linear A)*x,   excluding slacks.
!
!     On exit, if  iExit .lt. 30  it is safe to save the final
!     basis files and print the solution.  Otherwise, a fatal error
!     condition exists and numerous items will be undefined.
!     The last basis map saved (if any) retains the only useful
!     information.
!
!     30 Dec 1991: First version based on npsol routine npcore.
!     23 Oct 1993: Proximal point FP added.
!     29 Oct 1993: Crash on LG rows moved outside s5QP.
!     24 Apr 1994: Nx columns no longer in Q.
!     26 May 1995: Column order of R defined by kBS.
!     04 Aug 1995: Limited memory update
!     11 Aug 1995: tolg changed from 0.1 to 1.0d-4.
!     09 Nov 1995: Updated multipliers used to define Lagrangian.
!     19 Dec 1995: Finite-differences added.
!     09 Oct 1996: First Min Sum version.
!     16 Jul 1997: First thread-safe version.
!     09 Jul 1998: Quasi-Newton updates implemented correctly.
!     24 Aug 1998: Fixed bug in s8x1 found by Alan Brown at Nag.
!     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added.
!     16 Jan 1999: Name changed from s8core.
!     06 Apr 2001: For Hot Starts, don't mess with Utol.
!     27 Apr 2001: wtMax introduced as parameter to s8getWeights.
!     15 Jan 2003: CG and QN  QP solvers added.
!     03 Aug 2003: snEXIT and snPRNT adopted.
!     04 Jul 2005: Switched to vanilla CG for QN with nS > maxR.
!     16 Jun 2008: Call-status implemented correctly.
!     18 Jun 2008: pBS added as argument.
!     07 Oct 2014: infoTags added.
!     19 Oct 2014: Added user-defined Hessian.
!     29 Dec 2014: s8iQP and s8iQN merged to form s8solveQP.
!     26 Apr 2015: blQP, buQP added as arguments of s6linesearch.
!     10 Jun 2015: Removed slack reset after entering elastic mode.
!     ==================================================================
      character
     &     str*132
      external
     &     ddot, dnormi, dnrm1s, dnrm2, s8Hwrapper, s8Hx
      logical
     &     KTcond(2), Boosted, Done, DualFeas,
     &     FDObj, FDCon, Feasible, FPonly,
     &     FirstQP, GoodG, GotNewx, MaxIts, MaxnS, NearOpt,
     &     NeedDerivs, NeedLU, NewB, NewLU, NewTol, NonlinearCon,
     &     NonlinearObj, Optimal, PrimalFeas, PrintLog, PrintSum,
     &     Restart, UseFD, FonlyLS
      integer
     &     cdItns, HvCalls, iAbort, inform,
     &     iPrint, iSumm, itQP, j, jObj, jprimalInf, jdualInf, klog,
     &     kSumm, LUrequest, lvlDer, lvlDif, lvlHess,
     &     lvlPiv, lvlPre, lvlSrch, maxR, minmax, mMajor,
     &     modefg, mjrPrint, mnrPrint, mjrHdP, mjrHdS, mnrHdP,
     &     mRestarts, elastics, nInfQP, nInfEQP, nnJac, nSwap,
     &     outerItn, QNskips, restarts,
     &     PreCon, QPmode, QPsolver, RtRmods, typeLU,
     &     Utol1, Utol2
      integer
     &     infoTag,
     &     QNInfo, MdInfo, LSInfo, FPInfo, QPInfo, FDInfo, HDInfo
      double precision
     &     condZHZ, ddot, dnormi, dnrm1s, dnrm2, dxHdx, eps0,
     &     eps1, eps5, fObj1, fObj2, fObjQP, gMerit,
     &     gNorm0, gNorm, UD0, dRzmax, dRzmin,
     &     primalInf, pHpMrt, sInfQP, sInfEQP, sInfE1, signObj, step,
     &     tolOptFP, tolOptQP, tolOptQPk, tolOptNP,
     &     tolCon, tolx, Utol1s, Utol2s, Utolmn, weight, wtMax,
     &     wtFactor, wtScale, dxNorm, xNorm, xPen0
!     ------------------------------------------------------------------
      integer            SetWeight,       BoostWeight
      parameter         (SetWeight   = 0, BoostWeight   = 1)

      integer            QPChol,     CG
      parameter         (QPChol = 0, CG     = 1)
      integer            HOT
      parameter         (HOT    = 3)
      integer            Unit
      parameter         (Unit   = 2)
      integer            NO
      parameter         (NO     = 0)
      integer            Normal
      parameter         (Normal = 0)
      integer            LM   ,      FM
      parameter         (LM     = 0, FM     = 1)
      integer            B,          BS        , BT
      parameter         (B      = 0, BS     = 2, BT       = 3)
      integer            RedTol,     MinTol,     ResetTol
      parameter         (RedTol = 1, MinTol = 2, ResetTol = 3)
      integer            UnLim ,     VioLim,     UsrLim
      parameter         (UnLim  = 0, VioLim = 1, UsrLim   = 2)

      double precision   zero,          half,           one
      parameter         (zero =0.0d+0,  half   =0.5d+0, one =  1.0d+0)
      double precision   ten
      parameter         (ten  =10.0d+0)

      parameter         (Utol1   = 154) ! abs tol for small diag of U.
      parameter         (Utol2   = 155) ! rel tol for small diag of U.
      parameter         (lvlDif  = 181) ! =1(2) forwd (cntrl) diffs
      parameter         (HvCalls = 188) ! number of Hx products
      parameter         (QPmode  = 208) ! Current QP solver
      parameter         (PreCon  = 209) ! Current precon mode
      parameter         (mnrHdP  = 223) ! >0 => Mnr heading for iPrint
      parameter         (mjrHdP  = 224) ! >0 => Mjr heading for iPrint
      parameter         (mjrHdS  = 226) ! >0 => Mjr heading for iSumm

      parameter         (infoTag = 237) ! infoTag(1:7) status info array
      parameter         (QNInfo  = 237) ! (1): QN update type
      parameter         (MdInfo  = 238) ! (2): QN mod type (A or A+B)
      parameter         (LSInfo  = 239) ! (3): Line search result
      parameter         (FPInfo  = 240) ! (4): QP Feasibility status
      parameter         (QPInfo  = 241) ! (5)  QP Optimality  status
      parameter         (FDInfo  = 242) ! (6)
      parameter         (HDInfo  = 243) ! (7): Approx Hess type
!     ------------------------------------------------------------------
      character          line*4
      data               line/'----'/
!     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file
      iSumm     = iw( 13) ! Summary file

      nnJac     = iw( 21) ! # nonlinear Jacobian variables
      maxR      = iw( 52) ! max columns of R.
      QPsolver  = iw( 55) ! = 0:1:2   => QPChol:CG:QN QP solver
      klog      = iw( 61) ! log/print frequency
      kSumm     = iw( 62) ! Summary print frequency
      lvlHess   = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian
      lvlSrch   = iw( 76) ! >0     => use derivatives in the line search
      lvlPre    = iw( 77) ! >0     => QN preconditioned CG
      lvlPiv    = iw( 80) ! 0/1 Threshold partial/complete pivoting: user

      minmax    = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
      mMajor    = iw( 90) ! limit on major iterations
      mjrPrint  = iw( 92) ! Major print level
      mnrPrint  = iw( 93) ! Minor print level
      lvlDer    = iw( 70) ! = 0, 1, 2 or 3, the derivative level

!     Constants

      eps0      = rw(  2) ! eps**(4/5)          IEEE DP  3.00e-13
      eps1      = rw(  3) ! eps**(2/3)          IEEE DP  3.67e-11
      eps5      = rw(  7) ! eps**(1/5)          IEEE DP  7.40e-04

      tolOptFP  = rw( 51) ! Minor Phase 1 Opt tol
      tolOptQP  = rw( 52) ! Minor Phase 2 Opt tol
      tolOptNP  = rw( 53) ! Major Optimality tolerance
      tolx      = rw( 56) ! Minor feasibility tolerance.
      tolCon    = rw( 57) ! Major feasibility tolerance.
      xPen0     = rw( 89) ! initial penalty parameter.
      wtMax     = rw( 90) ! max     infeasibility weight

      NonlinearCon = nnCon  .gt. 0
      NonlinearObj = nnObj  .gt. 0
      FPonly       = minmax .eq. 0

      iw(HvCalls)  = 0
      iw(mnrHdP)   = 0
      iw(mjrHdP)   = 0
      iw(mjrHdS)   = 0

      iw(QPmode)   = QPsolver    ! Current QP solver
      iw(PreCon)   = lvlPre      ! Current precon mode

!     ------------------------------------------------------------------
!     s8SQP  operates in either ``Normal'' or ``Elastic'' mode.
!     In elastic mode, the nonlinear slacks are allowed to be infeasible
!     while a weighted sum of the slack infeasibilities is minimized.
!     ------------------------------------------------------------------
      Elastic        = .false.
      call iload ( nb, 0, eState, 1 )
      elastics       = 0

      nInfE          = 0
      sInfE          = zero
      sInfE1         = zero

      call dcopy ( nb, bl, 1, blQP, 1)
      call dcopy ( nb, bu, 1, buQP, 1)

      nInf      = 0
      sInf      = zero

      iExit     = 0
      LUrequest = 0
      QNskips   = 0
      restarts  = 0
      if (nnH .gt. 0) then
         mRestarts = 2
      else
         mRestarts = 0
      end if
      Restart  = .false.
      RtRmods  = 0

      call iload ( 7, 0, iw(infoTag), 1 )

      signObj = minimize

      gNorm   = zero
      if (iObj .eq. 0) then
         gNorm0 = zero
      else
         gNorm0 = scaleObj
      end if

      primalInf = zero
      dualInf   = zero
      wtInf     = wtInf0
      tolOptQPk = ten*tolOptQP

      gMerit    = zero
      step      = zero

      KTcond(1) =  .false.
      KTcond(2) =  .false.

      FirstQP   = .true.

      condZHZ   = one

      FDObj   = (lvlDer .eq. 0 .or. lvlDer .eq. 2) .and. (nnObj .gt. 0)
      FDCon   = (lvlDer .eq. 0 .or. lvlDer .eq. 1) .and. (nnJac .gt. 0)
      UseFD   =  FDObj         .or. FDCon
      FonlyLS =  UseFD         .or. lvlSrch .eq. 0

      if (mjrPrint .ge. 10  .or.  mnrPrint .ge. 10) then
         PrintLog = iPrint .gt. 0  .and.  klog  .eq. 1
         PrintSum = iSumm  .gt. 0  .and.  kSumm .eq. 1
         if (PrintLog) then
            write(str, 1000) (line, j=1,29)
            call snPRNT( 1, str, iw, leniw )
            write(str, 1010) majors
            call snPRNT( 1, str, iw, leniw )
         end if
         if (PrintSum  .and.  mnrPrint .ge. 10) then
            write(str, 1000) (line, j=1,19)
            call snPRNT( 2, str, iw, leniw )
            write(str, 1010) majors
            call snPRNT( 2, str, iw, leniw )
         end if
      end if

      jObj   = n + iObj

      if (NonlinearCon) then
         call s8InitPen
     &      ( nnCon, penParm, xPen0, xPen, rw, lenrw )
      end if

      if (nS .gt. maxR) then
         iw(QPmode) = CG      ! Use CG
         iw(PreCon) = NO      ! with no preconditioning
         GotR       = .false.
      end if

      call dcopy ( nb, x, 1, xQP, 1 )

      cdItns     = -1
      NeedDerivs = .false.

!     ======================Start of main loop==========================
!     Start of a Major Iteration.
!     ==================================================================
      do outerItn = 0, mMajor

         GotNewx = .false.

!        ===============================================================
!        Repeat                                          (until GotNewx)
!        ===============================================================
  100       minors = 0

!           ============================================================
!           Repeat                 (until an accurate gradient is found)

  200          if (NeedDerivs) then
                  if (UseFD) then
!                    ---------------------------------------------------
!                    Compute any missing derivatives.
!                    ---------------------------------------------------
                     call s6getMissing
     &                  ( iExit, n, negCon,
     &                    nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &                    funwrapper, funcon, funobj, userHv,
     &                    bl, bu, x, yCon,
     &                    neJ, nlocJ, locJ, indJ,
     &                    neH, nlocH, locH, indH, Hcol,
     &                    fCon, fObj, gCon, gObj, y,
     &                    cu, lencu, iu, leniu, ru, lenru,
     &                    cw, lencw, iw, leniw, rw, lenrw )
                     if (iExit .ne. 0) then
                        go to 999 ! break major iteration loop
                     end if
                  end if ! UseFD
                  NeedDerivs = .false.
               end if

               if (NonlinearObj) then
!                 gNorm = dnormi( nnObj, gObj, 1 )
                  gNorm = dnrm1s( nnObj, gObj, 1 ) ! Approximate 2-norm
               end if

               if (NonlinearCon) then
!                 ------------------------------------------------------
!                 Load the scaled Jacobian in J.
!                 Compute the QP right-hand side   QPrhs  =  Jx - fCon.
!                 Find Fx the nonlinear constraint values.
!                 ------------------------------------------------------
                  call s8Gcopy
     &               ( nnCon, nnJac,
     &                 neJ   , nlocJ, locJ, indJ,
     &                 negCon, nlocG, locG, gCon,
     &                 neJ   , nlocJ, locJ, Jcol )
                  call dcopy ( nnCon, fCon, 1, QPrhs, 1 )
                  call s2Aprod
     &               ( Normal, eps0,
     &                 neJ, nlocJ, locJ, indJ, Jcol,
     &                 one, x, nnJac, (-one), QPrhs, nnCon )
!                 ------------------------------------------------------
!                 s8optimizeSlacks finds the nonlinear slacks sN that
!                 minimize the merit function with x(1:n) and yCon held
!                 fixed.  The optimal slacks are loaded into  x(n+1:nb)
!                 and the violations are calculated:
!                  Fv = fCon  + A(linear)x - nonlinear slacks
!                     = Fx                 - sN
!                 ------------------------------------------------------
                  if (.not. Elastic) then
                     call s8getWeights
     &                  ( SetWeight, Boosted, itn, gNorm0, gNorm,
     &                    wtInf0, wtInf, wtMax,
     &                    weight, wtFactor, wtScale, iw, leniw )
                  end if

                  call s8optimizeSlacks
     &               ( n, nb, nnCon, nInfE, sInfE, tolx, wtInf,
     &                 bl, bu, Fv, x, yCon, xPen, Fx )
               end if

!              ---------------------------------------------------------
!              Prepare to (re-)solve the QP subproblem (possibly after
!              the elastic weight has been increased).
!              ---------------------------------------------------------
!              Factorize the basis at x.
!              Compute xQP such that (J -I)*xQP = rhs.

  300          if (FirstQP) then
!                 ------------------------------------------------------
!                 First QP subproblem.
!                 ------------------------------------------------------
                  call s8InitH
     &               ( itn, nnH0, nnH, gNorm0, gNorm, UD0, HD,
     &                 iw, leniw, rw, lenrw )

                  NeedLU = .true.
                  GotR   = .false.
                  nSwap  = 0
                  if (nS .eq. 0) then
                     typeLU = B
                  else
                     typeLU = BS
                  end if

                  Utol1s = rw(Utol1)
                  Utol2s = rw(Utol2)

!                 To avoid an unnecessarily ill-conditioned starting
!                 basis for the first QP, use big singularity tols
!                 (except if it's a Hot Start!).

                  if (startType .eq. HOT) then
                     Utolmn = eps1
                  else
                     Utolmn = eps5
                  end if

                  rw(Utol1) = max( Utol1s, Utolmn )
                  rw(Utol2) = max( Utol2s, Utolmn )

               else
!                 ------------------------------------------------------
!                 Subsequent factorizations.
!                 ------------------------------------------------------
!                 For linearly constrained problems, the factors L, U
!                 and R can be saved as long as a poor x does not force
!                 a new factorization. (Even in this case, R can be
!                 saved if there are no swaps.)

                  NeedLU = NonlinearCon

                  if (restarts .eq. 0) then
                     typeLU = BT

!                    Reset the factor and update tolerances if changed
!                    during the previous major iteration.

                     call s2tols
     &                  ( ResetTol, NewTol, itn, iw, leniw, rw, lenrw )
                  else
                     typeLU = BT
                     needLU = .true.
                     call s2tols
     &                  ( RedTol, NewTol, itn, iw, leniw, rw, lenrw )
                  end if
               end if

               call s2Bfac
     &            ( iExit, typeLU, NeedLU, NewLU, NewB,
     &              iObj, itn, mjrPrint, LUrequest,
     &              m, mBS, n, nb, nnH, nS, nSwap,
     &              neJ, nlocJ, locJ, indJ, Jcol,
     &              kBS, hs, blQP, buQP, blBS, buBS,
     &              nnCon0, nnCon, QPrhs, xQP, xBS,
     &              iy, iy1, y, y1, iw, leniw, rw, lenrw )

               if (iExit .ne. 0) then
                  go to 999     ! break major iteration loop
               end if

               if (iw(QPmode) .eq. QPChol) then
                  GotR = GotR  .and.  .not. NewLU
               end if

               NeedLU  = .false.
               if (mjrPrint .ge. 10) iw(mjrHdP) = 1

!              xQP satisfies the general constraints.

               if (FirstQP) then
                  iw(80)    = lvlPiv ! Reset original TPP or TCP
                  rw(Utol1) = Utol1s
                  rw(Utol2) = Utol2s
               end if

!              ---------------------------------------------------------
!              Solve the QP subproblem to obtain kBS, xQP and pi.
!              The search direction will be dx = xQP - x.
!              Use x1 to store the first feasible point.
!              ---------------------------------------------------------
               inform = 0

               call s8solveQP
     &            ( inform,
     &              Mnrlog, s8Hwrapper, s8Hx, iw(HvCalls),
     &              Elastic, GotR,
     &              itn, itQP, lenR, m, maxS, mBS, n, nb,
     &              nnCon0, nnCon, nnObj0, nnObj,
     &              nnH0, nnH, nS, nDegen,
     &              mjrPrint, mnrPrint, minimize,
     &              iObj, scaleObj, (objAdd+fObj), fObjQP,
     &              condZHZ, tolOptFP, tolOptQPk, tolx,
     &              nInfQP, sInfQP, elastics, nInfEQP, sInfEQP, wtInf,
     &              UD0, piNorm,
     &              neJ, nlocJ, locJ, indJ, Jcol,
     &              neH, nlocH, locH, indH, Hcol,
     &              eType, eState, feasType, hs, kBS,
     &              bl, bu, blQP, buQP, blBS, buBS,
     &              gBS, gQP, gObj, HD, Hdx,
     &              pBS, pi, R, rc, rg, rg2, QPrhs, scales,
     &              x, xBS, xQP0, xQP,
     &              iy, iy1, y, y1, y2,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )

!              inform    Status
!              ------    ------
!               >0       Fatal error
!                0       QP solution found
!               -1       Too many iterations
!               -2       Too many superbasics

               minors = minors + itQP
               if (inform .gt. 0) then
                  iExit = inform
                  go to 999     ! break major iteration loop
               end if

!              QP inform values are saved until after printing

               MaxnS  = inform .eq. -2
               MaxIts = inform .eq. -1

               FirstQP = .false.

!              ---------------------------------------------------------
!              In normal mode, the routine s8optimizeSlacks defines
!              feasible nonlinear slacks x(n+1:n+nnCon).
!              If the QP has just switched to elastic mode, the slacks
!              are no longer inside their bounds and x(n+1:n+nnCon)
!              can be defined as the true infeasibilities.
!              ---------------------------------------------------------
               if (NonlinearCon) then

!                 Compute the dual search direction.
!                 Set        dyCon = yCon - pi(qp)

                  call dcopy ( nnCon,            pi, 1, dyCon, 1 )
                  call daxpy ( nnCon, (-one),  yCon, 1, dyCon, 1 )
               end if

!              Compute the search direction dx.

               call dcopy ( nb,         xQP, 1, dx, 1 )
               call daxpy ( nb, (-one), x  , 1, dx, 1 )

               xNorm  = dnormi( n, x , 1 )
               dxNorm = dnormi( n, dx, 1 )

!              ---------------------------------------------------------
!              Test for convergence.
!              ---------------------------------------------------------
!              Compute all the QP reduced costs.
!              (We could use yCon for the nonlinear pi's).
!              Compute the maximum dual infeasibility.

               call s8rc
     &            ( scaleObj, minimize, iObj, m, n, nb,
     &              nnObj0, nnObj, nnCon, nnJac, negCon,
     &              neJ, nlocJ, locJ, indJ, Jcol,
     &              gObj, gCon, pi, rc )
               call s8Infs
     &            ( Elastic, n, nb, nnCon0, nnCon, tolx, wtInf,
     &              primalInf, dualInf, jprimalInf, jdualInf,
     &              bl, bu, Fx, rc, x )

               if (GotR) then
                  call s6Rcnd
     &               ( maxR, nS, lenR, R, dRzmax, dRzmin, condZHZ )
               end if

               primalInf  = primalInf/max(xNorm ,   one)
               dualInf    = dualInf  /max(piNorm, gNorm)

               PrimalFeas = primalInf  .le. tolCon
               DualFeas   = dualInf    .le. tolOptNP
               KTcond(1)  = PrimalFeas
               KTcond(2)  = DualFeas

               Optimal  = DualFeas
               Feasible = PrimalFeas
               Done     = (Feasible   .and. Optimal     ) .or.
     &                    (Optimal    .and. Elastic     ) .or.
     &                    (Feasible   .and. FPonly      )
               NearOpt  =  .not. Done .and.
     &                    ((Feasible  .and. dualInf   .lt. ten*tolOptNP)
     &                .or. (Optimal   .and. primalInf .lt. ten*tolCon)
     &                .or. ((FPonly   .or.  dualInf   .lt. ten*tolOptNP)
     &                                .and. primalInf .lt. ten*tolCon ))

               if (Elastic  .and.  Done) then
                  call s8getWeights
     &               ( BoostWeight, Boosted, itn, gNorm0, gNorm,
     &                 wtInf0, wtInf, wtMax,
     &                 weight, wtFactor, wtScale, iw, leniw )
                  if (Boosted) then ! Solve the QP again
                     NeedDerivs = .false.
                     go to 300
                  end if
               end if

               if (majors .eq. 0) then
                  tolOptQPk = min( dualInf, 1.0d-3 )
               end if
               if (minors .eq. 0  .and.  tolOptQPk .gt. tolOptQP) then
                  tolOptQPk = 0.2d+0*tolOptQPk
               end if
               tolOptQPk = min( 0.5d+0*tolOptQPk, 0.1d+0*dualInf )
               tolOptQPk = max( tolOptQPk, tolOptQP )

!              ---------------------------------------------------------
!              Compute the current augmented Lagrangian merit function.
!              objAdd is added in the log routine.
!              ---------------------------------------------------------
               if (iObj .eq. 0) then
                  fMerit = zero
               else
                  fMerit = signObj*x(jObj)*scaleObj
               end if

               if (NonlinearObj) then
                  fMerit = fMerit + signObj*fObj
               end if

               if (NonlinearCon) then
                  call dcopy ( nnCon, Fv  , 1, y, 1 )
                  call ddscl ( nnCon, xPen, 1, y, 1 )
                  fMerit = fMerit -      ddot( nnCon, yCon, 1, Fv, 1 )
     &                            + half*ddot( nnCon,    y, 1, Fv, 1 )

                  if (Elastic) then
                     fMerit = fMerit + wtInf*sInfE
                  end if
               end if

!              ---------------------------------------------------------
!              If the forward-difference estimate of the reduced gradient
!              of the Lagrangian is small,  prepare to: (i) switch to
!              central differences; (ii)  recompute the derivatives,  and
!              (iii) solve the QP again.
!
!              On the other hand, if central differences give a large
!              reduced-gradient norm, switch back to forward differences.
!              ---------------------------------------------------------
               call s8fdSwitch
     &            ( nnCon0, nnCon, nnObj, itn, cdItns,
     &              GoodG, NeedDerivs, UseFD, dualInf,
     &              fCon, fObj, iw, leniw, rw, lenrw )

!              ---------------------------------------------------------
!              Print the details of this iteration.
!              ---------------------------------------------------------
               call mjrLog
     &            ( iAbort,
     &              KTcond, mjrPrint, minimize,
     &              n, nb, nnCon0, nnObj,
     &              nS, itn, majors, minors, nSwap,
     &              condZHZ, iObj, scaleObj, objAdd,
     &              fMerit, penParm, step,
     &              primalInf, dualInf, maxVi, maxViRel, hs,
     &              neJ, nlocJ, locJ, indJ, Jcol,
     &              scales, bl, bu, fCon, yCon, x,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )

               call snSTOP
     &            ( iAbort,
     &              KTcond, mjrPrint, minimize,
     &              m, maxS, n, nb, nnCon0, nnCon, nnObj0, nnObj, nS,
     &              itn, majors, minors, nSwap,
     &              condZHZ, iObj, scaleObj, objAdd,
     &              fMerit, penParm, step,
     &              primalInf, dualInf, maxVi, maxViRel, hs,
     &              neJ, nlocJ, locJ, indJ, Jcol, negCon,
     &              scales, bl, bu, fCon, gCon, gObj,
     &              yCon, pi, rc, rg, x,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )

               if (iAbort .ne. 0) then
                  iExit = 74    ! User has aborted the run via snSTOP
                  go to 999     ! break major iteration loop
               end if

               iw(LSInfo) = UnLim

!+          until (.not. (UseFD  .and.  .not.GoodG))
            if (          UseFD  .and.  .not.GoodG ) go to 200
!           ============================================================

            if (MaxIts                                   ) iExit = 31
            if (                      majors .ge. mMajor ) iExit = 32
            if (MaxnS     .and.       majors .ge. mMajor ) iExit = 33
            if (Feasible  .and.   iw(QPInfo) .eq.      3 ) iExit = 21
            if (NearOpt   .and.        iExit .ne.      0
     &                    .and. (             .not.UseFD
     &                     .or.   iw(lvlDif) .eq.      1)) iExit =  3
!           if (Optimal   .and.        nInfE .gt.      0 ) iExit = 13
            if (Optimal   .and.                  Elastic ) iExit = 13
            if (Feasible  .and.                  Optimal ) iExit =  1
            if (Feasible  .and.                   FPonly ) iExit =  2

            if (iExit .ne. 0) then
               go to 999        ! break major iteration loop
            end if

            step  = zero
            nSwap = 0

!           ============================================================
!           Take a step in the right direction.
!           ============================================================
!           Compute  dxHdx = s'Hs  and other directional derivatives.
!           Be prepared to fix up pHpMrt if there are linear variables.
!           ------------------------------------------------------------
            if (nnH .gt. 0) then
               call s8xHx
     &            ( nnH, dx, Udx, Hdx, dxHdx, iw, leniw, rw, lenrw )
            else
               dxHdx = zero
            end if

            if (nnH .eq. n) then
               if (dxHdx .eq. zero  .and.  iw(HDInfo) .ne. Unit) then
                  call s8ResetH
     &               ( itn, nnH, UD0, HD, iw, leniw, rw, lenrw )
                  needDerivs = .false.
                  go to 100
               end if
               pHpMrt = dxHdx
            else
               pHpMrt = max( eps1*dxNorm*dxNorm, dxHdx )
            end if

!           ------------------------------------------------------------
!           Compute the contributions to the merit function and its
!           directional derivative from the nonlinear constraints.
!           The penalty parameters  xPen(j)  are increased if the
!           directional derivative is not sufficiently negative.
!           ------------------------------------------------------------
!           First, compute the value and directional derivative of the
!           Lagrangian with respect to x and the multipliers.

            if (iObj .eq. 0) then
               fMerit = zero
               gMerit = zero
            else
               fMerit = signObj* x(jObj)*scaleObj
               gMerit = signObj*dx(jObj)*scaleObj
            end if

            if (NonlinearObj) then
               fMerit = fMerit + signObj*fObj
               gMerit = gMerit + signObj*ddot( nnObj, gObj, 1, dx, 1 )
            end if

            if (Elastic) then
               fMerit = fMerit +            sInfE *wtInf
               gMerit = gMerit + (sInfEQP - sInfE)*wtInf
            end if

!           ------------------------------------------------------------
!           Compute the search direction for the multipliers and nonlinear
!           slacks, and the contributions to the merit function and its
!           directional derivative from the nonlinear constraints.
!           The penalty parameters  xPen(j)  are increased if the
!           directional derivative is not sufficiently negative.
!           ------------------------------------------------------------
            if (NonlinearCon) then
               fMerit = fMerit  - ddot( nnCon,  yCon, 1, Fv, 1 )
               gMerit = gMerit  + ddot( nnCon,  yCon, 1, Fv, 1 )
               gMerit = gMerit  - ddot( nnCon, dyCon, 1, Fv, 1 )

               call s8Merit
     &            ( nnCon, fMerit, gMerit, pHpMrt, penParm,
     &              Fv, xPen, y, rw, lenrw )
            end if

!           ============================================================
!           Compute the scalar step, the step length from x along dx.
!           Implicitly,  dx defines the change sInfEQP - sInfE in the
!           elastic variables.
!
!           x   is the base point, with associated  sInfE
!           x1  is x + step*dx,    with associated  sInfE1
!
!           pBS == x2 is workspace.
!           ============================================================
            call s6lineSearch
     &         ( inform, funwrapper, funcon, funobj, userHv,
     &           Elastic, FonlyLS, PrimalFeas,
     &           iObj, signObj, scaleObj,
     &           n, nb, nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &           itn, majors, QNskips, maxVi, maxViRel, supVi,
     &           step, dxNorm, xNorm, fMerit, gMerit,
     &           sInfE, sInfEQP, sInfE1, wtInf,
     &           bl, bu, dx, dyCon,
     &           neJ   , nlocJ, locJ, indJ, Jcol,
     &           neH   , nlocH, locH, indH, Hcol,
     &           negCon, nlocG, locG,
     &           fObj1, fCon1, gCon1, gObj1,
     &           fObj2, fCon2, gCon2, gObj2, Fx,
     &           x, x1, pBS, xQP, pi,
     &           yCon, yCon1, yCon2, xPen,
     &           y, y1, y2, cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

            Restart = restarts .lt. mRestarts ! permission to restart

            if (inform .eq. 0) then

               GotNewx = .true.

            else if (inform .gt. 0) then
               iExit = inform   ! The user wants to stop or violn limit
               go to 999        ! break major iteration loop

            else
!              =========================================================
!              No acceptable step
!              =========================================================
               if (NonlinearCon  .and.  sInfE .gt. zero) then
                  call s8getWeights
     &               ( BoostWeight, Boosted, itn, gNorm0, gNorm,
     &                 wtInf0, wtInf, wtMax,
     &                 weight, wtFactor, wtScale, iw, leniw )
                  if (Boosted) then
                     Elastic = .true.
                  end if
               end if

!              =========================================================
!              As yet, unknown reason for no acceptable step,
!              Deal with some obvious cases.
!              =========================================================
               if (UseFD  .and.  iw(lvlDif) .ne. 2) then
!                 ------------------------------------------------------
!                 Switch to central differences and solve the QP again.
!                 ------------------------------------------------------
                  cdItns     = 0
                  write(str, 3020) itn
                  call snPRNT( 23, str, iw, leniw )
                  iw(lvlDif) = 2
                  iw(FDInfo) = 2
               else
                  if      (MaxnS  ) then
                     iExit = 33 ! Superbasics limit
                  else if (NearOpt) then
                     iExit =  3 ! Requested accuracy could not be ...
                  end if

                  if (iExit .gt. 0) then
                     go to 999  ! break major iteration loop
                  end if

                  if (Restart) then
                     if (iw(HDInfo) .ne. Unit) then

!                       Discard the off-diagonals of H.
!                       If H is already diagonal, H is set to the identity.

                        call s8ResetH
     &                     ( itn, nnH, UD0, HD, iw, leniw, rw, lenrw )
                        GotR  = .false.
                     end if

                     if (NonlinearCon) then
                        call s8InitPen
     &                     ( nnCon, penParm, xPen0, xPen, rw, lenrw )
                     end if
                     restarts = restarts + 1

                  else
!                    ---------------------------------------------------
!                    We have run out of things to try. Bummer.
!                    ---------------------------------------------------
                     iExit = 41 ! Current point cannot be improved...
                     go to 999  ! break major iteration loop
                  end if
               end if
            end if

            if (.not. GotNewx) goto 100
!+       until       (GotNewx)

!        ===============================================================
!        The new point  x1  has been computed.
!        ===============================================================
         restarts = 0
         inform   = 0

!        ---------------------------------------------------------------
!        Some unknown derivatives may need to be calculated at x1.
!        ---------------------------------------------------------------
         if (FonlyLS  .and.  nnH .gt. 0)  then
            modefg = 1
            call funwrapper
     &         ( iExit,
     &           modefg, NonlinearCon, NonlinearObj,
     &           n, negCon, nnCon0, nnCon,
     &           nnJac, nnH, nnObj0, nnObj,
     &           funcon, funobj, userHv,
     &           x1, yCon1,
     &           neJ, nlocJ, locJ, indJ,
     &           neH, nlocH, locH, indH, Hcol,
     &           fCon2, fObj2, gCon1, gObj1,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

            if (iExit .ne. 0) then
               go to 999        ! break major iteration loop
            end if

            if (UseFD) then
               call s6getMissing
     &            ( iExit, n, negCon,
     &              nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &              funwrapper, funcon, funobj, userHv,
     &              bl, bu, x1, yCon1,
     &              neJ, nlocJ, locJ, indJ,
     &              neH, nlocH, locH, indH, Hcol,
     &              fCon1, fObj1, gCon1, gObj1, y,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )

               if (iExit .ne. 0) then
                  go to 999     ! break major iteration loop
               end if
            end if
         end if

         inform = 0
         majors = majors + 1
         if (iw(lvlDif) .eq. 2) then
            cdItns = cdItns + 1
         end if

         if (mjrPrint .ge. 10  .or.  mnrPrint .ge. 10) then
            PrintLog = iPrint .gt. 0  .and.  klog  .eq. 1
            PrintSum = iSumm  .gt. 0  .and.  kSumm .eq. 1

            if (PrintLog) then
               call s1page( 0, iw, leniw )
               write(str, 1000) (line, j=1,29)
               call snPRNT( 1, str, iw, leniw )
               write(str, 1010) majors
               call snPRNT( 1, str, iw, leniw )
            end if
            if (PrintSum  .and.  mnrPrint .ge. 10) then
               write(str, 1000) (line, j=1,19)
               call snPRNT( 2, str, iw, leniw )
               write(str, 1010) majors
               call snPRNT( 2, str, iw, leniw )
            end if
         end if

!        ===============================================================
!        The problem functions have been defined at the new x.
!        ===============================================================
         if ( nnH     .gt. 0  .and.
     &       (lvlHess .eq. LM .or. lvlHess .eq. FM)) then
                                ! Update a QN approximate Hessian.
            call s8HQN
     &         ( inform,
     &           funwrapper, funcon, funobj, userHv,
     &           UseFD, iw(QPmode), startType,
     &           itn, lenR, m, mBS, n, nb,
     &           nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &           nS, majors, QNskips, UD0,
     &           minimize, step, dxHdx,
     &           RtRmods, GotR, penParm,
     &           fObj, fCon, gCon, gObj, fCon1, gCon1, gObj1,
     &           neJ   , nlocJ, locJ, indJ, Jcol,
     &           neH   , nlocH, locH, indH, Hcol,
     &           negCon, nlocG, locG,
     &           kBS, bl, bu, dx, dg, Udx, Hdx, HD,
     &           yCon1, R, x, x1, xQP0, xPen, y, y1, y2,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

            if (inform .ne. 0) then
               iExit = inform
               go to 999        ! break,  major iteration loop
            end if
         end if

!        ---------------------------------------------------------------
!        Update the variables.
!        The QP solution, saved in xQP, is used to start the next QP.
!        (If a unit step was not taken last iteration, some more
!        nonbasics may be between their bounds.
!        Nov 10, 1994. Tried leaving the nonbasics between their
!        bounds after short step. In some cases, the number of minor
!        iterations increased dramatically with a very short step.)
!        ---------------------------------------------------------------
         call dcopy ( nb, x1, 1, x, 1 )

         if (NonlinearCon) then
            call dcopy ( negCon, gCon1, 1, gCon, 1 )
            call dcopy ( nnCon , yCon1, 1, yCon, 1 )
            call dcopy ( nnCon , fCon1, 1, fCon, 1 )
         end if

         if (NonlinearObj) then
            fObj  = fObj1
            call dcopy ( nnObj, gObj1, 1, gObj, 1 )
         end if

         sInfE = sInfE1

      end do
!     ======================end of main loop============================

  999 return

 1000 format(1x, 29a4)
 1010 format(' Start of major itn', i6)
 3020 format(' Itn', i7, ' -- Central differences invoked.',
     &       ' Small step length.' )

      end ! subroutine s8SQP

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8callStatus
     &   ( Status, iw, leniw)

      implicit
     &     none

      integer
     &     Status, leniw, iw(leniw)

!     ==================================================================
!     s8callStatus fetches the call-status for the snOpt user-defined
!     functions.
!
!     16 Jun 2008: First version of s8callStatus.
!     ==================================================================
      character
     &     str*80
      integer
     &     npStatus
!     ------------------------------------------------------------------
      parameter         (npStatus = 236) ! NP user-routine call-status
!     ------------------------------------------------------------------

      if (     iw(npStatus) .eq. 0) then
         ! Standard call

         Status       =  0
      else if (iw(npStatus) .lt. 0) then
         ! First call

         Status       =  1
         iw(npStatus) =  0
      else if (iw(npStatus) .ge. 2) then
         ! Last orders please

         Status       = iw(npStatus)
         iw(npStatus) = -1
      else
         Status       = iw(npStatus)
         write(str, 9999) Status
         call snPRNT( 3, str, iw, leniw )
      end if

      return

 9999 format(' XXX  user-function call-status not recognized.',
     &       ' Requested status =', i6)

      end ! subroutine s8callStatus
