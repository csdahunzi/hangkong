!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn80ncon.f
!
!     s8fdSwitch        s8Fv         s8Fx
!     s8Gcopy           s8getFeasLC  s8getR      s8getWeights
!     s8Gloc            s8Gprod      s8Gsize
!     s8Infs            s8InitH      s8InitPen
!     s8Merit
!     s8optimizeSlacks
!     s8HxLP            s8HxPP       s8HxQP      s8HxNull
!     s8rand            s8rc
!     s8scaleG          s8scaleJ     s8sInf      s8solveQP
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8fdSwitch
     &   ( nnCon0, nnCon, nnObj, itn, cdItns,
     &     GoodG, NeedDerivs, UseFD, dualInf,
     &     fCon, fObj, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     GoodG, NeedDerivs, UseFD
      integer
     &     nnCon0, nnCon, nnObj, itn, cdItns, leniw, lenrw, iw(leniw)
      double precision
     &     dualInf, fCon(nnCon0), fObj, rw(lenrw)

!     ==================================================================
!     s8fdSwitch   controls the switch between forward and central
!     differences.
!
!     If the forward-difference estimate of the reduced gradient of the
!     Lagrangian is small,  a switch is made to central differences.
!     In this case, the derivatives are recomputed and the QP is solved
!     again.
!
!     On the other hand, if central differences have produced a large
!     reduced-gradient norm, switch back to forward differences.
!
!     31 Mar 2000: First version of s8fdSwitch written for SNOPT 6.1.
!     03 Aug 2003: snPRNT adopted.
!     03 Aug 2003: Current version of s8fdSwitch.
!     ==================================================================
      character
     &     str*80
      external
     &     dnrm1s
      logical
     &     Central
      integer
     &     FDInfo, lvlDif
      double precision
     &     epsrf, cNorm, fdint1, ObjSize, rgNorm, rgTest, dnrm1s
!     ------------------------------------------------------------------
      parameter         (lvlDif = 181) ! forwd diffs or cntrl diffs
      parameter         (FDInfo = 242) ! infoTag(6)

      double precision   zero,          one,          ten
      parameter         (zero = 0.0d+0, one = 1.0d+0, ten = 10.0d+0)
!     ------------------------------------------------------------------
      epsrf      = rw( 73) ! relative function precision.
      fdint1     = rw( 76) ! (1) forwrd diff. interval

      Central    = iw(lvlDif) .eq. 2

      if (nnCon .eq. 0) then
         cNorm   = zero
      else
         cNorm   = dnrm1s( nnCon, fCon, 1 )
      end if

      if (nnObj .eq. 0) then
         ObjSize = zero
      else
         ObjSize = abs(fObj)
      end if

      GoodG  = .true.
      rgTest = (one + ObjSize + cNorm)*epsrf/fdint1
      rgNorm = dualInf

      if (Central) then
         if (rgNorm .gt. ten*rgTest  .and.  cdItns .gt. 0) then
            iw(lvlDif) =  1
            Central    = .false.
            if (UseFD) then
               iw(FDInfo) = 1  ! Forward differences
            end if
         end if
      else
         if (rgNorm .le.     rgTest) then
            cdItns     = 0
            iw(lvlDif) = 2
            if (UseFD) then
               GoodG      = .false.
               NeedDerivs = .true.
               iw(FDInfo) = 2   ! Central differences
               write(str, 1000) itn
               call snPRNT( 23, str, iw, leniw )
            end if
         end if
      end if

 1000 format( ' Itn', i7, ' -- Central differences invoked.',
     &       '  Small reduced gradient.' )

      end ! subroutine s8fdSwitch

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Fv
     &   ( Elastic, n, nnCon, tolz, wtInf,
     &     bl, bu, Fv, x, yCon, Fx )

      implicit
     &     none
      logical
     &     Elastic
      integer
     &     n, nnCon
      double precision
     &     tolz, wtInf, bl(n+nnCon), bu(n+nnCon), x(n+nnCon),
     &     Fx(nnCon), Fv(nnCon), yCon(nnCon)

!     ==================================================================
!     s8Fv  computes the vector of nonlinear constraint violations:
!        Fv = fCon + A(linear)*x - (nonlinear slacks)
!
!     If the Lagrange multiplier is zero, the violation can be set to
!     any value without changing the merit function.  In this case we
!     try and set the slack so that Fv is zero (subject to the slack
!     being feasible).
!
!     In elastic mode we implicitly adjust the variables v and w such
!     that   c - s(feas) + v - w = 0,  with  v >= 0  and  w >= 0.
!
!     On entry,
!        x   =  the current x.
!        Fx  =  fCon + A(linear)*x,   defined in s8Fx.
!
!     On exit,
!        x   =  x containing the modified slacks.
!        Fv  =  fCon + A(linear)*x -  slacks.
!        Fx  =  unaltered.
!
!     19 Apr 2001: First version.
!     19 Apr 2001: Current version.
!     ==================================================================
      integer
     &     i, j
      double precision
     &     blj, buj, Fxi, Fvi, FvL, FvU, xj, yConi, yConv, yConw
!     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
!     ------------------------------------------------------------------
      do i = 1, nnCon
         j     = n + i
         xj    = x (j)
         Fxi   = Fx(i)
         Fvi   = Fxi - xj

         yConi = yCon(i)

         blj   = bl(j)
         buj   = bu(j)

         FvU    = Fxi - buj
         FvL    = Fxi - blj

         yConv = abs( wtInf - yConi ) ! Multiplier for v in elastic mode
         yConw = abs( wtInf + yConi ) ! Multiplier for w in elastic mode

         if (     Elastic .and. xj .le. blj .and. yConv .le. tolz) then
            if (Fvi .gt. zero) then
               Fvi = max( zero, FvL )
            else
               Fvi = zero
            end if
         else if (Elastic .and. xj .ge. buj .and. yConw .le. tolz) then
            if (Fvi .lt. zero) then
               Fvi = min( zero, FvU )
            else
               Fvi = zero
            end if
         else
            if (     yConi .le.  tolz  .and.  Fvi .gt. zero) then
               Fvi = max( zero, FvU )
            else if (yConi .ge. -tolz  .and.  Fvi .lt. zero) then
               Fvi = min( zero, FvL )
            end if
         end if

         xj    = Fxi - Fvi
         Fv(i) = Fvi
         x(j)  = xj

      end do

      end ! subroutine s8Fv

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Fx
     &   ( n, nnCon, nnJac, tolz,
     &     neJ, nlocJ, locJ, indJ, Jcol, fCon, x, Fx )

      implicit
     &     none
      integer
     &     n, nnCon, nnJac, neJ, nlocJ, indJ(neJ), locJ(nlocJ)
      double precision
     &     tolz, Jcol(neJ), x(n+nnCon), fCon(nnCon), Fx(nnCon)

!     ==================================================================
!     s8Fx  defines the nonlinear constraint values
!       Fx  =  true nonlinear slack = fCon + A(linear)*x,
!
!     09 Jan 1992: First version based on Minos routine m8viol.
!     16 Nov 1998: Norm x changed to include only columns.
!     21 Oct 2000: Made compatible with SNOPT 6.1
!     21 Oct 2000: Current version of s8Fx
!     ==================================================================
      integer
     &     nlin
!     ------------------------------------------------------------------
      integer            Normal
      parameter         (Normal = 0)
      double precision   one
      parameter         (one = 1.0d+0)
!     ------------------------------------------------------------------
!     Compute the nonlinear constraint value.
!     Set  Fx  =  fCon + (linear A)*x,   excluding slacks.

      call dcopy ( nnCon, fCon, 1, Fx, 1 )

      nlin = n - nnJac
      if (nlin .gt. 0) then
         call s2Aprod
     &      ( Normal, tolz,
     &        neJ, nlin+1, locJ(nnJac+1), indJ, Jcol,
     &        one, x(nnJac+1), nlin, one, Fx, nnCon )
      end if

      end ! subroutine s8Fx

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Gcopy
     &   ( nnCon, nnJac, neJ, nlocJ, locJ, indJ,
     &     neG1, nlocG1, locG1, G1,
     &     neG2, nlocG2, locG2, G2 )

      implicit
     &     none
      integer
     &     nnCon, nnJac, neG1, neG2, nlocG1, nlocG2, neJ, nlocJ,
     &     indJ(neJ), locJ(nlocJ), locG1(nlocG1), locG2(nlocG2)
      double precision
     &     G1(neG1), G2(neG2)

!     ==================================================================
!     s8Gcopy  copies G1 into G2 when either  G1 or  G2
!     is stored in the upper-left hand corner of J.
!
!     16 Sep 1993: First version.
!     26 Oct 2000: Current version.
!     ==================================================================
      integer
     &     ir, j, k, l1, l2
!     ------------------------------------------------------------------
      do j  = 1, nnJac
         l1 = locG1(j)
         l2 = locG2(j)
         do k  = locJ(j), locJ(j+1)-1
            ir = indJ(k)
            if (ir .gt. nnCon) go to 100
            G2(l2) = G1(l1)
            l1  = l1 + 1
            l2  = l2 + 1
         end do
  100    continue
      end do

      end ! subroutine s8Gcopy

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8getFeasLC
     &   ( iExit,
     &     startType, mnrLog, lenR, m, maxS, mBS,
     &     n, nb, nnCon0, nnCon, nnH0, nnH, nDegen, nS,
     &     numLC, numLIQ, itn, itnlim, itQP, mnrPrtlvl,
     &     scaleObj, tolOptQP, tolx,
     &     nInf, sInf, nInfE, sInfE, wtInf, piNorm, rgNorm,
     &     neJ, nlocJ, locJ, indJ, Jcol,
     &     neH, nlocH, locH, indH, Hcol,
     &     eType, eState, feasType, hs, kBS,
     &     bl, bu, blQP, buQP, blBS, buBS, blSave, buSave,
     &     gBS, gQP, Hdx, pBS, pi,
     &     R, rc, rg, QPrhs, scales,
     &     x0, x, xBS,
     &     iy, iy1, y, y1, y2,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     mnrLog
      integer
     &     iExit, lenR, m, maxS, mBS, n, nb, nnCon0, nnCon, neJ, neH,
     &     nlocJ, nlocH, nnH0, nnH, nDegen, nS, numLC, numLIQ,
     &     nInf, nInfE, itn, itnlim, itQP, mnrPrtlvl, startType,
     &     lencw, leniw, lenrw, locJ(nlocJ), locH(nlocH), indJ(neJ),
     &     indH(neH), kBS(mBS), feasType(mBS),
     &     eState(nb), eType(nb), hs(nb), iy(nb), iy1(nb), iw(leniw)
      double precision
     &     scaleObj, tolOptQP, tolx, sInf, sInfE, wtInf, piNorm,
     &     rgNorm, bl(nb), bu(nb), blQP(nb), buQP(nb),
     &     blSave(nb), buSave(nb), blBS(mBS), buBS(mBS), gBS(mBS),
     &     gQP(nnH0), Hdx(nnH0), Hcol(neH), Jcol(neJ), pBS(nb), pi(m),
     &     QPrhs(nnCon0), R(lenR), rc(nb), rg(maxS), scales(nb),
     &     x0(nb), x(nb), xBS(mBS), y(nb), y1(nb), y2(nb), rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     s8getfeasLC finds a feasible point for the linear constraints.
!
!     On entry:
!     ---------
!     startType     is 0,1,2,3: an integer version of the solver's Start
!                   (character for sqopt, snoptb, snoptc, npopt
!                   integer    for snopta).
!
!     nInf          is the number of infeasible equality rows.
!
!     nInfE         is zero (elastic mode is turned off in s5getB).
!
!     kBS(m+1:m+nS) contains a basis specified by nS, hs(*) and x(*).
!                   In particular, there must be nS values hs(j) = 2,
!                   and the corresponding j's must be listed in
!                   kBS(m+1:m+nS). If nInf = 0, the basis is infeasible.
!
!     bl, bu        contain the true bounds (optionally scaled).
!
!     blSave, buSave are copies of the (optionally scaled) upper
!                   and bounds set in s5getB. Entries for the nonlinear
!                   rows are relaxed.
!
!     blQP, buQP    contain the working lower and upper bounds.
!
!
!      iExit       Result
!      -----       ------
!       >0         Fatal error
!        0         Feasible point found
!
!     11 May 1994: First version of s8feasLC.
!     19 Aug 1996: First minsum version.
!     05 Feb 1998: Proximal point norm changed to one-norm.
!     23 Dec 1999: Optional Proximal Point methods 0 and  2 added.
!     03 Aug 2003: snPRNT and snEXIT adopted.
!     16 May 2006: Explicit target itQP added.
!     18 Jun 2008: Hdx, pBS and rg added as arguments.
!     22 Feb 2015: First FP call to s5QP instead of s5LP (for nS > 0).
!     03 Mar 2015: Set nnH = 0 in first FP call.
!     ==================================================================
      character
     &     probTag*20, str*80
      logical
     &     Elastic, GotR, NeedLU, Needx
      integer
     &     eigH, HvCalls, inform, iObjPP, itQPmax, itQPtarget,
     &     j, lvlObjE, lvlPPm, eMode, minimize, mNewSB, mnrHdP, mnrHdS,
     &     mSBsave, elastics, nnHFP, nObjFP0, nObjFP, nObjPP0, nObjPP,
     &     nviol, probType, subOptimize, typeLU
      double precision
     &     condZmax0, condZmax, eps0, eps2, objA, objPP, x0j,
     &     tolOptFP, tolOptPP, Hcondbnd
      external
     &     s8HxPP, s8HxQP, s8HxNull
!     ------------------------------------------------------------------
      integer            COLD
      parameter         (COLD     = 0)
      integer            FP,           LP,          QPP
      parameter         (FP       = 0, LP      = 1, QPP = 6)
      integer            SEMDEF,       POSDEF
      parameter         (SEMDEF   = 0, POSDEF  = 1)
      integer            BT
      parameter         (BT       = 3)
      integer            Normal
      parameter         (Normal   = 0)
      integer            NoSubOpt
      parameter         (NoSubOpt =-1)

      parameter         (mNewSB =  95) ! max # of new superbasics
      parameter         (mnrHdP = 223) ! >0 => Minor heading for iPrint
      parameter         (mnrHdS = 225) ! >0 => Minor heading for iSumm

      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      lvlPPm    = iw( 79) ! 1(2)-norm proximal point method for x0

      eps0      = rw(  2) ! eps**(4/5)
      eps2      = rw(  4) ! eps**(1/2)
      Hcondbnd  = rw( 85) ! bound on the condition of ZHZ
      condZmax0 = rw( 86) ! initial bound on the condition est of Z

      iExit     = 0

      objA      = zero
      iObjPP    = 0
      minimize  = 1       ! Local value

      NeedLU    = .true.
      Needx     =  NeedLU

      condZmax  = condZmax0

!     Set the QP rhs to make x satisfy the (relaxed) nonlinear rows.
!     The array  QPrhs  contains the rhs.
!     Use a fairly tight optimality tolerance for phase 1.

      if (nnCon .gt. 0) then
         call dcopy ( nnCon, x(n+1), 1, QPrhs, 1 )
         call s2Aprod
     &      ( Normal, eps0, neJ, nlocJ, locJ, indJ, Jcol,
     &        one, x, n, (-one), QPrhs, nnCon )
      end if

      if (numLIQ .gt. 0  .or.  nInf .gt. 0) then
!        ---------------------------------------------------------------
!        s5getB has found a basis for the linear constraints that may
!        or may not be feasible.  Find a feasible point for all the
!        linear constraints.  If none exists, minimize the sum of
!        infeasibilities of the linear rows, subject to the bounds.
!        ---------------------------------------------------------------
!        Make the bounds on the linear rows elastic.

!        The bounds on the linear rows are elastic

         call iload ( numLC, 3, eType(n+nnCon+1), 1 )

         if (nInf .eq. 0) then

!           s5getB found feasible linear equality constraints.

            Elastic  = .false.  ! Start in normal mode
            eMode    = 1        ! Enter elastic mode if infeasible
            probType = FP
         else

!           s5getB found infeasible linear equality constraints.
!           Minimize the sum of infeasibilities.

            Elastic  = .true.   ! Start in elastic mode
            eMode    = 2        ! elastic mode throughout
            elastics = 0
            probType = LP
         end if

         lvlObjE     = 2        ! In elastic mode, use:
                                ! W1 = 0, W2 = 1   W1*trueObj + W2*sInfE
         probTag     = 'linear rows'
         subOptimize = NoSubOpt
         tolOptFP    = max( 1.0d-6, eps2 ) ! Default optimality tol.

         eigH        =  SEMDEF

         iw(mnrHdP)  = 1        ! Switch to QP print   heading
         iw(mnrHdS)  = 1        ! Switch to QP summary heading

         nObjFP      = 0        ! No explicit gradient in proximal point
         nObjFP0     = 1
         nnHFP       = 0        ! Ignore any Hessian
         GotR        = .false.
         typeLU      = BT
         HvCalls     = 0

         call s5QP
     &      ( inform,
     &        probType, probTag, subOptimize,
     &        mnrLog, s8HxQP, s8HxNull, HvCalls, eigH,
     &        Elastic, GotR, NeedLU, typeLU, Needx,
     &        lenR, m, maxS, mBS, n, nb, nDegen,
     &        nnH0, nnHFP, nObjFP0, nObjFP, nnH0, nnHFP, nS,
     &        itQP, itnlim, itnlim, itn,
     &        eMode, lvlObjE, mnrPrtlvl,
     &        minimize, iObjPP, scaleObj, objA, objPP,
     &        condZmax, Hcondbnd, tolOptFP, tolOptQP, tolx,
     &        nInf, sInf, elastics, nInfE, sInfE, wtInf,
     &        piNorm, rgNorm,
     &        neJ, nlocJ, locJ, indJ, Jcol,
     &        neH, nlocH, locH, indH, Hcol,
     &        eType, eState, feasType, hs, kBS,
     &        bl, bu, blQP, buQP, blBS, buBS,
     &        gBS, gQP, gQP, Hdx, pBS, pi, R, rc, rg,
     &        nnCon0, nnCon, QPrhs, scales,
     &        nnH0, nnH, x0, x, xBS, x0,
     &        iy, iy1, y, y1, y2,
     &        cw, lencw, iw, leniw, rw, lenrw,
     &        cw, lencw, iw, leniw, rw, lenrw )

!        Check for trouble in s5QP.
!        iExit        Status
!        -----        ------
!         -3          Too many iterations
!         -2          Phase 1 is unbounded
!         -1          infeasible nonelastics
!          0          infeasibilities minimized
!         >0          Fatal error
!
!        Time to stop if the linear constraints are infeasible.
!        If inform = 0, the sum of infeasibilities will have been
!        minimized.

         if (inform .ne. 0  .or.  nInf .gt. 0  .or.  nInfE .gt. 0) then
            if (inform .gt. 0) then
               iExit = inform   ! Fatal error
            else if (inform .eq. -3) then
               iExit = 31       ! iterations limit
            else if (nInf   .gt. 0) then
               iExit = 11       ! infeasible linear constraints
            else if (nInfE  .gt. 0) then
               iExit = 14       ! infeasible linear constraints
            end if
            if (iExit .ne. 0) go to 800
         end if

!        Print something brief if s5LP didn't already do so.

         if (mnrPrtlvl .ge. 1) then
            write(str, 8000) itn
            call snPRNT( 31, str, iw, leniw )
            call snPRNT( 22, str, iw, leniw )
         end if

         needLU = .false.

      end if

!     ------------------------------------------------------------------
!     x is feasible for the linear constraints.
!     The linear rows are not allowed to go infeasible again.
!     ------------------------------------------------------------------
      call iload ( numLC, 0, eType(n+nnCon+1), 1 )

      if (lvlPPm .gt. 0 .and. nnH .gt. 0 .and. startType .eq. COLD) then
!        ===============================================================
!        Find a feasible point closest to x0.
!        Minimize norm(x - x0).
!        ===============================================================
         if (mnrPrtlvl .ge. 1) then
            write(str, 8100) itn, lvlPPm
            call snPRNT( 23, str, iw, leniw )
         end if

         if (lvlPPm .eq. 1) then
!           ------------------------------------------------------------
!           Minimize the one-norm of (x-x0) by fixing the nonlinear
!           variables so that bl = x0 = bu.  Any bl or bu that is moved
!           to  x0  is made elastic.
!           ------------------------------------------------------------
            do j = 1, nnH
               if (bl(j) .eq. bu(j)) then
!                 Relax
               else
                  x0j      = x0(j)
                  bl(j)    = x0j
                  bu(j)    = x0j
                  blQP(j)  = bl(j)
                  buQP(j)  = bu(j)
                  eType(j) = 3

                  if (hs(j) .le. 1) then
                     x(j) = x0j
                  end if
               end if
            end do

            Elastic    = .false.     ! Start in normal mode
            eMode      = 1           ! Enter elastic mode if infeasible
            lvlObjE    = 2           ! In elastic mode, use:
                                     ! 0*trueObj + 1*sInfE
            probTag     = 'norm(x-x0) problem  '
            iw(mnrHdP)  = 1          ! New LP print   header
            iw(mnrHdS)  = 1          ! New LP summary header
            Needx       = .true.
            subOptimize = NoSubOpt
            tolOptFP    = 1.0d-2     ! Sloppy phase1 optimality tol for PP
            tolOptPP    = 1.0d-2     ! Sloppy phase 2 opt tol
            itQPmax     = itQP + 200 ! Limit the minor iterations

            call s5LP
     &         ( inform,
     &           LP, probTag, Elastic,
     &           subOptimize, mnrLog, NeedLU, Needx,
     &           m, n, nb, nDegen, itQP, itQPmax, itn,
     &           eMode, lvlObjE, mnrPrtlvl,
     &           minimize, iObjPP, scaleObj, objA,
     &           condZmax, tolOptFP, tolOptPP, tolx,
     &           nInf, sInf, elastics, nInfE, sInfE, wtInf,
     &           piNorm, rgNorm,
     &           neJ, nlocJ, locJ, indJ, Jcol,
     &           eType, eState, feasType, hs, kBS,
     &           bl, bu, blQP, buQP, blBS, buBS,
     &           gBS, pi, rc, nnCon0, nnCon, QPrhs, scales,
     &           x, xBS, x0,
     &           iy, iy1, y, y1,
     &           cw, lencw, iw, leniw, rw, lenrw )

!           Some elastic variables may have moved outside their bounds.
!           Count them.  Reset the true bounds.
!           If necessary,  get feasible again with the normal tolOptQP.

            nviol = 0
            do j  = 1, nnH
               bl(j)    = blSave(j)
               bu(j)    = buSave(j)
               blQP(j)  = bl(j)
               buQP(j)  = bu(j)
               eType(j) = 0

               if (x(j) .lt. bl(j) - tolx  .or.
     &             x(j) .gt. bu(j) + tolx      ) then
                  nviol = nviol + 1
               end if
            end do

!           Check for errors in s5LP.
!           inform values are = -3,-2,-1, 0, >0

            if (inform .ne. 0) then
               if (inform .gt. 0) then
                  iExit = inform ! Fatal error
               else if (inform .eq. -3  .and. itn .ge. itnlim) then
                  iExit = 31     ! iterations limit
               end if
               if (iExit .ne. 0) go to 800
            end if

            if (inform .eq. 0  .and.  mnrPrtlvl .ge. 1) then
               write(str, 8200) itn, lvlPPm, sInf
               call snPRNT( 33, str, iw, leniw )
               if (nviol .gt. 0) then
                  write(str, 8300) itn
                  call snPRNT( 23, str, iw, leniw )
               end if
            end if

            if (nviol .gt. 0) then

!              Return the infeasible variables to feasibility.

               probTag     = 'linear rows again   '
               Elastic     = .false.
               NeedLU      = .false.
               Needx       = .true.
               subOptimize = NoSubOpt
               tolOptFP    = eps2   ! Revert to accurate phase 1 opt tol

               if (inform .ne. 0) NeedLU = .true.

               call s5LP
     &            ( inform,
     &              FP, probTag, Elastic,
     &              subOptimize, mnrLog, NeedLU, Needx,
     &              m, n, nb, nDegen, itQP, itnlim, itn,
     &              eMode, lvlObjE, mnrPrtlvl,
     &              minimize, iObjPP, scaleObj, objA,
     &              condZmax, tolOptFP, tolOptQP, tolx,
     &              nInf, sInf, elastics, nInfE, sInfE, wtInf,
     &              piNorm, rgNorm,
     &              neJ, nlocJ, locJ, indJ, Jcol,
     &              eType, eState, feasType, hs, kBS,
     &              bl, bu, blQP, buQP, blBS, buBS,
     &              gBS, pi, rc, nnCon0, nnCon, QPrhs, scales,
     &              x, xBS, x0,
     &              iy, iy1, y, y1,
     &              cw, lencw, iw, leniw, rw, lenrw )

!              Possible inform values are = -3,-2,-1, 0, >0

               if (inform .ne. 0) then
                  if (inform .gt. 0) then
                     iExit = inform ! Fatal error
                  else if (inform .eq. -3) then
                     iExit = 31 ! iterations limit
                  else if (nInf .gt. 0) then
                     iExit = 11 ! infeasible (should not happen here)
                  end if
                  if (iExit .ne. 0) go to 800
               end if

               if (inform .eq. 0  .and.  mnrPrtlvl .ge. 1) then
                  write(str, 8400) itn, nviol
                  call snPRNT( 23, str, iw, leniw )
               end if
            end if

            nInf  = 0
            sInf  = zero
            nInfE = 0
            sInfE = zero

!           x(1:nnH) are feasible and must remain so.

            call iload ( nnH, 0, eType, 1 )

         else if (lvlPPm .eq. 2) then
!           ------------------------------------------------------------
!           Minimize the two-norm of (x-x0).
!           ------------------------------------------------------------
            nObjPP     = 0      ! No explicit gradient in proximal point
            nObjPP0    = 1
            GotR       = .false.
            typeLU     = BT
            HvCalls    = 0

            probTag    = 'norm(x-x0) problem  '
            if (nnH .lt. n) then
               eigH    =  SEMDEF
            else
               eigH    =  POSDEF
            end if
            iw(mnrHdP) = 1         ! Switch to QP print   heading
            iw(mnrHdS) = 1         ! Switch to QP summary heading

            Elastic    = .false.   ! Start in normal mode
            eMode      = 0         ! elastic mode not needed.
            lvlObjE    = 0         ! No objective in elastic mode

            Needx       = .false.
            itQPmax     = itQP + 100 ! Limit the minor iterations
            itQPtarget  = itQP + 100
            mSBsave     = iw(mNewSB)
            iw(mNewSB)  = 100       ! and the number of new superbasics
            subOptimize = NoSubOpt
            tolOptFP    = eps2
            tolOptPP    = 1.0d-2    ! Sloppy phase 2 opt tol

            call s5QP
     &         ( inform,
     &           QPP, probTag, subOptimize,
     &           mnrlog, s8HxPP, s8HxQP, HvCalls, eigH,
     &           Elastic, GotR, NeedLU, typeLU, Needx,
     &           lenR, m, maxS, mBS, n, nb, nDegen,
     &           nnH0, nnH, nObjPP0, nObjPP, nnH0, nnH, nS,
     &           itQP, itQPmax, itQPtarget, itn,
     &           eMode, lvlObjE, mnrPrtlvl,
     &           minimize, iObjPP, scaleObj, objA, objPP,
     &           condZmax, Hcondbnd, tolOptFP, tolOptPP, tolx,
     &           nInf, sInf, elastics, nInfE, sInfE, wtInf,
     &           piNorm, rgNorm,
     &           neJ, nlocJ, locJ, indJ, Jcol,
     &           neH, nlocH, locH, indH, Hcol,
     &           eType, eState, feasType, hs, kBS,
     &           bl, bu, blQP, buQP, blBS, buBS,
     &           gBS, gQP, gQP, Hdx, pBS, pi, R, rc, rg,
     &           nnCon0, nnCon, QPrhs, scales,
     &           nnH0, nnH, x0, x, xBS, x0,
     &           iy, iy1, y, y1, y2,
     &           cw, lencw, iw, leniw, rw, lenrw,
     &           cw, lencw, iw, leniw, rw, lenrw )
            iw(mNewSB) = mSBsave

!           Check for trouble.
!           Possible inform values are = -9(-1)-1, 0, >0

            if (inform .ne. 0  .or.  nInf .gt. 0) then
               if (inform .gt. 0) then
                  iExit = inform ! Fatal LU error or time limit
               else if (inform .eq. -3  .and.  itn  .ge. itnlim) then
                  iExit = 31    ! iterations limit
               else if (inform .eq. -1  .or.   nInf .gt.      0) then
                  iExit = 11    ! infeasible (should not happen here)
               end if
               if (iExit .ne. 0) go to 800
            end if

!           Note: objQP is an updated quantity that may be slightly
!           negative.

            if (mnrPrtlvl .ge. 1) then
               write(str, 8200) itn, lvlPPm, abs(objPP)
               call snPRNT( 31, str, iw, leniw )
               call snPRNT( 22, str, iw, leniw )
            end if
         end if ! Proximal Point method 1
      end if ! nnH > 0

  800 return

 8000 format(' Itn', i7, ': Feasible linear rows')
 8100 format(' Itn', i7, ': PP', i1, '.  Minimizing  Norm(x-x0)')
 8200 format(' Itn', i7,
     &       ': PP', i1, '.  Norm(x-x0) approximately minimized  (',
     &               1p, e8.2, ')')
 8300 format(' Itn', i7,
     &       ': PP1.  Making nonlinear variables feasible')
 8400 format(' Itn', i7, ': PP1. ',
     &               i7, ' nonlinear variables made feasible')

      end ! subroutine s8getFeasLC

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8getWeights
     &   ( job, Boosted, itn, g0Norm, gNorm,
     &     wtInf0, wtInf, wtMax,
     &     weight, wtFactor, wtScale, iw, leniw )

      implicit
     &     none
      logical
     &     Boosted
      integer
     &     job, itn, leniw, iw(leniw)
      double precision
     &     g0Norm, gNorm, wtInf0, wtInf, wtMax,
     &     weight, wtFactor, wtScale

!     ==================================================================
!     s8getWeights  initializes or updates the elastic weight  wtInf.
!     The elastic weight is given by  wtInf = wtScale*weight,
!     where wtScale is some scale-dependent quantity (fObj here).
!     wtInf is increased by redefining weight as weight*wtFactor, where
!     wtFactor is a constant factor.
!
!     weight, wtFactor and wtScale are 'saved' local variables.
!
!     20 Feb 1997: First version of s8getWeights.
!     27 Apr 2001: wtMax introduced as parameter instead of local.
!     03 Aug 2003: snPRNT adopted.
!     03 Aug 2003: Current version of s8getWeights.
!     ==================================================================
      character
     &     str*80
      double precision
     &     newWeight
!     ------------------------------------------------------------------
      double precision   one,           ten
      parameter         (one  = 1.0d+0, ten   = 10.0d+0)
      integer            SetWeight,     BoostWeight
      parameter         (SetWeight = 0, BoostWeight = 1)
!     ------------------------------------------------------------------
      if (job .eq. SetWeight) then

!        Set the weight.
!        weight is the ``unscaled'' weight on the infeasibilities.
!        wtScale is a scale factor based on the current gradient.

         wtScale  = max( one, g0Norm + gNorm)
!        wtScale  = max( 1.0d+2, g0Norm + gNorm)
!        wtScale  = max( 1.0d+0, g0Norm + gNorm)
         wtFactor = ten
         weight   = wtInf0
         wtInf    = wtScale*weight

      else if (job .eq. BoostWeight) then

!        If possible, boost the weight.

         newWeight = min( wtFactor*weight, wtMax )
         Boosted   = newWeight .gt. weight

         if (Boosted) then
            weight   = newWeight
            wtInf    = weight*wtScale
            wtFactor = ten*wtFactor
            write(str, 1000) itn, wtInf
            call snPRNT( 23, str, iw, leniw )
         end if
      end if

      return

 1000 format(' Itn', i7, ': Elastic weight increased to ', 1p, e11.3)

      end ! subroutine s8getWeights

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Gloc
     &   ( nnCon, nnJac, neJ, nlocJ, locJ, indJ, negCon, nlocG, locG )

      implicit
     &     none
      integer
     &     neJ, negCon, nlocG, nlocJ, nnCon, nnJac, indJ(neJ)
      integer
     &     locJ(nlocJ), locG(nlocG)

!     ==================================================================
!     s8Gloc  counts the number of nonlinear Jacobian elements and
!     assembles their column pointers in locG.
!
!     29 Oct 2000: First version of s8Gloc.
!     31 Aug 2008: Local variable used for negCon.
!     ==================================================================
      integer
     &     ir, j, k, neg
!     ------------------------------------------------------------------
      neg     = 0
      locG(1) = 1
      do j = 1, nnJac
         do  k = locJ(j), locJ(j+1)-1
            ir = indJ(k)
            if (ir .gt. nnCon) go to 100
            neg = neg + 1
         end do
  100    locG(j+1) = neg + 1
      end do

      end ! subroutine s8Gloc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Gprod
     &   ( Task, tolz,
     &     neJ, nlocJ, locJ, indJ, negCon, nlocG, locG, gCon,
     &     alpha, x, lenx, beta, y, leny )

      implicit
     &     none
      integer
     &     Task, neJ, negCon, nlocG, nlocJ, lenx, leny, indJ(neJ),
     &     locJ(nlocJ), locG(nlocG)
      double precision
     &     tolz, alpha, beta, gCon(negCon), x(lenx), y(leny)

!     ==================================================================
!     s8Gprod computes matrix-vector products involving J and x.  The
!     variable task specifies the operation to be performed as follows:
!       task = 'N' (normal)          y := alpha*J *x + beta*y,
!       task = 'T' (transpose)       y := alpha*J'*x + beta*y,
!     where alpha and beta are scalars, x and y are vectors and J is a
!     sparse matrix whose columns are in natural order.
!
!     26 Oct 2000: Current version.
!     ==================================================================
      integer
     &     i, ig, iJ, ir, j
      double precision
     &     alphxj, sum, xj
!     ------------------------------------------------------------------
      integer            Normal,        Transp
      parameter         (Normal = 0,    Transp = 1)
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      if (alpha .eq. zero  .and.  beta .eq. one)
     &   return

!     First form  y := beta*y.

      if (beta .ne. one) then
         if (beta .eq. zero) then
            do i = 1, leny
               y(i) = zero
            end do
         else
            do i = 1, leny
               y(i) = beta*y(i)
            end do
         end if
      end if

      if (alpha .eq. zero) then

!        Relax

      else if (alpha .eq. (-one)) then

         if (Task .eq. Normal) then
            do  j = 1, lenx
               xj = x(j)
               if (abs( xj ) .gt. tolz) then
                  ig = locG(j)
                  do iJ = locJ(j), locJ(j+1)-1
                     ir = indJ(iJ)
                     if (ir .gt. leny) go to 100
                     y(ir) = y(ir) - gCon(ig)*xj
                     ig    = ig + 1
                  end do
               end if
  100          continue
            end do

         else if (Task .eq. Transp) then

            do j   = 1, leny
               sum = y(j)
               ig  = locG(j)
               do iJ = locJ(j), locJ(j+1)-1
                  ir = indJ(iJ)
                  if (ir .gt. lenx) go to 200
                  sum = sum - gCon(ig)*x(ir)
                  ig  = ig + 1
               end do
  200          y(j) = sum
            end do
         end if

      else ! General alpha

         if (Task .eq. Normal) then
            do j = 1, lenx
               alphxj = alpha*x(j)
               if (abs( alphxj ) .gt. tolz) then
                  ig  = locG(j)
                  do iJ = locJ(j), locJ(j+1)-1
                     ir = indJ(iJ)
                     if (ir .gt. leny) go to 300
                     y(ir) = y(ir) + gCon(ig)*alphxj
                     ig    = ig + 1
                  end do
               end if
  300          continue
            end do
         else if (Task .eq. Transp) then
            do j   = 1, leny
               sum = zero
               ig  = locG(j)
               do iJ = locJ(j), locJ(j+1)-1
                  ir = indJ(iJ)
                  if (ir .gt. lenx) go to 400
                  sum = sum + gCon(ig)*x(ir)
                  ig  = ig + 1
               end do
  400          y(j) = y(j) + alpha*sum
            end do
         end if ! task .eq. Normal
      end if ! general alpha

      end ! subroutine s8Gprod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Gsize
     &   ( m, nnCon, nnJac, neJ, nlocJ, locJ, indJ, negCon )

      implicit
     &     none
      integer
     &     m, neJ, negCon, nlocJ, nnCon, nnJac, indJ(neJ)
      integer
     &     locJ(nlocJ)

      !=================================================================
      ! s8Gsize  counts the number of nonlinear Jacobian elements.
      !
      ! 04 Nov 2000: First version of s8Gsize
      ! 31 Aug 2008: Local variable used for negCon.
      !=================================================================
      integer
     &     ir, k, last, neg, nlocG
      !-----------------------------------------------------------------
      neg   = 0
      nlocG = nnJac + 1

      if (nnCon .gt. 0) then
         last = locJ(nlocG) - 1
         if (nnCon .eq. m) then
            neg = last
         else
            do  k = 1, last
               ir = indJ(k)
               if (ir .le. nnCon) neg = neg + 1
            end do
         end if
      end if
      negCon = max( 1, neg )

      end ! subroutine s8Gsize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Infs
     &   ( Elastic, n, nb, nnCon0, nnCon, tolx, wtInf,
     &     prInf, dualInf, jprInf, jdualInf, bl, bu, Fx, rc, x )

      implicit
     &     none
      logical
     &     Elastic
      integer
     &     n, nb, nnCon0, nnCon, jprInf, jdualInf
      double precision
     &     tolx, wtInf, dualInf, prInf, bl(nb), bu(nb), rc(nb), x(nb),
     &     Fx(nnCon0)

!     ==================================================================
!     s8Infs computes the maximum primal and dual infeasibilities,
!     using bl, bu, rc, x and the true nonlinear slacks Fxslk.
!     The linear constraints and bounds are assumed to be satisfied.
!     The primal infeasibility is therefore the maximum violation of
!     the nonlinear constraints.
!     The dual infeasibility is the maximum complementarity gap
!     for the bound constraints (with bounds assumed to be no further
!     than 1.0 from each x(j)).
!
!     prInf, dualInf   return the max primal and dual infeas.
!
!     20 Feb 1994: First version based on Minos 5.5 routine m8infs.
!     25 Oct 1996: Elastic mode added.
!     11 Sep 2014: Nonlinear Constraint violations replace slack values.
!     29 May 2015: Elastic variables computed from Fx instead of slacks.
!     ==================================================================
      integer
     &     i, j
      double precision
     &     dj, slack, tol, viol, v, w, xj
!     ------------------------------------------------------------------
      double precision   zero,           one
      parameter        ( zero = 0.0d+0,  one = 1.0d+0 )
!     ------------------------------------------------------------------
      jprInf = 0
      prInf  = zero
      tol    = tolx

!     See how much  Fx  violates the bounds on the nonlinear slacks.
!     prInf is the maximum violation.

      do i = 1, nnCon
         j     = n + i
         slack = Fx(i)
         viol  = max( zero, bl(j) - slack, slack - bu(j) )
         if (prInf .lt. viol) then
            prInf  = viol
            jprInf = j
         end if
      end do

!     ------------------------------------------------------------------
!     + rc(j)  is the multiplier for lower bound constraints.
!     - rc(j)  is the multiplier for upper bound constraints.
!     dualInf is the maximum component-wise complementarity measure.
!     ------------------------------------------------------------------
      jdualInf = 0
      dualInf  = zero
      do  j  = 1, nb
         dj  = rc(j)
         if (dj .ne. zero) then
            if (j .le. n  .or. j .gt. n+nnCon) then
               xj = x(j)
            else
               xj = Fx(j-n)
            end if

            if (     dj .gt. zero) then
               dj =   dj * min(   xj  - bl(j), one )
            else if (dj .lt. zero) then
               dj = - dj * min( bu(j) -   xj , one )
            end if

            if (dualInf .lt. dj) then
               dualInf   =  dj
               jdualInf  =  j
            end if
         end if ! dj nonzero
      end do

!     ------------------------------------------------------------------
!     Include contributions from the elastic variables.
!     ------------------------------------------------------------------
      if (Elastic) then
         do i = 1, nnCon
            j     = n + i
            dj    = rc(j)
            slack = Fx(i)
            v  = bl(j) - slack
            w  = slack - bu(j)

            if      (v .gt. tol) then
               dj = abs(wtInf - dj) * min( v, one )
            else if (w .gt. tol) then
               dj = abs(wtInf + dj) * min( w, one )
            else
               dj = zero
            end if

            if (dualInf .lt. dj) then
               dualInf   =  dj
               jdualInf  =  j
            end if
         end do
      end if

      end ! subroutine s8Infs

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8InitH
     &   ( itn, nnH0, nnH, gNorm0, gNorm, UD0, HD,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     itn, leniw, lenrw, nnH0, nnH, iw(leniw)
      double precision
     &     gNorm0, gNorm, UD0, HD(nnH0), rw(lenrw)

!     ==================================================================
!     s8InitH  defines the initial approximate Hessian H.
!
!     On entry
!     --------
!     gNorm0  is the term from the scaled constant objective row.
!     gNorm   approximates the two-norm of the nonlinear obj. gradient.
!
!     On exit
!     --------
!     UD0     is a positive scalar such that H = UD0*UD0*I
!     HD      is the diagonal marix such that diag(HD) = H.
!
!     19 Aug 2013: First version based on dnInit.
!     20 Aug 2013: Last update.
!     ==================================================================
      integer
     &     HDInfo
!     ------------------------------------------------------------------
      double precision   one
      parameter         (one     = 1.0d+0)

      integer            UnSet
      parameter         (UnSet   =-1)

      double precision   Hweight
      parameter         (Hweight = 1.0d+0)

      double precision   U0max,            U0min
      parameter         (U0max   = 1.0d+1, U0min = 1.0d-2)

      parameter         (HDInfo  = 243) ! Approximate Hessian type
!     ------------------------------------------------------------------

      if (nnH .eq. 0) then
         UD0  = one             ! UD0 is not used
      else
         UD0  = sqrt(Hweight*max(one,(gNorm0 + gNorm)))
      end if

      UD0  = min( max( UD0, U0min ), U0max )

      if (nnH .gt. 0) then
         iw(HDInfo) = UnSet
         call s8ResetH
     &      ( itn, nnH, UD0, HD, iw, leniw, rw, lenrw )
      end if

      end ! subroutine s8InitH

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8InitPen
     &   ( nnCon, penParm, xPen0, xPen, rw, lenrw )

      implicit
     &     none
      integer
     &     lenrw, nnCon
      double precision
     &     penParm(4), xPen0, xPen(nnCon), rw(lenrw)

!     ==================================================================
!     s8InitPen  defines initial values for the penalty parameters.
!
!     19 Aug 2013: First version based on dnInit.
!     20 Aug 2013: Last update.
!     ==================================================================
      double precision
     &     eps
!     ------------------------------------------------------------------
      integer            penDamp,   penMax,   penNorm,   IncRun
      parameter         (penDamp=1, penMax=2, penNorm=3, IncRun=4)

      double precision   one
      parameter         (one   = 1.0d+0)
!     ------------------------------------------------------------------
      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16

!     ---------------------------------------------
!     Initialize the penalty parameters.
!     ---------------------------------------------
      penParm(IncRun ) = one
      penParm(penDamp) = one
      penParm(penMax ) = one / eps
      penParm(penNorm) = xPen0

      call dload ( nnCon, xPen0, xPen, 1 )

      end ! subroutine s8InitPen

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Merit
     &   ( nnCon, fMerit, gMerit, HMerit, penParm,
     &     Fv, xPen, w, rw, lenrw )

      implicit
     &     none
      integer
     &     nnCon, lenrw
      double precision
     &     fMerit, gMerit, HMerit, penParm(4),
     &     Fv(nnCon), xPen(nnCon), w(nnCon), rw(lenrw)

!     ==================================================================
!     s8Merit  computes the contributions to the merit function and its
!     directional derivative from the nonlinear constraints.
!     The penalty parameters  xPen(j)  are increased if
!     the directional derivative is not sufficiently negative.
!
!     On entry:
!         Fv     is the violation c(x) + A(linear)x - s,  where
!                s  minimizes the merit function with respect to the
!                nonlinear slacks only.
!
!     30 Dec 1991: First version based on Npsol 4.0 routine npmrt.
!     02 Nov 1996: Multipliers no longer updated here.
!     19 Jul 1997: Thread-safe version.
!     21 Oct 2000: Made compatible with SNOPT 6.1
!     21 Oct 2000: Current version of s8Merit.
!     ==================================================================
      external
     &     ddiv, ddot, dnrm2
      logical
     &     Boost, OverFlow
      integer
     &     i
      double precision
     &     ddiv, ddot, dnrm2, eps0, incRun, penDamp, penMax, penNorm,
     &     ppscl, penlty, penMin, penNew,
     &     penOld, rtUndf, xPen0, xPeni, wnorm
!     ------------------------------------------------------------------
      double precision   zero,          half,          two
      parameter         (zero = 0.0d+0, half = 0.5d+0, two = 2.0d+0)
!     ------------------------------------------------------------------
      eps0     = rw(  2)
      rtUndf   = rw( 10)
      xPen0    = rw( 89)

      penDamp  = penParm(1)
      penMax   = penParm(2)
      penNorm  = penParm(3)
      incRun   = penParm(4)

      OverFlow = .false.

!     Find the quantities that define  penMin, the vector of minimum
!     two-norm such that the directional derivative is one half of
!     approximate curvature   - (p)'H(p).
!     The factor  rtUndf  tends to keep  xPen  sparse.

      do i = 1, nnCon
         if (abs( Fv(i) ) .le. rtUndf) then
            w(i) = zero
         else
            w(i) = Fv(i)**2
         end if
      end do

      wnorm  = dnrm2 ( nnCon, w, 1 )
      ppscl  = ddiv  ( gMerit + half*HMerit, wnorm, OverFlow )
      if (abs(ppscl) .le. penMax  .and.  .not. OverFlow) then
!        ---------------------------------------------------------------
!        Bounded  penMin  found.  The final value of  xPen(i)  will
!        never be less than  penMin(i).  A trial value  penNew  is
!        computed that is equal to the geometric mean of the previous
!        xPen  and a damped value of penMin.  The new  xPen  is defined
!        as  penNew  if it is less than half the previous  xPen  and
!        greater than  penMin.
!        ---------------------------------------------------------------
         do i = 1, nnCon
            penMin = max( (w(i)/wnorm)*ppscl, zero )
            xPeni  = xPen(i)

            penNew = sqrt( xPeni*(penDamp + penMin) )
            if (penNew .lt. half*xPeni) xPeni = penNew
            xPeni   = max (xPeni, penMin)
            xPen(i) = max (xPeni, xPen0 )
         end do

         penOld  = penNorm
         penNorm = dnrm2( nnCon, xPen, 1 )

!        ---------------------------------------------------------------
!        If  IncRun = true,  there has been a run of iterations in
!        which the norm of  xPen  has not decreased.  Conversely,
!        IncRun = false  implies that there has been a run of
!        iterations in which the norm of xPen has not increased.  If
!        IncRun changes during this iteration the damping parameter
!        penDamp is increased by a factor of two.  This ensures that
!        xPen(j) will oscillate only a finite number of times.
!        ---------------------------------------------------------------
         Boost  = .false.
         if (incRun .gt. zero .and. penNorm .lt. penOld) then
            Boost = .true.
         end if
         if (incRun .lt. zero .and. penNorm .gt. penOld) then
            Boost = .true.
         end if

         if (Boost) then
            penDamp = min( 1/eps0, two*penDamp )
            incRun  = -incRun
         end if
      end if

!     ------------------------------------------------------------------
!     Compute the new value and directional derivative of the
!     merit function.
!     ------------------------------------------------------------------
      call dcopy ( nnCon, Fv  , 1, w, 1 )
      call ddscl ( nnCon, xPen, 1, w, 1 )

      penlty = ddot  ( nnCon, w, 1, Fv, 1 )
      fMerit = fMerit  + half*penlty
      gMerit = gMerit  -      penlty

      penParm(1) = penDamp
      penParm(2) = penMax
      penParm(3) = penNorm
      penParm(4) = incRun

      end ! subroutine  s8Merit

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8optimizeSlacks
     &   ( n, nb, nnCon, nInfE, sInfE, featol, wtInf,
     &     bl, bu, Fv, x, yCon, xPen, Fx )

      implicit
     &     none
      integer
     &     n, nb, nnCon, nInfE
      double precision
     &     featol, sInfE, wtInf, bl(nb), bu(nb),
     &     x(nb), xPen(nnCon), Fx(nnCon), Fv(nnCon), yCon(nnCon)

!     ==================================================================
!     s8optimizeSlacks computes the nonlinear constraint violations:
!        Fv = fCon + A(linear)*x - (optimal nonlinear slacks)
!           =        Fx          - (optimal nonlinear slacks)
!
!     The nonlinear slacks are adjusted so that they minimize the merit
!     function with  x  and  yCon  held constant.
!
!     On entry,
!        bl     =  the array of lower bounds.
!        bu     =  the array of upper bounds.
!        x      =  the current x.
!        Fx     =  fCon + A(linear)*x,   defined in s8Fx.
!
!     On exit,
!        x      =  contains the optimal slacks in x(n+1:n+nnCon).
!                  Other elements are left unchanged.
!        Fv     =  fCon + A(linear)*x - optimal slacks.
!        Fx     =  unaltered.
!        nInfE  =  number of optimal slack infeasibilities.
!        sInfE  =  sum of optimal slack infeasibilities.
!
!     19 Nov 2012: First version based on SNOPT routine s8sOpt.
!     21 Sep 2013: Added the sum and number of slack infeasibilities.
!     01 Jan 2014: Reworked to handle elastic mode better.
!     20 Jan 2014: Reworked and simplified to compute the correct sInfE.
!     30 Apr 2015: Major revision guarantees slacks are bounded.
!     24 May 2015: Piecewise differentiability implemented correctly.
!     ==================================================================
      integer
     &     i, j
      double precision
     &     blj, buj, ci, dM, yi, vL, vU, sj, peni
!     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
!     ------------------------------------------------------------------
      nInfE = 0
      sInfE = zero

      do i = 1, nnCon
         ci    =   Fx(i)
         peni  = xPen(i)
         yi    = yCon(i)

         j     = n + i
         sj    =  x(j)
         blj   = bl(j)
         buj   = bu(j)

         if (sj .ge. blj - featol  .and.  sj .le. buj + featol ) then
!           ------------------------------------------------------------
!           Feasible slack.
!           Reset with the optimal slack in  [blj, buj].
!           ------------------------------------------------------------
            dM    = yi -  peni*(ci - sj)

            if      (dM .eq. zero) then
!              Relax
            else if (dM .lt. zero) then

!              s must increase to reduce the merit function.
!
!              sMin = c - y/xPen
!              Mul >  0 and dM <  0   =>  s0 <= sMin <= c, Pen ne 0.
!              Mul <= 0 and c  >= s0  =>  s0 <= s    <= c.

               if (yi .gt. zero) then
                  sj = min(buj, ci - yi/peni)
               else if (ci .ge. sj) then
                  sj = min(buj, ci)
               end if
            else if (dM .gt. zero) then

!              s must decrease to reduce the merit function.
!
!              sMin = c - y/Pen
!              Mul <  0 and dM >   0  =>  c  <  sMin <= s0, Pen ne 0.
!              Mul >= 0 and c  <= s0  =>  s0 <= s    <= c.

               if (yi .lt. zero) then
                  sj = max( blj, ci - yi/peni )
               else if (ci .le. sj) then
                  sj = max( blj, ci )
               end if
            end if

         else
!           ------------------------------------------------------------
!           Infeasible slack in elastic mode.
!           ------------------------------------------------------------
            if (sj .lt. blj - featol) then

!              Reset with the optimal slack in  (-inf, blj].

               yi    = yi - wtInf
               dM    = yi - peni*(ci - sj)

               if (dM .eq. zero) then
!                 Relax
               else if (dM .lt. zero) then

!                 s must increase to reduce the merit function.

                  if (yi .gt. zero) then
                     sj = min(ci - yi/peni, blj)
                  else if (ci .ge. sj) then
                     sj = min(ci,           blj)
                  end if
               else if (dM .gt. zero) then

!                 s must decrease to reduce the merit function.

                  if (yi .lt. zero) then
                     sj = ci - yi/peni
                  else if (ci .le. sj) then
                     sj = ci
                  end if
               end if

            else !if (sj .gt. buj + featol) then

!              Reset with the optimal slack in  [buj, +inf)

               yi  = yi + wtInf
               dM  = yi - peni*(ci - sj)

               if      (dM .eq. zero) then
!                 Relax
               else if (dM .lt. zero) then

!                 s must increase to reduce the merit function.

                  if (yi .gt. zero) then
                     sj = ci - yi/peni
                  else if (ci .ge. sj) then
                     sj = ci
                  end if
               else if (dM .gt. zero) then

!                 s must decrease to reduce the merit function.

                  if (yi .lt. zero) then
                     sj = max(buj, ci - yi/peni)
                  else if (ci .le. sj) then
                     sj = max(buj, ci)
                  end if
               end if
            end if

            vL    = blj -  sj
            vU    = sj  - buj
            if (vL .gt. featol  .or.  vU .gt. featol) then
               nInfE = nInfE + 1
               sInfE = sInfE + max (vL, vU )
            end if
         end if

         Fv(i) = ci - sj
         x(j)  = sj

      end do

      end ! subroutine s8optimizeSlacks

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8HxLP
     &   ( nnH, x, Hx, status, cu, lencu, iu, leniu, ru, lenru)

      implicit
     &     none
      integer
     &     nnH, status, lencu, leniu, lenru, iu(leniu)
      double precision
     &     x(nnH), Hx(nnH), ru(lenru)
      character
     &     cu(lencu)*8

!     ==================================================================
!     s8HxLP is the argument qpHx for s5solve when s5solve is called
!     from one of the snOpt wrappers.
!
!     04 Dec 2004: First version of s8HxLP.
!     04 Dec 2004: Current version of s8HxLP.
!     ==================================================================

      ! Relax

      end ! subroutine s8HxLP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8HxNP
     &   ( userHv, nnH,
     &     neH, nlocH, locH, indH, Hcol,
     &     v, Hv, status,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     userHv
      integer
     &     lencu, leniu, lenru, lencw, leniw, lenrw, neH, nlocH,
     &     nnH, status, locH(nlocH), indH(neH), iu(leniu), iw(leniw)
      double precision
     &     Hcol(neH), Hv(nnH), v(nnH), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

      !=================================================================
      ! s8HxNP (aka Hprod)  is called by the QP solver.
      ! It is a wrapper that calls userHv  (aka Hprod1) to compute Hv,
      ! the product of the Hessian and v and scales it.
      !
      ! 30 Apr 2011: First version of s8npHx
      ! 07 May 2011: HvMode added to the argument list of userHv.
      !=================================================================
      logical
     &     scaled
      integer
     &     HvMode, lvlScale, lscales, lx, lxScaled, lyCon, nnCon
      !-----------------------------------------------------------------
      nnCon     = iw( 23) ! # of nonlinear constraints
      lvlScale  = iw( 75) ! scale option

      lx        = iw(300)       ! x(nb) at which H is defined
      lscales   = iw(296) ! scales(nb)  = row and column scales
      lxScaled  = iw(302) ! xScaled(n)  = copy of the scaled x
      lyCon     = iw(348) ! yCon (nnCon) = multipliers for fCon

      scaled    = lvlScale .gt. 0

      ! Determine the call-status for the NP subproblem.

      call s8callStatus( status, iw, leniw )

      if ( scaled ) then
         call dcopy ( nnH, v          , 1, rw(lxScaled), 1 )
         call ddscl ( nnH, rw(lscales), 1, v           , 1 )
         ! Scale the base x and the multipliers
      end if

      HvMode = 1                ! Do NOT compute a new Hessian.
      call userHv
     &   ( HvMode, nnCon, nnH,
     &     neH, nlocH, locH, indH, Hcol,
     &     rw(lyCon), rw(lx),
     &     v, Hv, status,
     &     cu, lencu, iu, leniu, ru, lenru )

      if (scaled) then
         call dcopy ( nnH, rw(lxScaled), 1, v , 1 )
         call ddscl ( nnH, rw(lscales) , 1, Hv, 1 )
         ! Unscale the base x and the multipliers!
      end if

      end ! subroutine s8HxNP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8HxPP
     &   ( Hprod, nnH,
     &     neH, nlocH, locH, indH, Hcol,
     &     x, Hx, status,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod
      integer
     &     lencu, leniu, lenru, lencw, leniw, lenrw, neH, nnH, nlocH,
     &     indH(neH), locH(nlocH), status, iu(leniu), iw(leniw)
      double precision
     &     Hcol(neH), Hx(nnH), x(nnH), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s8HxPP  defines the product  H*x  for the proximal-point QP
!     subproblem of snopt.
!
!     On exit,    Hx   = x.
!
!     23 Oct 1993: First version of s8HxPP.
!     02 Aug 2000: Current version.
!     ==================================================================
      call dcopy ( nnH, x, 1, Hx, 1 )

      end ! subroutine s8HxPP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8HxQP
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
!     s8HxQP  computes the user-defined product  Hx  and scales it.
!
!     07 Nov 2014: First   version of s8Hxqp.
!     ==================================================================
      logical
     &     Scaled
      integer
     &     lvlScale, lscales, lxScaled, status
!     ------------------------------------------------------------------
      lvlScale  = iw( 75) ! scale option
      lscales   = iw(296) ! scales(nb)  = row and column scales
      lxScaled  = iw(302) ! xScaled(n)  = copy of scaled x(nnL)

      Scaled    = lvlScale .gt. 0

      ! Determine the user-function call-status.

      call s5CallStatus( status, iw, leniw )

      if (Scaled) then
         call dcopy ( nnH,           x, 1, rw(lxScaled), 1 )
         call ddscl ( nnH, rw(lscales), 1,            x, 1 )
      end if

      call usrHx
     &   ( nnH, x, Hx, status,
     &     cu, lencu, iu, leniu, ru, lenru )

      if (Scaled) then
         call dcopy ( nnH, rw(lxScaled), 1, x , 1 )
         call ddscl ( nnH, rw(lscales) , 1, Hx, 1 )
      end if

      end ! subroutine s8HxQP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8HxNull
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
!     It should never be called.
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
     &     '     has been called from SNOPT. ')

      end ! subroutine s8HxNull

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8rand( leng, neg, g )

      implicit
     &     none
      integer
     &     leng, neg
      double precision
     &     g(leng)

!     ==================================================================
!     s8rand  fills the array g with random numbers.
!
!     15 Nov 1991: First version of s8rand in s8aux.
!     30 Jun 1999: Current version.
!     ==================================================================
      integer
     &     seeds(3)
!     ------------------------------------------------------------------
      if (neg .le. 0) return

      seeds(1) = 1547
      seeds(2) = 2671
      seeds(3) = 3770

      call ddrand( neg, g, 1, seeds )

      end ! subroutine s8rand

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8rc
     &   ( scaleObj, minimize, iObj, m, n, nb,
     &     nnObj0, nnObj, nnCon, nnJac, negCon,
     &     neJ, nlocJ, locJ, indJ, Jcol,
     &     gObj, gCon, pi, rc )

      implicit
     &     none
      integer
     &     minimize, iObj, m, n, nb, nnObj0, nnObj, nnCon, nnJac,
     &     negCon, neJ, nlocJ, indJ(neJ), locJ(nlocJ)
      double precision
     &     scaleObj, Jcol(neJ), gObj(nnObj0), gCon(negCon),
     &     pi(m), rc(nb)

!     ==================================================================
!     s8rc   computes reduced costs rc = gObj - ( A  -I )'*pi,
!     using  gCon  as the top left-hand corner of A.
!     gCon, gObj and pi are assumed to exist.
!
!     s8rc   is called by s8SQP.
!
!     28 Sep 1993: First version, derived from m4rc.
!     31 Oct 1996: Min sum option added.
!     30 Oct 2000: Current version of s8rc.
!     ==================================================================
      integer
     &     ir, j, k, l
      double precision
     &     dj, signObj
!     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )
!     ------------------------------------------------------------------
      l     = 0

      do j  = 1, nnJac
         dj = zero
         do k  = locJ(j), locJ(j+1) - 1
            ir = indJ(k)
            if (ir .le. nnCon) then
               l  = l  + 1
               dj = dj + pi(ir)*gCon(l)
            else
               dj = dj + pi(ir)*Jcol(k)
            end if
         end do
         rc(j) = -dj
      end do

      do j  = nnJac+1, n
         dj = zero
         do k  = locJ(j), locJ(j+1) - 1
            ir = indJ(k)
            dj = dj  +  pi(ir) * Jcol(k)
         end do
         rc(j) = -dj
      end do

      call dcopy ( m, pi, 1, rc(n+1), 1 )

!     Include the nonlinear objective gradient.

      signObj = minimize
      if (nnObj .gt. 0) then
         call daxpy ( nnObj, signObj, gObj, 1, rc, 1 )
      end if

      if (iObj .gt. 0) rc(n+iObj) =  rc(n+iObj) + signObj*scaleObj

      end ! subroutine s8rc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8scaleG
     &   ( nnObj, scales, gObj, rw, lenrw )

      implicit
     &     none
      integer
     &     nnObj, lenrw
      double precision
     &     scales(nnObj), gObj(nnObj), rw(lenrw)

!     ==================================================================
!     s8scaleG  scales the objective gradient.
!     s8scaleG is called by funwrapper only if modefg = 2.
!     Hence, it is used to scale known gradient elements (if any),
!     but is not called when missing gradients are being estimated
!     by s6dobj.
!
!     17 Feb 1992: First version.
!     16 Jul 1997: Thread-safe version.
!     02 Jan 2001: Current version of s8scaleG.
!     ==================================================================
      integer
     &     j
      double precision
     &     gdummy, grad
!     ------------------------------------------------------------------
      gdummy = rw( 69) ! definition of 'unset' value

      do j = 1, nnObj
         grad = gObj(j)
         if (grad .ne. gdummy) gObj(j) = grad*scales(j)
      end do

      end ! subroutine s8scaleG

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8scaleJ
     &   ( nnCon, nnJac, negCon, n, scales,
     &     neJ, nlocJ, locJ, indJ, gCon, rw, lenrw )

      implicit
     &     none
      integer
     &     n, neJ, negCon, nnCon, nnJac, nlocJ, lenrw, indJ(neJ),
     &     locJ(nlocJ)
      double precision
     &     scales(n+nnCon), gCon(negCon), rw(lenrw)

!     ==================================================================
!     s8scaleJ  scales the Jacobian.
!     s8scaleJ is called by funwrapper only if modefg = 2.
!     Hence, it is used to scale known gradient elements (if any),
!     but is not called when missing gradients are being estimated
!     by s6dcon.
!
!     17 Feb 1992: First version based on Minos routine m8sclj.
!     16 Jul 1997: Thread-safe version.
!     02 Dec 2001: Current version of s8scaleJ.
!     ==================================================================
      integer
     &    ir, j, k, l
      double precision
     &     Cscale, gdummy, grad
!     ------------------------------------------------------------------
      gdummy = rw( 69) ! definition of 'unset' value

      l    = 0
      do j = 1, nnJac
         Cscale = scales(j)

         do k = locJ(j), locJ(j+1)-1
            ir     = indJ(k)
            if (ir .gt. nnCon) go to 300
            l      = l + 1
            grad   = gCon(l)
            if (grad .ne. gdummy)
     &         gCon(l) = grad*cscale/scales(n+ir)
         end do
  300    continue
      end do

      end ! subroutine s8scaleJ

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8sInf
     &   ( n, nb, nnCon, tolx, nInf, sInf, bl, bu, x )

      implicit
     &     none
      integer
     &      n, nb, nnCon, nInf
      double precision
     &     tolx, sInf, bl(nb), bu(nb), x(nb)

!     ==================================================================
!     s8sInf computes the sum of infeasibilities of the nonlinear slacks
!     using bl, bu and x.
!
!     10 Jan 1997: First version of s8sInf.
!     30 Oct 2000: Current version.
!     ==================================================================
      integer
     &     i, j
       double precision
     &     slack, tol, violL, violU
!     ------------------------------------------------------------------
      double precision   zero
      parameter        ( zero = 0.0d+0 )
!     ------------------------------------------------------------------
      nInf   = 0
      sInf   = zero
      tol    = tolx

!     See how much  x(n+1:n+nnCon) violates its bounds.

      do i = 1, nnCon
         j     = n + i
         slack = x(j)
         violL = bl(j) - slack
         violU = slack - bu(j)
         if (violL .gt. tol  .or.  violU .gt. tol) then
            nInf = nInf + 1
            sInf = sInf + max (violL, violU )
         end if
      end do

      end ! subroutine s8sInf

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8solveQP
     &   ( iExit,
     &     mnrlog, Hprod, Hprod1, HvCalls,
     &     Elastic, GotR,
     &     itn, itQP, lenR, m, maxS, mBS, n, nb,
     &     nnCon0, nnCon, nnObj0, nnObj,
     &     nnH0, nnH, nS, nDegen,
     &     mjrPrtlvl, mnrPrtlvl, minimize,
     &     iObj, scaleObj, objAdd, objQP,
     &     condZHZ, tolOptFP, tolOptQPk, tolx,
     &     nInf, sInf, elastics, nInfE, sInfE, wtInf,
     &     UD0, piNorm,
     &     neJ, nlocJ, locJ, indJ, Jcol,
     &     neH, nlocH, locH, indH, Hcol,
     &     eType, eState, feasType, hs, kBS,
     &     bl, bu, blQP, buQP, blBS, buBS,
     &     gBS, gQP, gObj, HD, Hdx,
     &     pBS, pi, R, rc, rg, rg2, QPrhs, scales,
     &     x, xBS, xQP0, xQP,
     &     iy, iy1, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod, Hprod1, mnrlog
      logical
     &     Elastic, GotR
      integer
     &     elastics, HvCalls, iExit, iObj, itn, itQP, lenR,
     &     lencu, leniu, lenru, lencw, leniw, lenrw, m, maxS, mBS,
     &     mjrPrtlvl, mnrPrtlvl, minimize, n, nb, nDegen, neJ, neH,
     &     nlocJ, nlocH, nnCon0, nnCon, nnObj0, nnObj, nnH0, nnH,
     &     nInf, nInfE, nS, locJ(nlocJ), locH(nlocH),
     &     indJ(neJ), indH(neH), eType(nb), eState(nb), hs(nb),
     &     feasType(mBS), kBS(mBS), iy(nb), iy1(nb),
     &     iu(leniu), iw(leniw)
      double precision
     &     condZHZ, objAdd, objQP, piNorm, scaleObj, sInf, sInfE,
     &     tolOptFP, tolOptQPk, tolx, UD0, wtInf, bl(nb), bu(nb),
     &     blQP(nb), buQP(nb), blBS(mBS), buBS(mBS), gBS(mBS),
     &     gQP(nnH0), gObj(nnObj0), HD(nnH0), Hdx(nnH0), Hcol(neH),
     &     Jcol(neJ), pBS(mBS), pi(m),
     &     QPrhs(nnCon0), R(lenR), rc(nb), rg(maxS), rg2(maxS),
     &     scales(nb), x(nb), xBS(mBS), xQP0(nb), xQP(nb),
     &     y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s8solveQP  computes  xQP, the solution of the QP subproblem.
!     By construction, the problem has  nnH  nonlinear variables,
!
!     The SQP base point  x  is not altered.
!
!     On entry:
!     ---------
!     eType    contains the allowed type of elastic violation.
!     eState   contains the elastic states.

!     The LU factorization is assumed to be known.
!     The arrays  xBS, blBS and buBS are defined.
!
!     On exit:
!     --------
!     xQP0    contains the first QPfeasible point for the nonelastics.
!     xQP     contains the QP solution.
!
!     iExit     Status
!     -----     ------
!      >0         Fatal error
!       0         QP solution found
!      -1         Too many iterations
!      -2         Too many superbasics
!
!     LUrequest =  1  Frequency
!     LUrequest =  2  LU nonzeros increased
!     LUrequest =  3
!     LUrequest =  4
!     LUrequest =  5  Singular after LU mod
!     LUrequest =  6  Unstable LU mod (growth in new column of U)
!     LUrequest =  7  Not enough memory
!     LUrequest =  8
!     LUrequest =  9
!     LUrequest = 10  Row error in setx
!     LUrequest = 11  Big  dx   in setx
!
!     LUrequest = 20  Infeasible nonelastics in phase 2.
!     LUrequest = 21  Iterative refinement failed in QP
!     LUrequest = 22  Unbounded QP
!     LUrequest = 23
!     LUrequest = 24  Small directional derivative in QP
!     LUrequest = 25  Ill-conditioned Z
!     LUrequest = 26  Indefinite Z'HZ in QP
!     LUrequest = 27  R singular after bound swap in QP
!
!     On output,
!     QPerr points to ' ', 't', 'u' or 'w'.
!     QPfea points to ' '  or 'i'.
!
!     30 Dec 1991: First version of s8solveQP (as s8iQP).
!     19 Jul 1997: Thread-safe version.
!     31 Jul 2003: dnormj used for norm of the nonlinear pis.
!     03 Aug 2003: snEXIT and snPRNT  adopted.
!     19 Jun 2008: Hprod, Hprod1 added as arguments.
!     07 Oct 2014: infoTags added.
!     18 Oct 2014: Switch to elastic mode based on piNorm.
!     03 Nov 2014: neH, indH, locH, Hcol added as arguments.
!     29 Dec 2014: Nonzero nnH argument for FP phase.
!     29 Dec 2014: Merged with s8iQN.
!     07 May 2015: Check nS le maxR when limiting maxS = maxR.
!     12 May 2015: condZ estimate checked in s5ZHZ.
!     ==================================================================
      character
     &     probTag*20, str*80
      external
     &     dnormj
      logical
     &     Feasible, LUok, NeedLU, Needx, NonlinearCon,
     &     NormalIn, ReSolve, Solved, Terminate
      integer
     &     eigH, inform, itnlim, itQPmax, itQPtarget, eMode,
     &     linesL, linesS, LUrequest, LUitn, lvlObjE, lvlObjFP, lvlPre,
     &     maxR, maxSB, mMinor, mnrHdP, mnrHdS, ngQP0, ngQP,
     &     preCon, QPmode, QPSolver, subOptimize, typeLU
      integer
     &     QNInfo, MdInfo, LSInfo, FPInfo, QPInfo, FDInfo, HDInfo
      double precision
     &     condZmax0, condZmax, dnormj, flmax, condZHZbnd,
     &     objFP, piNormNLN, plInfy, rgNorm, condZHZmax
!     ------------------------------------------------------------------
      integer            BT
      parameter         (BT       = 3)
      integer            FPS,          QPS
      parameter         (FPS      = 4, QPS    = 5)
      integer            SEMDEF,       POSDEF
      parameter         (SEMDEF   = 0, POSDEF = 1)
      integer            CHOL
      parameter         (CHOL     = 0)
      integer            NO
      parameter         (NO       = 0)
      integer            NoSubOpt,     SubOpt
      parameter         (NoSubOpt =-1, SubOpt = 0)
      integer            Unit
      parameter         (Unit     = 2)
      integer            QPChol,       CG,         QN
      parameter         (QPChol   = 0, CG     = 1, QN = 2)

      parameter         (QPmode = 208) ! Current QP solver
      parameter         (PreCon = 209) ! Current precon mode
      parameter         (LUitn  = 215) ! itns since last factorize
      parameter         (linesL = 220) ! # lines in log     file
      parameter         (linesS = 221) ! # lines in summary file
      parameter         (mnrHdP = 223) ! >0 => Minor heading for iPrint
      parameter         (mnrHdS = 225) ! >0 => Minor heading for iSumm

      parameter         (QNInfo = 237) ! TagInfo(1): QN update type
      parameter         (MdInfo = 238) ! TagInfo(2):
      parameter         (LSInfo = 239) ! TagInfo(3): Line search outcome
      parameter         (FPInfo = 240) ! TagInfo(4): QP Feasibility
      parameter         (QPInfo = 241) ! TagInfo(5)  QP Optimality
      parameter         (FDInfo = 242) ! TagInfo(6)
      parameter         (HDInfo = 243) ! TagInfo(7): Approx Hessian type

      double precision   zero
      parameter         (zero       = 0.0d+0)
!     ------------------------------------------------------------------
      maxR       = iw( 52) ! max columns of R
      QPsolver   = iw( 55) ! = 0:1:2   => QPChol:CG:QN QP solver
      lvlPre     = iw( 77) ! >0     => QN preconditioned CG
      itnlim     = iw( 89) ! limit on total iterations
      mMinor     = iw( 91) ! limit on minor iterations

      flmax      = rw(  8) ! est. of the largest pos. real
      condZHZbnd = rw( 85) ! bound on the condition of ZHZ
      condZmax0  = rw( 86) ! bound on the condition of Z

      iExit      = 0
      itQP       = 0
      probTag    = 'QP subproblem'

      plInfy     = flmax

      condZHZmax = condZHZbnd   ! condH > condZHZmax  => Bad cond Z'HZ
      condZmax   = condZmax0    ! condZ > condZmax    => Bad cond Z
      condZHZ    = zero
      NormalIn   = .not. Elastic

      ! eigH    encodes the inertia of the Hessian H
      !         eigH        Hessian
      !         -----   ---------------------
      !          -1                indefinite
      !           0     positive semidefinite
      !           1     positive semidefinite

      if (nnH .lt. n) then
         eigH    = SEMDEF
      else
         eigH    = POSDEF
      end if

      iw(linesL) = 0
      iw(linesS) = 0
      iw(mnrHdP) = 1
      iw(mnrHdS) = 1

      iw(QPInfo) = 0              ! Optimal  QP subproblem
      iw(FPInfo) = 0              ! Feasible QP subproblem

!     If nS dips below maxR, go back to using the default solver.

      if (nS .le. maxR) then
         if (QPsolver .eq. QPChol  .or.  QPsolver .eq. QN) then
            iw(QPmode) = QPsolver
         end if
      end if

      iw(PreCon) = lvlPre         ! Current precon mode

      ngQP       = nnH
      ngQP0      = max( ngQP, 1 )

      typeLU     = BT
      LUrequest  = 0              ! The first LU has been computed.

      NonlinearCon = nnCon .gt. 0
      Feasible     = nnCon .eq. 0 ! Status of the non-elastic variables

      !-----------------------------------------------------------------
      ! Find a feasible point.
      ! If the constraints are linear, x is already feasible.
      ! Otherwise, find a feasible x for this linearization.
      ! Minimize the sum of the elastic variables
      ! subject to keeping the non-elastic variables feasible.
      ! Elastic variables can move outside their bounds.
      !-----------------------------------------------------------------
      ! Loop back here in the unlikely event that the nonelastic bounds
      ! become infeasible while solving the QP subproblem.

  100 if (.not. Feasible) then

         itQPmax     = itnlim
         itQPtarget  = itnlim
         subOptimize = NoSubOpt

         if (NormalIn) then
            ! Set eMode to switch to Elastic mode on infeasibility.
            eMode = 1
         else
            ! Already in elastic mode.
            eMode = 2
         end if

         ! When getting feasible in elastic mode, set lvlObjFP to use
         ! the composite objective:
         !  w1*Obj + w2*sInf,  with w1 = 0, w2 = wtInf.
         ! This minimizes the sum of the infeasibilities of the
         ! elastic variables subject to the nonelastic constraints.

         lvlObjFP   = 2

         if (iw(QPmode) .eq. QPChol) then
            GotR    = .false.
         end if

         LUok       = .true.
         Terminate  = .false.

!        ===============================================================
       ! while (.not. Terminate  .and.  LUok) do
  500    if    (.not. Terminate  .and.  LUok) then

            NeedLU = LUrequest .gt. 0
            Needx  = NeedLU

            call s5QN
     &         ( inform,
     &           FPS, probTag, subOptimize,
     &           mnrlog, Hprod, Hprod1, HvCalls,
     &           Elastic, GotR, NeedLU, typeLU, Needx,
     &           lenR, m, maxS, mBS, n, nb, nDegen,
     &           ngQP0, ngQP, nnObj0, nnObj, nnH0, nnH, nS,
     &           itQP, itQPmax, itQPtarget, itn,
     &           eMode, lvlObjFP, mnrPrtlvl,
     &           minimize, iObj, scaleObj, objAdd, objFP,
     &           condZHZ, condZmax, tolOptFP, tolOptQPk, tolx,
     &           nInf, sInf, elastics, nInfE, sInfE, wtInf, piNorm,
     &           neJ, nlocJ, locJ, indJ, Jcol,
     &           neH, nlocH, locH, indH, Hcol,
     &           eType, eState, feasType, hs, kBS,
     &           bl, bu, blQP, buQP, blBS, buBS,
     &           gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg, rg2,
     &           nnCon0, nnCon, QPrhs, scales,
     &           nnH0, nnH, x, xQP, xBS, x,
     &           iy, iy1, y, y1, y2,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

            ! Check for trouble.  Here are the possibilities:
            !
            ! inform      Result
            ! ------      ------
            !  >0         Fatal LU error
            !   0         Found a feasible point for the nonelastics
            !  -1         The nonelastics are infeasible
            !  -2         Phase 1 is unbounded
            !  -3         Too many iterations
            !  -4         Void
            !  -5         Superbasic limit exceeded

            if (inform .gt. 0) then
               iExit = inform   ! Fatal LU error or time limit
               go to 900
            end if

            Terminate =       inform .eq.  0  .or.  inform .eq. -3
     &                  .or.  inform .eq. -5

            if (.not. Terminate) then
               !========================================================
               ! Trouble.
               ! inform = -2 implies that phase 1 was unbounded,
               !             which can only occur if a bad basis gives
               !             a large search direction
               ! inform = -1 implies the nonelastics are infeasible,
               !             which should not happen since we already
               !             know a feasible point for the nonelastics.
               !========================================================
               ! Treat both cases as infeasible. Repeatedly refactorize
               ! with tighter tols before declaring LC infeasibility.

               inform = -1
               write(str, 1500) itn
               call snPRNT( 23, str, iw, leniw )
               call s2tryLU
     &            ( itn, 22, nS, LUrequest, LUok, typeLU,
     &              iw, leniw, rw, lenrw )
               if (.not. LUok) then
                  iExit = 15    ! infeasible linear constraints in QP.
                  go to 900
               end if
            end if
            go to 500
         end if
       ! end while
!        ---------------------------------------------------------------

         if (inform .lt. 0) go to 800 ! Itns or time limit exceeded

         if (Elastic  .and.  NormalIn) then

            ! The QP switched to elastic mode.
            ! The linearized constraints are infeasible.

            if (mjrPrtlvl .ge. 1  .or.  mnrPrtlvl .ge. 10) then
               write(str, 1100) itn, wtInf
               call snPRNT( 23, str, iw, leniw )
               iw(mnrHdP) = 1
               iw(mnrHdS) = 1
            end if
         end if
      end if ! .not. Feasible

!     ------------------------------------------------------------------
!     The inelastic variables (x's and linear slacks) are now feasible.
!     Save them in xQP0 for use with the BFGS update.
!
!     Solve the QP subproblem.
!     Loop back sometimes if we need a BS factorize.
!     ------------------------------------------------------------------
      call dcopy ( nb, xQP, 1, xQP0, 1 )

      Feasible   = .true.       ! the nonelastics are feasible

      ! For nonlinear constraints.
      ! Set lvlObjE to use the composite objective  Obj + wtInf*sInf
      ! after any switch to elastic mode.
      ! Linear constraints.  In theory, the subproblem should be
      ! feasible. If it is not, do not switch to Elastic mode.

      if (NonlinearCon) then
         lvlObjE = 1
      else
         eMode   = 0
         lvlObjE = 0
      end if

      if (nnH .gt. 0) then
         subOptimize = SubOpt
      else
         subOptimize = NoSubOpt
      end if

      itQPmax    = itnlim
      itQPtarget = itQP + mMinor
      LUrequest  = 0
      typeLU     = BT

      LUok       = .true.
      Terminate  = .false.
      Solved     = .false.

!     ==================================================================
!     while (.not. (Solved  .or.  Terminate)  .and.  LUok) do
  600 if    (.not. (Solved  .or.  Terminate)  .and.  LUok) then

         inform = 0
         NeedLU = LUrequest .gt. 0
         Needx  = NeedLU

         if (mnrPrtlvl .ge. 1) then
            iw(mnrHdP) = 1      ! QP print   header
            iw(mnrHdS) = 1      ! QP summary header
         end if

         if (iw(QPmode) .eq. QPChol) then
!           -----------------------------------------------------------
!           Solve the QP subproblem using the Cholesky QP solver..
!           maxS = maxR to force termination at SB limit.
!           -----------------------------------------------------------
            if (nS .le. maxR) then
               maxSB  = maxR
            else
               maxSB  = maxS
            end if

            call s5QP
     &         ( inform,
     &           QPS, probTag, subOptimize,
     &           mnrlog, Hprod, Hprod1, HvCalls, eigH,
     &           Elastic, GotR, NeedLU, typeLU, Needx,
     &           lenR, m, maxSB, mBS, n, nb, nDegen,
     &           ngQP0, ngQP, nnObj0, nnObj, nnH0, nnH, nS,
     &           itQP, itQPmax, itQPtarget, itn,
     &           eMode, lvlObjE, mnrPrtlvl,
     &           minimize, iObj, scaleObj, objAdd, objQP,
     &           condZmax, condZHZmax, tolOptFP, tolOptQPk, tolx,
     &           nInf, sInf, elastics, nInfE, sInfE, wtInf,
     &           piNorm, rgNorm,
     &           neJ, nlocJ, locJ, indJ, Jcol,
     &           neH, nlocH, locH, indH, Hcol,
     &           eType, eState, feasType, hs, kBS,
     &           bl, bu, blQP, buQP, blBS, buBS,
     &           gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg,
     &           nnCon0, nnCon, QPrhs, scales,
     &           nnH0, nnH, x, xQP, xBS, x,
     &           iy, iy1, y, y1, y2,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

         else if (iw(QPmode) .eq. QN  .or.  iw(QPmode) .eq. CG) then
!           ------------------------------------------------------------
!           Solve the QP using the quasi-Newton or CG solver..
!           maxS = maxR to avoid automatic switch to CG mode
!           ------------------------------------------------------------
            if (iw(QPmode) .eq. QN  .and.  nS .lt. maxR) then
               maxSB = maxR
            else
               maxSB = maxS
            end if

            call s5QN
     &         ( inform,
     &           QPS, probTag, subOptimize,
     &           mnrlog, Hprod, Hprod1, HvCalls,
     &           Elastic, GotR, NeedLU, typeLU, Needx,
     &           lenR, m, maxSB, mBS, n, nb, nDegen,
     &           ngQP0, ngQP, nnObj0, nnObj, nnH0, nnH, nS,
     &           itQP, itQPmax, itQPtarget, itn,
     &           eMode, lvlObjE, mnrPrtlvl,
     &           minimize, iObj, scaleObj, objAdd, objQP,
     &           condZHZ, condZmax, tolOptFP, tolOptQPk, tolx,
     &           nInf, sInf, elastics, nInfE, sInfE, wtInf, piNorm,
     &           neJ, nlocJ, locJ, indJ, Jcol,
     &           neH, nlocH, locH, indH, Hcol,
     &           eType, eState, feasType, hs, kBS,
     &           bl, bu, blQP, buQP, blBS, buBS,
     &           gBS, gObj, gQP, Hdx, pBS, pi, R, rc, rg, rg2,
     &           nnCon0, nnCon, QPrhs, scales,
     &           nnH0, nnH, x, xQP, xBS, x,
     &           iy, iy1, y, y1, y2,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )

         end if ! QPSolver = QN or CG

!        inform  Meaning
!        ------  -------
!         >0     Fatal LU error or time limit
!          0     QP solution found
!         -1     The nonelastics are infeasible
!         -2     The QP subproblem is unbounded
!         -3     Too many iterations
!         -4     The QP subproblem has a weak minimizer   (s5QP)
!         -5     Too many superbasics
!         -6     Reduced Hessian not psd with new column  (s5QP)
!         -7     Z'g could not be made sufficiently small (s5QP)
!         -8     Ill-conditioned Z                        (s5QN)
!        -10     Too many CG subspace iterations          (s5QN)

         if (inform .gt. 0) then
            iExit = inform      ! Fatal LU error or time limit
            go to 900
         end if

         Solved    = inform .eq.  0
         Terminate = inform .eq. -3

         ReSolve   = inform .eq. -1  .or.
     &               inform .eq. -2  .or.
     &               inform .eq. -4  .or.
     &               inform .eq. -5  .or.
     &               inform .eq. -6  .or.
     &               inform .eq. -7  .or.
     &               inform .eq. -8  .or.
     &               inform .eq.-10

         if (Solved) then

            ! Finish if there are no large multipliers.
            ! Otherwise, set elastic mode and solve the QP again,

            if (.not. Elastic .and. NonlinearCon) then
               piNormNLN = dnormj( nnCon, pi, 1 )
               if (piNormNLN .gt. wtInf) then
                  Elastic    = .true.
                  write(str, 1400) itn, wtInf
                  call snPRNT( 23, str, iw, leniw )
                  Solved     = .false.
               end if
            end if

         else if (Terminate) then

            ! Relax

         else if (ReSolve) then

            !===========================================================
            ! Neither solved nor terminated.
            ! inform = -1, -2, -4, -5, -6, -7, -8, -10
            ! Try to solve again with one or more of the following:
            !  (1) with a different H
            !  (2) a different QP solver
            !  (3) or more accurate initial LU.
            !===========================================================

            if (inform .eq. -1) then
               !--------------------------------------------------------
               ! The nonelastics are infeasible. This should not happen
               ! because Phase 1 has already found a feasible point for
               ! the nonelastics. The basis must be ill-conditioned.
               ! Refactorize with tighter tols and restart at the known
               ! feasible point.
               !--------------------------------------------------------
               call dcopy ( nb, xQP0, 1, xQP, 1 )

               write(str, 2001) itn
               call snPRNT( 23, str, iw, leniw )
               call s2tryLU
     &            ( itn, 20, nS, LUrequest, LUok, typeLU,
     &              iw, leniw, rw, lenrw )
               if (.not. LUok) then
                  iExit = 15    ! infeasible linear constraints in QP.
                  go to 900
               end if

               Elastic   = .false.
               Feasible  = .false.
               go to 100

            else if (inform .eq. -2) then
!              ---------------------------------------------------------
!              The QP is unbounded.
!              ------------------------------------------------------
               write(str, 2002) itn
               call snPRNT( 23, str, iw, leniw )

               if (iw(QPmode) .eq. QPChol) then

!                 An unbounded QP => Z'HZ is positive semidefinite.
!                 Refactor for safety.

                  if (iw(LUitn) .gt. 0) then
                     call s2tryLU
     &                  ( itn, 22, nS, LUrequest, LUok, typeLU,
     &                    iw, leniw, rw, lenrw )
                     if (.not. LUok) then
                        Terminate = .true.
                     end if
                  else
                     Terminate = .true.
                  end if
               else if (iw(QPmode) .eq. QN) then

!                 The approximate Hessian is positive definite.
!                 Reset both the full and reduced Hessian.

                  if (iw(HDInfo) .ne. Unit) then
                     if (nnH .gt. 0) then
                        iw(HDInfo) = Unit
                        call s8ResetH
     &                     ( itn, nnH, UD0, HD, iw, leniw, rw, lenrw )
                     end if
                     GotR = .false.
                  end if

                  call s2tryLU
     &               ( itn, 25, nS, LUrequest, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )
                  if (.not. LUok) then
                     Terminate = .true.
                  end if
               else if (iw(QPmode) .eq. CG) then

!                 No Hessian in this case.
!                 Bad LU is likely giving a large reduced gradient

                  call s2tryLU
     &               ( itn, 25, nS, LUrequest, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )
                  if (.not. LUok) then
                     Terminate = .true.
                  end if
               end if

            else if (inform .eq. -4) then
!              ---------------------------------------------------------
!              Weak QP minimizer.
!              ------------------------------------------------------
               if (iw(LUitn) .gt. 0) then
                  call s2tryLU
     &               ( itn, 27, nS, LUrequest, LUok, typeLU,
     &                 iw, leniw, rw, lenrw )
                  if (.not. LUok) then
                     Terminate = .true.
                  end if
               else
                  Terminate = .true.
               end if

            else if (inform .eq. -5) then
!              ---------------------------------------------------------
!              Too many superbasics.
!              Switch to CG mode if possible.
!              ------------------------------------------------------
               if (maxR .lt. maxS) then
                  iw(QPmode) = CG ! Switch to CG
                  iw(PreCon) = NO ! with no preconditioning
                  GotR       = .false.
                  iw(QPInfo) = 0
               else
                  Terminate  = .true.
               end if

            else if (inform .eq. -8) then
!              ---------------------------------------------------------
!              condZ > condZmax0  while computing Z'HZ.
!              Refactorize B, possibly with a reduced factor tol.
!              ---------------------------------------------------------
               write(str, 2008) itn, condZmax
               call s5setCondZmax( condZmax, rw, lenrw )
               call snPRNT( 23, str, iw, leniw )
               call s2tryLU
     &            ( itn, 25, nS, LUrequest, LUok, typeLU,
     &              iw, leniw, rw, lenrw )
               if (.not. LUok) then
                  iExit = 44    ! ill-conditioned null-space basis
                  go to 900
               end if

            else if (inform .eq. -6  .or.  inform .eq. -7 .or.
     &               inform .eq.-10) then
!              ---------------------------------------------------------
!              Either Z'HZ is not positive definite or Z'g is large.
!              Most likely, Z'HZ is ill-conditioned.
!              Refactorize B, possibly with a reduced factor tol. If
!              the factor tol is already tight, and H is not a multiple
!              of the identity, accept Z, however bad it is.
!
!              Discard the off-diagonals of H, or, if H is already
!              diagonal, reset it to a multiple of the identity.
!              ----------------------------------------------------------
               if (     inform .eq. -6) then
                  write(str, 2006) itn
               else if (inform .eq. -7) then
                  write(str, 2007) itn
               else if (inform .eq.-10) then
                  write(str, 2010) itn
               end if
               call snPRNT( 23, str, iw, leniw )

               if (     inform .eq. -6) then
                  LUrequest = 25
               else if (inform .eq. -7) then
                  LUrequest = 21
               end if

               call s2tryLU
     &            ( itn, LUrequest, nS, LUrequest, LUok, typeLU,
     &              iw, leniw, rw, lenrw )

               if (.not. LUok) then
                  if (iw(HDInfo) .eq. Unit) then
!                    H = I, so Z'Z must not be positive definite.
                     iExit = 44 ! ill-conditioned null-space basis
                     go to 900
                  else
!                    Relax, and assume that resetting H will fix things.
                  end if
               end if

               if (iw(HDInfo) .ne. Unit) then

!                 Discard the off-diagonals of H.
!                 If H is already diagonal, H is set to the identity.

                  call s8ResetH
     &               ( itn, nnH, UD0, HD, iw, leniw, rw, lenrw )
                  GotR  = .false.
               end if
            end if ! inform ne 0
         end if

         go to 600
      end if
    ! end while not solved or terminated
!     ------------------------------------------------------------------
  800 if (nInfE .gt. 0) iw(FPInfo) =  1

      if (inform .eq. 0) then
         iw(QPInfo) = max(subOptimize,0) ! = 0,1 or 2
      else if (inform .eq. -1) then
         iExit      = 15       ! infeasible nonelastics
      else if (inform .eq. -2) then
         iw(QPInfo) =  3       ! unbounded subproblem
      else if (inform .eq. -3) then
         iw(QPInfo) =  2
         iExit      = -1       ! too many iterations
      else if (inform .eq. -4) then
         iw(QPInfo) =  4       ! weak QP solution
      else if (inform .eq. -5  .and.  Feasible) then
         iw(QPInfo) =  5
         iExit      = -2       ! superbasic limit
      else if (inform .eq. -5) then
         iExit      = 33       ! superbasic limit
      end if

  900 return

 1100 format(' Itn', i7, ': Infeasible subproblem.',
     &       ' Elastic mode started with weight = ', 1p, e8.1)
 1400 format(' Itn', i7, ': Large multipliers.',
     &       ' Elastic mode started with weight = ', 1p, e8.1)
 1500 format(' Itn', i7, ': Infeasible nonelastics in QP feasibility',
     &                    ' phase')
 2001 format(' Itn', i7, ': Infeasible nonelastics in QP optimality',
     &                    ' phase')
 2002 format(' Itn', i7, ': Unbounded QP subproblem')
 2006 format(' Itn', i7, ': Indefinite QP reduced Hessian')
 2007 format(' Itn', i7, ': Large QP reduced gradient')
 2008 format(' Itn', i7, ': Ill-conditioned QP null-space basis.',
     &                    ' Cond Z = ', 1p, e8.1)
 2010 format(' Itn', i7, ': Too many CG subspace iterations.')

      end ! subroutine s8solveQP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8checkLC
     &   ( labelstring,
     &     n, nb, nnCon, bl, bu, x, iw, leniw, rw, lenrw )

      implicit
     &     none
      character*(*)
     &     labelstring
      integer
     &     n, nb, nnCon, leniw, lenrw, iw(leniw)
      double precision
     &     bl(nb), bu(nb), x(nb), rw(lenrw)

!     ==================================================================
!     Checks that the linear constraints and bounds are satisfied.
!     ==================================================================
      integer
     &     i, iPrint, iSumm, j, length, nLinear, npCon
      double precision
     &     tolx
!     ------------------------------------------------------------------
      character          form*4
      parameter        ( form  = '(/a)')
!     ------------------------------------------------------------------
      tolx     = rw( 56) ! Minor feasibility tolerance
      iPrint   = iw( 12) ! Print (Log) file
      iSumm    = iw( 13) ! Summary file

!     length = len_trim(labelstring)     ! An F90 intrinsic
      call s1trim( labelstring, length ) ! The F77 equivalent

      if (iPrint .gt. 0)
     &write(iPrint, form) labelstring(1:length)
      if (iSumm  .gt. 0)
     &write(iSumm , form) labelstring(1:length)

      do j = 1, n
         if (x(j) .lt. bl(j) - tolx .or.
     &       x(j) .gt. bu(j) + tolx) then
            write (iPrint,1000) '    variable outside its bound',
     &                        j, x(j), bl(j), bu(j)
            write (iSumm ,1000) '    variable outside its bound',
     &                        j, x(j), bl(j), bu(j)

         end if
      end do

      npCon     = n  + nnCon
      nLinear   = nb - npCon
      do j      = 1, nLinear
         i      = npCon + j
         if (x(i) .lt. bl(i) - tolx .or.
     &       x(i) .gt. bu(i) + tolx) then
            write (iPrint,1000) 'linear slack outside its bound',
     &                        i, x(i), bl(i), bu(i)
            write (iSumm ,1000) 'linear slack outside its bound',
     &                        i, x(i), bl(i), bu(i)
         end if
      end do

 1000 format( 1p, a, i10, 3e24.14 )

      end ! subroutine s8CheckLC
