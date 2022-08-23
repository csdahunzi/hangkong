!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn82qn.f
!
!     s8getH     s8ResetH   s8HQN    s8Hwrapper   s8Hx   s8xHx
!     s8Hupdate  s8HmodA    s8HmodB
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8getH
     &   ( nnH, lenH, U, H, y, y1 )

      implicit
     &     none
      integer
     &     lenH, nnH
      double precision
     &     H(lenH), U(lenH), y(nnH), y1(nnH)

!     ==================================================================
!     s8getH  computes the product H = U'U, where  U is the Cholesky
!     factor of the approximate Hessian of the Lagrangian.  The matrix
!     U is stored by rows in the one-dimensional array  U.
!     lenH defines the length of U.  lenH must be at least
!     nnH*(nnH + 1)/2.  The result is stored by columns in the upper
!     triangular array H.
!
!     03 Sep 2006: First version of s8getH.
!     03 Sep 2006: Current version.
!     ==================================================================
      integer
     &     j, jthcol
!     ------------------------------------------------------------------
      integer            WithU,      WithUt
      parameter         (WithU  = 0, WithUt = 1)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one   = 1.0d+0)
!     ------------------------------------------------------------------
!     Compute the product y1 = U'Uy, where  U is an upper-
!     ------------------------------------------------------------------
      jthcol = 1

      do     j = 1, nnH
         jthcol = jthcol + j - 1

         call dload ( nnH, zero, y, 1 )
         y(j) = one

         call s6Rprod( WithU , nnH, nnH, lenH, U,  y, y1 )
         call s6Rprod( WithUt, nnH, nnH, lenH, U, y1,  y )
         call dcopy  ( j, y, 1, H(jthcol), 1 )

      end do

      end ! subroutine s8getH

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8ResetH
     &   ( itn, nnH, UD0, HD, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     itn, leniw, lenrw, nnH, iw(leniw)
      double precision
     &     UD0, HD(nnH), rw(lenrw)

!     ==================================================================
!     s8ResetH  initializes or resets  the BFGS approximate Hessian.
!
!     The BFGS Hessian is  H = U'U,  where the form of U is defined
!     according to the definition of lvlHess. Copies of the diagonals of
!     H are held in HD.
!
!     At a reset,  H is replaced by a diagonal matrix.
!
!     If HD is well-conditioned, H is reset to HD (i.e., the off
!     diagonals of H are discarded). This strategy implies that HD is
!     always the diagonal of H,  regardless of resets.
!
!     If HD is ill-conditioned, H and HD are reset to the scaled
!     identity matrix.
!
!     On entry,
!     ---------
!      UD0    is the new diagonal of U if U is reset to a diagonal.
!             Its value reflects the scale of the gradient and is set
!             in s8SQP.
!
!      HD     are the diagonal elements of the BFGS QN Hessian.
!             i.e.,  HD = diag(U*U)
!
!      HDInfo defines the input status of H.
!             UnSet  (-1)  => U not set
!             Normal ( 0)  => U'U is a non-diagonal BFGS approx
!             Diag   ( 1)  => U'U is general diagonal
!             Unit   ( 2)  => U'U is a multiple of the identity
!
!     On exit
!     -------
!      HD     is defined according to the input value of HDInfo
!
!             HDInfo
!             -----
!             UnSet   =>    U =     U0*I,   HD = (U0*U0)*I
!             Normal  =>    U = sqrt(HD),   HD   unchanged if HD is good
!                        or U =     U0*I,   HD = (U0*U0)*I if HD is bad
!             Diag    =>    U =     U0*I,   HD = (U0*U0)*I
!             Unit    =>    U =     U0*I,   HD = (U0*U0)*I
!
!     HDInfo  indicates what happened.
!             Diag    =>    normal reset,   HD not changed.
!             Unit    =>    Unit   reset or HD was bad.
!
!             In all cases, diag(U'U) = HD
!
!     iw(QNmods), the number of BFGS updates since the last reset,
!                 is set to 0.
!
!     19 Nov 2012: First version based on SNOPT routine s8H0.
!     20 May 2015: Reset info printed based on itn.
!     20 May 2015: itn added as argument.
!     ==================================================================
      character
     &     str*132
      logical
     &     OverFlow
      integer
     &     j, lenU, lvlHess, lU0, lU, HDInfo, HDtype, QNmods
      double precision
     &     condHD, condUmax, ddiv, HDmax, HDmin
!     ------------------------------------------------------------------
      integer            Unset,      Normal,     Diag,      Unit
      parameter         (Unset = -1, Normal = 0, Diag  = 1, Unit  = 2)

      double precision   zero
      parameter         (zero   = 0.0d+0)

      integer            LM   ,      FM,     Exact
      parameter         (LM     = 0, FM = 1, Exact = 2)

      parameter         (HDInfo = 243) ! Approximate Hessian type
      parameter         (QNmods = 381) ! BFGS updates since last reset
!     ------------------------------------------------------------------
      condUmax = rw( 87) ! max cond estimator for U with H = U'U
      lvlHess  = iw( 72) ! LM, FM or Exact Hessian

      HDtype   = iw(HDInfo)

      if (HDtype .eq. Normal) then
!        ---------------------------------------------------------------
!        Try and set U so that U'U = HD, where HD is the diagonal of the
!        Hessian.  First, check the condition of HD and reset it to
!        UD0*UD0*I  if its ill-conditioned or not positive definite.
!        ---------------------------------------------------------------
         OverFlow = .false.

         HDmin  = HD(1)         ! strictly positive in exact arithmetic
         HDmax  = HDmin
         do j = 2, nnH
            HDmin = min( HD(j), HDmin)
            HDmax = max( HD(j), HDmax)
         end do

         condHD = ddiv( HDmax, HDmin, OverFlow )

         if (HDmin .le. zero  .or.  condHD .ge. condUmax*condUmax) then
            HDtype = Unit
         end if
      end if

      if (HDtype .eq. Normal) then
         write(str, 1000) itn
         iw(HDInfo) = Diag      ! Set U to sqrt(HD)
      else                      ! HDInfo == Unset, Diag or Unit
         write(str, 2000) itn
         iw(HDInfo) = Unit      ! Set diag(U) to UD0*I
         call dload ( nnH, (UD0*UD0), HD, 1 )
      end if
      call snPRNT( 23, str, iw, leniw )

!     ------------------------------------------------------------------
!     Zero the off-diagonal elements of U.
!     How this is done depends on the way that U is stored.
!     ------------------------------------------------------------------
      if      (lvlHess .eq. LM) then
!        -----------------------
!        Limited-memory Hessian.
!        -----------------------
         lU0  = iw(346) ! Square root of initial BFGS diagonal

         call s8LMH0( nnH, HD, rw(lU0) )

      else if (lvlHess .eq. FM) then
!        -----------------------
!        Full-memory Hessian.
!        -----------------------
         lU   = iw(391) ! U(lenU), full-memory BFGS H = U'U
         lenU = iw(392) !

         call s8FMH0( nnH, HD, lenU, rw(lU) )

      end if

      iw(QNmods) = 0

 1000 format(' Itn', i7, ': Hessian off-diagonals discarded')
 2000 format(' Itn', i7, ': Hessian set to a scaled identity matrix')

      end ! subroutine s8ResetH

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8HQN
     &   ( iExit,
     &     funwrapper, funcon, funobj, userHv,
     &     useFD, QPtype, startType,
     &     itn, lenR, m, mBS, n, nb,
     &     nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &     nS, nMajor, QNskips, UD0,
     &     minimize, step, dxHdx,
     &     RtRmods, GotR, penParm,
     &     fObj, fCon, gCon, gObj, fCon1, gCon1, gObj1,
     &     neJ, nlocJ, locJ, indJ, Jcol,
     &     neH, nlocH, locH, indH, Hcol,
     &     negCon, nlocG, locG,
     &     kBS, bl, bu, dx, dg, Udx, Hdx, HD,
     &     yCon1, R, x, x1, xQP0, xPen, y, y1, y2,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     funwrapper, funcon, funobj, userHv
      logical
     &     GotR, useFD
      integer
     &     iExit, itn, lencu, lencw,
     &     leniu, leniw, lenru, lenrw, lenR, mBS,
     &     minimize, m, n, nb, neJ, neH, negCon, nlocG, nlocH,
     &     nlocJ, nMajor, nnCon0, nnCon, nnJac, nnH, nnObj0,
     &     nnObj, nS, QNskips, QPtype, RtRmods, startType,
     &     kBS(mBS), locG(nlocG), locH(nlocH), locJ(nlocJ),
     &     indH(neH), indJ(neJ), iu(leniu), iw(leniw)
      double precision
     &     dxHdx, UD0, penParm(4), step, bl(nb), bu(nb),
     &     dg(nnH), dx(nnH), HD(nnH), Hdx(nnH),
     &     Hcol(neH), Jcol(neJ), yCon1(nnCon0), fObj, fCon(nnCon0),
     &     gCon(negCon), gObj(nnObj0), fCon1(nnCon0), gCon1(negCon),
     &     gObj1(nnObj0), R(lenR), Udx(nnH),
     &     x(nb), x1(nnH), xPen(nnCon0), xQP0(nb),
     &     y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s8HQN  does the quasi-Newton update with vectors
!        dx = x1 - x   and   dg = gL(x1) - gL(x).
!
!     On entry:
!       xQP0    holds the first QP feasible point for the nonelastics.
!
!     23 Apr 1999: First version of s8HQN,
!     18 Feb 2001: LM H stored in product form.
!     12 Oct 2003: snEXIT and SNPRNT adopted
!     10 Jan 2005: FM H stored in product form.
!     23 Jun 2008: y3 no longer an argument.
!     19 Oct 2014: Added user-defined Hessian.
!     27 Apr 2015: s8HmodA checks the new base point for feasibility.
!     27 Apr 2015: x and xQP0 now of length nb
!     20 May 2015: itn added as argument.
!     ==================================================================
      external
     &     ddot, ddiv, dnrm2
      logical
     &     NonlinearCon, OverFlow, Updated
      integer
     &     inform, kFac, maxR,
     &     mQNskips, nBS, nzero, QNInfo, MdInfo, HDInfo
      double precision
     &     ddot, ddiv, dnrm2, eps0, eps1, gLnorm,
     &     PenUnm, rdxHdx, rydx, rnnH, signObj, sNorm, U0scale,
     &     xPen0, ydx, ydxmin
!     ------------------------------------------------------------------
      integer            Normal
      parameter         (Normal = 0)
      integer            QPChol
      parameter         (QPChol = 0)
      integer            Transp
      parameter         (Transp = 1)
      integer            HOT
      parameter         (HOT    = 3)
      integer            mGap
      parameter         (mGap   = 2)
      double precision   tolg,             tolg2
      parameter         (tolg   =  1.0d-3, tolg2  = 1.0d-1)
      double precision   U0max,            U0min
      parameter         (U0max  =  1.0d+1, U0min  = 1.0d-2)
      double precision   zero,             one
      parameter         (zero   =  0.0d+0, one    = 1.0d+0)

      parameter         (QNInfo = 237) ! Tags(1)
      parameter         (MdInfo = 238) ! Hessian mod type
      parameter         (HDInfo = 243) ! Approximate Hessian type
!     ------------------------------------------------------------------
      maxR     = iw( 52) ! max columns of R
      kFac     = iw( 59) ! factorization frequency
      mQNskips = iw( 67) ! # largest allowable  QNskips

      eps0     = rw(  2) ! eps**(4/5)
      eps1     = rw(  3) ! eps**(2/3)
      xPen0    = rw( 89) ! initial penalty parameter.

      iExit    = 0

      OverFlow     = .false.
      NonlinearCon = nnCon  .gt. 0

      nBS        = m + nS
      signObj    = minimize

      iw(QNInfo) = Normal
      iw(MdInfo) = Normal

      ydx        = zero

!     ---------------------------------------------------------------
!     Compute  dx = x1 - x  and  dg = gL1 - gL.
!     Compute the approx. curvature ydx and new scale factor U0.
!     ---------------------------------------------------------------
      call dcopy ( nnH,          x1, 1, dx, 1 )
      call daxpy ( nnH, (-one),  x , 1, dx, 1 )
      call dscal ( nnH,   step, Hdx, 1 )
      call dscal ( nnH,   step, Udx, 1 )
      dxHdx   = dxHdx*step*step

      if (nnObj .gt. 0) then
         call dcopy ( nnObj, gObj1, 1, dg, 1 )
         if (minimize .lt. 0) then
            call dscal ( nnObj, signObj, dg, 1 )
         end if
      end if

      nzero = nnH - nnObj
      if (nzero .gt. 0) then
         call dload ( nzero, zero, dg(nnObj+1), 1 )
      end if

      if (nnCon  .gt. 0) then
         call s8Gprod
     &      ( Transp, eps0,
     &        neJ, nlocJ, locJ, indJ, negCon, nlocG, locG, gCon1,
     &        (-one), yCon1, nnCon, one, dg, nnJac )
      end if

!     gLnorm = dnormi( nnH, dg, 1 )
      gLnorm = dnrm2 ( nnH, dg, 1 )

      if (nnObj .gt. 0) then
         call daxpy ( nnObj, (-signObj), gObj , 1, dg, 1 )
      end if

      if (nnCon  .gt. 0) then
         call s8Gprod
     &      ( Transp, eps0,
     &        neJ, nlocJ, locJ, indJ, negCon, nlocG, locG, gCon,
     &        one , yCon1, nnCon, one, dg, nnJac )
      end if

      ydx  = ddot ( nnH, dg, 1, dx, 1 )

      if (nMajor .eq. 1  .and.  startType .ne. HOT) then
!        ===============================================================
!        The BFGS update is not applied on the first iteration, but the
!        latest curvature information is used to get a better scaled H.
!        ===============================================================
         if (gLnorm .gt. zero) then
            rnnH = nnH
            UD0  = sqrt(gLnorm/sqrt(rnnH))
         else
            UD0  = one
         end if
         UD0    = min( max( UD0, U0min ), U0max )
         call s8ResetH
     &      ( itn, nnH, UD0, HD, iw, leniw, rw, lenrw )
         GotR   = .false.

      else
!        ===============================================================
!        Except on the first iteration, attempt a BFGS update.
!        Compute the smallest allowable curvature.
!        If the update cannot be done, s8HmodA attempts to find a
!        modified update using  dx = x1 - x defined with a new x.
!        Arrays fCon, gCon and gObj must be redefined at the new x.
!        ===============================================================
         sNorm  = dnrm2 ( nnH, dx, 1 )
         UD0    = ddiv  ( sqrt(abs(ydx)), sNorm, OverFlow )
         UD0    = min   ( max( UD0, U0min ), U0max )

         PenUnm  = zero
         ydxmin  = tolg*dxHdx
         Updated =  dxHdx .gt. zero   .and.
     &              ( ydx .ge. ydxmin  .or.   ydx .ge. eps1)

         if (NonlinearCon  .and. .not. Updated) then
!           ------------------------------------------------------------
!           Redefine the base point x and get new  dx, Hdx and dg.
!           The problem functions are recomputed at x.
!           ------------------------------------------------------------
            call s8HmodA
     &         ( inform, funwrapper, funcon, funobj, userHv, useFD,
     &           n, nb, nnCon0, nnCon, nnJac, nnObj0, nnObj, nnH,
     &           minimize, step, dxHdx, ydx,
     &           fObj, fCon, gCon, gObj, gCon1, gObj1,
     &           neJ   , nlocJ, locJ, indJ,
     &           neH   , nlocH, locH, indH, Hcol,
     &           negCon, nlocG, locG,
     &           bl, bu, dx, dg, Udx, Hdx, yCon1,
     &           y1, x, xQP0, y,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (inform .gt. 0) then
               iExit = inform   ! User wants to stop
               go to 999
            end if

            ydxmin  = tolg*dxHdx
            Updated = dxHdx .gt. zero    .and.
     &                ( ydx .ge. ydxmin  .or.   ydx .ge. eps1)

            if (Updated) then
               iw(MdInfo) = 1 ! Mod A  succeeded
            end if

            if (.not. Updated  .and.  dxHdx .gt. zero  ) then
!              ---------------------------------------------------------
!              If all else fails, attempt to update the Hessian of
!              the augmented Lagrangian.
!              The target ydx is defined via tolg2.
!              ---------------------------------------------------------
               ydxmin = tolg2*dxHdx
               call s8HmodB
     &            ( nnCon, nnJac, eps0,
     &              neJ, nlocJ, locJ, indJ, negCon, nlocG, locG,
     &              ydx, ydxmin, PenUnm, fCon, fCon1, gCon, gCon1,
     &              dx, dg, y, y1, y2 )

               Updated = ydx .ge. ydxmin

               if (Updated) then
                  iw(MdInfo) = 2 ! Mod A + B succeeded
               end if
            end if
         end if ! NonlinearCon

         if (Updated) then
!           ------------------------------------------------------------
!           Update the approximate Hessian using (dg,Hdx).
!           If there are no nonlinear constraints,  apply the update
!           to the reduced Hessian.
!           ------------------------------------------------------------
            QNskips = 0

            if (ydx .ge. ydxmin  .and.  iw(HDInfo) .eq. Normal) then
               iw(QNInfo) = 0 ! conventional BFGS
            else
               iw(QNInfo) = 1 ! self-scaled  BFGS
            end if

            rydx   = sqrt(   ydx )
            rdxHdx = sqrt( dxHdx )

            call s8Hupdate
     &         ( iw(QNInfo), itn, nnH,
     &           UD0, U0scale, rydx, rdxHdx, HD, dx, Hdx, dg,
     &           iw, leniw, rw, lenrw )

            if (     QPtype .eq. QPChol) then
               GotR = GotR  .and.  nnCon      .eq. 0
     &                      .and.  nS         .gt. 0
     &                      .and.  iw(HDInfo) .eq. Normal
     &                      .and.  RtRmods    .lt. kfac
            else
               GotR = GotR  .and.  nS         .gt. 0
     &                      .and.  iw(HDInfo) .eq. Normal
            end if

            if (GotR) then
               call s6Rupdate
     &            ( iw(QNInfo), maxR, lenR, m, n, nBS, nnH, nS,
     &              U0scale, rdxHdx, neJ, nlocJ, locJ, indJ, Jcol,
     &              kBS, dg, Hdx, R, y, y1, y2,
     &              iw, leniw, rw, lenrw )
               RtRmods = RtRmods + 1
            else
               RtRmods = 0
            end if
         else
!           ------------------------------------------------------------
!           No suitable update pair (dg,Hdx) could be found.
!           Skip the update.  Too many skips and we reset.
!           ------------------------------------------------------------
            QNskips = QNskips  + 1

!           Apply all updates to H and discard the off-diagonals.

            if (mod( QNskips, mQNskips  ) .eq. 0) then
               call s8ResetH
     &            ( itn, nnH, UD0, HD, iw, leniw, rw, lenrw )

               if (mod( QNskips, mGap*mQNskips ) .eq. 0) then
!                 ------------------------------------------------------
!                 Reset the multipliers and penalty parameters
!                 ------------------------------------------------------
                  call s8InitPen
     &               ( nnCon, penParm, xPen0, xPen, rw, lenrw )
                  call dload ( nnCon, zero , yCon1, 1 )
               end if
            end if
         end if
      end if ! nMajor > 1

  999 return

      end ! subroutine s8HQN

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Hwrapper
     &   ( Hprod, nnH,
     &     neH, nlocH, locH, indH, Hcol,
     &     x, Hx, Status,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     Hprod
      integer
     &     lencu, leniu, lenru, lencw, leniw, lenrw, neH, nlocH, nnH,
     &     indH(neH), locH(nlocH), Status, iu(leniu), iw(leniw)
      double precision
     &     Hcol(neH), Hx(nnH), x(nnH), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s8Hwrapper wraps Hprod, which multiplies the QP Hessian H by the
!     vector  x.   It is called by the QP solver.
!
!     On entry:
!        Status  = 0  => a normal call for H*x.
!        Status  = 1  => the first entry for a given QP.
!        Status ge 2  => last call for a given QP. Status = 2+iExit.
!
!     On exit:
!        Status lt 0   the user wants to stop.
!
!     03 Nov 2000: First version of s8Hwrapper.
!     03 Nov 2014: neH, indH, locH, Hcol added as arguments.
!     ==================================================================
      call Hprod
     &   ( nnH, x, Hx, cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine s8Hwrapper

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Hx
     &   ( nnH, x, Hx, cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     nnH, lencw, leniw, lenrw, iw(leniw)
      double precision
     &     Hx(nnH), x(nnH), rw(lenrw)
      character
     &     cw(lencw)*8

!     ==================================================================
!     s8Hx  multiplies the QP Hessian  H by the vector  x.
!     It is used to define Hx for the QP subproblem.
!
!     This routine is called by a general QP solver, which will rescale
!     Hx by  signObj  when maximizing.
!
!     s8Hx calls one of the Hessian routines s8LMH, s8FMH, s8SDH, ...
!     according to the value of the options lvlDer and lvlHess.
!     Each of these routines defines a particular form of the Hessian.
!     At the moment the options are:
!
!        lvlHess = LM      Limited-Memory (LM) BFGS  (the default).
!        lvlHess = FM      Full-Memory    (FM) BFGS
!        lvlHess = Exact   FD or exact Hessian
!
!     30 Dec 1991: First version of s8Hx.
!     12 Jan 1996: Full memory Hessian option added.
!     04 Apr 1999: Exact and FD Hessian option added.
!     18 Feb 2001: LM H stored in product form.
!     10 Jan 2005: FM H stored in product form.
!     23 Jun 2008: Exact option added.
!     ==================================================================
      integer
     &     lenU, lS, lU, lU0, lUx, lV, lvlHess, minimize, mQNmods
      double precision
     &     signObj
!     ------------------------------------------------------------------
      integer            LM   ,      FM,     Exact
      parameter         (LM     = 0, FM = 1, Exact = 2)
      integer            QNmods
      parameter         (QNmods = 381)
!     ------------------------------------------------------------------
      lvlHess   = iw( 72) ! LM, FM or Exact Hessian
      minimize  = iw(199) ! (-1)(+1)    => (max)(min)
      lUx       = iw(345) ! Ux(nnH)     = product of U with x

      signObj   = minimize

      if      (lvlHess .eq. LM) then
!        -----------------------
!        Limited-memory Hessian.
!        -----------------------
         mQNmods   = iw( 54) ! (ge 0) max # of BFGS updates
         lU0       = iw(346) ! Square root of initial BFGS diagonal
         lS        = iw(401) ! sk's for BFGS products: (I + sk*vk')
         lV        = iw(402) ! vk's for BFGS products: (I + sk*vk')

         call s8LMHx
     &      ( nnH, x, rw(lUx), Hx,
     &        mQNmods, iw(QNmods), rw(lU0), rw(lS), rw(lV) )

      else if (lvlHess .eq. FM) then
!        -----------------------
!        Full-memory Hessian.
!        -----------------------
         lU        = iw(391) !
         lenU      = iw(392) !

         call s8FMHx
     &      ( nnH, x, rw(lUx), Hx, lenU, rw(lU) )
      end if

      if (minimize .lt. 0)
     &   call dscal
     &     ( nnH, signObj, Hx, 1 )

      end ! subroutine s8Hx

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8xHx
     &   ( nnH, x, Ux, Hx, xHx, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     nnH, leniw, lenrw, iw(leniw)
      double precision
     &     xHx, Hx(nnH), Ux(nnH), x(nnH), rw(lenrw)

!     ==================================================================
!     s8xHx  computes x'Hx and Hx, where H = U'U.
!
!     s8xHx calls one of the Hessian routines s8LMH, s8FMH, s8SDH, ...
!     according to the value of the options lvlDer and lvlHess.
!     Each of these routines defines a particular form of the Hessian.
!     At the moment the options are:
!
!        lvlHess = LM      Limited-Memory (LM) BFGS  (the default).
!        lvlHess = FM      Full-Memory    (FM) BFGS
!        lvlHess = Exact   FD or exact Hessian
!
!     10 Jan 2005: First version of s8Hx based on s8Hx
!     18 Feb 2001: LM H stored in product form.
!     10 Jan 2005: FM H stored in product form.
!     23 Jun 2008: Exact option added.
!     ==================================================================
      external
     &     ddot
      integer
     &     lenU, lvlHess, lU, lU0, lS, lV, mQNmods
      double precision
     &     ddot
!     ------------------------------------------------------------------
      integer            LM,         FM,     Exact
      parameter         (LM     = 0, FM = 1, Exact = 2)
      integer            QNmods
      parameter         (QNmods = 381)
!     ------------------------------------------------------------------
      lvlHess = iw( 72) ! LM, FM or Exact Hessian

      if      (lvlHess .eq. LM) then
!        -----------------------
!        Limited memory Hessian.
!        -----------------------
         mQNmods   = iw( 54) ! (ge 0) max # of BFGS updates
         lU0       = iw(346) ! Square root of initial BFGS diagonal
         lS        = iw(401) ! sk's for BFGS products: (I + sk*vk')
         lV        = iw(402) ! vk's for BFGS products: (I + sk*vk')

         call s8LMHx
     &      ( nnH, x, Ux, Hx,
     &        mQNmods, iw(QNmods), rw(lU0), rw(lS), rw(lV) )

      else if (lvlHess .eq. FM) then
!        -----------------------
!        Full memory Hessian.
!        -----------------------
         lU        = iw(391) ! U(lenU), dense Hessian factor
         lenU      = iw(392) !

         call s8FMHx
     &      ( nnH, x, Ux, Hx, lenU, rw(lU) )
      end if

      xHx = ddot  ( nnH, Ux, 1, Ux, 1 )

      end ! subroutine s8xHx

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8Hupdate
     &   ( updateType, itn, nnH,
     &     UD0, U0scale, rydx, rdxHdx, HD, dx, Hdx, y,
     &     iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     itn, updateType, nnH, leniw, lenrw, iw(leniw)
      double precision
     &     UD0, U0scale, rydx, rdxHdx, dx(nnH), HD(nnH),
     &     Hdx(nnH), y(nnH), rw(lenrw)

!     ==================================================================
!     s8Hupdate  applies the pair of vectors that define the BFGS update
!     or self-scaled BFGS update.
!
!     On entry:
!     ---------
!     Hdx   contains  H  times the difference x1 - x.
!     y     contains the gradient  difference g1 - g.
!
!     On exit:
!     ---------
!     Hdx   contains  H  times the difference x1 - x(new).
!     y     contains the gradient  difference g1 - g(new).
!
!     s8Hupdate calls one of the Hessian routines s8LMH, s8FMH, s8SDH, ...
!     according to the value of the option lvlHess.
!     At the moment the options are:
!
!        lvlHess = LM      Limited-Memory (LM) BFGS  (the default).
!        lvlHess = FM      Full-Memory    (FM) BFGS
!        lvlHess = Exact   FD or exact Hessian
!
!     19 Jul 1995: First version (s8Hupd) of s8Hupdate.
!     12 Jan 1996: Full-memory Hessian option added.
!     18 Feb 2001: LM H stored in product form.
!     12 Jan 2005: FM H stored in product form.
!     20 May 2015: itn added as argument.
!     ==================================================================
      integer
     &     i, lvlHess, lenU, mQNmods, lU0, lU, lUdx, lS, lV
      double precision
     &     H0scale, Hdxi, yi
!     ------------------------------------------------------------------
      integer            Normal
      parameter         (Normal = 0)
      integer            ssBFGS
      parameter         (ssBFGS = 1)
      integer            LM        , FM
      parameter         (LM     = 0, FM = 1)
      integer            QNmods
      parameter         (QNmods = 381) ! # of updates since last reset
      integer            HDInfo
      parameter         (HDInfo = 243) ! Tags(7): Approx Hessian type
!     ------------------------------------------------------------------
      mQNmods = iw( 54) ! (ge 0) max # of BFGS updates
      lvlHess = iw( 72) ! LM, FM or Exact Hessian

      if (updateType .eq. ssBFGS) then
         U0scale =    rydx/rdxHdx
         H0scale = U0scale*U0scale
         rdxHdx  =  rdxHdx*U0scale   ! Used later for the LC update.

         call dscal ( nnH, H0scale, HD , 1 )
         call dscal ( nnH, H0scale, Hdx, 1 )
      end if

!     Include the latest update in the Hessian diagonal.

      do     i = 1, nnH
         Hdxi  = Hdx(i)/rdxHdx
         yi    =   y(i)/rydx
         HD(i) = HD(i) - Hdxi**2 + yi**2
      end do

      if (iw(QNmods) .ge. mQNmods) then
!        ---------------------------------------------------------------
!        Too many modifications.
!        Reset U to be the root of the diagonal of the current H.
!        ---------------------------------------------------------------
         call s8ResetH
     &      ( itn, nnH, UD0, HD, iw, leniw, rw, lenrw )
      else
!        ---------------------------------------------------------------
!        Apply the update to H = U'U.
!        The form of U depends on the selected implementation type.
!        ---------------------------------------------------------------
         iw(QNmods) = iw(QNmods) + 1
         iw(HDInfo) = Normal

         if (lvlHess .eq. LM) then
            lU0       = iw(346) ! Square root of initial BFGS diagonal
            lS        = iw(401) ! sk's for BFGS products: (I + sk*vk')
            lV        = iw(402) ! vk's for BFGS products: (I + sk*vk')

            call s8LMupdate
     &         ( updateType, iw(QNmods), mQNmods, nnH, U0scale,
     &           rydx, rdxHdx, Hdx, y, dx, rw(lU0), rw(lS), rw(lV) )

         else if (lvlHess .eq. FM) then
            lUdx      = iw(345) ! Ux(nnH)      = product of U with x
            lU        = iw(391) ! U(lenU), full-memory BFGS Hessian H = U'U
            lenU      = iw(392) !

            call s8FMupdate
     &         ( updateType, nnH, U0scale,
     &           rydx, rdxHdx, Hdx, y, rw(lUdx), lenU, rw(lU) )
         end if
      end if

      end ! subroutine s8HUpdate

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8HmodA
     &   ( iExit,
     &     funwrapper, funcon, funobj, userHv, useFD,
     &     n, nb, nnCon0, nnCon, nnJac, nnObj0, nnObj, nnH,
     &     minimize, step, dxHdx, ydx,
     &     fObj, fCon, gCon, gObj, gCon1, gObj1,
     &     neJ   , nlocJ, locJ, indJ,
     &     neH   , nlocH, locH, indH, Hcol,
     &     negCon, nlocG, locG,
     &     bl, bu, dx, dg, Udx, Hdx, yCon1,
     &     tdx, x, xQP0, y,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     funwrapper, funcon, funobj, userHv
      logical
     &     useFD
      integer
     &     iExit, lencu, lencw, leniu, leniw, lenru, lenrw, minimize,
     &     n, nb, negCon, neH, neJ, nlocG, nlocH, nlocJ, nnCon0,
     &     nnCon, nnJac, nnH, nnObj0, nnObj, locG(nlocG), locH(nlocH),
     &     locJ(nlocJ), indH(neH), indJ(neJ), iu(leniu), iw(leniw)
      double precision
     &     dxHdx, fObj, step, ydx,
     &     bl(nb), bu(nb), dg(nnH), dx(nnH),
     &     yCon1(nnCon0), fCon(nnCon0), gCon(negCon), gObj(nnObj0),
     &     gCon1(negCon), gObj1(nnObj0), Hcol(neH), Hdx(nnH),
     &     tdx(nnH), Udx(nnH), x(nb), xQP0(nb), y(nb),
     &     ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s8HmodA   redefines the quantities x and dx, Hdx and dg  used for the
!     quasi-Newton update.  The problem functions are recomputed at x.
!
!     The recomputed  x  is  x + step*(xQP0 - x1),  where xQP0 is a
!     (nonelastic) feasible point from the QP subproblem.
!
!     s8HmodA is always called with nnH > 0.
!
!     02 Dec 1994: First version (s8x1) of s8HmodA.
!     20 Jul 1998: s8HmodA made self-contained
!     24 Aug 1998: Fixed bug found by Alan Brown at Nag.
!                  FD derivatives now computed correctly.
!                  Parameter useFD added.
!     11 Oct 1998: Facility to combine funobj and funcon added.
!     12 Oct 2003: snEXIT and snPRNT adopted.
!     14 Jan 2005: Argument Udx added for call to s8xHx.
!     16 Jun 2008: Call-status implemented correctly
!     19 Oct 2014: Added user-defined Hessian.
!     28 Apr 2015: Check the redefined base point for LC feasibility.
!     28 Apr 2015: x and xQP0 declared length nb
!     ==================================================================
      external
     &     ddot
      logical
     &     LCcheck, NonlinearCon, NonlinearObj
      integer
     &     modefg, nzero
      double precision
     &     ddot, dnormi, dxNorm, eps, eps0, signObj, stepMax
!     ------------------------------------------------------------------
      integer            Transp
      parameter         (Transp = 1)
      double precision   zero,            one
      parameter         (zero  =  0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      eps     = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      eps0    = rw(  2) ! eps**(4/5)

      iExit   = 0
      modefg  = 2
      signObj = minimize

      NonlinearCon  = nnCon  .gt. 0
      NonlinearObj  = nnObj  .gt. 0

!     Save dx in case a better update pair cannot be found.

      call dcopy ( nnH, dx, 1, tdx, 1 )

!     ------------------------------------------------------------------
!     Define a new dx such that dx = x1 - z, where z = x + step*(xF - x)
!     with xF the first QP feasible point (stored in xQP0). The step is
!     chosen as the largest value such that z = x + step*(xF - x) is
!     feasible for the linear constraints and bounds.  The new dx is
!     dx = x1 - z =  dx - step*(xQP0 - x). The vector xQP0 is
!     overwritten by xQP0 - x. The base x is overwritten with the new x.
!     ------------------------------------------------------------------
      call daxpy ( nb, (-one), x, 1, xQP0, 1 )
      dxNorm  = dnormi( n, xQP0, 1 )
      stepMax = step

      call s6LCstepMax
     &   ( nb, n, nnCon, step, stepMax, dxNorm,
     &     bl, bu, x, xQP0, rw, lenrw )

!     Once a feasible z is known, only z(1:nnH) is used for the update.
!     Store the new dx = dx - step*y in xQP0

      call daxpy ( nnH, (-step), xQP0, 1,   dx, 1 )

!     -------------------------------------------------
!     Compute the minimum curvature.
!     If nnH < n, dxHdx may be zero (or negative fuzz).
!     -------------------------------------------------
      call s8xHx
     &   ( nnH, dx, Udx, Hdx, dxHdx, iw, leniw, rw, lenrw )

      if (dxHdx .ge. eps) then

         call daxpy ( nb, step, xQP0, 1, x, 1 )

         LCcheck = .false.
         if (LCcheck) then
            call s8checkLC
     &         ( 's8HmodA',
     &           n, nb, nnCon, bl, bu, x, iw, leniw, rw, lenrw )
         end if

!        -----------------------------------------------
!        Evaluate the functions at the new x.
!        -----------------------------------------------
         call funwrapper
     &      ( iExit,
     &        modefg, NonlinearCon, NonlinearObj,
     &        n, negCon, nnCon0, nnCon,
     &        nnJac, nnH, nnObj0, nnObj,
     &        funcon, funobj, userHv,
     &        x, yCon1,
     &        neJ, nlocJ, locJ, indJ,
     &        neH, nlocH, locH, indH, Hcol,
     &        fCon, fObj, gCon, gObj,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         if (iExit .gt. 0) return

         if (iExit .eq. 0  .and.  useFD) then
            call s6GetMissing
     &         ( iExit,
     &           n, negCon,
     &           nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &           funwrapper, funcon, funobj, userHv,
     &           bl, bu, x, yCon1,
     &           neJ, nlocJ, locJ, indJ,
     &           neH, nlocH, locH, indH, Hcol,
     &           fCon, fObj, gCon, gObj, y,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (iExit .gt. 0) return
         end if

         if (iExit .eq. 0) then
!           ------------------------------------------------------------
!           The functions have been computed at x.
!           ------------------------------------------------------------
            if (nnObj .gt. 0) then
               call dcopy ( nnObj,         gObj1, 1, dg, 1 )
               call daxpy ( nnObj, (-one), gObj , 1, dg, 1 )
               if (minimize .lt. 0) then
                  call dscal ( nnObj, signObj, dg, 1 )
               end if
            end if

            nzero = nnH - nnObj
            if (nzero .gt. 0) call dload ( nzero, zero, dg(nnObj+1), 1 )

            if (nnCon  .gt. 0) then
               call s8Gprod
     &            ( Transp, eps0,
     &              neJ, nlocJ, locJ, indJ, negCon, nlocG, locG, gCon1,
     &              (-one), yCon1, nnCon, one, dg, nnJac )
               call s8Gprod
     &            ( Transp, eps0,
     &              neJ, nlocJ, locJ, indJ, negCon, nlocG, locG, gCon,
     &              one   , yCon1, nnCon, one, dg, nnJac )
            end if
            ydx   = ddot ( nnH, dg, 1, dx, 1 )
         end if
      end if

      if (dxHdx .lt. eps  .or.  iExit .ne. 0) then
         call dcopy
     &      ( nnH, tdx, 1, dx, 1 )
         call s8xHx
     &      ( nnH, dx, Udx, Hdx, dxHdx, iw, leniw, rw, lenrw )
      end if

      end ! subroutine s8HmodA

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8HmodB
     &   ( nnCon, nnJac, tolz,
     &     neJ, nlocJ, locJ, indJ, negCon, nlocG, locG,
     &     ydx, ydxmin, PenUnm, fCon, fCon1, gCon, gCon1,
     &     dx, gd, PenU, v, w )

      implicit
     &     none
      integer
     &     nnCon, nnJac, neJ, negCon, nlocG, nlocJ, indJ(neJ),
     &     locJ(nlocJ), locG(nlocG)
      double precision
     &     ydx, ydxmin, PenUnm, tolz, gCon(negCon), gCon1(negCon),
     &     fCon(nnCon), fCon1(nnCon), PenU(nnCon), v(nnCon), w(nnCon),
     &     dx(nnJac), gd(nnJac)

!     ==================================================================
!     s8HmodB  attempts to find the a vector xPen  of minimum two-norm
!     such that there exists a BFGS update for the modified Lagrangian
!       La   = f(x) - lambda'(fCon1 - LfCon)
!                   + 1/2 (fCon1 - LfCon)'*diag(PenU)*(fCon1 - LfCon),
!
!     where  LfCon = fCon + J(x1)*dx.
!
!     On entry:
!     ---------
!     dx      is the nonlinear part of the search direction x2 - x1.
!     gd      is the Lagrangian gradient difference.
!     gCon    is the Jacobian at the old x.
!     gCon1   is the Jacobian at the new x.
!     ydx     is the approximate curvature of the Lagrangian.
!     ydxmin  (ydx < ydxmin) is the smallest acceptable approximate
!             curvature.
!
!     On exit,
!     ---------
!     gd      is the augmented Lagrangian gradient difference.
!     PenU    are the penalty parameters.
!     ydx     is unchanged unless gotPen is true, in which case
!              ydx = ydxmin.
!
!     08 Dec 1991: First version based on  Npsol  routine npupdt.
!     26 Oct 2000: Current version of s8HmodB.
!     ==================================================================
      external
     &     ddiv, dnrm2
      logical
     &     gotPen, OverFlow
      integer
     &     i
      double precision
     &    beta, ddiv, diff, dnrm2, Peni, wi, wmax, wnorm
!     ------------------------------------------------------------------
      integer            Normal,     Transp
      parameter         (Normal = 0, Transp = 1)
      double precision   PenMax
!-->  parameter         (PenMax = 1.0d+5)
!-->  parameter         (PenMax = 1.0d+16)
      parameter         (PenMax = 1.0d+5)
      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      OverFlow = .false.

!     Try an augmented Lagrangian term to increase ydx.

      PenUnm = zero

!     Compute  v = J1*dx and w = (J2 - J1)*dx = J2*dx - v.

      call s8Gprod
     &   ( Normal, tolz,
     &     neJ, nlocJ, locJ, indJ,
     &     negCon, nlocG, locG, gCon,
     &     one, dx, nnJac, zero, v, nnCon )
      call s8Gprod
     &   ( Normal, tolz,
     &     neJ, nlocJ, locJ, indJ,
     &     negCon, nlocG, locG, gCon1,
     &     one, dx, nnJac, zero, w, nnCon )

      call daxpy ( nnCon, (-one), v, 1, w, 1 )

!     Compute the difference between c and its linearization.
!     v  =  c - cL = fCon1 - (fCon + J1*s) = fCon1 - fCon - J1*s.

      call daxpy ( nnCon, (-one), fCon1, 1, v, 1 )
      call daxpy ( nnCon,   one , fCon , 1, v, 1 )
      call dscal ( nnCon, (-one), v    , 1 )

!     ---------------------------------------------------------
!     Compute the minimum-length vector of penalty parameters
!     that makes the approximate curvature equal to  ydxmin.
!     ---------------------------------------------------------
!     Use w to hold the constraint on PenU.
!     Minimize            norm(PenU)
!     subject to   ( Sum( w(i)*PenU(i) )  =   const,
!                  (           PenU(i)   .ge. 0.

      wmax = zero
      do i    = 1, nnCon
         wi   = w(i)*v(i)
         wmax = max( wmax, wi )
         w(i) = max( zero, wi )
      end do

      wnorm  = dnrm2 ( nnCon, w, 1 )
      diff   = ydxmin - ydx
      beta   = ddiv  ( wmax*diff, wnorm**2, OverFlow )
      gotPen = .not. OverFlow  .and.  wmax .gt. zero
     &                         .and.  beta .lt. PenMax

      if ( gotPen ) then
         beta   = diff/wnorm**2

         do    i = 1, nnCon
            wi   = w(i)
            Peni = beta*wi
            v(i) =       Peni*v(i)
            ydx  = ydx + Peni*wi
            PenU(i) = Peni
         end do
         ydx    = max   ( ydx, ydxmin )
         PenUnm = dnrm2 ( nnCon, PenU, 1 )

!        Update  gd  by the term  (J2' - J1')*v,
!        with v = diag(PenU)*(fCon1 - fCon - J1*s) from above.

         call s8Gprod
     &      ( Transp, tolz,
     &        neJ, nlocJ, locJ, indJ,
     &        negCon, nlocG, locG, gCon1,
     &          one , v, nnCon, one, gd, nnJac )
         call s8Gprod
     &      ( Transp, tolz,
     &        neJ, nlocJ, locJ, indJ,
     &        negCon, nlocG, locG, gCon,
     &        (-one), v, nnCon, one, gd, nnJac )
      end if ! gotPen

      end ! subroutine s8HmodB
