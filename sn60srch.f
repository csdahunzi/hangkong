*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     file  sn60srch.f
*
*     s6getMissing  s6getMissing1  s6LCstepMax
*     s6lineSearch  s6search       s6stepLimits   s6tols   s6Userf
*     srchc         srchq
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6getMissing
     &   ( iExit,
     &     n, negCon, nnCon0, nnCon, nnJac,
     &     nnH, nnObj0, nnObj,
     &     funwrapper, fgcon, fgobj, userHv,
     &     bl, bu, x, yCon,
     &     neJ, nlocJ, locJ, indJ,
     &     neH, nlocH, locH, indH, Hcol,
     &     fCon, fObj, gCon, gObj, y,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     funwrapper, fgcon, fgobj, userHv
      integer
     &     iExit, lencu, leniu, lenru, lencw, leniw, lenrw,
     &     n, neJ, neH, negCon, nnCon0, nnCon, nnJac, nnH,
     &     nnObj0, nnObj, nlocH, nlocJ, locH(nlocH), locJ(nlocJ),
     &     indH(neH), indJ(neJ), iu(leniu), iw(leniw)
      double precision
     &     fObj, bl(n), bu(n), fCon(nnCon0), gObj(nnObj0),
     &     gCon(negCon), Hcol(neH), x(n), yCon(nnCon0), y(nnCon0),
     &     ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s6getMissing   computes any missing objective and constraint gradients
*     for the current value of the variables in  x.
*
*     NOTE --  s6dcon overwrites the first  nnCon  elements of  y
*     if central differences are needed.
*
*     30 Dec 1991: First version (s6fdG) based on Minos routine m6grd.
*     17 Jul 1997: First thread-safe version.
*     11 Oct 1998: s6dcon and s6dobj merged.
*     24 Oct 2000: Updated for SNOPT 6.1
*     12 Oct 2003: snEXIT and snPRNT added.
*     14 Oct 2003: Enforced feasible perturbation.
*     01 Apr 2005: Current version of s6getMissing.
*     19 Oct 2014: Added user-defined Hessian.
*     ==================================================================
      integer
     &     gotFD, llocG, lfCon2, lgConU, lgObjU, nlocG
*     ------------------------------------------------------------------
      gotFD     = iw(183) ! > 0 => some differences needed
      llocG     = iw(260) ! locG(nlocG) = column pointers for indG
      lfCon2    = iw(318) ! fCon2(nnCon) work vector
      lgConU    = iw(319) ! record of unknown derivatives and constants
      lgObjU    = iw(323) ! record of unknown derivatives

      iExit     = 0
      nlocG     = nnJac + 1

      if (gotFD .gt. 0) then
         call s6getMissing1
     &      ( iExit,
     &        n, nnCon0, nnCon, nnJac,
     &        nnH, nnObj0, nnObj,
     &        funwrapper, fgcon, fgobj, userHv,
     &        bl, bu, x, yCon,
     &        neJ   , nlocJ, locJ, indJ,
     &        neH   , nlocH, locH, indH, Hcol,
     &        negCon, nlocG, iw(llocG),
     &        fObj, gObj, fCon, gCon,
     &        rw(lgObjU), rw(lfCon2), rw(lgConU), y,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
         if (iExit .eq. 63) then ! The user didn't like some x's
            call snPRNT(  3,
     &         ' XXX  Unable to apply reversion when differencing',
     &         iw, leniw )
         end if
      end if

      end ! subroutine s6getMissing

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6getMissing1
     &   ( iExit,
     &     n, nnCon0, nnCon, nnJac,
     &     nnH, nnObj0, nnObj,
     &     funwrapper, fgcon, fgobj, userHv,
     &     bl, bu, x, yCon,
     &     neJ, nlocJ, locJ, indJ,
     &     neH, nlocH, locH, indH, Hcol,
     &     negCon, nlocG, locG,
     &     fObj, gObj, fCon, gCon,
     &     gObjU, fConU, gConU, y,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     funwrapper, fgcon, fgobj, userHv
      integer
     &     iExit, lencu, leniu, lenru, lencw, leniw, lenrw,
     &     n, negCon, neH, neJ, nlocG, nlocH, nlocJ, nnCon0, nnCon,
     &     nnJac, nnH, nnObj0, nnObj, indH(neH), indJ(neJ),
     &     locG(nlocG), locH(nlocH), locJ(nlocJ), iu(leniu), iw(leniw)
      double precision
     &     fObj, bl(n), bu(n), fCon(nnCon0), fConU(nnCon0),
     &     gObj(nnObj0), gObjU(nnObj0), gCon(negCon), gConU(negCon),
     &     Hcol(neH), x(n), yCon(nnCon0), y(nnCon0),
     &     ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s6getMissing1  estimates missing elements in the objective
*     gradient and constraint Jacobian using finite differences of the
*     problem functions fObj and fCon.
*
*     The arrays y, fConU, gConU  and gObjU are used as workspace.
*     Dummy elements of gObjU and gConU define the unknown derivatives.
*
*     11 Oct 1998: First version based on combining s6dobj and s6dcon.
*     24 Oct 2000: Updated for SNOPT 6.1
*     12 Oct 2003: snEXIT and snPRNT added.
*     14 Oct 2003: Implemented feasible perturbation.
*     16 Jun 2008: Call-status implemented correctly.
*     15 Nov 2010: Call-status removed from argument list.
*     19 Oct 2014: Added user-defined Hessian.
*     18 Feb 2015: Function evaluations counted correctly.
*     ==================================================================
      integer
     &     inform, ir, j, k, l, kmax, lvlDer, lvlDif, modefg,
     &     nfCon3, nfCon4, nfObj3, nfObj4, numcd, numfd
      double precision
     &     infBnd, buj, delta, dxj, Fback, Fforwd, fdint(2), gdummy,
     &     tolx, xj
      logical
     &     Centrl, Done, Found, NonlinearCon , NonlinearObj,
     &     SomeG, SomeJ
*     ------------------------------------------------------------------
      logical            yes  ,           no
      parameter         (yes    = .true., no   = .false.)
      double precision   one
      parameter         (one    = 1.0d+0)
      double precision   three,           four
      parameter         (three  = 3.0d+0, four = 4.0d+0)

      parameter         (nfCon3 = 191) ! calls to fCon: (forward differencing)
      parameter         (nfCon4 = 192) ! calls to fCon: (central differencing)
      parameter         (nfObj3 = 196) ! calls to fObj: (forward differencing)
      parameter         (nfObj4 = 197) ! calls to fObj: (central differencing)
*     ------------------------------------------------------------------
      lvlDer    = iw( 70) ! = 0, 1, 2 or 3, the derivative level
      lvlDif    = iw(181) !    =1 (2) for forwd (cntrl) diffs

      tolx      = rw( 56) ! Minor feasibility tolerance
      gdummy    = rw( 69) ! definition of an 'unset' value
      infBnd    = rw( 70) ! definition of an infinite bound
      fdint(1)  = rw( 76) ! (1) forwrd diff. interval
      fdint(2)  = rw( 77) ! (2) cntrl  diff. interval

      iExit  = 0

*     The problem functions are called to provide functions only.

      modefg = 0

      SomeG  = lvlDer .eq. 0  .or.  lvlDer .eq. 2
      SomeJ  = lvlDer .eq. 0  .or.  lvlDer .eq. 1

      Centrl = lvlDif .eq. 2
      delta  = fdint(lvlDif)

      numcd  = 0
      numfd  = 0

      do j = 1, nnH

*        Look for the first missing element in this column.

         Found  = no

         NonlinearCon  = j .le. nnJac  .and.  SomeJ
         NonlinearObj  = j .le. nnObj  .and.  SomeG

         if (NonlinearObj) then
            if (gObjU(j) .eq. gdummy) Found = yes
         end if

         if (NonlinearCon) then
            l      = locG(j)
            k      = locJ(j)
            kmax   = locJ(j+1) - 1
            Done   = no

*+          ------------------------------------------------------------
*+          while (k .le. kmax  .and. .not.(Found  .or.  Done)) do
  120       if    (k .le. kmax  .and. .not.(Found  .or.  Done)) then
               ir = indJ(k)
               if (ir .gt. nnCon) then
                  Done = yes
               else
                  if (gConU(l) .eq. gdummy) Found = yes
                  l    = l + 1
               end if
               k  = k + 1
               go to 120
*+          end while
*+          ------------------------------------------------------------
            end if
         end if

         if (Found) then
*           ------------------------------------------------------------
*           Some missing derivatives for this variable.
*           A finite difference is needed.
*           ------------------------------------------------------------
            xj     =  x(j)
            dxj    = delta*(one + abs(  xj ))
            buj    = bu(j)

            if (buj .lt. infBnd  .and.  xj+dxj+dxj .gt. buj+tolx) then
               dxj = -dxj
            end if

            x(j)   = xj    + dxj

            numfd  = numfd + 1

            call funwrapper
     &         ( inform,
     &           modefg, NonlinearCon, NonlinearObj,
     &           n, negCon, nnCon0, nnCon,
     &           nnJac, nnH, nnObj0, nnObj,
     &           fgcon, fgobj, userHv,
     &           x, yCon,
     &           neJ, nlocJ, locJ, indJ,
     &           neH, nlocH, locH, indH, Hcol,
     &           fConu, Fforwd, gConu, gObjU,
     &           cu, lencu, iu, leniu, ru, lenru,
     &           cw, lencw, iw, leniw, rw, lenrw )
            if (inform .ne. 0) then
               if (inform .gt. 0) then
                  iExit = inform
               else
                  iExit = 63   ! unable to move into undefined region
               end if
               go to 999
            end if

            if (Centrl) then
               dxj   =   dxj + dxj
               x(j)  =    xj + dxj
               numcd = numcd + 1

               call funwrapper
     &            ( inform,
     &              modefg, NonlinearCon, NonlinearObj,
     &              n, negCon, nnCon0, nnCon,
     &              nnJac, nnH, nnObj0, nnObj,
     &              fgcon, fgobj, userHv,
     &              x, yCon,
     &              neJ, nlocJ, locJ, indJ,
     &              neH, nlocH, locH, indH, Hcol,
     &              y, Fback, gConu, gObjU,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )
               if (inform .ne. 0) then
                  if (inform .gt. 0) then
                     iExit = inform
                  else
                     iExit = 63 ! unable to move into undefined region
                  end if
                  go to 999
               end if
            end if

            if ( NonlinearObj ) then
               if (gObjU(j) .eq. gdummy) then
                  if (Centrl) then
                     gObj(j) = (four*Fforwd - three*fObj - Fback)
     &                                         / dxj
                  else
                     gObj(j) = (Fforwd - fObj) / dxj
                  end if
               end if
            end if

            if ( NonlinearCon ) then
               l      = locG(j)
               k      = locJ(j)
               Done   = no

*+             ---------------------------------------------------------
*+             while (k .le. kmax  .and.  .not. Done) do
  140          if    (k .le. kmax  .and.  .not. Done) then
                  ir = indJ(k)
                  if (ir .gt. nnCon) then
                     Done = yes
                  else
                     if (gConu(l) .eq. gdummy) then
                        if (Centrl) then
                           gCon(l) = (four*fConu(ir) - three*fCon(ir)
     &                                               - y(ir))/dxj
                        else
                           gCon(l) = (fConu(ir) - fCon(ir))/dxj
                        end if
                     end if
                     l   = l + 1
                  end if

                  k  = k + 1
                  go to 140
*+             end while
*+             ---------------------------------------------------------
               end if
            end if ! j <= nnJac

            x(j)   = xj
         end if ! Found
      end do

*     ------------------------------------------------------------------
*     The missing derivatives have been estimated.
*     Finish up with some housekeeping.
*     ------------------------------------------------------------------
      iw(nfCon3) = iw(nfCon3) + numfd
      iw(nfCon4) = iw(nfCon4) + numcd
      iw(nfObj3) = iw(nfObj3) + numfd
      iw(nfObj4) = iw(nfObj4) + numcd

  999 return

      end ! subroutine s6getMissing1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6LCstepMax
     &   ( nb, n, nnCon, step, stepMax, dxNorm,
     &     bl, bu, x, dx, rw, lenrw )

      implicit
     &     none
      integer
     &     nb, n, nnCon, lenrw
      double precision
     &     step, stepMax, dxNorm,
     &     bl(nb), bu(nb), x(nb), dx(nb), rw(lenrw)

!     ==================================================================
!     s6LCstepMax  finds the maximum feasible step for the bounds and
!     lineear constraints.
!
!     On entry
!     --------
!       x  is feasible for the linear constraints.
!
!     26 Apr 2015: First version.
!     ==================================================================
      integer
     &     i, j, nLinear, npCon
      double precision
     &     eps, pivot, pivabs, res, tolp, tolx
!     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
!     ------------------------------------------------------------------
      tolx   = rw( 56) ! Minor feasibility tolerance
      eps    = rw(  1) ! machine precision.  IEEE DP  2.22e-16

!     ==================================================================
!     stepMax  is the largest feasible steplength.
!     step     is initialized as stepMax.
!     ==================================================================
      tolp   = eps*dxNorm
      step   = stepMax

      do j      = 1, n
         pivot  = dx(j)
         pivabs = abs( pivot )
         if (pivabs .gt. tolp) then
            if (pivot  .le. zero  ) then
               res = x(j) - bl(j) + tolx
               if (step*pivabs .gt. res) step = res / pivabs
            else
               res = bu(j) + tolx - x(j)
               if (step*pivabs .gt. res) step = res / pivabs
            end if
         end if
      end do

      npCon     = n  + nnCon
      nLinear   = nb - npCon
      do i      = 1, nLinear
         j      = npCon + i
         pivot  = dx(j)
         pivabs = abs( pivot )
         if (pivabs .gt. tolp) then
            if (pivot  .le. zero  ) then
               res = x(j) - bl(j) + tolx
               if (step*pivabs .gt. res) step = res / pivabs
            else
               res = bu(j) + tolx - x(j)
               if (step*pivabs .gt. res) step = res / pivabs
            end if
         end if
      end do

      end ! subroutine s6LCstepMax

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6lineSearch
     &   ( iExit,
     &     funwrapper, fgcon, fgobj, userHv,
     &     Elastic, FonlyLS, PrimalFeasible,
     &     iObj, signObj, scaleObj,
     &     n, nb, nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &     itn, nMajor, QNskips, maxVi, maxViRel, supVi,
     &     step, dxNorm, xNorm, fMerit, gMerit,
     &     sInfE, sInfEQP, sInfE1, wtInf,
     &     bl, bu, dx, dyCon,
     &     neJ, nlocJ, locJ, indJ, Jcol,
     &     neH, nlocH, locH, indH, Hcol,
     &     negCon, nlocG, locG,
     &     fObj1, fCon1, gCon1, gObj1,
     &     fObj2, fCon2, gCon2, gObj2, Fx,
     &     x, x1, x2, xQP, piQP,
     &     yCon, yCon1, yCon2, xPen,
     &     y, y1, y2, cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     funwrapper, fgcon, fgobj, userHv
      logical
     &     Elastic, FonlyLS, PrimalFeasible
      integer
     &     iExit, iObj, itn, lencu, lencw, leniu, leniw, lenru, lenrw,
     &     n, nb, neJ, negCon, neH, nlocG, nlocH, nlocJ, nnCon0, nnCon,
     &     nnJac, nnH, nnObj0, nnObj, nMajor, QNskips, indH(neH),
     &     indJ(neJ), locG(nlocG), locH(nlocH), locJ(nlocJ),
     &     iu(leniu), iw(leniw)
      double precision
     &     dxNorm, fMerit, gMerit, maxVi, maxViRel, supVi, sInfE,
     &     sInfEQP, sInfE1, scaleObj, signObj, step, wtInf, xNorm,
     &     bl(nb), bu(nb), Fx(nnCon0),
     &     fObj1, fCon1(nnCon0), gCon1(negCon), gObj1(nnObj0),
     &     fObj2, fCon2(nnCon0), gCon2(negCon), gObj2(nnObj0),
     &     dx(nb), dyCon(nnCon0), Hcol(neH), Jcol(neJ), piQP(nnCon0),
     &     yCon(nnCon0), yCon1(nnCon0), yCon2(nnCon0),
     &     x(nb), xQP(nb), x1(nb), x2(nb), xPen(nnCon0),
     &     y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s6lineSearch  finds the step length.
!
!     The set length defines points
!        x1 = x + step*dx  and  fmul1 = fmul + step*dfmul.
!     where, on entry,  x1 = xQP and  fmul1 = pi.
!
!     fCon1, JCon1, gObj1 and fmul1 are defined at the new point x1.
!     fCon2, JCon2, gObj2 and fmul2 are temporary work arrays.
!
!        iExit
!        -----
!         >0     Fatal error
!                = 22  Small step backing away from a violation limit.
!                = 63  User rejected the step
!          0     acceptable step found
!         -1     No step, no reason
!
!
!     20 Aug 2013: First version of s6lineSearch.
!     30 Aug 2013: Last update.
!     19 Oct 2014: Added user-defined Hessian.
!     17 Feb 2015: iExit = 22 added.
!     28 Apr 2015: s6stepLimits for linear constraints and bounds only.
!     ==================================================================
      character
     &     str*132
      logical
     &     StepFound, BackTrack, NoStep, Debug, NonlinearCon
      integer
     &     iMsg, inform, lprSrch, imaxVi,
     &     userBackTracks, violBackTracks, LSInfo
      double precision
     &     backFactor, eps0, fMerit1, gMerit1,
     &     stepLimit, stepMin, stepMax, wolfeG
!     ------------------------------------------------------------------
      integer            violLim,     userLim,     userRej
      parameter         (violLim = 1, userLim = 2, userRej = 3)
      parameter         (LSInfo  = 239) !  Line search result

      double precision   one
      parameter         (one     = 1.0d+0)
!     ------------------------------------------------------------------
      character          msg(4:8)*19
      data               msg /'max step too small.', ! inform 4
     &                        'step too small.    ', ! inform 5
     &                        'no minimizer.      ', ! inform 6
     &                        'too many functions.', ! inform 7
     &                        'uphill direction.  '/ ! inform 8
!     ------------------------------------------------------------------
      lprSrch      = iw( 82) ! line search debug starting itn

!     Constants

      eps0         = rw(  2) ! eps**(4/5)          IEEE DP  3.00e-13
      wolfeG       = rw( 84) ! Wolfe line search gradient tol

      NonlinearCon = nnCon .gt. 0

      Debug        = nMajor .ge. lprSrch ! turn on lines search o/p
      backFactor   = 0.1d+0              ! backtracking factor

!     ==================================================================
!     Find stepMin, stepMax and step,  the maximum, minimum and initial
!     values for the line search step length.
!     ==================================================================
      call s6stepLimits
     &   ( FonlyLS,
     &     nb, n, nnCon, nnObj, nMajor, QNskips,
     &     step, stepMin, stepLimit, stepMax, eps0, dxNorm, xNorm,
     &     bl, bu, x, dx, iw, leniw, rw, lenrw )

      userBackTracks = 0
      violBackTracks = 0

      StepFound   = .false.         ! Acceptable step found
      NoStep      = .false.         ! No step could be found.

!     ==================================================================
!     Repeat                               (until StepFound  or  NoStep)

  100    call s6search
     &      ( inform,
     &        funwrapper, fgcon, fgobj, userHv,
     &        Debug, Elastic, FonlyLS, PrimalFeasible,
     &        iObj, signObj, scaleObj,
     &        n, nb, nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &        itn, wolfeG, step, stepMin, stepMax,
     &        dxNorm, xNorm, fMerit, fMerit1, gMerit, gMerit1,
     &        sInfE, sInfEQP, sInfE1, wtInf,
     &        neJ   , nlocJ, locJ, indJ, Jcol,
     &        neH   , nlocH, locH, indH, Hcol,
     &        negCon, nlocG, locG,
     &        fObj1, fCon1, gCon1, gObj1, fObj2, fCon2, gCon2, gObj2,
     &        dx, dyCon, x, xQP, piQP, x1, x2,
     &        yCon, yCon1, yCon2, xPen,
     &        y, y1, y2, cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )

         if (inform .gt. 0) then
            iExit = inform      ! Fatal error
            go to 999
         end if

!     srSearch returns the following values:
!
!     inform    Result
!     ------    ------
!       >0      Fatal error
!        0      The user rejected the step.
!
!       -1      The search is successful and step < stepMax.
!       -2      The search is successful and step = stepMax.
!       -3      A better point was found but no sufficient decrease.
!               Most likely, the merit function is decreasing at the
!               boundary, but there could be too many function calls.
!
!       -4      stepMax < tolabs (too small to start a search).
!       -5      step   < stepMin  (srchq only)
!               Central differences may get a better direction.
!
!       -6      No useful step.
!               The interval of uncertainty is less than 2*tolabs.
!               The minimizer is very close to step = zero
!               or the gradients are not sufficiently accurate.
!       -7      Too many function calls with no better point.
!       -8      Input parameters were bad
!               (stepMax le toltny  or  oldg ge 0).

         StepFound  = inform .lt.    0 .and.
     &                inform .gt.   -4     ! step gives improvement
         BackTrack  = inform .eq.    0     ! User rejected the step
         NoStep     = inform .le.   -4     ! no step computed

         if (BackTrack) then
            userBackTracks = userBackTracks + 1
         end if

         if (StepFound .and. NonlinearCon) then
!           ------------------------------------------------------------
!           Acceptable step, but it may give a big constraint violation.
!           If so, the search is redone with a smaller stepMax.
!           ------------------------------------------------------------
            call s8Fx
     &         ( n, nnCon, nnJac, eps0,
     &           neJ, nlocJ, locJ, indJ, Jcol, fCon1, x1, Fx )
            call s2vmax
     &         ( n, nnCon, imaxVi, maxVi, bl, bu, Fx )
            maxViRel = maxVi / (one + xNorm)

            if (maxVi .gt. supVi) then
               violBackTracks = violBackTracks + 1
               BackTrack      = .true.
               StepFound      = .false.
            end if
         end if

         if (BackTrack) then
            stepMax    = backFactor*step
            step       = stepMax
            backFactor = backFactor*backFactor
         end if

!+    Until    (StepFound  or  NoStep)
!+    ------------------------------------------------------------------
      if (.not.(StepFound .or. NoStep)) go to 100

!     Record what happened and set the exit code.

      if (     userBackTracks .gt. 0) then
         iw(LSInfo) = userRej
      else if (violBackTracks .gt. 0) then
         iw(LSInfo) = violLim
      else if (step .ge. stepLimit) then
         iw(LSInfo) = userLim
      end if

      if (StepFound) then
         iExit = 0
      else if (violBackTracks .gt. 0 .and. inform .eq. -4) then
         iExit =  22           ! No step, caused by violation limit
      else if (userBackTracks .gt. 0) then
         iExit =  63           ! No step, user wants to stop
      else
         iExit = -1            ! No step, no reason
         iMsg  = -inform
         write(str, 1050) iMsg, msg(iMsg), nMajor
         call snPRNT( 23, str, iw, leniw )
      end if

  999 return

 1050 format(' Search exit', i3, ' -- ', a, ' Major itn =', i7)

      end ! subroutine s6lineSearch

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6search
     &   ( iExit,
     &     funwrapper, fgcon, fgobj, userHv,
     &     Debug, Elastic, FonlyLS, PrimalFeasible,
     &     iObj, signObj, scaleObj,
     &     n, nb, nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &     itn, wolfeG, step, stepMin, stepMax,
     &     dxNorm, xNorm, fMerit, fMerit1, gMerit, gMerit1,
     &     sInfE, sInfEQP, sInfE1, wtInf,
     &     neJ   , nlocJ, locJ, indJ, Jcol,
     &     neH   , nlocH, locH, indH, Hcol,
     &     negCon, nlocG, locG,
     &     fObj1, fCon1, gCon1, gObj1, fObj2, fCon2, gCon2, gObj2,
     &     dx, dyCon,  x, xQP, piQP, x1, x2,
     &     yCon, yCon1, yCon2, xPen,
     &     y, y1, y2, cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     funwrapper, fgcon, fgobj, userHv
      logical
     &     Debug, Elastic, FonlyLS, PrimalFeasible
      integer
     &     iExit, iObj, itn, lencu, lencw, leniu, leniw, lenru, lenrw,
     &     n, nb, negCon, neH, neJ, nlocG, nlocH, nlocJ, nnCon0, nnCon,
     &     nnJac, nnH, nnObj0, nnObj, indH(neH), indJ(neJ), locG(nlocG),
     &     locH(nlocH), locJ(nlocJ), iu(leniu), iw(leniw)
      double precision
     &     wolfeG, fMerit, fMerit1, gMerit, gMerit1, dxNorm,
     &     sInfE, sInfEQP, sInfE1, scaleObj, signObj,
     &     step, stepMin, stepMax, wtInf, xNorm,
     &     fObj1, fCon1(nnCon0), gCon1(negCon), gObj1(nnObj0),
     &     fObj2, fCon2(nnCon0), gCon2(negCon), gObj2(nnObj0),
     &     dx(nb), dyCon(nnCon0), Hcol(neH), Jcol(neJ), piQP(nnCon0),
     &     yCon(nnCon0), yCon1(nnCon0), yCon2(nnCon0),
     &     x(nb), xQP(nb), x1(nb), x2(nb), xPen(nnCon0),
     &     y(nb), y1(nb), y2(nb), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

*     ==================================================================
*     s6search  finds a step along the search direction  p,  such that
*     the function  fMerit  is sufficiently reduced, i.e.,
*               fMerit(x + step*p)  <  fMerit(x).
*
*     On entry
!     --------
*     step     is an initial estimate of the step length.
*     FonlyLS  is true if function-only search should be used.
*     fMerit   is the value of fMerit at the base point x.
*     gMerit   is the directional derivative of the merit function.
!     sInfE    is the sum of elastic infeasibilities at x.
!     sInfEQP  is the sum of elastic infeasibilities at xQP.
*     x(nb)    is the vector of variables   at the base point.
*     yCon     is the vector of multipliers at the base point.
*     dx(nb)   is the search direction.
*     piQP     contains the  QP multipliers.
*     dx(nb)   The search direction, xQP - x1.
*
*     On exit
!     --------
*     step     is the final step.
*     fMerit1  is the final value of fMerit.
*     x1       is the final value of the variables.
*     fCon1,gCon1,gObj1,yCon1 are defined at the new point x1.
*
*     fCon2,gCon2,gObj2,yCon2 and x2 are work arrays.
*
*     iExit    Result
*     -----    ------
*      >0      Fatal error
*       0      Repeat the search with smaller stpmax.
*      -1      The search is successful and step < stpmax.
*      -2      The search is successful and step = stpmax.
*      -3      A better point was found but too many functions
*              were needed (not sufficient decrease).
*      -4      stpmax < tolabs (too small to do a search).
*      -5      step   < stepMin (srchq only -- maybe want to switch
*              to central differences to get a better direction).
*      -6      No useful step.
*              The interval of uncertainty is less than 2*tolabs.
*              The minimizer is very close to step = zero
*              or the gradients are not sufficiently accurate.
*      -7      Too many function calls.
*      -8      Bad input parameters
*              (stpmax le toltny  or  oldg ge 0).
*
*     30 Dec 1991: First version based on NPSOL 4.6 routine npsrch.
*     28 Sep 1993: Allow functions to say "undefined at this point".
*                  (Back up and try again.)
*     18 Feb 1994: Back up in a similar way if vimax increases a lot.
*                  Deleted after first limited-memory version.
*     29 Dec 1994: Merit function calculations included explicitly.
*     06 Apr 1996: Special coding for the unit step.
*                  On entry, x2 is the QP solution.
*     18 Oct 1996: First Min-sum version.
*     17 Jul 1997: First thread-safe version.
*     11 Oct 1998: Facility to combine funobj and funcon added.
*     11 Jun 2000: Tolerances computed in a subroutine.
*     24 Oct 2000: Updated for SNOPT 6.1
*     04 Aug 2003: snEXIT and snPRNT adopted.
*     15 Jun 2008: Call-status implemented correctly.
*     15 Nov 2010: Call-status removed from argument list.
*     19 Oct 2014: Added user-defined Hessian.
*     25 May 2015: Defn of ftry rearranged to reduce cancellation error.
*     ==================================================================
      character
     &     str*80
      external
     &     ddot
      logical
     &     UsingG, Done, First, Improved, NonlinearCon, NonlinearObj,
     &     Brakted, Cramped, Extrap, Moved, vset, wset
      integer
     &     inform, iPrint, jObj, maxf, modefg, nLin, nnJ1,
     &     nout, nsamea, nsameb, numf
      double precision
     &     eps0,  aLow, bUpp, bigFx, dsInfE, epsaf, fMerit2, g0,
     &     gMerit2, oldf, oldg,
     &     fa, factor, ftry, fv, fw, gtry, gw,
     &     sbest, fbest, gbest, sInfE2, stpmax,
     &     targetg, tolAbs, tolRel, tolTiny,
     &     tolMax, xtry, xv, xw, ddot
*     ------------------------------------------------------------------
      integer            Normal
      parameter         (Normal = 0)
      double precision   zero,          half,          one
      parameter         (zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0)
      double precision   wolfeF
      parameter         (wolfeF = 1.0d-4)
*     ------------------------------------------------------------------
      iPrint    = iw( 12) ! Print file

      eps0      = rw(  2) ! eps**(4/5)
      bigFx     = rw( 71) ! unbounded objective.

      iExit     = 0

*     ------------------------------------------------------------------
*     Set the input parameters for srchc or srchq.
*
*     stepMin is used by srchq.  If  step  would be less than  stepMin,
*             the search will be terminated early.
*             If p was found using forward or backward differences,
*             stepMin  should be positive (and related to the difference
*             interval used).
*             If p was found using central differences (lvlDif = 2)
*             stepMin  should be zero.
*
*     epsaf   is the absolute function precision. If f(x1) and f(x2) are
*             as close as  epsaf,  we cannot safely conclude from the
*             function values alone which of x1 or x2 is a better point.
*
*     tolAbs  is an estimate of the absolute spacing between points
*             along  p.  This step should produce a perturbation
*             of  epsaf  in the merit function.
*
*     tolRel  is an estimate of the relative spacing between points
*             along  p.
*
*     tolTiny is the minimum allowable absolute spacing between points
*             along  p.
*     ------------------------------------------------------------------
      stpmax = stepMax
      UsingG = .not. FonlyLS

      if (FonlyLS) then
         maxf   = 15
         modefg =  0
      else
         maxf   = 10
         modefg =  2
      end if

      nlin         = n - nnJac
      nnJ1         = nnJac + 1

      nout         = iPrint
      jObj         = n + iObj
      NonlinearCon = nnCon .gt. 0
      NonlinearObj = nnObj .gt. 0

*     Define the line search tolerances.

      call s6tols
     &   ( nb, epsaf, stepMax, tolAbs, tolRel, tolTiny,
     &     dxNorm, xNorm, fMerit, dx, x, rw, lenrw )

      oldf    = fMerit
      oldg    = gMerit
      fMerit1 = fMerit
      gMerit1 = gMerit
      call dcopy ( nb   ,    x, 1,    x1, 1 )

      if (NonlinearCon)
     &call dcopy ( nnCon, yCon, 1, yCon1, 1 )

      if (Elastic) then
         sInfE1  = sInfE
         dsInfE  = sInfEQP - sInfE
      end if

      fObj2   = zero            ! keeps ftnchek quiet
      fMerit2 = zero            ! keeps ftnchek quiet
      ftry    = zero            ! keeps ftnchek quiet
      gtry    = zero            ! keeps ftnchek quiet
      fv      = zero            ! keeps ftnchek quiet
      fw      = zero            ! keeps ftnchek quiet

      First   = .true.
      sbest   = zero
      fbest   = zero
      gbest   = (one     - wolfeF)*oldg
      targetg = (wolfeF  - wolfeG)*oldg
      g0      = gbest

      if (Debug) write(nout, 1000) itn, dxNorm

*     ------------------------------------------------------------------
*     Commence main loop, entering srchc or srchq two or more times.
*     First = true for the First entry,  false for subsequent entries.
*     Done  = true indicates termination, in which case inform gives
*     the result of the search (with inform = iExit as above).
*     ------------------------------------------------------------------
*+    repeat
  200    if (UsingG) then
            call srchc
     &         ( inform , First  , Debug  , Done  , Improved,
     &           maxf   , numf   , nout   ,
     &           stpmax ,          epsaf  ,
     &           g0     , targetg, ftry   , gtry  ,
     &           tolAbs , tolRel , tolTiny,
     &           step   , sbest  , fbest  , gbest ,
     &           Brakted, Cramped, Extrap , Moved , wset  ,
     &           nsamea , nsameb ,
     &           aLow   , bUpp   , factor ,
     &           xtry   , xw     , fw     , gw    , tolMax )
         else
            call srchq
     &         ( inform , First  , Debug  , Done , Improved,
     &           maxf   , numf   , nout   ,
     &           stpmax , stepMin, epsaf  ,
     &           g0     , targetg, ftry   ,
     &           tolAbs , tolRel , tolTiny,
     &           step   , sbest  , fbest  ,
     &           Brakted, Cramped, Extrap , Moved , vset  , wset  ,
     &           nsamea , nsameb ,
     &           aLow   , bUpp   , fa     , factor,
     &           xtry   , xw     , fw     , xv    , fv    , tolMax)
         end if

         if (Improved) then
            fMerit1 = fMerit2

            if (NonlinearCon) then
               call dcopy ( nnCon, fCon2, 1, fCon1, 1 )
               if (UsingG) then
                  call dcopy ( negCon, gCon2, 1, gCon1, 1 )
               end if
            end if

            fObj1 = fObj2
            if (NonlinearObj) then
               if (UsingG) then
                  call dcopy ( nnObj , gObj2, 1, gObj1, 1 )
               end if

*              Terminate if the objective is unbounded below in the
*              feasible region.

               if (signObj*fObj1 .lt. -bigFx .and. PrimalFeasible) then
                  iExit = 21
                  go to 900
               end if
            end if
         end if ! Improved

*        ---------------------------------------------------------------
*           Done = false  First time through.
*        If Done = false, the functions must be computed for the next
*                  entry to srchc or srchq.
*        If Done = true,  this is the last time through and inform ge 1.
*        ---------------------------------------------------------------
         if (.not. Done) then

            if (step .eq. one) then
               call dcopy ( nb, xQP, 1, x2, 1 )

               if (Elastic) then
                  sInfE2 = sInfEQP
               end if

               if (NonlinearCon) then
                  call dcopy ( nnCon, piQP, 1, yCon2, 1 )
               end if
            else
               call dcopy ( nb,       x1 , 1, x2, 1 )
               call daxpy ( nb, step, dx , 1, x2, 1 )

               if (Elastic) then
                  sInfE2 = sInfE1 + step*dsInfE
               end if
            end if

            if (nnH .gt. 0) then
               call funwrapper
     &            ( inform,
     &              modefg, NonlinearCon, NonlinearObj,
     &              n, negCon, nnCon0, nnCon,
     &              nnJac, nnH, nnObj0, nnObj,
     &              fgcon, fgobj, userHv,
     &              x2, yCon,
     &              neJ, nlocJ, locJ, indJ,
     &              neH, nlocH, locH, indH, Hcol,
     &              fCon2, fObj2, gCon2, gObj2,
     &              cu, lencu, iu, leniu, ru, lenru,
     &              cw, lencw, iw, leniw, rw, lenrw )
               if (inform .ne. 0) then
                  if (inform .gt. 0) then
                     iExit  = inform
                  else
                     inform = 0 ! Redo the search
                  end if
                  go to 900
               end if
            end if

            if (iObj .eq. 0) then
               fMerit2 = zero
            else
               fMerit2 = signObj*x2(jObj)*scaleObj
            end if

            if (NonlinearObj) then
               fMerit2 = fMerit2 + signObj*fObj2
            end if

            if (Elastic) then
               fMerit2 = fMerit2 + wtInf*sInfE2
            end if

            if (NonlinearCon) then
*              ---------------------------------------------------------
*              Compute w and y, the constraint violations and the
*              negative of the gradient of the merit function with
*              respect to the nonlinear slacks. These quantities define
*              the directional derivative of the merit function.
*              Finally, add the nonlinear terms to the merit function.
*              ---------------------------------------------------------
               if (step .le. one) then
                  call dcopy ( nnCon,       yCon1, 1, yCon2, 1 )
                  call daxpy ( nnCon, step, dyCon, 1, yCon2, 1 )
               end if

*              Compute the constraint violations and aux. multipliers:
*              y1 = fCon  + A(linear) x - nonlinear slacks.
*              y  = yCon2 - xPen*viol.

               call dcopy ( nnCon,           fCon2, 1, y1, 1 )
               call daxpy ( nnCon, (-one), x2(n+1), 1, y1, 1 )
               if (nlin .gt. 0) then
                  call s2Aprod
     &               ( Normal, eps0,
     &                 neJ, nlin+1, locJ(nnJ1), indJ, Jcol,
     &                 one, x2(nnJ1), nlin, one, y1, nnCon )
               end if

               call dcopy ( nnCon,   y1, 1, y, 1 )
               call ddscl ( nnCon, xPen, 1, y, 1 )
               fMerit2 = fMerit2 -      ddot( nnCon, yCon2, 1, y1, 1 )
     &                           + half*ddot( nnCon, y    , 1, y1, 1 )

               if (UsingG) then
                  call daxpy ( nnCon, (-one), yCon2, 1, y, 1 )
                  call dscal ( nnCon, (-one), y    , 1 )
               end if
            end if ! NonlinearCon

            ftry  = fMerit2 - (oldf + wolfeF*oldg*step)

            if (UsingG) then
*              ---------------------------------------------------------
*              A gradient search is requested.
*              Compute the directional derivative gtry.
*              ---------------------------------------------------------
               if (iObj .eq. 0) then
                  gMerit2 = zero
               else
                  gMerit2 = signObj*dx(jObj)*scaleObj
               end if

               if (NonlinearObj) then
                  gMerit2 = gMerit2
     &                          + signObj*ddot( nnObj, gObj2, 1, dx, 1 )
               end if

               if (Elastic) then
                  gMerit2 = gMerit2 + wtInf*dsInfE
               end if

               if (NonlinearCon) then
*                 ------------------------------------------------------
*                 Form  J*dx (including linear columns).
*                 Set  y2 = J*dx + A(linear) dx - dx(slacks).
*                 ------------------------------------------------------
                  call s8Gprod
     &               ( Normal, eps0,
     &                 neJ   , nlocJ, locJ, indJ,
     &                 negCon, nlocG, locG, gCon2,
     &                 one, dx, nnJac, zero, y2, nnCon )

                  if (nlin .gt. 0) then
                     call s2Aprod
     &                  ( Normal, eps0,
     &                    neJ, nlin+1, locJ(nnJ1), indJ, Jcol,
     &                    one, dx(nnJ1), nlin, one, y2, nnCon )
                  end if
                  call daxpy ( nnCon, (-one), dx(n+1), 1, y2, 1 )

                  gMerit2 = gMerit2 - ddot  ( nnCon, y , 1, y2   , 1 )
     &                              - ddot  ( nnCon, y1, 1, dyCon, 1 )

               end if ! NonlinearCon

               gtry   = gMerit2 - wolfeF*oldg

            end if ! UsingG
         end if ! not Done
*+    until (      Done)
      if    (.not. Done) go to 200

*     ==================================================================
*     The search is Done.
*     Finish with  x1 = the best point found so far.
*     ==================================================================
      step = sbest

      if (Improved) then
         call dcopy ( nb   ,    x2, 1,    x1, 1 )
         if (NonlinearCon) then
            call dcopy ( nnCon, yCon2, 1, yCon1, 1 )
         end if
         if (Elastic) then
            sInfE1 = sInfE2
         end if

      else if (step .gt. zero) then
         call daxpy ( nb   ,  step, dx   , 1, x1   , 1 )
         if (NonlinearCon) then
            call daxpy ( nnCon,  step, dyCon, 1, yCon1, 1 )
         end if
         if (Elastic) then
            sInfE1 = sInfE1 + step*dsInfE
         end if
      end if

*     ------------------------------------------------------------------
*     Print any warning messages.
*     ------------------------------------------------------------------
      if (inform .eq. 7) then
         write(str, 1700) numf
         call snPRNT( 23, str, iw, leniw )

      else if (inform .eq. 8) then
         if (oldg .ge. zero) then
            write(str, 1600) oldg
            call snPRNT( 21, str, iw, leniw )
         else
            write(str, 1800) stepMax, dxNorm, oldg, numf
            call snPRNT( 21, str, iw, leniw )
         end if
      end if

      iExit = -inform

  900 return

 1000 format(// ' --------------------------------------------'
     &       /  ' Output from s6search following iteration', i7,
     &      5x, ' Norm p =', 1p, e11.2 )
 1600 format(   ' XXX  The search direction is uphill.  gMerit  =',
     &          1p, e9.1)
 1700 format(   ' XXX  The line search has evaluated the functions', i5,
     &          '  times')
 1800 format(   ' stepMax =', 1p, e11.2, '   dxNorm =', e11.2,
     &          '   gMerit =',     e11.2, '    numf  =', i3)

      end ! subroutine s6search

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6stepLimits
     &   ( FonlyLS,
     &     nb, n, nnCon, nnObj, majors, QNskips,
     &     step, stepMin, stepLimit, stepMax, tolz, dxNorm, xNorm,
     &     bl, bu, x, dx, iw, leniw, rw, lenrw )

      implicit
     &     none
      logical
     &     FonlyLS
      integer
     &     nb, n, nnCon, nnObj, majors, QNskips, leniw, lenrw,
     &     iw(leniw)
      double precision
     &     step, stepMin, stepLimit, stepMax, tolz, dxNorm, xNorm,
     &     bl(nb), bu(nb), x(nb), dx(nb), rw(lenrw)

!     ==================================================================
!     s6stepLimits  finds the maximum, minimum and initial values for
!     the line search step.
!
!     If the problems has nonlinear constraints and the QN update is
!     has been successful in prior iterations,  then the maximum step
!     stepMax is one.  Otherwise, the maximum step is the largest step
!     such that x + step*dx reaches one of the linear constraint bounds.
!
!     All step sizes are subject to the user-specified limit  stepLimit.
!
!     04 Dec 1992: First version (s6step) based on npsol routine npalf.
!     31 Mar 2000: Updated for SNOPT 6.1.
!     26 Apr 2015: Renamed from s8stepLimits.
!     ==================================================================
      external
     &     ddiv
      logical
     &     Switch, OverFlow
      integer
     &     HDInfo,  lvlDif, gotFD
      double precision
     &     bigdx, fdint1, stepQP, dxlimit
      double precision
     &     ddiv
!     ------------------------------------------------------------------
      integer            Unit
      parameter         (Unit = 2)
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      bigdx     = rw( 72) ! unbounded step.
      fdint1    = rw( 76) ! (1) forwrd diff. interval
      dxlimit   = rw( 80) ! Step limit

      lvlDif    = iw(181) ! 1(2) for forwd (cntrl) diffs
      gotFD     = iw(183) ! > 0 => some differences needed
      HDInfo    = iw(243) ! infoTag(7): Approximate Hessian type

      OverFlow    = .false.

!     ==================================================================
!     Switch   indicates if there is an option to switch to
!              central differences to get a better search direction.
!     stepQP   is the step predicted by the QP subproblem (usually 1).
!     stepMax  is the largest feasible steplength subject to a
!              user-defined limit, bigdx, on the change in  x.
!     step     is initialized subject to a user-defined limit, dxlim.
!     ==================================================================
      if (nnCon .eq. 0  .and.  nnObj .eq. 0) then  ! LP !!
!        Linear program
         step      = one
         stepMin   = one
         stepLimit = one
         stepMax   = one
      else
!        Nolinear program
         Switch = gotFD .gt. 0  .and.  lvlDif .ne. 2

         stepMin = zero
         if (FonlyLS  .and.  Switch) then
            stepMin = fdint1*(one + xNorm)/dxNorm
         end if

         stepLimit = ddiv( (one+xNorm)*dxlimit, dxNorm, OverFlow )
         if (majors .le. 1)
     &   stepLimit = min (stepLimit, ddiv( one, dxNorm, OverFlow ))

         stepQP    = one

         if (nnCon .gt. 0  .and.  (QNskips .eq. 0  .or.
     &                             HDInfo  .ne. Unit)) then
            stepMax = one
         else
            stepMax = ddiv  ( bigdx, dxNorm, OverFlow )
            if (nnCon .gt. 0)
     &      stepMax = min( stepMax, stepLimit )
            step    = stepMax

            call s6LCstepMax
     &         ( nb, n, nnCon, step, stepMax, dxNorm,
     &           bl, bu, x, dx, rw, lenrw )
            step   = max( step, stepQP )
            if (step .lt. stepQP + tolz) step = stepQP
            stepMax = step
         end if

         stepMax   = min ( stepLimit, stepMax )
         step      = min ( stepLimit, one     )
      end if

      end ! subroutine s6stepLimits

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6tols
     &   ( nb, epsaf, stepMax, tolAbs, tolRel, tolTiny,
     &     dxNorm, xNorm, f, dx, x, rw, lenrw )

      implicit
     &     none
      integer
     &     nb, lenrw
      double precision
     &     epsaf, f, stepMax, tolAbs, tolRel, tolTiny, dxNorm, xNorm,
     &     dx(nb), x(nb), rw(lenrw)

*     ==================================================================
*     s6tols defines various tolerances for the line search.
*
*     11 Jun 2000: First   version of s6tols.
*     11 Jun 2000: Current version of s6tols.
*     ==================================================================
      integer
     &     j
      double precision
     &     eps, eps0, epsrf, s, q, t, tolax, tolrx
*     ------------------------------------------------------------------
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
*     ------------------------------------------------------------------
      eps       = rw(  1) ! machine precision.  IEEE DP  2.22e-16
      eps0      = rw(  2) ! eps**(4/5)
      epsrf     = rw( 73) ! relative function precision.
*     ------------------------------------------------------------------
      epsaf  = max( epsrf, eps )*(one + abs(f))
      tolax  = eps0
      tolrx  = eps0
      t      = xNorm*tolrx  +  tolax
      if (t .lt. dxNorm*stepMax) then
         tolAbs = t/dxNorm
      else
         tolAbs = stepMax
      end if
      tolRel = max( tolrx, eps )

      t      = zero
      do j = 1, nb
         s     = abs(dx(j))
         q     = abs( x(j))*tolrx + tolax
         if (s .gt. q*t) t = s / q
      end do

      if (t*tolAbs .gt. one) then
         tolTiny = one / t
      else
         tolTiny = tolAbs
      end if

      end ! subroutine s6tols

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s6Userf
     &   ( mode, n, x, xPert, damper, userfg,
     &     needF, nF  , FPert,
     &     NonlinearObj, lenG, GPert,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     iw, leniw )

      implicit
     &     none
      external
     &     userfg
      integer
     &     mode, n, needF, NonlinearObj, nF, lenG, lencu, leniu, lenru,
     &     leniw, iu(leniu), iw(leniw)
      double precision
     &     damper, FPert(nF), GPert(lenG), x(n), xPert(n), ru(lenru)
      character
     &     cu(lencu)*8

*     ==================================================================
*     s6Userf  attempts to compute the snoptA problem functions at a
*     point xPert.  If this computation is undefined, then the
*     evaluation takes place at a point on the ray joining xPert and
*     a point  x  at which the functions are known to be well-defined.
*
*     On entry:
*        xPert  is the point at which the functions are required.
*
*        x      is a point at which the problem functions have been
*               computed successfully.
*
*     On exit:
*        mode   is nonnegative if the problem function were evaluated
*               successfully.   Otherwise, five attempts to evaluate
*               the functions at points closer to x were unsuccessful.
*
*        If mode ge 0,  then the output values of damper and xPert are
*        defined as follows:
*
*        damper (0 lt damper le 1) is 1 if the functions were
*               evaluated successfuly at the input value of xPert.
*               Otherwise the problem functions were evaluated
*               successfuly at xPert + damper*(xPert - x).
*
*        xPert  is  xPert(in) + damper*(xPert(in) - x),
*
*     26 Oct 2002: First version.
*     22 Apr 2007: damper added as an output argument.
*     15 Jun 2008: Call-status implemented correctly.
*     ==================================================================
      integer
     &     j, Status, tries
*     ------------------------------------------------------------------
      double precision   one,            ten
      parameter         (one   = 1.0d+0, ten  =10.0d+0)
*     ------------------------------------------------------------------

      mode = 0

      ! Determine the status of this call.

      call s8callStatus( Status, iw, leniw )

      damper = one
      tries  = 0
*     ==================================================================
*     Repeat                       (until problem functions are defined)

  100    tries  = tries + 1
         call userfg
     &      ( Status, n, xPert,
     &        needF, nF, FPert,
     &        NonlinearObj, lenG, GPert,
     &        cu, lencu, iu, leniu, ru, lenru )

         mode   = Status
         Status = 0

         if (mode .eq. -1) then
            damper = damper*damper/ten
            do   j = 1, n
               xPert(j) = damper*xPert(j) + (one - damper)*x(j)
            end do
         end if

*+    until (.not. (mode .eq. -1  .and.  tries .lt. 5))
      if (          mode .eq. -1  .and.  tries .lt. 5 ) go to 100
*     ==================================================================

      end ! subroutine s6Userf

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine srchc
     &   ( iExit  , First  , Debug  , Done  , Improved,
     &     maxf   , numf   , nout   ,
     &     alfmax ,          epsaf  ,
     &     g0     , targetg, ftry   , gtry  ,
     &     tolAbs , tolRel , tolTiny,
     &     alfa   , alfbst , fbest  , gbest ,
     &     Brakted, Cramped, Extrap , Moved , wset  ,
     &     nsamea , nsameb ,
     &     a      , b      , factor ,
     &     xtry   , xw     , fw     , gw    , tolMax )

      implicit
     &     none
      logical
     &     First, Debug, Done, Improved, Brakted, Cramped, Extrap,
     &     Moved, wset
      integer
     &     iExit, maxf, numf, nout, nsamea, nsameb
      double precision
     &     alfmax, epsaf, g0, targetg, ftry, gtry, tolAbs, tolRel,
     &     tolTiny, alfa, alfbst, fbest, gbest, a, b, factor, xtry,
     &     xw, fw, gw, tolMax

*     ==================================================================
*     srchc  finds a sequence of improving estimates of a minimizer of
*     the univariate function f(alpha) in the interval (0,alfmax].
*     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
*     srchc requires both  f(alpha)  and  f'(alpha)  to be evaluated at
*     points in the interval.  Estimates of the minimizer are computed
*     using safeguarded cubic interpolation.
*
*     Reverse communication is used to allow the calling program to
*     evaluate f and f'.  Some of the parameters must be set or tested
*     by the calling program.  The remainder would ordinarily be local
*     variables.
*
*     Input parameters (relevant to the calling program)
*     --------------------------------------------------
*
*     First         must be true on the First entry. It is subsequently
*                   altered by srchc.
*
*     Debug         specifies whether detailed output is wanted.
*
*     maxf          is an upper limit on the number of times srchc is
*                   to be entered consecutively with Done = false
*                   (following an initial entry with First = true).
*
*     alfa          is the first estimate of a minimizer.  alfa is
*                   subsequently altered by srchc (see below).
*
*     alfmax        is the upper limit of the interval to be searched.
*
*     epsaf         is an estimate of the absolute precision in the
*                   computed value of f(0).
*
*     ftry, gtry    are the values of f, f'  at the new point
*                   alfa = alfbst + xtry.
*
*     g0            is the value of f'(0).  g0 must be negative.
*
*     tolAbs,tolRel define a function tol(alfa) = tolRel*alfa + tolAbs
*                   such that if f has already been evaluated at alfa,
*                   it will not be evaluated closer than tol(alfa).
*                   These values may be reduced by srchc.
*
*     targetg       is the target value of abs(f'(alfa)). The search
*                   is terminated when
*                    abs(f'(alfa)) le targetg and f(alfa) lt 0.
*
*     tolTiny       is the smallest value that tolAbs is allowed to be
*                   reduced to.
*
*     Output parameters (relevant to the calling program)
*     ---------------------------------------------------
*
*     Improved      is true if the previous alfa was the best point so
*                   far.  Any related quantities should be saved by the
*                   calling program (e.g., gradient arrays) before
*                   paying attention to the variable Done.
*
*     Done = false  means the calling program should evaluate
*                      ftry = f(alfa),  gtry = f'(alfa)
*                   for the new trial alfa, and re-enter srchc.
*
*     Done = true   means that no new alfa was calculated.  The value
*                   of iExit gives the result of the search as follows
*
*                   iExit = 1 means the search has terminated
*                             successfully with alfbst < alfmax.
*
*                   iExit = 2 means the search has terminated
*                             successfully with alfbst = alfmax.
*
*                   iExit = 3 means that the search failed to find a
*                             point of sufficient decrease.
*                             The function is either decreasing at
*                             alfmax or maxf function evaluations
*                             have been exceeded.
*
*                   iExit = 4 means alfmax is so small that a search
*                             should not have been attempted.
*
*                   iExit = 5 is never set by srchc.
*
*                   iExit = 6 means the search has failed to find a
*                             useful step.  The interval of uncertainty
*                             is [0,b] with b < 2*tolAbs. A minimizer
*                             lies very close to alfa = 0, or f'(0) is
*                             not sufficiently accurate.
*
*                   iExit = 7 if no better point could be found after
*                             maxf  function calls.
*
*                   iExit = 8 means the input parameters were bad.
*                             alfmax le tolTiny  or g0 ge zero.
*                             No function evaluations were made.
*
*     numf          counts the number of times srchc has been entered
*                   consecutively with Done = false (i.e., with a new
*                   function value ftry).
*
*     alfa          is the point at which the next function ftry and
*                   derivative gtry must be computed.
*
*     alfbst        should be accepted by the calling program as the
*                   approximate minimizer, whenever srchc returns
*                   iExit = 1 or 2 (and possibly 3).
*
*     fbest, gbest  will be the corresponding values of f, f'.
*
*
*     The following parameters retain information between entries
*     -----------------------------------------------------------
*
*     Brakted       is false if f and f' have not been evaluated at
*                   the far end of the interval of uncertainty.  In this
*                   case, the point b will be at alfmax + tol(alfmax).
*
*     Cramped       is true if alfmax is very small (le tolAbs).  If the
*                   search fails, this indicates that a zero step should
*                   be taken.
*
*     Extrap        is true if xw lies outside the interval of
*                   uncertainty.  In this case, extra safeguards are
*                   applied to allow for instability in the polynomial
*                   fit.
*
*     Moved         is true if a better point has been found, i.e.,
*                   alfbst gt 0.
*
*     wset          records whether a second-best point has been
*                   determined it will always be true when convergence
*                   is tested.
*
*     nsamea        is the number of consecutive times that the
*                   left-hand end point of the interval of uncertainty
*                   has remained the same.
*
*     nsameb        similarly for the right-hand end.
*
*     a, b, alfbst  define the current interval of uncertainty.
*                   A minimizer lies somewhere in the interval
*                   [alfbst + a, alfbst + b].
*
*     alfbst        is the best point so far.  It is always at one end
*                   of the interval of uncertainty.  hence we have
*                   either  a lt 0,  b = 0  or  a = 0,  b gt 0.
*
*     fbest, gbest  are the values of f, f' at the point alfbst.
*
*     factor        controls the rate at which extrapolated estimates
*                   of alfa may expand into the interval of uncertainty.
*                   factor is not used if a minimizer has been bracketed
*                   (i.e., when the variable Brakted is true).
*
*     fw, gw        are the values of f, f' at the point alfbst + xw.
*                   they are not defined until wset is true.
*
*     xtry          is the trial point in the shifted interval (a, b).
*
*     xw            is such that  alfbst + xw  is the second-best point.
*                   it is not defined until  wset  is true.
*                   in some cases,  xw  will replace a previous  xw
*                   that has a lower function but has just been excluded
*                   from the interval of uncertainty.
*
*
*     Systems Optimization Laboratory, Stanford University, California.
*     Original version February 1982.  Rev. May 1983.
*     Original f77 version 22-August-1985.
*     14 Sep 1992: Introduced QuitI, QuitF, etc.
*     22 Nov 1995: Altered criterion for reducing the step below tolAbs.
*     17 Jul 1997: Removed saved variables for thread-safe version.
*     19 Apr 2000: QuitF only allowed after a move.
*     19 Apr 2000: Current version.
*     ==================================================================
      logical
     &     BadFun, closef, Found, QuitF, QuitI, Fitok, setxw
      double precision
     &     absr, artifa, artifb, daux, dtry,
     &     q, r, s, scale, tol, truea, trueb, xmidpt
*     ------------------------------------------------------------------
      double precision   zero,          point1,          half
      parameter         (zero = 0.0d+0, point1 = 0.1d+0, half = 0.5d+0)
      double precision   one,           three,           five
      parameter         (one  = 1.0d+0, three  = 3.0d+0, five = 5.0d+0)
      double precision   ten,           eleven
      parameter         (ten  = 1.0d+1, eleven = 1.1d+1               )
*     ------------------------------------------------------------------
*     Local variables
*     ===============
*
*     closef     is true if the new function ftry is within epsaf of
*                fbest (up or down).
*
*     Found      is true if the sufficient decrease conditions hold at
*                alfbst.
*
*     QuitF      is true when  maxf  function calls have been made.
*
*     QuitI      is true when the interval of uncertainty is less than
*                2*tol.
*  ---------------------------------------------------------------------

      BadFun   = .false.
      QuitF    = .false.
      QuitI    = .false.
      Improved = .false.

      if (First) then
*        ---------------------------------------------------------------
*        First entry.  Initialize various quantities, check input data
*        and prepare to evaluate the function at the initial alfa.
*        ---------------------------------------------------------------
         First  = .false.
         numf   = 0
         alfbst = zero
         BadFun = alfmax .le. tolTiny  .or.  g0 .ge. zero
         Done   = BadFun
         Moved  = .false.

         if (.not. Done) then
            Brakted = .false.
            Cramped = alfmax .le. tolAbs
            Extrap  = .false.
            wset    = .false.
            nsamea  = 0
            nsameb  = 0

            tolMax = tolAbs + tolRel*alfmax
            a      = zero
            b      = alfmax + tolMax
            factor = five
            tol    = tolAbs
            xtry   = alfa
            if (Debug) then
               write(nout, 1000)
     &              g0     , tolAbs, alfmax,
     &              targetg, tolRel, epsaf , Cramped
            end if
         end if
      else
*        ---------------------------------------------------------------
*        Subsequent entries. The function has just been evaluated at
*        alfa = alfbst + xtry,  giving ftry and gtry.
*        ---------------------------------------------------------------
         if (Debug) write(nout, 1100) alfa, ftry, gtry

         numf   = numf   + 1
         nsamea = nsamea + 1
         nsameb = nsameb + 1

         if (.not. Brakted) then
            tolMax = tolAbs + tolRel*alfmax
            b      = alfmax - alfbst + tolMax
         end if

*        See if the new step is better.  If alfa is large enough that
*        ftry can be distinguished numerically from zero,  the function
*        is required to be sufficiently negative.

         closef = abs( ftry - fbest ) .le.  epsaf
         if (closef) then
            Improved =  abs( gtry ) .le. abs( gbest )
         else
            Improved = ftry .lt. fbest
         end if

         if (Improved) then

*           We seem to have an improvement.  The new point becomes the
*           origin and other points are shifted accordingly.

            fw     = fbest
            fbest  = ftry
            gw     = gbest
            gbest  = gtry
            alfbst = alfa
            Moved  = .true.

            a      = a    - xtry
            b      = b    - xtry
            xw     = zero - xtry
            wset   = .true.
            Extrap =       (xw .lt. zero  .and.  gbest .lt. zero)
     &               .or.  (xw .gt. zero  .and.  gbest .gt. zero)

*           Decrease the length of the interval of uncertainty.

            if (gtry .le. zero) then
               a       = zero
               nsamea  = 0
            else
               b       = zero
               nsameb  = 0
               Brakted = .true.
            end if
         else

*           The new function value is not better than the best point so
*           far.  The origin remains unchanged but the new point may
*           qualify as xw.  xtry must be a new bound on the best point.

            if (xtry .le. zero) then
               a       = xtry
               nsamea  = 0
            else
               b       = xtry
               nsameb  = 0
               Brakted = .true.
            end if

*           If xw has not been set or ftry is better than fw, update the
*           points accordingly.

            if (wset) then
               setxw = ftry .lt. fw  .or.  .not. Extrap
            else
               setxw = .true.
            end if

            if (setxw) then
               xw     = xtry
               fw     = ftry
               gw     = gtry
               wset   = .true.
               Extrap = .false.
            end if
         end if

*        ---------------------------------------------------------------
*        Check the termination criteria.  wset will always be true.
*        ---------------------------------------------------------------
         tol    = tolAbs + tolRel*alfbst
         truea  = alfbst + a
         trueb  = alfbst + b

         Found  = abs(gbest) .le. targetg
         QuitF  = numf       .ge. maxf    .and.  Moved
         QuitI  = b - a      .le. tol + tol

         if (QuitI  .and. .not. Moved) then

*           The interval of uncertainty appears to be small enough,
*           but no better point has been found.  Check that changing
*           alfa by b-a changes f by less than epsaf.

            tol    = tol/ten
            tolAbs = tol
            QuitI  =      tol .le. tolTiny .or.
     &               (abs(fw) .le. epsaf   .and.  gw .le. epsaf)
         end if

         Done  = QuitF  .or.  QuitI  .or.  Found

         if (Debug) then
            write(nout, 1200)
     &           truea    , trueb , b - a , tol   ,
     &           nsamea   , nsameb, numf  ,
     &           Brakted  , Extrap, closef, Improved,
     &           Found    , QuitI ,
     &           alfbst   , fbest , gbest ,
     &           alfbst+xw, fw    , gw
         end if

*        ---------------------------------------------------------------
*        Proceed with the computation of a trial steplength.
*        The choices are...
*        1. Parabolic fit using derivatives only, if the f values are
*           close.
*        2. Cubic fit for a minimizer, using both f and f'.
*        3. Damped cubic or parabolic fit if the regular fit appears to
*           be consistently overestimating the distance to a minimizer.
*        4. Bisection, geometric bisection, or a step of  tol  if
*           choices 2 or 3 are unsatisfactory.
*        ---------------------------------------------------------------
         if (.not. Done) then
            xmidpt = half*(a + b)
            s      = zero
            q      = zero

            if (closef) then
*              ---------------------------------------------------------
*              Fit a parabola to the two best gradient values.
*              ---------------------------------------------------------
               s      = gbest
               q      = gbest - gw
               if (Debug) write(nout, 2200)
            else
*              ---------------------------------------------------------
*              Fit cubic through  fbest  and  fw.
*              ---------------------------------------------------------
               if (Debug) write(nout, 2100)
               Fitok  = .true.
               r      = three*(fbest - fw)/xw + gbest + gw
               absr   = abs( r )
               s      = sqrt( abs( gbest ) ) * sqrt( abs( gw ) )

*              Compute  q =  the square root of  r*r - gbest*gw.
*              The method avoids unnecessary underflow and overflow.

               if ((gw .lt. zero  .and.  gbest .gt. zero) .or.
     &             (gw .gt. zero  .and.  gbest .lt. zero)) then
                  scale  = absr + s
                  if (scale .eq. zero) then
                     q  = zero
                  else
                     q  = scale*sqrt( (absr/scale)**2 + (s/scale)**2 )
                  end if
               else if (absr .ge. s) then
                  q     = sqrt(absr + s)*sqrt(absr - s)
               else
                  Fitok = .false.
               end if

               if (Fitok) then

*                 Compute a minimizer of the fitted cubic.

                  if (xw .lt. zero) q = - q
                  s  = gbest -  r - q
                  q  = gbest - gw - q - q
               end if
            end if
*           ------------------------------------------------------------
*           Construct an artificial interval  (artifa, artifb)  in which
*           the new estimate of a minimizer must lie.  Set a default
*           value of xtry that will be used if the polynomial fit fails.
*           ------------------------------------------------------------
            artifa = a
            artifb = b
            if (.not. Brakted) then

*              A minimizer has not been bracketed.  Set an artificial
*              upper bound by expanding the interval  xw  by a suitable
*              factor.

               xtry   = - factor*xw
               artifb =   xtry
               if (alfbst + xtry .lt. alfmax) factor = five*factor

            else if (Extrap) then

*              The points are configured for an extrapolation.
*              Set a default value of  xtry  in the interval  (a, b)
*              that will be used if the polynomial fit is rejected.  In
*              the following,  dtry  and  daux  denote the lengths of
*              the intervals  (a, b)  and  (0, xw)  (or  (xw, 0),  if
*              appropriate).  The value of  xtry is the point at which
*              the exponents of  dtry  and  daux  are approximately
*              bisected.

               daux = abs( xw )
               dtry = b - a
               if (daux .ge. dtry) then
                  xtry = five*dtry*(point1 + dtry/daux)/eleven
               else
                  xtry = half * sqrt( daux ) * sqrt( dtry )
               end if
               if (xw .gt. zero)   xtry = - xtry
               if (Debug) write(nout, 2400) xtry, daux, dtry

*              Reset the artificial bounds.  If the point computed by
*              extrapolation is rejected,  xtry will remain at the
*              relevant artificial bound.

               if (xtry .le. zero) artifa = xtry
               if (xtry .gt. zero) artifb = xtry
            else

*              The points are configured for an interpolation.  The
*              default value xtry bisects the interval of uncertainty.
*              the artificial interval is just (a, b).

               xtry   = xmidpt
               if (Debug) write(nout, 2300) xtry
               if (nsamea .ge. 3  .or.  nsameb .ge. 3) then

*                 If the interpolation appears to be overestimating the
*                 distance to a minimizer,  damp the interpolation.

                  factor = factor / five
                  s      = factor * s
               else
                  factor = one
               end if
            end if
*           ------------------------------------------------------------
*           The polynomial fits give  (s/q)*xw  as the new step.
*           Reject this step if it lies outside  (artifa, artifb).
*           ------------------------------------------------------------
            if (q .ne. zero) then
               if (q .lt. zero) s = - s
               if (q .lt. zero) q = - q
               if (s*xw .ge. q*artifa  .and.  s*xw .le. q*artifb) then

*                 Accept the polynomial fit.

                  if (abs( s*xw ) .ge. q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = zero
                  end if
                  if (Debug) write(nout, 2500) xtry
               end if
            end if
         end if
      end if

*     ==================================================================

      if (.not. Done) then
         alfa  = alfbst + xtry
         if (Brakted  .or.  alfa .lt. alfmax - tolMax) then

*           The function must not be evaluated too close to a or b.
*           (It has already been evaluated at both those points.)

            if (xtry .le. a + tol  .or.  xtry .ge. b - tol) then
               if (half*(a + b) .le. zero) then
                  xtry = - tol
               else
                  xtry =   tol
               end if
               alfa = alfbst + xtry
            end if
         else

*           The step is close to, or larger than alfmax, replace it by
*           alfmax to force evaluation of  f  at the boundary.

            Brakted = .true.
            xtry    = alfmax - alfbst
            alfa    = alfmax
         end if
      end if

*     ------------------------------------------------------------------
*     Exit.
*     ------------------------------------------------------------------
      if (Done) then
         if      (BadFun) then
            iExit = 8           ! bad arguments
         else if (Found) then
            if (alfbst .lt. alfmax) then
               iExit = 1        ! Sufficient decrease
            else
               iExit = 2        ! Suff. Decrease on the boundary
            end if
         else if (Moved ) then
            iExit = 3           ! Decr at boundary or max funs
         else if (QuitF) then
            iExit = 7           ! No new point after max funs
         else if (Cramped) then
            iExit = 4           ! alfmax too mall
         else
            iExit = 6           ! [a,b] too small
         end if
      end if

      if (Debug) write(nout, 3000)
      return

 1000 format(/'     g0  tolAbs  alfmax        ', 1p, 2e22.14,   e16.8
     &       /' targetg tolRel   epsaf        ', 1p, 2e22.14,   e16.8
     &       /' Cramped                       ',  l3)
 1100 format(/' alfa    ftry    gtry          ', 1p, 2e22.14,   e16.8)
 1200 format(/' a       b       b - a   tol   ', 1p, 2e22.14,  2e16.8
     &       /' nsamea  nsameb  numf          ', 3i3
     &       /' Brakted Extrap  closef  Imprvd', 4l3
     &       /' Found   QuitI                 ', 2l3
     &       /' alfbst  fbest   gbest         ', 1p, 3e22.14
     &       /' alfaw   fw      gw            ', 1p, 3e22.14)
 2100 format( ' Cubic.')
 2200 format( ' Parabola.')
 2300 format( ' Bisection.              xmidpt', 1p,  e22.14)
 2400 format( ' Geo. bisection. xtry,daux,dtry', 1p, 3e22.14)
 2500 format( ' Polynomial fit accepted.  xtry', 1p,  e22.14)
 3000 format( ' ----------------------------------------------------'/)

      end ! subroutine srchc

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine srchq
     &   ( iExit  , First  , Debug  , Done  , Improved,
     &     maxf   , numf   , nout   ,
     &     alfmax , alfsml , epsaf  ,
     &     g0     , targetg, ftry   ,
     &     tolAbs , tolRel , tolTiny,
     &     alfa   , alfbst , fbest  ,
     &     Brakted, Cramped, Extrap , Moved , vset  , wset,
     &     nsamea , nsameb ,
     &     a      , b      , fa     , factor,
     &     xtry   , xw     , fw     , xv    , fv    , tolMax )

      implicit
     &     none
      logical
     &     First, Debug, Done, Improved, Brakted, Cramped, Extrap,
     &     Moved, vset , wset
      integer
     &     iExit, maxf, numf, nout, nsamea, nsameb
      double precision
     &     alfmax, alfsml, epsaf, g0, targetg, ftry, tolAbs, tolRel,
     &     tolTiny, alfa, alfbst, fbest, a, b, factor, xtry,
     &     xw, fw, gw, tolMax

*     ==================================================================
*     srchq  finds a sequence of improving estimates of a minimizer of
*     the univariate function f(alpha) in the interval (0,alfmax].
*     f(alpha) is a smooth function such that  f(0) = 0  and  f'(0) < 0.
*     srchq  requires  f(alpha) (but not f'(alpha)) to be evaluated
*     in the interval.  New estimates of a minimizer are computed using
*     safeguarded quadratic interpolation.
*
*     Reverse communication is used to allow the calling program to
*     evaluate f.  Some of the parameters must be set or tested by the
*     calling program.  The remainder would ordinarily be local
*     variables.
*
*     Input parameters (relevant to the calling program)
*     --------------------------------------------------
*
*     First         must be true on the first entry.  It is subsequently
*                   altered by srchq.
*
*     Debug         specifies whether detailed output is wanted.
*
*     maxf          is an upper limit on the number of times srchq is
*                   to be entered consecutively with Done = false
*                   (following an initial entry with First = true).
*
*     alfa          is the first estimate of a minimizer.  alfa is
*                   subsequently altered by srchq (see below).
*
*     alfmax        is the upper limit of the interval to be searched.
*
*     alfsml        is intended to prevent inefficiency when a minimizer
*                   is very small, for cases where the calling program
*                   would prefer to redefine f'(alfa).  alfsml is
*                   allowed to be zero.  Early termination will occur if
*                   srchq determines that a minimizer lies somewhere in
*                   the interval [0, alfsml) (but not if alfmax is
*                   smaller that alfsml).
*
*     epsaf         is an estimate of the absolute precision in the
*                   computed value of f(0).
*
*     ftry          the value of f at the new point
*                   alfa = alfbst + xtry.
*
*     g0            is the value of f'(0).  g0 must be negative.
*
*     tolAbs,tolRel define a function tol(alfa) = tolRel*alfa + tolAbs
*                   such that if f has already been evaluated at alfa,
*                   it will not be evaluated closer than tol(alfa).
*                   These values may be reduced by srchq.
*
*     targetg       is the target value of abs(f'(alfa)). The search
*                   is terminated when
*                    abs(f'(alfa)) le targetg and f(alfa) lt 0.
*
*     tolTiny       is the smallest value that tolAbs is allowed to be
*                   reduced to.
*
*     Output parameters (relevant to the calling program)
*     ---------------------------------------------------
*
*     Improved      is true if the previous alfa was the best point so
*                   far.  Any related quantities should be saved by the
*                   calling program (e.g., arrays) before paying
*                   attention to the variable Done.
*
*     Done = false  means the calling program should evaluate ftry
*                   for the new trial step alfa, and reenter srchq.
*
*     Done = true   means that no new alfa was calculated.  The value
*                   of iExit gives the result of the search as follows
*
*                   iExit = 1 means the search has terminated
*                             successfully with alfbst < alfmax.
*
*                   iExit = 2 means the search has terminated
*                             successfully with alfbst = alfmax.
*
*                   iExit = 3 means that the search failed to find a
*                             point of sufficient decrease in maxf
*                             functions, but a lower point was found.
*
*                   iExit = 4 means alfmax is so small that a search
*                             should not have been attempted.
*
*                   iExit = 5 means that the search was terminated
*                             because of alfsml (see above).
*
*                   iExit = 6 means the search has failed to find a
*                             useful step.  The interval of uncertainty
*                             is [0,b] with b < 2*tolAbs. A minimizer
*                             lies very close to alfa = 0, or f'(0) is
*                             not sufficiently accurate.
*
*                   iExit = 7 if no better point could be found after
*                             maxf  function calls.
*
*                   iExit = 8 means the input parameters were bad.
*                             alfmax le tolTiny  or  g0 ge zero.
*                             No function evaluations were made.
*
*     numf          counts the number of times srchq has been entered
*                   consecutively with Done = false (i.e., with a new
*                   function value ftry).
*
*     alfa          is the point at which the next function ftry must
*                   be computed.
*
*     alfbst        should be accepted by the calling program as the
*                   approximate minimizer, whenever srchq returns
*                   iExit = 1, 2 or 3.
*
*     fbest         will be the corresponding value of f.
*
*     The following parameters retain information between entries
*     -----------------------------------------------------------
*
*     Brakted       is false if f has not been evaluated at the far end
*                   of the interval of uncertainty.  In this case, the
*                   point b will be at alfmax + tol(alfmax).
*
*     Cramped       is true if alfmax is very small (le tolAbs).  If the
*                   search fails, this indicates that a zero step should
*                   be taken.
*
*     Extrap        is true if alfbst has moved at least once and xv
*                   lies outside the interval of uncertainty.  In this
*                   case, extra safeguards are applied to allow for
*                   instability in the polynomial fit.
*
*     Moved         is true if a better point has been found, i.e.,
*                   alfbst gt 0.
*
*     vset          records whether a third-best point has been defined.
*
*     wset          records whether a second-best point has been
*                   defined.  It will always be true by the time the
*                   convergence test is applied.
*
*     nsamea        is the number of consecutive times that the
*                   left-hand end point of the interval of uncertainty
*                   has remained the same.
*
*     nsameb        similarly for the right-hand end.
*
*     a, b, alfbst  define the current interval of uncertainty.
*                   A minimizer lies somewhere in the  interval
*                   [alfbst + a, alfbst + b].
*
*     alfbst        is the best point so far.  It lies strictly within
*                   [atrue,btrue]  (except when alfbst has not been
*                   moved, in which case it lies at the left-hand end
*                   point).  Hence we have a .le. 0 and b .gt. 0.
*
*     fbest         is the value of f at the point alfbst.
*
*     fa            is the value of f at the point alfbst + a.
*
*     factor        controls the rate at which extrapolated estimates of
*                   alfa  may expand into the interval of uncertainty.
*                   Factor is not used if a minimizer has been bracketed
*                   (i.e., when the variable Brakted is true).
*
*     fv, fw        are the values of f at the points alfbst + xv  and
*                   alfbst + xw.  They are not defined until  vset  or
*                   wset  are true.
*
*     xtry          is the trial point within the shifted interval
*                   (a, b).  The new trial function value must be
*                   computed at the point alfa = alfbst + xtry.
*
*     xv            is such that alfbst + xv is the third-best point.
*                   It is not defined until vset is true.
*
*     xw            is such that alfbst + xw is the second-best point.
*                   It is not defined until wset is true.  In some
*                   cases,  xw will replace a previous xw that has a
*                   lower function but has just been excluded from
*                   (a,b).
*
*     Systems Optimization Laboratory, Stanford University, California.
*     Original version February 1982.  Rev. May 1983.
*     Original F77 version 22-August-1985.
*     17 Jul 1997: Removed saved variables for thread-safe version.
*     31 Jul 1999: Current version.
*     ==================================================================
      logical
     &     BadFun, CloseF, Found, QuitF, QuitFZ, QuitI , QuitS,
     &     setxv , xinxw
      double precision
     &     artifa, artifb, daux, dtry, endpnt, fa, fv, gv,
     &     q, s, tol, truea, trueb, xmidpt, xv
*     ------------------------------------------------------------------
      double precision   zero,          point1,          half
      parameter         (zero = 0.0d+0, point1 = 0.1d+0, half = 0.5d+0)
      double precision   one,           two,             five
      parameter         (one  = 1.0d+0, two    = 2.0d+0, five = 5.0d+0)
      double precision   ten,           eleven
      parameter         (ten  = 1.0d+1, eleven = 1.1d+1               )
*     ------------------------------------------------------------------
*     Local variables
*     ===============
*
*     CloseF     is true if the worst function fv is within epsaf of
*                fbest (up or down).
*
*     Found      is true if the sufficient decrease conditions holds at
*                alfbst.
*
*     QuitF      is true when  maxf  function calls have been made.
*
*     QuitFZ     is true when the three best function values are within
*                epsaf of each other, and the new point satisfies
*                fbest le ftry le fbest+epsaf.
*
*     QuitI      is true when the interval of uncertainty is less than
*                2*tol.
*
*     QuitS      is true as soon as alfa is too small to be useful;
*                i.e., btrue le alfsml.
*
*     xinxw      is true if xtry is in (xw,0) or (0,xw).
*     ------------------------------------------------------------------

      Improved = .false.
      BadFun   = .false.
      QuitF    = .false.
      QuitFZ   = .false.
      QuitS    = .false.
      QuitI    = .false.

      if (First) then
*        ---------------------------------------------------------------
*        First entry.  Initialize various quantities, check input data
*        and prepare to evaluate the function at the initial step alfa.
*        ---------------------------------------------------------------
         First  = .false.
         numf   = 0
         alfbst = zero
         BadFun = alfmax .le. tolTiny  .or.  g0 .ge. zero
         Done   = BadFun
         Moved  = .false.

         if (.not. Done) then
            Brakted = .false.
            Cramped = alfmax .le. tolAbs
            Extrap  = .false.
            vset    = .false.
            wset    = .false.
            nsamea  = 0
            nsameb  = 0

            tolMax  = tolRel*alfmax + tolAbs
            a       = zero
            b       = alfmax + tolMax
            fa      = zero
            factor  = five
            tol     = tolAbs
            xtry    = alfa
            if (Debug) then
               write(nout, 1000)
     &              g0     , tolAbs, alfmax,
     &              targetg, tolRel, epsaf , Cramped
            end if
         end if
      else
*        ---------------------------------------------------------------
*        Subsequent entries.  The function has just been evaluated at
*        alfa = alfbst + xtry,  giving ftry.
*        ---------------------------------------------------------------
         if (Debug)
     &   write(nout, 1100)
     &        alfa, ftry

         numf   = numf   + 1
         nsamea = nsamea + 1
         nsameb = nsameb + 1

         if (.not. Brakted) then
            tolMax = tolAbs + tolRel*alfmax
            b      = alfmax - alfbst + tolMax
         end if

*        Check if xtry is in the interval (xw,0) or (0,xw).

         if (wset) then
            xinxw =       (zero .lt. xtry  .and.  xtry .le. xw  )
     &               .or. (  xw .le. xtry  .and.  xtry .lt. zero)
         else
            xinxw = .false.
         end if

         Improved = ftry .lt. fbest
         if (vset) then
            CloseF = abs( fbest - fv ) .le. epsaf
         else
            CloseF = .false.
         end if

         if (Improved) then

*           We seem to have an improvement.  The new point becomes the
*           origin and other points are shifted accordingly.

            if (wset) then
               xv     = xw - xtry
               fv     = fw
               vset   = .true.
            end if

            xw     = zero - xtry
            fw     = fbest
            wset   = .true.
            fbest  = ftry
            alfbst = alfa
            Moved  = .true.

            a      = a    - xtry
            b      = b    - xtry
            Extrap = .not. xinxw

*           Decrease the length of (a,b).

            if (xtry .ge. zero) then
               a       = xw
               fa      = fw
               nsamea  = 0
            else
               b       = xw
               nsameb  = 0
               Brakted = .true.
            end if
         else if (CloseF  .and.  ftry - fbest .lt. epsaf) then

*           Quit if there has been no progress and ftry, fbest, fw
*           and fv are all within epsaf of each other.

            QuitFZ = .true.
         else

*           The new function value is no better than the current best
*           point.  xtry must an end point of the new (a,b).

            if (xtry .lt. zero) then
               a       = xtry
               fa      = ftry
               nsamea  = 0
            else
               b       = xtry
               nsameb  = 0
               Brakted = .true.
            end if

*           The origin remains unchanged but xtry may qualify as xw.

            if (wset) then
               if (ftry .lt. fw) then
                  xv     = xw
                  fv     = fw
                  vset   = .true.

                  xw     = xtry
                  fw     = ftry
                  if (Moved) Extrap = xinxw
               else if (Moved) then
                  if (vset) then
                     setxv = ftry .lt. fv  .or.  .not. Extrap
                  else
                     setxv = .true.
                  end if

                  if (setxv) then
                     if (vset  .and.  xinxw) then
                        xw = xv
                        fw = fv
                     end if
                     xv   = xtry
                     fv   = ftry
                     vset = .true.
                  end if
               else
                  xw  = xtry
                  fw  = ftry
               end if
            else
               xw     = xtry
               fw     = ftry
               wset   = .true.
            end if
         end if

*        ---------------------------------------------------------------
*        Check the termination criteria.
*        ---------------------------------------------------------------
         tol    = tolAbs + tolRel*alfbst
         truea  = alfbst + a
         trueb  = alfbst + b

         Found  = Moved  .and.  abs(fa - fbest) .le. -a*targetg
         QuitF  = numf  .ge. maxf
         QuitI  = b - a .le. tol + tol
         QuitS  = trueb .le. alfsml

         if (QuitI  .and.  .not. Moved) then

*           The interval of uncertainty appears to be small enough,
*           but no better point has been found.  Check that changing
*           alfa by b-a changes f by less than epsaf.

            tol    = tol/ten
            tolAbs = tol
            QuitI  = abs(fw) .le. epsaf  .or.  tol .le. tolTiny
         end if

         Done  = QuitF  .or.  QuitFZ  .or.  QuitS  .or.  QuitI
     &                  .or.  Found

         if (Debug) then
            write(nout, 1200)
     &           truea    , trueb , b-a   , tol     ,
     &           nsamea   , nsameb, numf  ,
     &           Brakted  , Extrap, CloseF, Improved,
     &           Found    , QuitI , QuitFZ, QuitS   ,
     &           alfbst   , fbest ,
     &           alfbst+xw, fw
            if (vset) then
               write(nout, 1300) alfbst + xv, fv
            end if
         end if

*        ---------------------------------------------------------------
*        Proceed with the computation of an estimate of a minimizer.
*        The choices are...
*        1. Parabolic fit using function values only.
*        2. Damped parabolic fit if the regular fit appears to be
*           consistently overestimating the distance to a minimizer.
*        3. Bisection, geometric bisection, or a step of tol if the
*           parabolic fit is unsatisfactory.
*        ---------------------------------------------------------------
         if (.not. Done) then
            xmidpt = half*(a + b)
            s      = zero
            q      = zero

*           ============================================================
*           Fit a parabola.
*           ============================================================
*           See if there are two or three points for the parabolic fit.

            gw = (fw - fbest)/xw
            if (vset  .and.  Moved) then

*              Three points available.  Use fbest, fw and fv.

               gv = (fv - fbest)/xv
               s  = gv - (xv/xw)*gw
               q  = two*(gv - gw)
               if (Debug) write(nout, 2200)
            else

*              Only two points available.  Use fbest, fw and g0.

               if (Moved) then
                  s  = g0 - two*gw
               else
                  s  = g0
               end if
               q = two*(g0 - gw)
               if (Debug) write(nout, 2100)
            end if

*           ------------------------------------------------------------
*           Construct an artificial interval (artifa, artifb) in which
*           the new estimate of the steplength must lie.  Set a default
*           value of  xtry  that will be used if the polynomial fit is
*           rejected. In the following, the interval (a,b) is considered
*           the sum of two intervals of lengths  dtry  and  daux, with
*           common end point the best point (zero).  dtry is the length
*           of the interval into which the default xtry will be placed
*           and endpnt denotes its non-zero end point.  The magnitude of
*           xtry is computed so that the exponents of dtry and daux are
*           approximately bisected.
*           ------------------------------------------------------------
            artifa = a
            artifb = b
            if (.not. Brakted) then

*              A minimizer has not yet been bracketed.
*              Set an artificial upper bound by expanding the interval
*              xw  by a suitable factor.

               xtry   = - factor*xw
               artifb =   xtry
               if (alfbst + xtry .lt. alfmax) factor = five*factor
            else if (vset .and. Moved) then

*              Three points exist in the interval of uncertainty.
*              Check if the points are configured for an extrapolation
*              or an interpolation.

               if (Extrap) then

*                 The points are configured for an extrapolation.

                  if (xw .lt. zero) endpnt = b
                  if (xw .gt. zero) endpnt = a
               else

*                 If the interpolation appears to be overestimating the
*                 distance to a minimizer,  damp the interpolation step.

                  if (nsamea .ge. 3  .or.   nsameb .ge. 3) then
                     factor = factor / five
                     s      = factor * s
                  else
                     factor = one
                  end if

*                 The points are configured for an interpolation.  The
*                 artificial interval will be just (a,b).  Set endpnt so
*                 that xtry lies in the larger of the intervals (a,b)
*                 and  (0,b).

                  if (xmidpt .gt. zero) then
                     endpnt = b
                  else
                     endpnt = a
                  end if

*                 If a bound has remained the same for three iterations,
*                 set endpnt so that  xtry  is likely to replace the
*                 offending bound.

                  if (nsamea .ge. 3) endpnt = a
                  if (nsameb .ge. 3) endpnt = b
               end if

*              Compute the default value of  xtry.

               dtry = abs( endpnt )
               daux = b - a - dtry
               if (daux .ge. dtry) then
                  xtry = five*dtry*(point1 + dtry/daux)/eleven
               else
                  xtry = half*sqrt( daux )*sqrt( dtry )
               end if
               if (endpnt .lt. zero) xtry = - xtry
               if (Debug) write(nout, 2500) xtry, daux, dtry

*              If the points are configured for an extrapolation set the
*              artificial bounds so that the artificial interval lies
*              within (a,b).  If the polynomial fit is rejected,  xtry
*              will remain at the relevant artificial bound.

               if (Extrap) then
                  if (xtry .le. zero) then
                     artifa = xtry
                  else
                     artifb = xtry
                  end if
               end if
            else

*              The gradient at the origin is being used for the
*              polynomial fit.  Set the default xtry to one tenth xw.

               if (Extrap) then
                  xtry = - xw
               else
                  xtry   = xw/ten
               end if
               if (Debug) write(nout, 2400) xtry
            end if

*           ------------------------------------------------------------
*           The polynomial fits give (s/q)*xw as the new step.  Reject
*           this step if it lies outside (artifa, artifb).
*           ------------------------------------------------------------
            if (q .ne. zero) then
               if (q .lt. zero) s = - s
               if (q .lt. zero) q = - q
               if (s*xw .ge. q*artifa   .and.   s*xw .le. q*artifb) then

*                 Accept the polynomial fit.

                  if (abs( s*xw ) .ge. q*tol) then
                     xtry = (s/q)*xw
                  else
                     xtry = zero
                  end if
                  if (Debug) write(nout, 2600) xtry
               end if
            end if
         end if
      end if
*     ==================================================================

      if (.not. Done) then
         alfa  = alfbst + xtry
         if (Brakted  .or.  alfa .lt. alfmax - tolMax) then

*           The function must not be evaluated too close to a or b.
*           (It has already been evaluated at both those points.)

            xmidpt = half*(a + b)
            if (xtry .le. a + tol  .or.  xtry .ge. b - tol) then
               if (xmidpt .le. zero) then
                  xtry = - tol
               else
                  xtry =   tol
               end if
            end if

            if (abs( xtry ) .lt. tol) then
               if (xmidpt .le. zero) then
                  xtry = - tol
               else
                  xtry =   tol
               end if
            end if
            alfa  = alfbst + xtry
         else

*           The step is close to or larger than alfmax, replace it by
*           alfmax to force evaluation of the function at the boundary.

            Brakted = .true.
            xtry    = alfmax - alfbst
            alfa    = alfmax
         end if
      end if
*     ------------------------------------------------------------------
*     Exit.
*     ------------------------------------------------------------------
      if (Done) then
         if      (BadFun) then
            iExit = 8
         else if (QuitS ) then
            iExit = 5
         else if (Found) then
            if (alfbst .lt. alfmax) then
               iExit = 1
            else
               iExit = 2
            end if
         else if (Moved ) then
            iExit = 3
         else if (QuitF) then
            iExit = 7
         else if (Cramped) then
            iExit = 4
         else
            iExit = 6
         end if
      end if

      if (Debug) write(nout, 3000)
      return

 1000 format(/'     g0  tolAbs  alfmax        ', 1p, 2e22.14,   e16.8
     &       /' targetg tolRel   epsaf        ', 1p, 2e22.14,   e16.8
     &       /' Cramped                       ',  l3)
 1100 format(/' alfa    ftry                  ', 1p,2e22.14          )
 1200 format(/' a       b       b - a   tol   ', 1p,2e22.14,   2e16.8
     &       /' nsamea  nsameb  numf          ', 3i3
     &       /' Brakted Extrap  CloseF  Imprvd', 4l3
     &       /' Found   QuitI   QuitFZ  QuitS ', 4l3
     &       /' alfbst  fbest                 ', 1p,2e22.14
     &       /' alfaw   fw                    ', 1p,2e22.14)
 1300 format( ' alfav   fv                    ', 1p,2e22.14 /)
 2100 format( ' Parabolic fit,    two points. ')
 2200 format( ' Parabolic fit,  three points. ')
 2400 format( ' Exponent reduced.  Trial point', 1p,  e22.14)
 2500 format( ' Geo. bisection. xtry,daux,dtry', 1p, 3e22.14)
 2600 format( ' Polynomial fit accepted.  xtry', 1p,  e22.14)
 3000 format( ' ----------------------------------------------------'/)

      end ! subroutine srchq
