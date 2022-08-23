!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn05wrpa.f  --- user-function wrapper for  snOptA.
!
!     s0fgA    s0fgA1
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s0fgA
     &   ( iExit,
     &     modefg, NonlinearCon, NonlinearObj,
     &     n, negCon, nnCon0, nnCon,
     &     nnJac, nnH, nnObj0, nnObj,
     &     userfg, dummyf, dummyH,
     &     x, yCon,
     &     neJ, nlocJ, locJ, indJ,
     &     neH, nlocH, locH, indH, Hcol,
     &     fCon, fObj, gCon, gObj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     userfg, dummyf, dummyH
      logical
     &     NonlinearCon, NonlinearObj
      integer
     &     iExit, lencu, leniu, lenru, lencw, leniw, lenrw, modefg,
     &     n, negCon, neH, neJ, nlocH, nlocJ, nnCon0, nnCon, nnJac,
     &     nnH, nnObj0, nnObj, indH(neH), indJ(neJ), locH(nlocH),
     &     locJ(nlocJ), iu(leniu), iw(leniw)
      double precision
     &     fObj, fCon(nnCon0), gObj(nnObj0), gCon(negCon),
     &     Hcol(neH), x(n), yCon(nnCon0), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s0fgA  is an instance of funwrapper that calls the user-written
!     routine  userfg  to evaluate the problem functions and possibly
!     their gradients.
!
!     Subroutine  userfg  is called using modefg to control
!     the gradients as follows:
!
!     modefg        Task
!     ------        ----
!       2     Assign fCon, fObj and all known elements of gCon and gObj.
!       1     Assign all known elements of gCon and gObj.
!             (fObj and fCon are ignored).
!       0     Assign fObj, fCon.  (gCon and gObj are ignored).
!
!     Since objective and constraints are computed simultaneously,
!     the input variables  NonlinearCon  and  NonlinearObj are ignored.
!
!     31 Oct 1998: First version based on s0fg in SNOPT 5.3-4.
!     03 Aug 2003: snEXIT and snPRNT adopted.
!     16 Jun 2008: Call-status implemented correctly.
!     15 Nov 2010: Call-status removed from the argument list.
!     07 Nov 2014: yCon added as argument.
!     13 Nov 2014: Scales applied for lvlScale > 0.
!     14 Feb 2015: proxWeight added.
!     ==================================================================
      character
     &     str*80
      logical
     &     Scaled
      integer
     &     gotFD, gotG, k, kG, knownG, lscales,
     &     lenG, lgConU, lgObjU, lindG, llocG, liGfun, ljGvar, lkxN,
     &     lyCon, lF, lFmul, lG, lnGlin, lvlDer, lvlScale, lxN,
     &     lxScaled, mode, nF, neG, nkx, nfCon1, nfCon2, nfObj1,
     &     nfObj2, nlocG, Status
      double precision
     &     Gdummy
!     ------------------------------------------------------------------
      parameter         (lvlDer =  70) ! = 0,1,2 or 3, deriv. level
      parameter         (gotFD  = 183) ! > 0 => some differences needed
      parameter         (gotG   = 184) ! > 0 => some exact derivs
      parameter         (nfCon1 = 189) ! calls to fCon: mode = 0
      parameter         (nfCon2 = 190) ! calls to fCon  mode > 0
      parameter         (nfObj1 = 194) ! calls to fObj: mode = 0
      parameter         (nfObj2 = 195) ! calls to fObj: mode > 0
!     ------------------------------------------------------------------
      lvlScale  = iw( 75) ! scale option
      nF        = iw(248) ! # of components of user-defined F
      neG       = iw(249) ! # of components of user-defined G
      lkxN      = iw(252) ! jN = kxN(j ) => col j of Jcol is variable jN
      llocG     = iw(260) ! locG(nlocG) = column pointers for indG
      lindG     = iw(261) ! indG(neG) holds the row indices for gij
      lnGlin    = iw(262) ! nGlin(j) = # linear elems in col j of gCon
      liGfun    = iw(266) ! iGfun(neG) row list of reordered G nonzeros
      ljGvar    = iw(267) ! iGvar(neG) col list of reordered G nonzeros
      lscales   = iw(296) ! scales(nb)  = row and column scales
      lxScaled  = iw(302) ! xScaled(n)  = copy of scaled  x
      lgConU    = iw(319) ! record of unknown derivatives and constants
      lgObjU    = iw(323) ! record of unknown derivatives
      lxN       = iw(327) ! xN(n)       = variables in natural order
      lF        = iw(328) ! F(nF)       = user-defined F
      lFmul     = iw(329) ! Fmul(nF)    = user-defined multipliers
      lG        = iw(330) ! G (lenG)    = problem derivatives
      lyCon     = iw(348) ! yCon(nnCon) = multipliers for F

      Gdummy    = rw( 69) ! definition of an 'unset' value

      iExit     = 0
      Scaled    = lvlScale .gt. 0

      nlocG     = nnJac + 1
      lenG      = max( neG, 1 )
      nkx       = n     + nF

      ! Determine the user-function call-status.

      call s8callStatus( Status, iw, leniw )

      if (Status .eq. 1) then
         !--------------------------------------------------------------
         ! First evaluation of the problem functions in snOptA
         ! On entry, lvlScale = 0.
         !--------------------------------------------------------------
         iw(gotFD) =  0 ! > 0 => some differences needed
         iw(gotG)  =  0 ! > 0 => some exact derivatives provided
         call snPRNT( 13, ' ', iw, leniw )
         call dload ( neG, Gdummy, rw(lG), 1 )
      end if

      !-----------------------------------------------------------------
      ! Unscale x.
      !-----------------------------------------------------------------
      if (Scaled) then
         call dcopy ( n,           x, 1, rw(lxScaled), 1 )
         call ddscl ( n, rw(lscales), 1,            x, 1 )
      end if

      !-----------------------------------------------------------------
      ! Compute the user-defined functions and derivatives.
      !-----------------------------------------------------------------
      mode = Status
      call s0fgA1
     &   ( modefg, mode,
     &     nF, n, nkx, nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &     iw(lkxN), userfg, rw(lxN), x,
     &     rw(lF), rw(lFmul), fCon, rw(lyCon), fObj, gObj,
     &     negCon, nlocG, iw(llocG), iw(lindG), gCon, rw(lgConU),
     &     iw(lnGlin), iw(liGfun), iw(ljGvar), lenG, neG, rw(lG),
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      iw(nfCon1) = iw(nfCon1) + 1
      iw(nfObj1) = iw(nfObj1) + 1
      if (modefg .gt. 0) then
         iw(nfCon2) = iw(nfCon2) + 1
         iw(nfObj2) = iw(nfObj2) + 1
      end if

      !-----------------------------------------------------------------
      ! Scale  x and the derivatives.
      !-----------------------------------------------------------------
      if (Scaled) then
         call dcopy
     &      ( n    , rw(lxScaled) , 1, x   , 1 )
         call dddiv
     &      ( nnCon, rw(lscales+n), 1, fCon, 1 )
         if (modefg .gt. 0  .and.  iw(gotG) .gt. 0) then
            call s8scaleg
     &         ( nnObj, rw(lscales), gObj, rw, lenrw )
            call s8scaleJ
     &         ( nnCon, nnJac, negCon, n, rw(lscales),
     &           neJ, nlocJ, locJ, indJ, gCon, rw, lenrw )
         end if
      end if

      if (mode .lt. 0) then
         !--------------------------------------------------------------
         ! The user may be saying the function is undefined (mode = -1)
         ! or may just want to stop                         (mode < -1).
         !--------------------------------------------------------------
         if (mode .eq. -1) then
            iExit = -1
         else
            iExit = 71
         end if
      end if

      !=================================================================
      ! Do some housekeeping after the first call.
      !=================================================================
      if (Status .eq. 1  .and.  iExit .eq. 0) then
         !--------------------------------------------------------------
         ! Count how many Jacobian elements are provided.
         !--------------------------------------------------------------
         knownG = 0
         kG     = lG
         do k = 1, neG
            if (rw(kG) .ne. Gdummy) knownG  = knownG + 1
            kG = kG + 1
         end do

         write(str, 1100) knownG, neG
         call snPRNT( 3, str, iw, leniw )

         if (knownG .lt. neG) iw(gotFD) = 1
         if (knownG .gt.   0) iw(gotG ) = 1

         if (knownG .lt. neG) then
            !-----------------------------------------------------------
            ! Missing derivatives.
            ! Keep a record of them in gObjU and gConU.
            ! Reduce the derivative level if necessary.
            !-----------------------------------------------------------
            call dcopy ( nnObj , gObj, 1, rw(lgObjU), 1 )
            call dcopy ( negCon, gCon, 1, rw(lgConU), 1 )

            if (iw(lvlDer) .eq. 3) then
               iw(lvlDer) = 0
               write(str, 2100)
               call snPRNT( 3, str, iw, leniw )
            end if
         end if
      end if

      return

 1100 format(' The user has defined', i8, '   out of', i8,
     &       '   first  derivatives')
 2100 format(' XXX  Some first  derivatives are missing ---',
     &       ' derivative level reduced to 0')

      end ! subroutine s0fgA

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s0fgA1
     &   ( modefg, Status,
     &     nF, n, nkx, nnCon0, nnCon, nnJac, nnH, nnObj0, nnObj,
     &     kxN, userfg, xN, x,
     &     F, Fmul, fCon, yCon, fObj, gObj,
     &     negCon, nlocG, locG, indG, gCon, gConU, nGlin,
     &     iGfun, jGvar, lenG, neG, G,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     userfg
      integer
     &     lenG, lencu, leniu, lenru, lencw, leniw, lenrw, modefg,
     &     n, nF, neG, negCon, nkx, nlocG, nnCon, nnCon0, nnJac, nnH,
     &     nnObj0, nnObj, Status, kxN(nkx), locG(nlocG), indG(negCon),
     &     nGlin(nLocG), iGfun(lenG), jGvar(lenG), iu(leniu), iw(leniw)
      double precision
     &     fObj, F(nF), Fmul(nF), fCon(nnCon0), gCon(negCon),
     &     gConU(negCon), gObj(nnObj0),
     &     G(lenG), yCon(nnCon0), x(n), xN(n), ru(lenru),
     &     rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s0fgA1   does the work for s0fgA.
!
!     08 Nov 1998: First version of s0fgA1.
!     09 Apr 1999: Updated for SnoptA.
!     01 Apr 2005: Cosmetic changes.
!     14 Feb 2015: proxWeight added.
!     ==================================================================
      external
     &     ddot
      integer
     &     i, iN, iObj, j, jN, k, lvlTim, lx0,
     &     minimize, needF, needG, needH, nextG, ObjRow
      double precision
     &     ddot, signObj, proxWeight
!     ------------------------------------------------------------------
      double precision   half
      parameter         (half = 0.5d+0)
      double precision   zero,          one
      parameter         (zero = 0.0d+0, one = 1.0d+0)
!     ------------------------------------------------------------------
      proxWeight = rw( 91) ! Proximal-point weight

      minimize   = iw( 87)      ! 1, -1  => MIN, MAX
      lvlTim     = iw(182) ! Timing level
      iObj       = iw(204) ! position of the objective row in J
      lx0        = iw(298) ! x0(nnL) = FP base point

      ObjRow     = iw(103) ! Objective row of user-defined F

      signObj    = minimize

      needF      = 0
      needG      = 0
      needH      = 0
      if (modefg .eq. 0  .or.  modefg .eq. 2) needF = 1
      if (modefg .eq. 1  .or.  modefg .eq. 2) needG = 1
      if (modefg .eq. 4                     ) needH = 1

      !-----------------------------------------------------------------
      ! Save x in natural order.
      ! For safety, zero out the linear components of x, just in case
      ! they are used by  userfg.
      !-----------------------------------------------------------------
      do j  = 1, n
         jN = kxN(j)
         if (j .le. nnH) then
            xN(jN) = x(j)
         else
            xN(jN) = zero
         end if
      end do

      if (needH .gt. 0) then
         !--------------------------------------------------------------
         ! The Hessian of the Lagrangian is requested.
         ! Expand the multipliers into Fmul.
         !--------------------------------------------------------------
         call dload ( nF, zero, Fmul, 1 )

         if (ObjRow .gt. 0) then
            Fmul(ObjRow) = one
         end if

         do i  = 1, nnCon
            j  = n + i
            iN = kxN(j)
            Fmul(iN) = -signObj*yCon(i)
         end do
      end if

      !-----------------------------------------------------------------
      ! Compute the user-defined functions.
      !-----------------------------------------------------------------
      if (lvlTim .ge. 2) call s1time( 4, 0, iw, leniw, rw, lenrw )
      call userfg
     &   ( Status, n, xN,
     &     needF, nF, F,
     &     needG, lenG, G,
     &     cu, lencu, iu, leniu, ru, lenru )
      if (lvlTim .ge. 2) call s1time(-4, 0, iw, leniw, rw, lenrw )

      !-----------------------------------------------------------------
      ! Gather the nonlinear elements of F.
      !-----------------------------------------------------------------
      if (needF .gt. 0) then
         do i  = 1, nnCon
            j  = n + i
            iN = kxN(j)
            fCon(i) = F(iN)
         end do

         if (nnObj .eq. 0) then
            fObj = zero
         else if (ObjRow .eq. 0) then
            call dcopy ( nnObj, x, 1, gObj, 1 )
            call daxpy ( nnObj, (-one), rw(lx0), 1, gObj, 1 )
            fObj = half*proxWeight*ddot ( nnObj, gObj, 1, gObj, 1 )
            call dscal ( nnObj, proxWeight, gObj, 1 )
         else
            fObj = F(ObjRow)
         end if
      end if

      if (needG .gt. 0) then

         ! Extract the derivatives from G.

         do k = 1, neG
            i = iGfun(k)
            j = jGvar(k)

            if (i .eq. iObj) then
               gObj(j) = G(k)
            else
               nextG       = locG(j)
               indG(nextG) = i
               gCon(nextG) = G(k)
               locG(j)     = nextG + 1
            end if
         end do
      else

         ! Keep locG in synch ready for any constant elements.

         do k = 1, neG
            i = iGfun(k)
            j = jGvar(k)
            if (i .ne. iObj) locG(j) = locG(j) + 1
         end do
      end if

      if (needF .gt. 0) then

         ! Add to fCon the linear term associated with every constant
         ! element in the nonlinear part of Jcol.

         do j  = 1, nnJac
            nextG = locG(j)
            do k  = 1, nGlin(j)
               i  = indG(nextG)
               fCon(i) = fCon(i) + gConU(nextG)*x(j)
               nextG   = nextG + 1
            end do
         end do
      end if

      ! locG(j) points to the first linear element in column j.
      ! Update it so that it points to the start of column j.

      do j = nnJac, 2, -1
         locG(j-1) = locG(j-1) + nGlin(j-1)
         locG(j)   = locG(j-1)
      end do
      locG(1) = 1

      end ! subroutine s0fgA1

