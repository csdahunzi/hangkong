!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn05wrpn.f  --- user-function interfaces for npOpt.
!
!     s0fgN
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s0fgN
     &   ( iExit,
     &     modefg, NonlinearCon, NonlinearObj,
     &     n, negCon, nnCon0, nnCon,
     &     nnJac, nnH, nnObj0, nnObj,
     &     fgcon, fgobj, dummyH,
     &     x, yCon,
     &     neJ, nlocJ, locJ, indJ,
     &     neH, nlocH, locH, indH, Hcol,
     &     fCon, fObj, gCon, gObj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     fgcon, fgobj, dummyH
      logical
     &     NonlinearCon, NonlinearObj
      integer
     &     iExit, lencu, leniu, lenru, lencw, leniw, lenrw, modefg,
     &     n, negCon, neH, neJ, nlocH, nlocJ, nnCon0, nnCon, nnJac,
     &     nnH, nnObj0, nnObj, indH(neH), indJ(neJ), locH(nlocH),
     &     locJ(nlocJ), iu(leniu), iw(leniw)
      double precision
     &     fObj, fCon(nnCon0),  gObj(nnObj0), gCon(negCon),
     &     Hcol(neH), x(n), yCon(nnCon0), ru(lenru), rw(lenrw)
      character
     &     cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     s0fgN   is a version of funwrapper that calls the user-written
!     routines  fgcon  and  fgobj  to evaluate the problem functions
!     and possibly their gradients.
!
!     Arguments  fgcon  and  fgobj  are called using modefg to control
!     the gradients as follows:
!
!     Version for the NPSOL interface snOptN.
!
!     modefg        Task
!     ------        ----
!       2     Assign fCon, fObj and all known elements of gCon and gObj.
!       1     Assign all known elements of gCon and gObj.
!             (fObj and fCon are ignored).
!       0     Assign fObj, fCon.  (gCon and gObj are ignored).
!
!     If s0fgN is called with minmax = 0 (feasible point only) then
!     nnObj = nnJac and the user objective is not used.
!
!     09-Jan 1992: First version of s0fgN  based on snwrap.
!     03 Aug 2003: snPRNT adopted.
!     16 Jun 2008: Call-status implemented correctly.
!     15 Nov 2010: Call-status removed from the argument list.
!     04 Nov 2014: neH, indH, locH, Hcol added as arguments.
!     07 Nov 2014: yCon added as argument.
!     14 Feb 2015: proxWeight added.
!     ==================================================================
      character
     &     str*80
      external
     &     ddot
      logical
     &     FPonly, Scaled
      integer
     &     gotFD, gotG, conJij, l, lG, lgConU, lgObjU, lgSave,
     &     liy1, lvlScale, lvlTim, lvlDer, lscales, lx0, lxScaled,
     &     minmax, modeC, modeF, nfCon1, nfCon2, nfObj1, nfObj2,
     &     supplied, nnGlin, nnObjU, Status, StatusUser
      double precision
     &     ddot, Gdummy, proxWeight
!     ------------------------------------------------------------------
      parameter         (lvlDer =  70) ! = 0,1,2 or 3, deriv. level
      parameter         (gotFD  = 183) ! > 0 => some differences needed
      parameter         (gotG   = 184) ! > 0 => some exact derivs
      parameter         (conJij = 185) ! > 0 => constant Jacob elements
      parameter         (nfCon1 = 189) ! calls to fCon: mode = 0
      parameter         (nfCon2 = 190) ! calls to fCon  mode > 0
      parameter         (nfObj1 = 194) ! calls to fObj: mode = 0
      parameter         (nfObj2 = 195) ! calls to fObj: mode > 0

      double precision   half,            one
      parameter         (half   = 0.5d+0, one   = 1.0d+0)
!     ------------------------------------------------------------------
      nnObjU     = iw( 22) ! # of objective variables
      lvlScale   = iw( 75) ! scale option
      minmax     = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
      lvlTim     = iw(182) ! Timing level

      lscales    = iw(296) ! scales(nb) = row and column scales
      lx0        = iw(298) ! x0(nnH)    = Feasible starting point
      lxScaled   = iw(302) ! xScaled(n) = copy of scaled x
      liy1       = iw(309) ! iy1(nb)    =  integer work vector
      lgConU     = iw(319) ! record of unknown derivatives and constants
      lgObjU     = iw(323) ! record of unknown derivatives
      lgSave     = iw(339) ! gSav(nnObjU) holds user-defined gObj

      Gdummy     = rw( 69) ! definition of an 'unset' value
      proxWeight = rw( 91) ! Proximal-point weight

      iExit      = 0

      FPonly     = minmax   .eq. 0
      Scaled     = lvlScale .eq. 2

      modeC      = modefg
      modeF      = modefg

      ! Determine the status of this call.

      call s8callStatus( Status, iw, leniw )

      if (Status .eq. 1) then
         !--------------------------------------------------------------
         ! First evaluation of the problem functions in npOpt
         ! On entry, lvlScale = 0.
         !--------------------------------------------------------------
         iw(gotFD)  =  0
         iw(gotG)   =  0
         iw(conJij) =  0
         call snPRNT( 13, ' ', iw, leniw )
         call dload ( negCon, Gdummy, gCon, 1 )
         call dload ( nnObj , Gdummy, gObj, 1 )
      end if

      !-----------------------------------------------------------------
      ! Unscale x.
      !-----------------------------------------------------------------
      if (Scaled) then
         call dcopy ( nnH,           x, 1, rw(lxScaled), 1 )
         call ddscl ( nnH, rw(lscales), 1,            x, 1 )

         ! If the Jacobian has some constant elements, they are wrecked
         ! by the scaling.  Restore them from gConU.

         if (NonlinearCon) then
            if (modefg .gt. 0  .and.  iw(conJij) .gt. 0) then
               call dcopy ( negCon, rw(lgConU), 1, gCon, 1 )
            end if
         end if
      end if

      !-----------------------------------------------------------------
      ! Compute the constraints.
      !-----------------------------------------------------------------
      ! To incorporate user workspace in fgcon, replace the next
      ! call to fgcon with:
      ! call fgcon ( modeC, nnCon, nnJac, nnCon,
      !&             iw(liy1), x, fCon, gCon, Status,
      !&             cu, lencu, iu, leniu, ru, lenru )

      if (NonlinearCon) then
         StatusUser = Status    ! In case fgCon alters Status

         if (lvlTim .ge. 2) call s1time( 4, 0, iw, leniw, rw, lenrw )
         call iload ( nnCon, (1), iw(liy1), 1 )
         call fgcon
     &      ( modeC, nnCon, nnJac, nnCon,
     &        iw(liy1), x, fCon, gCon, StatusUser )
         if (lvltim .ge. 2) call s1time(-4, 0, iw, leniw, rw, lenrw )

         iw(nfCon1) = iw(nfCon1) + 1
         if (modefg .gt. 0)
     &        iw(nfCon2) = iw(nfCon2) + 1
      end if

      !-----------------------------------------------------------------
      ! Compute the objective.
      !-----------------------------------------------------------------
      ! To incorporate user workspace in fgobj, replace the next
      ! call to fgobj with:
      ! call fgobj ( modeF, nnObjU, x, fObj, gObj, Status,
      !&             cu, lencu, iu, leniu, ru, lenru )

      if (NonlinearObj  .and.  modeC .ge. 0) then
         StatusUser = Status    ! In case fgObj alters Status

         if (lvlTim .ge. 2) call s1time( 5, 0, iw, leniw, rw, lenrw )

         if (FPonly) then
            call fgobj ( modeF, nnObjU, x, fObj, rw(lgSave), Status )
            call dcopy ( nnObj, x, 1, gObj, 1 )
            call daxpy ( nnObj, (-one), rw(lx0), 1, gObj, 1 )
            fObj = half*proxWeight*ddot ( nnObj, gObj, 1, gObj, 1 )
            call dscal ( nnObj, proxWeight, gObj, 1 )

         else ! nnObj = nnObjU
            call fgobj ( modeF, nnObjU, x, fObj, gObj, Status )
         end if

         if (lvlTim .ge. 2) call s1time(-5, 0, iw, leniw, rw, lenrw )
         iw(nfObj1) = iw(nfObj1) + 1
         if (modefg .gt. 0) iw(nfObj2) = iw(nfObj2) + 1
      end if

      !-----------------------------------------------------------------
      ! Scale  x and the derivatives.
      !-----------------------------------------------------------------
      if ( Scaled ) then
         call dcopy ( nnH, rw(lxScaled), 1, x, 1 )

         if (NonlinearCon) then
            call dddiv ( nnCon, rw(lscales+n), 1, fCon, 1 )
            if (modefg .gt. 0  .and.  iw(gotG) .gt. 0) then
               call s8scaleJ
     &            ( nnCon, nnJac, negCon, n, rw(lscales),
     &              neJ, nlocJ, locJ, indJ, gCon, rw, lenrw )
            end if
         end if

         if (NonlinearObj  .and.  modeC .ge. 0) then
            if (modefg .gt. 0  .and.  iw(gotG) .gt. 0) then
               call s8scaleg
     &            ( nnObj, rw(lscales), gObj, rw, lenrw )
            end if
         end if
      end if

      if (modeC .lt. 0  .or.  modeF .lt. 0) then
         !--------------------------------------------------------------
         ! The user may be saying the function is undefined (mode = -1)
         ! or may just want to stop                         (mode < -1).
         !--------------------------------------------------------------
         if (modeC .eq. -1  .or.  modeF .eq. -1) then
            iExit = -1
         else
            if (modeC .lt. 0) then
               iExit = 72
            else
               iExit = 73
            end if
         end if
      end if

      !=================================================================
      ! Do some housekeeping on the first entry.
      !=================================================================
      if (Status .eq. 1  .and.  iExit .eq. 0) then
         if (NonlinearCon) then
            !-----------------------------------------------------------
            ! Count how many Jacobian elements are provided.
            !-----------------------------------------------------------
            nnGlin   = 0
            supplied = 0
            do l = 1, negCon
               if (gCon(l) .ne. Gdummy) supplied = supplied + 1
            end do

            write(str, 1100) supplied, negCon
            call snPRNT( 3, str, iw, leniw )

            if (supplied .lt. negCon) then

               ! Some Jacobian elements are missing.

               if (iw(lvlDer) .ge. 2) then
                  !-----------------------------------------------------
                  ! All the Jacobian is known.  Any undefined elements
                  ! are assumed constant, and are restored from gConU.
                  !-----------------------------------------------------
                  write(str, 3100)
                  call snPRNT( 3, str, iw, leniw )

                  lG  = lgConU
                  do l  = 1, negCon
                     if (gCon(l) .eq. Gdummy) then
                        gCon(l) = rw(lG)
                        nnGlin  = nnGlin + 1
                     end if
                     lG = lG + 1
                  end do
               else
                  !-----------------------------------------------------
                  ! Save a permanent copy of gCon in gConU so that we
                  ! know which derivatives must be estimated.
                  !-----------------------------------------------------
                  call dcopy ( negCon, gCon, 1, rw(lgConU), 1 )
               end if
            end if ! supplied < negCon
            if (supplied + nnGlin .lt. negCon) iw(gotFD)  = 1
            if (supplied          .gt.      0) iw(gotG )  = 1
            if (nnGlin            .gt.      0) iw(conJij) = 1
         end if

         if (NonlinearObj) then
            !-----------------------------------------------------------
            ! Count how many gradient elements are known.
            ! (These may be the gradients of the FP objective.)
            !-----------------------------------------------------------
            if (FPonly) then
               supplied  = nnObj
               iw(gotG ) = 1
               write(str, 2010) nnObj
               call snPRNT( 3, str, iw, leniw )

            else
               supplied = 0
               do l = 1, nnObjU
                  if (gObj(l) .ne. Gdummy) then
                     supplied = supplied + 1
                  end if
               end do

               write(str, 2000) supplied, nnObjU
               call snPRNT( 3, str, iw, leniw )

               if (supplied .gt. 0) then
                  iw(gotG) = 1
               end if

               if (supplied .lt. nnObjU) then

!                 Some objective gradients are missing.

                  iw(gotFD) = 1

                  if (iw(lvlDer) .eq. 1  .or.  iw(lvlDer) .eq. 3) then
                     !--------------------------------------------------
                     ! The objective gradient was meant to be known.
                     !--------------------------------------------------
                     iw(lvlDer) = iw(lvlDer) - 1
                     write(str, 2100) iw(lvlDer)
                     call snPRNT( 3, str, iw, leniw )
                  end if
               end if
            end if

            !-----------------------------------------------------------
            ! Copy gObj into gObjU.
            !-----------------------------------------------------------
            call dcopy ( nnObj, gObj, 1, rw(lgObjU), 1 )

         end if
      end if

      return

 1100 format(  ' The user has defined', i8, '   out of', i8,
     &         '   constraint gradients.')
 2000 format(  ' The user has defined', i8, '   out of', i8,
     &         '   objective  gradients.')
 2010 format(  ' NpOpt   will define ', i8, '   gradients for the ',
     &         ' FP objective.')
 2100 format(  ' XXX  Some objective  derivatives are missing ---',
     &         ' derivative level reduced to', i3)
 3100 format(  ' ==>  Some constraint derivatives are missing, ',
     &         ' assumed constant.'/ )

      end ! subroutine s0fgN
