!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  snctrl.f  --- the Basic interface for SNCTRL.
!
!     snKerCT
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snKerCT
     &   ( Start, m, n, neJ, nNames,
     &     nnCon, nnObjU, nnJac,
     &     iObjU, objUAdd, Prob,
     &     s0fgctrl, odecon, algcon,
     &     snLog, snLog2, sqLog, snSTOP,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, obj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     s0fgctrl, odecon, algcon,
     &     snLog, snLog2, sqLog, snSTOP
      integer
     &     INFO, iObjU, lencu, lencw, leniu, leniw, lenru, lenrw, m,
     &     mincw, miniw, minrw, n, neJ, nInf, nNames, nnCon, nnJac,
     &     nnObjU, nS, hs(n+m), indJ(neJ), iu(leniu), iw(leniw),
     &     locJ(n+1)
      double precision
     &     obj, objUAdd, sInf, Jcol(neJ), bl(n+m), bu(n+m), pi(m),
     &     rc(n+m), ru(lenru), rw(lenrw), x(n+m)
      character*(*)
     &     Start
      character
     &     Prob*8, cu(lencu)*8, cw(lencw)*8, Names(nNames)*8
!     ==================================================================
!     snKerCT does the work for SNCTRL. (Kernel for snctrl)
!
!     Developers can call this version with customized versions of
!     snLog, snLog2  and  snSTOP.
!
!     17 Oct 2004: First version of snKerB.
!     01 Sep 2007: Sticky parameters added.
!     04 Jul 2010: mincw, miniw, minrw added to workspace.
!     11 Sep 2014: nnObjU and nnObj separated for FP mode.
!     29 Apr 2015: Control version taken from SNOPTB.
!     ==================================================================
      character
     &     Solver*6, str*80, str2*80
      integer
     &     Useriw(130)
      double precision
     &     Userrw(130)
      character
     &     Usercw(130)*8
      logical
     &     FPonly, gotR, PrintMem
      integer
     &     Errors, HDInfo, inform, iObj, lenR, lenx0, lgObj, leType,
     &     lkx, llocG, liwEst, lrwEst, lvlHes, lvlStart, lx0, maxcw,
     &     maxiw, maxR, maxrw, maxS, minmax, mProb, mQNmod, nb, negCon,
     &     nextcw, nextiw, nextrw, ngQP, nInfE, nkx, nlocG, nlocJ,
     &     nMajor, nnObj, nnObj0, nnH, nnH0, nrhs, nrhs0, nx0,
     &     startType, stickyOp
      integer
     &     neH, nlocH, indH(1), locH(1)
      double precision
     &     Hcol(1)
       double precision
     &     fObj, objAdd, objTrue, sInfE, rhs(1), x0(1)
      external
     &     s8HxLP, s8HxNP, s8HxQP
!     ------------------------------------------------------------------
      integer            COLD,       BASIS,      WARM,       HOT
      parameter         (COLD   = 0, BASIS  = 1, WARM   = 2, HOT    = 3)
      integer            StdIn
      parameter         (StdIn  = 2)
      integer            Unset,      Normal
      parameter         (Unset  =-1, Normal = 0)

      parameter         (mProb     =  51) ! Problem name
      parameter         (lvlStart  =  69) ! cold:warm:basis:hot start
      parameter         (HDInfo    = 243) ! Approximate Hessian type

      double precision   zero
      parameter         (zero      = 0.0d+0)
!     ------------------------------------------------------------------
      Solver = 'SNCTRL'
      INFO   = 0

!     ------------------------------------------------------------------
!     Check memory limits and fetch the workspace starting positions.
!     ------------------------------------------------------------------
      call s2Mem0
     &   ( INFO, Solver, lencw, leniw, lenrw, iw,
     &     mincw, miniw, minrw, maxcw, maxiw, maxrw,
     &     nextcw, nextiw, nextrw )
      if (INFO .gt. 0) go to 999 ! Exit without printing

      cw(1)  = Solver//'  '

!     Save the user's option choices  (weird choices get overwritten).
!     Initialize timers and the standard input file.

      call chcopy( 130, cw(51), 1, Usercw, 1 )
      call icopy ( 130, iw(51), 1, Useriw, 1 )
      call dcopy ( 130, rw(51), 1, Userrw, 1 )

      call s1time( 0, 0, iw, leniw, rw, lenrw  )
      call s1file( StdIn, iw, leniw )

!     Check the input arguments.

      startType = iw(lvlStart)
      call s3chkArgsB
     &   ( inform, Start, m, n, neJ, nNames, nS,
     &     nnCon, nnObjU, nnJac, iObjU,
     &     indJ, locJ, bl, bu, Names, hs, pi, startType, Errors,
     &     iw, leniw, rw, lenrw )
      if (inform .gt. 0) then
         INFO = inform
         go to 800
      end if

!     ------------------------------------------------------------------
!     The obligatory call to snInit has already set the defaults.
!     Check that the optional parameters have sensible values.
!     Print the options.
!     ------------------------------------------------------------------
      cw(mProb) = Prob

      call s8Defaults
     &   ( m, n, nnCon, nnJac, nnObjU, iObjU,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call s3printB
     &   ( m, n, nnCon, nnJac, nnObjU, startType, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Compute the storage requirements for SNOPT  from the following
!     variables:
!         m,      n,     ne
!         lenR  , maxS
!         nnCon , nnJac, nnObjU,
!         negCon
!     All have to be known before calling s2Mem.
!     The only one in doubt is negCon, the number of Jacobian elements.
!     Count them in s8Gsize.
!     ------------------------------------------------------------------
      nb      = n + m
      nlocJ   = n + 1
      nkx     = nb

      call s8Gsize
     &   ( m, nnCon, nnJac, neJ, nlocJ, locJ, indJ, negCon )

      maxR    = iw( 52) ! max columns of R.
      maxS    = iw( 53) ! max # of superbasics
      mQNmod  = iw( 54) ! (ge 0) max # of BFGS updates
      lvlHes  = iw( 72) ! 0, 1,  2  => LM, FM, Exact Hessian
      minmax  = iw( 87) ! 1, 0, -1  => MIN, FP, MAX

      lenR    = maxR*(maxR + 1)/2  +  (maxS - maxR)

      iw( 20) = negCon  ! # of nonzeros in gCon
      iw( 28) = lenR    ! R(lenR) is the reduced Hessian factor

      nnH     = max( nnJac, nnObjU )

      neH     = 1               ! Placeholders
      nlocH   = 1

!     Load iw with various problem dimensions.

      iw( 15) = n       ! copy of the number of columns
      iw( 16) = m       ! copy of the number of rows
      iw( 17) = neJ     ! copy of the number of nonzeros in Jcol
      iw( 21) = nnJac   ! # of Jacobian  variables
      iw( 22) = nnObjU  ! # of user-defined objective variables
      iw( 23) = nnCon   ! # of nonlinear constraints
      iw( 24) = nnH     !   max( nnObj, nnJac )
      iw(204) = iObjU   ! position of the objective row in J

!     ------------------------------------------------------------------
!     If only a feasible point is requested, save the base point for the
!     objective function:  1/2 || x - x0 ||^2
!
!     Set the working objective gradient dimensions.
!     ------------------------------------------------------------------
      FPonly  = minmax .eq. 0
      if (FPonly) then
         nnObj  = nnH
         iObj   = 0
         objAdd = zero
      else
         nnObj  = nnObjU       ! working # nonlinear objective vars
         iObj   = iObjU
         objAdd = objUAdd
      end if

      nnObj0 = max( nnObj, 1 )

!     ------------------------------------------------------------------
!     Allocate the local arrays for snOpt.
!     s8Map  maps snOpt integer and double arrays.
!     s2BMap maps the arrays for the LU routines.
!     s2Mem  checks what space is available and prints any messages.
!     ------------------------------------------------------------------
      call s8Map
     &   ( m, n, negCon, nkx, nnCon, nnJac, nnObjU, nnObj, nnH,
     &     lenR, maxR, maxS,  mQNmod, lvlHes,
     &     nextcw, nextiw, nextrw, iw, leniw )
      call s2Bmap
     &   ( m, n, neJ, maxS,
     &     nextiw, nextrw, maxiw, maxrw, liwEst, lrwEst, iw, leniw )
      PrintMem = .true.         ! Print all messages in s2Mem
      call s2Mem
     &   ( inform, PrintMem, liwEst, lrwEst,
     &     nextcw, nextiw, nextrw,
     &     maxcw, maxiw, maxrw, lencw, leniw, lenrw,
     &     mincw, miniw, minrw, iw )
      if (inform .ne. 0) then
         INFO = inform
         go to 800
      end if

!     ------------------------------------------------------------------
!     If only a feasible point is requested, save the base point for the
!     objective function:  1/2 || x - x0 ||^2
!     ------------------------------------------------------------------
      if (FPonly) then
         lx0 = iw(298)
         call dcopy ( nnH, x, 1, rw(lx0), 1 )
      end if

!     Define the row and column ordering for J.
!     SNOPT  uses natural order throughout, so kx = kxN.

      iw(247) = nkx     ! dimension of kx and its inverse, kxN
      lkx     = iw(251) ! j  = kx (jN) => col j of Jcol is variable jN
      iw(252) = lkx     ! jN = kxN(j ) => col j of Jcol is variable jN

      call s1perm( n, iw(lkx) )
      call s1perm( m, iw(lkx+n) )

!     ------------------------------------------------------------------
!     Construct column pointers for the nonlinear part of the  Jacobian.
!     ------------------------------------------------------------------
      if (nnCon .gt. 0) then
         llocG = iw(260) ! locG(nlocG) = column pointers for indG
         nlocG = nnJac + 1

         call s8Gloc
     &      ( nnCon, nnJac,
     &        neJ, nlocJ, locJ, indJ, negCon, nlocG, iw(llocG) )
      end if

!     ------------------------------------------------------------------
!     Solve the problem.
!     ------------------------------------------------------------------
      if (nnH .eq. 0) then

!        The problem is a linear program.
!        nnObjU = nnJac = 0

         nrhs   = 0             ! No constraint rhs vector.
         nx0    = 0             ! No constant shift for x.
         lenx0  = 1
         nrhs0  = 1
         nnH0   = 1
         ngQP   = nnObj

         leType = iw(283)       ! eType(nb) definition of elastic vars
         lgObj  = iw(297)       ! gObj(ngObj) = Objective gradient

         call iload ( nb, 3, iw(leType), 1 )

         call s5solve
     &      ( INFO,
     &        Solver, startType,
     &        sqLog, s8HxQP, s8HxLP, gotR,
     &        m, n, nb, nnH0, nnH,
     &        nNames, ngQP, nnObj0, nnObj,
     &        iObj, objadd, fObj, objTrue,
     &        nInf, sInf, nInfE, sInfE,
     &        neJ, nlocJ, locJ, indJ, Jcol,
     &        neH, nlocH, locH, indH, Hcol,
     &        bl, bu, rw(lgObj), Names,
     &        nrhs0, nrhs, rhs, lenx0, nx0, x0,
     &        iw(leType), hs, x, pi, rc, nS,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )

      else

!        The problem is nonlinear.
!        Define the type of initial Hessian.

         if      (startType .eq. COLD ) then
            iw(HDInfo) = Unset
         else if (startType .eq. BASIS) then
            iw(HDInfo) = Unset
         else if (startType .eq. WARM ) then
            iw(HDInfo) = Unset
         else if (startType .eq. HOT  ) then
            iw(HDInfo) = Normal
         end if

         call s8solve
     &      ( INFO, Solver, startType,
     &        s0fgctrl, odecon, algcon, s8HxNP,
     &        snLog, snLog2, snSTOP, gotR,
     &        m, n, nb, nnH, nnCon, nnJac, nnObj,
     &        nNames, iObj, objadd, fObj, objTrue,
     &        nInf, sInf, nInfE, sInfE,
     &        neJ, nlocJ, locJ, indJ, Jcol,
     &        neH, nlocH, locH, indH, Hcol,
     &        bl, bu, Names,
     &        hs, x, pi, rc, nMajor, nS,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
      end if

      nInf   = nInf + nInfE
      sInf   = sInf + sInfE
      obj    = fObj

      mincw  = iw(47)            ! minimum length of cw
      miniw  = iw(48)            ! minimum length of iw
      minrw  = iw(49)            ! minimum length of rw

!     If "sticky parameters no",  restore the user-defined options

      stickyOp = iw(116)

      if (stickyOp .le. 0) then
         call chcopy( 130, Usercw, 1, cw(51), 1 )
         call icopy ( 130, Useriw, 1, iw(51), 1 )
         call dcopy ( 130, Userrw, 1, rw(51), 1 )
      end if

!     Print times for all clocks (if lvlTim > 0).

      call s1time( 0, 2, iw, leniw, rw, lenrw )

      return

!     Local exit messages.

  800 call snWRAP( INFO, Solver, str, str2, iw, leniw )

  999 return

      end ! subroutine snKerB

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

