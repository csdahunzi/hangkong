!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn04stats.f
!
!     s4qpGetstats   s4npGetStats
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4qpGetStats
     &   ( m, n, nnH,
     &     probName,
     &     INFO,
     &     nS, nInf, sInf, iObj, objAdd, fobj, x,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     INFO, iObj, lencw, leniw, lenrw, m, n, nInf,
     &     nnH, nS, iw(leniw)
      double precision
     &     fobj, objAdd, sInf, x(n+m), rw(lenrw)
      character
     &     ProbName*8, cw(lencw)*8

!     ==================================================================
!     s4qpGetStats prints statistics associated with an SQOPT run.
!
!     04 May 2015: Print overflow fixed for  itns, HvCalls.
!     ==================================================================
      character
     &     flag*3, msg*43,  String*132
      logical
     &     FileExists, FPonly
      integer
     &     iAll, iPP, iStats, iTeX
      integer
     &     HvCalls, iflag, itns, length, minmax
      double precision
     &     objTrue, rHvCalls, runTime, zero
!     ------------------------------------------------------------------
      parameter         (zero    = 0.0d+0)

      parameter         (minmax  = 87) ! 1, 0, -1  => MIN, FP, MAX
      parameter         (iAll    = 7, iTeX =  8, iStats = 9, iPP = 11)
!     ------------------------------------------------------------------

      runTime  = rw(462) ! Solve time
      itns     = iw(421) ! Total iteration count
      HvCalls  = iw(188) ! Hessian- vector products

      FPonly   = iw(minmax) .eq. 0

!     Calls for qpStatus ge 2 are not included in the totals.

      HvCalls = max(HvCalls,0)
      if (HvCalls .gt. 0) then
         HvCalls =  HvCalls - 1
      end if

      itns    = mod( itns   , 1000000 )
      HvCalls = mod( HvCalls, 1000000 )

      if (FPonly) then
         objTrue = zero
      else
         objTrue = objAdd + fobj
         if (iObj .gt. 0) then
            objTrue = objTrue + x(n+iObj)
         end if
      end if

      flag = ' '
      msg  = 'Unknown exit'

      if (     INFO .eq.  1) then
         flag  = '   '
         iflag = 1
         msg   = 'Optimal solution found'
      else if (INFO .eq.  2) then
         flag  = 'fp '
         iflag = 2
         msg   = 'Feasible point found'

      else if (INFO .eq.  3) then
         flag  = 'acc'
         iflag = 3
         msg   = 'Optimal, but accuracy not achvd'

      else if (INFO .eq.  6) then
         msg   = 'Elastic infeas. minimized'
         flag  = 'Inf'
         iflag = 11

      else if (INFO .eq. 11) then
         flag  = 'Inf'
         iflag = 11
         msg   = 'infeasible linear constraints'

      else if (INFO .eq. 12) then
         flag  = 'Lnf'
         iflag = 11
         msg   = 'infeasible linear equalities'

      else if (INFO .eq. 13) then
         flag  = 'Inf'
         iflag = 13
         msg   = 'nonlinear infeasibilities minimized'

      else if (INFO .eq. 14) then
         flag  = 'Lnf'
         iflag = 11
         msg   = 'linear infeasibilities minimized'

      else if (INFO .eq. 15) then
         flag  = 'Lnf'
         iflag = 11
         msg   = 'infeasible QP subproblem'

      else if (INFO .eq. 16) then
         flag  = 'Lnf'
         iflag = 11
         msg   = 'infeasible nonelastics'

      else if (INFO .eq. 21) then
         flag  = 'Unb'
         iflag = 21
         msg   = 'Unbounded problem'

      else if (INFO .eq. 22) then
         flag   = 'vlm'
         iflag = 22
         msg    = 'Violation limit exceeded'

      else if (INFO .eq. 31) then
         flag  = 'Itr'
         iflag = 31
         msg   = 'Iteration limit exceeded'

      else if (INFO .eq. 32) then
         flag  = 'Itr'
         iflag = 32
         msg   = 'Major iterations limit'

      else if (INFO .eq. 33) then
         flag  = 'SB '
         iflag = 33
         msg   = 'Superbasics limit too small'

      else if (INFO .eq. 34) then
         flag  = 'Cpu'
         iflag = 34
         msg   = 'Cpu time limit exceeded'

      else if (INFO .eq. 41) then
         flag  = 'Cbi'
         iflag = 41
         msg   = 'Current point cannot be improved'
*                 1234567890123456789012345678901234567890123
*                           1         2         3         4
      else if (INFO .eq. 42) then
         flag  = 'Gnf'
         iflag = 42
         msg   = 'Singular basis'

      else if (INFO .eq. 43) then
         flag  = 'Gnf'
         iflag = 43
         msg   = 'Constraints cannot be satisfied'

      else if (INFO .eq. 44) then
         flag  = 'Gnf'
         iflag = 44
         msg   = 'Ill-conditioned null-space basis'

      else if (INFO .eq. 53) then
         flag  = 'Ind'
         iflag = 53
         msg   = 'QP Hessian is indefinite'

      else if (INFO .gt. 0) then
         flag  = '???'
         iflag = 99
         msg   = ' '
      end if

!     ------------------------------------------------------------------
!     Write the summary of the problem statistics
!     ------------------------------------------------------------------
      inquire( FILE='temp.stats', EXIST=FileExists )
      if (FileExists) then
         open (iStats,  file='temp.stats',  status='old',
     &                                    position='append')
      else
         open (iStats,  file='temp.stats',  status='new')
         write(iStats, 2000)
      end if

      write(String, 2100) probName, n, m, nnH, nS, objTrue
      call s1trim( String, length )
      write(iStats,7000) String(1:length)

      close(iStats)

!     ------------------------------------------------------------------
!     Tex summary
!     ------------------------------------------------------------------
      inquire( FILE='temp.tex', EXIST=FileExists )
      if (FileExists) then
         open (iTeX,  file='temp.tex',      status='old',
     &                                    position='append')
      else
         open (iTeX,  file='temp.tex',      status='new')
         write(iTeX, 3000)
      end if

      write(String, 3100) probName, itns, HvCalls, objTrue,
     &                    runTime, flag
      call s1trim( String, length )
      write(iTeX,7000) String(1:length)

      close(iTeX)

!     ------------------------------------------------------------------
!     All data
!     ------------------------------------------------------------------
      inquire( FILE='temp.all', EXIST=FileExists )
      if (FileExists) then
         open (iAll,  file='temp.all',  status='old',
     &                                  position='append')
      else
         open (iAll,  file='temp.all',  status='new')
         write(iAll, 4000)
      end if

      write(String, 4100) probName, n, m, itns, HvCalls,
     &                    runTime, objTrue, nInf, msg
      call s1trim( String, length )
      write(iAll,7000) String(1:length)

      close(iAll)

!     ------------------------------------------------------------------
!     Performance profile
!     ------------------------------------------------------------------
      inquire( FILE='temp.pp', EXIST=FileExists )
      if (FileExists) then
         open (iPP,  file='temp.pp',    status='old',
     &                                position='append')
      else
         open (iPP,  file='temp.pp',    status='new')
         write(iPP, 6000)
      end if

      RunTime  = max(1.0d-4, RunTime)
      rHvCalls = max(1.0d-3, dble(HvCalls) )
      write(String, 6100) probName, iflag, RunTime, rHvCalls
      if (INFO .le. 14 .or. INFO .eq. 21 .or. INFO .eq. 13) then
!        Relax
      else
         String(16:27) = ' '
         String(31:42) = ' '
         String(25:27) = 'NaN'
         String(40:42) = 'NaN'
      end if

      call s1trim( String, length )
      write(iPP,7000) String(1:length)

      close (iPP)

      return

 2000 format( 28x, 'Linear Nonlinear'
     &       /'Problem Name', 3x, 'Variables constrnts variables',
     &       7x, 'SBs', 3x, 'Final objective' )
 2100 format(1x, a8, 5x, 4i10, 1p, e18.6)

 3000 format(' \\begin{tabular}{lrrrl}' /
     &   ' Problem Name&      Itns',
     &   ' &  H prods&',
     &   '    Objective   &     Time& Status \\\\[2ex]')
 3100 format( 1x, a12, '&', i10, ' &', i8, ' &$',
     &        1p, e13.6, ' $&$', 0p, f8.2, '&{\em ',
     &        a3, '} \\' )

!3200 format(' \\end{verbatim}'
!    &    // ' \\end{document}')

 4000 format(' Name', 11x, 'n', 6x, 'm', 3x, 'Itns', 5x, 'Hv', 4x,
     &       'cpu(s)', 6x, 'Objective', 6x, 'inf', 2x,  'Result')
 4100 format( a8, 2x, 2i7, 2i7, 1x, f9.2, 1p, 2x, e16.9, i6, 2x, a)

 6000 format('ProName', 4x, 'Info', 7x, 'Time', 8x, 'Hv calls' )
 6100 format(a8, 3x, i2, 2x, f12.4, 3x, f12.3)

 7000 format( a )

      end ! s4qpGetStats

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s4npGetStats
     &   ( m, n,
     &     nnCon, nnObj, nnJac,
     &     probName,
     &     INFO,
     &     nS, nInf, sInf, iObj, objAdd, fobj, x,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     INFO, iObj, lencw, leniw, lenrw, m, n, nInf, nnCon, nnJac,
     &     nnObj, nS, iw(leniw)
      double precision
     &     fobj, objAdd, sInf, x(n+m), rw(lenrw)
      character
     &     ProbName*8, cw(lencw)*8

!     ==================================================================
!     s4npGetStats prints statistics associated with an SNOPT run.
!
!     04 May 2015: Print overflow fixed for  itns, majors, nFuns, nCons.
!     ==================================================================
      character
     &     flag*3, msg*43,  String*132
      logical
     &     FileExists, FPonly
      integer
     &     iAll, iPP, iStats, iTeX
      integer
     &     iflag, itns, length, linCon, majors, minmax, nnL,
     &     nCons, nFuns, nObjs
      double precision
     &     maxVi, objTrue, rnObjs, runTime, zero
!     ------------------------------------------------------------------
      parameter         (zero    = 0.0d+0)

      parameter         (minmax  = 87) ! 1, 0, -1  => MIN, FP, MAX
      parameter         (iAll    = 7, iTeX =  8, iStats = 9, iPP = 11)
!     ------------------------------------------------------------------

      maxVi    = rw(432) ! Inf norm of the constraint violation
      runTime  = rw(462) ! Solve time

      nCons    = iw(190) ! calls to fCon: mode > 0
      nFuns    = iw(195) ! calls to fObj: mode > 0

      itns     = iw(421) ! Total iteration count
      majors   = iw(422) ! Major iterations

      FPonly   = iw(minmax) .eq. 0

      nnL      = max( nnObj, nnJac )

      linCon   = m - nnCon
!     linVar   = n - nnL ! Not used right now

      if (nnL .eq. 0) then
         majors = 0
         maxVi  = zero
      end if

!     nFuns, nCons are undefined for an LP.
!     Calls for npStatus ge 2 are not included in the totals.

      nFuns = max(nFuns,0)
      nCons = max(nCons,0)
      if (nFuns .gt. 0) then
         nFuns = nFuns - 1
      end if

      if (nCons .gt. 0) then
         nCons = nCons - 1
      end if

      itns    = mod( itns  , 1000000 )
      majors  = mod( majors, 1000000 )

      nFuns   = mod( nFuns , 1000000 )
      nCons   = mod( nCons , 1000000 )

      nObjs    = max(nCons,nFuns)

      if (FPonly) then
         objTrue = zero
      else
         objTrue = objAdd + fobj
         if (iObj .gt. 0) then
            objTrue = objTrue + x(n+iObj)
         end if
      end if

      flag = ' '
      msg  = 'Unknown exit'

      if (     INFO .eq.  1) then
         flag  = '   '
         iflag = 1
         msg   = 'Optimal solution found'
      else if (INFO .eq.  2) then
         flag  = 'fp '
         iflag = 2
         msg   = 'Feasible point found'

      else if (INFO .eq.  3) then
         flag  = 'acc'
         iflag = 3
         msg   = 'Optimal, but accuracy not achvd'

      else if (INFO .eq.  6) then
         msg   = 'Elastic infeas. minimized'
         if (nnCon .eq. 0) then
            flag  = 'Lnf'
            iflag = 11
         else
            flag  = 'Inf'
            iflag = 13
         end if

      else if (INFO .eq. 11) then
         flag  = 'Lnf'
         iflag = 11
         msg   = 'infeasible linear constraints'

      else if (INFO .eq. 12) then
         flag  = 'Lnf'
         iflag = 11
         msg   = 'infeasible linear equalities'

      else if (INFO .eq. 13) then
         flag  = 'Inf'
         iflag = 13
         msg   = 'nonlinear infeasibilities minimized'

      else if (INFO .eq. 14) then
         flag  = 'Lnf'
         iflag = 11
         msg   = 'linear infeasibilities minimized'

      else if (INFO .eq. 15) then
         flag  = 'Lnf'
         iflag = 11
         msg   = 'infeasible QP subproblem'

      else if (INFO .eq. 16) then
         flag  = 'Lnf'
         iflag = 11
         msg   = 'infeasible nonelastics'

      else if (INFO .eq. 21) then
         flag  = 'Unb'
         iflag = 21
         msg   = 'Unbounded problem'

      else if (INFO .eq. 22) then
         flag   = 'vlm'
         iflag = 22
         msg    = 'Violation limit exceeded'

      else if (INFO .eq. 31) then
         flag  = 'Itr'
         iflag = 31
         msg   = 'Iteration limit exceeded'

      else if (INFO .eq. 32) then
         flag  = 'Itr'
         iflag = 32
         msg   = 'Major iterations limit'

      else if (INFO .eq. 33) then
         flag  = 'SB '
         iflag = 33
         msg   = 'Superbasics limit too small'

      else if (INFO .eq. 34) then
         flag  = 'Cpu'
         iflag = 34
         msg   = 'Cpu time limit exceeded'

      else if (INFO .eq. 41) then
         flag  = 'Cbi'
         iflag = 41
         msg   = 'Current point cannot be improved'

      else if (INFO .eq. 42) then
         flag  = 'Gnf'
         iflag = 42
         msg   = 'Singular basis'

      else if (INFO .eq. 43) then
         flag  = 'Gnf'
         iflag = 43
         msg   = 'Constraints cannot be satisfied'

      else if (INFO .eq. 44) then
         flag  = 'Gnf'
         iflag = 44
         msg   = 'Ill-conditioned null-space basis'

      else if (INFO .gt. 0) then
         flag  = '???'
         iflag = 99
         msg   = ' '
      end if

!     ------------------------------------------------------------------
!     Write the summary of the problem statistics
!     ------------------------------------------------------------------
      inquire( FILE='temp.stats', EXIST=FileExists )
      if (FileExists) then
         open (iStats,  file='temp.stats',  status='old',
     &                                    position='append')
      else
         open (iStats,  file='temp.stats',  status='new')
         write(iStats, 2000)
      end if

      write(String, 2100) probName, n, linCon, nnCon, nS, objTrue
      call s1trim( String, length )
      write(iStats,7000) String(1:length)

      close(iStats)

!     ------------------------------------------------------------------
!     Tex summary
!     ------------------------------------------------------------------
      inquire( FILE='temp.tex', EXIST=FileExists )
      if (FileExists) then
         open (iTeX,  file='temp.tex',      status='old',
     &                                    position='append')
      else
         open (iTeX,  file='temp.tex',      status='new')
         write(iTeX, 3000)
      end if

      write(String, 3100) probName, itns, majors, nObjs, objTrue, maxVi,
     &                    runTime, flag
      call s1trim( String, length )
      write(iTeX,7000) String(1:length)

      close(iTeX)

!     ------------------------------------------------------------------
!     All data
!     ------------------------------------------------------------------
      inquire( FILE='temp.all', EXIST=FileExists )
      if (FileExists) then
         open (iAll,  file='temp.all',  status='old',
     &                                  position='append')
      else
         open (iAll,  file='temp.all',  status='new')
         write(iAll, 4000)
      end if

      write(String, 4100) probName, n, m, itns, majors, nFuns, nCons,
     &                    runTime, objTrue, nInf, msg
      call s1trim( String, length )
      write(iAll,7000) String(1:length)

      close(iAll)

!     ------------------------------------------------------------------
!     Performance profile
!
!     Assume that a locally infeasible point counts as a "success".
!     ------------------------------------------------------------------
      inquire( FILE='temp.pp', EXIST=FileExists )
      if (FileExists) then
         open (iPP,  file='temp.pp',    status='old',
     &                                position='append')
      else
         open (iPP,  file='temp.pp',    status='new')
         write(iPP, 6000)
      end if

      RunTime = max(1.0d-4, RunTime)
      rnObjs  = dble(nObjs)
      rnObjs  = max(1.0d-3,  rnObjs )
      write(String, 6100) probName, iflag, RunTime, rnObjs
      if (INFO .le. 14 .or. INFO .eq. 21 .or. INFO .eq. 13) then
!        Relax
      else
         String(16:27) = ' '
         String(31:42) = ' '
         String(25:27) = 'NaN'
         String(40:42) = 'NaN'
      end if

      call s1trim( String, length )
      write(iPP,7000) String(1:length)

      close (iPP)

      return

 2000 format( 28x, 'Linear Nonlinear'
     &       /'Problem Name', 3x, 'Variables constrnts constrnts',
     &       7x, 'SBs', 3x, 'Final objective' )
 2100 format(1x, a8, 5x, 4i10, 1p, e18.6)

 3000 format(' \\begin{tabular}{lrrrrrrl}' /
     &   ' Problem Name&       Min',
     &   ' &     Maj &     Fun &',
     &   '    Objective   &  norm(c)   &    Time& Status \\\\[2ex]')
 3100 format( 1x, a12, '&', i10, 2(' &',i8), ' &$',
     &        1p, e13.6, ' $&$', e9.2 , ' $&', 0p, f8.2, '&{\em ',
     &        a3, '} \\' )

!3200 format(' \\end{tabular}'
!    &    // ' \\begin{verbatim}')

 4000 format(' Name', 11x, 'n', 6x, 'm', 3x, 'Minr', 3x, 'Majr', 3x,
     &       'nFun', 3x, 'Ncon', 4x, 'cpu(s)', 6x, 'Objective', 6x,
     &       'inf', 2x,  'Result')
 4100 format( a8, 2x, 2i7, 2i7, 2i7, 1x, f9.2, 1p, 2x, e16.9, i6, 2x, a)

 6000 format('ProName', 4x, 'Info', 7x, 'Time', 11x, 'Funs' )
 6100 format(a8, 3x, i2, 2x, f12.4, 3x, f12.3)
 7000 format( a )

      end ! s4npGetStats

