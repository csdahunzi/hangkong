!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  snoptc.f --- snOpt interface with Combined obj/constraints.
!
!     snOptC   snKerC
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snOptC
     &   ( Start, m, n, neJ, nNames,
     &     nnCon, nnObjU, nnJac,
     &     iObjU, objUAdd, Prob,
     &     userfun,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     userfun
      integer
     &     INFO, iObjU, lencu, leniu, lenru, lencw, leniw, lenrw,
     &     mincw, miniw, minrw, m, n, neJ, nNames, nS, nInf, nnCon,
     &     nnObjU, nnJac, indJ(neJ), hs(n+m), locJ(n+1), iu(leniu),
     &     iw(leniw)
      double precision
     &     sInf, Obj, objUAdd, Jcol(neJ), bl(n+m), bu(n+m), x(n+m),
     &     pi(m), rc(n+m), ru(lenru), rw(lenrw)
      character*(*)
     &     Start
      character
     &     Prob*8, Names(nNames)*8, cu(lencu)*8, cw(lencw)*8

!     ------------------------------------------------------------------
!     snOptC  is a Fortran subroutine for constrained nonlinear
!     optimization.  The constraints take the form
!
!                            (   x  )
!                      bl <= (      ) <= bu,
!                            ( F(x) )
!
!     where bl and bu are constant lower and upper bounds.
!
!     o If all constraints are linear, F = J x for some sparse matrix J.
!
!     o If all constraints are nonlinear, F = fCon(x) for some vector
!       fCon of smooth functions.
!
!     o In general, there is a mixture of constraints of the form
!                      ( fCon(x1) +  J2 x2 ),
!                      (   J3 x1  +  J4 x2 )
!       where the nonlinear variables x1 must appear first as shown.
!
!     o fCon(x1) and (optionally) its partial derivatives J1(x) are set
!       in subroutine userfun (see below).
!
!     o The matrices J2, J3, J4 and the sparsity pattern of J1(x) are
!       entered column-wise in the arrays Jcol, indJ, locJ (below).
!
!     o Internally, the constraints are converted into the form
!
!           fCon(x1) +  J2 x2  - s1      = 0,     bl <= ( x ) <= bu
!             J3 x1  +  J4 x2       - s2 = 0            ( s )
!
!       where s = (s1,s2)  and the components of (x,s) are the
!       variables and slacks respectively.
!
!     ------------------------------------------------------------------
!     NOTE: Before calling snOptC, your calling program must call:
!     call snInit( iPrint, iSumm,
!    &             cw, lencw, iw, leniw, rw, lenrw )
!     This sets the default values of the optional parameters. You can
!     also alter the default values of iPrint and iSumm before snOptC
!     is used.  iPrint = 0, etc, is OK.
!     ------------------------------------------------------------------
!
!     ON ENTRY:
!
!     Start   specifies how a starting basis (and certain other items)
!             are to be obtained.
!             Start = 'Cold' means that Crash should be used to choose
!                      an initial basis, unless a basis file is given
!                      via Old basis, Insert or Load in the Specs file.
!             Start = 'Basis file' means the same (but is more
!                      meaningful in the latter case).
!             Start = 'Warm' means that a basis is already defined in hs
!                      (probably from an earlier call).
!
!     m       is the number of slacks (i.e., general constraints).
!             For LP, QP or LC  problems this means the number of rows
!             in the constraint matrix J.
!             m > 0.
!
!             For problems with no general constraints, set m = 1 and
!             impose the constraint that the sum of the variables
!             must lie between plus and minus infinity. This gives
!             J  one ``free row'' that will not alter the solution.
!
!     n       is the number of variables, excluding slacks.
!             For LP problems, this is the number of columns in J.
!             n > 0.
!
!     neJ     is the number of nonzero entries in J (including the
!             Jacobian for any nonlinear constraints).
!             neJ gt 0.
!
!     nNames  is the number of column and row names provided in the
!             array  Names.  If nNames = 1, there are NO names.
!             Generic names will be used in the printed solution.
!             Otherwise, nNames = n+m and all names must be provided.
!
!     nnCon   is the number of nonlinear constraints.
!             nnCon ge 0.
!
!             a nonzero nnCon defines the row dimension of the
!             constraint Jacobian J1(x) defined in subroutine userfun.
!
!     nnObjU   is the number of nonlinear Objective variables.
!             nnObjU ge 0.
!
!     nnJac   is the number of nonlinear Jacobian variables.
!             If nnCon = 0, nnJac = 0.
!             if nnCon > 0, nnJac > 0.
!
!             a nonzero nnJac defines the column dimension of the
!             constraint Jacobian J1(x) defined in subroutine userfun.
!
!     iObjU    says which row of J is a free row containing a linear
!             objective vector  c  (iObjU = 0 if none).
!             iObjU = 0  or  nnCon < iObjU le m.
!
!     objUAdd  is a constant that will be added to the objective.
!             Typically objUAdd = 0.0d+0.
!
!     Prob    is an 8-character name for the problem, used in the
!             output.  A blank name can be assigned if necessary.
!
!     Jcol(neJ) is the constraint Jacobian J, stored column-wise.  Every
!             element of Jcol(*) must be assigned a value.  Elements
!             in the nonlinear part (see NOTE 2 below) can be any dummy
!             value (e.g., zero) since they are initialized by snOptC at
!             the first point that is feasible with respect to the
!             linear constraints.  The linear part of Jcol(*) must
!             contain the constant Jacobian elements.
!
!     indJ(neJ)  is the list of row indices for each nonzero in Jcol(*).
!
!     locJ(n+1) is a set of pointers to the beginning of each column of
!             the constraint matrix within Jcol(*) and indJ(*).
!             Must have locJ(1) = 1 and locJ(n+1) = ne+1.
!
!  NOTES:  1. If the problem has a nonlinear objective, the first nnObjU
!             columns of Jcol and indJ belong to the nonlinear objective
!             variables. Subroutine userfun deals with these variables.
!
!          2. If the problem has nonlinear constraints, the first nnJac
!             columns of Jcol and indJ belong to the nonlinear Jacobian
!             variables, and the first nnCon rows of Jcol and indJ
!             belong to the nonlinear constraints.  Subroutine userfun
!             deals with these variables and constraints.
!
!          3. If nnObjU > 0 and nnJac > 0, the two sets of
!             nonlinear variables overlap.  The total number of
!             nonlinear variables is nnL = max( nnObjU, nnJac ).
!
!          4. The Jacobian forms the top left corner of Jcol and indJ.
!             If a Jacobian column j (1 le j le nnJac) contains
!             any entries Jcol(k), indJ(k) associated with nonlinear
!             constraints (1 le indJ(k) le nnCon), those entries must
!             come before any other (linear) entries.
!
!          5. The row indices indJ(k) for a column may be in any order
!             (subject to Jacobian entries appearing first).
!             Subroutine userfun must define Jacobian entries in the
!             same order.
!
!     bl(n+m) is the lower bounds on each variable (x,s).
!
!     bu(n+m) is the upper bounds on each variable (x,s).
!
!     Names(nNames) is an character*8 array.
!             If nNames =  1, Names is not used.  The printed solution
!             will use generic names for the columns and row.
!             If nNames = n+m, Names(j) should contain an 8 character
!             name of the jth variable (j = 1, n+m).
!             If j = n+i, the jth variable is the ith row.
!
!     hs(n+m) sometimes contains a set of initial states for each
!             variable (x, s).  See the following NOTES.
!
!     x(n+m)  is a set of initial values for each variable (x, s).
!
!  NOTES:  1. If Start = 'Cold' or 'Basis file' and a BASIS file
!             of some sort is to be input
!             (an OLD BASIS file, INSERT file or LOAD file),
!             hs and x need not be set at all.
!
!          2. Otherwise, hs(1:n) must be defined for a cold start.
!             If nothing special is known about the problem, or if
!             there is no wish to provide special information,
!             you may set hs(j) = 0, x(j) = 0.0d+0 for all j=1:n.
!             All variables will be eligible for the initial basis.
!
!             Less trivially, to say that variable j will probably
!             be equal to one of its bounds,
!             set hs(j) = 4 and x(j) = bl(j)
!             or  hs(j) = 5 and x(j) = bu(j) as appropriate.
!
!          3. For Cold starts with no basis file, a Crash procedure
!             is used to select an initial basis.  The initial basis
!             matrix will be triangular (ignoring certain small
!             entries in each column).
!             The values hs(j) = 0, 1, 2, 3, 4, 5 have the following
!             meaning:
!
!             hs(j)    State of variable j during Crash
!
!             0, 1, 3  Eligible for the basis.  3 is given preference.
!             2, 4, 5  Ignored.
!
!             After Crash, hs(j) = 2 entries are made superbasic.
!             Other entries not selected for the basis are made
!             nonbasic at the value x(j) if bl(j) <= x(j) <= bu(j),
!             or at the value bl(j) or bu(j) closest to x(j).
!
!          4. For Warm starts, all of hs(1:n+m) is assumed to be
!             set to the values 0, 1, 2 or 3 from some previous call.
!
!     pi(m)   contains an estimate of the vector of Lagrange multipliers
!             (shadow prices) for the NONLINEAR constraints.  The first
!             nnCon components must be defined.  They will be used as
!             lambda in the subproblem objective function for the first
!             major iteration.  If nothing is known about lambda,
!             set pi(i) = 0.0d+0, i = 1 to nnCon.
!
!     nS      need not be specified for Cold starts,
!             but should retain its value from a previous call
!             when a Warm start is used.
!
!
!     ON EXIT:
!
!     hs(n+m) is the final state vector:
!
!                hs(j)    State of variable j    Normal value of x(j)
!
!                  0      nonbasic               bl(j)
!                  1      nonbasic               bu(j)
!                  2      superbasic             Between bl(j) and bu(j)
!                  3      basic                  ditto
!
!             Very occasionally there may be nonbasic variables for
!             which x(j) lies strictly between its bounds.
!             If nInf = 0, basic and superbasic variables may be outside
!             their bounds by as much as the Feasibility tolerance.
!             Note that if Scale is specified, the Feasibility tolerance
!             applies to the variables of the SCALED problem.
!             In this case, the variables of the original problem may be
!             as much as 0.1 outside their bounds, but this is unlikely
!             unless the problem is very badly scaled.
!
!     x(n+m)  contains the final variables and slacks (x, s).
!
!     pi(m)   is the vector of Lagrange multipliers (shadow prices)
!             for the general constraints.
!
!     rc(n+m) is a vector of reduced costs: rc = g - (J -I)'pi, where g
!             is the gradient of the objective if x is feasible
!             (or the gradient of the Phase-1 objective otherwise).
!             If nInf = 0, the last m entries are pi.
!
!     INFO    says what happened; see the User's Guide.
!             The possible values are as follows:
!
!             INFO    Meaning
!
!                0    finished successfully
!                1       optimality conditions satisfied
!                2       feasible point found
!                3       requested accuracy could not be achieved
!
!               10    the problem appears to be infeasible
!               11       infeasible linear constraints
!               12       infeasible linear equalities
!               13       nonlinear infeasibilities minimized
!               14       infeasibilities minimized
!
!               20    the problem appears to be unbounded
!               21       unbounded objective
!               22       constraint violation limit reached
!
!               30    resource limit error
!               31       iteration limit reached
!               32       major iteration limit reached
!               33       the superbasics limit is too small
!
!               40    terminated after numerical difficulties
!               41       current point cannot be improved
!               42       singular basis
!               43       cannot satisfy the general constraints
!               44       ill-conditioned null-space basis
!
!               50    error in the user-supplied functions
!               51       incorrect objective  derivatives
!               52       incorrect constraint derivatives
!
!               60    undefined user-supplied functions
!               61       undefined function at the first feasible point
!               62       undefined function at the initial point
!               63       unable to proceed into undefined region
!
!               70    user requested termination
!               71       terminated during function evaluation
!               72       terminated during constraint evaluation
!               73       terminated during objective evaluation
!               74       terminated from monitor routine
!
!               80    insufficient storage allocated
!               81       work arrays must have at least 500 elements
!               82       not enough character storage
!               83       not enough integer storage
!               84       not enough real storage
!
!               90    input arguments out of range
!               91       invalid input argument
!               92       basis file dimensions do not match this problem
!
!              140    system error
!              141       wrong no of basic variables
!              142       error in basis package
!
!     mincw   says how much character storage is needed to solve the
!             problem.  If INFO = 82, the work array cw(lencw) was
!             too small.  snOptA may be called again with lencw suitably
!             larger than mincw.
!
!     miniw   says how much integer storage is needed to solve the
!             problem.  If INFO = 83, the work array iw(leniw) was too
!             small.  snOptA may be called again with leniw suitably
!             larger than miniw.  (The bigger the better, since it is
!             not certain how much storage the basis factors need.)
!
!     minrw   says how much real storage is needed to solve the
!             problem.  If INFO = 84, the work array rw(lenrw) was too
!             small.  (See the comments above for miniw.)
!
!     nS      is the final number of superbasics.
!
!     nInf    is the number of infeasibilities.
!
!     sInf    is the sum    of infeasibilities.
!
!     Obj     is the value of the nonlinear part of the objective.
!             If nInf = 0, Obj includes the nonlinear objective if any.
!             If nInf > 0, Obj is just the linear objective if any.
!
!     cu(lencu), iu(leniu), ru(lenru)  are character, integer and real
!             arrays of USER workspace.  These arrays are available to
!             pass data to the user-defined routine userfun.
!             If no workspace is required, you can either use dummy
!             arrays for cu, iu and ru, or use cw, iw and rw
!             (see below).
!
!     cw(lencw), iw(leniw), rw(lenrw)  are character*8, integer and real
!             arrays of workspace used by snOptC.
!             lencw  should be about at least 500.
!             leniw  should be about max( 500, 20(m+n) ) or larger.
!             lenrw  should be about max( 500, 40(m+n) ) or larger.
!
!     snOptC is maintained by Philip E. Gill,
!     Dept of Mathematics, University of California, San Diego.
!
!     LUSOL is maintained by Michael A. Saunders,
!     Systems Optimization Laboratory,
!     Dept of Management Science & Engineering, Stanford University.
!
!     31 Oct 1998: First snoptC based on the snopt in SNOPT 5.3-4.
!     22 Dec 2002: Added input argument checking.
!     02 Jan 2003: Call sqOpt for LP's.
!     20 Jan 2003: QP-QN added.
!     31 Jul 2003: snEXIT and snPRNT adopted.
!     15 Oct 2004: snSTOP adopted.
!     01 Sep 2007: Sticky parameters added.
!     11 Sep 2014: nnObjU and nnObj separated for FP mode.
!     ==================================================================
      external
     &     snLog, snLog2, sqLog, snSTOP
!     ------------------------------------------------------------------
      call snKerC
     &   ( Start, m, n, neJ, nNames,
     &     nnCon, nnObjU, nnJac,
     &     iObjU, objUAdd, Prob,
     &     userfun, snLog, snLog2, sqLog, snSTOP,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snOptc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snKerC
     &   ( Start, m, n, neJ, nNames,
     &     nnCon, nnObjU, nnJac,
     &     iObjU, objUAdd, Prob,
     &     userfun, snLog, snLog2, sqLog, snSTOP,
     &     Jcol, indJ, locJ, bl, bu, Names,
     &     hs, x, pi, rc,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf, Obj,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     userfun, snLog, snLog2, sqLog, snSTOP
      integer
     &     INFO, iObjU, lencu, leniu, lenru, lencw, leniw, lenrw,
     &     mincw, miniw, minrw, m, n, neJ, nNames, nS, nInf, nnCon,
     &     nnObjU, nnJac, indJ(neJ), hs(n+m), locJ(n+1), iu(leniu),
     &     iw(leniw)
      double precision
     &     sInf, Obj, objUAdd, Jcol(neJ), bl(n+m), bu(n+m), x(n+m),
     &     pi(m), rc(n+m), ru(lenru), rw(lenrw)
      character*(*)
     &     Start
      character
     &     Prob*8, Names(nNames)*8, cu(lencu)*8, cw(lencw)*8

!     ==================================================================
!     snKerC does the work for snOptC. (Kernel for snoptC)
!
!     Developers can call this version with customized versions of
!     snLog, snLog2  and  snSTOP.
!
!     17 Oct 2004: First version of snKerC.
!     01 Sep 2007: Sticky parameters added.
!     04 Jul 2010: mincw, miniw, minrw added to workspace.
!     21 Aug 2013: Ad hoc allocation of x0 and gSave made safe.
!     11 Sep 2014: nnObjU and nnObj separated for FP mode.
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
     &     Errors, inform, iObj, HDInfo, lenR, lenx0, lgObj,
     &     leType, lkx, llocG, liwEst, lrwEst, lvlHes, lvlStart, lx0,
     &     maxcw, maxiw, maxR, maxrw, maxS, minmax, mProb, mQNmod, nb,
     &     negCon, nextcw, nextiw, nextrw, nInfE, nkx, ngQP, nnObj0,
     &     nnObj, nlocG, nlocJ, nMajor, nnH, nnH0, nrhs, nrhs0, nx0,
     &     startType, stickyOp
      integer
     &     neH, nlocH, indH(1), locH(1)
      double precision
     &     Hcol(1)
       double precision
     &     fObj, objAdd, objTrue, sInfE, rhs(1), x0(1)
      external
     &     s0fgC, s8HxLP, s8HxNP, s8HxQP
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
      Solver = 'SNOPTC'
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

!     Check the arguments of snOpt.

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
!     Print the options if iPrint > 0, Print level > 0 and lvlPrm > 0.
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
      lvlHes  = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian
      minmax  = iw( 87) ! 1, 0, -1  => MIN, FP, MAX

      lenR    = maxR*(maxR + 1)/2  +  (maxS - maxR)

      iw( 20) = negCon
      iw( 28) = lenR    ! R(lenR) is the reduced Hessian factor

!     Load the iw array with various problem dimensions.

      nnH     = max( nnJac, nnObjU )

      neH     = 1       ! Placeholders
      nlocH   = 1

      iw( 15) = n       ! copy of the number of columns
      iw( 16) = m       ! copy of the number of rows
      iw( 17) = neJ     ! copy of the number of nonzeros in Jcol
      iw( 21) = nnJac   ! # of Jacobian variables
      iw( 22) = nnObjU  ! # of user-defined objective variables
      iw( 23) = nnCon   ! # of nonlinear constraints
      iw( 24) = nnH     !   max( nnObjU, nnJac )
      iw(204) = iObjU   ! position of the objective row in J

!     ------------------------------------------------------------------
!     If only a feasible point is requested, save the base point for the
!     objective function:  1/2 || x - x0 ||^2
!
!     Set the working objective gradient dimensions.
!     ------------------------------------------------------------------
      FPonly  = minmax .eq. 0
      if ( FPonly ) then
         nnObj  = nnH
         iObj   = 0
         objAdd = zero
      else
         nnObj  = nnObjU       ! working # nonlinear objective vars
         iObj   = iObjU
         objAdd = objUAdd
      end if

      nnObj0 = max( nnObj, 1   )

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

*     ------------------------------------------------------------------
*     If only a feasible point is requested, save the base point for the
*     objective function:  1/2 || x - x0 ||^2
*     ------------------------------------------------------------------
      if (FPonly) then
         lx0 = iw(298)
         call dcopy ( nnH, x, 1, rw(lx0), 1 )
      end if

!     Define the row and column ordering for J.
!     snOptC  uses natural order throughout, so kx = kxN.

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

         nrhs   = 0             ! No constraint rhs vector.
         nx0    = 0             ! No constant shift for x.
         lenx0  = max( nx0  , 1   )
         nrhs0  = 1
         nnH0   = 1
         ngQP   = nnObj

         leType = iw(283)      ! eType(nb) definition of elastic vars
         lgObj  = iw(297)      ! gObj(nnObjU) = Objective gradient

         call iload ( nb, 3, iw(leType), 1 )

         call s5solve
     &      ( INFO,
     &        Solver, startType,
     &        sqLog, s8HxQP, s8HxLP, GotR,
     &        m, n, nb, nnH0, nnH, nNames,
     &        ngQP, nnObj0, nnObj,
     &        iObj, objAdd, fObj, objTrue,
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
     &        s0fgC, userfun, userfun, s8HxNP,
     &        snLog, snLog2, snSTOP, gotR,
     &        m, n, nb, nnH, nnCon, nnJac, nnObj,
     &        nNames, iObj, objAdd, fObj, objTrue,
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
      Obj    = fObj

      mincw  = iw(47)           ! minimum length of cw
      miniw  = iw(48)           ! minimum length of iw
      minrw  = iw(49)           ! minimum length of rw

!     If "sticky parameters no",  restore the user-defined options

      stickyOp = iw(116)

      if (stickyOp .le. 0) then
         call chcopy
     &      ( 130, Usercw, 1, cw(51), 1 )
         call icopy
     &      ( 130, Useriw, 1, iw(51), 1 )
         call dcopy
     &      ( 130, Userrw, 1, rw(51), 1 )
      end if

!     Print times for all clocks (if lvlTim > 0).

      call s1time( 0, 2, iw, leniw, rw, lenrw )

      return

!     Local exit messages.

  800 call snWRAP( INFO, Solver, str, str2, iw, leniw )

  999 return

      end ! subroutine snKerC
