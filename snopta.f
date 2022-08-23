!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     file  snopta.f --- the free format interface for snOpt
!
!     snOptA    snKerA
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snOptA
     &   ( start, nF, n, nxname, nFname,
     &     objUAdd, objRow, Prob, usrfun,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     usrfun
      integer
     &     INFO, lenA, lencu, lencw, lenG, leniu, leniw, lenru, lenrw,
     &     mincw, miniw, minrw, n, neA, neG, nF, nFname, nInf, nS,
     &     nxname, objRow, start, iAfun(lenA), iGfun(lenG), iu(leniu),
     &     iw(leniw), jAvar(lenA), jGvar(lenG), xstate(n), Fstate(nF)
      double precision
     &     objUAdd, sInf, A(lenA), F(nF), Fmul(nF), Flow(nF), Fupp(nF),
     &     ru(lenru), rw(lenrw), x(n), xlow(n), xmul(n), xupp(n)
      character
     &     Prob*8, cu(lencu)*8, cw(lencw)*8, Fnames(nFname)*8,
     &     xnames(nxname)*8

!     ==================================================================
!     snOptA  is a Fortran wrappper for the SNOPT solver.
!     snOptA   is a subroutine for constrained nonlinear
!     optimization.  The optimization problem involves m  functions
!     F(1), F(2), ... , F(nF), each a function of n variables
!     x(1), x(2), ... , x(n).  The  problem has the form:
!
!           minimize/maximize    objUAdd + F(objRow)
!
!                            ( xlow <=  x  <= xupp,
!                 subject to (
!                            ( Flow <=  F  <= Fupp,
!
!     where objUAdd is a constant, objRow is a user-specified row of  F,
!     xlow, Flow, xupp and Fupp are constant lower and upper bounds.
!
!     ------------------------------------------------------------------
!     NOTE: Before calling SNOPTA, your calling program MUST call the
!     initialization routine using the call:
!     call snInit( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )
!     This sets the default values of the optional parameters. You can
!     also alter the default values of iPrint and iSumm before snOptA
!     is used.  iPrint = 0, etc, is OK.
!     ------------------------------------------------------------------
!
!     o If objRow = 0, then snOptA will find a point satisfying the
!       constraints.
!
!     o If all functions are linear, F = A x for some sparse matrix A.
!       This defines a linear program (LP).  In this case,  the nonzero
!       elements of A can be input in coordinate form (i,j,A_ij) (see
!       below).
!
!     o If all functions are nonlinear, F = F(x) for some vector
!       F(x) of smooth functions.  In this case, the elements of  F  and
!       (optionally) their first and second partial derivatives must be
!       coded by the user in the subroutine usrfun  (see below).
!
!     o If some functions are linear and some are nonlinear, the user
!       can choose to set every component in usrfun.  It is usually more
!       efficient, however,  to supply the coefficients of the linear
!       functions via the sparse array  A (see below).   In this case,
!       the linear elements of  F  need not be assigned (SNOPTA will
!       figure out which components of  F  are needed).
!
!     o In the most general situation, the ith component of F(x) is the
!       sum of linear and nonlinear terms.  In this case, if F(x) can be
!       defined as a sum of "non-overlapping" linear and nonlinear
!       functions, then the nonlinear part of F can be defined in usrfun
!       and the linear part can be defined via the array A.
!
!       Suppose that the ith component of F(x) is of the form
!            F_i(x) = f_i(x) + sum (over j)  A_ij x_j,
!       where f_i(x) is a nonlinear function and the elements A_ij
!       are constant.   It is convenient to write  F_i(x)  in the
!       compact form  F_i(x) = f_i(x) + A_i' x, where A_i denotes a
!       column vector with components ( A_i1, A_i2, ..., A_in ), and
!       "A_i'" denotes the transpose of A_i.
!
!       Functions f_i and A_i are said to be "non-overlapping" if any
!       variable x_j  appearing explicitly in f_i(x) does not appear
!       explicitly in A_i'x, i.e., A_ij = 0.  (Equivalently, any
!       variable with a nonzero A_ij must not appear explicitly in
!       f_i(x).)  For example, the function
!         F_i(x) = 3x_1 + exp(x_2)x_4 + x_2^2 + 4x_4 - x_3 + x_5
!       can be written as the sum of non-overlapping functions f_i and
!       A_i'x, such that
!           f_i(x) = exp(x_2)x_4 + x_2^2  + 4x_4  and
!           A_i'x  = 3x_1 - x_3 + x_5.
!
!       Given a non-overlapping sum for each component of F, we can
!       write  F(x) = f(x) + Ax, where f(x) is a vector-valued function
!       of x and A is a sparse matrix whose ith row is A_i'.
!
!       The nF by n  Jacobian of  F(x)  is the sum of two  nF by n
!       sparse matrices G and A,  i.e.,  J = G + A,  where G and A
!       contain the nonlinear and constant elements of J respectively.
!       The important property of non-overlapping functions is that
!       a nonzero entry of J is either an element of A, or an element
!       of G, but NOT BOTH (i.e., the nonzeros of  A  and  G  do not
!       overlap.
!
!       The nonzero elements of A and G must be provided in coordinate
!       form.  In coordinate form, a nonzero element G_ij of a matrix
!       G  is stored as the triple (i,j,G_ij).  The kth coordinate is
!       defined by iGfun(k) and jGvar(k)  (i.e., if i=iGfun(k) and
!       j=jGvar(k), then G(k) is the ijth element of G.)  Any known
!       values of G(k) must be assigned by the user in the routine
!       usrfun.
!
!       RESTRICTIONS:
!        1.  If the elements of G cannot be provided because they are
!            either too expensive or too complicated to evaluate,  it
!            is still necessary to specify the position of the nonzeros
!            as specified by the arrays iGfun and jGvar.
!
!        2.  If an element of G happens to be zero at a given point,
!            it must still be loaded in usrfun. (The order of the
!            list of coordinates (triples) is meaningful in snOptA.)
!
!       The elements of A and G can be stored in any order, (e.g., by
!       rows, by columns, or mixed). Duplicate entries are ignored.
!
!     ON ENTRY:
!
!     start   specifies how a starting basis (and certain other items)
!             are to be obtained.
!             start =  0 (Cold) means that Crash should be used to
!                      choose an initial basis, unless a basis file is
!                      given by reference in the Specs file to an
!                      Old basis file.
!             start =  1 (Basis file) means the same (but is more
!                      meaningful in the latter case).
!             start =  2 (Warm) means that a basis is already defined
!                      in xstate and Fstate (probably from an earlier
!                      call).
!
!     nF      is the number  of problem functions in F, including the
!             objective function (if any) and the linear
!             and nonlinear constraints.  Simple upper and lower bound
!             constraints on the variables should not be included in  F.
!             nF > 0.
!
!     n       is the number of variables.
!             n > 0.
!
!     neA     is the number of nonzero entries in A.
!             neA >= 0.
!
!     nxname  is the number of 8-character column (i.e., variable) names
!             provided in the array xnames.  If nxname = 1,  then there
!             are NO column names (generic names will be used in the
!             printed solution).  Otherwise, nxname = n and every
!             column name must be provided.
!
!     nFname  is the number of 8-character row (i.e., constraint and
!             objective) names provided in the array Fnames.
!             If nFname = 1,  then there are NO row names (generic
!             names will be used in the printed solution).  Otherwise,
!             nFname = nF and every row name must be provided.
!
!     objUAdd is a constant that will be added to the objective.
!             Typically objUAdd = 0.0d+0.
!
!     Prob    is an 8-character name for the problem, used in the
!             output.  A blank name can be assigned if necessary.
!
!     xlow(n) are the lower bounds on x.
!
!     xupp(n) are the upper bounds on x.
!
!     xnames(nxname) is an character*8 array of names for each x(j).
!             If nxname =  1, xnames is not used.  The printed solution
!             will use generic names for the variables.
!             If nxname = n, xnames(j) should contain an 8 character
!             name of the jth variable (j = 1:n).
!
!     Flow(n) are the lower bounds on F.  If component F(objRow)
!             is being optimized,  Flow(objRow) is ignored.
!
!     Fupp(n) are the upper bounds on F.  If component F(objRow)
!             is being optimized,  Fupp(objRow) is ignored.
!
!     Fnames(nFname) is an character*8 array of names for each F(i).
!             If nFname =  1, Fnames is not used.  The printed solution
!             will use generic names for the objective and constraints.
!             If nNames = nF, Fnames(j) should contain an 8 character
!             name of the jth constraint (j=1:nF).
!
!     xstate(n) sometimes contains a set of initial states for each
!             variable x.  See the following NOTES.
!
!     x(n)    is a set of initial values for each variable  x.
!
!  NOTES:  1. If start = 0 (Cold) or 1 (Basis file) and an OLD BASIS
!             file is to be input, xstate and x need not be set at all.
!
!          2. Otherwise, xstate(1:n) must be defined for a cold start.
!             If nothing special is known about the problem, or if
!             there is no wish to provide special information,
!             you may set xstate(j) = 0, x(j) = 0.0d+0 for all j=1:n.
!             All variables will be eligible for the initial basis.
!
!             Less trivially, to say that variable j will probably
!             be equal to one of its bounds,
!             set xstate(j) = 4 and x(j) = bl(j)
!             or  xstate(j) = 5 and x(j) = bu(j) as appropriate.
!
!          3. For Cold starts with no basis file, a Crash procedure
!             is used to select an initial basis.  The initial basis
!             matrix will be triangular (ignoring certain small
!             entries in each column).
!             The values xstate(j) = 0, 1, 2, 3, 4, 5 have the following
!             meaning:
!
!             xstate(j)  State of variable j during Crash
!             ---------  --------------------------------
!             0, 1, 3    Eligible for the basis.  3 is given preference.
!             2, 4, 5    Ignored.
!
!             After Crash, xstate(j) = 2 entries are made superbasic.
!             Other entries not selected for the basis are made
!             nonbasic at the value x(j) if bl(j) <= x(j) <= bu(j),
!             or at the value bl(j) or bu(j) closest to x(j).
!
!          4. For Warm starts, all of Fstate(1:nF) is assumed to be
!             set to the values 0, 1, 2 or 3 from some previous call.
!
!     Fmul(nF) contains an estimate of the Lagrange multipliers
!             (shadow prices) for the F- constraints.  They are used
!             to define the Lagrangian for the first major iteration.
!             If nothing is known about Fmul, set
!             Fmul(i) = 0.0d+0, i = 1:nF
!
!     ON EXIT:
!
!     xstate(n) is the final state vector for x:
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
!     x(n)    contains the final variables.
!
!     F(nF)   contains the final values of F.
!
!     xmul(nF) is the vector of Lagrange multipliers (shadow prices)
!             for the variables constraints.
!
!     Fmul(nF) is the vector of Lagrange multipliers (shadow prices)
!             for the general constraints.
!
!     INFO    says what happened; see the User's Guide.
!             The possible values are as follows:
!
!             INFO       Meaning
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
!               93       the QP Hessian is indefinite
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
!     cu(lencu), iu(leniu), ru(lenru)  are character, integer and real
!             arrays of USER workspace.  These arrays are available to
!             pass data to the user-defined routine usrfun.
!             If no workspace is required, you can either use dummy
!             arrays for cu, iu and ru, or use cw, iw and rw
!             (see below).
!
!     cw(lencw), iw(leniw), rw(lenrw)  are character*8, integer and real
!             arrays of workspace used by snOptA.
!             lencw should be at least 500, or nF+n if names are given.
!                              +.
!             leniw should be about max( 500, 20(nF+n) ) or larger.
!             lenrw should be about max( 500, 40(nF+n) ) or larger.
!
!     SNOPT package maintained by Philip E. Gill,
!     Dept of Mathematics, University of California, San Diego.
!
!     08 Nov 1998: First version based on the snopt of SNOPT 5.3-4.
!     25 Aug 1999: for SNOPT Version 6.0.
!     04 Nov 2001: LP's solved explicitly.
!     31 Jul 2003: snEXIT and snPRNT adopted.
!     02 May 2004: Call to base routine added.
!     01 Sep 2007: Sticky parameters added.
!     11 Sep 2014: Added nnObjU and nnObj for FP mode.
!     08 Mar 2015: objRow may be a parameter. Passed local value.
!     ==================================================================
      integer
     &     objRowU
      external
     &     snLog, snLog2, sqLog, snSTOP
!     ------------------------------------------------------------------
      objRowU = objRow
      call snKerA
     &   ( start, nF, n, nxname, nFname,
     &     objUAdd, objRow, Prob,
     &     usrfun, snLog, snLog2, sqLog, snSTOP,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      end ! subroutine snOptA

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snKerA
     &   ( start, nF, n, nxname, nFname,
     &     objUAdd, objRow, Prob,
     &     usrfun, snLog, snLog2, sqLog, snSTOP,
     &     iAfun, jAvar, lenA, neA, A,
     &     iGfun, jGvar, lenG, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     x, xstate, xmul, F, Fstate, Fmul,
     &     INFO, mincw, miniw, minrw,
     &     nS, nInf, sInf,
     &     cu, lencu, iu, leniu, ru, lenru,
     &     cw, lencw, iw, leniw, rw, lenrw )

      implicit
     &     none
      external
     &     usrfun, snLog, snLog2, sqLog, snSTOP
      integer
     &     INFO, lenA, lencu, lencw, lenG, leniu, leniw, lenru, lenrw,
     &     mincw, miniw, minrw, n, neA, neG, nF, nFname, nInf, nS,
     &     nxname, objRow, start, iAfun(lenA), iGfun(lenG), iu(leniu),
     &     iw(leniw), jAvar(lenA), jGvar(lenG), xstate(n), Fstate(nF)
      double precision
     &     objUAdd, sInf, A(lenA), F(nF), Fmul(nF), Flow(nF), Fupp(nF),
     &     ru(lenru), rw(lenrw), x(n), xlow(n), xmul(n), xupp(n)
      character
     &     Prob*8, cu(lencu)*8, cw(lencw)*8, Fnames(nFname)*8,
     &     xnames(nxname)*8

!     ==================================================================
!     snKerA does the work for snOptA. (Kernel for snoptA)
!
!     Developers can call this version with customized versions of
!     snLog, snLog2  and  snSTOP.
!
!     17 Oct 2004: First version of snKerA
!     01 Sep 2007: Sticky parameters added
!     24 Nov 2007: s3defaultsA sets options dependent on prob dimensions
!     04 Jul 2010: mincw, miniw, minrw added to workspace
!     11 Sep 2014: Added nnObjU and nnObj for FP mode
!     ==================================================================
      character
     &     Solver*6, str*80, str2*80, Usercw(130)*8
      integer
     &     Useriw(130)
      double precision
     &     Userrw(130)
      logical
     &     gotR, PrintMem
      integer
     &     Errors, HDInfo, inform, iObj, iObjU, lbl, lbu, lenR,
     &     lenx0, leType, lFx, lgObj, lgObj1, lgObj2, liGfun, ljGvar,
     &     lJcol, lkx, lkxN, lhs, liwEst, lrwEst, liy, llocJ, lindJ,
     &     llocG, lindG, lNames, lnGlin, lpi, lrc, lvlHes, lvlStart,
     &     lx0, lx, maxcw, maxiw, maxrw, mProb, m, maxR, maxS,
     &     minBld, mQNmod, nb, neJ, negCon, nextcw, nextiw, nextrw,
     &     ngQP, nInfE, nkx, nNames, nnCon, nnH0, nnH, nnJac, nnObj0,
     &     nnObj, nnObjU, nlocJ, nlocG, nMajor, nrhs0, npStatus, nrhs,
     &     nx0, objSave, ObjSpc, qpStatus, startType, stickyOp
      integer
     &      lenRTmp, lvlHesTmp, maxRTmp, maxSTmp, mQNmodTmp
      integer
     &     neH, nlocH, indH(1), locH(1)
      double precision
     &     Hcol(1)
      double precision
     &     fObj, objAdd, objTrue, infBnd, sInfE, rhs(1), x0(1)
      external
     &     s0fgA, s8HxLP, s8HxNP, s8HxQP
!     ------------------------------------------------------------------
      integer            COLD,       BASIS,      WARM,       HOT
      parameter         (COLD   = 0, BASIS  = 1, WARM   = 2, HOT    = 3)
      integer            StdIn
      parameter         (StdIn     = 2)
      integer            Unset,         Normal
      parameter         (Unset     =-1, Normal = 0)

      parameter         (mProb     =  51) ! Problem name
      parameter         (lvlStart  =  69) ! cold:warm:basis:hot start
      parameter         (HDInfo    = 243) ! Approximate Hessian type

      parameter         (qpStatus  = 235) ! QP user-routine call-status
      parameter         (npStatus  = 236) ! NP user-routine call-status

      double precision   zero
      parameter         (zero      = 0.0d+0)
      integer            idummy
      parameter         (idummy    =  -11111)
!     ------------------------------------------------------------------
      Solver = 'SNOPTA'
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

!     Check the arguments of snOptA.

      startType = iw(lvlStart)
      call s3chkArgsA
     &   ( inform, start, nF, n, nS, nxname, nFname,
     &     objRow, neA, neG,
     &     xlow, xupp, xnames, Flow, Fupp, Fnames,
     &     xstate, xmul, Fstate, Fmul, startType, Errors,
     &     iw, leniw, rw, lenrw )
      if (inform .gt. 0) then
         INFO = inform
         go to 800
      end if

!     ------------------------------------------------------------------
!     Any objective row specified in the specs file overrides  objRow.
!     ------------------------------------------------------------------
!     There is always an objective function, even if the user didn't
!     specify one.

      objSave   = objRow
      ObjSpc    = iw(103) ! Optional parameter: objective row.
      if (ObjSpc .ne. idummy) then
         objRow  = ObjSpc
      else
         iw(103) = objRow
      end if

!     Allocate temporary work arrays for s3sizeA.

      nkx    = n + nF
      nlocJ  = n + 1

!     Permanent addresses first.

      lkx    = nextiw
      minbld = lkx + nkx

      if (minbld .gt. maxiw) then
!        ---------------------------------------------------------------
!        Not enough memory to build the problem.
!        Provide the user an (over) estimate of what is needed. The
!        problem dimensions have not been computed yet, so s3defaultsA
!        assigns temporary estimates of  lvlHes, maxR, maxS  and  mQNmod
!        if they have not yet been specified by the user.
!        ---------------------------------------------------------------
         neJ   = neA  + neG
         m     = nF
         nlocG = n    + 1

         if (nxname .eq. 1  .and.  nFname .eq. 1) then
            nNames = 1
         else
            nNames = n + m
         end if

         nnCon  = m
         nnJac  = n
         nnObjU = n
         nnObj  = n
         negCon = neJ
         nnH    = max( nnJac, nnObjU )

         call s3defaultsA
     &      ( n, lvlHesTmp, maxRTmp, maxSTmp, mQNmodTmp, iw, leniw )

         lenRTmp = maxRTmp*(maxRTmp+1)/2 + (maxSTmp - maxRTmp)

         call s8Map
     &      ( m, n, negCon, nkx, nnCon, nnJac, nnObjU, nnObj, nnH,
     &        lenRTmp, maxRTmp, maxSTmp, mQNmodTmp, lvlHesTmp,
     &        nextcw, nextiw, nextrw, iw, leniw )
         call s3mapA
     &      ( m, n, neJ, nF, neG, negCon, nkx, nnJac, nNames,
     &        nextcw, nextiw, nextrw, iw, leniw )
         call s2Bmap
     &      ( m, n, neJ, maxSTmp,
     &        nextiw, nextrw, maxiw, maxrw, liwEst, lrwEst, iw, leniw )
         PrintMem = .true.        ! Print all messages in s2Mem
         call s2Mem
     &      ( inform, PrintMem, liwEst, lrwEst,
     &        nextcw, nextiw, nextrw,
     &        maxcw, maxiw, maxrw, lencw, leniw, lenrw,
     &        mincw, miniw, minrw, iw )
         INFO = inform
         go to 800
      end if

!     Compute  m, negCon, neJ, nnCon, nnJac, nnObjU and iObjU.
!     The integer array kx defines the order of the variables
!     and constraints given to SNOPT.

      call s3sizeA
     &   ( INFO, n, nF, nkx, objRow,
     &     iAfun, jAvar, lenA, neA, iGfun, jGvar, lenG, neG,
     &     m, negCon, neJ, nnCon, nnJac, nnObjU, iObjU,
     &     iw(lkx), leniw, iw )
      if (INFO .gt. 0) then
         go to 800
      end if

!     ------------------------------------------------------------------
!     Values of  neJ,  nnCon,  nnJac  and  nnObj  are now known exactly.
!     Load the iw array with various problem dimensions.
!     ------------------------------------------------------------------
      nnH     = max( nnJac, nnObjU )
      neH     = 1                ! Placeholders
      nlocH   = 1


      iw( 15) = n      ! copy of the number of columns
      iw( 16) = m      ! copy of the number of rows
      iw( 17) = neJ    ! copy of the number of nonzeros in Jcol
      iw( 20) = negCon ! # of nonzero elems in J
      iw( 21) = nnJac  ! # of Jacobian  variables
      iw( 22) = nnObjU ! # of objective variables
      iw( 23) = nnCon  ! # of nonlinear constraints
      iw( 24) = nnH    !   max( nnObjU, nnJac )
      iw(204) = iObjU  ! position of the objective row in J

      iw(248) = nF     ! # of components of user-defined F
      iw(249) = neG    ! # of components of user-defined G

!     ------------------------------------------------------------------
!     The obligatory call to snInit has already set the defaults.
!     All problem dimensions have been computed.
!     Check that the optional parameters have sensible values.
!     Print the options.
!     ------------------------------------------------------------------
      cw(mProb)  = Prob

      call s8Defaults
     &   ( m, n, nnCon, nnJac, nnObjU, iObjU,
     &     cw, lencw, iw, leniw, rw, lenrw )
      call s3printA
     &   ( m, n, nnCon, nnJac, nnObjU, startType, iw, leniw, rw, lenrw )

!     ------------------------------------------------------------------
!     Compute the addresses of all work arrays.
!     ------------------------------------------------------------------
      nb     = n     + m
      nlocG  = nnJac + 1

      if (nxname .eq. 1  .and.  nFname .eq. 1) then
         nNames = 1
      else
         nNames = nb
      end if

      maxR    = iw( 52) ! max columns of R.
      maxS    = iw( 53) ! max # of superbasics
      mQNmod  = iw( 54) ! (ge 0) max # of BFGS updates
      lvlHes  = iw( 72) ! 0,1,2  => LM, FM, Exact Hessian

      lenR    = maxR*(maxR + 1)/2  +  (maxS - maxR)
      iw( 28) = lenR

!     ------------------------------------------------------------------
!     If only a feasible point is requested, save the base point for the
!     objective function:  1/2 || x - x0 ||^2
!
!     Set the working objective gradient dimensions.
!     ------------------------------------------------------------------
      if (objRow .eq. 0) then
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
!     Allocate the local arrays for snOptA.
!     s8Map  maps snOptA integer and double arrays.
!     s3mapA maps additional arrays for snOptA.
!     s2BMap maps the arrays for the LU routines.
!     s2Mem  checks what space is available and prints any messages.
!     ------------------------------------------------------------------
      call s8Map
     &   ( m, n, negCon, nkx, nnCon, nnJac, nnObjU, nnObj, nnH,
     &     lenR, maxR, maxS, mQNmod, lvlHes,
     &     nextcw, nextiw, nextrw, iw, leniw )
      call s3mapA
     &   ( m, n, neJ, nF, neG, negCon, nkx, nnJac, nNames,
     &     nextcw, nextiw, nextrw, iw, leniw )
      call s2Bmap
     &   ( m, n, neJ, maxS,
     &     nextiw, nextrw, maxiw, maxrw, liwEst, lrwEst, iw, leniw )
      PrintMem = .true.         ! OK to print messages in s2Mem
      call s2Mem
     &   ( inform, PrintMem, liwEst, lrwEst,
     &     nextcw, nextiw, nextrw,
     &     maxcw, maxiw, maxrw, lencw, leniw, lenrw,
     &     mincw, miniw, minrw, iw )
      if (inform .ne. 0) then
         INFO = inform
         go to 800
      end if

!     Allocate local work arrays.

      lkxN   = iw(252) ! jN = kxN(j ) => col j of Jcol is variable jN
      lJcol  = iw(256) ! Jcol(ne)    = Constraint Jacobian by columns
      llocJ  = iw(257) ! locJ(n+1)   = column pointers for indJ
      lindJ  = iw(258) ! indJ(ne) holds the row indices for Jij

      llocG  = iw(260) ! locG(nlocG) = column pointers for indG
      lindG  = iw(261) ! indG(neG) holds the row indices for gij
      lnGlin = iw(262) ! nGlin(j) = # linear elems in col j of gCon

      liGfun = iw(266) ! iGfun(neG) row list of reordered G nonzeros
      ljGvar = iw(267) ! iGvar(neG) col list of reordered G nonzeros

      lgObj  = iw(297) ! gObj(nnObj) = Objective gradient
      lgObj1 = iw(324) ! gObj1(nnObj) objective gradients at x1
      lgObj2 = iw(325) ! gObj2(nnObj) work gObj

      lbl    = iw(269) ! bl(nb) = user-defined lower bounds for SNOPTA
      lbu    = iw(270) ! bu(nb) = user-defined upper-bounds for SNOPTA
      lpi    = iw(279) ! pi(m)  = the pi-vector
      lrc    = iw(280) ! rc(nb) = the reduced costs
      lhs    = iw(282) ! the column state vector
      leType = iw(283) ! eType(nb) definition of elastic vars
      lx     = iw(299) ! x(nb)       = the solution (x,s)

      liy    = iw(308) ! iy (nb)     =  integer work vector
      lFx    = iw(336) ! Fx (nnCon)  = F(x) + A(linear)x

      lNames = iw(359) ! Names(nNames), row and column names

!     ------------------------------------------------------------------
!     Build the column-wise data structure for the Jacobian.
!     ------------------------------------------------------------------
      call s3buildA
     &   ( objRow, n, nkx, nnCon, nnJac, iw(lkx), iw(lnGlin),
     &     iAfun, jAvar, lenA, neA, A, iGfun, jGvar, lenG, neG,
     &     neJ   , nlocJ, iw(llocJ), iw(lindJ), rw(lJcol),
     &     negCon, nlocG, iw(llocG), iw(lindG), iw(liy) )

!     ------------------------------------------------------------------
!     Re-order the input data and invert the row and column orderings.
!     ------------------------------------------------------------------
      call s3permuteA
     &   ( n, nF, nkx,
     &     iGfun, jGvar, iw(liGfun), iw(ljGvar), lenG, neG,
     &     iw(lkx), iw(lkxN) )

      iw(247) = nkx     ! dimension of kx and its inverse, kxN

!     ------------------------------------------------------------------
!     Load the arrays used by SNOPTA.
!     These are for the data,
!              Jcol, indJ, locJ, bl, bu
!     and for the solution
!              hs, x, pi, rc, hs.
!     ------------------------------------------------------------------
      infBnd  = rw( 70) ! definition of an infinite bound.

      call s3inA
     &   ( startType, iObjU,
     &     m, n, nb, nnCon, nF, nkx, iw(lkxN), infBnd,
     &     xnames, Fnames, cw(lNames), nxname, nFname, nNames,
     &     xlow, xupp, Flow, Fupp, rw(lbl), rw(lbu), xstate, Fstate,
     &     iw(lhs), x, F, rw(lx), rw(lFx), Fmul, rw(lpi) )

*     ------------------------------------------------------------------
*     If only a feasible point is requested, save the base point for the
*     objective function:  1/2 || x - x0 ||^2
*     ------------------------------------------------------------------
      if (objRow .eq. 0) then
         lx0 = iw(298)
         call dcopy ( nb, rw(lx), 1,  rw(lx0), 1 )
      end if

!     ------------------------------------------------------------------
!     Sparse obj. gradients are scattered into gObj, gObj1 and gObj2.
!     ------------------------------------------------------------------
      call dload ( nnObj, zero, rw(lgObj) , 1 )
      call dload ( nnObj, zero, rw(lgObj1), 1 )
      call dload ( nnObj, zero, rw(lgObj2), 1 )

!     ------------------------------------------------------------------
!     The problem has been built.  The call-status variables npstat and
!     qpstat are reset to ensure that the necessary housekeeping is done
!     on the first call to s0fgA.
!
!     In future versions, the function build will be separated from the
!     snoptA. This will allow the function build to be called from the
!     function wrapper s0fgA.
!     ------------------------------------------------------------------
      iw(qpStatus) = idummy
      iw(npStatus) = idummy

!     ------------------------------------------------------------------
!     Solve the problem.
!     ------------------------------------------------------------------
      if (nnH .eq. 0) then

!        The problem is a linear program.

         nrhs   = 0             ! No constraint rhs vector.
         nx0    = 0             ! No constant shift for x.
         lenx0  = 1
         nrhs0  = 1
         nnH0   = 1
         ngQP   = nnObj

         call iload
     &      ( nb, 3, iw(leType), 1 )

         call s5solve
     &      ( INFO,
     &        Solver, startType,
     &        sqLog, s8HxQP, s8HxLP, gotR,
     &        m, n, nb, nnH0, nnH,
     &        nNames, ngQP, nnObj0, nnObj,
     &        iObj, objadd, fObj, objTrue,
     &        nInf, sInf, nInfE, sInfE,
     &        neJ, nlocJ, iw(llocJ), iw(lindJ), rw(lJcol),
     &        neH, nlocH,     locH,      indH,      Hcol,
     &        rw(lbl), rw(lbu), rw(lgObj), cw(lNames),
     &        nrhs0, nrhs, rhs, lenx0, nx0, x0,
     &        iw(leType), iw(lhs), rw(lx), rw(lpi), rw(lrc), nS,
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
     &      ( INFO,
     &        Solver, startType,
     &        s0fgA, usrfun, usrfun, s8HxNP,
     &        snLog, snLog2, snSTOP, gotR,
     &        m, n, nb, nnH, nnCon, nnJac, nnObj,
     &        nNames, iObj, objadd, fObj, objTrue,
     &        nInf, sInf, nInfE, sInfE,
     &        neJ, nlocJ, iw(llocJ), iw(lindJ), rw(lJcol),
     &        neH, nlocH,     locH,      indH,      Hcol,
     &        rw(lbl), rw(lbu), cw(lNames),
     &        iw(lhs), rw(lx), rw(lpi), rw(lrc), nMajor, nS,
     &        cu, lencu, iu, leniu, ru, lenru,
     &        cw, lencw, iw, leniw, rw, lenrw )
      end if

!     ------------------------------------------------------------------
!     Restore and update the user data.
!     ------------------------------------------------------------------
      call s3outA
     &   ( n, nb, nF, nnCon, nkx, iw(lkxN), objRow, iObjU,
     &     fObj, xstate, Fstate, iw(lhs), x, F, rw(lx),
     &     rw(lFx), xmul, Fmul, rw(lrc) )

      nInf   = nInf + nInfE
      sInf   = sInf + sInfE
      objRow = objSave

      mincw  = iw(47)            ! minimum length of cw
      miniw  = iw(48)            ! minimum length of iw
      minrw  = iw(49)            ! minimum length of rw

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

      end ! subroutine snKerA
