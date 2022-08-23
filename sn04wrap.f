!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn04wrap.f
!
!     snEXIT   snWRAP   snSolF
!
!     13 Dec 2013: Three subroutines separated from sn02lib to build
!                  standalone SQOPT library.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snEXIT
     &   ( iExit, Solver, string, string2 )

      implicit
     &     none
      character*(*)
     &     string, string2
      character
     &     Solver*6
      integer
     &     iExit

!     ==================================================================
!     snEXIT  returns the strings associated with EXIT condition iExit.
!
!     On exit, string1 and string2 are trimmed of trailing blanks and
!              each have a maximum length of 76 chars.
!
!     25 Sep 2002: First version of snEXIT.
!     12 Mar 2004: Trimmed message strings to avoid trouble on SGI's.
!     07 May 2006: Second derivative exits added (Exit 54)
!     01 Jun 2006: Stand-alone derivative checker added (Exits 55,105)
!     01 Jun 2006: Student edition exit added (Exit 66)
!     22 Apr 2007: Unrecognized options flagged (Exit 135)
!     22 Apr 2007: No derivatives checked (Exit 106)
!     16 Mar 2014: Added time limit (Exit 34)
!     17 Sep 2014: Added explicit elastic mode (Exits 5, 6)
!     02 Nov 2014: Added no acceptable LU (Exit 45)
!     ==================================================================
      integer
     &     i1, i2, mjr, mnr
!     ------------------------------------------------------------------
      integer
     &     length
      integer            Msgs
      parameter         (Msgs = 74)
      integer
     &     indc(0:14)
      character
     &     c(0:Msgs)*80
      data
     &     indc
     &     / 0, 7, 14, 17, 22, 28, 34, 38, 43, 48, 52, 60, 64, 66, 71 /
      data
     & c( 0)/'finished successfully'/,                             ! EXIT  0
     & c( 1)   /'optimality conditions satisfied'/,                ! EXIT  1
     & c( 2)   /'feasible point found'/,                           ! EXIT  2
     & c( 3)   /'requested accuracy could not be achieved'/,       ! EXIT  3
     & c( 4)   /'weak QP minimizer'/,                              ! EXIT  4
     & c( 5)   /'elastic objective minimized'/,                    ! EXIT  5
     & c( 6)   /'elastic infeasibilities minimized'/,              ! EXIT  6
     & c( 7)/'the problem appears to be infeasible'/,              ! EXIT 10
     & c( 8)   /'infeasible linear constraints'/,                  ! EXIT 11
     & c( 9)   /'infeasible linear equality constraints'/,         ! EXIT 12
     & c(10)   /'nonlinear infeasibilities minimized'/,            ! EXIT 13
     & c(11)   /'linear infeasibilities minimized'/,               ! EXIT 14
     & c(12)   /'infeasible linear constraints in QP subproblem'/, ! EXIT 15
     & c(13)   /'infeasible nonelastic constraints'/,              ! EXIT 16
     & c(14)/'the problem appears to be unbounded'/,               ! EXIT 20
     & c(15)   /'unbounded objective'/,                            ! EXIT 21
     & c(16)   /'constraint violation limit reached'/,             ! EXIT 22
     & c(17)/'resource limit error'/,                              ! EXIT 30
     & c(18)   /'iteration limit reached'/,                        ! EXIT 31
     & c(19)   /'major iteration limit reached'/,                  ! EXIT 32
     & c(20)   /'the superbasics limit is too small'/,             ! EXIT 33
     & c(21)   /'time limit reached'/,                             ! EXIT 34
     & c(22)/'terminated after numerical difficulties'/,           ! EXIT 40
     & c(23)   /'current point cannot be improved'/,               ! EXIT 41
     & c(24)   /'singular basis'/,                                 ! EXIT 42
     & c(25)   /'cannot satisfy the general constraints'/,         ! EXIT 43
     & c(26)   /'ill-conditioned null-space basis'/,               ! EXIT 44
     & c(27)   /'unable to compute acceptable LU factors'/,        ! EXIT 45
     & c(28)/'error in the user-supplied functions'/,              ! EXIT 50
     & c(29)   /'incorrect objective  derivatives'/,               ! EXIT 51
     & c(30)   /'incorrect constraint derivatives'/,               ! EXIT 52
     & c(31)   /'the QP Hessian is indefinite'/,                   ! EXIT 53
     & c(32)   /'incorrect second derivatives'/,                   ! EXIT 54
     & c(33)   /'incorrect derivatives'/,                          ! EXIT 55
     & c(34)/'undefined user-supplied functions'/,                 ! EXIT 60
     & c(35)   /'undefined function at the first feasible point'/, ! EXIT 61
     & c(36)   /'undefined function at the initial point'/,        ! EXIT 62
     & c(37)   /'unable to proceed into undefined region'/,        ! EXIT 63
     & c(38)/'user requested termination'/,                        ! EXIT 70
     & c(39)   /'terminated during function evaluation'/,          ! EXIT 71
     & c(40)   /'terminated during constraint evaluation'/,        ! EXIT 72
     & c(41)   /'terminated during objective evaluation'/,         ! EXIT 73
     & c(42)   /'terminated from monitor routine'/,                ! EXIT 74
     & c(43)/'insufficient storage allocated'/,                    ! EXIT 80
     & c(44)   /'work arrays must have at least 500 elements'/,    ! EXIT 81
     & c(45)   /'not enough character storage'/,                   ! EXIT 82
     & c(46)   /'not enough integer storage'/,                     ! EXIT 83
     & c(47)   /'not enough real storage'/,                        ! EXIT 84
     & c(48)/'input arguments out of range'/,                      ! EXIT 90
     & c(49)   /'invalid input argument'/,                         ! EXIT 91
     & c(50)   /'basis file dimensions do not match this problem'/,! EXIT 92
     & c(51)   /'the QP Hessian is indefinite'/,                   ! EXIT 93
     & c(52)/'finished successfully'/,                             ! EXIT100
     & c(53)   /'SPECS file read'/,                                ! EXIT101
     & c(54)   /'Jacobian structure estimated'/,                   ! EXIT102
     & c(55)   /'MPS file read'/,                                  ! EXIT103
     & c(56)   /'memory requirements estimated'/,                  ! EXIT104
     & c(57)   /'user-supplied derivatives appear to be correct'/, ! EXIT105
     & c(58)   /'no derivatives were checked'/,                    ! EXIT106
     & c(59)   /'some SPECS keywords were not recognized'/,        ! EXIT107
     & c(60)/'errors while processing MPS data'/,                  ! EXIT110
     & c(61)   /'no MPS file specified'/,                          ! EXIT111
     & c(62)   /'problem-size estimates too small'/,               ! EXIT112
     & c(63)   /'fatal error in the MPS file'/,                    ! EXIT113
     & c(64)/'errors while estimating Jacobian structure'/,        ! EXIT120
     & c(65)   /'cannot find Jacobian structure at given point'/,  ! EXIT121
     & c(66)/'fatal errors while reading the SPECS'/,              ! EXIT130
     & c(67)   /'no SPECS file (iSpecs le 0 or iSpecs gt 99)'/,    ! EXIT131
     & c(68)   /'End-of-file while looking for a BEGIN'/,          ! EXIT132
     & c(69)   /'End-of-file while reading SPECS file'/,           ! EXIT133
     & c(70)   /'ENDRUN found before any valid SPECS'/,            ! EXIT134
     & c(71)/'system error'/,                                      ! EXIT140
     & c(72)   /'wrong no of basic variables'/,                    ! EXIT141
     & c(73)   /'error in basis package'/                          ! EXIT142
     & c(74)   /'Problem dimensions are too large'/                ! EXIT143
!     ------------------------------------------------------------------
!     Find the "major" and "minor" iExit modes

      mjr = iExit/10
      mnr = iExit - 10*mjr

      i1  = indc(mjr)
      i2  = i1 + mnr

      call s1trim( c(i1), length )
      write(string , '(1x,2a,i4,a,(a))')
     &     Solver, ' EXIT', 10*mjr, ' -- ', c(i1)(1:length)

      call s1trim( c(i2), length )
      write(string2, '(1x,2a,i4,a,(a))')
     &     Solver, ' INFO',  iExit, ' -- ', c(i2)(1:length)

      end ! subroutine snEXIT

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snWRAP
     &   ( iExit, Solver, string, string2, iw, leniw )

      implicit
     &     none
      character*(*)
     &     string, string2
      character
     &     Solver*6
      integer
     &     iExit, leniw, iw(leniw)

!     ==================================================================
!     snWRAP  it's a wrap!
!
!     18 Oct 2003: First version of snWRAP.
!     21 Jun 2004: iExit saved for use in snSolF.
!     21 Jun 2004: Current version of snWRAP.
!     ==================================================================
      integer
     &     iExit0
!     ------------------------------------------------------------------

      iw(424) = iExit   ! INFO code from all solvers

      call snEXIT( iExit, Solver, string, string2 )

      if (iExit .eq. 81) then   ! Print without using accessing iw, etc.
         call snPRNT( 15, string , iw, leniw )
         call snPRNT(  5, string2, iw, leniw )
      else
         iExit0 = iExit/10
         if (iExit0 .eq. 0) then ! Normal exit
            call s1page( 1, iw, leniw )
         else
            call s1page( 2, iw, leniw )
         end if

         call snPRNT( 3, string , iw, leniw )
         call snPRNT( 3, string2, iw, leniw )
      end if

      end ! subroutine snWRAP

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine snSolF
     &   ( m, n, nb, ninf, j, jkey, jstate,
     &     hs, bl, bu, rc, xs, iw, leniw, rw, lenrw )

      implicit
     &     none
      integer
     &     m, n, nb, ninf, j, jkey, jstate, leniw, lenrw
      integer
     &     hs(nb), iw(leniw)
      double precision
     &     bl(nb), bu(nb), rc(nb), xs(nb), rw(lenrw)

!     ==================================================================
!     snSolF sets the solution flags for the j-th variable:
!                      ' ' A  D  I  N    and   LL  UL SBS  BS  EQ  FR
!     by returning jkey=0  1  2  3  4,  jstate= 0   1   2   3   4   5
!
!     snSolF is called by SNOPT from s4soln.
!     snSolF may also be called externally (e.g. by GAMS)
!     following a normal call of SNOPT.
!     At this stage the solution will be UNSCALED!!
!     Hence, SNOPT (via m4soln) now outputs flags for the UNSCALED soln.
!
!     Input parameters m, n, nb, hs, bl, bu, rc, xs
!     are the same as for SNOPT.
!
!     j      (input ) is column j if j <= n;  otherwise row i = j - n.
!     jkey   (output) is one of 0 1 2 3 4.
!     jstate (output) is one of 0 1 2 3 4 5.
!
!     09 Mar 2004: First version of snSolF, derived from misolf.
!     18 Jun 2004: If the scaled problem was infeasible
!                  (with max inf at j = jbInf1), always flag that j
!                  as infeasible in the unscaled solution.
!     21 Jun 2004: Similarly, if the scaled problem wasn't optimal,
!                  (with max dual inf at j = jdInf1), always flag that j
!                  as nonoptimal in the unscaled solution.
!     23 Jun 2004: Suppress nonoptimal flag if iExit <= 2 (not 0).
!     ==================================================================

      logical
     &     Feasible, maximz
      integer
     &     iExit, jbInf1, jdInf1, js, minimize
      double precision
     &     b1, b2, d1, d2, dj, djtest, piNorm,
     &     tolfea, tolNLP, tolopt, tolx, xj


      minimize = iw(199) ! (-1)(+1)    => (max)(min)
      iExit    = iw(424) ! INFO code from all solvers
      jbInf1   = iw(427) ! Largest bound infeasibility (  scaled)
      jdInf1   = iw(428) ! Largest dual  infeasibility (  scaled)

      tolNLP   = rw( 53) ! Major Optimality tolerance
      tolx     = rw( 56) ! Minor feasibility tolerance
      piNorm   = rw(422) ! Lagrange multiplier norm

      tolfea   = tolx
      tolopt   = tolNLP * piNorm
      Feasible = ninf     .eq. 0
      maximz   = minimize .lt. 0

      js       = hs(j)
      b1       =  bl(j)
      b2       = bu(j)
      xj       = xs(j)
      dj       = rc(j)
      d1       = b1 - xj
      d2       = xj - b2
      djtest   = - dj

      if (Feasible) then
         if (maximz) djtest = - djtest
         jbInf1 = 0
      end if

      if (iExit .le. 2) then
         jdInf1 = 0
      end if

      ! Set keys and states.

      jkey   = 0   ! blank
      jstate = js  ! 0, 1, 2, 3

      if (js .le. 1) then            ! Nonbasic variables.
         if (b1 .eq. b2) jstate = 4
         if (- d1 .gt. tolfea  .and.     - d2 .gt. tolfea) jstate = 5
         if (jstate .eq. 1 ) djtest = - djtest
         if (jstate .ge. 4 ) djtest =   abs(djtest)
         if (                     abs(djtest) .le. tolopt) jkey = 1  ! A
         if (jstate .ne. 4     .and.  djtest  .gt. tolopt) jkey = 4  ! N

      else                           ! Basic and superbasic variables.
         if (abs(d1).le. tolfea  .or. abs(d2) .le. tolfea) jkey = 2  ! D
         if (jstate .eq. 2  .and. abs(djtest) .gt. tolopt) jkey = 4  ! N
         if (    d1 .gt. tolfea  .or.     d2  .gt. tolfea) jkey = 3  ! I
         if (     j .eq. jbInf1                          ) jkey = 3  ! I
      end if

      if (j .eq. jdInf1) jkey = 4  ! N

      end ! subroutine snSolF
