*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*
*     File  sn83bfgs.f.   Full and limited memory BFGS routines.
*
*     s8FMH0   s8FMupdate   s8FMHx
*     s8LMH0   s8LMupdate   s8LMHx
*
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8FMH0
     &   ( nnH, HD, lenU, U )

      implicit
     &     none
      integer
     &     lenU, nnH
      double precision
     &     HD(nnH), U(lenU)

*     ==================================================================
!     s8FMH0 zeros the off-diagonal elements of H such that H = U'U.
!     The diagonals of U are set to the square roots of the diagonal HD.
*
*     19 Jul 1995: First version of s8FMH0.
*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added.
*     13 Jan 2005: HD always positive semidefinite.
*     05 Oct 2014: Reorganized to reflect qnInit in dnopt.
*     ==================================================================
      integer
     &     incr, j, l, nzeros
*     ------------------------------------------------------------------
      double precision   zero
      parameter         (zero = 0.0d+0)
*     ------------------------------------------------------------------
*     Zero the off-diagonal elements of U.
*     ------------------------------------------------------------
      incR   = nnH
      nzeros = nnH - 1
      l      = 1

      do  j   = 1, nnH
         U(l) = sqrt(HD(j))
         if (j .lt. nnH) then
            call dload ( nzeros, zero, U(l+1), 1 )
            l      = l      + incR
            incR   = incR   - 1
            nzeros = nzeros - 1
         end if
      end do

      end ! subroutine s8FMH0

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8FMupdate
     &   ( updateType, nnH, U0scale,
     &     rydx, rdxHdx, Hdx, y, Udx, lenU, U )

      implicit
     &     none
      integer
     &     updateType, lenU, nnH
      double precision
     &     rydx, rdxHdx, U0scale, Hdx(nnH), y(nnH), Udx(nnH), U(lenU)

*     ==================================================================
*     s8FMupdate applies the full-memory BFGS update to H = U'U.
*     If defined, the self-scaling BFGS update parameter is saved.
*     It is needed to update the reduced Hessian when there are only
*     linear constraints.
*
*     19 Jul 1995: First version of s8FMup.
*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added.
*     18 Feb 2001: LM H stored in product form.
*     13 Jan 2005: FM H stored in product form.
*     15 Jan 2005: Current version.
*     ==================================================================
      external
     &     ddot
      integer
     &     iExit, lastnz, numU
      double precision
     &     ddot, t, told, tolz, Ulast
*     ------------------------------------------------------------------
      integer            ssBFGS
      parameter         (ssBFGS = 1)

      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one    = 1.0d+0)
*     ------------------------------------------------------------------
      told   = zero
      tolz   = zero

      if (updateType .eq. ssBFGS) then
         numU   = nnH*(nnH + 1)/2
         call dscal ( numU, U0scale, U  , 1 ) ! multiplies U by U0scale.
         call dscal ( nnH , U0scale, Udx, 1 )
      end if

*     ------------------------------------------------------------------
*     Overwrite (Udx,y) with the vectors (Us,v)  such that
*       Us  = Udx / rdxHdx,    v = (1/rydx) gdif - (1/rdxHdx) Hdx.
*
*     Then, U(new) = U + Us v',  with H = U'U.
*
*     Hdx and v  are saved to update R for LC problems.
*     ------------------------------------------------------------------
      t      = ddot ( nnH, y, 1, Hdx, 1 )
      if (t .ge. zero) then
         call dscal ( nnH, ( one/rydx), y, 1 )
      else
         call dscal ( nnH, (-one/rydx), y, 1 )
      end if

      call daxpy ( nnH, (-one/rdxHdx), Hdx, 1, y, 1 )
      call dscal ( nnH, ( one/rdxHdx), Udx, 1 )

*     ------------------------------------------------------------------
*     Restore  U + Us y' to triangular form  (overwriting Udx).
*     ------------------------------------------------------------------
      Ulast  = zero
      lastnz = nnH

      call s6Rmod
     &   ( iExit, nnH, nnH, lenU, U, Udx, y, lastnz, Ulast,
     &     told, tolz )

      end ! subroutine s8FMupdate

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8FMHx
     &   ( nnH, x, Ux, Hx, lenU, U )

      implicit
     &     none
      integer
     &     nnH, lenU
      double precision
     &     U(lenU), Hx(nnH), Ux(nnH), x(nnH)

*     ==================================================================
*     s8FMHx  computes the product Hx = U'Ux, where  U is an upper-
*     triangular matrix stored by rows in the one-dimensional array  U.
*     lenU defines the length of U.  lenU must be at least
*     nnH*(nnH + 1)/2.
*
*     12 Jan 1996: First version of s8FMHx
*     12 Jan 2005: H held as U'U.
*     15 Jan 2005: Current version.
*     ==================================================================
      integer            WithU,      WithUt
      parameter         (WithU  = 0, WithUt = 1)
*     ------------------------------------------------------------------
      call s6Rprod
     &   ( WithU , nnH, nnH, lenU, U,  x, Ux )
      call s6Rprod
     &   ( WithUt, nnH, nnH, lenU, U, Ux, Hx )

      end ! subroutine s8FMHx

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8LMH0
     &   ( nnH, HD, U0 )

      implicit
     &     none
      integer
     &     nnH
      double precision
     &     HD(nnH), U0(nnH)

*     ==================================================================
*     s8LMH0  zeros the off-diagonal elements of H such that H = U'U.
!     The diagonals of U are set to the square roots of the diagonal HD.
!
*     19 Jul 1995: First version of s8LMH0.
*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added.
*     18 Feb 2001: H stored in product form.
*     13 Jan 2005: Hd always positive semidefinite.
*     05 Oct 2014: Reorganized to reflect qnInit in dnopt.
*     ==================================================================
      integer
     &     i
*     ------------------------------------------------------------------
      do     i = 1, nnH
         U0(i) = sqrt(Hd(i))
      end do

      end ! subroutine s8LMH0

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8LMupdate
     &   ( updateType, QNmods, mQNmods, nnH, U0scale,
     &     rydx, rdxHdx, Hdx, y, dx, U0, S, V )

      implicit
     &     none
      integer
     &     updateType, QNmods, mQNmods, nnH
      double precision
     &     rydx, rdxHdx, U0scale, dx(nnH),
     &     Hdx(nnH), y(nnH), U0(nnH), S(nnH,mQNmods),
     &     V(nnH,mQNmods)

*     ==================================================================
*     s8LMupdate computes the limited-memory BFGS update.
*
*     If defined, the self-scaling BFGS parameter U0scalee is saved.
*     It is needed to update the reduced Hessian when there are only
*     linear constraints.
*
*     19 Jul 1995: First version of s8LMup.
*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added.
*     18 Feb 2001: LM H stored in product form.
*     13 Jan 2005: FM H stored in product form.
*     16 Jan 2005: Current version.
*     ==================================================================
      external
     &     ddot
      double precision
     &     ddot, t
*     ------------------------------------------------------------------
      integer            ssBFGS
      parameter         (ssBFGS = 1)

      double precision   zero,            one
      parameter         (zero   = 0.0d+0, one    = 1.0d+0)
*     ------------------------------------------------------------------
      if (updateType .eq. ssBFGS) then
         call dscal ( nnH, U0scale, U0 , 1 ) ! multiplies U0 by U0scale.
      end if

*     ------------------------------------------------------------------
*     Space remains. Store s and v, where
*     U(new) = U(I + sv'), with H = U'U.
*       S  =  dx / rdxHdx,    V =  (1/rydx) gdif - (1/rdxHdx) Hdx.
*
*     Hdx and the modified y (= v) are used to update the reduced
*     Hessian for LC problems.
*     ------------------------------------------------------------------
      call dcopy ( nnH,         dx, 1, S(1,QNmods), 1 )
      call dscal ( nnH, ( one/rdxHdx), S(1,QNmods), 1 )

      t = ddot ( nnH, y, 1, Hdx, 1 )
      if (t .ge. zero) then
         call dscal ( nnH, ( one/rydx), y, 1 )
      else
         call dscal ( nnH, (-one/rydx), y, 1 )
      end if

      call daxpy ( nnH, (-one/rdxHdx), Hdx, 1,           y, 1 )
      call dcopy ( nnH,                  y, 1, V(1,QNmods), 1 )

      end ! subroutine s8LMupdate

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8LMHx
     &   ( nnH, x, Ux, Hx, mQNmods, QNmods, U0, S, V )

      implicit
     &     none
      integer
     &     nnH, mQNmods, QNmods
      double precision
     &     Hx(nnH), U0(nnH), Ux(nnH), x(nnH),
     &     S(nnH,mQNmods), V(nnH,mQNmods)

*     ==================================================================
*     s8LMHx forms the product  Hx  for the limited-memory s8Hx.
*     H = U'U, where U = U0*(I + s1*v1')*(I + s2*v2')...(I + sk*vk').
*     with  k = QNmods
*
*     19 Jul 1995: First version of s8LMHx
*     18 Feb 2001: H stored in product form.
*     12 Jan 2005: Ux added as argument.
*     16 Jan 2005: Current version.
*     ==================================================================
      external
     &     ddot
      integer
     &     k
      double precision
     &     c, ddot
*     ------------------------------------------------------------------
      call dcopy ( nnH,  x, 1, Ux, 1 )

*     Multiply by U.

      do k = QNmods, 1, -1
         c = ddot   ( nnH,    V(1,k), 1, Ux, 1 )
         call daxpy ( nnH, c, S(1,k), 1, Ux, 1 )
      end do

      call ddscl ( nnH, U0, 1, Ux, 1 )

*     Multiply by U'.

      call dcopy ( nnH, Ux, 1, Hx, 1 )
      call ddscl ( nnH, U0, 1, Hx, 1 )

      do k = 1, QNmods
         c = ddot   ( nnH,    S(1,k), 1, Hx, 1 )
         call daxpy ( nnH, c, V(1,k), 1, Hx, 1 )
      end do

      end ! subroutine s8LMHx
