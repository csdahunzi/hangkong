!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     File  sn85Hess.f
!
*     s8SDH0   s8SDIx
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8SDH0
     &   ( HQNType, nnH, U0pre, U0 )

      implicit
     &     none
      integer
     &     HQNType, nnH
      double precision
     &     U0pre, U0(nnH)

*     ==================================================================
*     s8SDH0 resets the approximate Hessian H to a diagonal matrix.
*     On entry, the value of HQNType is as follows:
*
*       HQNType
*       -------
*       HUnset (-1)      H not set.
*       HNorml ( 0)      H is a Hessian of the form defined by  lvlHes.
*       HDiag  ( 1)      H is a diagonal matrix.
*       HUnit  ( 2)      H is an identity matrix.
*
*     19 Jul 1995: First version of s8SDH0.
*     06 Sep 1998: Pre- and post-QP diagonal Hessian scaling added.
*     25 Mar 2005: Resurrected for snopt8.
*     25 Mar 2005: Current version.
*     ==================================================================
      integer            HUnit
      parameter         (HUnit  = 2)
*     ------------------------------------------------------------------
*     Set H0 to a multiple of the identity.
*     ------------------------------------------------------------------
      HQNType = HUnit
      call dload ( nnH, U0pre, U0, 1 )

      end ! subroutine s8SDH0

*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine s8SDIx
     &   ( nnH, U0, x, Hx )

      implicit
     &     none
      integer
     &     nnH
      double precision
     &     U0(nnH), Hx(nnH), x(nnH)

*     ==================================================================
*     s8SDIx  multiplies the QP Hessian  H by the vector  x.
*     It is used to define Hx for the first QP subproblem.
*
*     19 Jul 1995: First version of s8SDIx.
*     13 Apr 2005: Current version.
*     ==================================================================

      call dcopy ( nnH,  x, 1, Hx, 1 )
      call ddscl ( nnH, U0, 1, Hx, 1 )
      call ddscl ( nnH, U0, 1, Hx, 1 )

      end ! subroutine s8SDIx

