C12345678901234567890123456789012345678901234567890123456789012345678901
C
C     $Id: xmxm.f,v 8.1 2015/04/20 11:14:20 hmb Exp $
C
C     Matrix-matrix, matrix-vector multiply routines,
C     designed to be called from C.
C     E.g. where C = A * B; the FORTRAN equivalent is C' = B' * A'.
C
C     Copyright (c) 1998 <--> $Date: 2015/04/20 11:14:20 $, Hugh Blackburn
C
C     This file is part of Semtex.
C 
C     Semtex is free software; you can redistribute it and/or modify it
C     under the terms of the GNU General Public License as published by
C     the Free Software Foundation; either version 2 of the License, or
C     (at your option) any later version.
C 
C     Semtex is distributed in the hope that it will be useful, but
C     WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C     General Public License for more details.
C 
C     You should have received a copy of the GNU General Public License
C     along with Semtex (see the file COPYING); if not, write to the
C     Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
C     Boston, MA 02110-1301 USA
C
C
C     ------------------------------------------------------------------
C     C = A * B.  (As written in C, with row-major storage.)
C
      SUBROUTINE DMXM (A, NRA, B, NCA, C, NCB)
      IMPLICIT NONE
      INTEGER          NRA, NCA, NCB, I, J, K
      DOUBLE PRECISION A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
      DO J = 1, NCB
         DO I = 1, NRA
            C(J, I) = 0.0D0
            DO K = 1, NCA
               C(J, I) = C(J, I) + B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE SMXM (A, NRA, B, NCA, C, NCB)
      IMPLICIT NONE
      INTEGER  NRA, NCA, NCB, I, J, K
      REAL     A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
      DO J = 1, NCB
         DO I = 1, NRA
            C(J, I) = 0.0
            DO K = 1, NCA
               C(J, I) = C(J, I) + B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
C
C     ------------------------------------------------------------------
C     C += A * B.
C
      SUBROUTINE DMXMA (A, NRA, B, NCA, C, NCB)
      IMPLICIT NONE
      INTEGER          NRA, NCA, NCB, I, J, K
      DOUBLE PRECISION A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
      DO J = 1, NCB
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) + B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE SMXMA (A, NRA, B, NCA, C, NCB)
      IMPLICIT NONE
      INTEGER  NRA, NCA, NCB, I, J, K
      REAL     A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
      DO J = 1, NCB
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) + B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
C
C     ------------------------------------------------------------------
C     C -= A * B.
C
      SUBROUTINE DMXMS (A, NRA, B, NCA, C, NCB)
      IMPLICIT NONE
      INTEGER          NRA, NCA, NCB, I, J, K
      DOUBLE PRECISION A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
      DO J = 1, NCB
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) - B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE SMXMS (A, NRA, B, NCA, C, NCB)
      IMPLICIT NONE
      INTEGER  NRA, NCA, NCB, I, J, K
      REAL     A(NCA, NRA), B(NCB, NCA), C(NCB, NRA)
      DO J = 1, NCB
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) - B(J, K) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
C
C     ------------------------------------------------------------------
C               t
C     C -= A * B.
C
      SUBROUTINE DMXMTS (A, NRA, B, NCA, C, NCBT)
      IMPLICIT NONE
      INTEGER          NRA, NCA, NCBT, I, J, K
      DOUBLE PRECISION A(NCA, NRA), B(NCA, NCBT), C(NCBT, NRA)
      DO J = 1, NCBT
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) - B(K, J) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE SMXMTS (A, NRA, B, NCA, C, NCBT)
      IMPLICIT NONE
      INTEGER  NRA, NCA, NCBT, I, J, K
      REAL     A(NCA, NRA), B(NCA, NCBT), C(NCBT, NRA)
      DO J = 1, NCBT
         DO I = 1, NRA
            DO K = 1, NCA
               C(J, I) = C(J, I) - B(K, J) * A(K, I)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
C
C     ------------------------------------------------------------------
C     Y = A * X.  Matrix-vector multiply.
C
      SUBROUTINE DMXV (A, NRA, X, NCA, Y)
      IMPLICIT NONE
      INTEGER          NRA, NCA, I, J
      DOUBLE PRECISION A(NCA, NRA), X(NRA), Y(NCA)
C     
C     -- ALTERNATIVE ORDERING FOR TESTING:
C
C      DO I = 1, NCA
C         Y(I) = 0.0D0
C      ENDDO
C      
C      DO J = 1, NRA
C         DO I = 1, NCA
C            Y(I) = Y(I) + A(J, I) * X(J)
C         ENDDO
C      ENDDO
      DO I = 1, NCA
         Y(I) = 0.0D0
         DO J = 1, NRA
            Y(I) = Y(I) + A(J, I) * X(J)
         ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE SMXV (A, NRA, X, NCA, Y)
      IMPLICIT NONE
      INTEGER NRA, NCA, I, J
      REAL    A(NCA, NRA), X(NRA), Y(NCA)
      DO I = 1, NCA
         Y(I) = 0.0
         DO J = 1, NRA
            Y(I) = Y(I) + A(J, I) * X(J)
         ENDDO
      ENDDO
      RETURN
      END
