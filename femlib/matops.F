C12345678901234567890123456789012345678901234567890123456789012345678901
C
C     -- Matrix operators for gradients, etc.
C
C     $Id: matops.F,v 8.1 2015/04/20 11:14:14 hmb Exp $
C
C     Copyright (c) 2002 <--> $Date: 2015/04/20 11:14:14 $, Hugh Blackburn
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
      SUBROUTINE DGRAD2 (X, Y, UR, US, DV, DT, NR, NS, NEL)
C
C     This operator can be used for any field operation that can be 
C     described by the tensor-product forms
C
C     UR =  X    [DT] (This is a derivative in the "r" direction.)
C     US =  [DV]  Y   (This is a derivative in the "s" direction.)
C
C     Often for tensor-product forms, and when NR=NS, DT = transpose
C     (DV), but this restriction is not assumed here.  A typical
C     application is for taking derivatives: Ur = dX/dr, Us = dY/ds over
C     all nel elements each np*np.  The routine is designed for use with
C     row-major (C-like) matrices.
C
C     UR & US are updated by addition, so should be zero on entry.
C     Conceivably UR & US can be aliased to one another, so
C     the operations are broken out into two loops.
C
C     Matrix X:  NR x NS (x NEL)
C     Matrix Y:  NR x NS (x NEL)
C     Matrix DV: NS x NS
C     Matrix DT: NR x NR
C
      IMPLICIT NONE
C
      INTEGER          NR, NS, NEL, I, J, K, L
      DOUBLE PRECISION  X(NR,NS,NEL),  Y(NR,NS,NEL)
      DOUBLE PRECISION UR(NR,NS,NEL), US(NR,NS,NEL)
      DOUBLE PRECISION DV(NS,NS),     DT(NR,NR)
C
#if defined (_SX) || defined (__uxp__) || defined (_VECTOR_ARCH)
      IF (NEL.EQ.1) THEN
         CALL DGEMM ('N','N',NR,NS,NR,1.0D0,DT,NR, X,NR,1.0D0,UR,NR)
         CALL DGEMM ('N','N',NR,NS,NS,1.0D0,Y, NR,DV,NS,1.0D0,US,NR)
      ELSE
         DO K = 1, NR
            DO J = 1, NS
               DO L = 1, NR
                  DO I = 1, NEL
                     UR(K,J,I) = UR(K,J,I) + X(L,J,I) * DT(K,L)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
#if defined (_SX)
!CDIR NOVECTOR
#endif
         DO J = 1, NS
#if defined (_SX)
!CDIR NOVECTOR
#endif
            DO L = 1, NS
#if defined (_SX)
!CDIR NOVECTOR
#endif
               DO K = 1, NR
                  DO I = 1, NEL
                     US(K,J,I) = US(K,J,I) + Y(K,L,I) * DV(L,J)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
#else
      CALL DGEMM ('N','N',NR,NS*NEL,NR,1.0D0,DT,NR,X,NR,1.0D0,UR,NR)
      DO I = 1, NEL
         CALL DGEMM ('N','N',NR,NS,NS,1.0D0,Y(1,1,I),
     +        NR,DV,NS,1.0D0,US(1,1,I),NR)
      ENDDO
#endif
      RETURN
      END
C
C
       SUBROUTINE DTPR2D (X, Y, T, DV, DT, NR, NS, NEL)
C
C     This operator can be used for any field operation that can be 
C     described (in C, row-major ordering) by the tensor-product form
C
C     Y  = [DV] X [DT].
C
C     X and Y are both square: X is NS x NS, Y is NR x NR.
C     X and Y are not required to be distinct if NR = NS (i.e. in-place).
C
C     Typically for tensor-product forms, DT = transpose (DV), but this
C     is not enforced here.  A typical application is for projections of
C     one element order to another.  The routine is designed for use
C     with row-major (C-like) matrices.
C     
C     Matrix X:  NS x NS (x NEL)
C     Matrix Y:  NR x NR (x NEL)
C     Matrix T:  NS x NR (x NEL, workspace)
C     Matrix DV: NR x NS
C     Matrix DT: NS x NR
C
      IMPLICIT NONE
C
      INTEGER          NTOT, NR, NS, NEL, I, J, K, L
      DOUBLE PRECISION  X(NS,NS,NEL),  Y(NR,NR,NEL), T(NR,NS,NEL)
      DOUBLE PRECISION DV(NS,NR),     DT(NR,NS)
C
#if defined (_SX) || defined (__uxp__) || defined (_VECTOR_ARCH)
      IF (NEL.EQ.1) THEN
         CALL DGEMM ('N','N',NR,NS,NS,1.0D0,DT,NR, X,NS,0.0D0,T,NR)
         CALL DGEMM ('N','N',NR,NR,NS,1.0D0, T,NR,DV,NS,0.0D0,Y,NR)
      ELSE
         NTOT = NR*NS*NEL
         DO I = 1, NTOT
            T(I,1,1) = 0.0D0
         ENDDO
#if defined (_SX)
!CDIR NOVECTOR
#endif
         DO J = 1, NS
#if defined (_SX)
!CDIR NOVECTOR
#endif
            DO L = 1, NS
#if defined (_SX)
!CDIR NOVECTOR
#endif
               DO K = 1, NR
                  DO I = 1, NEL
                     T(K,J,I) = T(K,J,I) + X(L,J,I) * DT(K,L)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         NTOT = NR*NR*NEL
         DO I = 1, NTOT
            Y(I,1,1) = 0.0D0
         ENDDO
         DO K = 1, NR
            DO J = 1, NR
               DO L = 1, NS
                  DO I = 1, NEL
                     Y(K,J,I) = Y(K,J,I) + T(K,L,I) * DV(L,J)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
#else
      CALL DGEMM ('N','N',NR,NS*NEL,NS,1.0D0,DT,NR,X,NS,0.0D0,T,NR)
      DO I = 1, NEL
         CALL DGEMM ('N','N',NR,NR,NS,1.0D0,T(1,1,I),
     +        NR,DV,NS,0.0D0,Y(1,1,I),NR)
      ENDDO
#endif
      RETURN
      END
