C12345678901234567890123456789012345678901234567890123456789012345678901
C     ==================================================================
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

C     Routines for generalized multiple 1D prime factor complex--complex
C     self-sorting in-place FFT, and wrapper routine  for  real--complex
C     FFT.
C
C     CALL PRF235 (N, IP, IQ, IR, IPQR2)
C     CALL DSETPF (TRIG, IP, IQ, IR)
C     CALL DMPFFT (V, W, NP, NZ, IP, IQ, IR, TRIG, ISIGN)
C
C     Direct call to Temperton's underlying complex--complex FFT:
C     CALL DGPFA  (A,B,TRIG,INC,JUMP,N,IP,IQ,IR,LOT,ISIGN) 
C     
C     Reference:
C     C. Temperton, "A generalized prime factor FFT algorithm
C       for any N = (2**P)(3**Q)(5**R)",
C       SIAM J. Sci. Stat. Comp., May 1992.
C
C     The complex--complex FFT routines are the original GPFA code
C     written by Clive Temperton, Recherche en Prevision Numerique
C     Atmospheric Environment Service, Canada.
C
C     Note that the definition of the DFT employed in DGPFA is
C
C        X(J) = SUM(K=0,...,N-1)(C(K)*EXP(ISIGN*2*I*J*K*PI/N)).
C
C     This file contains DOUBLE PRECISION versions of routines.
C
C     $Id: temfftd.F,v 8.1 2015/04/20 11:14:14 hmb Exp $
C     ==================================================================
C
C
      SUBROUTINE PRF235 (N, IP, IQ, IR, IPQR2)
C     ------------------------------------------------------------------
C     Compute prime factors of N (factors of 2, 3, and 5 are allowed).
C     Return the factors in the form IP, IQ, IR, where
C       N = (2**IP) * (3**IQ) * (5**IR)
C
C     Also return IPQR2 (storage size for trigonometric array)
C       IPQR2 = 2 * ((2**IP) + (3**IQ) + (5**IR))
C
C     If N does not have suitable prime factors, it is returned as zero.
C     ------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER N, IP, IQ, IR, IPQR2
      INTEGER IFAC, KK, LL, NN, NJ(3)
C
      NN   = N
      IFAC = 2
C     
      DO 30 LL = 1, 3
         KK = 0
 10      CONTINUE
         IF (MOD(NN,IFAC).NE.0) GO TO 20
         KK = KK + 1
         NN = NN / IFAC
         GO TO 10
 20      CONTINUE
         NJ(LL) = KK
         IFAC   = IFAC + LL
 30   CONTINUE
C
      IF (NN.NE.1) THEN
         N = 0
         RETURN
      ENDIF
C
      IP = NJ(1)
      IQ = NJ(2)
      IR = NJ(3)
      IPQR2 = 2**IP + 3**IQ + 5**IR
      IPQR2 = IPQR2 + IPQR2
C
      RETURN
      END
C
C
      SUBROUTINE DSETPF (TRIG,N,IP,IQ,IR)
C     ------------------------------------------------------------------
C     Setup routine for self-sorting in-place
C     generalized prime factor (complex) FFT [GPFA]
C
C     Usage: CALL DSETPF (TRIG,N,IP,IQ,IR)
C
C     Input :
C     -----
C     N is the length of the transforms. N must be of the form:
C       -------------------------------
C       N = (2**IP) * (3**IQ) * (5**IR)
C       -------------------------------
C     IP, IQ, IR are also input (call PRF235 to compute).
C
C     Output:
C     ------
C     Trig is a table of twiddle factors,
C     of length 2*IPQR (real) words, where:
C       --------------------------------------
C       IPQR = (2**IP) + (3**IQ) + (5**IR)
C       --------------------------------------
C     ------------------------------------------------------------------
      IMPLICIT NONE
C
      DOUBLE PRECISION TRIG(*)
      INTEGER          N, IP, IQ, IR
      DOUBLE PRECISION DEL, TWOPI, ANGLE
      INTEGER          NJ(3), LL, KK
      INTEGER          I, K, NI, IROT, KINK
      PARAMETER        ( TWOPI = 6.28318530717958647688D0 )
C
      NJ(1) = 2**IP
      NJ(2) = 3**IQ
      NJ(3) = 5**IR
C
      I = 1
C
      DO 60 LL = 1, 3
         NI = NJ(LL)
         IF (NI.EQ.1) GO TO 60
C
         DEL  = TWOPI / DFLOAT(NI)
         IROT = N / NI
         KINK = MOD(IROT,NI)
         KK   = 0
C
         DO 50 K = 1, NI
            ANGLE     = DFLOAT(KK) * DEL
            TRIG(I)   = DCOS(ANGLE)
            TRIG(I+1) = DSIN(ANGLE)
            I  = I + 2
            KK = KK + KINK
            IF (KK.GT.NI) KK = KK - NI
 50      CONTINUE
 60   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE DMPFFT (V,W,NP,NZ,IP,IQ,IR,TRIG,ISIGN)
C     ------------------------------------------------------------------
C     Compute multiple 1D real--complex transforms of V, 
C     return in place, using generalized prime factor FFT.
C
C     V     is the 2D array to be transformed, size NP * NZ;
C     W     is a work array, size NP * NZ;
C     NP    is the number of transforms (must be even);
C     NZ    is the number of real data in each transform, 
C           factorizable by 2, 3, 5:
C           --------------------------------
C           NZ = (2**IP) * (3**IQ) * (5**IR);
C           --------------------------------
C     TRIG  is an array of trigonometric factors, length
C           2*IPQR (real) words, where:
C           ----------------------------------
C           IPQR = (2**IP) + (3**IQ) + (5**IR);
C           ----------------------------------
C     ISIGN specifies the direction of the transform.
C
C     IP, IQ, IR and TRIG should be precomputed by DSETPF.
C     
C     For the real-->complex transform (ISIGN=+1), V is interpreted on 
C     input as a 2D array with each row as a distinct data set: DFTs are
C     to be carried out along each of the NP rows, each of which has NZ
C     points in it.
C
C     After the transform, the real and imaginary parts of the V~ are
C     stored alternating in the rows of V.  For NZ even, the Nyquist
C     frequency datum is stored in the imaginary storage location of
C     the zeroth mode Fourier component.
C
C     Work array W could be reduced in size to 2*NZ if the inner loops
C     of scrambling/unscrambling were over NZ rather than NPH.
C     This would reduce vector performance on suitable machines, but
C     might be preferable on scalar ones.
C     ------------------------------------------------------------------
      IMPLICIT NONE
C
      INTEGER          NP, NPH, NZ, NZH, NZHM, IP, IQ, IR, ISIGN
      INTEGER          I, J, II, JE, JO, MJ
      DOUBLE PRECISION V    (NP,0:NZ-1)
      DOUBLE PRECISION W    (NP/2,2,0:NZ-1)
      DOUBLE PRECISION TRIG (*)
      LOGICAL          NZODD
C 
      NPH  = NP / 2
      NZH  = NZ / 2
      NZHM = NZH - 1
      IF (NZH+NZH .EQ. NZ) THEN
         NZODD = .FALSE.
      ELSE
         NZODD = .TRUE.
      ENDIF
C         
      IF (ISIGN .EQ. +1) THEN
C
C     -- Forward transform
C
         CALL DGPFA (V,V(1+NPH,0),TRIG,NP,1,NZ,IP,IQ,IR,NPH,-1)
         CALL DCOPY (NP*NZ, V, 1, W, 1)
C
         DO 40 J = 1, NZHM
            JE = J  + J
            JO = JE + 1
            MJ = NZ - J
            DO 50 I = 1, NPH
               II = I + NPH
               V (I,  JE) = 0.5D0 * (W (I, 1,  J) + W (I, 1, MJ))
               V (I,  JO) = 0.5D0 * (W (I, 2,  J) - W (I, 2, MJ))
               V (II, JE) = 0.5D0 * (W (I, 2,  J) + W (I, 2, MJ))
               V (II, JO) = 0.5D0 * (W (I, 1, MJ) - W (I, 1,  J))
 50         CONTINUE
 40      CONTINUE
         DO 60 I = 1, NPH
            II = I + NPH
            V (I,  0) = W (I, 1,   0)
            V (I,  1) = W (I, 1, NZH)
            V (II, 0) = W (I, 2,   0)
            V (II, 1) = W (I, 2, NZH)
 60      CONTINUE
C     
      ELSE
C
C     -- Inverse transform
C
         DO 100 J = 1, NZHM
            JE = J  + J
            JO = JE + 1
            MJ = NZ - J
            DO 110 I = 1, NPH
               II = I + NPH
               W (I, 1,  J) = V (I,  JE) - V (II, JO)
               W (I, 2,  J) = V (I,  JO) + V (II, JE)
               W (I, 1, MJ) = V (I,  JE) + V (II, JO)
               W (I, 2, MJ) = V (II, JE) - V (I,  JO)
 110        CONTINUE
 100     CONTINUE
         DO 120 I = 1, NPH
            II = I + NPH
            W (I, 1,   0) = V (I,  0)
            W (I, 2,   0) = V (II, 0)
            W (I, 1, NZH) = V (I,  1)
            W (I, 2, NZH) = V (II, 1)
 120     CONTINUE
         CALL DGPFA (W,W(1,2,0),TRIG,NP,1,NZ,IP,IQ,IR,NPH,+1)
         CALL DCOPY (NZ*NP, W, 1, V, 1)
C     
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE DGPFA (A,B,TRIG,INC,JUMP,N,IP,IQ,IR,LOT,ISIGN)
C     ------------------------------------------------------------------
C     Self-sorting in-place generalized prime factor (complex) FFT.
C
C     A     is first real      input/output vector
C     B     is first imaginary input/output vector
C     TRIG  is a table of twiddle factors, precalculated
C             by calling subroutine 'DSETPF'
C     INC   is the increment within each data vector
C     JUMP  is the increment between data vectors
C     N     is the length of the transforms:
C     IP    is the power of 2 in the prime factor decomposition of N
C     IQ    is the power of 3 in the prime factor decomposition of N
C     IR    is the power of 5 in the prime factor decomposition of N
C           -----------------------------------
C           N = (2**IP) * (3**IQ) * (5**IR)
C           -----------------------------------
C     LOT   is the number of transforms
C     ISIGN = +1 for 'forward' transform  --- NB the definition below
C           = -1 for 'inverse' transform
C     
C     Definition of transform
C     -----------------------
C
C     X(J) = SUM(K=0,...,N-1)(C(K)*EXP(ISIGN*2*I*J*K*PI/N))
C     ------------------------------------------------------------------
      IMPLICIT NONE
C
      DOUBLE PRECISION A(*), B(*), TRIG(*)
      INTEGER          INC, JUMP, N, IP, IQ, IR, LOT, ISIGN
      INTEGER          I
C
C     COMPUTE THE TRANSFORM
C     ---------------------
      I = 1
      IF (IP.GT.0) THEN
         CALL DPFA2F (A,B,TRIG(I),INC,JUMP,N,IP,LOT,ISIGN)
         I = I + 2 * (2**IP)
      ENDIF
      IF (IQ.GT.0) THEN
         CALL DPFA3F (A,B,TRIG(I),INC,JUMP,N,IQ,LOT,ISIGN)
         I = I + 2 * (3**IQ)
      ENDIF
      IF (IR.GT.0) THEN
         CALL DPFA5F (A,B,TRIG(I),INC,JUMP,N,IR,LOT,ISIGN)
      ENDIF
C     
      RETURN
      END


      SUBROUTINE DPFA2F(A,B,TRIG,INC,JUMP,N,MM,LOT,ISIGN)
C     ------------------------------------------------------------------
C     Radix-2 FFT.
C
C     ******************************************************************
C     *                                                                *
C     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90.    *
C     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE    *
C     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.         *
C     *                                                                *
C     ******************************************************************
C     ------------------------------------------------------------------
      IMPLICIT NONE
C
      DOUBLE PRECISION A(*), B(*), TRIG(*)
      INTEGER          INC, JUMP, N, MM, LOT, ISIGN
C
      DOUBLE PRECISION C1, C2, C3, CO1, CO2, CO3, CO4, CO5, CO6, CO7
      DOUBLE PRECISION S,  SS,     SI1, SI2, SI3, SI4, SI5, SI6, SI7
      DOUBLE PRECISION AJA, AJB, AJC, AJD, AJE, AJF, AJG, AJH
      DOUBLE PRECISION AJI, AJJ, AJK, AJL, AJM, AJN, AJO, AJP
      DOUBLE PRECISION BJA, BJB, BJC, BJD, BJE, BJF, BJG, BJH 
      DOUBLE PRECISION BJI, BJJ, BJK, BJL, BJM, BJN, BJO, BJP
      DOUBLE PRECISION T0, T1, T2, T3, U0, U1, U2, U3
C
      INTEGER          ISTART, IPASS, INQ, INK
      INTEGER          JA, JB, JC, JD, JE, JF, JG, JH
      INTEGER          JI, JJ, JK, JL, JM, JN, JO, JP
      INTEGER          J, JJJ, JSTEP, JSTEPL, JSTEPX
      INTEGER          K, KK
      INTEGER          L, LL, LA, LAINCL, LEFT, LVR
      INTEGER          M, MH, MU, M2, M8
      INTEGER          N2, NINC, NBLOX, NB, NVEX, NU
C
#if defined(_CRAY)
      PARAMETER        ( LVR = 64 )
#elif defined(_SX)
      PARAMETER        ( LVR = 256 )
#elif defined(__uxp__)
      PARAMETER        ( LVR = 2048 )
#else
      PARAMETER        ( LVR = 1024 )
#endif
C
      N2     = 2**MM
      INQ    = N/N2
      JSTEPX = (N2-N) * INC
      NINC   = N * INC
      INK    = INC * INQ
C
      M2 = 0
      M8 = 0
      IF (MOD(MM,2).EQ.0) THEN
         M = MM/2
      ELSE IF (MOD(MM,4).EQ.1) THEN
         M  = (MM-1)/2
         M2 = 1
      ELSE IF (MOD(MM,4).EQ.3) THEN
         M  = (MM-3)/2
         M8 = 1
      ENDIF
      MH = (M+1)/2
C
      NBLOX  = 1 + (LOT-1)/LVR
      LEFT   = LOT
      S      = DFLOAT(ISIGN)
      ISTART = 1
C
C     Loop on blocks of LVR transforms
C     --------------------------------
      DO 500 NB = 1 , NBLOX
C
         IF (LEFT.LE.LVR) THEN
            NVEX = LEFT
         ELSE IF (LEFT.LT.(2*LVR)) THEN
            NVEX = LEFT/2
            NVEX = NVEX + MOD(NVEX,2)
         ELSE
            NVEX = LVR
         ENDIF
         LEFT = LEFT - NVEX
C     
         LA = 1
C     
C     Loop on type I radix-4 passes
C     -----------------------------
         MU = MOD(INQ,4)
         IF (ISIGN.EQ.-1) MU = 4 - MU
         SS = 1.0
         IF (MU.EQ.3) SS = -1.0
C
         IF (MH.EQ.0) GO TO 200
C
         DO 160 IPASS = 1 , MH
            JSTEP  = (N*INC) / (4*LA)
            JSTEPL = JSTEP - NINC
C
C     K = 0 loop (no twiddle factors)
C     -------------------------------
            DO 120 JJJ = 0 , (N-1)*INC , 4*JSTEP
               JA = ISTART + JJJ
C
C     "Transverse" loop
C     -----------------
               DO 115 NU = 1 , INQ
                  JB = JA + JSTEPL
                  IF (JB.LT.ISTART) JB = JB + NINC
                  JC = JB + JSTEPL
                  IF (JC.LT.ISTART) JC = JC + NINC
                  JD = JC + JSTEPL
                  IF (JD.LT.ISTART) JD = JD + NINC
                  J = 0
C     
C     Loop across transforms
C     ----------------------
#if   defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
                  DO 110 L = 1 , NVEX
                     AJA = A(JA+J)
                     AJC = A(JC+J)
                     T0  = AJA + AJC
                     T2  = AJA - AJC
                     AJB = A(JB+J)
                     AJD = A(JD+J)
                     T1  = AJB + AJD
                     T3  = SS * ( AJB - AJD )
                     BJA = B(JA+J)
                     BJC = B(JC+J)
                     U0  = BJA + BJC
                     U2  = BJA - BJC
                     BJB = B(JB+J)
                     BJD = B(JD+J)
                     U1  = BJB + BJD
                     U3  = SS * ( BJB - BJD )
                     A(JA+J) = T0 + T1
                     A(JC+J) = T0 - T1
                     B(JA+J) = U0 + U1
                     B(JC+J) = U0 - U1
                     A(JB+J) = T2 - U3
                     A(JD+J) = T2 + U3
                     B(JB+J) = U2 + T3
                     B(JD+J) = U2 - T3
                     J = J + JUMP
 110              CONTINUE
                  JA = JA + JSTEPX
                  IF (JA.LT.ISTART) JA = JA + NINC
 115           CONTINUE
 120        CONTINUE
C
C     Finished if N2 = 4
C     ------------------
            IF (N2.EQ.4) GO TO 490
            KK = 2 * LA
C     
C     Loop on nonzero K
C     -----------------
            DO 150 K = INK , JSTEP-INK , INK
               CO1 = TRIG(KK+1)
               SI1 = S*TRIG(KK+2)
               CO2 = TRIG(2*KK+1)
               SI2 = S*TRIG(2*KK+2)
               CO3 = TRIG(3*KK+1)
               SI3 = S*TRIG(3*KK+2)
C
C     Loop along transform
C     --------------------
               DO 140 JJJ = K , (N-1)*INC , 4*JSTEP
                  JA = ISTART + JJJ
C
C     "Transverse" loop
C     -----------------
                  DO 135 NU = 1 , INQ
                     JB = JA + JSTEPL
                     IF (JB.LT.ISTART) JB = JB + NINC
                     JC = JB + JSTEPL
                     IF (JC.LT.ISTART) JC = JC + NINC
                     JD = JC + JSTEPL
                     IF (JD.LT.ISTART) JD = JD + NINC
                     J = 0
C
C     Loop across transforms
C     ----------------------
#if defined (_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
                     DO 130 L = 1 , NVEX
                        AJA = A(JA+J)
                        AJC = A(JC+J)
                        T0  = AJA + AJC
                        T2  = AJA - AJC
                        AJB = A(JB+J)
                        AJD = A(JD+J)
                        T1  = AJB + AJD
                        T3  = SS * ( AJB - AJD )
                        BJA = B(JA+J)
                        BJC = B(JC+J)
                        U0  = BJA + BJC
                        U2  = BJA - BJC
                        BJB = B(JB+J)
                        BJD = B(JD+J)
                        U1  = BJB + BJD
                        U3  = SS * ( BJB - BJD )
                        A(JA+J) = T0 + T1
                        B(JA+J) = U0 + U1
                        A(JB+J) = CO1*(T2-U3) - SI1*(U2+T3)
                        B(JB+J) = SI1*(T2-U3) + CO1*(U2+T3)
                        A(JC+J) = CO2*(T0-T1) - SI2*(U0-U1)
                        B(JC+J) = SI2*(T0-T1) + CO2*(U0-U1)
                        A(JD+J) = CO3*(T2+U3) - SI3*(U2-T3)
                        B(JD+J) = SI3*(T2+U3) + CO3*(U2-T3)
                        J = J + JUMP
 130                 CONTINUE
C     -----(End of loop across transforms)
                     JA = JA + JSTEPX
                     IF (JA.LT.ISTART) JA = JA + NINC
 135              CONTINUE
 140           CONTINUE
C     -----(End of loop along transforms)
               KK = KK + 2*LA
 150        CONTINUE
C     -----(End of loop on nonzero K)
            LA = 4*LA
 160     CONTINUE
C     -----(End of loop on type I radix-4 passes)
C
C     Central radix-2 pass
C     --------------------
 200     CONTINUE
         IF (M2.EQ.0) GO TO 300
C
         JSTEP  = (N*INC) / (2*LA)
         JSTEPL = JSTEP - NINC
C
C     K=0 loop (no twiddle factors)
C     -----------------------------
         DO 220 JJJ = 0 , (N-1)*INC , 2*JSTEP
            JA = ISTART + JJJ
C
C     "Transverse" loop
C     -----------------
            DO 215 NU = 1 , INQ
               JB = JA + JSTEPL
               IF (JB.LT.ISTART) JB = JB + NINC
               J = 0
C
C     Loop across transforms
C     ----------------------
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
               DO 210 L = 1 , NVEX
                  AJA     = A(JA+J)
                  AJB     = A(JB+J)
                  T0      = AJA - AJB
                  A(JA+J) = AJA + AJB
                  A(JB+J) = T0
                  BJA     = B(JA+J)
                  BJB     = B(JB+J)
                  U0      = BJA - BJB
                  B(JA+J) = BJA + BJB
                  B(JB+J) = U0
                  J = J + JUMP
 210           CONTINUE
C     -----(End of loop across transforms)
               JA = JA + JSTEPX
               IF (JA.LT.ISTART) JA = JA + NINC
 215        CONTINUE
 220     CONTINUE
C
C     Finished if N2=2
C     ----------------
         IF (N2.EQ.2) GO TO 490
C
         KK = 2 * LA
C
C     Loop on nonzero K
C     -----------------
         DO 260 K = INK , JSTEP - INK , INK
            CO1 = TRIG(KK+1)
            SI1 = S*TRIG(KK+2)
C
C     Loop along transforms
C     ---------------------
            DO 250 JJJ = K , (N-1)*INC , 2*JSTEP
               JA = ISTART + JJJ
C     
C     "Transverse" loop
C     -----------------
               DO 245 NU = 1 , INQ
                  JB = JA + JSTEPL
                  IF (JB.LT.ISTART) JB = JB + NINC
                  J = 0
C
C     Loop across transforms
C     ----------------------
                  IF (KK.EQ.N2/2) THEN
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
                     DO 230 L = 1 , NVEX
                        AJA     = A(JA+J)
                        AJB     = A(JB+J)
                        T0      = SS * ( AJA - AJB )
                        A(JA+J) = AJA + AJB
                        BJB     = B(JB+J)
                        BJA     = B(JA+J)
                        A(JB+J) = SS * ( BJB - BJA )
                        B(JA+J) = BJA + BJB
                        B(JB+J) = T0
                        J = J + JUMP
 230                 CONTINUE
C
      ELSE
C
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
         DO 240 L = 1 , NVEX
            AJA     = A(JA+J)
            AJB     = A(JB+J)
            T0      = AJA - AJB
            A(JA+J) = AJA + AJB
            BJA     = B(JA+J)
            BJB     = B(JB+J)
            U0      = BJA - BJB
            B(JA+J) = BJA + BJB
            A(JB+J) = CO1*T0 - SI1*U0
            B(JB+J) = SI1*T0 + CO1*U0
            J = J + JUMP
 240     CONTINUE
C
      ENDIF
C
C     -----(End of loop across transforms)
      JA = JA + JSTEPX
      IF (JA.LT.ISTART) JA = JA + NINC
 245  CONTINUE
 250  CONTINUE
C     -----(End of loop along transforms)
      KK = KK + 2 * LA
 260  CONTINUE
C     -----(End of loop on nonzero K)
C     -----(End of radix-2 pass)
C
      LA = 2 * LA
      GO TO 400
C
C     Central radix-8 pass
C     --------------------
 300  CONTINUE
      IF (M8.EQ.0) GO TO 400
      JSTEP  = (N*INC) / (8*LA)
      JSTEPL = JSTEP - NINC
      MU     = MOD(INQ,8)
      IF (ISIGN.EQ.-1) MU = 8 - MU
      C1 = 1.0
      IF (MU.EQ.3.OR.MU.EQ.7) C1 = -1.0
      C2 = DSQRT(0.5D0)
      IF (MU.EQ.3.OR.MU.EQ.5) C2 = -C2
      C3 = C1 * C2
C
C     Stage 1
C     -------
      DO 320 K = 0 , JSTEP - INK , INK
         DO 315 JJJ = K , (N-1)*INC , 8*JSTEP
            JA = ISTART + JJJ
C
C     "Transverse" loop
C     -----------------
            DO 312 NU = 1 , INQ
               JB = JA + JSTEPL
               IF (JB.LT.ISTART) JB = JB + NINC
               JC = JB + JSTEPL
               IF (JC.LT.ISTART) JC = JC + NINC
               JD = JC + JSTEPL
               IF (JD.LT.ISTART) JD = JD + NINC
               JE = JD + JSTEPL
               IF (JE.LT.ISTART) JE = JE + NINC
               JF = JE + JSTEPL
               IF (JF.LT.ISTART) JF = JF + NINC
               JG = JF + JSTEPL
               IF (JG.LT.ISTART) JG = JG + NINC
               JH = JG + JSTEPL
               IF (JH.LT.ISTART) JH = JH + NINC
               J = 0
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
               DO 310 L = 1 , NVEX
                  AJA     = A(JA+J)
                  AJE     = A(JE+J)
                  T0      = AJA - AJE
                  A(JA+J) = AJA + AJE
                  AJC     = A(JC+J)
                  AJG     = A(JG+J)
                  T1      = C1 * ( AJC - AJG )
                  A(JE+J) = AJC + AJG
                  AJB     = A(JB+J)
                  AJF     = A(JF+J)
                  T2      = AJB - AJF
                  A(JC+J) = AJB + AJF
                  AJD     = A(JD+J)
                  AJH     = A(JH+J)
                  T3      = AJD - AJH
                  A(JG+J) = AJD + AJH
                  A(JB+J) = T0
                  A(JF+J) = T1
                  A(JD+J) = C2 * ( T2 - T3 )
                  A(JH+J) = C3 * ( T2 + T3 )
                  BJA     = B(JA+J)
                  BJE     = B(JE+J)
                  U0      = BJA - BJE
                  B(JA+J) = BJA + BJE
                  BJC     = B(JC+J)
                  BJG     = B(JG+J)
                  U1      = C1 * ( BJC - BJG )
                  B(JE+J) = BJC + BJG
                  BJB     = B(JB+J)
                  BJF     = B(JF+J)
                  U2      = BJB - BJF
                  B(JC+J) = BJB + BJF
                  BJD     = B(JD+J)
                  BJH     = B(JH+J)
                  U3      = BJD - BJH
                  B(JG+J) = BJD + BJH
                  B(JB+J) = U0
                  B(JF+J) = U1
                  B(JD+J) = C2 * ( U2 - U3 )
                  B(JH+J) = C3 * ( U2 + U3 )
                  J = J + JUMP
 310           CONTINUE
               JA = JA + JSTEPX
               IF (JA.LT.ISTART) JA = JA + NINC
 312        CONTINUE
 315     CONTINUE
 320  CONTINUE
C
C     Stage 2
C     -------
C
C     K=0 (No twiddle factors)
C     ------------------------
      DO 330 JJJ = 0 , (N-1)*INC , 8*JSTEP
         JA = ISTART + JJJ
C
C     "Transverse" loop
C     -----------------
         DO 328 NU = 1 , INQ
            JB = JA + JSTEPL
            IF (JB.LT.ISTART) JB = JB + NINC
            JC = JB + JSTEPL
            IF (JC.LT.ISTART) JC = JC + NINC
            JD = JC + JSTEPL
            IF (JD.LT.ISTART) JD = JD + NINC
            JE = JD + JSTEPL
            IF (JE.LT.ISTART) JE = JE + NINC
            JF = JE + JSTEPL
            IF (JF.LT.ISTART) JF = JF + NINC
            JG = JF + JSTEPL
            IF (JG.LT.ISTART) JG = JG + NINC
            JH = JG + JSTEPL
            IF (JH.LT.ISTART) JH = JH + NINC
            J = 0
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
            DO 325 L = 1 , NVEX
               AJA     = A(JA+J)
               AJE     = A(JE+J)
               T0      = AJA + AJE
               T2      = AJA - AJE
               AJC     = A(JC+J)
               AJG     = A(JG+J)
               T1      = AJC + AJG
               T3      = C1 * ( AJC - AJG )
               BJA     = B(JA+J)
               BJE     = B(JE+J)
               U0      = BJA + BJE
               U2      = BJA - BJE
               BJC     = B(JC+J)
               BJG     = B(JG+J)
               U1      = BJC + BJG
               U3      = C1 * ( BJC - BJG )
               A(JA+J) = T0 + T1
               A(JE+J) = T0 - T1
               B(JA+J) = U0 + U1
               B(JE+J) = U0 - U1
               A(JC+J) = T2 - U3
               A(JG+J) = T2 + U3
               B(JC+J) = U2 + T3
               B(JG+J) = U2 - T3
               AJB     = A(JB+J)
               AJD     = A(JD+J)
               T0      = AJB + AJD
               T2      = AJB - AJD
               AJF     = A(JF+J)
               AJH     = A(JH+J)
               T1      = AJF - AJH
               T3      = AJF + AJH
               BJB     = B(JB+J)
               BJD     = B(JD+J)
               U0      = BJB + BJD
               U2      = BJB - BJD
               BJF     = B(JF+J)
               BJH     = B(JH+J)
               U1      = BJF - BJH
               U3      = BJF + BJH
               A(JB+J) = T0 - U3
               A(JH+J) = T0 + U3
               B(JB+J) = U0 + T3
               B(JH+J) = U0 - T3
               A(JD+J) = T2 + U1
               A(JF+J) = T2 - U1
               B(JD+J) = U2 - T1
               B(JF+J) = U2 + T1
               J = J + JUMP
 325        CONTINUE
            JA = JA + JSTEPX
            IF (JA.LT.ISTART) JA = JA + NINC
 328     CONTINUE
 330  CONTINUE
C
      IF (N2.EQ.8) GO TO 490
C
C     Loop on nonzero K
C     -----------------
      KK = 2 * LA
C
      DO 350 K = INK , JSTEP - INK , INK
C     
         CO1 = TRIG(KK+1)
         SI1 = S * TRIG(KK+2)
         CO2 = TRIG(2*KK+1)
         SI2 = S * TRIG(2*KK+2)
         CO3 = TRIG(3*KK+1)
         SI3 = S * TRIG(3*KK+2)
         CO4 = TRIG(4*KK+1)
         SI4 = S * TRIG(4*KK+2)
         CO5 = TRIG(5*KK+1)
         SI5 = S * TRIG(5*KK+2)
         CO6 = TRIG(6*KK+1)
         SI6 = S * TRIG(6*KK+2)
         CO7 = TRIG(7*KK+1)
         SI7 = S * TRIG(7*KK+2)
C
         DO 345 JJJ = K , (N-1)*INC , 8*JSTEP
            JA = ISTART + JJJ
C
C     "Transverse" loop
C     -----------------
            DO 342 NU = 1 , INQ
               JB = JA + JSTEPL
               IF (JB.LT.ISTART) JB = JB + NINC
               JC = JB + JSTEPL
               IF (JC.LT.ISTART) JC = JC + NINC
               JD = JC + JSTEPL
               IF (JD.LT.ISTART) JD = JD + NINC
               JE = JD + JSTEPL
               IF (JE.LT.ISTART) JE = JE + NINC
               JF = JE + JSTEPL
               IF (JF.LT.ISTART) JF = JF + NINC
               JG = JF + JSTEPL
               IF (JG.LT.ISTART) JG = JG + NINC
               JH = JG + JSTEPL
               IF (JH.LT.ISTART) JH = JH + NINC
               J = 0
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
               DO 340 L = 1 , NVEX
                  AJA     = A(JA+J)
                  AJE     = A(JE+J)
                  T0      = AJA + AJE
                  T2      = AJA - AJE
                  AJC     = A(JC+J)
                  AJG     = A(JG+J)
                  T1      = AJC + AJG
                  T3      = C1 * ( AJC - AJG )
                  BJA     = B(JA+J)
                  BJE     = B(JE+J)
                  U0      = BJA + BJE
                  U2      = BJA - BJE
                  BJC     = B(JC+J)
                  BJG     = B(JG+J)
                  U1      = BJC + BJG
                  U3      = C1 * ( BJC - BJG )
                  A(JA+J) = T0 + T1
                  B(JA+J) = U0 + U1
                  A(JE+J) = CO4*(T0-T1) - SI4*(U0-U1)
                  B(JE+J) = SI4*(T0-T1) + CO4*(U0-U1)
                  A(JC+J) = CO2*(T2-U3) - SI2*(U2+T3)
                  B(JC+J) = SI2*(T2-U3) + CO2*(U2+T3)
                  A(JG+J) = CO6*(T2+U3) - SI6*(U2-T3)
                  B(JG+J) = SI6*(T2+U3) + CO6*(U2-T3)
                  AJB     = A(JB+J)
                  AJD     = A(JD+J)
                  T0      = AJB + AJD
                  T2      = AJB - AJD
                  AJF     = A(JF+J)
                  AJH     = A(JH+J)
                  T1      = AJF - AJH
                  T3      = AJF + AJH
                  BJB     = B(JB+J)
                  BJD     = B(JD+J)
                  U0      = BJB + BJD
                  U2      = BJB - BJD
                  BJF     = B(JF+J)
                  BJH     = B(JH+J)
                  U1      = BJF - BJH
                  U3      = BJF + BJH
                  A(JB+J) = CO1*(T0-U3) - SI1*(U0+T3)
                  B(JB+J) = SI1*(T0-U3) + CO1*(U0+T3)
                  A(JH+J) = CO7*(T0+U3) - SI7*(U0-T3)
                  B(JH+J) = SI7*(T0+U3) + CO7*(U0-T3)
                  A(JD+J) = CO3*(T2+U1) - SI3*(U2-T1)
                  B(JD+J) = SI3*(T2+U1) + CO3*(U2-T1)
                  A(JF+J) = CO5*(T2-U1) - SI5*(U2+T1)
                  B(JF+J) = SI5*(T2-U1) + CO5*(U2+T1)
                  J = J + JUMP
 340           CONTINUE
               JA = JA + JSTEPX
               IF (JA.LT.ISTART) JA = JA + NINC
 342        CONTINUE
 345     CONTINUE
         KK = KK + 2 * LA
 350  CONTINUE
C
      LA = 8 * LA
C
C     Loop on type II radix-4 passes
C     ------------------------------
 400  CONTINUE
      MU = MOD(INQ,4)
      IF (ISIGN.EQ.-1) MU = 4 - MU
      SS = 1.0
      IF (MU.EQ.3) SS = -1.0
C
      DO 480 IPASS = MH+1 , M
         JSTEP = (N*INC) / (4*LA)
         JSTEPL = JSTEP - NINC
         LAINCL = LA * INK - NINC
C
C     K=0 loop (no twiddle factors)
C     -----------------------------
         DO 430 LL = 0 , (LA-1)*INK , 4*JSTEP
C
            DO 420 JJJ = LL , (N-1)*INC , 4*LA*INK
               JA = ISTART + JJJ
C
C     "Transverse" loop
C     -----------------
               DO 415 NU = 1 , INQ
                  JB = JA + JSTEPL
                  IF (JB.LT.ISTART) JB = JB + NINC
                  JC = JB + JSTEPL
                  IF (JC.LT.ISTART) JC = JC + NINC
                  JD = JC + JSTEPL
                  IF (JD.LT.ISTART) JD = JD + NINC
                  JE = JA + LAINCL
                  IF (JE.LT.ISTART) JE = JE + NINC
                  JF = JE + JSTEPL
                  IF (JF.LT.ISTART) JF = JF + NINC
                  JG = JF + JSTEPL
                  IF (JG.LT.ISTART) JG = JG + NINC
                  JH = JG + JSTEPL
                  IF (JH.LT.ISTART) JH = JH + NINC
                  JI = JE + LAINCL
                  IF (JI.LT.ISTART) JI = JI + NINC
                  JJ = JI + JSTEPL
                  IF (JJ.LT.ISTART) JJ = JJ + NINC
                  JK = JJ + JSTEPL
                  IF (JK.LT.ISTART) JK = JK + NINC
                  JL = JK + JSTEPL
                  IF (JL.LT.ISTART) JL = JL + NINC
                  JM = JI + LAINCL
                  IF (JM.LT.ISTART) JM = JM + NINC
                  JN = JM + JSTEPL
                  IF (JN.LT.ISTART) JN = JN + NINC
                  JO = JN + JSTEPL
                  IF (JO.LT.ISTART) JO = JO + NINC
                  JP = JO + JSTEPL
                  IF (JP.LT.ISTART) JP = JP + NINC
                  J = 0
C
C     Loop across transforms
C     ----------------------
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
                  DO 410 L = 1 , NVEX
                     AJA     = A(JA+J)
                     AJC     = A(JC+J)
                     T0      = AJA + AJC
                     T2      = AJA - AJC
                     AJB     = A(JB+J)
                     AJD     = A(JD+J)
                     T1      = AJB + AJD
                     T3      = SS * ( AJB - AJD )
                     AJI     = A(JI+J)
                     AJC     = AJI
                     BJA     = B(JA+J)
                     BJC     = B(JC+J)
                     U0      = BJA + BJC
                     U2      = BJA - BJC
                     BJB     = B(JB+J)
                     BJD     = B(JD+J)
                     U1      = BJB + BJD
                     U3      = SS * ( BJB - BJD )
                     AJE     = A(JE+J)
                     AJB     =  AJE
                     A(JA+J) = T0 + T1
                     A(JI+J) = T0 - T1
                     B(JA+J) = U0 + U1
                     BJC     =  U0 - U1
                     BJM     = B(JM+J)
                     BJD     =  BJM
                     A(JE+J) = T2 - U3
                     AJD     =  T2 + U3
                     BJB     =  U2 + T3
                     B(JM+J) = U2 - T3
C     ----------------------
                     AJG     = A(JG+J)
                     T0      = AJB + AJG
                     T2      = AJB - AJG
                     AJF     = A(JF+J)
                     AJH     = A(JH+J)
                     T1      = AJF + AJH
                     T3      = SS * ( AJF - AJH )
                     AJJ     = A(JJ+J)
                     AJG     =  AJJ
                     BJE     = B(JE+J)
                     BJG     = B(JG+J)
                     U0      = BJE + BJG
                     U2      = BJE - BJG
                     BJF     = B(JF+J)
                     BJH     = B(JH+J)
                     U1      = BJF + BJH
                     U3      = SS * ( BJF - BJH )
                     B(JE+J) = BJB
                     A(JB+J) = T0 + T1
                     A(JJ+J) = T0 - T1
                     BJJ     = B(JJ+J)
                     BJG     = BJJ
                     B(JB+J) = U0 + U1
                     B(JJ+J) = U0 - U1
                     A(JF+J) = T2 - U3
                     AJH     = T2 + U3
                     B(JF+J) = U2 + T3
                     BJH     =  U2 - T3
C     ----------------------
                     AJK     = A(JK+J)
                     T0      = AJC + AJK
                     T2      = AJC - AJK
                     AJL     = A(JL+J)
                     T1      = AJG + AJL
                     T3      = SS * ( AJG - AJL )
                     BJI     = B(JI+J)
                     BJK     = B(JK+J)
                     U0      = BJI + BJK
                     U2      = BJI - BJK
                     AJO     = A(JO+J)
                     AJL     =  AJO
                     BJL     = B(JL+J)
                     U1      = BJG + BJL
                     U3      = SS * ( BJG - BJL )
                     B(JI+J) = BJC
                     A(JC+J) = T0 + T1
                     A(JK+J) = T0 - T1
                     BJO     = B(JO+J)
                     BJL     =  BJO
                     B(JC+J) = U0 + U1
                     B(JK+J) = U0 - U1
                     A(JG+J) = T2 - U3
                     A(JO+J) = T2 + U3
                     B(JG+J) = U2 + T3
                     B(JO+J) = U2 - T3
C     ----------------------
                     AJM     = A(JM+J)
                     T0      = AJM + AJL
                     T2      = AJM - AJL
                     AJN     = A(JN+J)
                     AJP     = A(JP+J)
                     T1      = AJN + AJP
                     T3      = SS * ( AJN - AJP )
                     A(JM+J) = AJD
                     U0      = BJD + BJL
                     U2      = BJD - BJL
                     BJN     = B(JN+J)
                     BJP     = B(JP+J)
                     U1      = BJN + BJP
                     U3      = SS * ( BJN - BJP )
                     A(JN+J) = AJH
                     A(JD+J) = T0 + T1
                     A(JL+J) = T0 - T1
                     B(JD+J) = U0 + U1
                     B(JL+J) = U0 - U1
                     B(JN+J) = BJH
                     A(JH+J) = T2 - U3
                     A(JP+J) = T2 + U3
                     B(JH+J) = U2 + T3
                     B(JP+J) = U2 - T3
                     J = J + JUMP
 410              CONTINUE
C     -----(End of loop across transforms)
                  JA = JA + JSTEPX
                  IF (JA.LT.ISTART) JA = JA + NINC
 415           CONTINUE
 420        CONTINUE
 430     CONTINUE
C     -----(End of double loop for K=0 )
C
C     Finished if last pass
C     ---------------------
         IF (IPASS.EQ.M) GO TO 490
C
         KK = 2*LA
C
C     Loop on nonzero K
C     -----------------
         DO 470 K = INK , JSTEP-INK , INK
            CO1 = TRIG(KK+1)
            SI1 = S*TRIG(KK+2)
            CO2 = TRIG(2*KK+1)
            SI2 = S*TRIG(2*KK+2)
            CO3 = TRIG(3*KK+1)
            SI3 = S*TRIG(3*KK+2)
C
C     Double loop along first transform in block
C  ------------------------------------------
            DO 460 LL = K , (LA-1)*INK , 4*JSTEP
C
               DO 450 JJJ = LL , (N-1)*INC , 4*LA*INK
                  JA = ISTART + JJJ
C
C     "Transverse" loop
C     -----------------
                  DO 445 NU = 1 , INQ
                     JB = JA + JSTEPL
                     IF (JB.LT.ISTART) JB = JB + NINC
                     JC = JB + JSTEPL
                     IF (JC.LT.ISTART) JC = JC + NINC
                     JD = JC + JSTEPL
                     IF (JD.LT.ISTART) JD = JD + NINC
                     JE = JA + LAINCL
                     IF (JE.LT.ISTART) JE = JE + NINC
                     JF = JE + JSTEPL
                     IF (JF.LT.ISTART) JF = JF + NINC
                     JG = JF + JSTEPL
                     IF (JG.LT.ISTART) JG = JG + NINC
                     JH = JG + JSTEPL
                     IF (JH.LT.ISTART) JH = JH + NINC
                     JI = JE + LAINCL
                     IF (JI.LT.ISTART) JI = JI + NINC
                     JJ = JI + JSTEPL
                     IF (JJ.LT.ISTART) JJ = JJ + NINC
                     JK = JJ + JSTEPL
                     IF (JK.LT.ISTART) JK = JK + NINC
                     JL = JK + JSTEPL
                     IF (JL.LT.ISTART) JL = JL + NINC
                     JM = JI + LAINCL
                     IF (JM.LT.ISTART) JM = JM + NINC
                     JN = JM + JSTEPL
                     IF (JN.LT.ISTART) JN = JN + NINC
                     JO = JN + JSTEPL
                     IF (JO.LT.ISTART) JO = JO + NINC
                     JP = JO + JSTEPL
                     IF (JP.LT.ISTART) JP = JP + NINC
                     J = 0
C
C     Loop across transforms
C     ----------------------
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
                     DO 440 L = 1 , NVEX
                        AJA     = A(JA+J)
                        AJC     = A(JC+J)
                        T0      = AJA + AJC
                        T2      = AJA - AJC
                        AJB     = A(JB+J)
                        AJD     = A(JD+J)
                        T1      = AJB + AJD
                        T3      = SS * ( AJB - AJD )
                        AJI     = A(JI+J)
                        AJC     = AJI
                        BJA     = B(JA+J)
                        BJC     = B(JC+J)
                        U0      = BJA + BJC
                        U2      = BJA - BJC
                        BJB     = B(JB+J)
                        BJD     = B(JD+J)
                        U1      = BJB + BJD
                        U3      = SS * ( BJB - BJD )
                        AJE     = A(JE+J)
                        AJB     = AJE
                        A(JA+J) = T0 + T1
                        B(JA+J) = U0 + U1
                        A(JE+J) = CO1*(T2-U3) - SI1*(U2+T3)
                        BJB     = SI1*(T2-U3) + CO1*(U2+T3)
                        BJM     = B(JM+J)
                        BJD     = BJM
                        A(JI+J) = CO2*(T0-T1) - SI2*(U0-U1)
                        BJC     = SI2*(T0-T1) + CO2*(U0-U1)
                        AJD     = CO3*(T2+U3) - SI3*(U2-T3)
                        B(JM+J) = SI3*(T2+U3) + CO3*(U2-T3)
C     ----------------------------------------
                        AJG     = A(JG+J)
                        T0      = AJB + AJG
                        T2      = AJB - AJG
                        AJF     = A(JF+J)
                        AJH     = A(JH+J)
                        T1      = AJF + AJH
                        T3      = SS * ( AJF - AJH )
                        AJJ     = A(JJ+J)
                        AJG     = AJJ
                        BJE     = B(JE+J)
                        BJG     = B(JG+J)
                        U0      = BJE + BJG
                        U2      = BJE - BJG
                        BJF     = B(JF+J)
                        BJH     = B(JH+J)
                        U1      = BJF + BJH
                        U3      = SS * ( BJF - BJH )
                        B(JE+J) = BJB
                        A(JB+J) = T0 + T1
                        B(JB+J) = U0 + U1
                        BJJ     = B(JJ+J)
                        BJG     = BJJ
                        A(JF+J) = CO1*(T2-U3) - SI1*(U2+T3)
                        B(JF+J) = SI1*(T2-U3) + CO1*(U2+T3)
                        A(JJ+J) = CO2*(T0-T1) - SI2*(U0-U1)
                        B(JJ+J) = SI2*(T0-T1) + CO2*(U0-U1)
                        AJH     = CO3*(T2+U3) - SI3*(U2-T3)
                        BJH     = SI3*(T2+U3) + CO3*(U2-T3)
C     ----------------------------------------
                        AJK     = A(JK+J)
                        T0      = AJC + AJK
                        T2      = AJC - AJK
                        AJL     = A(JL+J)
                        T1      = AJG + AJL
                        T3      = SS * ( AJG - AJL )
                        BJI     = B(JI+J)
                        BJK     = B(JK+J)
                        U0      = BJI + BJK
                        U2      = BJI - BJK
                        AJO     = A(JO+J)
                        AJL     = AJO
                        BJL     = B(JL+J)
                        U1      = BJG + BJL
                        U3      = SS * ( BJG - BJL )
                        B(JI+J) = BJC
                        A(JC+J) = T0 + T1
                        B(JC+J) = U0 + U1
                        BJO     = B(JO+J)
                        BJL     = BJO
                        A(JG+J) = CO1*(T2-U3) - SI1*(U2+T3)
                        B(JG+J) = SI1*(T2-U3) + CO1*(U2+T3)
                        A(JK+J) = CO2*(T0-T1) - SI2*(U0-U1)
                        B(JK+J) = SI2*(T0-T1) + CO2*(U0-U1)
                        A(JO+J) = CO3*(T2+U3) - SI3*(U2-T3)
                        B(JO+J) = SI3*(T2+U3) + CO3*(U2-T3)
C     ----------------------------------------
                        AJM     = A(JM+J)
                        T0      = AJM + AJL
                        T2      = AJM - AJL
                        AJN     = A(JN+J)
                        AJP     = A(JP+J)
                        T1      = AJN + AJP
                        T3      = SS * ( AJN - AJP )
                        A(JM+J) = AJD
                        U0      = BJD + BJL
                        U2      = BJD - BJL
                        A(JN+J) = AJH
                        BJN     = B(JN+J)
                        BJP     = B(JP+J)
                        U1      = BJN + BJP
                        U3      = SS * ( BJN - BJP )
                        B(JN+J) = BJH
                        A(JD+J) = T0 + T1
                        B(JD+J) = U0 + U1
                        A(JH+J) = CO1*(T2-U3) - SI1*(U2+T3)
                        B(JH+J) = SI1*(T2-U3) + CO1*(U2+T3)
                        A(JL+J) = CO2*(T0-T1) - SI2*(U0-U1)
                        B(JL+J) = SI2*(T0-T1) + CO2*(U0-U1)
                        A(JP+J) = CO3*(T2+U3) - SI3*(U2-T3)
                        B(JP+J) = SI3*(T2+U3) + CO3*(U2-T3)
                        J = J + JUMP
 440                 CONTINUE
C     -----(End of loop across transforms)
                     JA = JA + JSTEPX
                     IF (JA.LT.ISTART) JA = JA + NINC
 445              CONTINUE
 450           CONTINUE
 460        CONTINUE
C     -----(End of double loop for this K)
            KK = KK + 2*LA
 470     CONTINUE
C     -----(End of loop over values of K)
         LA = 4*LA
 480  CONTINUE
C     -----(End of loop on type II radix-4 passes)
C     -----(NVEX transforms completed)
 490  CONTINUE
      ISTART = ISTART + NVEX * JUMP
 500  CONTINUE
C     -----(End of loop on blocks of transforms)
C
      RETURN
      END


      SUBROUTINE DPFA3F (A,B,TRIG,INC,JUMP,N,MM,LOT,ISIGN)
C     ------------------------------------------------------------------
C     Radix-3 FFT.
C
C     ******************************************************************
C     *                                                                *
C     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90.    *
C     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE    *
C     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.         *
C     *                                                                *
C     ******************************************************************
C     ------------------------------------------------------------------
      IMPLICIT NONE
C
      DOUBLE PRECISION A(*), B(*), TRIG(*)
      INTEGER          INC, JUMP, N, MM, LOT, ISIGN
C
      DOUBLE PRECISION AJA, AJB, AJC, AJD, AJE, AJF, AJG, AJH, AJI
      DOUBLE PRECISION BJA, BJB, BJC, BJD, BJE, BJF, BJG, BJH, BJI
      DOUBLE PRECISION S, C1, SIN60, CO1, CO2, SI1, SI2
      DOUBLE PRECISION T1, T2, T3, U1, U2, U3
C
      INTEGER          INK, INQ, ISTART, IPASS
      INTEGER          J, JJJ, JSTEP, JSTEPL, JSTEPX
      INTEGER          JA, JB, JC, JD, JE, JF, JG, JH, JI
      INTEGER          K, KK
      INTEGER          L, LA, LL, LEFT, LVR, LAINCL
      INTEGER          M, MU, MH
      INTEGER          N3, NB, NU, NBLOX, NVEX, NINC
C
      PARAMETER        ( SIN60 = 0.866025403784439D0 )
#if defined(_CRAY)
      PARAMETER        ( LVR = 64 )
#elif defined(_SX)
      PARAMETER        ( LVR = 256 )
#elif defined(__uxp__)
      PARAMETER        ( LVR = 2048 )
#else
      PARAMETER        ( LVR = 1024 )
#endif
C
      N3     = 3**MM
      INQ    = N/N3
      JSTEPX = (N3-N) * INC
      NINC   = N * INC
      INK    = INC * INQ
      MU     = MOD(INQ,3)
      IF (ISIGN.EQ.-1) MU = 3-MU
      M      = MM
      MH     = (M+1)/2
      S      = DFLOAT(ISIGN)
      C1     = SIN60
      IF (MU.EQ.2) C1 = -C1
C
      NBLOX  = 1 + (LOT-1)/LVR
      LEFT   = LOT
      S      = DFLOAT(ISIGN)
      ISTART = 1
C
C     Loop on blocks of LVR transforms
C     --------------------------------
      DO 500 NB = 1 , NBLOX
C
         IF (LEFT.LE.LVR) THEN
            NVEX = LEFT
         ELSE IF (LEFT.LT.(2*LVR)) THEN
            NVEX = LEFT/2
            NVEX = NVEX + MOD(NVEX,2)
         ELSE
            NVEX = LVR
         ENDIF
         LEFT = LEFT - NVEX
C     
         LA = 1
C     
C     Loop on type I radix-3 passes
C     -----------------------------
         DO 160 IPASS = 1 , MH
            JSTEP = (N*INC) / (3*LA)
            JSTEPL = JSTEP - NINC
C     
C     K = 0 loop (no twiddle factors)
C     -------------------------------
            DO 120 JJJ = 0 , (N-1)*INC , 3*JSTEP
               JA = ISTART + JJJ
C
C     "Transverse" loop
C     -----------------
               DO 115 NU = 1 , INQ
                  JB = JA + JSTEPL
                  IF (JB.LT.ISTART) JB = JB + NINC
                  JC = JB + JSTEPL
                  IF (JC.LT.ISTART) JC = JC + NINC
                  J = 0
C
C     Loop across transforms
C     ----------------------
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
                  DO 110 L = 1 , NVEX
                     AJB     = A(JB+J)
                     AJC     = A(JC+J)
                     T1      = AJB + AJC
                     AJA     = A(JA+J)
                     T2      = AJA - 0.5D0 * T1
                     T3      = C1 * ( AJB - AJC )
                     BJB     = B(JB+J)
                     BJC     = B(JC+J)
                     U1      = BJB + BJC
                     BJA     = B(JA+J)
                     U2      = BJA - 0.5D0 * U1
                     U3      = C1 * ( BJB - BJC )
                     A(JA+J) = AJA + T1
                     B(JA+J) = BJA + U1
                     A(JB+J) = T2 - U3
                     B(JB+J) = U2 + T3
                     A(JC+J) = T2 + U3
                     B(JC+J) = U2 - T3
                     J = J + JUMP
 110              CONTINUE
                  JA = JA + JSTEPX
                  IF (JA.LT.ISTART) JA = JA + NINC
 115           CONTINUE
 120        CONTINUE
C     
C     Finished if N3 = 3
C     ------------------
            IF (N3.EQ.3) GO TO 490
            KK = 2 * LA
C
C     Loop on nonzero K
C     -----------------
            DO 150 K = INK , JSTEP-INK , INK
               CO1 = TRIG(KK+1)
               SI1 = S*TRIG(KK+2)
               CO2 = TRIG(2*KK+1)
               SI2 = S*TRIG(2*KK+2)
C
C     Loop along transform
C     --------------------
               DO 140 JJJ = K , (N-1)*INC , 3*JSTEP
                  JA = ISTART + JJJ
C
C     "Transverse" loop
C     -----------------
                  DO 135 NU = 1 , INQ
                     JB = JA + JSTEPL
                     IF (JB.LT.ISTART) JB = JB + NINC
                     JC = JB + JSTEPL
                     IF (JC.LT.ISTART) JC = JC + NINC
                     J = 0
C
C     Loop across transforms
C     ----------------------
#if defined(_CRAY)
CDIR$ IVDEP,SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
                     DO 130 L = 1 , NVEX
                        AJB     = A(JB+J)
                        AJC     = A(JC+J)
                        T1      = AJB + AJC
                        AJA     = A(JA+J)
                        T2      = AJA - 0.5D0 * T1
                        T3      = C1 * ( AJB - AJC )
                        BJB     = B(JB+J)
                        BJC     = B(JC+J)
                        U1      = BJB + BJC
                        BJA     = B(JA+J)
                        U2      = BJA - 0.5D0 * U1
                        U3      = C1 * ( BJB - BJC )
                        A(JA+J) = AJA + T1
                        B(JA+J) = BJA + U1
                        A(JB+J) = CO1*(T2-U3) - SI1*(U2+T3)
                        B(JB+J) = SI1*(T2-U3) + CO1*(U2+T3)
                        A(JC+J) = CO2*(T2+U3) - SI2*(U2-T3)
                        B(JC+J) = SI2*(T2+U3) + CO2*(U2-T3)
                        J = J + JUMP
 130                 CONTINUE
C     -----(End of loop across transforms)
                     JA = JA + JSTEPX
                     IF (JA.LT.ISTART) JA = JA + NINC
 135              CONTINUE
 140           CONTINUE
C     -----(End of loop along transforms)
               KK = KK + 2*LA
 150        CONTINUE
C     -----(End of loop on nonzero K)
            LA = 3*LA
 160     CONTINUE
C     -----(End of loop on type I radix-3 passes)
C
C     Loop on type II radix-3 passes
C     ------------------------------
 400     CONTINUE
C
         DO 480 IPASS = MH+1 , M
            JSTEP  = (N*INC) / (3*LA)
            JSTEPL = JSTEP - NINC
            LAINCL = LA*INK - NINC
C
C     K=0 loop (no twiddle factors)
C     -----------------------------
            DO 430 LL = 0 , (LA-1)*INK , 3*JSTEP
C
               DO 420 JJJ = LL , (N-1)*INC , 3*LA*INK
                  JA = ISTART + JJJ
C
C     "Transverse" loop
C     -----------------
                  DO 415 NU = 1 , INQ
                     JB = JA + JSTEPL
                     IF (JB.LT.ISTART) JB = JB + NINC
                     JC = JB + JSTEPL
                     IF (JC.LT.ISTART) JC = JC + NINC
                     JD = JA + LAINCL
                     IF (JD.LT.ISTART) JD = JD + NINC
                     JE = JD + JSTEPL
                     IF (JE.LT.ISTART) JE = JE + NINC
                     JF = JE + JSTEPL
                     IF (JF.LT.ISTART) JF = JF + NINC
                     JG = JD + LAINCL
                     IF (JG.LT.ISTART) JG = JG + NINC
                     JH = JG + JSTEPL
                     IF (JH.LT.ISTART) JH = JH + NINC
                     JI = JH + JSTEPL
                     IF (JI.LT.ISTART) JI = JI + NINC
                     J = 0
C
C     Loop across transforms
C     ----------------------
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
                     DO 410 L = 1 , NVEX
                        AJB     = A(JB+J)
                        AJC     = A(JC+J)
                        T1      = AJB + AJC
                        AJA     = A(JA+J)
                        T2      = AJA - 0.5D0 * T1
                        T3      = C1 * ( AJB - AJC )
                        AJD     = A(JD+J)
                        AJB     = AJD
                        BJB     = B(JB+J)
                        BJC     = B(JC+J)
                        U1      = BJB + BJC
                        BJA     = B(JA+J)
                        U2      = BJA - 0.5D0 * U1
                        U3      = C1 * ( BJB - BJC )
                        BJD     = B(JD+J)
                        BJB     = BJD
                        A(JA+J) = AJA + T1
                        B(JA+J) = BJA + U1
                        A(JD+J) = T2 - U3
                        B(JD+J) = U2 + T3
                        AJC     =  T2 + U3
                        BJC     =  U2 - T3
C     ----------------------
                        AJE     = A(JE+J)
                        AJF     = A(JF+J)
                        T1      = AJE + AJF
                        T2      = AJB - 0.5D0 * T1
                        T3      = C1 * ( AJE - AJF )
                        AJH     = A(JH+J)
                        AJF     = AJH
                        BJE     = B(JE+J)
                        BJF     = B(JF+J)
                        U1      = BJE + BJF
                        U2      = BJB - 0.5D0 * U1
                        U3      = C1 * ( BJE - BJF )
                        BJH     = B(JH+J)
                        BJF     =  BJH
                        A(JB+J) = AJB + T1
                        B(JB+J) = BJB + U1
                        A(JE+J) = T2 - U3
                        B(JE+J) = U2 + T3
                        A(JH+J) = T2 + U3
                        B(JH+J) = U2 - T3
C     ----------------------
                        AJI     = A(JI+J)
                        T1      = AJF + AJI
                        AJG     = A(JG+J)
                        T2      = AJG - 0.5D0 * T1
                        T3      = C1 * ( AJF - AJI )
                        T1      = AJG + T1
                        A(JG+J) = AJC
                        BJI     = B(JI+J)
                        U1      = BJF + BJI
                        BJG     = B(JG+J)
                        U2      = BJG - 0.5D0 * U1
                        U3      = C1 * ( BJF - BJI )
                        U1      = BJG + U1
                        B(JG+J) = BJC
                        A(JC+J) = T1
                        B(JC+J) = U1
                        A(JF+J) = T2 - U3
                        B(JF+J) = U2 + T3
                        A(JI+J) = T2 + U3
                        B(JI+J) = U2 - T3
                        J = J + JUMP
 410                 CONTINUE
C     -----(End of loop across transforms)
                     JA = JA + JSTEPX
                     IF (JA.LT.ISTART) JA = JA + NINC
 415              CONTINUE
 420           CONTINUE
 430        CONTINUE
C     -----(End of double loop for K=0)
C
C     Finished if last pass
C     ---------------------
            IF (IPASS.EQ.M) GO TO 490
C
            KK = 2*LA
C
C     Loop on nonzero K
C     -----------------
            DO 470 K = INK , JSTEP-INK , INK
               CO1 = TRIG(KK+1)
               SI1 = S*TRIG(KK+2)
               CO2 = TRIG(2*KK+1)
               SI2 = S*TRIG(2*KK+2)
C
C     Double loop along first transform in block
C     ------------------------------------------
               DO 460 LL = K , (LA-1)*INK , 3*JSTEP
C
                  DO 450 JJJ = LL , (N-1)*INC , 3*LA*INK
                     JA = ISTART + JJJ
C
C     "Transverse" loop
C     -----------------
                     DO 445 NU = 1 , INQ
                        JB = JA + JSTEPL
                        IF (JB.LT.ISTART) JB = JB + NINC
                        JC = JB + JSTEPL
                        IF (JC.LT.ISTART) JC = JC + NINC
                        JD = JA + LAINCL
                        IF (JD.LT.ISTART) JD = JD + NINC
                        JE = JD + JSTEPL
                        IF (JE.LT.ISTART) JE = JE + NINC
                        JF = JE + JSTEPL
                        IF (JF.LT.ISTART) JF = JF + NINC
                        JG = JD + LAINCL
                        IF (JG.LT.ISTART) JG = JG + NINC
                        JH = JG + JSTEPL
                        IF (JH.LT.ISTART) JH = JH + NINC
                        JI = JH + JSTEPL
                        IF (JI.LT.ISTART) JI = JI + NINC
                        J = 0
C
C     Loop across transforms
C     ----------------------
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
                        DO 440 L = 1 , NVEX
                           AJB     = A(JB+J)
                           AJC     = A(JC+J)
                           T1      = AJB + AJC
                           AJA     = A(JA+J)
                           T2      = AJA - 0.5D0 * T1
                           T3      = C1 * ( AJB - AJC )
                           AJD     = A(JD+J)
                           AJB     = AJD
                           BJB     = B(JB+J)
                           BJC     = B(JC+J)
                           U1      = BJB + BJC
                           BJA     = B(JA+J)
                           U2      = BJA - 0.5D0 * U1
                           U3      = C1 * ( BJB - BJC )
                           BJD     = B(JD+J)
                           BJB     = BJD
                           A(JA+J) = AJA + T1
                           B(JA+J) = BJA + U1
                           A(JD+J) = CO1*(T2-U3) - SI1*(U2+T3)
                           B(JD+J) = SI1*(T2-U3) + CO1*(U2+T3)
                           AJC     = CO2*(T2+U3) - SI2*(U2-T3)
                           BJC     = SI2*(T2+U3) + CO2*(U2-T3)
C     ----------------------
                           AJE     = A(JE+J)
                           AJF     = A(JF+J)
                           T1      = AJE + AJF
                           T2      = AJB - 0.5D0 * T1
                           T3      = C1 * ( AJE - AJF )
                           AJH     = A(JH+J)
                           AJF     = AJH
                           BJE     = B(JE+J)
                           BJF     = B(JF+J)
                           U1      = BJE + BJF
                           U2      = BJB - 0.5D0 * U1
                           U3      = C1 * ( BJE - BJF )
                           BJH     = B(JH+J)
                           BJF     = BJH
                           A(JB+J) = AJB + T1
                           B(JB+J) = BJB + U1
                           A(JE+J) = CO1*(T2-U3) - SI1*(U2+T3)
                           B(JE+J) = SI1*(T2-U3) + CO1*(U2+T3)
                           A(JH+J) = CO2*(T2+U3) - SI2*(U2-T3)
                           B(JH+J) = SI2*(T2+U3) + CO2*(U2-T3)
C     ----------------------
                           AJI     = A(JI+J)
                           T1      = AJF + AJI
                           AJG     = A(JG+J)
                           T2      = AJG - 0.5D0 * T1
                           T3      = C1 * ( AJF - AJI )
                           T1      = AJG + T1
                           A(JG+J) = AJC
                           BJI     = B(JI+J)
                           U1      = BJF + BJI
                           BJG     = B(JG+J)
                           U2      = BJG - 0.5D0 * U1
                           U3      = C1 * ( BJF - BJI )
                           U1      = BJG + U1
                           B(JG+J) = BJC
                           A(JC+J) = T1
                           B(JC+J) = U1
                           A(JF+J) = CO1*(T2-U3) - SI1*(U2+T3)
                           B(JF+J) = SI1*(T2-U3) + CO1*(U2+T3)
                           A(JI+J) = CO2*(T2+U3) - SI2*(U2-T3)
                           B(JI+J) = SI2*(T2+U3) + CO2*(U2-T3)
                           J = J + JUMP
 440                    CONTINUE
C     -----(End of loop across transforms)
                        JA = JA + JSTEPX
                        IF (JA.LT.ISTART) JA = JA + NINC
 445                 CONTINUE
 450              CONTINUE
 460           CONTINUE
C     -----(End of double loop for this K)
               KK = KK + 2*LA
 470        CONTINUE
C     -----(End of loop over values of K)
            LA = 3*LA
 480     CONTINUE
C     -----(End of loop on type II radix-3 passes)
C     -----(NVEX transforms completed)
 490     CONTINUE
         ISTART = ISTART + NVEX * JUMP
 500  CONTINUE
C     -----(End of loop on blocks of transforms)
C
      RETURN
      END


      SUBROUTINE DPFA5F (A,B,TRIG,INC,JUMP,N,MM,LOT,ISIGN)
C     ------------------------------------------------------------------
C     Radix-5 FFT.
C
C     ******************************************************************
C     *                                                                *
C     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90.    *
C     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE    *
C     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.         *
C     *                                                                *
C     ******************************************************************
C     ------------------------------------------------------------------
      IMPLICIT NONE
C
      DOUBLE PRECISION A(*), B(*), TRIG(*)
      INTEGER          INC, JUMP, N, MM, LOT, ISIGN
C
      DOUBLE PRECISION AJA, AJB, AJC, AJD, AJE, AJF, AJG, AJH, AJI, AJJ
      DOUBLE PRECISION AJK, AJL, AJM, AJN, AJO, AJP, AJQ, AJR, AJS, AJT
      DOUBLE PRECISION AJU, AJW, AJV, AJX, AJY, AX
      DOUBLE PRECISION BJA, BJB, BJC, BJD, BJE, BJF, BJG, BJH, BJI, BJJ
      DOUBLE PRECISION BJK, BJL, BJM, BJN, BJO, BJP, BJQ, BJR, BJS, BJT
      DOUBLE PRECISION BJU, BJW, BJV, BJX, BJY, BX
      DOUBLE PRECISION S, C1, C2, C3
      DOUBLE PRECISION CO1, CO2, CO3, CO4
      DOUBLE PRECISION SI1, SI2, SI3, SI4
      DOUBLE PRECISION SIN36, SIN72, QRT5
      DOUBLE PRECISION T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11
      DOUBLE PRECISION U1, U2, U3, U4, U5, U6, U7, U8, U9, U10, U11
C
      INTEGER          INK, INQ, ISTART, IPASS
      INTEGER          JA, JB, JC, JD, JE, JF, JG, JH, JI, JJ, JK, JL
      INTEGER          JM, JN, JO, JP, JQ, JR, JS, JT, JU, JV, JW, JX
      INTEGER          JY
      INTEGER          J, JJJ, JSTEP, JSTEPL, JSTEPX
      INTEGER          K, KK
      INTEGER          L, LA, LL, LEFT, LAINCL, LVR
      INTEGER          M, MH, MU
      INTEGER          NB, NU, N5, NINC, NBLOX, NVEX
C     
      PARAMETER        ( SIN36 = 0.587785252292473D0 )
      PARAMETER        ( SIN72 = 0.951056516295154D0 )
      PARAMETER        ( QRT5  = 0.559016994374947D0 )
#if defined(_CRAY)
      PARAMETER        ( LVR = 64 )
#elif defined(_SX)
      PARAMETER        ( LVR = 256 )
#elif defined(__uxp__)
      PARAMETER        ( LVR = 2048 )
#else
      PARAMETER        ( LVR = 1024 )
#endif
C
      N5     = 5 ** MM
      INQ    = N / N5
      JSTEPX = (N5-N) * INC
      NINC   = N * INC
      INK    = INC * INQ
      MU     = MOD(INQ,5)
      IF (ISIGN.EQ.-1) MU = 5 - MU
C
      M      = MM
      MH     = (M+1)/2
      S      = DFLOAT(ISIGN)
      C1     = QRT5
      C2     = SIN72
      C3     = SIN36
      IF (MU.EQ.2.OR.MU.EQ.3) THEN
         C1 = -C1
         C2 = SIN36
         C3 = SIN72
      ENDIF
      IF (MU.EQ.3.OR.MU.EQ.4) C2 = -C2
      IF (MU.EQ.2.OR.MU.EQ.4) C3 = -C3
C
      NBLOX  = 1 + (LOT-1)/LVR
      LEFT   = LOT
      S      = DFLOAT(ISIGN)
      ISTART = 1
C
C     Loop on blocks of LVR transforms
C     --------------------------------
      DO 500 NB = 1 , NBLOX
C     
         IF (LEFT.LE.LVR) THEN
            NVEX = LEFT
         ELSE IF (LEFT.LT.(2*LVR)) THEN
            NVEX = LEFT/2
            NVEX = NVEX + MOD(NVEX,2)
         ELSE
            NVEX = LVR
         ENDIF
         LEFT = LEFT - NVEX
C     
         LA = 1
C
C     Loop on type I radix-5 passes
C     -----------------------------
         DO 160 IPASS = 1 , MH
            JSTEP = (N*INC) / (5*LA)
            JSTEPL = JSTEP - NINC
            KK = 0
C     
C     Loop on K
C     ---------
            DO 150 K = 0 , JSTEP-INK , INK
C     
               IF (K.GT.0) THEN
                  CO1 = TRIG(KK+1)
                  SI1 = S*TRIG(KK+2)
                  CO2 = TRIG(2*KK+1)
                  SI2 = S*TRIG(2*KK+2)
                  CO3 = TRIG(3*KK+1)
                  SI3 = S*TRIG(3*KK+2)
                  CO4 = TRIG(4*KK+1)
                  SI4 = S*TRIG(4*KK+2)
               ENDIF
C
C     Loop along transform
C     --------------------
               DO 140 JJJ = K , (N-1)*INC , 5*JSTEP
                  JA = ISTART + JJJ
C
C     "Transverse" loop
C     -----------------
                  DO 135 NU = 1 , INQ
                     JB = JA + JSTEPL
                     IF (JB.LT.ISTART) JB = JB + NINC
                     JC = JB + JSTEPL
                     IF (JC.LT.ISTART) JC = JC + NINC
                     JD = JC + JSTEPL
                     IF (JD.LT.ISTART) JD = JD + NINC
                     JE = JD + JSTEPL
                     IF (JE.LT.ISTART) JE = JE + NINC
                     J = 0
C
C     Loop across transforms
C     ----------------------
                     IF (K.EQ.0) THEN
C
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
                        DO 110 L = 1 , NVEX
                           AJB     = A(JB+J)
                           AJE     = A(JE+J)
                           T1      = AJB + AJE
                           AJC     = A(JC+J)
                           AJD     = A(JD+J)
                           T2      = AJC + AJD
                           T3      = AJB - AJE
                           T4      = AJC - AJD
                           T5      = T1 + T2
                           T6      = C1 * ( T1 - T2 )
                           AJA     = A(JA+J)
                           T7      = AJA - 0.25D0 * T5
                           A(JA+J) = AJA + T5
                           T8      = T7 + T6
                           T9      = T7 - T6
                           T10     = C3 * T3 - C2 * T4
                           T11     = C2 * T3 + C3 * T4
                           BJB     = B(JB+J)
                           BJE     = B(JE+J)
                           U1      = BJB + BJE
                           BJC     = B(JC+J)
                           BJD     = B(JD+J)
                           U2      = BJC + BJD
                           U3      = BJB - BJE
                           U4      = BJC - BJD
                           U5      = U1 + U2
                           U6      = C1 * ( U1 - U2 )
                           BJA     = B(JA+J)
                           U7      = BJA - 0.25D0 * U5
                           B(JA+J) = BJA + U5
                           U8      = U7 + U6
                           U9      = U7 - U6
                           U10     = C3 * U3 - C2 * U4
                           U11     = C2 * U3 + C3 * U4
                           A(JB+J) = T8 - U11
                           B(JB+J) = U8 + T11
                           A(JE+J) = T8 + U11
                           B(JE+J) = U8 - T11
                           A(JC+J) = T9 - U10
                           B(JC+J) = U9 + T10
                           A(JD+J) = T9 + U10
                           B(JD+J) = U9 - T10
                           J = J + JUMP
 110                    CONTINUE
C
                     ELSE
C     
#if defined(_CRAY)
CDIR$ IVDEP,SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
                        DO 130 L = 1 , NVEX
                           AJB     = A(JB+J)
                           AJE     = A(JE+J)
                           T1      = AJB + AJE
                           AJC     = A(JC+J)
                           AJD     = A(JD+J)
                           T2      = AJC + AJD
                           T3      = AJB - AJE
                           T4      = AJC - AJD
                           T5      = T1 + T2
                           T6      = C1 * ( T1 - T2 )
                           AJA     = A(JA+J)
                           T7      = AJA - 0.25D0 * T5
                           A(JA+J) = AJA + T5
                           T8      = T7 + T6
                           T9      = T7 - T6
                           T10     = C3 * T3 - C2 * T4
                           T11     = C2 * T3 + C3 * T4
                           BJB     = B(JB+J)
                           BJE     = B(JE+J)
                           U1      = BJB + BJE
                           BJC     = B(JC+J)
                           BJD     = B(JD+J)
                           U2      = BJC + BJD
                           U3      = BJB - BJE
                           U4      = BJC - BJD
                           U5      = U1 + U2
                           U6      = C1 * ( U1 - U2 )
                           BJA     = B(JA+J)
                           U7      = BJA - 0.25D0 * U5
                           B(JA+J) = BJA + U5
                           U8      = U7 + U6
                           U9      = U7 - U6
                           U10     = C3 * U3 - C2 * U4
                           U11     = C2 * U3 + C3 * U4
                           A(JB+J) = CO1*(T8-U11) - SI1*(U8+T11)
                           B(JB+J) = SI1*(T8-U11) + CO1*(U8+T11)
                           A(JE+J) = CO4*(T8+U11) - SI4*(U8-T11)
                           B(JE+J) = SI4*(T8+U11) + CO4*(U8-T11)
                           A(JC+J) = CO2*(T9-U10) - SI2*(U9+T10)
                           B(JC+J) = SI2*(T9-U10) + CO2*(U9+T10)
                           A(JD+J) = CO3*(T9+U10) - SI3*(U9-T10)
                           B(JD+J) = SI3*(T9+U10) + CO3*(U9-T10)
                           J = J + JUMP
 130                    CONTINUE
C     
                     ENDIF
C     
C     -----(End of loop across transforms)
C     
                     JA = JA + JSTEPX
                     IF (JA.LT.ISTART) JA = JA + NINC
 135              CONTINUE
 140           CONTINUE
C     -----(End of loop along transforms)
               KK = KK + 2*LA
 150        CONTINUE
C     -----(End of loop on nonzero K)
            LA = 5*LA
 160     CONTINUE
C     -----(End of loop on type I radix-5 passes)
C     
         IF (N.EQ.5) GO TO 490
C     
C     Loop on type II radix-5 passes
C     ------------------------------
 400     CONTINUE
C     
         DO 480 IPASS = MH+1 , M
            JSTEP = (N*INC) / (5*LA)
            JSTEPL = JSTEP - NINC
            LAINCL = LA * INK - NINC
            KK = 0
C     
C     Loop on K
C     ---------
            DO 470 K = 0 , JSTEP-INK , INK
C     
               IF (K.GT.0) THEN
                  CO1 = TRIG(KK+1)
                  SI1 = S*TRIG(KK+2)
                  CO2 = TRIG(2*KK+1)
                  SI2 = S*TRIG(2*KK+2)
                  CO3 = TRIG(3*KK+1)
                  SI3 = S*TRIG(3*KK+2)
                  CO4 = TRIG(4*KK+1)
                  SI4 = S*TRIG(4*KK+2)
               ENDIF
C     
C     Double loop along first transform in block
C     ------------------------------------------
               DO 460 LL = K , (LA-1)*INK , 5*JSTEP
C     
                  DO 450 JJJ = LL , (N-1)*INC , 5*LA*INK
                     JA = ISTART + JJJ
C     
C     "Transverse" loop
C     -----------------
                     DO 445 NU = 1 , INQ
                        JB = JA + JSTEPL
                        IF (JB.LT.ISTART) JB = JB + NINC
                        JC = JB + JSTEPL
                        IF (JC.LT.ISTART) JC = JC + NINC
                        JD = JC + JSTEPL
                        IF (JD.LT.ISTART) JD = JD + NINC
                        JE = JD + JSTEPL
                        IF (JE.LT.ISTART) JE = JE + NINC
                        JF = JA + LAINCL
                        IF (JF.LT.ISTART) JF = JF + NINC
                        JG = JF + JSTEPL
                        IF (JG.LT.ISTART) JG = JG + NINC
                        JH = JG + JSTEPL
                        IF (JH.LT.ISTART) JH = JH + NINC
                        JI = JH + JSTEPL
                        IF (JI.LT.ISTART) JI = JI + NINC
                        JJ = JI + JSTEPL
                        IF (JJ.LT.ISTART) JJ = JJ + NINC
                        JK = JF + LAINCL
                        IF (JK.LT.ISTART) JK = JK + NINC
                        JL = JK + JSTEPL
                        IF (JL.LT.ISTART) JL = JL + NINC
                        JM = JL + JSTEPL
                        IF (JM.LT.ISTART) JM = JM + NINC
                        JN = JM + JSTEPL
                        IF (JN.LT.ISTART) JN = JN + NINC
                        JO = JN + JSTEPL
                        IF (JO.LT.ISTART) JO = JO + NINC
                        JP = JK + LAINCL
                        IF (JP.LT.ISTART) JP = JP + NINC
                        JQ = JP + JSTEPL
                        IF (JQ.LT.ISTART) JQ = JQ + NINC
                        JR = JQ + JSTEPL
                        IF (JR.LT.ISTART) JR = JR + NINC
                        JS = JR + JSTEPL
                        IF (JS.LT.ISTART) JS = JS + NINC
                        JT = JS + JSTEPL
                        IF (JT.LT.ISTART) JT = JT + NINC
                        JU = JP + LAINCL
                        IF (JU.LT.ISTART) JU = JU + NINC
                        JV = JU + JSTEPL
                        IF (JV.LT.ISTART) JV = JV + NINC
                        JW = JV + JSTEPL
                        IF (JW.LT.ISTART) JW = JW + NINC
                        JX = JW + JSTEPL
                        IF (JX.LT.ISTART) JX = JX + NINC
                        JY = JX + JSTEPL
                        IF (JY.LT.ISTART) JY = JY + NINC
                        J = 0
C     
C     Loop across transforms
C     ----------------------
                        IF (K.EQ.0) THEN
C     
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
                           DO 410 L = 1 , NVEX
                              AJB     = A(JB+J)
                              AJE     = A(JE+J)
                              T1      = AJB + AJE
                              AJC     = A(JC+J)
                              AJD     = A(JD+J)
                              T2      = AJC + AJD
                              T3      = AJB - AJE
                              T4      = AJC - AJD
                              AJF     = A(JF+J)
                              AJB     =  AJF
                              T5      = T1 + T2
                              T6      = C1 * ( T1 - T2 )
                              AJA     = A(JA+J)
                              T7      = AJA - 0.25D0 * T5
                              A(JA+J) = AJA + T5
                              T8      = T7 + T6
                              T9      = T7 - T6
                              AJK     = A(JK+J)
                              AJC     = AJK
                              T10     = C3 * T3 - C2 * T4
                              T11     = C2 * T3 + C3 * T4
                              BJB     = B(JB+J)
                              BJE     = B(JE+J)
                              U1      = BJB + BJE
                              BJC     = B(JC+J)
                              BJD     = B(JD+J)
                              U2      = BJC + BJD
                              U3      = BJB - BJE
                              U4      = BJC - BJD
                              BJF     = B(JF+J)
                              BJB     = BJF
                              U5      = U1 + U2
                              U6      = C1 * ( U1 - U2 )
                              BJA     = B(JA+J)
                              U7      = BJA - 0.25D0 * U5
                              B(JA+J) = BJA + U5
                              U8      = U7 + U6
                              U9      = U7 - U6
                              BJK     = B(JK+J)
                              BJC     =  BJK
                              U10     = C3 * U3 - C2 * U4
                              U11     = C2 * U3 + C3 * U4
                              A(JF+J) = T8 - U11
                              B(JF+J) = U8 + T11
                              AJE     = T8 + U11
                              BJE     = U8 - T11
                              A(JK+J) = T9 - U10
                              B(JK+J) = U9 + T10
                              AJD     = T9 + U10
                              BJD     = U9 - T10
C     ----------------------
                              AJG     = A(JG+J)
                              AJJ     = A(JJ+J)
                              T1      = AJG + AJJ
                              AJH     = A(JH+J)
                              AJI     = A(JI+J)
                              T2      = AJH + AJI
                              T3      = AJG - AJJ
                              T4      = AJH - AJI
                              AJL     = A(JL+J)
                              AJH     = AJL
                              T5      = T1 + T2
                              T6      = C1 * ( T1 - T2 )
                              T7      = AJB - 0.25D0 * T5
                              A(JB+J) = AJB + T5
                              T8      = T7 + T6
                              T9      = T7 - T6
                              AJQ     = A(JQ+J)
                              AJI     = AJQ
                              T10     = C3 * T3 - C2 * T4
                              T11     = C2 * T3 + C3 * T4
                              BJG     = B(JG+J)
                              BJJ     = B(JJ+J)
                              U1      = BJG + BJJ
                              BJH     = B(JH+J)
                              BJI     = B(JI+J)
                              U2      = BJH + BJI
                              U3      = BJG - BJJ
                              U4      = BJH - BJI
                              BJL     = B(JL+J)
                              BJH     = BJL
                              U5      = U1 + U2
                              U6      = C1 * ( U1 - U2 )
                              U7      = BJB - 0.25D0 * U5
                              B(JB+J) = BJB + U5
                              U8      = U7 + U6
                              U9      = U7 - U6
                              BJQ     = B(JQ+J)
                              BJI     = BJQ
                              U10     = C3 * U3 - C2 * U4
                              U11     = C2 * U3 + C3 * U4
                              A(JG+J) = T8 - U11
                              B(JG+J) = U8 + T11
                              AJJ     = T8 + U11
                              BJJ     = U8 - T11
                              A(JL+J) = T9 - U10
                              B(JL+J) = U9 + T10
                              A(JQ+J) = T9 + U10
                              B(JQ+J) = U9 - T10
C     ----------------------
                              AJO     = A(JO+J)
                              T1      = AJH + AJO
                              AJM     = A(JM+J)
                              AJN     = A(JN+J)
                              T2      = AJM + AJN
                              T3      = AJH - AJO
                              T4      = AJM - AJN
                              AJR     = A(JR+J)
                              AJN     = AJR
                              T5      = T1 + T2
                              T6      = C1 * ( T1 - T2 )
                              T7      = AJC - 0.25D0 * T5
                              A(JC+J) = AJC + T5
                              T8      = T7 + T6
                              T9      = T7 - T6
                              AJW     = A(JW+J)
                              AJO     = AJW
                              T10     = C3 * T3 - C2 * T4
                              T11     = C2 * T3 + C3 * T4
                              BJO     = B(JO+J)
                              U1      = BJH + BJO
                              BJM     = B(JM+J)
                              BJN     = B(JN+J)
                              U2      = BJM + BJN
                              U3      = BJH - BJO
                              U4      = BJM - BJN
                              BJR     = B(JR+J)
                              BJN     = BJR
                              U5      = U1 + U2
                              U6      = C1 * ( U1 - U2 )
                              U7      = BJC - 0.25D0 * U5
                              B(JC+J) = BJC + U5
                              U8      = U7 + U6
                              U9      = U7 - U6
                              BJW     = B(JW+J)
                              BJO     = BJW
                              U10     = C3 * U3 - C2 * U4
                              U11     = C2 * U3 + C3 * U4
                              A(JH+J) = T8 - U11
                              B(JH+J) = U8 + T11
                              A(JW+J) = T8 + U11
                              B(JW+J) = U8 - T11
                              A(JM+J) = T9 - U10
                              B(JM+J) = U9 + T10
                              A(JR+J) = T9 + U10
                              B(JR+J) = U9 - T10
C     ----------------------
                              AJT     = A(JT+J)
                              T1      = AJI + AJT
                              AJS     = A(JS+J)
                              T2      = AJN + AJS
                              T3      = AJI - AJT
                              T4      = AJN - AJS
                              AJX     = A(JX+J)
                              AJT     = AJX
                              T5      = T1 + T2
                              T6      = C1 * ( T1 - T2 )
                              AJP     = A(JP+J)
                              T7      = AJP - 0.25D0 * T5
                              AX      = AJP + T5
                              T8      = T7 + T6
                              T9      = T7 - T6
                              A(JP+J) = AJD
                              T10     = C3 * T3 - C2 * T4
                              T11     = C2 * T3 + C3 * T4
                              A(JD+J) = AX
                              BJT     = B(JT+J)
                              U1      = BJI + BJT
                              BJS     = B(JS+J)
                              U2      = BJN + BJS
                              U3      = BJI - BJT
                              U4      = BJN - BJS
                              BJX     = B(JX+J)
                              BJT     = BJX
                              U5      = U1 + U2
                              U6      = C1 * ( U1 - U2 )
                              BJP     = B(JP+J)
                              U7      = BJP - 0.25D0 * U5
                              BX      = BJP + U5
                              U8      = U7 + U6
                              U9      = U7 - U6
                              B(JP+J) = BJD
                              U10     = C3 * U3 - C2 * U4
                              U11     = C2 * U3 + C3 * U4
                              B(JD+J) = BX
                              A(JI+J) = T8 - U11
                              B(JI+J) = U8 + T11
                              A(JX+J) = T8 + U11
                              B(JX+J) = U8 - T11
                              A(JN+J) = T9 - U10
                              B(JN+J) = U9 + T10
                              A(JS+J) = T9 + U10
                              B(JS+J) = U9 - T10
C     ----------------------
                              AJV     = A(JV+J)
                              AJY     = A(JY+J)
                              T1      = AJV + AJY
                              T2      = AJO + AJT
                              T3      = AJV - AJY
                              T4      = AJO - AJT
                              A(JV+J) = AJJ
                              T5      = T1 + T2
                              T6      = C1 * ( T1 - T2 )
                              AJU     = A(JU+J)
                              T7      = AJU - 0.25D0 * T5
                              AX      = AJU + T5
                              T8      = T7 + T6
                              T9      = T7 - T6
                              A(JU+J) = AJE
                              T10     = C3 * T3 - C2 * T4
                              T11     = C2 * T3 + C3 * T4
                              A(JE+J) = AX
                              BJV     = B(JV+J)
                              BJY     = B(JY+J)
                              U1      = BJV + BJY
                              U2      = BJO + BJT
                              U3      = BJV - BJY
                              U4      = BJO - BJT
                              B(JV+J) = BJJ
                              U5      = U1 + U2
                              U6      = C1 * ( U1 - U2 )
                              BJU     = B(JU+J)
                              U7      = BJU - 0.25D0 * U5
                              BX      = BJU + U5
                              U8      = U7 + U6
                              U9      = U7 - U6
                              B(JU+J) = BJE
                              U10     = C3 * U3 - C2 * U4
                              U11     = C2 * U3 + C3 * U4
                              B(JE+J) = BX
                              A(JJ+J) = T8 - U11
                              B(JJ+J) = U8 + T11
                              A(JY+J) = T8 + U11
                              B(JY+J) = U8 - T11
                              A(JO+J) = T9 - U10
                              B(JO+J) = U9 + T10
                              A(JT+J) = T9 + U10
                              B(JT+J) = U9 - T10
                              J = J + JUMP
 410                       CONTINUE
C     
                        ELSE
C     
#if defined(_CRAY)
CDIR$ IVDEP, SHORTLOOP
#elif defined(_SX)
*VDIR NODEP, SHORTLOOP
#elif defined(__uxp__)
!OCL NOVREC
#endif
                           DO 440 L = 1 , NVEX
                              AJB     = A(JB+J)
                              AJE     = A(JE+J)
                              T1      = AJB + AJE
                              AJC     = A(JC+J)
                              AJD     = A(JD+J)
                              T2      = AJC + AJD
                              T3      = AJB - AJE
                              T4      = AJC - AJD
                              AJF     = A(JF+J)
                              AJB     = AJF
                              T5      = T1 + T2
                              T6      = C1 * ( T1 - T2 )
                              AJA     = A(JA+J)
                              T7      = AJA - 0.25D0 * T5
                              A(JA+J) = AJA + T5
                              T8      = T7 + T6
                              T9      = T7 - T6
                              AJK     = A(JK+J)
                              AJC     = AJK
                              T10     = C3 * T3 - C2 * T4
                              T11     = C2 * T3 + C3 * T4
                              BJB     = B(JB+J)
                              BJE     = B(JE+J)
                              U1      = BJB + BJE
                              BJC     = B(JC+J)
                              BJD     = B(JD+J)
                              U2      = BJC + BJD
                              U3      = BJB - BJE
                              U4      = BJC - BJD
                              BJF     = B(JF+J)
                              BJB     =  BJF
                              U5      = U1 + U2
                              U6      = C1 * ( U1 - U2 )
                              BJA     = B(JA+J)
                              U7      = BJA - 0.25D0 * U5
                              B(JA+J) = BJA + U5
                              U8      = U7 + U6
                              U9      = U7 - U6
                              BJK     = B(JK+J)
                              BJC     = BJK
                              U10     = C3 * U3 - C2 * U4
                              U11     = C2 * U3 + C3 * U4
                              A(JF+J) = CO1*(T8-U11) - SI1*(U8+T11)
                              B(JF+J) = SI1*(T8-U11) + CO1*(U8+T11)
                              AJE     = CO4*(T8+U11) - SI4*(U8-T11)
                              BJE     = SI4*(T8+U11) + CO4*(U8-T11)
                              A(JK+J) = CO2*(T9-U10) - SI2*(U9+T10)
                              B(JK+J) = SI2*(T9-U10) + CO2*(U9+T10)
                              AJD     = CO3*(T9+U10) - SI3*(U9-T10)
                              BJD     = SI3*(T9+U10) + CO3*(U9-T10)
C     ----------------------
                              AJG     = A(JG+J)
                              AJJ     = A(JJ+J)
                              T1      = AJG + AJJ
                              AJH     = A(JH+J)
                              AJI     = A(JI+J)
                              T2      = AJH + AJI
                              T3      = AJG - AJJ
                              T4      = AJH - AJI
                              AJL     = A(JL+J)
                              AJH     = AJL
                              T5      = T1 + T2
                              T6      = C1 * ( T1 - T2 )
                              T7      = AJB - 0.25D0 * T5
                              A(JB+J) = AJB + T5
                              T8      = T7 + T6
                              T9      = T7 - T6
                              AJQ     = A(JQ+J)
                              AJI     = AJQ
                              T10     = C3 * T3 - C2 * T4
                              T11     = C2 * T3 + C3 * T4
                              BJG     = B(JG+J)
                              BJJ     = B(JJ+J)
                              U1      = BJG + BJJ
                              BJH     = B(JH+J)
                              BJI     = B(JI+J)
                              U2      = BJH + BJI
                              U3      = BJG - BJJ
                              U4      = BJH - BJI
                              BJL     = B(JL+J)
                              BJH     = BJL
                              U5      = U1 + U2
                              U6      = C1 * ( U1 - U2 )
                              U7      = BJB - 0.25D0 * U5
                              B(JB+J) = BJB + U5
                              U8      = U7 + U6
                              U9      = U7 - U6
                              BJQ     = B(JQ+J)
                              BJI     = BJQ
                              U10     = C3 * U3 - C2 * U4
                              U11     = C2 * U3 + C3 * U4
                              A(JG+J) = CO1*(T8-U11) - SI1*(U8+T11)
                              B(JG+J) = SI1*(T8-U11) + CO1*(U8+T11)
                              AJJ     = CO4*(T8+U11) - SI4*(U8-T11)
                              BJJ     = SI4*(T8+U11) + CO4*(U8-T11)
                              A(JL+J) = CO2*(T9-U10) - SI2*(U9+T10)
                              B(JL+J) = SI2*(T9-U10) + CO2*(U9+T10)
                              A(JQ+J) = CO3*(T9+U10) - SI3*(U9-T10)
                              B(JQ+J) = SI3*(T9+U10) + CO3*(U9-T10)
C     ----------------------
                              AJO     = A(JO+J)
                              T1      = AJH + AJO
                              AJM     = A(JM+J)
                              AJN     = A(JN+J)
                              T2      = AJM + AJN
                              T3      = AJH - AJO
                              T4      = AJM - AJN
                              AJR     = A(JR+J)
                              AJN     = AJR
                              T5      = T1 + T2
                              T6      = C1 * ( T1 - T2 )
                              T7      = AJC - 0.25D0 * T5
                              A(JC+J) = AJC + T5
                              T8      = T7 + T6
                              T9      = T7 - T6
                              AJW     = A(JW+J)
                              AJO     = AJW
                              T10     = C3 * T3 - C2 * T4
                              T11     = C2 * T3 + C3 * T4
                              BJO     = B(JO+J)
                              U1      = BJH + BJO
                              BJM     = B(JM+J)
                              BJN     = B(JN+J)
                              U2      = BJM + BJN
                              U3      = BJH - BJO
                              U4      = BJM - BJN
                              BJR     = B(JR+J)
                              BJN     = BJR
                              U5      = U1 + U2
                              U6      = C1 * ( U1 - U2 )
                              U7      = BJC - 0.25D0 * U5
                              B(JC+J) = BJC + U5
                              U8      = U7 + U6
                              U9      = U7 - U6
                              BJW     = B(JW+J)
                              BJO     = BJW
                              U10     = C3 * U3 - C2 * U4
                              U11     = C2 * U3 + C3 * U4
                              A(JH+J) = CO1*(T8-U11) - SI1*(U8+T11)
                              B(JH+J) = SI1*(T8-U11) + CO1*(U8+T11)
                              A(JW+J) = CO4*(T8+U11) - SI4*(U8-T11)
                              B(JW+J) = SI4*(T8+U11) + CO4*(U8-T11)
                              A(JM+J) = CO2*(T9-U10) - SI2*(U9+T10)
                              B(JM+J) = SI2*(T9-U10) + CO2*(U9+T10)
                              A(JR+J) = CO3*(T9+U10) - SI3*(U9-T10)
                              B(JR+J) = SI3*(T9+U10) + CO3*(U9-T10)
C     ----------------------
                              AJT     = A(JT+J)
                              T1      = AJI + AJT
                              AJS     = A(JS+J)
                              T2      = AJN + AJS
                              T3      = AJI - AJT
                              T4      = AJN - AJS
                              AJX     = A(JX+J)
                              AJT     = AJX
                              T5      = T1 + T2
                              T6      = C1 * ( T1 - T2 )
                              AJP     = A(JP+J)
                              T7      = AJP - 0.25D0 * T5
                              AX      = AJP + T5
                              T8      = T7 + T6
                              T9      = T7 - T6
                              A(JP+J) = AJD
                              T10     = C3 * T3 - C2 * T4
                              T11     = C2 * T3 + C3 * T4
                              A(JD+J) = AX
                              BJT     = B(JT+J)
                              U1      = BJI + BJT
                              BJS     = B(JS+J)
                              U2      = BJN + BJS
                              U3      = BJI - BJT
                              U4      = BJN - BJS
                              BJX     = B(JX+J)
                              BJT     = BJX
                              U5      = U1 + U2
                              U6      = C1 * ( U1 - U2 )
                              BJP     = B(JP+J)
                              U7      = BJP - 0.25D0 * U5
                              BX      = BJP + U5
                              U8      = U7 + U6
                              U9      = U7 - U6
                              B(JP+J) = BJD
                              U10     = C3 * U3 - C2 * U4
                              U11     = C2 * U3 + C3 * U4
                              B(JD+J) = BX
                              A(JI+J) = CO1*(T8-U11) - SI1*(U8+T11)
                              B(JI+J) = SI1*(T8-U11) + CO1*(U8+T11)
                              A(JX+J) = CO4*(T8+U11) - SI4*(U8-T11)
                              B(JX+J) = SI4*(T8+U11) + CO4*(U8-T11)
                              A(JN+J) = CO2*(T9-U10) - SI2*(U9+T10)
                              B(JN+J) = SI2*(T9-U10) + CO2*(U9+T10)
                              A(JS+J) = CO3*(T9+U10) - SI3*(U9-T10)
                              B(JS+J) = SI3*(T9+U10) + CO3*(U9-T10)
C     ----------------------
                              AJV     = A(JV+J)
                              AJY     = A(JY+J)
                              T1      = AJV + AJY
                              T2      = AJO + AJT
                              T3      = AJV - AJY
                              T4      = AJO - AJT
                              A(JV+J) = AJJ
                              T5      = T1 + T2
                              T6      = C1 * ( T1 - T2 )
                              AJU     = A(JU+J)
                              T7      = AJU - 0.25D0 * T5
                              AX      = AJU + T5
                              T8      = T7 + T6
                              T9      = T7 - T6
                              A(JU+J) = AJE
                              T10     = C3 * T3 - C2 * T4
                              T11     = C2 * T3 + C3 * T4
                              A(JE+J) = AX
                              BJV     = B(JV+J)
                              BJY     = B(JY+J)
                              U1      = BJV + BJY
                              U2      = BJO + BJT
                              U3      = BJV - BJY
                              U4      = BJO - BJT
                              B(JV+J) = BJJ
                              U5      = U1 + U2
                              U6      = C1 * ( U1 - U2 )
                              BJU     = B(JU+J)
                              U7      = BJU - 0.25D0 * U5
                              BX      = BJU + U5
                              U8      = U7 + U6
                              U9      = U7 - U6
                              B(JU+J) = BJE
                              U10     = C3 * U3 - C2 * U4
                              U11     = C2 * U3 + C3 * U4
                              B(JE+J) = BX
                              A(JJ+J) = CO1*(T8-U11) - SI1*(U8+T11)
                              B(JJ+J) = SI1*(T8-U11) + CO1*(U8+T11)
                              A(JY+J) = CO4*(T8+U11) - SI4*(U8-T11)
                              B(JY+J) = SI4*(T8+U11) + CO4*(U8-T11)
                              A(JO+J) = CO2*(T9-U10) - SI2*(U9+T10)
                              B(JO+J) = SI2*(T9-U10) + CO2*(U9+T10)
                              A(JT+J) = CO3*(T9+U10) - SI3*(U9-T10)
                              B(JT+J) = SI3*(T9+U10) + CO3*(U9-T10)
                              J = J + JUMP
 440                       CONTINUE
C     
                        ENDIF
C     
C     -----(End of loop across transforms)
C     
                        JA = JA + JSTEPX
                        IF (JA.LT.ISTART) JA = JA + NINC
 445                 CONTINUE
 450              CONTINUE
 460           CONTINUE
C     -----(End of double loop for this K)
               KK = KK + 2*LA
 470        CONTINUE
C     -----(End of loop over values of K)
            LA = 5*LA
 480     CONTINUE
C     -----(End of loop on type II radix-5 passes)
C     -----(Nvex transforms completed)
 490     CONTINUE
         ISTART = ISTART + NVEX * JUMP
 500  CONTINUE
C     -----(End of loop on blocks of transforms)
C     
      RETURN
      END
