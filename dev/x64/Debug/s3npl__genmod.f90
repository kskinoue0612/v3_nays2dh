        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul 13 12:37:19 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE S3NPL__genmod
          INTERFACE 
            SUBROUTINE S3NPL(X,Y,S,XX,YY,SS,N,NN)
              INTEGER(KIND=4) :: NN
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(0:N)
              REAL(KIND=8) :: Y(0:N)
              REAL(KIND=8) :: S(0:N)
              REAL(KIND=8) :: XX(0:NN)
              REAL(KIND=8) :: YY(0:NN)
              REAL(KIND=8) :: SS(0:NN)
            END SUBROUTINE S3NPL
          END INTERFACE 
        END MODULE S3NPL__genmod
