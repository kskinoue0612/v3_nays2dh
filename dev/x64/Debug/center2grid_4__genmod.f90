        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul 13 12:37:18 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CENTER2GRID_4__genmod
          INTERFACE 
            SUBROUTINE CENTER2GRID_4(F_CEN,F_GRID,NK)
              USE COMMON_HH
              INTEGER(KIND=4), INTENT(IN) :: NK
              REAL(KIND=8), INTENT(IN) :: F_CEN(0:IM,0:JM,NK)
              REAL(KIND=8), INTENT(OUT) :: F_GRID(0:IM,0:JM,NK)
            END SUBROUTINE CENTER2GRID_4
          END INTERFACE 
        END MODULE CENTER2GRID_4__genmod
