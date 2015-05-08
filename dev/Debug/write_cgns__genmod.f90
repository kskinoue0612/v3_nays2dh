        !COMPILER-GENERATED INTERFACE MODULE: Thu Apr 23 12:20:38 2015
        MODULE WRITE_CGNS__genmod
          INTERFACE 
            SUBROUTINE WRITE_CGNS(INPUTFILE,TIME,DISCH,IM,JM,X,Y,U,V,HS,&
     &Z,Z0,ZB,VORT,C,DMN,PHI,FR,TS,ZAVE,ZMIN,HAVE,QBX,QBY,CC_M,NK,J_MIX)
              INTEGER(KIND=4), INTENT(IN) :: NK
              INTEGER(KIND=4), INTENT(IN) :: JM
              INTEGER(KIND=4), INTENT(IN) :: IM
              CHARACTER(*), INTENT(IN) :: INPUTFILE
              REAL(KIND=8), INTENT(IN) :: TIME
              REAL(KIND=8), INTENT(IN) :: DISCH
              REAL(KIND=8), INTENT(IN) :: X(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: Y(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: U(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: V(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: HS(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: Z(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: Z0(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: ZB(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: VORT(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: C(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: DMN(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: PHI(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: FR(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: TS(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: ZAVE(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: ZMIN(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: HAVE(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: QBX(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: QBY(0:IM,0:JM)
              REAL(KIND=8), INTENT(IN) :: CC_M(0:IM,0:JM,NK)
              INTEGER(KIND=4), INTENT(IN) :: J_MIX
            END SUBROUTINE WRITE_CGNS
          END INTERFACE 
        END MODULE WRITE_CGNS__genmod
