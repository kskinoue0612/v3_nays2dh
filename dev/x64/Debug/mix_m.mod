	  �%  _   k820309    �          15.0        3:U                                                                                                           
       D:\01_iRIC\v3.x\solvers\nays2dh\src\Nays2DH.f90 MIX_M                                                    
                                                          
                                                                                                                                                                                           @ @                                   
                  @ @                                   
                                                                                                             	                                                         
                                                                                                                                                                               
                                                        
                                                        
                                                        
                                                                    
                &                   &                                                                                                                                                     
                                                        
                                                                    
                &                   &                                                                                                             @                                                   
                &                   &                                                    @                                                   
                &                   &                                                    @                                                   
                &                   &                                                    @                                                   
                &                   &                                                    @                                                   
                &                   &                   &                                                    @                                                   
                &                   &                   &                                                    @                                                   
                &                   &                   &                                                    @                                                   
                &                   &                   &                                                    @                                                   
                &                   &                   &                                                    @                                                    
                &                   &                   &                                           #         @                                   !                    #NK "             
                                 "           #         @                                   #                   #MIXINI%DSQRT $   #MIXINI%DLOG10 %   #SNU00 &   #DM0 '                                             $     DSQRT                                           %     DLOG10           
@ @                              &     
                D                                '     
       #         @                                   (                    #SNU00 )             
  @                              )     
      #         @                                  *                   #DMTSCM%DSQRT +   #SNU00 ,                                             +     DSQRT           
@ @                              ,     
      #         @                                   -                     #         @                                   .                   #QBCAL_W_MIX%DSIN /   #QBCAL_W_MIX%DATAN2 0   #QBCAL_W_MIX%DCOS 1   #QBCAL_W_MIX%DTAN 2   #QBCAL_W_MIX%DMAX1 3   #QBCAL_W_MIX%DATAN 4   #QBCAL_W_MIX%DABS 5   #QBCAL_W_MIX%DSQRT 6   #UX0 7   #UY0 8   #HS0 9   #GAMMA_M :   #DSMT ;   #PI_BED <   #TANTC =   #J_BANK >   #I_EROSION_START ?   #I_EROSION_END @   #BHEIGHT A                                             /     DSIN                                           0     DATAN2                                           1     DCOS                                           2     DTAN                                           3     DMAX1                                           4     DATAN                                           5     DABS                                           6     DSQRT          
      �                           7                    
      p           5 r    p         p        p           & p         5 r      & p         5 r          5 r    p         p            5 r    p         p                                   
      �                           8                    
      p           5 r    p         p        p           & p         5 r      & p         5 r          5 r    p         p            5 r    p         p                                   
  @   �                           9                    
      p           5 r    p         p        p           & p         5 r      & p         5 r          5 r    p         p            5 r    p         p                                    
                                 :     
                
                                 ;     
                
                                 <     
                
                                 =     
                
                                  >                     
                                  ?                     
                                  @                     
                                 A     
      #         @                                   B                   #ETACAL_MIX%NK C   #DSMT D                                               C                      
                                 D     
      #         @                                  E                   #SORTING_FIXED%NK F   #I G   #J H   #NB_NEW I   #E_T_NEW J   #P_M_NEW K   #P_T_NEW L   #P_D_NEW M                                               F                      
                                  G                     
                                  H                     D                                 I                      D                                 J     
                D                                 K                    
     p          5 r F       5 r F                              D                                 L                    
     p          5 r F       5 r F                              D                                 M                    
     p          5 r F       5 r F                     #         @                                  N                   #SORTING_MOVABLE%NK O   #I P   #J Q   #NB_NEW R   #E_T_NEW S   #P_M_NEW T   #P_T_NEW U   #P_D_NEW V                                               O                      
                                  P                     
                                  Q                     D                                 R                      D                                 S     
                D                                 T                    
     p          5 r O       5 r O                              D                                 U                    
     p          5 r O       5 r O                              D                                 V                    
     p          5 r O       5 r O                     #         @                                   W                   #ETACAL_MIX_C%NK X   #DSMT Y                                               X                      
                                 Y     
      #         @                                   Z                   #C_TRANSPORT_MIX%MAX [   #C_TRANSPORT_MIX%DABS \   #DSMT ]                                             [     MAX                                           \     DABS           
                                 ]     
      #         @                                   ^                        �   >      fn#fn    �   @   j   COMMON_HH      @   J   FIXED_BED    ^  @       IM+COMMON_HH    �  @       JM+COMMON_HH $   �  @       J_MIX_DIS+COMMON_HH      @       SPEC+COMMON_HH    ^  @       G+COMMON_HH "   �  @       NM_CELL+COMMON_HH (   �  @       J_MIX_DIS_DEP+COMMON_HH      @       NY+COMMON_HH    ^  @       NX+COMMON_HH    �  @       JREP+COMMON_HH    �  @       HMIN+COMMON_HH       @       R_DXI+COMMON_HH     ^  @       R_DET+COMMON_HH    �  @       MU_S+COMMON_HH    �  �       PHI+FIXED_BED #   �  @       J_QB_VEC+COMMON_HH    �  @       DT+COMMON_HH      @       CSM+COMMON_HH    B  �       EMB+FIXED_BED !   �  @       J_QBUP+COMMON_HH    &  �       UBXI    �  �       UBET    n  �       UANG    	  �       DZDN    �	  �       QBTI_MIX    r
  �       QUCK    .  �       QVCK    �  �       DCDXI_K    �  �       DCDET_K    b  �       SOURCEK )     P       ALLOC_MIX_TEMP_VARIABLES ,   n  @   a   ALLOC_MIX_TEMP_VARIABLES%NK    �  �       MIXINI    /  >      MIXINI%DSQRT    m  ?      MIXINI%DLOG10    �  @   a   MIXINI%SNU00    �  @   a   MIXINI%DM0    ,  S       INI_LAYER       @   a   INI_LAYER%SNU00    �  e       DMTSCM    $  >      DMTSCM%DSQRT    b  @   a   DMTSCM%SNU00    �  H       DMCAL    �  �      QBCAL_W_MIX !   q  =      QBCAL_W_MIX%DSIN #   �  ?      QBCAL_W_MIX%DATAN2 !   �  =      QBCAL_W_MIX%DCOS !   *  =      QBCAL_W_MIX%DTAN "   g  >      QBCAL_W_MIX%DMAX1 "   �  >      QBCAL_W_MIX%DATAN !   �  =      QBCAL_W_MIX%DABS "      >      QBCAL_W_MIX%DSQRT     ^  T  a   QBCAL_W_MIX%UX0     �  T  a   QBCAL_W_MIX%UY0       T  a   QBCAL_W_MIX%HS0 $   Z  @   a   QBCAL_W_MIX%GAMMA_M !   �  @   a   QBCAL_W_MIX%DSMT #   �  @   a   QBCAL_W_MIX%PI_BED "     @   a   QBCAL_W_MIX%TANTC #   Z  @   a   QBCAL_W_MIX%J_BANK ,   �  @   a   QBCAL_W_MIX%I_EROSION_START *   �  @   a   QBCAL_W_MIX%I_EROSION_END $     @   a   QBCAL_W_MIX%BHEIGHT    Z  e       ETACAL_MIX "   �  @     ETACAL_MIX%NK+MIX     �  @   a   ETACAL_MIX%DSMT    ?  �       SORTING_FIXED %   �  @     SORTING_FIXED%NK+MIX     +  @   a   SORTING_FIXED%I     k  @   a   SORTING_FIXED%J %   �  @   a   SORTING_FIXED%NB_NEW &   �  @   a   SORTING_FIXED%E_T_NEW &   +  �   a   SORTING_FIXED%P_M_NEW &   �  �   a   SORTING_FIXED%P_T_NEW &   S  �   a   SORTING_FIXED%P_D_NEW     �  �       SORTING_MOVABLE '   �   @     SORTING_MOVABLE%NK+MIX "   �   @   a   SORTING_MOVABLE%I "   !  @   a   SORTING_MOVABLE%J '   U!  @   a   SORTING_MOVABLE%NB_NEW (   �!  @   a   SORTING_MOVABLE%E_T_NEW (   �!  �   a   SORTING_MOVABLE%P_M_NEW (   i"  �   a   SORTING_MOVABLE%P_T_NEW (   �"  �   a   SORTING_MOVABLE%P_D_NEW    �#  g       ETACAL_MIX_C $   �#  @     ETACAL_MIX_C%NK+MIX "   8$  @   a   ETACAL_MIX_C%DSMT     x$  �       C_TRANSPORT_MIX $   �$  <      C_TRANSPORT_MIX%MAX %   9%  =      C_TRANSPORT_MIX%DABS %   v%  @   a   C_TRANSPORT_MIX%DSMT    �%  H       BOUND_C_MIX 