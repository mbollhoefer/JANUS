#include <ilupack_fortran.h>
C----- SUBROUTINE GENRCM
C***************************************************************         
C***************************************************************         
C********   GENRCM ..... GENERAL REVERSE CUTHILL MCKEE   *******         
C***************************************************************         
C***************************************************************         
C                                                                        
C     PURPOSE - GENRCM FINDS THE REVERSE CUTHILL-MCKEE                   
C        ORDERING FOR A GENERAL GRAPH. FOR EACH CONNECTED                
C        COMPONENT IN THE GRAPH, GENRCM OBTAINS THE ORDERING             
C        BY CALLING THE SUBROUTINE RCM.                                  
C                                                                        
C     INPUT PARAMETERS -                                                 
C        NEQNS - NUMBER OF EQUATIONS                                     
C        (XADJ, ADJNCY) - ARRAY PAIR CONTAINING THE ADJACENCY            
C               STRUCTURE OF THE GRAPH OF THE MATRIX.                    
C                                                                        
C     OUTPUT PARAMETER -                                                 
C        PERM - VECTOR THAT CONTAINS THE RCM ORDERING.                   
C                                                                        
C     WORKING PARAMETERS -                                               
C        MASK - IS USED TO MARK VARIABLES THAT HAVE BEEN                 
C               NUMBERED DURING THE ORDERING PROCESS. IT IS              
C               INITIALIZED TO 1, AND SET TO ZERO AS EACH NODE           
C               IS NUMBERED.                                             
C        XLS - THE INDEX VECTOR FOR A LEVEL STRUCTURE.  THE              
C               LEVEL STRUCTURE IS STORED IN THE CURRENTLY               
C               UNUSED SPACES IN THE PERMUTATION VECTOR PERM.            
C                                                                        
C     PROGRAM SUBROUTINES -                                              
C        FNROOT, RCM.                                                    
C                                                                        
C***************************************************************         
C                                                                        
      SUBROUTINE  GENRCM ( NEQNS, XADJ, ADJNCY, PERM, MASK, XLS )        
C                                                                        
C***************************************************************         
C                                                                        
         integer ADJNCY(1), MASK(1), PERM(1), XLS(1)                     
         integer XADJ(1), CCSIZE, I, NEQNS, NLVL,                        
     1           NUM, ROOT                                               
C                                                                        
C***************************************************************         
C                                                                        
         DO 100 I = 1, NEQNS                                             
            MASK(I) = 1                                                  
  100    CONTINUE                                                        
         NUM = 1                                                         
         DO 200 I = 1, NEQNS                                             
C           ---------------------------------------                      
C           FOR EACH MASKED CONNECTED COMPONENT ...                      
C           ---------------------------------------                      
            IF (MASK(I) .EQ. 0) GO TO 200                                
               ROOT = I                                                  
C              -----------------------------------------                 
C              FIRST FIND A PSEUDO-PERIPHERAL NODE ROOT.                 
C              NOTE THAT THE LEVEL STRUCTURE FOUND BY                    
C              FNROOT IS STORED STARTING AT PERM(NUM).                   
C              THEN RCM IS CALLED TO ORDER THE COMPONENT                 
C              USING ROOT AS THE STARTING NODE.                          
C              -----------------------------------------                 
               CALL  FNROOT ( ROOT, XADJ, ADJNCY, MASK,                  
     1                        NLVL, XLS, PERM(NUM) )                     
               CALL     RCM ( ROOT, XADJ, ADJNCY, MASK,                  
     1                        PERM(NUM), CCSIZE, XLS )                   
               NUM = NUM + CCSIZE                                        
               IF (NUM .GT. NEQNS) RETURN                                
  200    CONTINUE                                                        
         RETURN                                                          
      END                                                                
