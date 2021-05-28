#include <ilupack_fortran.h>
C----- SUBROUTINE GENND
C***************************************************************         
C***************************************************************         
C********     GENND ..... GENERAL NESTED DISSECTION     ********         
C***************************************************************         
C***************************************************************         
C                                                                        
C     PURPOSE - SUBROUTINE GENND FINDS A NESTED DISSECTION               
C        ORDERING FOR A GENERAL GRAPH.                                   
C                                                                        
C                                                                        
C     INPUT PARAMETERS -                                                 
C        NEQNS - NUMBER OF EQUATIONS.                                    
C        (XADJ, ADJNCY) - ADJACENCY STRUCTURE PAIR.                      
C                                                                        
C     OUTPUT PARAMETERS -                                                
C        PERM - THE NESTED DISSECTION ORDERING.                          
C                                                                        
C     WORKING PARAMETERS -                                               
C        MASK - IS USED TO MASK OFF VARIABLES THAT HAVE                  
C               BEEN NUMBERED DURING THE ORDERNG PROCESS.                
C        (XLS, LS) - THIS LEVEL STRUCTURE PAIR IS USED AS                
C               TEMPORARY STORAGE BY FNROOT.                             
C                                                                        
C     PROGRAM SUBROUTINES -                                              
C        FNDSEP, REVRSE.                                                 
C                                                                        
C***************************************************************         
C                                                                        
       SUBROUTINE  GENND ( NEQNS, XADJ, ADJNCY, MASK,                    
     1                     PERM, XLS, LS )                               
C                                                                        
C***************************************************************         
C                                                                        
          integer ADJNCY(1), MASK(1), LS(1), PERM(1),                    
     1            XLS(1)                                                 
          integer XADJ(1), I, NEQNS, NSEP, NUM, ROOT                     
C                                                                        
C***************************************************************         
C                                                                        
          DO 100 I = 1, NEQNS                                            
             MASK(I) = 1                                                 
  100     CONTINUE                                                       
          NUM   = 0                                                      
          DO 300 I = 1, NEQNS                                            
C            -----------------------------                               
C            FOR EACH MASKED COMPONENT ...                               
C            -----------------------------                               
  200        IF ( MASK(I) .EQ. 0 )  GO TO 300                            
                ROOT = I                                                 
C               -------------------------------------------              
C               FIND A SEPARATOR AND NUMBER THE NODES NEXT.              
C               -------------------------------------------              
                CALL  FNDSEP ( ROOT, XADJ, ADJNCY, MASK,                 
     1                         NSEP, PERM(NUM+1), XLS, LS )              
                NUM  = NUM + NSEP                                        
                IF ( NUM .GE. NEQNS )  GO TO 400                         
                GO TO 200                                                
  300     CONTINUE                                                       
C         ----------------------------------------------                 
C         SINCE SEPARATORS FOUND FIRST SHOULD BE ORDERED                 
C         LAST, ROUTINE REVRSE IS CALLED TO ADJUST THE                   
C         ORDERING VECTOR.                                               
C         ----------------------------------------------                 
  400     CALL  REVRSE ( NEQNS, PERM )                                   
          RETURN                                                         
       END                                                               
