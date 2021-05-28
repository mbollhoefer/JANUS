#include <ilupack_fortran.h>
C----- SUBROUTINE FNDSEP
C***************************************************************         
C***************************************************************         
C************     FNDSEP ..... FIND SEPARATOR       ************         
C***************************************************************         
C***************************************************************         
C                                                                        
C     PURPOSE - THIS ROUTINE IS USED TO FIND A SMALL                     
C               SEPARATOR FOR A CONNECTED COMPONENT SPECIFIED            
C               BY MASK IN THE GIVEN GRAPH.                              
C                                                                        
C     INPUT PARAMETERS -                                                 
C        ROOT - IS THE NODE THAT DETERMINES THE MASKED                   
C               COMPONENT.                                               
C        (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE PAIR.                  
C                                                                        
C     OUTPUT PARAMETERS -                                                
C        NSEP - NUMBER OF VARIABLES IN THE SEPARATOR.                    
C        SEP - VECTOR CONTAINING THE SEPARATOR NODES.                    
C                                                                        
C     UPDATED PARAMETER -                                                
C        MASK - NODES IN THE SEPARATOR HAVE THEIR MASK                   
C               VALUES SET TO ZERO.                                      
C                                                                        
C     WORKING PARAMETERS -                                               
C        (XLS, LS) - LEVEL STRUCTURE PAIR FOR LEVEL STRUCTURE            
C               FOUND BY FNROOT.                                         
C                                                                        
C     PROGRAM SUBROUTINES -                                              
C        FNROOT.                                                         
C                                                                        
C***************************************************************         
C                                                                        
      SUBROUTINE  FNDSEP ( ROOT, XADJ, ADJNCY, MASK,                     
     1                     NSEP, SEP, XLS , LS )                         
C                                                                        
C***************************************************************         
C                                                                        
         integer ADJNCY(1), LS(1), MASK(1), SEP(1), XLS(1)               
         integer XADJ(1), I, J, JSTOP, JSTRT, MIDBEG,                    
     1           MIDEND, MIDLVL, MP1BEG, MP1END,                         
     1           NBR, NLVL, NODE, NSEP, ROOT                             
C                                                                        
C***************************************************************         
C                                                                        
         CALL  FNROOT ( ROOT, XADJ, ADJNCY, MASK,                        
     1                  NLVL, XLS, LS )                                  
C        ----------------------------------------------                  
C        IF THE NUMBER OF LEVELS IS LESS THAN 3, RETURN                  
C        THE WHOLE COMPONENT AS THE SEPARATOR.                           
C        ----------------------------------------------                  
         IF ( NLVL .GE. 3 )  GO TO 200                                   
            NSEP = XLS(NLVL+1) - 1                                       
            DO 100 I = 1, NSEP                                           
               NODE = LS(I)                                              
               SEP(I) = NODE                                             
               MASK(NODE) = 0                                            
  100       CONTINUE                                                     
            RETURN                                                       
C        ----------------------------------------------------            
C        FIND THE MIDDLE LEVEL OF THE ROOTED LEVEL STRUCTURE.            
C        ----------------------------------------------------            
  200    MIDLVL = (NLVL + 2)/2                                           
         MIDBEG = XLS(MIDLVL)                                            
         MP1BEG = XLS(MIDLVL + 1)                                        
         MIDEND = MP1BEG - 1                                             
         MP1END = XLS(MIDLVL+2) - 1                                      
C        -------------------------------------------------               
C        THE SEPARATOR IS OBTAINED BY INCLUDING ONLY THOSE               
C        MIDDLE-LEVEL NODES WITH NEIGHBORS IN THE MIDDLE+1               
C        LEVEL. XADJ IS USED TEMPORARILY TO MARK THOSE                   
C        NODES IN THE MIDDLE+1 LEVEL.                                    
C        -------------------------------------------------               
         DO 300 I = MP1BEG, MP1END                                       
            NODE = LS(I)                                                 
            XADJ(NODE) = - XADJ(NODE)                                    
  300    CONTINUE                                                        
         NSEP  = 0                                                       
         DO 500 I = MIDBEG, MIDEND                                       
            NODE = LS(I)                                                 
            JSTRT = XADJ(NODE)                                           
            JSTOP = IABS(XADJ(NODE+1)) - 1                               
            DO 400 J = JSTRT, JSTOP                                      
               NBR = ADJNCY(J)                                           
               IF ( XADJ(NBR) .GT. 0 )  GO TO 400                        
                  NSEP = NSEP + 1                                        
                  SEP(NSEP) = NODE                                       
                  MASK(NODE) = 0                                         
                  GO TO 500                                              
  400       CONTINUE                                                     
  500    CONTINUE                                                        
C        -------------------------------                                 
C        RESET XADJ TO ITS CORRECT SIGN.                                 
C        -------------------------------                                 
         DO 600 I = MP1BEG, MP1END                                       
            NODE = LS(I)                                                 
            XADJ(NODE) = - XADJ(NODE)                                    
  600    CONTINUE                                                        
         RETURN                                                          
      END                                                                
