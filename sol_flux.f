! This file is particularly tailored for application to MFiX simulations 
! Output is in format [frame#, cell#, voidage] where cell# is proxy for
! spatial coordinates written in order r-, y- and then theta-directions

! Data is only written for cells with (a) voidage > 0.65 and (b) y-
! coordinate < average bed height 
! jbed = integer(jmax*ybed/ylength). Also use jbed and ybed as y-
! variables in geometry.xlsx 

! For MFiX users: compile this file with postmfix! Next, invoke postmfix, enter option 4, then 4 again and then enter tstart, tend and jbed 


      SUBROUTINE SOL_FLUX

      Use param
      Use param1
      Use fldvar
      Use run
      Use geometry
      Use indices
      Use post3d
      Use physprop
      Use compar
	Use functions 

      IMPLICIT NONE
      INCLUDE 'xforms.inc'

      DOUBLE PRECISION  TAVG(DIMENSION_3,7)
      REAL              TIME_REAL(N_SPX)
      REAL              TIME_FOUND, TIME_NOW
      DOUBLE PRECISION  EPGCUTOFF
      INTEGER           FILE1_INDEX , FILE2_INDEX , NSTEP_1
      INTEGER           NTT, CTR, JBED
      INTEGER           REC_POINTER(N_SPX) , L, NT 
      LOGICAL           READ_SPX(N_SPX) , AT_EOF(N_SPX)
      INTEGER           I, J, K, IJK 

      IF (.NOT.DO_XFORMS) THEN
         WRITE (*,'(A,$)')&
                 ' Enter start time, end time & jbed > '
         READ  (*,*) TIME_START, TIME_END, JBED 
         CALL GET_FILE_NAME(TEMP_FILE)
      END IF


    
     EPGCUTOFF = 0.65


      OPEN (UNIT=40,FILE=TEMP_FILE,STATUS='UNKNOWN')

      DO L = 1,N_SPX
         READ_SPX(L)    = .FALSE.
         REC_POINTER(L) = 4
         AT_EOF(L)      = .FALSE.
      END DO
      READ_SPX(1) = .TRUE.         !  EP_g

      CALL SEEK_TIME(READ_SPX, TIME_START, REC_POINTER, TIME_FOUND)
      IF(TIME_FOUND .LT. ZERO) THEN
        WRITE(*,*)' Could not find record for TIME_START'
        RETURN
      ENDIF
      NTT = 0
      NT  = 0

100   CONTINUE
      CALL GET_SAME_TIME (READ_SPX, REC_POINTER, AT_EOF,&
                          TIME_NOW, TIME_REAL, NSTEP_1)
      IF (TIME_NOW .LT. ZERO .OR. TIME_NOW .GT. TIME_END) GOTO 200
      IF (TIME_NOW .LT. TIME_START) GOTO 100
      NTT = NTT + 1
      NT  = NT  +1

      CTR = 0
      DO K = KMIN1, KMAX1
         DO J = JMIN1,JBED+1  
             DO I=IMIN1, IMAX1
                IJK = FUNIJK(I,J,K)
                CTR = CTR + 1
                 IF(EP_g(IJK)>EPGCUTOFF) THEN  
                    WRITE (40,'(2I10, 1p1E15.4)') NTT, CTR, EP_g(IJK)
                 ENDIF
            ENDDO
         ENDDO
      ENDDO

      GOTO 100

200   IF (NT.EQ.0) THEN
         WRITE (*,*) ' No times found in common'
         RETURN
      END IF

      CLOSE (UNIT=40)
      WRITE(*,*)' Number of sampling points = ', NT

      RETURN
      END
