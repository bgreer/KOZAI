
! Contains global variables, settings, etc

MODULE ParseInput
	USE omp_lib

	IMPLICIT NONE

	INTEGER :: nsystems
	PARAMETER(nsystems=100)
	LOGICAL :: verbose

CONTAINS

	! Note the lack of IMPLICIT NONE
	! Poor coding? Maybe.
	SUBROUTINE ReadCommandLine()
		INTEGER :: ii, argcount, currarg, tempint, tempint2
		CHARACTER(LEN=200) :: strbuffer

		argcount = IARGC()

		! Read command line arguments
		DO ii=1,argcount
			CALL getarg(ii,strbuffer)

			! help
			IF (strbuffer .EQ. "--help") THEN
				CALL Usage()
				STOP
			ENDIF

			! verbose mode
			IF (strbuffer .EQ. "-v") THEN
				verbose = .TRUE.
			ENDIF

		ENDDO
	END SUBROUTINE ReadCommandLine



	! Set the default values of various run parameters
	! It is best to call this before reading in any conflicting information..
	SUBROUTINE SetDefaults()
		verbose = .FALSE.
	END SUBROUTINE SetDefaults

	! Print the run details to stdout
	SUBROUTINE PrintDetails()
		INTEGER :: id
		WRITE(*,'(A)') "Run Parameters:"
		!$OMP PARALLEL
		id = OMP_GET_THREAD_NUM()
		IF (id .EQ. 0) THEN
			WRITE(*,'(A,I4)') "  Number of threads = ",OMP_GET_NUM_THREADS()
		ENDIF
		!$OMP END PARALLEL
		WRITE(*,'(A,I7,A)') "  ",nsystems," Systems"
	END SUBROUTINE PrintDetails

	SUBROUTINE Usage()
		PRINT*, "Minimal Usage: ./orbit"
		PRINT*, "Other Options:"
		PRINT*, " -v      Verbose Mode"
		PRINT*, " -l [#]  Number of timesteps to track for"
		PRINT*, " -r [#]  Tracking rate (0=car, 1=snod, 2=custom)"
		PRINT*, " -ml [file]  Master list of dopplergrams"
		PRINT*, " -and others"
	END SUBROUTINE Usage

END MODULE
