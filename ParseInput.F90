
! Contains global variables, settings, etc

MODULE ParseInput
	USE omp_lib

	IMPLICIT NONE

	INTEGER :: nsystems
	LOGICAL :: verbose
	DOUBLE PRECISION :: driver_mass_min, driver_mass_max, test_mass_min, test_mass_max
	DOUBLE PRECISION :: driver_sma_min, driver_sma_max, test_sma_min, test_sma_max
	DOUBLE PRECISION :: driver_inc_min, driver_inc_max, test_inc_min, test_inc_max
	DOUBLE PRECISION :: driver_ecc_min, driver_ecc_max, test_ecc_min, test_ecc_max

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

			! number of systems
			IF (strbuffer .EQ. "-n") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) nsystems
				IF (nsystems .LT. 1) THEN
					PRINT*, "ERROR: nsystems must be > 0  : ", nsystems
					STOP
				ENDIF
			ENDIF

			! driver/test mass range
			IF (strbuffer .EQ. "--driver_mass_max") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) driver_mass_max
			ENDIF
			IF (strbuffer .EQ. "--driver_mass_min") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) driver_mass_min
			ENDIF
			IF (strbuffer .EQ. "--test_mass_max") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) test_mass_max
			ENDIF
			IF (strbuffer .EQ. "--test_mass_min") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) test_mass_min
			ENDIF

			! semi-major axes ranges
			IF (strbuffer .EQ. "--driver_sma_max") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) driver_sma_max
			ENDIF
			IF (strbuffer .EQ. "--driver_sma_min") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) driver_sma_min
			ENDIF
			IF (strbuffer .EQ. "--test_sma_max") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) test_sma_max
			ENDIF
			IF (strbuffer .EQ. "--test_sma_min") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) test_sma_min
			ENDIF

			! inclination ranges
			IF (strbuffer .EQ. "--driver_inc_max") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) driver_inc_max
			ENDIF
			IF (strbuffer .EQ. "--driver_inc_min") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) driver_inc_min
			ENDIF
			IF (strbuffer .EQ. "--test_inc_max") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) test_inc_max
			ENDIF
			IF (strbuffer .EQ. "--test_inc_min") THEN
				CALL getarg(ii+1,strbuffer)
				READ(strbuffer,*) test_inc_min
			ENDIF

		ENDDO
	END SUBROUTINE ReadCommandLine



	! Set the default values of various run parameters
	! It is best to call this before reading in any conflicting information..
	SUBROUTINE SetDefaults()
		verbose = .FALSE.
		nsystems = 1
		driver_mass_min = 1D-3
		driver_mass_max = 1D-3
		test_mass_min = 1D-6
		test_mass_max = 1D-6
		driver_sma_min = 1D2
		driver_sma_max = 1D2
		test_sma_min = 1D0
		test_sma_max = 1D0
		driver_inc_min = 1D0
		driver_inc_max = 1D0
		test_inc_min = 0D0
		test_inc_max = 0D0
		driver_ecc_min = 0D0
		driver_ecc_max = 0D0
		test_ecc_min = 5D-1
		test_ecc_max = 5D-1
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
		PRINT*, "Minimal Usage: ./orbit -n #"
		PRINT*, "Other Options:"
		PRINT*, " -v      Verbose Mode"
		PRINT*, " -n #    Number of Random Systems"
	END SUBROUTINE Usage

END MODULE
