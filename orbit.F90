PROGRAM ORBIT
	USE omp_lib
	USE BSINT
	USE PHYSICS
	IMPLICIT NONE

	! misc variables
	INTEGER :: ii, ij, ik, id
	DOUBLE PRECISION :: ecc

	! number of systems to evolve simultaneously
	INTEGER :: nsystems
	PARAMETER(nsystems=10)
	
	! body parameters
	! mass in solar mass
	! position in AU
	! velocity in AU/yr
	DOUBLE PRECISION, DIMENSION(nsystems) :: mass_a, mass_b, mass_c
	DOUBLE PRECISION, DIMENSION(3,nsystems) :: pos_a, pos_b, pos_c
	DOUBLE PRECISION, DIMENSION(3,nsystems) :: vel_a, vel_b, vel_c

	! integration parameters
	! time in years
	DOUBLE PRECISION, DIMENSION(nsystems) :: currtime, dt, endtime
	! 0 = ok, 1 = escaped, 2 = collided
	INTEGER, DIMENSION(nsystems) :: stat

	WRITE(*,'(A)') "Run Parameters:"
	!$OMP PARALLEL
	id = OMP_GET_THREAD_NUM()
	IF (id .EQ. 0) THEN
		WRITE(*,'(A,I4)') "  Number of threads = ",OMP_GET_NUM_THREADS()
	ENDIF
	!$OMP END PARALLEL
	WRITE(*,'(A,I7,A)') "  ",nsystems," Systems"

	! initialize a bunch of systems
	DO ii=1,nsystems
		mass_a(ii) = 1D0 ! sun
		pos_a(1,ii) = 0D0
		pos_a(2,ii) = 0D0
		pos_a(3,ii) = 0D0
		vel_a(1,ii) = 0D0
		vel_a(2,ii) = 0D0
		vel_a(3,ii) = 0D0

		mass_b(ii) = 1D-3 ! driver
		pos_b(1,ii) = 4.5D0
		pos_b(2,ii) = 0D0
		pos_b(3,ii) = 2.6D0
		vel_b(1,ii) = 0D0
		vel_b(2,ii) = 2.7D0
		vel_b(3,ii) = 0D0

		mass_c(ii) = 1D-6 ! test particle
		pos_c(1,ii) = 0D0
		pos_c(2,ii) = 1D0
		pos_c(3,ii) = 0D0
		vel_c(1,ii) = 6.32D0
		vel_c(2,ii) = 0D0
		vel_c(3,ii) = 0D0

		currtime(ii) = 0D0
		dt(ii) = 0.0001D0
		endtime(ii) = 1e5
		stat(ii) = 0
	ENDDO

	WRITE(*,'(A)') "Systems Initialized, Beginning Integration."

	! evolve the systems in parallel
	!$OMP PARALLEL DO PRIVATE(ii,ij) SCHEDULE(DYNAMIC,1)
	DO ii=1,nsystems
		ij = 0
		DO WHILE (currtime(ii) .LT. endtime(ii) .AND. stat(ii) .EQ. 0)

			! decide timestep?

			CALL DriverPos(pos_b(:,ii),currtime(ii),11.86D0)
			! integrate
			CALL Euler2(pos_a(:,ii), vel_a(:,ii), mass_a(ii), &
					pos_b(:,ii), vel_b(:,ii), mass_b(ii), &
					pos_c(:,ii), vel_c(:,ii), mass_c(ii), dt(ii))

			! check for escape
			IF (Escaping(pos_c(:,ii), vel_c(:,ii), mass_c(ii), &
					pos_a(:,ii), vel_a(:,ii), mass_a(ii))) THEN
				WRITE(*,'(A,I6,A)') "--System ",ii," Escaped"
				stat(ii) = 1
			ENDIF

			! check for collision
			IF (HasCollided(pos_c(:,ii), pos_a(:,ii))) THEN
				stat(ii) = 2
			ENDIF

			IF (ii .EQ. 1 .AND. ij .EQ. 100000) THEN
				CALL Eccentricity(pos_a(:,ii), vel_a(:,ii), mass_a(ii), &
					pos_c(:,ii), vel_c(:,ii), mass_c(ii), ecc)
				!WRITE(*,'(E15.6,7E15.6)') currtime(ii), pos_b(:,ii), pos_c(:,ii), ecc
				WRITE(*,'(2E15.6)') currtime(ii), ecc
				ij = 0
			ENDIF
			currtime(ii) = currtime(ii) + dt(ii)
			ij = ij + 1
		ENDDO
	ENDDO
	!$OMP END PARALLEL DO

	WRITE(*,'(A)') "Integration Complete"

	! output stats
	

END PROGRAM

