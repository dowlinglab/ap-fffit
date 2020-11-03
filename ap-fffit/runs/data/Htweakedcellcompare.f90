PROGRAM UCC

IMPLICIT NONE

INTEGER, PARAMETER :: ik8 = SELECTED_INT_KIND (8)
INTEGER, PARAMETER :: ik4 = SELECTED_INT_KIND (4)
INTEGER, PARAMETER :: rk8 = SELECTED_REAL_KIND (8)

LOGICAL check
INTEGER (KIND = ik4) :: dt, id, n_atom, n_header, mol, typ, n_atom_pc
INTEGER (KIND = ik8) :: i, j, k, m, l, q, time_start, time_stop
INTEGER (KIND = ik8) :: n_step, counter, n_a, n_b, n_c, n_atom_uc
INTEGER (KIND = ik4), DIMENSION (8) :: PC_Sym
REAL (KIND = rk8) :: dx, dy, dz, vx1, vx2, vy1, vy2, vz1, vz2, v1mag, v2mag, r
REAL (KIND = rk8) :: bx, by, bz, x, y, z, charge, ux, uy, uz, fac
REAL (KIND = rk8) :: x1, x2, y1, y2, z1, z2
REAL (KIND = rk8), DIMENSION (3) :: NH
REAL (KIND = rk8), DIMENSION (40, 3) :: UC_Expt
REAL (KIND = rk8), DIMENSION (40, 3) :: UC_Sim
REAL (KIND = rk8), DIMENSION (8, 3) :: PC_Expt
REAL (KIND = rk8), DIMENSION (8, 3) :: PC_Sim
REAL (KIND = rk8), ALLOCATABLE, DIMENSION (:, :, :) :: Coords !atom info
REAL (KIND = rk8), ALLOCATABLE, DIMENSION (:, :, :) :: Box !box size info
CHARACTER filename1*40, filename2*40
CHARACTER element*2

NH (1) = 1.028
NH (2) = 1.058
NH (3) = 1.031
n_atom_pc = 8
n_atom_uc = 40
n_a = 6
n_b = 9
n_c = 7
n_header = 5 !header lines to skip in trajectory
n_atom = 15120 !number of atoms
time_start = 100000 !starting timestep
time_stop = 200000 !ending timestep
dt = 10000 !frequency of trajectory dump
n_step = (time_stop - time_start) / dt + 1 !snapshots in trajectory
ALLOCATE (Coords (5, n_step, n_atom), Box (2, 3, n_step)) !atom + box info

OPEN (UNIT = 10, FILE = 'in.data', STATUS = 'OLD')

DO i = 1, 41

  READ (10, *)

END DO

DO i = 1, n_atom_uc

  READ (10, *) id, mol, Coords (5, 1, i), charge, UC_Expt (i, 1:3)

END DO

!Hardcoded tweaking of Hydrogen atom coordinates
!Scale N-H lengths to the values listed in the NH array
x1 = UC_Expt (5, 1)
y1 = UC_Expt (5, 2)
z1 = UC_Expt (5, 3)

DO j = 1, 3

  x2 = UC_Expt (5 + j, 1)
  y2 = UC_Expt (5 + j, 2)
  z2 = UC_Expt (5 + j, 3)

  dx = x2 - x1
  dy = y2 - y1
  dz = z2 - z1

  r = SQRT (dx ** 2 + dy ** 2 + dz ** 2)
  fac = NH (j) / r
  dx = dx * fac
  dy = dy * fac
  dz = dz * fac

  UC_Expt (5 + j, 1) = UC_Expt (5, 1) + dx
  UC_Expt (5 + j, 2) = UC_Expt (5, 2) + dy
  UC_Expt (5 + j, 3) = UC_Expt (5, 3) + dz

END DO

x2 = UC_Expt (18, 1)
y2 = UC_Expt (18, 2)
z2 = UC_Expt (18, 3)

dx = x2 - x1
dy = y2 - y1
dz = z2 - z1

r = SQRT (dx ** 2 + dy ** 2 + dz ** 2)
fac = NH (3) / r
dx = dx * fac
dy = dy * fac
dz = dz * fac

UC_Expt (18, 1) = UC_Expt (5, 1) + dx
UC_Expt (18, 2) = UC_Expt (5, 2) + dy
UC_Expt (18, 3) = UC_Expt (5, 3) + dz

x1 = UC_Expt (13, 1)
y1 = UC_Expt (13, 2)
z1 = UC_Expt (13, 3)

DO j = 1, 3

  x2 = UC_Expt (13 + j, 1)
  y2 = UC_Expt (13 + j, 2)
  z2 = UC_Expt (13 + j, 3)

  dx = x2 - x1
  dy = y2 - y1
  dz = z2 - z1

  r = SQRT (dx ** 2 + dy ** 2 + dz ** 2)
  fac = NH (j) / r
  dx = dx * fac
  dy = dy * fac
  dz = dz * fac

  UC_Expt (13 + j, 1) = UC_Expt (13, 1) + dx
  UC_Expt (13 + j, 2) = UC_Expt (13, 2) + dy
  UC_Expt (13 + j, 3) = UC_Expt (13, 3) + dz

END DO

x2 = UC_Expt (40, 1)
y2 = UC_Expt (40, 2)
z2 = UC_Expt (40, 3)

dx = x2 - x1
dy = y2 - y1
dz = z2 - z1

r = SQRT (dx ** 2 + dy ** 2 + dz ** 2)
fac = NH (3) / r
dx = dx * fac
dy = dy * fac
dz = dz * fac

UC_Expt (40, 1) = UC_Expt (13, 1) + dx
UC_Expt (40, 2) = UC_Expt (13, 2) + dy
UC_Expt (40, 3) = UC_Expt (13, 3) + dz

x1 = UC_Expt (23, 1)
y1 = UC_Expt (23, 2)
z1 = UC_Expt (23, 3)

DO j = 1, 3

  x2 = UC_Expt (23 + j, 1)
  y2 = UC_Expt (23 + j, 2)
  z2 = UC_Expt (23 + j, 3)

  dx = x2 - x1
  dy = y2 - y1
  dz = z2 - z1

  r = SQRT (dx ** 2 + dy ** 2 + dz ** 2)
  fac = NH (j) / r
  dx = dx * fac
  dy = dy * fac
  dz = dz * fac

  UC_Expt (23 + j, 1) = UC_Expt (23, 1) + dx
  UC_Expt (23 + j, 2) = UC_Expt (23, 2) + dy
  UC_Expt (23 + j, 3) = UC_Expt (23, 3) + dz

END DO

x2 = UC_Expt (36, 1)
y2 = UC_Expt (36, 2)
z2 = UC_Expt (36, 3)

dx = x2 - x1
dy = y2 - y1
dz = z2 - z1

r = SQRT (dx ** 2 + dy ** 2 + dz ** 2)
fac = NH (3) / r
dx = dx * fac
dy = dy * fac
dz = dz * fac

UC_Expt (36, 1) = UC_Expt (23, 1) + dx
UC_Expt (36, 2) = UC_Expt (23, 2) + dy
UC_Expt (36, 3) = UC_Expt (23, 3) + dz

x1 = UC_Expt (31, 1)
y1 = UC_Expt (31, 2)
z1 = UC_Expt (31, 3)

DO j = 1, 3

  x2 = UC_Expt (31 + j, 1)
  y2 = UC_Expt (31 + j, 2)
  z2 = UC_Expt (31 + j, 3)

  dx = x2 - x1
  dy = y2 - y1
  dz = z2 - z1

  r = SQRT (dx ** 2 + dy ** 2 + dz ** 2)
  fac = NH (j) / r
  dx = dx * fac
  dy = dy * fac
  dz = dz * fac

  UC_Expt (31 + j, 1) = UC_Expt (31, 1) + dx
  UC_Expt (31 + j, 2) = UC_Expt (31, 2) + dy
  UC_Expt (31 + j, 3) = UC_Expt (31, 3) + dz

END DO

x2 = UC_Expt (38, 1)
y2 = UC_Expt (38, 2)
z2 = UC_Expt (38, 3)

dx = x2 - x1
dy = y2 - y1
dz = z2 - z1

r = SQRT (dx ** 2 + dy ** 2 + dz ** 2)
fac = NH (3) / r
dx = dx * fac
dy = dy * fac
dz = dz * fac

UC_Expt (38, 1) = UC_Expt (31, 1) + dx
UC_Expt (38, 2) = UC_Expt (31, 2) + dy
UC_Expt (38, 3) = UC_Expt (31, 3) + dz

!End hardcoded H position tweaking
CLOSE (UNIT = 10, STATUS = 'KEEP')
PC_Expt (1:8, :) = UC_Expt (1:8, :)

OPEN (UNIT = 10, FILE = 'Tweaked_Expt_UC.xyz', STATUS = 'UNKNOWN')
WRITE (10, '(I0)') n_atom_uc
WRITE (10, *)

DO i = 1, n_atom_uc

  SELECT CASE (INT(Coords(5,1,i)))
  CASE(1)
    element = "Cl"
  CASE(2)
    element = "H"
  CASE(3)
    element = "N"
  CASE(4)
    element = "O"
  CASE DEFAULT
    WRITE(*, *) "Invalid atom type: ", Coords(5,1,i),". Exiting."
    CALL EXIT(1)
  END SELECT

  WRITE (10, *) element, UC_Expt (i, 1:3)

END DO

CLOSE (UNIT = 10, STATUS = 'KEEP')
OPEN (UNIT = 10, FILE = 'Tweaked_Expt_PC.xyz', STATUS = 'UNKNOWN')
WRITE (10, '(I0)') n_atom_pc
WRITE (10, *)

DO i = 1, n_atom_pc

  SELECT CASE (INT(Coords(5,1,i)))
  CASE(1)
    element = "Cl"
  CASE(2)
    element = "H"
  CASE(3)
    element = "N"
  CASE(4)
    element = "O"
  CASE DEFAULT
    WRITE(*, *) "Invalid atom type: ", Coords(5,1,i),". Exiting."
    CALL EXIT(1)
  END SELECT

  WRITE (10, *) element, PC_Expt (i, 1:3)

END DO

CLOSE (UNIT = 10, STATUS = 'KEEP')

PC_Sym (1) = 4
PC_Sym (2) = 4
PC_Sym (3) = 4
PC_Sym (4) = 8
PC_Sym (5) = 4
PC_Sym (6) = 4
PC_Sym (7) = 4
PC_Sym (8) = 8

!WRITE (*, *) 'Allocating'

!Initialize some arrays
Coords (:, :, :) = 0.0
Box (:, :, :) = 0.0
UC_Sim (:, :) = 0.0
PC_Sim (:, :) = 0.0

!Open trajectory file
OPEN (UNIT = 10, FILE = 'out.lammpstrj', STATUS = 'OLD')

!WRITE (*, *) 'Reading'

!Loop over snap shots of trajectory file
DO i = 1, n_step

  !Skip over header lines in trajectory file
  DO j = 1, n_header

    READ (10, *)

  ENDDO

  !Read box dimensions into an array
  DO j = 1, 3

    READ (10, *) Box (1, j, i), Box (2, j, i)

  ENDDO

  !Skip line describing elements listed per line in trajectory
  READ (10, *)

  !Loop over total number of atoms, saving per atom information
  DO j = 1, n_atom

    READ (10, *) id, mol, typ, x, y, z
    Coords (1, i, id) = x
    Coords (2, i, id) = y
    Coords (3, i, id) = z
    Coords (4, i, id) = mol
    Coords (5, i, id) = typ

  ENDDO

ENDDO

!WRITE (*, *) 'Done Reading'

CLOSE (UNIT = 10, STATUS = 'KEEP') !Close trajectory file

!Loop over all snapshots
DO i = 1, n_step

!  WRITE (*, '(A, I0)') 'snapshot=', i

  !Calculate box dimensions for snapshot
  bx = Box (2, 1, i) - Box (1, 1, i)
  by = Box (2, 2, i) - Box (1, 2, i)
  bz = Box (2, 3, i) - Box (1, 3, i)

  ux = bx / REAL (n_a) !lattice a
  uy = by / REAL (n_b) !lattice b
  uz = bz / REAL (n_c) !lattice c

  !Loop over all atoms, by unit cells
  DO j = 1, n_a !lattice a repeats

    DO k = 1, n_b !lattice b repeats

      DO m = 1, n_c !lattice c repeats

        DO l = 1, n_atom_uc !atoms per unit cell

          !find atom coordinates and reduce to (0,0,0), (ux, uy, uz)
          q = n_atom_uc * ((j - 1) * n_b * n_c + (k - 1) * n_c + (m - 1)) + l 
          x = Coords (1, i, q) - Box (1, 1, i)
          y = Coords (2, i, q) - Box (1, 2, i)
          z = Coords (3, i, q) - Box (1, 3, i)

          DO WHILE (x .LT. 0.0)

            x = x + ux

          END DO

          DO WHILE (x .GT. ux)

            x = x - ux

          END DO

          DO WHILE (y .LT. 0.0)

            y = y + uy

          END DO

          DO WHILE (y .GT. uy)

            y = y - uy

          END DO

          DO WHILE (z .LT. 0.0)

            z = z + uz

          END DO

          DO WHILE (z .GT. uz)

            z = z - uz

          END DO

          UC_Sim (l, 1) = UC_Sim (l, 1) + x
          UC_Sim (l, 2) = UC_Sim (l, 2) + y
          UC_Sim (l, 3) = UC_Sim (l, 3) + z

          CALL SYM (l)
          PC_Sim (q, 1) = PC_Sim (q, 1) + x
          PC_Sim (q, 2) = PC_Sim (q, 2) + y
          PC_Sim (q, 3) = PC_Sim (q, 3) + z

        END DO

      END DO

    END DO

  END DO

END DO

UC_Sim (:, :) = UC_Sim (:, :) / REAL (n_step * n_a * n_b * n_c)
PC_Sim (:, :) = PC_Sim (:, :) / REAL (n_step * n_a * n_b * n_c)
PC_Sim (:, 1) = PC_Sim (:, 1) / REAL (PC_Sym (:))
PC_Sim (:, 2) = PC_Sim (:, 2) / REAL (PC_Sym (:))
PC_Sim (:, 3) = PC_Sim (:, 3) / REAL (PC_Sym (:))

OPEN (UNIT = 20, FILE = 'Sim_UC.xyz', STATUS = 'UNKNOWN')
OPEN (UNIT = 30, FILE = 'Residuals_UC.txt', STATUS = 'UNKNOWN')
OPEN (UNIT = 40, FILE = 'Sim_PC.xyz', STATUS = 'UNKNOWN')
OPEN (UNIT = 50, FILE = 'Residuals_PC.txt', STATUS = 'UNKNOWN')

WRITE (20, *) n_atom_uc
WRITE (20, *)
WRITE (40, *) n_atom_pc
WRITE (40, *)

DO i = 1, n_atom_uc
  SELECT CASE (INT(Coords(5,1,i)))
  CASE(1)
    element = "Cl"
  CASE(2)
    element = "H"
  CASE(3)
    element = "N"
  CASE(4)
    element = "O"
  CASE DEFAULT
    WRITE(*, *) "Invalid atom type: ", Coords(5,1,i),". Exiting."
    CALL EXIT(1)
  END SELECT

  WRITE (20, *) element, UC_Sim (i, :)
  WRITE (30, *) element, (UC_Expt (i, :) - UC_Sim (i, :)) ** 2

END DO

DO i = 1, n_atom_pc
  SELECT CASE (INT(Coords(5,1,i)))
  CASE(1)
    element = "Cl"
  CASE(2)
    element = "H"
  CASE(3)
    element = "N"
  CASE(4)
    element = "O"
  CASE DEFAULT
    WRITE(*, *) "Invalid atom type: ", Coords(5,1,i),". Exiting."
    CALL EXIT(1)
  END SELECT

  WRITE (40, *) element, PC_Sim (i, :)
  WRITE (50, *) element, (PC_Expt (i, :) - PC_Sim (i, :)) ** 2

END DO

CONTAINS

SUBROUTINE SYM (a)

  IMPLICIT NONE

  INTEGER, INTENT (IN) :: a

  IF (a .EQ. 1) THEN

      x = x
      y = y
      z = z
      q = 1

    ELSE IF (a .EQ. 2) THEN

      x = x
      y = y
      z = z
      q = 2

    ELSE IF (a .EQ. 3) THEN

      x = x
      y = y
      z = z
      q = 3

    ELSE IF (a .EQ. 4) THEN

      x = x
      y = y
      z = z
      q = 4

    ELSE IF (a .EQ. 5) THEN

      x = x
      y = y
      z = z
      q = 5

    ELSE IF (a .EQ. 6) THEN

      x = x
      y = y
      z = z
      q = 6

    ELSE IF (a .EQ. 7) THEN

      x = x
      y = y
      z = z
      q = 7

    ELSE IF (a .EQ. 8) THEN

      x = x
      y = y
      z = z
      q = 8

    ELSE IF (a .EQ. 9) THEN

      x = 0.5 * ux - x
      y = y - 0.5 * uy
      z = z - 0.5 * uz
      q = 1

    ELSE IF (a .EQ. 10) THEN

      x = 0.5 * ux - x
      y = y - 0.5 * uy
      z = z - 0.5 * uz
      q = 2

    ELSE IF (a .EQ. 11) THEN

      x = x - ux
      x = 0.5 * ux - x
      y = y - 0.5 * uy
      z = z - 0.5 * uz
      q = 3

    ELSE IF (a .EQ. 12) THEN

      x = 0.5 * ux - x
      y = y - 0.5 * uy
      z = z - 0.5 * uz
      q = 4

    ELSE IF (a .EQ. 13) THEN

      x = 0.5 * ux - x
      y = y - 0.5 * uy
      z = z + uz
      z = z - 0.5 * uz
      q = 5

    ELSE IF (a .EQ. 14) THEN

      x = 0.5 * ux - x
      y = y - 0.5 * uy
      z = z + uz
      z = z - 0.5 * uz
      q = 6

    ELSE IF (a .EQ. 15) THEN

      x = 0.5 * ux - x
      y = y - 0.5 * uy
      z = z + uz
      z = z - 0.5 * uz
      q = 7

    ELSE IF (a .EQ. 16) THEN

      x = 0.5 * ux - x
      y = y - 0.5 * uy
      z = z + uz
      z = z - 0.5 * uz
      q = 8

    ELSE IF (a .EQ. 17) THEN

      x = x
      y = 0.5 * uy - y
      z = z
      q = 4

    ELSE IF (a .EQ. 18) THEN

      x = x
      y = 0.5 * uy - y
      z = z
      q = 8

    ELSE IF (a .EQ. 19) THEN

      x = x - 0.5 * ux
      y = y
      z = 0.5 * uz - z
      q = 1

    ELSE IF (a .EQ. 20) THEN

      x = x - 0.5 * ux
      y = y
      z = 0.5 * uz - z
      q = 2

    ELSE IF (a .EQ. 21) THEN

      x = x + ux
      x = x - 0.5 * ux
      y = y
      z = 0.5 * uz - z
      q = 3

    ELSE IF (a .EQ. 22) THEN

      x = x - 0.5 * ux
      y = y
      z = 0.5 * uz - z
      q = 4

    ELSE IF (a .EQ. 23) THEN

      x = x - 0.5 * ux
      y = y
      z = z - uz
      z = 0.5 * uz - z
      q = 5

    ELSE IF (a .EQ. 24) THEN

      x = x - 0.5 * ux
      y = y
      z = z - uz
      z = 0.5 * uz - z
      q = 6

    ELSE IF (a .EQ. 25) THEN

      x = x - 0.5 * ux
      y = y
      z = z - uz
      z = 0.5 * uz - z
      q = 7

    ELSE IF (a .EQ. 26) THEN

      x = x - 0.5 * ux
      y = y
      z = z - uz
      z = 0.5 * uz - z
      q = 8

    ELSE IF (a .EQ. 27) THEN

      x = x - ux
      x = -x
      y = y - uy
      y = -y
      z = z - uz
      z = -z
      q = 1

    ELSE IF (a .EQ. 28) THEN

      x = x - ux
      x = -x
      y = y - uy
      y = -y
      z = z - uz
      z = -z
      q = 2

    ELSE IF (a .EQ. 29) THEN

      x = x - ux
      x = -x
      y = y - uy
      y = -y
      z = z - uz
      z = -z
      q = 3

    ELSE IF (a .EQ. 30) THEN

      x = x - ux
      x = -x
      y = y - uy
      y = -y
      z = z - uz
      z = -z
      q = 4

    ELSE IF (a .EQ. 31) THEN

      x = x - ux
      x = -x
      y = y - uy
      y = -y
      z = z - uz
      z = -z
      q = 5

    ELSE IF (a .EQ. 32) THEN

      x = x - ux
      x = -x
      y = y - uy
      y = -y
      z = z - uz
      z = -z
      q = 6

    ELSE IF (a .EQ. 33) THEN

      x = x - ux
      x = -x
      y = y - uy
      y = -y
      z = z - uz
      z = -z
      q = 7

    ELSE IF (a .EQ. 34) THEN

      x = x - ux
      x = -x
      y = y - uy
      y = -y
      z = z - uz
      z = -z
      q = 8

    ELSE IF (a .EQ. 35) THEN

      x = x - 0.5 * ux
      y = 0.5 * uy - y
      z = 0.5 * uz - z
      q = 4

    ELSE IF (a .EQ. 36) THEN

      x = x - 0.5 * ux
      y = 0.5 * uy - y
      z = z - uz
      z = 0.5 * uz - z
      q = 8

    ELSE IF (a .EQ. 37) THEN

      x = x - ux
      x = -x
      y = y - 0.5 * uy
      z = z - uz
      z = -z
      q = 4

    ELSE IF (a .EQ. 38) THEN

      x = x - ux
      x = -x
      y = y - 0.5 * uy
      z = z - uz
      z = -z
      q = 8

    ELSE IF (a .EQ. 39) THEN

      x = 0.5 * ux - x
      y = y - uy
      y = -y
      z = z - 0.5 * uz
      q = 4

    ELSE IF (a .EQ. 40) THEN

      x = 0.5 * ux - x
      y = y - uy
      y = -y
      z = z + uz
      z = z - 0.5 * uz
      q = 8

  END IF

END SUBROUTINE SYM

END PROGRAM UCC
