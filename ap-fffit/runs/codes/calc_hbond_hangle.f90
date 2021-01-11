PROGRAM HBOND

IMPLICIT NONE

INTEGER, PARAMETER :: ik8 = SELECTED_INT_KIND (8)
INTEGER, PARAMETER :: ik4 = SELECTED_INT_KIND (4)
INTEGER, PARAMETER :: rk8 = SELECTED_REAL_KIND (8)
INTEGER, PARAMETER :: prec = SELECTED_REAL_KIND (8) !for the sorting routine

LOGICAL check
INTEGER (KIND = ik4) :: dt, id, n_ion, n_atom, n_header, mol, typ, nn
INTEGER (KIND = ik8) :: i, j, k, m, l, time_start, time_stop, mino
INTEGER (KIND = ik8) :: n_step, counter, ind, a_count, a_end
INTEGER (KIND = ik8) :: b_start, a_start, b_end
INTEGER (KIND = ik4), ALLOCATABLE, DIMENSION (:, :) :: Neigh_O !O neigh list
INTEGER (KIND = ik4), ALLOCATABLE, DIMENSION (:, :) :: NH4 !ammonium atom list
REAL (KIND = rk8) :: db, da, max_bond, min_bond, max_angle, min_angle, dist
REAL (KIND = rk8) :: dx, dy, dz, vx1, vx2, vy1, vy2, vz1, vz2, v1mag, v2mag
REAL (KIND = rk8) :: dot, angle, bx, by, bz, x, y, z, x1, y1, z1, x2, y2, z2
REAL (KIND = rk8) :: mindist, x3, y3, z3
REAL (KIND = rk8), DIMENSION (4) :: Bavg !Average hydrogen bond values
REAL (KIND = rk8), DIMENSION (4) :: Aavg !Average hydrogen bond angle values
REAL (KIND = rk8), DIMENSION (4) :: Bavgsort !sorted avg bond values
REAL (KIND = rk8), DIMENSION (4) :: Aavgsort !sorted avg hbond values
REAL (KIND = rk8), DIMENSION (4) :: Btemp1 !temp bond values
REAL (KIND = rk8), DIMENSION (4) :: Btemp2 !sorted temp bond values
REAL (KIND = rk8), DIMENSION (4) :: Atemp1 !temp angle values
REAL (KIND = rk8), DIMENSION (4) :: Atemp2 !sorted temp angle values
REAL (KIND = rk8), ALLOCATABLE, DIMENSION (:) :: Olist !sorted oxygen list
REAL (KIND = rk8), ALLOCATABLE, DIMENSION (:, :) :: Bonds !bond histograms
REAL (KIND = rk8), ALLOCATABLE, DIMENSION (:, :) :: Angles !angle histograms
REAL (KIND = rk8), ALLOCATABLE, DIMENSION (:, :, :) :: O_Dist !N-O distances
REAL (KIND = rk8), ALLOCATABLE, DIMENSION (:, :, :) :: Coords !atom info
REAL (KIND = rk8), ALLOCATABLE, DIMENSION (:, :, :) :: Box !box size info
CHARACTER filename1*40, filename2*40

nn = 10 !number of nearest oxygens for neighborlist
n_header = 5 !header lines to skip in trajectory
n_atom = 15120 !number of atoms
n_ion = 1512 !number of each type of ion
time_start = 100000 !starting timestep
time_stop = 200000 !ending timestep
dt = 10000 !frequency of trajectory dump
n_step = (time_stop - time_start) / dt + 1 !snapshots in trajectory
db = 0.00001 !resolution for bond histogram
da = 0.001 !resolution for angle histogram
max_bond = 5.0 !maximum hydrogen bond length for histogram - Changed from 3.5 on
!10/15
min_bond = 0.10 !minimum hydrogen bond length for histogram - Changed from 0.3
!on 10/15
max_angle = 180.0 !maximum hydrogen bond angle for histogram
min_angle = 1.0 !minimum hydrogen bond angle for histogram - Changed from 30 on
!10/15
b_start = INT (min_bond / db) - 1 !lower index for bond histogram
b_end = INT (max_bond / db) + 1 !upper index for bond histogram
a_start = INT (min_angle / da) - 1 !lower index for angle histogram
a_end = INT (max_angle / da) + 1 !upper index for angle histogram

!WRITE (*, *) 'Allocating'

ALLOCATE (Olist (n_ion * 4)) !total number of O atoms per snapshot
ALLOCATE (Neigh_O (n_ion, 10)) !number of O neigh per NH4
ALLOCATE (NH4 (n_ion, 5)) !number of atoms per NH4
ALLOCATE (Coords (5, n_step, n_atom), Box (2, 3, n_step)) !atom + box info
ALLOCATE (Bonds (4, b_start:b_end), Angles (4, a_start:b_end)) !histogramming
ALLOCATE (O_Dist (n_ion, 2, n_ion * 4)) !for each NH4, all O-N dist + O atom id

!Initialize some arrays
NH4 (:, :) = 0
Neigh_O (:, :) = 0
Coords (:, :, :) = 0.0
Box (:, :, :) = 0.0
O_Dist (:, :, :) = 0.0

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

!Create separate list for NH4 atom ids, sorted by mol id for easy reference
!Loop over total number of atoms in first snapshot
!NH4 array = (1:n_ion, 5) where (:, 1) are N and (:, 2:5) are H
DO j = 1, n_atom

  !Find Nitrogen atoms, type = 3
  IF (Coords (5, 1, j) .EQ. 3) THEN

      NH4 (Coords (4, 1, j) - n_ion, 1) = j

    !Find Hydrogen atoms, type = 2
    ELSE IF (Coords (5, 1, j) .EQ. 2) THEN

      IF (NH4 (Coords (4, 1, j) - n_ion, 2) .EQ. 0) THEN

          NH4 (Coords (4, 1, j) - n_ion, 2) = j

        ELSE IF (NH4 (Coords (4, 1, j) - n_ion, 3) .EQ. 0) THEN

          NH4 (Coords (4, 1, j) - n_ion, 3) = j

        ELSE IF (NH4 (Coords (4, 1, j) - n_ion, 4) .EQ. 0) THEN

          NH4 (Coords (4, 1, j) - n_ion, 4) = j

        ELSE

          NH4 (Coords (4, 1, j) - n_ion, 5) = j

      END IF

    ELSE

      CYCLE

  END IF

END DO

!WRITE (*, *) 'Assigned NH4 array'

!Calculate box lengths for first snapshot and initialize counter
bx = Box (2, 1, 1) - Box (1, 1, 1)
by = Box (2, 2, 1) - Box (1, 2, 1)
bz = Box (2, 3, 1) - Box (1, 3, 1)
counter = 0 !keeps track of number of oxygen atoms found per NH4

!Calculate and save all N-O distances
!Loop over all atoms
DO i = 1, n_atom

  !Search for Oxygen, type = 4, load coordinates, increment oxygen counter
  IF (Coords (5, 1, i) .EQ. 4) THEN

      counter = counter + 1
      x1 = Coords (1, 1, i)
      y1 = Coords (2, 1, i)
      z1 = Coords (3, 1, i)

      !Loop over all Nitrogen, calculate O-N distance, save dist and O atom id
      DO j = 1, n_ion

        x2 = Coords (1, 1, NH4 (j, 1))
        y2 = Coords (2, 1, NH4 (j, 1))
        z2 = Coords (3, 1, NH4 (j, 1))
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1
        dx = dx - bx * ANINT (dx / bx)
        dy = dy - by * ANINT (dy / by)
        dz = dz - bz * ANINT (dz / bz)
        O_Dist (j, 1, counter) = i !saving oxygen atom id
        O_Dist (j, 2, counter) = dx ** 2 + dy ** 2 + dz ** 2 !saving O distance

      END DO

    ELSE

      CYCLE

  END IF

END DO

!WRITE (*, *) 'Calculated Oxygen distances to Nitrogens'

!Sort O-N distances and find associated O atom ids
!Only save desired number of nearest neighboring O atoms "nn"
!Loop over all NH4 ions
DO i = 1, n_ion

  !Create temporary list to send to sorting subroutine
  Olist (:) = O_Dist (i, 2, :)
  CALL dpquicksort (Olist)

  !Only save O atom ids of the "nn" quantity of shortest N-O distances
  DO j = 1, nn

    !Initialize check and k
    check = .FALSE. !if TRUE, break out of do while loop
    k = 1 !atom id index, incremented in the do while loop

    !Search for oxygen atom corresponding to the jth shortest O-N distance
    DO WHILE (check .EQ. .FALSE.)

      !If jth shortest O-N distance found, save O atom id in neighborlist
      IF (Olist (j) .EQ. O_Dist (i, 2, k)) THEN

          check = .TRUE. !breaks out of do while loop
          Neigh_O (i, j) = O_dist (i, 1, k) !i=ion, j=O neigh, (i, 1, k)=O id

        ELSE

          !increment to next O-N distance comparison
          k = k + 1

      END IF

    END DO

  END DO

END DO

!WRITE (*, *) 'Created Neighborlists'

!No longer need temp Olist or full list of O-N distances of first snapshot
DEALLOCATE (Olist, O_Dist)
CLOSE (UNIT = 10, STATUS = 'KEEP') !Close trajectory file
OPEN (UNIT = 40, FILE = 'HBondstats.txt', STATUS = 'UNKNOWN') !Bond averages
OPEN (UNIT = 50, FILE = 'HAnglestats.txt', STATUS = 'UNKNOWN') !Angle averages

!WRITE (*, *) 'Starting Bond and Angle Calculations'

!Initialize arrays
Bonds (:, :) = 0.0
Angles (:, :) = 0.0
Bavg (:) = 0.0
Aavg (:) = 0.0
Bavgsort (:) = 0.0
Aavgsort (:) = 0.0

!Loop over all snapshots
DO i = 1, n_step

!  WRITE (*, '(A, I0)') 'snapshot=', i

  !Calculate box dimensions for snapshot
  bx = Box (2, 1, i) - Box (1, 1, i)
  by = Box (2, 2, i) - Box (1, 2, i)
  bz = Box (2, 3, i) - Box (1, 3, i)

  !Loop over all NH4 ions
  DO j = 1, n_ion

    !Nitrogen coordinates of the jth NH4 ion
    x1 = Coords (1, i, NH4 (j, 1))
    y1 = Coords (2, i, NH4 (j, 1))
    z1 = Coords (3, i, NH4 (j, 1))

    !Loop over the hydrogen atoms of the jth NH4 ion
    DO k = 2, 5

      !kth Hydrogen coordinates of the jth NH4 ion
      x2 = Coords (1, i, NH4 (j, k))
      y2 = Coords (2, i, NH4 (j, k))
      z2 = Coords (3, i, NH4 (j, k))
      mindist = 10.0 !initialize min. hydrogen bond dist (arb. large value)

      !Loop over the "nn" nearest oxygen neighbors of the jth NH4 ion
      !Attempt to find minimum O-H distance possible from provided neighlist
      DO l = 1, nn

        !Oxygen coordinates of nnth Oxygen in neighborlist for the jth NH4 ion
        !Calculate O-H distance
        x3 = Coords (1, i, Neigh_O (j, l))
        y3 = Coords (2, i, Neigh_O (j, l))
        z3 = Coords (3, i, Neigh_O (j, l))
        dx = x3 - x2
        dy = y3 - y2
        dz = z3 - z2
        dx = dx - bx * ANINT (dx / bx)
        dy = dy - by * ANINT (dy / by)
        dz = dz - bz * ANINT (dz / bz)
        dist = SQRT (dx ** 2 + dy ** 2 + dz ** 2)

        !Keep the information of the smallest O-H distance
        IF (dist .LT. mindist) THEN

            mindist = dist !update minimum O-H distance for kth H of jth NH4
            mino = l !save associated oxygen in neighborlist

          ELSE

            CYCLE

        END IF

      END DO

      !Add hydrogen bond information to histogram and average arrays
      ind = INT (mindist / db) + 1 !calculate index for bond histogram
      Bonds (k - 1, ind) = Bonds (k - 1, ind) + 1 !increment bond histogram
      Btemp1 (k - 1) = mindist !save hbond dist to temp array to be sorted

      !Oxygen coordinates of hydrogen bond + calculate angle
      x3 = Coords (1, i, Neigh_O (j, mino))
      y3 = Coords (2, i, Neigh_O (j, mino))
      z3 = Coords (3, i, Neigh_O (j, mino))
      vx1 = x2 - x1
      vy1 = y2 - y1
      vz1 = z2 - z1
      vx2 = x3 - x2
      vy2 = y3 - y2
      vz2 = z3 - z2
      vx1 = vx1 - bx * ANINT (vx1 / bx)
      vy1 = vy1 - by * ANINT (vy1 / by)
      vz1 = vz1 - bz * ANINT (vz1 / bz)
      vx2 = vx2 - bx * ANINT (vx2 / bx)
      vy2 = vy2 - by * ANINT (vy2 / by)
      vz2 = vz2 - bz * ANINT (vz2 / bz)
      v1mag = SQRT ((vx1 ** 2) + (vy1 ** 2) + (vz1 ** 2))
      v2mag = SQRT ((vx2 ** 2) + (vy2 ** 2) + (vz2 ** 2))
      vx1 = vx1 / v1mag
      vy1 = vy1 / v1mag
      vz1 = vz1 / v1mag
      vx2 = vx2 / v2mag
      vy2 = vy2 / v2mag
      vz2 = vz2 / v2mag
      dot = vx1 * vx2 + vy1 * vy2 + vz1 * vz2
      angle = 180.0 - (DACOS (dot)) * (180.0 / 3.141592654)

      !Add hydrogen bond angle to histogram and average arrays
      ind = INT (angle / da) + 1 !calculate index for angle histogram
      Angles (k - 1, ind) = Angles (k - 1, ind) + 1 !increment angle histogram
      Atemp1 (k - 1) = angle !save hbond angle to temp array to be sorted

    END DO

      !Create copy of hydrogen bond list for jth NH4 ion, sort, add to avg tot.
      Btemp2 (:) = Btemp1 (:) !copy hbond list
      CALL dpquicksort (Btemp2) !sort hbond list
      Bavgsort (:) = Bavgsort (:) + Btemp2 (:) !add hbond values for running total
      Bavg (:) = Bavg (:) + Btemp1 (:) !running total of unsorted hbonds

      !sort hydrogen bond angles to keep initial association with hbonds
      DO l = 1, 4

        !Compares original Hbond list to sorted Hbond list
        !Perform same sorting of the indicies for hydrogen bond angles
        IF (Btemp1 (l) .EQ. Btemp2 (1)) THEN

            Atemp2 (1) = Atemp1 (l)

          ELSE IF (Btemp1 (l) .EQ. Btemp2 (2)) THEN

            Atemp2 (2) = Atemp1 (l)

          ELSE IF (Btemp1 (l) .EQ. Btemp2 (3)) THEN

            Atemp2 (3) = Atemp1 (l)

          ELSE

            Atemp2 (4) = Atemp1 (l)

        END IF

      END DO

      !add hbond angle values to running total to be averaged at the end
      Aavg (:) = Aavg (:) + Atemp1 (:) !unsorted
      Aavgsort (:) = Aavgsort (:) + Atemp2 (:) !sorted

  END DO !end loop over jth NH4 ion

END DO !end loop over ith snapshort

!open file for Hbond histogram
!WRITE (filename1, '(A)') "HBondhist.txt"
!OPEN (UNIT = 20, FILE = filename1, STATUS = 'UNKNOWN')

!time average hbond histogram
!Bonds (:, :) = Bonds (:, :) / REAL (n_step)

!write hbond histogram to file
!DO i = b_start, b_end

!  WRITE (20, *) (i - 1) * db, Bonds (:, i)

!ENDDO

!CLOSE (UNIT = 20, STATUS = 'keep') !close Hbond histogram file

!open file for Hbond angle histogram
!WRITE (filename2, '(A)') "HAnglehist.txt"
!OPEN (UNIT = 30, FILE = filename2, STATUS = 'UNKNOWN')

!time average hbond angle histogram
!Angles (:, :) = Angles (:, :) / REAL (n_step)

!write hbond angle histogram to file
!DO i = a_start, a_end

!  WRITE (30, *) (i - 1) * da, Angles (:, i)

!ENDDO

!CLOSE (UNIT = 30, STATUS = 'KEEP') !close Hbond angle histogram file

!time and ion average the different hydrogen bonds and angles
Bavg (:) = Bavg (:) / REAL (n_step * n_ion) !unsorted bond averaging
Bavgsort (:) = Bavgsort (:) / REAL (n_step * n_ion) !sorted bond averaging
Aavg (:) = Aavg (:) / REAL (n_step * n_ion) !unsorted angle averaging
Aavgsort (:) = Aavgsort (:) / REAL (n_step * n_ion) !sorted angle averaging
WRITE (40, *) 'Unsorted', Bavg (:) !write unsorted average bonds
WRITE (40, *) 'Sorted', Bavgsort (:) !write sorted average bonds
WRITE (50, *) 'Unsorted', Aavg (:) !write unsorted average angles
WRITE (50, *) 'Sorted', Aavgsort (:) !write sorted average angles

CLOSE (UNIT = 40, STATUS = 'KEEP') !close average hbond file
CLOSE (UNIT = 50, STATUS = 'KEEP') !close average hbond angle file

!subroutine included for sorting algorithm
CONTAINS

  ! ripped from https://www.mjr19.org.uk/IT/sorts/sorts.f90
  ! dual pivot quicksort
  recursive subroutine dpquicksort(array)
    real(prec), intent(inout)::array(:)
    real(prec) :: temp,p1,p2
    integer :: i,j,last,l,k,g

    last=size(array)

    if (last.lt.40) then ! use insertion sort on small arrays
       do i=2,last
          temp=array(i)
          do j=i-1,1,-1
             if (array(j).le.temp) exit
             array(j+1)=array(j)
          enddo
          array(j+1)=temp
       enddo
       return
    endif
    p1=array(last/3)
    p2=array(2*last/3)
    if (p2.lt.p1) then
       temp=p1
       p1=p2
       p2=temp
    endif
    array(last/3)=array(1)
    array(1)=p1
    array(2*last/3)=array(last)
    array(last)=p2

    g=last
    l=2
    do while (array(l).lt.p1)
       l=l+1
    enddo
    k=l

    do while(k.lt.g)
       temp=array(k)
       if (temp.lt.p1) then
          array(k)=array(l)
          array(l)=temp
          l=l+1
       else if (temp.gt.p2) then
          do while(array(g-1).gt.p2)
             g=g-1
          enddo
          if (k.ge.g) exit
          g=g-1
          if (array(g).lt.p1) then
             array(k)=array(l)
             array(l)=array(g)
             array(g)=temp
             l=l+1
          else
             array(k)=array(g)
             array(g)=temp
          endif
       endif
       k=k+1
    enddo
    if (l.gt.2) then
       array(1)=array(l-1)
       array(l-1)=p1
       call dpquicksort(array(1:l-2))
    endif
    call dpquicksort(array(l:g-1))
    if (g.lt.last) then
       array(last)=array(g)
       array(g)=p2
       call dpquicksort(array(g+1:last))
    endif

  end subroutine dpquicksort

END PROGRAM HBOND
