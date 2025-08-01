module cartesianInitProcedures

  use universalVariables,           only : ZERO, INF
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use genericProcedures,            only : findCommon
  use numPrecision   
  use cartesianGenericProcedures,   only : testIntervalIntersection

  implicit none

contains

  !!
  !!
  !!
  pure subroutine testEdgeIntersection(vertices, edges, edgeIdx, circumscribedBallRadius, &
                                  targetDistance, centroid, currEdgeVector, currVertexIdxs, a, &
                                  phi, phiCapital)
    type(topologicalObjectShelf), intent(in)            :: vertices, edges
    integer(shortInt), intent(in)                       :: edgeIdx
    integer(shortInt), intent(inout)                    :: phi, phiCapital
    real(defReal), intent(in)                           :: circumscribedBallRadius, targetDistance, a
    real(defReal), dimension(3), intent(in)             :: centroid, currEdgeVector
    integer(shortInt), dimension(2), intent(in)         :: currVertexIdxs
    real(defReal), dimension(3)                         :: dummyVector
    real(defReal)                                       :: b, c, discriminant, edgeLength, &
                                                           sqrtDiscriminant
    real(defReal), dimension(2)                         :: t, targetDistanceRatio
    integer(shortInt)                                   :: i

    ! calculate constants
    edgeLength = edges % getEdgeLength(edgeIdx)

    ! calculate dummyVector = vertex1 coordinate - centroid of the cube
    dummyVector = vertices % getVertexCoordinates(currVertexIdxs(1)) - centroid

    ! calculate constants for quadratic formula
    b = dot_product(dummyVector, currEdgeVector)
    c = dot_product(dummyVector, dummyVector) - circumscribedBallRadius**2
    discriminant = b**2 - a*c
  
    ! calculate t and leave the subroutine if t is invalid
    if (discriminant < ZERO) then
        return
    elseif (discriminant == ZERO) then
        if (b > ZERO) then
            return
        else
            t(:) = -b/a
        end if
    else
      sqrtDiscriminant = sqrt(discriminant)
      t(1) = (-b - sqrtDiscriminant)/a
      t(2) = (-b + sqrtDiscriminant)/a

      if (t(1) > 1) return
      if (t(2) < 0) return
    end if

    ! If passed to this point, t is valid. Hence, construct phi mapping
    if (t(2) < 0.5) then
      phi = currVertexIdxs(1)
    else
      phi = currVertexIdxs(2)
    end if

    ! test and constuct phiCapital mapping
    targetDistanceRatio(1) = targetDistance/edgeLength
    targetDistanceRatio(2) = (edgeLength - targetDistance)/edgeLength
    do i = 1, 2
      if (targetDistanceRatio(1) < t(i) .AND. targetDistanceRatio(2) > t(i)) then
            phiCapital = edgeIdx
      end if

    end do

  end subroutine testEdgeIntersection

  !!
  !!
  !!
  subroutine testPolyhedronInclusion(faces, currElementFaceIdxs, centroid, &
                                     faceNormalSigns, elementIdx, chi)
    type(topologicalObjectShelf), intent(in)            :: faces
    integer(shortInt), dimension(:), intent(in)         :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)             :: centroid
    real(defReal), dimension(:,:), intent(in)           :: faceNormalSigns
    integer(shortInt), intent(in)                       :: elementIdx
    integer(shortInt), intent(inout)                    :: chi
    integer(shortInt)                                   :: i, j
    real(defReal), dimension(3)                         :: furthestVertexCoord

    do i = 1, size(currElementFaceIdxs)
      
      do j = 1, 3
        furthestVertexCoord(j) = centroid(j) + faceNormalSigns(j,i) 
      end do

      if (dot_product(faces % getFaceNormal(currElementFaceIdxs(i)), furthestVertexCoord) &
          + faces % getFaceConst(currElementFaceIdxs(i)) > 0) then
          ! if TRUE, then at least a part of this cell lies outside of the polyhedron

          !!!!!
          ! print*, "----------------------------------------------"
          ! print*, "fail"
          ! print*, currElementFaceIdxs(i)
          ! print*, dot_product(faces % getFaceNormal(currElementFaceIdxs(i)), furthestVertexCoord) &
          !         + faces % getFaceConst(currElementFaceIdxs(i))
          !!!!!

          return
      end if 

          !!!!!
          ! print*, "----------------------------------------------"
          ! print*, "passed"
          ! print*, currElementFaceIdxs(i)
          ! print*, furthestVertexCoord
          ! print*, dot_product(faces % getFaceNormal(currElementFaceIdxs(i)), furthestVertexCoord) &
          !         + faces % getFaceConst(currElementFaceIdxs(i))
          ! print*, "----------------------------------------------"
          !!!!!

    end do

    ! if survived to this point, then the cell is entired enclosed by the polyhedron. 
    ! Hence, set chi mapping to the current element index.
    chi = elementIdx

  end subroutine testPolyhedronInclusion

  !!
  !!
  !!
  subroutine testFaceIntersection(vertices, edges, faces, currVertexIdxs, extraDistance, currFaceNormal, centroid, &
                                  cellSpacing, faceIdx, currFaceEdgeIdxs, targetDistance, chi, phi, phiCapital, faceIdxs)
    type(topologicalObjectShelf), intent(in)            :: vertices, edges, faces
    integer(shortInt), dimension(:), intent(in)         :: currVertexIdxs, currFaceEdgeIdxs
    real(defReal), dimension(3), intent(in)             :: currFaceNormal, centroid
    real(defReal), intent(in)                           :: extraDistance, cellSpacing, targetDistance
    integer(shortInt), intent(in)                       :: faceIdx, chi
    integer(shortInt), intent(inout)                    :: phi, phiCapital
    integer(shortInt), dimension(2), intent(inout)      :: faceIdxs
    integer(shortInt)                                   :: i, j
    real(defReal), dimension(3)                         :: currVertexCoords, min1, max1, currEdgeUnitVector
    real(defReal)                                       :: faceConst, currValue, vectorDotCentroid, extraDistance2 

    !------------------------------------------------------------------------------------------------
    !if cell is contained within a polyhedron, the cell cannot intersect with the polyhedron's faces
    !Or, if phi and phiCapital mapping informations are assigned already from edge interesetion, no need to find new one.
    !------------------------------------------------------------------------------------------------
    if (chi /= 0) then
      return
    elseif (phi /= 0 .AND. phiCapital /= 0) then
      return
    end if

    !------------------------------------------------------------------------------------------------
    ! Testing along face normal
    !------------------------------------------------------------------------------------------------
    vectorDotCentroid = dot_product(currFaceNormal, centroid)
    faceConst = faces % getFaceConst(faceIdx)

    if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance, vectorDotCentroid + extraDistance, &
        -faceConst, -faceConst)) then
          !!!!!
          ! print*, "faceNormal"
          ! print*, vectorDotCentroid - extraDistance
          ! print*, vectorDotCentroid + extraDistance
          ! print*, -faceConst
          !!!!!
      return
    end if

    !------------------------------------------------------------------------------------------------
    ! Testing along crossProduct(each of edgeUnitVector, three coordinate basis)
    !------------------------------------------------------------------------------------------------
    ! loop over all edgesUnitVectors of the current face
    do i = 1, size(currFaceEdgeIdxs)
      currEdgeUnitVector = edges % getEdgeUnitvector(currFaceEdgeIdxs(i))
      
      ! loop over all vertices of the current face to find the min and max of polygon interval
      do j = 1, size(currVertexIdxs)
        currVertexCoords = vertices % getVertexCoordinates(currVertexIdxs(j))

        ! for coordinate basis = (1,0,0)
        currValue = currEdgeUnitVector(3)*currVertexCoords(2) - currEdgeUnitVector(2)*currVertexCoords(3)
        if (j == 1) then
          min1(1) = currValue
          max1(1) = currValue
        else
          if (currValue < min1(1)) then
            min1(1) = currValue
          elseif (currValue > max1(1)) then
            max1(1) = currValue
          end if 
        end if

        ! for coordinate basis = (0,1,0)
        currValue = - currEdgeUnitVector(3)*currVertexCoords(1) + currEdgeUnitVector(1)*currVertexCoords(3)
        if (j == 1) then
          min1(2) = currValue
          max1(2) = currValue
        else
          if (currValue < min1(2)) then
            min1(2) = currValue
          elseif (currValue > max1(2)) then
            max1(2) = currValue
          end if 
        end if

        ! for coordinate basis = (0,0,1)
        currValue =  currEdgeUnitVector(2)*currVertexCoords(1) - currEdgeUnitVector(1)*currVertexCoords(2)
        if (j == 1) then
          min1(3) = currValue
          max1(3) = currValue
        else
          if (currValue < min1(3)) then
            min1(3) = currValue
          elseif (currValue > max1(3)) then
            max1(3) = currValue
          end if 
        end if

      end do

      ! test intersections of interval
      vectorDotCentroid = currEdgeUnitVector(3)*centroid(2)-currEdgeUnitVector(2)*centroid(3)
      extraDistance2 = (abs(currEdgeUnitVector(3)) + abs(currEdgeUnitVector(2)))*cellSpacing*0.5
      if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
                                                min1(1), max1(1))) then
          !!!!!
          ! print*, "Crossx"
          !!!!!
          return
      end if

      vectorDotCentroid = currEdgeUnitVector(1)*centroid(3)-currEdgeUnitVector(3)*centroid(1)
      extraDistance2 = (abs(currEdgeUnitVector(3)) + abs(currEdgeUnitVector(1)))*cellSpacing*0.5
      if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
                                                min1(2), max1(2))) then
          !!!!!
          ! print*, "Crossy"
          !!!!!
          return
      end if

      vectorDotCentroid = currEdgeUnitVector(2)*centroid(1)-currEdgeUnitVector(1)*centroid(2)
      extraDistance2 = (abs(currEdgeUnitVector(2)) + abs(currEdgeUnitVector(1)))*cellSpacing*0.5
      if (.NOT. testIntervalIntersection(vectorDotCentroid - extraDistance2, vectorDotCentroid + extraDistance2, &
                                                min1(3), max1(3))) then
          !!!!!
          ! print*, "Crossz"
          !!!!!
          return
      end if

    end do

    !------------------------------------------------------------------------------------------------
    ! if survived to this point, then there is no separating axis. Hence, construct mapping accordingly
    !------------------------------------------------------------------------------------------------

    ! (needs to be changed) (due to memory, taking shortCut)
    ! (originally, we just have to add faceIdx to an array of intersected faces)
    if (faceIdxs(1) == 0) then
      faceIdxs(1) = faceIdx
    else
      faceIdxs(2) = faceIdx
      ! tests if the two faces have a common edge. If yes, assign phi and phiCapital mappings.
      ! if there are more than two faces intersected with the cell, and phi and phicaptial have been assigned already,
      ! then, this cell is complete in terms of mapping construction, and this subroutine is terminated at the beginning
      ! of this subroutine.
      call testTwoIntersectedFaces(vertices, edges, faces, centroid, targetDistance, phi, phiCapital, faceIdxs)
    end if
    
    !!!!!
    ! print*, "PASSESD"
    !!!!!

  end subroutine testFaceIntersection

  !!
  !!
  !! 
  subroutine testTwoIntersectedFaces(vertices, edges, faces, centroid, targetDistance, phi, phiCapital, faceIdxs)
    type(topologicalObjectShelf), intent(in)            :: vertices, edges, faces
    real(defReal), dimension(3), intent(in)             :: centroid
    real(defReal), intent(in)                           :: targetDistance
    integer(shortInt), intent(inout)                    :: phi, phiCapital
    integer(shortInt), dimension(2), intent(in)         :: faceIdxs
    integer(shortInt)                                   :: commonEdgeIdx
    integer(shortInt), dimension(:), allocatable        :: commonEdgeVertexIdxs
    real(defReal)                                       :: dummyConstant, dotProduct1, dotProduct2
    real(defReal), dimension(3)                         :: targetFaceNormal, edgeVertexCoord, c_1, &
                                                           edgeUnitVector, edgeVertex2Coord, &
                                                           tempVector1, tempVector2, chiCoord

    ! find the common edge index of the given two faces intersected by the cell 
    commonEdgeIdx = faces % findCommonedgeIdx(faceIdxs(1), faceIdxs(2))

    !!!!!
    !print*, "ENTERED TWO FACES"
    !!!!!

    ! if there is no common edge, exit the subroutine early (other combinations of faces to be tried later)
    if (commonEdgeIdx == 0) then
      !!!!!
      !print*, "Exiting two faces"
      !!!!!
      return
    end if 

    ! calculate centre of intersection between the circumscribed ball and the plane parallel to the polygon
    commonEdgeVertexIdxs = edges % getEdgeVertexIdxs(commonEdgeIdx)
    targetFaceNormal = faces % getFaceNormal(faceIdxs(1))
    edgeVertexCoord = vertices % getVertexCoordinates(commonEdgeVertexIdxs(1))
    dummyConstant = dot_product(centroid - edgeVertexCoord, targetFaceNormal)
    c_1 = centroid - dummyConstant*targetFaceNormal

    ! calculate chi
    edgeUnitVector = edges % getEdgeUnitVector(commonEdgeIdx)
    dummyConstant = dot_product(c_1 - edgeVertexCoord, edgeUnitVector)
    chiCoord = edgeVertexCoord + dummyConstant*edgeUnitVector

    ! calculate distances to be compared
    edgeVertex2Coord = vertices % getVertexCoordinates(commonEdgeVertexIdxs(2))
    tempVector1 = chiCoord - edgeVertexCoord
    tempVector2 = chiCoord - edgeVertex2Coord
    dotProduct1 = dot_product(tempVector1, tempVector1)
    dotProduct2 = dot_product(tempVector2, tempVector2)

    ! make comparison between the two distances (dotProduct1/2)
    if (dotProduct1 < dotProduct2) then
      phi = commonEdgeVertexIdxs(1)
    else
      phi = commonEdgeVertexIdxs(2)
    end if

    if (dotProduct1 > targetDistance .AND. dotProduct2 > targetDistance) then
      phiCapital = commonEdgeIdx
    end if
    
  end subroutine testTwoIntersectedFaces

  !!
  !!
  !!
  pure subroutine constructMapSingleFace(faces, faceIdxs, phiCapital)
    type(topologicalObjectShelf), intent(in)            :: faces
    integer(shortInt), dimension(2), intent(in)         :: faceIdxs
    integer(shortInt), intent(inout)                    :: phiCapital
    integer(shortInt), dimension(:), allocatable        :: currFaceEdgeIdxs

    ! testing if the cell intersects with a single face
    ! (necessary to test this because this subroutine is called for all cells)
    ! no need to test if chi or both (phi and phiCapital) have been assigned from edge intersection or 
    ! polyhedron inclusion because these cells are not tested against face intersection (self % faceIdxs(1) = 0)
    if (faceIdxs(1) /= 0 .AND. faceIdxs(2) == 0) then 

      ! (needs to be changed) (rather than saving candidate indices, write a subroutine that directly returns the first index of the array)
      ! (Does this subroutine has to save the array locally anyways?)
      currFaceEdgeIdxs = faces % getFaceEdgeIdxs(faceIdxs(1))

      ! (needs to be changed) (can find one in currFaceEdgeIdxs that are already been assigned to other cells)
      ! (Hence, we reduce the number of edges about which angles are calculated as well as memory if virtual grid is used)
      phiCapital = currFaceEdgeIdxs(1)

    end if

  end subroutine constructMapSingleFace

  !!
  !!
  !!
  subroutine coverFinitePrecision(faces, elements, centroid, elementIdx, chi)
    type(topologicalObjectShelf), intent(in)            :: faces, elements
    real(defReal), dimension(3), intent(in)             :: centroid
    integer(shortInt), intent(in)                       :: elementIdx
    integer(shortInt), intent(inout)                    :: chi
    integer(shortInt), dimension(:), allocatable        :: currElementFaceIdxs
    integer(shortInt)                                   :: i

    currElementFaceIdxs = elements % getElementFaceIdxs(elementIdx)

    do i = 1, size(currElementFaceIdxs)

      if (dot_product(faces % getFaceNormal(currElementFaceIdxs(i)), centroid) &
          + faces % getFaceConst(currElementFaceIdxs(i)) > 0) then
          ! if TRUE, then centroid of this cell lies outside of the polyhedron
          return
      end if 

    end do

    !print*, "TESTING", centroid
    ! If survived to this point, this cell is contained within a single mesh element, but
    ! chi mapping info was not assigned due to floating-point-finite-precision.
    chi = elementIdx

    !print*, "hahahoho"!!!!!

  end subroutine coverFinitePrecision




end module cartesianInitProcedures