module polyhedron_class
  
  use edgeShelf_class,    only : edgeShelf
  use element_inter,      only : element, elementBox, kill_super => kill, notch
  use face_inter,         only : faceBox
  use faceShelf_class,    only : faceShelf
  use genericProcedures,  only : append, areEqual, computePyramidCentre, computePyramidVolume, &
                                 computeTetrahedronCentre, computeTetrahedronVolume, findCommon, fatalError, numToChar
  use numPrecision
  use tetrahedron_class,  only : tetrahedron
  use triangle_class,     only : triangle
  use universalVariables, only : SURF_TOL, INF, ZERO
  use vertexShelf_class,  only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Element (cell) of an OpenFOAM mesh. Consists of a list of vertices and faces indices, as well
  !! as a list of tetrahedra indices into which the element is decomposed.
  !!
  !! Private members:
  !!   idx      -> Index of the element.
  !!   vertices -> Array of vertices indices making the element up.
  !!   faces    -> Array of faces indices making the element up.
  !!   Volume   -> Volume of the element.
  !!   Centroid -> Vector pointing to the centroid of the element.
  !!
  type, public, extends(element) :: polyhedron
    private
  contains
    ! Build procedures.
    procedure                    :: computeComponents
    procedure                    :: split
    procedure                    :: splitConcave
    ! Runtime procedures.
    procedure                    :: kill
  end type polyhedron

contains
  
  !! Subroutine 'computeVolumeAndCentroid'
  !!
  !! Basic description:
  !!   Computes the volume and volume-weighted centroid of the element.
  !!
  !! Detailed description:
  !!   First estimates the centroid of the element by taking an area-weighted average of the faces' 
  !!   centroids. Using this estimate, the element is decomposed into a number of pyramids whose 
  !!   apices are the estimated centroid of the element. Looping through all the faces, the volume 
  !!   and centroid of each pyramid is computed (see pyramid_class for more details); the centroid 
  !!   of the element is then obtained by taking a volume-weighted average of the pyramids' 
  !!   centroids.
  !!
  !! Arguments:
  !!   faces [in]    -> An array of face structures making the element up.
  !!   vertices [in] -> An array of vertex structures making the element up.
  !!
  !! Error:
  !!   fatalError if the volume is element is negative or infinite.
  !!
  pure subroutine computeComponents(self, faceIdxs, vertexIdxs, faces, vertices, centroid, volume)
    class(polyhedron), intent(inout)             :: self
    integer(shortInt), dimension(:), intent(in)  :: faceIdxs, vertexIdxs
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    real(defReal), dimension(3), intent(out)     :: centroid
    real(defReal), intent(out)                   :: volume
    integer(shortInt)                            :: i, nFaces, nVertices, absFaceIdx
    real(defReal)                                :: area, sumAreas, sumVolumes
    real(defReal), dimension(3)                  :: sumVolumesCentroid, C
    real(defReal), dimension(:, :), allocatable  :: array
    character(100), parameter                    :: Here = 'computeVolumeAndCentroid (element_class.f90)'

    ! Retrieve the number of faces and vertices in the element and initialise variables.
    nFaces = size(faceIdxs)
    nVertices = size(vertexIdxs)
    C = ZERO
    sumAreas = ZERO
    
    ! Loop through all faces and compute the element's approximate centroid
    ! by performing an area-weighted average of the different faces' centroids.
    do i = 1, nFaces
      ! Retrieve the area of the current face.
      absFaceIdx = abs(faceIdxs(i))
      area = faces % getFaceArea(absFaceIdx)
      ! Update the area-weighted centroid and the sum of faces' areas.
      C = C + faces % getFaceCentroid(absFaceIdx) * area
      sumAreas = sumAreas + area

    end do
    
    ! Using the approximate centroid, compute the actual centroid by performing a volume-weighted
    ! average of the different pyramids' centroids.
    sumVolumes = ZERO
    sumVolumesCentroid = ZERO
    allocate(array(2, 3))
    C = C / sumAreas
    
    ! Loop through all faces (pyramids).
    do i = 1, nFaces
      ! Retrieve the volume of the current pyramid and update the volume-weighted centroid and the sum of volumes.
      absFaceIdx = abs(faceIdxs(i))
      array(1, :) = faces % getFaceNormal(absFaceIdx) * faces % getFaceArea(absFaceIdx)
      array(2, :) = C - faces % getFaceCentroid(absFaceIdx)
      volume = computePyramidVolume(array)

      array(1, :) = 3.0_defReal * faces % getFaceCentroid(absFaceIdx)
      array(2, :) = C
      sumVolumesCentroid = sumVolumesCentroid + computePyramidCentre(array) * volume
      sumVolumes = sumVolumes + volume

    end do
    ! The volume of the element is simply the sum of volumes, while the centroid is the average of
    ! the volume-weighted sum.
    volume = sumVolumes
    centroid = sumVolumesCentroid / sumVolumes

  end subroutine computeComponents
  
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(polyhedron), intent(inout) :: self
    
    ! Element.
    call kill_super(self)

  end subroutine kill
  
  !! Subroutine 'split'
  !!
  !! Basic description:
  !!   Splits the element into a number of pyramids corresponding to the number of faces in the 
  !!   element.
  !!
  !! Detailed description:
  !!   The element is split into pyramids by subdividing it from its centroid: this then forms the
  !!   apex of each generated pyramid. For each face in the element, triangles are created by
  !!   joining each edge in the face to the common apex. During the creation of triangles a check is
  !!   made to ensure that their normal vectors point in the correct direction (from the pyramid of
  !!   lowest index to that of greatest index, just as OpenFOAM does for elements).
  !!
  !! Arguments:
  !!   faces [in]              -> A faceShelf.
  !!   edges [inout]           -> An edgeShelf.
  !!   vertices [inout]        -> A vertexShelf.
  !!   triangles [inout]       -> A triangleShelf.
  !!   pyramids [inout]        -> A pyramidShelf.
  !!   lastEdgeIdx [inout]     -> Index of the first free item in the edgeShelf.
  !!   lastTriangleIdx [inout] -> Index of the first free item in the triangleShelf.
  !!   lastPyramidIdx [inout]  -> Index of the first free item in the pyramidShelf.
  !!   lastVertexIdx [in]      -> Index of the first free item in the vertexShelf.
  !!
  subroutine split(self, faces, lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx, &
                        newEdges, newFaces, newVertices, tetrahedra, triangles)
    class(polyhedron), intent(inout)              :: self
    type(faceShelf), intent(inout)                :: faces, newFaces
    integer(shortInt), intent(inout)              :: lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx
    type(edgeShelf), intent(inout)                :: newEdges
    type(vertexShelf), intent(inout)              :: newVertices
    type(elementBox), dimension(:), intent(inout) :: tetrahedra
    type(faceBox), dimension(:), intent(inout)    :: triangles
    integer(shortInt)                             :: i, j, k, l, faceIdx, absFaceIdx, faceTriangleIdx, triangleVertexIdx, &
                                                     nFaces, nVertices, edgeIdx, commonTriangleIdx
    integer(shortInt), dimension(:), allocatable  :: faceIdxs, faceTriangleIdxs, vertexIdxs
    integer(shortInt), dimension(2)               :: edgeVertexIdxs
    integer(shortInt), dimension(3)               :: triangleVertexIdxs, testVertexIdxs
    integer(shortInt), dimension(4)               :: triangleIdxs
    integer(shortInt), dimension(6)               :: edgeIdxs
    real(defReal)                                 :: volume
    real(defReal), dimension(3)                   :: centroid
    real(defReal), dimension(4, 3)                :: array
    
    ! Initialise a new vertex corresponding to the centroid of the polyhedron.
    lastNewVertexIdx = lastNewVertexIdx + 1
    call newVertices % initVertex(lastNewVertexIdx, self % getCentroid())
    array(4, :) = self % getCentroid()

    ! Retrieve the indices of the vertices in the polyhedron and compute the number of vertices.
    vertexIdxs = self % getVertexIdxs()
    nVertices = size(vertexIdxs)

    ! Loop through all the vertices in the polyhedron and create new edges joining each vertex to the centroid.
    do i = 1, nVertices
      lastNewEdgeIdx = lastNewEdgeIdx + 1
      edgeVertexIdxs = [vertexIdxs(i), lastNewVertexIdx]
      call newEdges % initEdge(lastNewEdgeIdx, [vertexIdxs(i), lastNewVertexIdx])

      do j = 1, 2
        call newVertices % addEdgeIdxToVertex(edgeVertexIdxs(j), lastNewEdgeIdx)

      end do

    end do

    ! Retrieve the indices of the faces in the polyhedron and compute the number of faces.
    faceIdxs = self % getFaceIdxs()
    nFaces = size(faceIdxs)

    ! Loop through all the faces in the polyhedron.
    do i = 1, nFaces
      ! Retrieve the current face index and create its absolute value.
      faceIdx = faceIdxs(i)
      absFaceIdx = abs(faceIdx)

      ! Retrieve the indices of the triangles created from the current face.
      faceTriangleIdxs = faces % getFaceTriangleIdxs(absFaceIdx)
      do j = 1, size(faceTriangleIdxs)
        ! Increment lastNewElementIdx and retrieve the indices of the vertices in the triangle.
        lastNewElementIdx = lastNewElementIdx + 1
        faceTriangleIdx = faceTriangleIdxs(j)
        call newFaces % addElementIdxToFace(faceTriangleIdx, lastNewElementIdx)
        
        triangleVertexIdxs = triangles(faceTriangleIdx) % item % getVertexIdxs()
        edgeIdxs(1:3) = triangles(faceTriangleIdx) % item % getEdgeIdxs()

        ! Create the indices of the vertices in the new tetrahedron.
        vertexIdxs = [triangleVertexIdxs, lastNewVertexIdx]
        call newVertices % addElementIdxToVertex(lastNewVertexIdx, lastNewElementIdx)
        triangleIdxs(1) = sign(faceTriangleIdx, faceIdx)

        ! Create remaining array entries and update mesh connectivity.
        do k = 1, 3
          triangleVertexIdx = triangleVertexIdxs(k)
          array(k, :) = newVertices % getVertexCoordinates(triangleVertexIdx)
          call newVertices % addElementIdxToVertex(triangleVertexIdx, lastNewElementIdx)

        end do

        ! Compute tetrahedron centroid and volume.
        centroid = computeTetrahedronCentre(array)
        volume = computeTetrahedronVolume(array)

        ! Loop through all the remaining faces in the new tetrahedron.
        do k = 1, 3
          edgeIdxs(k + 3) = newVertices % findCommonEdgeIdx(triangleVertexIdxs(k), lastNewVertexIdx)
          call newEdges % addElementIdxToEdge(edgeIdxs(k), lastNewElementIdx)
          call newEdges % addElementIdxToEdge(edgeIdxs(k + 3), lastNewElementIdx)
          
          ! Check if a triangle containing the three vertices already exists.
          testVertexIdxs = [triangleVertexIdxs(k), triangleVertexIdxs(mod(k, 3) + 1), lastNewVertexIdx]
          commonTriangleIdx = newVertices % findCommonFaceIdx(testVertexIdxs)
          if (commonTriangleIdx > 0) then
            call newFaces % addElementIdxToFace(commonTriangleIdx, lastNewElementIdx)
            triangleIdxs(k + 1) = -commonTriangleIdx
            cycle

          end if

          ! Create a new internal triangle.
          lastNewFaceIdx = lastNewFaceIdx + 1
          allocate(triangle :: triangles(lastNewFaceIdx) % item)
          call triangles(lastNewFaceIdx) % item % build(lastNewFaceIdx, 0, .false., testVertexIdxs, newVertices, 'Triangle', &
                                                        centroid)
          call triangles(lastNewFaceIdx) % item % addElementIdx(lastNewElementIdx)
          triangleIdxs(k + 1) = lastNewFaceIdx

          ! Update mesh connectivity.
          do l = 1, 3
            edgeIdx = newVertices % findCommonEdgeIdx(testVertexIdxs(l), testVertexIdxs(mod(l, 3) + 1))
            call newEdges % addFaceIdxToEdge(edgeIdx, lastNewFaceIdx)
            call triangles(lastNewFaceIdx) % item % addEdgeIdx(edgeIdx)
            call newVertices % addFaceIdxToVertex(testVertexIdxs(l), lastNewFaceIdx)

          end do

          ! Set new triangle in the new faceShelf.
          call newFaces % addFace(lastNewFaceIdx, triangles(lastNewFaceIdx))

        end do

        ! Initialise new tetrahedron in the shelf.
        allocate(tetrahedron :: tetrahedra(lastNewElementIdx) % item)
        call tetrahedra(lastNewElementIdx) % item % init(lastNewElementIdx, self % getIdx(), triangleIdxs, vertexIdxs, &
                                                         centroid, volume, .true., 'Tetrahedron', edgeIdxs)

      end do

    end do
  
  end subroutine split

  !!
  !!
  !!
  subroutine splitConcave(self, edges, faces, vertices, newEdges, convexElements, newFaces, newVertices)
    class(polyhedron), intent(inout)              :: self
    type(edgeShelf), intent(inout)                :: edges, newEdges
    type(faceShelf), intent(inout)                :: faces, newFaces
    type(vertexShelf), intent(inout)              :: vertices, newVertices
    type(elementBox), dimension(:), allocatable, intent(out) :: convexElements
    integer(shortInt)                             :: i, j, k, l, edgeIdx, nEdges, nElements, nVertices, commonFaceIdx, &
                                                     firstVertexIdx, cutVertexIdx, idx, previousIdx, &
                                                     minPositiveIdx, minNegativeIdx, nFaces, nInitialVertices
    integer(shortInt), dimension(2)               :: edgeVertexIdxs, faceIdxs, verticesEdge, newEdgeVertexIdxs, &
                                                     edgeFaceIdxs, childEdgeIdxs
    real(defReal), dimension(3)                   :: u, v, normalDifference, CoorVertexNotch, coord1, coord2, &
                                                     newVertexCoords
    real(defReal)                                 :: t, denominator, const, newEdgeLength, testLength, dotProduct, &
                                                     currentVertexDotProduct, prev
    integer(shortInt), dimension(:), allocatable  :: edgeIdxs, edgeFaceA, edgeFaceB, newVertexIdxs, faceToSplitIdxs, &
                                                     faceVertexIdxs, PosNewFaceVertexIdxs, NegNewFaceVertexIdxs, &
                                                     minPositiveIdxs, minNegativeIdxs, newFaceVertexIdxs, vertexEdgeIdxs, &
                                                     faceEdgeIdxs, newFaceEdgeIdxs, positiveIdxs, negativeIdxs, newElementVertexIdxs
    real(defReal), dimension(:), allocatable      :: dotProducts
    character(:), allocatable                     :: type
    type(notch), dimension(:), allocatable        :: notches
    
    ! Initialise nEdges, nElements, nFaces, and nVertices.
    nEdges = 0
    nElements = 0
    nFaces = 0
    nVertices = newVertices % getSize()
    nInitialVertices = nVertices

    ! Allocate newVertexIdxs to zero-size.
    allocate(newVertexIdxs(0))
    allocate(faceToSplitIdxs(0))

    ! Loop through all notches.
    notches = self % getNotches()
    do i = 1, size(notches)

      ! Step 1: retrieve first direction vector from the edge of the current notch.
      edgeVertexIdxs = edges % getEdgeVertexIdxs(notches(i) % edgeIdx)
      u = vertices % getVertexCoordinates(edgeVertexIdxs(2)) - vertices % getVertexCoordinates(edgeVertexIdxs(1))
      u = u / norm2(u)

      ! Step 2: retrieve normal vectors for each face in the current notch, and find the normal of the bisector
      !         as a difference between the two
      faceIdxs = notches(i) % faceIdxs
      normalDifference = faces % getFaceNormal(faceIdxs(2)) - faces % getFaceNormal(faceIdxs(1))
      normalDifference = normalDifference / norm2(normalDifference)

      ! Get all the edge indices for this element
      edgeIdxs = self % getEdgeIdxs()
      edgeFaceA = faces % getFaceEdgeIdxs(notches(i) % faceIdxs(1))
      edgeFaceB = faces % getFaceEdgeIdxs(notches(i) % faceIdxs(2))

      ! Retrieve the coordinates of one of the vertices forming the current edge
      CoorVertexNotch = vertices % getVertexCoordinates(edgeVertexIdxs(1))

      ! Loop through all the edges in the element to find the intersection points with the cut plane
      do j = 1, size(edgeIdxs)
        edgeIdx = edgeIdxs(j)
        verticesEdge = edges % getEdgeVertexIdxs(edgeIdx)

        ! To find intersection points, skip edges which are already in the problematic faces
        if (any(edgeFaceA == edgeIdx) .or. any(edgeFaceB == edgeIdx)) then
          do k = 1, nEdges
            if (size(findCommon(newEdges % getEdgeVertexIdxs(k), verticesEdge)) == 2) cycle

          end do

          ! Copy current edge in the new edgeShelf.
          nEdges = nEdges + 1
          if (nEdges > newEdges % getSize()) call newEdges % expandShelf(1)
          call newEdges % initEdge(nEdges, verticesEdge)
          call edges % addChildIdxToEdge(edgeIdx, nEdges)

          ! Update connectivity information for the vertices in the new edge.
          do k = 1, 2
            call newVertices % addEdgeIdxToVertex(verticesEdge(k), nEdges)

          end do

          cycle

        end if

        ! Retrieve coordinates of the two vertices in the edge
        coord1 = vertices % getVertexCoordinates(verticesEdge(1))
        coord2 = vertices % getVertexCoordinates(verticesEdge(2))

        ! Compute denominator and cycle to the next edge if areEqual(denominator, ZERO).
        denominator = dot_product(normalDifference, coord2 - coord1)
        if (areEqual(denominator, ZERO)) then
          do k = 1, nEdges
            if (size(findCommon(newEdges % getEdgeVertexIdxs(k), verticesEdge)) == 2) cycle

          end do

          ! Copy current edge in the new edgeShelf.
          nEdges = nEdges + 1
          if (nEdges > newEdges % getSize()) call newEdges % expandShelf(1)
          call newEdges % initEdge(nEdges, verticesEdge)
          call edges % addChildIdxToEdge(edgeIdx, nEdges)

          ! Update connectivity information for the vertices in the new edge.
          do k = 1, 2
            call newVertices % addEdgeIdxToVertex(verticesEdge(k), nEdges)

          end do

          cycle

        end if
        
        ! Calculate the instersection point. Start by calculating t
        const = dot_product(normalDifference, CoorVertexNotch)
        t = abs((dot_product(normalDifference, coord1) + const) / denominator)
        
        ! If t calculated is not valid, skip the current edge.
        if (t <= ZERO .or. t >= ONE) then
          ! Create a new edge only if not already present in the new edgeShelf. Retrieve all
          ! edges linked to the first vertex of the current edge.
          vertexEdgeIdxs = newVertices % getVertexEdgeIdxs(verticesEdge(1))
          do k = 1, size(vertexEdgeIdxs)
            if (size(findCommon(newEdges % getEdgeVertexIdxs(vertexEdgeIdxs(k)), verticesEdge)) == 2) cycle
            nEdges = nEdges + 1
            if (nEdges > newEdges % getSize()) call newEdges % expandShelf(1)
            call newEdges % initEdge(nEdges, verticesEdge)
            call edges % addChildIdxToEdge(edgeIdx, nEdges)

            ! Update connectivity information for the vertices in the new edge.
            do l = 1, 2
              call newVertices % addEdgeIdxToVertex(verticesEdge(l), nEdges)

            end do

          end do

        !If t calculated is valid, append the corresponding intersection point to the list.
        else

          ! Update nVertices and create a new vertex.
          nVertices = nVertices + 1
          call append(newVertexIdxs, nVertices)
          call newVertices % expandShelf(1)
          newVertexCoords = (ONE - t) * coord1 + t * coord2
          call newVertices % initVertex(nVertices, newVertexCoords)

          ! Update nEdges and create a new edge connecting one vertex in the notch edge to
          ! the new vertex.
          nEdges = nEdges + 1
          newEdgeLength = INF
          newEdgeVertexIdxs(2) = nVertices
          do k = 1, 2
            ! Compute distance of edge from the current vertex.
            testLength = norm2(newVertices % getVertexCoordinates(edgeVertexIdxs(k)) &
                               - newVertices % getVertexCoordinates(nVertices))
            
            ! Pick first vertex index which minimises new edge length.
            if (testLength < newEdgeLength) then
              newEdgeLength = testLength
              newEdgeVertexIdxs(1) = edgeVertexIdxs(k)

            end if

          end do
          if (nEdges > newEdges % getSize()) call newEdges % expandShelf(1)
          call newEdges % initEdge(nEdges, newEdgeVertexIdxs)

          do k = 1, 2
            call newVertices % addEdgeIdxToVertex(newEdgeVertexIdxs(k), nEdges)

          end do

          ! Find the index of the face being split and append it to the list of faces to
          ! split if it is not already present.
          edgeFaceIdxs = findCommon(abs(self % getFaceIdxs()), edges % getEdgeFaceIdxs(edgeIdx))
          call append(faceToSplitIdxs, edgeFaceIdxs, .true.)

          ! Split the original edge into two.
          do k = 1, 2
            nEdges = nEdges + 1
            newEdgeVertexIdxs(1) = verticesEdge(k)
            newEdgeVertexIdxs(2) = nVertices
            if (nEdges > newEdges % getSize()) call newEdges % expandShelf(1)
            call newEdges % initEdge(nEdges, newEdgeVertexIdxs)
            call edges % addChildIdxToEdge(edgeIdx, nEdges)

            ! Add connectivity information for new vertices.
            do l = 1, 2
              call newVertices % addEdgeIdxToVertex(newEdgeVertexIdxs(l), nEdges)

            end do

          end do

          ! Set the index of the cut vertex for the edge being split.
          call edges % setEdgeCutVertexIdx(edgeIdx, nVertices)

        end if 

      end do

      ! Loop through all the newly created vertices and create edges joining them.
      do j = 1, size(newVertexIdxs) - 1
        ! Create a new edge.
        nEdges = nEdges + 1
        newEdgeVertexIdxs = newVertexIdxs(j:j + 1)
        if (nEdges > newEdges % getSize()) call newEdges % expandShelf(1)
        call newEdges % initEdge(nEdges, newEdgeVertexIdxs)

        do k = 1, 2
          call newVertices % addEdgeIdxToVertex(newEdgeVertexIdxs(k), nEdges)

        end do

      end do

      ! Loop through all initial faces and copy those not cut by the plane to the newFaces shelf.
      do j = 1, faces % getSize()
        if (any(faceToSplitIdxs == j)) cycle

        ! Retrieve the indices of the vertices and edges in the original face.
        faceVertexIdxs = faces % getFaceVertexIdxs(j)
        faceEdgeIdxs = faces % getFaceEdgeIdxs(j)
        if (allocated(newFaceEdgeIdxs)) deallocate(newFaceEdgeIdxs)
        do k = 1, size(faceEdgeIdxs)
          call append(newFaceEdgeIdxs, edges % getEdgeChildrenIdxs(faceEdgeIdxs(k)))

        end do

        nFaces = nFaces + 1
        call newFaces % expandShelf(1)
        call newFaces % initFace(nFaces, 0, faces % getFaceIsBoundary(j), faces % getFaceVertexIdxs(j), &
                                 faces % getFaceAB(j), faces % getFaceAC(j), faces % getFaceCentroid(j), &
                                 faces % getFaceNormal(j), faces % getFaceArea(j), faces % getFaceType(j), &
                                 newFaceEdgeIdxs)

        ! Update mesh connectivity information.
        do k = 1, size(faceVertexIdxs)
          call newVertices % addFaceIdxToVertex(faceVertexIdxs(k), nFaces)

        end do

        do k = 1, size(newFaceEdgeIdxs)
          call newEdges % addFaceIdxToEdge(newFaceEdgeIdxs(k), nFaces)

        end do

      end do

      ! Loop through all faces to be split and split them.
      do j = 1, size(faceToSplitIdxs)
        ! Initialise the arrays
        if (allocated(PosNewFaceVertexIdxs)) deallocate(PosNewFaceVertexIdxs)
        if (allocated(NegNewFaceVertexIdxs)) deallocate(NegNewFaceVertexIdxs)
        allocate(PosNewFaceVertexIdxs(0))
        allocate(NegNewFaceVertexIdxs(0))

        ! Retrieve the indices of the vertices in the current face.
        faceVertexIdxs = faces % getFaceVertexIdxs(faceToSplitIdxs(j))

        ! Calculate the dot product between the cut plane normal and a vector
        ! from the currect vertex to one of vertices in the reflex edge 
        if (allocated(dotProducts)) deallocate(dotProducts)
        allocate(dotProducts(size(faceVertexIdxs)))
        do k = 1, size(faceVertexIdxs)
          if (any(edgeVertexIdxs == faceVertexIdxs(k))) then
            dotProducts(k) = ZERO
            cycle

          end if
          dotProducts(k) = dot_product(normalDifference, CoorVertexNotch - vertices % getVertexCoordinates(faceVertexIdxs(k)))

        end do

        ! Start calculation for the first vertex
        if (dotProducts(1) > ZERO) then
          call append(PosNewFaceVertexIdxs, faceVertexIdxs(1))
          prev = dotProducts(1)
        elseif (dotProducts(1) < ZERO) then
          call append(NegNewFaceVertexIdxs, faceVertexIdxs(1))
          prev = dotProducts(1)
        else
          call append(NegNewFaceVertexIdxs, faceVertexIdxs(1))
          call append(PosNewFaceVertexIdxs, faceVertexIdxs(1))
          prev = dotProducts(2)

        end if

        ! Loop through from the second to the last vertex
        do k = 2, size(dotProducts)
          if (dotproducts(k) == ZERO) then
            call append(NegNewFaceVertexIdxs, faceVertexIdxs(k))
            call append(PosNewFaceVertexIdxs, faceVertexIdxs(k))
            prev = -1 * prev

          else
            if (dotproducts(k)*prev < ZERO) then
              cutVertexIdx = edges % getEdgeCutVertexIdx(vertices % findCommonEdgeIdx(faceVertexIdxs(k), &
                             faceVertexIdxs(k-1)))
              call append(NegNewFaceVertexIdxs, cutVertexIdx)
              call append(PosNewFaceVertexIdxs, cutVertexIdx)
              prev = -1 * prev

            end if

            if (dotproducts(k) > ZERO) then
              call append(PosNewFaceVertexIdxs, faceVertexIdxs(k))
            else
              call append(NegNewFaceVertexIdxs, faceVertexIdxs(k))
            end if

          end if 

        end do

        if (prev * dotproducts(1) < ZERO) then
          cutVertexIdx = edges % getEdgeCutVertexIdx(vertices % findCommonEdgeIdx(faceVertexIdxs(1), &
                             faceVertexIdxs(size(faceVertexIdxs))))
          call append(NegNewFaceVertexIdxs, cutVertexIdx)
          call append(PosNewFaceVertexIdxs, cutVertexIdx)

        end if

        ! Create the new faces in the newFaces shelf.
        nFaces = nFaces + 1
        type = 'Polygon'
        if (size(NegNewFaceVertexIdxs) == 3) type = 'Triangle'
        call newFaces % expandShelf(1)
        call newFaces % buildFace(nFaces, faceToSplitIdxs(j), faces % getFaceIsBoundary(faceToSplitIdxs(j)), &
                                  NegNewFaceVertexIdxs, newVertices, type)

        ! Update mesh connectivity information.
        do k = 1, size(NegNewFaceVertexIdxs)
          call newVertices % addFaceIdxToVertex(NegNewFaceVertexIdxs(k), nFaces)

        end do

        nFaces = nFaces + 1
        type = 'Polygon'
        if (size(PosNewFaceVertexIdxs) == 3) type = 'Triangle'
        call newFaces % expandShelf(1)
        call newFaces % buildFace(nFaces, faceToSplitIdxs(j), faces % getFaceIsBoundary(faceToSplitIdxs(j)), &
                                  PosNewFaceVertexIdxs, newVertices, type)

        ! Update mesh connectivity information.
        do k = 1, size(PosNewFaceVertexIdxs)
          call newVertices % addFaceIdxToVertex(PosNewFaceVertexIdxs(k), nFaces)

        end do

      end do

      ! Create a new face corresponding to the cut plane.
      nFaces = nFaces + 1
      allocate(newFaceVertexIdxs(nVertices - nInitialVertices + 2))
      do j = 1, size(newFaceVertexIdxs) - 2
        newFaceVertexIdxs(j) = nInitialVertices + j

      end do
      newFaceVertexIdxs(size(newFaceVertexIdxs) - 1:size(newFaceVertexIdxs)) = edgeVertexIdxs
      type = 'Polygon'
      if (size(newFaceVertexIdxs) == 3) type = 'Triangle'
      call newFaces % expandShelf(1)
      call newFaces % buildFace(nFaces, 0, .false., newFaceVertexIdxs, newVertices, type)

      ! Update mesh connectivity information.
      do j = 1, size(newFaceVertexIdxs)
        call newVertices % addFaceIdxToVertex(newFaceVertexIdxs(j), nFaces)

      end do

      ! Split the element into two. Calculate the dot product between the cut face normal and a vector
      ! from the current face centroid to the centroid of the cut face.
      if (allocated(dotProducts)) deallocate(dotProducts)
      allocate(dotProducts(nFaces - 1))
      do j = 1, nFaces - 1
        dotProducts(j) = dot_product(newFaces % getFaceNormal(nFaces), &
                                     newFaces % getFaceCentroid(nFaces) - newFaces % getFaceCentroid(j))
        if (dotProducts(j) > ZERO) call append(positiveIdxs, j)
        if (dotProducts(j) < ZERO) call append(negativeIdxs, j)

      end do
      call append(positiveIdxs, nFaces)
      call append(negativeIdxs, -nFaces)

      ! Create a new element, first for the positive indices (owner of cut face).
      allocate(convexElements(2))
      nElements = nElements + 1
      if (allocated(newElementVertexIdxs)) deallocate(newElementVertexIdxs)
      do j = 1, size(positiveIdxs)
        newFaceVertexIdxs = newFaces % getFaceVertexIdxs(abs(positiveIdxs(j)))
        do k = 1, size(newFaceVertexIdxs)
          call append(newElementVertexIdxs, newFaceVertexIdxs(k), .true.)

        end do

      end do
      if (size(newElementVertexIdxs) == 4) then
        type = 'Tetrahedron'
        allocate(tetrahedron :: convexElements(nElements) % item)

      else
        type = 'Polyhedron'
        allocate(polyhedron :: convexElements(nElements) % item)

      end if
      call convexElements(nElements) % item % build(nElements, self % getIdx(), positiveIdxs, newElementVertexIdxs, &
                                                    newFaces, newVertices, type)
      
      do j = 1, size(positiveIdxs)
        call newFaces % addElementIdxToFace(abs(positiveIdxs(j)), nElements)

      end do

      do j = 1, size(newElementVertexIdxs)
        call newVertices % addElementIdxToVertex(newElementVertexIdxs(j), nElements)

      end do

      nElements = nElements + 1
      if (allocated(newElementVertexIdxs)) deallocate(newElementVertexIdxs)
      do j = 1, size(negativeIdxs)
        newFaceVertexIdxs = newFaces % getFaceVertexIdxs(abs(negativeIdxs(j)))
        do k = 1, size(newFaceVertexIdxs)
          call append(newElementVertexIdxs, newFaceVertexIdxs(k), .true.)

        end do

      end do
      if (size(newElementVertexIdxs) == 4) then
        type = 'Tetrahedron'
        allocate(tetrahedron :: convexElements(nElements) % item)

      else
        type = 'Polyhedron'
        allocate(polyhedron :: convexElements(nElements) % item)

      end if
      call convexElements(nElements) % item % build(nElements, self % getIdx(), negativeIdxs, newElementVertexIdxs, &
                                                    newFaces, newVertices, type)
      
      do j = 1, size(negativeIdxs)
        call newFaces % addElementIdxToFace(abs(negativeIdxs(j)), nElements)

      end do

      do j = 1, size(newElementVertexIdxs)
        call newVertices % addElementIdxToVertex(newElementVertexIdxs(j), nElements)

      end do

    end do

  end subroutine splitConcave

end module polyhedron_class