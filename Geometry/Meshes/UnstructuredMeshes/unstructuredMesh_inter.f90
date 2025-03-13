module unstructuredMesh_inter

  use coord_class,         only : coord
  use cellZoneShelf_class, only : cellZoneShelf
  use dictionary_class,    only : dictionary
  use edgeShelf_class,     only : edgeShelf
  use element_inter,       only : elementBox
  use elementShelf_class,  only : elementShelf
  use face_inter,          only : faceBox
  use faceShelf_class,     only : faceShelf
  use genericProcedures,   only : append, fatalError, findDifferent, numToChar
  use mesh_inter,          only : mesh, kill_super => kill
  use numPrecision
  use universalVariables
  use vertexShelf_class,   only : vertexShelf
  use kdTree_class,        only : kdTree

  implicit none
  private

  ! Extendable procedures.
  public :: distanceToBoundaryFace, distanceToNextFace, findElementAndParentIdxs, kill
  
  !! Abstract interface to group all unstructured meshes. An unstructured mesh uses a vertex -> face 
  !! -> element representation of space. Each element is composed by a set of faces which are themselves 
  !! composed by a number of vertices. Elements can be grouped together into zones. This is useful to 
  !! assign material filling to mesh elements. Local ids are assigned in the order of the cell zone 
  !! definition.
  !!
  !! Public members:
  !!   cellZones                -> Shelf that stores cell zones.
  !!   edges                    -> Shelf that stores edges.
  !!   elements                 -> Shelf that stores elements.
  !!   faces                    -> Shelf that stores faces.
  !!   vertices                 -> Shelf that stores vertices.
  !!   nVertices                -> Number of vertices in the mesh.
  !!   nFaces                   -> Number of faces in the mesh.
  !!   nEdges                   -> Number of edges in the mesh.
  !!   nElements                -> Number of elements in the mesh.
  !!   nInternalFaces           -> Number of internal faces in the mesh.
  !!   tree                     -> kd-tree used for nearest-neighbour searches and entry checks.
  !!
  !! Interface:
  !!   kill                     -> Returns to an unitialised state.
  !!   printComposition         -> Displays mesh composition to the user.
  !!   distanceToBoundaryFace   -> Checks if a particle enters the mesh and returns distance to entry 
  !!                               intersection.
  !!   distanceToNextFace       -> Returns the distance to the next mesh face.
  !!   findElementAndParentIdxs -> Returns the index of the mesh element occupied by a particle. Also
  !!                               returns the index of the parent mesh element containing the occupied
  !!                               element.
  !!
  type, public, abstract, extends(mesh) :: unstructuredMesh
    private
    integer(shortInt), public           :: nVertices = 0, nFaces = 0, nEdges = 0, &
                                           nElements = 0, nInternalFaces = 0
    type(edgeShelf), public             :: edges
    type(elementShelf), public          :: elements
    type(faceShelf), public             :: faces
    type(vertexShelf), public           :: centroids, vertices
    type(kdTree), public                :: tree, centroidTree
  contains
    ! Build procedures.
    procedure                           :: computePrimitives
    procedure(importMesh), deferred     :: importMesh
    procedure                           :: init
    procedure                           :: kill
    procedure, non_overridable          :: printComposition
    procedure                           :: setCentroidShelf
    procedure                           :: setEdgeShelf
    procedure                           :: setElementShelf
    procedure                           :: setFaceShelf
    procedure                           :: setVertexShelf
    procedure                           :: split
    procedure                           :: splitConcaveElements
    procedure                           :: splitElements
    procedure                           :: splitFaces
    ! Runtime procedures.
    procedure                           :: distanceToBoundaryFace
    procedure                           :: distanceToNextFace
    procedure                           :: findElementAndParentIdxs
    procedure                           :: findElementFromDirection
    procedure                           :: getAllCentroidCoordinates
    procedure                           :: getAllVertexCoordinates
  end type unstructuredMesh

  abstract interface

    !! Subroutine 'distanceToNextFace'
    !!
    !! Basic description:
    !!   Returns the distance to the next intersected mesh face.
    !!
    !! Arguments:
    !!   d [out]        -> Distance to the next intersected face.
    !!   coords [inout] -> Particle's coordinates.
    !!
    subroutine importMesh(self, folderPath, centroids, edges, elements, elementZones, faces, vertices, concaveElementIdxs)
      import                                                    :: unstructuredMesh, cellZoneShelf, edgeShelf, &
                                                                   elementShelf, faceShelf, shortInt, vertexShelf
      class(unstructuredMesh), intent(inout)                    :: self
      character(*), intent(in)                                  :: folderPath
      type(edgeShelf), intent(out)                              :: edges
      type(elementShelf), intent(out)                           :: elements
      type(cellZoneShelf), intent(out)                          :: elementZones
      type(faceShelf), intent(out)                              :: faces
      type(vertexShelf), intent(out)                            :: centroids, vertices
      integer(shortInt), dimension(:), allocatable, intent(out) :: concaveElementIdxs

    end subroutine importMesh

  end interface

contains

  !! Subroutine 'computePrimitives'
  !!
  !! Basic description:
  !!   Computes the number of pyramids, triangles and tetrahedra to be created during the mesh
  !!   splitting process.
  !!
  !! Detailed description:
  !!   The number of pyramids is simply given by the sum of the number of faces in each element in
  !!   the original element. The number of triangles is more complex: each pyramid created during
  !!   the splitting process also creates a number of triangles equal to the number of edges (or
  !!   vertices) in the current face. However, since all these triangles are internal they are
  !!   always shared between two pyramids; therefore, the number of triangles to be generated during
  !!   the pyramid creation process is, for a given element, equal to the sum of the number of
  !!   vertices in each of the element's face divided by two. Triangles are also created during the
  !!   splitting of the original mesh's faces: for a given face, the number of triangles to be
  !!   created is simply equal to the number of vertices in the face less two. Lastly, during the
  !!   splitting of pyramids into tetrahedra, additional internal triangles are created, given by
  !!   the number of vertices in a given pyramid's base less three. The number of tetrahedra to be
  !!   generated simply is, for a given face, the number of triangles it is decomposed into.
  !!
  !! Arguments:
  !!   nEdges [out]      -> Number of edges to be generated.
  !!   nTriangles [out]  -> Number of triangles to be generated.
  !!   nTetrahedra [out] -> Number of tetrahedra to be generated.
  !!   nVertices [out]   -> Number of vertices to be generated.
  !!
  elemental subroutine computePrimitives(self, elements, faces, nEdges, nInternalTriangles, nTetrahedra, nTriangles, nVertices)
    class(unstructuredMesh), intent(in)          :: self
    type(elementShelf), intent(in)               :: elements
    type(faceShelf), intent(in)                  :: faces
    integer(shortInt), intent(out)               :: nEdges, nInternalTriangles, nTetrahedra, nTriangles, nVertices
    integer(shortInt)                            :: i, j, nVerticesInElement, nFaces, nVerticesInFace, &
                                                    absFaceIdx
    integer(shortInt), dimension(:), allocatable :: faceIdxs

    ! Initialise nEdges = 0, nInternalTriangles = 0, nTetrahedra = 0, nTriangles = 0 and nVertices = 0.
    nEdges = 0
    nInternalTriangles = 0
    nTetrahedra = 0
    nTriangles = 0
    nVertices = 0

    ! Loop through all elements.
    do i = 1, self % nElements
      ! Retrieve the number of vertices and indices of the faces in the current element.
      nVerticesInElement = size(elements % getElementVertexIdxs(i))
      faceIdxs = elements % getElementFaceIdxs(i)
      
      ! Check if the current element is already a tetrahedron. If yes, increment nTetrahedra by 1
      ! and nTriangles by the number of triangles owned by the tetrahedron then cycle.
      if (nVerticesInElement == 4) then
        nTetrahedra = nTetrahedra + 1
        nTriangles = nTriangles + count(faceIdxs > 0)

        do j = 1, 4
          if (faceIdxs(j) > 0) then
            if (.not. faces % getFaceIsBoundary(faceIdxs(j))) nInternalTriangles = nInternalTriangles + 1

          end if

        end do
        cycle
      
      end if

      ! If the current element is not a tetrahedron it will be split from its centroid so we need
      ! to add the current element's centroid to the list of vertices.
      nVertices = nVertices + 1

      ! Increment nEdges.
      nEdges = nEdges + nVerticesInElement

      ! Compute the number of faces in the current element.
      nFaces = size(faceIdxs)
      
      ! Initialise nVertices and loop through all faces.
      nVerticesInElement = 0
      do j = 1, nFaces
        absFaceIdx = abs(faceIdxs(j))
        ! Retrieve the number of vertices in the current face and increase the total 
        ! number of vertices by the number of vertices in the current face.
        nVerticesInFace = size(faces % getFaceVertexIdxs(absFaceIdx))
        nVerticesInElement = nVerticesInElement + nVerticesInFace
        
        ! Increase the number of triangles corresponding to new internal faces by nVerticesInFace - 3.
        nTriangles = nTriangles + nVerticesInFace - 3
        nInternalTriangles = nInternalTriangles + nVerticesInFace - 3

        ! If the element owns the current face, increase the number of triangles by nVerticesInFace - 2.
        if (faceIdxs(j) > 0) then
          nEdges = nEdges + nVerticesInFace - 3
          nTriangles = nTriangles + nVerticesInFace - 2
          if (.not. faces % getFaceIsBoundary(faceIdxs(j))) nInternalTriangles = nInternalTriangles + nVerticesInFace - 2

        end if
        
        ! There will be as many tetrahedra as the number of triangles in each face, which is given
        ! by nVerticesInFace - 2.
        nTetrahedra = nTetrahedra + nVerticesInFace - 2
      
      end do
      
      ! The number of pyramids' faces is given by half the total number of vertices.
      nTriangles = nTriangles + nVerticesInElement / 2
      nInternalTriangles = nInternalTriangles + nVerticesInElement / 2

    end do

  end subroutine computePrimitives

  !! Subroutine 'distanceToBoundaryFace'
  !!
  !! Basic description:
  !!   Returns the distance to the mesh boundary face intersected by a particle's path. Also returns the index
  !!   of the parent element containing the intersected boundary face.
  !!
  !! See mesh_inter for details.
  !!
  elemental subroutine distanceToBoundaryFace(self, d, coords, parentIdx)
    class(unstructuredMesh), intent(in)          :: self
    real(defReal), intent(out)                   :: d
    type(coord), intent(inout)                   :: coords
    integer(shortInt), intent(out)               :: parentIdx
    integer(shortInt)                            :: edgeIdx, vertexIdx
    integer(shortInt), dimension(:), allocatable :: testFaceIdxs, elementIdxs
    
    ! Initialise parentIdx = 0, edgeIdx = 0 and vertexIdx = 0 then search the tree for the intersected boundary face.
    parentIdx = 0
    edgeIdx = 0
    vertexIdx = 0
    call self % tree % findIntersectedFace(self % vertices, self % faces, d, coords, edgeIdx, vertexIdx)

    ! If edgeIdx > 0 or vertexIdx > 0 then the particle intersects a boundary face through and edge or a vertex. In this
    ! case assign element occupation from particle direction.
    if (edgeIdx > 0 .or. vertexIdx > 0) then
      if (edgeIdx > 0) then
        testFaceIdxs = self % edges % getEdgeFaceIdxs(edgeIdx)
        elementIdxs = self % edges % getEdgeElementIdxs(edgeIdx)

      elseif (vertexIdx > 0) then
        testFaceIdxs = self % vertices % getVertexFaceIdxs(vertexIdx)
        elementIdxs = self % vertices % getVertexElementIdxs(vertexIdx)

      end if
      coords % elementIdx = self % findElementFromDirection(coords % dir, elementIdxs, testFaceIdxs)
      
      ! If the particle's elementIdx is 0 (e.g., particle enters the mesh via a boundary vertex but still points outside)
      ! then update d = INF and return.
      if (coords % elementIdx == 0) then
        d = INF
        return

      end if

    end if

    ! Update parentIdx only if particle enters the mesh.
    if (coords % elementIdx > 0) parentIdx = self % elements % getElementParentIdx(coords % elementIdx)

  end subroutine distanceToBoundaryFace

  !! Subroutine 'distanceToNextFace'
  !!
  !! Basic description:
  !!   Returns the distance to the next face intersected by the particle's path. Returns INF if the particle
  !!   does not intersect any face (i.e., if its path is entirely contained in the element the particle 
  !!   currently is). Algorithm adapted from Macpherson, et al. (2009). DOI: 10.1002/cnm.1128.
  !!
  !! See mesh_inter for details.
  !!
  elemental subroutine distanceToNextFace(self, d, coords)
    class(unstructuredMesh), intent(in)          :: self
    real(defReal), intent(out)                   :: d
    type(coord), intent(inout)                   :: coords
    real(defReal), dimension(3)                  :: r, rEnd
    integer(shortInt), dimension(:), allocatable :: potentialFaces, faceToElements
    integer(shortInt)                            :: elementIdx, intersectedFaceIdx
    real(defReal)                                :: lambda
    
    ! Initialise d = INF, retrieve the element currently occupied by the particle and compute potential 
    ! face intersections.
    d = INF
    elementIdx = coords % elementIdx
    rEnd = coords % rEnd
    potentialFaces = self % elements % computePotentialFaceIdxs(elementIdx, rEnd, self % faces)

    ! If no potential intersections are detected return early.
    if (size(potentialFaces) == 0) return
    
    ! If reached here, compute which face is actually intersected and update d.
    r = coords % r
    call self % elements % computeFaceIntersection(elementIdx, r, rEnd, potentialFaces, self % faces, &
                                                   intersectedFaceIdx, lambda)
    d = norm2(min(ONE, max(ZERO, lambda)) * (rEnd - r))
    
    ! If the intersected face is a boundary face then the particle is leaving the mesh.
    if (self % faces % getFaceIsBoundary(intersectedFaceIdx)) then
      coords % elementIdx = 0
      coords % localId = 1
      return

    end if

    ! Else, retrieve the elements sharing the intersected face from mesh connectivity then
    ! update elementIdx and localId.
    faceToElements = self % faces % getFaceElementIdxs(intersectedFaceIdx)
    coords % elementIdx = findDifferent(faceToElements, elementIdx)
    coords % localId = self % findElementZoneIdx(self % elements % getElementParentIdx(coords % elementIdx))

  end subroutine distanceToNextFace

  !! Subroutine 'findElementAndParentIdxs'
  !!
  !! Basic description:
  !!   Returns the index of the mesh element occupied by a particle. Also returns the index of the parent mesh
  !!   element containing the occupied element.
  !!
  !! See mesh_inter for details.
  !!
  pure subroutine findElementAndParentIdxs(self, r, u, elementIdx, parentIdx)
    class(unstructuredMesh), intent(in)          :: self
    real(defReal), dimension(3), intent(in)      :: r, u
    integer(shortInt), intent(out)               :: elementIdx, parentIdx
    real(defReal), dimension(6)                  :: boundingBox
    integer(shortInt), dimension(:), allocatable :: potentialElements, faceToElements, zeroDotProductFaceIdxs, &
                                                    testFaceIdxs, elementFaceIdxs, elementVertexIdxs, vertexFaceIdxs
    integer(shortInt)                            :: i, nearestVertexIdx, failedFaceIdx, commonEdgeIdx, commonVertexIdx
    logical(defBool), dimension(self % nElements) :: isVisited
    
    ! Initialise elementIdx = 0 and parentIdx = 0. Retrieve the mesh's bounding box. If the particle is outside the
    ! bounding box we can return early.
    elementIdx = 0
    parentIdx = 0
    boundingBox = self % getBoundingBox()
    do i = 1, 3
      if (r(i) <= boundingBox(i) .or. boundingBox(i + 3) <= r(i)) return

    end do

    ! If the point is inside the bounding box then we need to determine if the point is inside a mesh element.
    ! First find the vertex which is nearest to the coordinates. Query the tree with the lowest number of points.
    if (self % nElements < self % nVertices) then
      nearestVertexIdx = self % centroidTree % findNearestVertex(r)

    else
      nearestVertexIdx = self % tree % findNearestVertex(r)

    end if

    ! Retrieve potential elements occupied by the particle from mesh connectivity information.
    potentialElements = self % vertices % getVertexElementIdxs(nearestVertexIdx)

    ! Initialise isVisited = .false. then do a first pass over all elements sharing the nearest vertex.
    isVisited = .false.
    outerLoop: do i = 1, size(potentialElements)
      elementIdx = potentialElements(i)
      searchLoop: do
        if (isVisited(elementIdx)) cycle outerLoop
        call self % elements % testForInclusion(elementIdx, r, self % faces, failedFaceIdx, zeroDotProductFaceIdxs)
        isVisited(elementIdx) = .true.

        if (failedFaceIdx > 0) then
          if (self % faces % getFaceIsBoundary(failedFaceIdx)) cycle outerLoop
          elementIdx = findDifferent(self % faces % getFaceElementIdxs(failedFaceIdx), elementIdx)
          cycle searchLoop

        end if

        ! If there are faces on which the particle lies we need to employ some more specific
        ! procedures to correctly determine which element is actually occupied.
        if (allocated(zeroDotProductFaceIdxs)) then
          select case (size(zeroDotProductFaceIdxs))
            ! Particle is on a element's face.
            case(1)
              potentialElements = self % faces % getFaceElementIdxs(zeroDotProductFaceIdxs(1))
              testFaceIdxs = zeroDotProductFaceIdxs
            ! Particle is on a element's edge.
            case(2)
              ! Find the common edge and retrieve elements sharing this edge.
              commonEdgeIdx = self % faces % findCommonEdgeIdx(zeroDotProductFaceIdxs(1), zeroDotProductFaceIdxs(2))
              potentialElements = self % edges % getEdgeElementIdxs(commonEdgeIdx)
              testFaceIdxs = self % edges % getEdgeFaceIdxs(commonEdgeIdx)
            ! Particle is on a element's vertex.
            case default
              ! Find the common vertex and retrieve all elements sharing this vertex.
              commonVertexIdx = self % faces % findCommonVertexIdx(zeroDotProductFaceIdxs)
              potentialElements = self % vertices % getVertexElementIdxs(commonVertexIdx)
              testFaceIdxs = self % vertices % getVertexFaceIdxs(commonVertexIdx)
  
          end select
          elementIdx = self % findElementFromDirection(u, potentialElements, testFaceIdxs)
          if (elementIdx > 0) parentIdx = self % elements % getElementParentIdx(elementIdx)
          return
        
        end if

        parentIdx = self % elements % getElementParentIdx(elementIdx)
        return

      end do searchLoop

    end do outerLoop

    ! If host element has not been found yet, do a search on all remaining elements.
    do i = 1, self % nElements
      if (isVisited(i)) cycle
      elementIdx = i
      call self % elements % testForInclusion(elementIdx, r, self % faces, failedFaceIdx, zeroDotProductFaceIdxs)
      isVisited(elementIdx) = .true.
      
      if (failedFaceIdx > 0) cycle

      ! If there are faces on which the particle lies we need to employ some more specific
      ! procedures to correctly determine which element is actually occupied.
      if (allocated(zeroDotProductFaceIdxs)) then
        select case (size(zeroDotProductFaceIdxs))
          ! Particle is on a element's face.
          case(1)
            potentialElements = self % faces % getFaceElementIdxs(zeroDotProductFaceIdxs(1))
            testFaceIdxs = zeroDotProductFaceIdxs
          ! Particle is on a element's edge.
          case(2)
            ! Find the common edge and retrieve elements sharing this edge.
            commonEdgeIdx = self % faces % findCommonEdgeIdx(zeroDotProductFaceIdxs(1), zeroDotProductFaceIdxs(2))
            potentialElements = self % edges % getEdgeElementIdxs(commonEdgeIdx)
            testFaceIdxs = self % edges % getEdgeFaceIdxs(commonEdgeIdx)
          ! Particle is on a element's vertex.
          case default
            ! Find the common vertex and retrieve all elements sharing this vertex.
            commonVertexIdx = self % faces % findCommonVertexIdx(zeroDotProductFaceIdxs)
            potentialElements = self % vertices % getVertexElementIdxs(commonVertexIdx)
            testFaceIdxs = self % vertices % getVertexFaceIdxs(commonVertexIdx)

        end select
        elementIdx = self % findElementFromDirection(u, potentialElements, testFaceIdxs)
        if (elementIdx > 0) parentIdx = self % elements % getElementParentIdx(elementIdx)
        return
      
      end if

      parentIdx = self % elements % getElementParentIdx(elementIdx)
      return

    end do

    ! If reached here, particle is outside the mesh.
    elementIdx = 0

  end subroutine findElementAndParentIdxs

  !! Function 'findElementFromDirection'
  !!
  !! Basic description:
  !!   Returns the index of the element occupied by a particle in case the particle lies on one or more
  !!   element face(s).
  !!
  !! Detailed description:
  !!   When the particle lies on one or more element face(s), it is technically in more than one element
  !!   at once. In this case element occupation can be assigned by using the particle's direction vector.
  !!   The idea is to take all the faces on which the particle lies and perform the dot product of the
  !!   particle's direction with each of the face's (signed) normal vector. The element actually occupied
  !!   by the partice is then that for which the number of negative dot products is the greatest.
  !!
  !! Arguments:
  !!   u [in]            -> Particle's direction.
  !!   elementIdxs [in]  -> Indices of the elements to be tested.
  !!   testFaceIdxs [in] -> Indices of the faces to be tested.
  !!
  !! Result:
  !!   elementIdx        -> Index of the element occupied by the particle.
  !!
  pure function findElementFromDirection(self, u, elementIdxs, testFaceIdxs) result(elementIdx)
    class(unstructuredMesh), intent(in)          :: self
    real(defReal), dimension(3), intent(in)      :: u
    integer(shortInt), dimension(:), intent(in)  :: elementIdxs, testFaceIdxs
    integer(shortInt)                            :: elementIdx, i, j, count, maxCount, maxIdx
    integer(shortInt), dimension(:), allocatable :: faceIdxs, absFaceIdxs
    real(defReal), dimension(3)                  :: normal

    ! Initialise elementIdx = 0 and maxCount = 0 then loop over all elements.
    elementIdx = 0
    maxCount = 0
    do i = 1, size(elementIdxs)
      ! Retrieve the faces making the current element and convert them to absolute indices.
      faceIdxs = self % elements % getElementFaceIdxs(elementIdxs(i))
      absFaceIdxs = abs(faceIdxs)
      
      ! Initialise count = 0 then loop through all faces in the current element.
      count = 0
      do j = 1, size(faceIdxs)

        ! Cycle to the next face if it is not a face on which the particle is.
        if (.not. any(testFaceIdxs == absFaceIdxs(j))) cycle

        ! Retrieve the current face's signed normal vector.
        normal = self % faces % getFaceNormal(faceIdxs(j))

        ! If the dot product is negative increment count and cycle to the next face.
        if (dot_product(u, normal) <= ZERO) then
          count = count + 1
          cycle

        end if
        
        ! If the test fails and the current face is a boundary face, the particle is outside the mesh
        ! and we can return early.
        if (self % faces % getFaceIsBoundary(absFaceIdxs(j))) return

      end do

      ! If count > maxCount, update maxCount and maxIdx.
      if (count > maxCount) then
        maxCount = count
        maxIdx = i

      end if

    end do

    ! If reached here, update elementIdx.
    elementIdx = elementIdxs(maxIdx)

  end function findElementFromDirection

  !! Function 'getAllCentroidCoordinates'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of all the centroids in the mesh.
  !!
  !! Result:
  !!   coords -> 3-D coordinates of all the centroids in the mesh.
  !!
  pure function getAllCentroidCoordinates(self) result(coords)
    class(unstructuredMesh), intent(in)           :: self
    real(defReal), dimension(self % nElements, 3) :: coords

    coords = self % centroids % getAllCoordinates()

  end function getAllCentroidCoordinates

  !! Function 'getAllVertexCoordinates'
  !!
  !! Basic description:
  !!   Returns the 3-D coordinates of all the vertices in the mesh.
  !!
  !! Result:
  !!   coords -> 3-D coordinates of all the vertices in the mesh.
  !!
  pure function getAllVertexCoordinates(self) result(coords)
    class(unstructuredMesh), intent(in)           :: self
    real(defReal), dimension(self % nVertices, 3) :: coords

    coords = self % vertices % getAllCoordinates()

  end function getAllVertexCoordinates

  !!
  !!
  !!
  subroutine init(self, folderPath, dict)
    class(unstructuredMesh), intent(inout)       :: self
    character(*), intent(in)                     :: folderPath
    class(dictionary), intent(in)                :: dict
    type(edgeShelf)                              :: edges, newEdges
    type(elementShelf)                           :: elements, newElements
    type(cellZoneShelf)                          :: elementZones
    type(faceShelf)                              :: faces, newFaces
    type(vertexShelf)                            :: centroids, newCentroids, newVertices, vertices
    integer(shortInt), dimension(:), allocatable :: concaveElementIdxs
    integer(shortInt)                            :: i, nConcaveElements
    logical(defBool)                             :: triangulate
    character(*), parameter                      :: Here = 'init (unstructuredMesh_inter.f90)'

    ! Set up base components.
    call self % setupBase(dict)
    
    ! Import mesh from files.
    call self % importMesh(folderPath, centroids, edges, elements, elementZones, faces, vertices, concaveElementIdxs)

    ! Check if any elements are concave.
    nConcaveElements = size(concaveElementIdxs)
    if (nConcaveElements > 0) then 
      call self % splitConcaveElements(concaveElementIdxs, edges, elements, faces, vertices, &
                                       newCentroids, newEdges, newElements, newFaces, newVertices)
      call fatalError(Here, 'Temporary error message.')
      
    end if

    ! Check if triangulation was requested.
    call dict % getOrDefault(triangulate, 'triangulate', .false.)
    if (triangulate) then
      call self % split(edges, elements, faces, vertices, newCentroids, newEdges, newElements, newFaces, newVertices)
      call self % setCentroidShelf(newCentroids)
      call self % setEdgeShelf(newEdges)
      call self % setElementShelf(newElements)
      call self % setFaceShelf(newFaces)
      call self % setVertexShelf(newVertices)

    else
      call self % setCentroidShelf(centroids)
      call self % setEdgeShelf(edges)
      call self % setElementShelf(elements)
      call self % setFaceShelf(faces)
      call self % setVertexShelf(vertices)

    end if

    ! Set elements zones and initialise kd-tree for the mesh.
    call self % setElementZones(elementZones)
    call self % tree % init(self % getAllVertexCoordinates(), .true.)
    call self % centroidTree % init(self % getAllCentroidCoordinates(), .true.)

  end subroutine init

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  elemental subroutine kill(self)
    class(unstructuredMesh), intent(inout) :: self

    ! Superclass.
    call kill_super(self)
    
    ! Local.
    self % nVertices = 0
    self % nFaces = 0
    self % nInternalFaces = 0
    self % nElements = 0
    self % nEdges = 0
    call self % centroids % kill()
    call self % edges % kill()
    call self % elements % kill()
    call self % faces % kill()
    call self % centroidTree % kill()
    call self % tree % kill()
    call self % vertices % kill()

  end subroutine kill

  !! Subroutine 'printComposition'
  !!
  !! Basic description:
  !!   Prints the initial polyhedral composition of the mesh.
  !!
  !! Arguments:
  !!   nTetrahedra [out] -> Number of tetrahedra in the mesh.
  !!
  subroutine printComposition(self, nTetrahedra)
    class(unstructuredMesh), intent(in) :: self
    integer(shortInt), intent(out)      :: nTetrahedra
    integer(shortInt)                   :: nFaces, nPentahedra, nHexahedra, nOthers, i
    
    ! Initialise the numbers of various polyhedra to zero.
    nTetrahedra = 0
    nPentahedra = 0
    nHexahedra = 0
    nOthers = 0
    
    ! Loop over all elements in the mesh.
    do i = 1, self % nElements
      ! Retrieve the number of faces in the current element and increment specific polyhedra
      ! accordingly.
      nFaces = size(self % elements % getElementFaceIdxs(i))
      select case (nFaces)
        case (4)
          nTetrahedra = nTetrahedra + 1
        case (5)
          nPentahedra = nPentahedra + 1
        case (6)
          nHexahedra = nHexahedra + 1
        case default
          nOthers = nOthers + 1

      end select

    end do
    
    ! Print to screen.
    print *, 'Displaying unstructured mesh composition:'
    print *, '  Number of tetrahedra     : '//numToChar(nTetrahedra)//'.'
    print *, '  Number of pentahedra     : '//numToChar(nPentahedra)//'.'
    print *, '  Number of hexahedra      : '//numToChar(nHexahedra)//'.'
    print *, '  Number of other polyhedra: '//numToChar(nOthers)//'.'

  end subroutine printComposition

  !! Subroutine 'setCentroidShelf'
  !!
  !! Basic description:
  !!   Sets the centroidShelf of the unstructured mesh.
  !!
  !! Arguments:
  !!   centroids [in] -> A vertexShelf.
  !!
  elemental subroutine setCentroidShelf(self, centroids)
    class(unstructuredMesh), intent(inout) :: self
    type(vertexShelf), intent(in)          :: centroids

    self % centroids = centroids

  end subroutine setCentroidShelf

  !! Subroutine 'setEdgeShelf'
  !!
  !! Basic description:
  !!   Sets the edgeShelf of the unstructured mesh.
  !!
  !! Arguments:
  !!   edges [in] -> An edgeShelf.
  !!
  elemental subroutine setEdgeShelf(self, edges)
    class(unstructuredMesh), intent(inout) :: self
    type(edgeShelf), intent(in)            :: edges

    self % edges = edges

  end subroutine setEdgeShelf

  !! Subroutine 'setElementShelf'
  !!
  !! Basic description:
  !!   Sets the elementShelf of the unstructured mesh.
  !!
  !! Arguments:
  !!   elements [in] -> An elementShelf.
  !!
  elemental subroutine setElementShelf(self, elements)
    class(unstructuredMesh), intent(inout) :: self
    type(elementShelf), intent(in)         :: elements

    self % elements = elements

  end subroutine setElementShelf

  !! Subroutine 'setFaceShelf'
  !!
  !! Basic description:
  !!   Sets the faceShelf of the unstructured mesh.
  !!
  !! Arguments:
  !!   faces [in] -> A faceShelf.
  !!
  elemental subroutine setFaceShelf(self, faces)
    class(unstructuredMesh), intent(inout) :: self
    type(faceShelf), intent(in)            :: faces

    self % faces = faces

  end subroutine setFaceShelf

  !! Subroutine 'setVertexShelf'
  !!
  !! Basic description:
  !!   Sets the vertexShelf of the unstructured mesh.
  !!
  !! Arguments:
  !!   vertices [in] -> A vertexShelf.
  !!
  elemental subroutine setVertexShelf(self, vertices)
    class(unstructuredMesh), intent(inout) :: self
    type(vertexShelf), intent(in)          :: vertices

    self % vertices = vertices

  end subroutine setVertexShelf

!! Subroutine 'split'
  !!
  !! Basic description:
  !!   Splits a mesh into tetrahedral elements. If a given element is already a tetrahedron it is
  !!   not split but simply added to the shelf of tetrahedra in the mesh.
  !!
  !! Detailed description:
  !!   'split' starts by computing the number of pyramids, tetrahedra and triangles that will be
  !!   generated in the resulting mesh. Then, each element is split into a set of pyramids, whose
  !!   bases are each of the element's face and whose (common) apex is the element's centroid. This
  !!   apex is also appended to the list of vertices in the mesh in the process. Once this is done,
  !!   each face in the original mesh is subdivided into triangles. Lastly, each pyramid previously
  !!   created is further split into tetrahedra.
  !!
  !! Arguments:
  !!   lastVertexIdx [out] -> Index of the last vertex in the resulting mesh.
  !!
  subroutine split(self, edges, elements, faces, vertices, newCentroids, newEdges, newElements, newFaces, newVertices)
    class(unstructuredMesh), intent(inout)    :: self
    type(edgeShelf), intent(in)               :: edges
    type(elementShelf), intent(inout)         :: elements
    type(faceShelf), intent(inout)            :: faces
    type(vertexShelf), intent(in)             :: vertices
    type(edgeShelf), intent(out)              :: newEdges
    type(elementShelf), intent(out)           :: newElements
    type(faceShelf), intent(out)              :: newFaces
    type(vertexShelf), intent(out)            :: newCentroids, newVertices
    integer(shortInt)                         :: i, j, nEdges, nInternalTriangles, nNewEdges, nTetrahedra, nTriangles, &
                                                 nVertices, nNewVertices, lastEdgeIdx, lastFaceIdx, lastElementIdx, lastVertexIdx
    integer(shortInt), dimension(:), allocatable :: edgeIdxs
    type(elementBox), dimension(:), allocatable  :: tetrahedra
    type(faceBox), dimension(:), allocatable     :: triangles
    
    ! Retrieve sizes of the original shelves.
    nEdges = self % nEdges
    nVertices = self % nVertices
    
    ! Compute the number of edges, pyramids, triangles and tetrahedra to be created and
    ! allocate memory to the corresponding structures.
    call self % computePrimitives(elements, faces, nNewEdges, nInternalTriangles, nTetrahedra, nTriangles, nNewVertices)
    
    ! Allocate memory in the new shelves.
    call newCentroids % allocateShelf(nTetrahedra)
    call newEdges % allocateShelf(nEdges + nNewEdges)
    call newElements % allocateShelf(nTetrahedra)
    call newFaces % allocateShelf(nTriangles)
    call newVertices % allocateShelf(nVertices + nNewVertices)

    ! Copy original edges and vertices into the new shelves.
    do i = 1, nEdges
      call newEdges % initEdge(i, edges % getEdgeVertexIdxs(i))

    end do

    call newVertices % setExtremalCoordinates(vertices % getExtremalCoordinates())
    call newVertices % setOffset(vertices % getOffset())
    do i = 1, nVertices
      call newVertices % initVertex(i, vertices % getVertexCoordinates(i))
      edgeIdxs = vertices % getVertexEdgeIdxs(i)

      do j = 1, size(edgeIdxs)
        call newVertices % addEdgeIdxToVertex(i, edgeIdxs(j))

      end do

    end do

    ! Initialise new triangles and tetrahedra to be generated.
    allocate(triangles(nTriangles))
    allocate(tetrahedra(nTetrahedra))
    
    ! Initialise lastVertexIdx, lastPyramidIdx, lastTetrahedronIdx and lastTriangleIdx then
    ! split all elements into pyramids and all pyramids into tetrahedra.
    lastEdgeIdx = nEdges
    lastElementIdx = 0
    lastFaceIdx = 0
    lastVertexIdx = nVertices
    call self % splitFaces(faces, newEdges, newFaces, newVertices, lastEdgeIdx, lastFaceIdx, triangles)
    call self % splitElements(elements, faces, lastEdgeIdx, lastElementIdx, lastFaceIdx, lastVertexIdx, &
                              newCentroids, newEdges, newElements, newFaces, newVertices, tetrahedra, triangles)

    ! Update the number of edges, faces, elements and vertices in the mesh.
    self % nEdges = nEdges + nNewEdges
    self % nElements = nTetrahedra
    self % nFaces = nTriangles
    self % nInternalFaces = nInternalTriangles
    self % nVertices = nVertices + nNewVertices

  end subroutine split

  !!
  !!
  !!
  subroutine splitConcaveElements(self, concaveElementIdxs, edges, elements, faces, vertices, newCentroids, &
                                  newEdges, newElements, newFaces, newVertices)
    class(unstructuredMesh), intent(in)         :: self
    integer(shortInt), dimension(:), intent(in) :: concaveElementIdxs
    type(edgeShelf), intent(inout)              :: edges, newEdges
    type(elementShelf), intent(inout)           :: elements, newElements
    type(faceShelf), intent(inout)              :: faces, newFaces
    type(vertexShelf), intent(inout)            :: vertices, newCentroids, newVertices
    integer(shortInt)                           :: i
    type(elementBox), dimension(:), allocatable :: convexElements

    ! Loop through all concave elements.
    do i = 1, size(concaveElementIdxs)
      ! Build all concave notches in current concave element.
      call elements % buildElementNotches(concaveElementIdxs(i), edges, faces, vertices)

      ! Now split current element.
      call elements % splitConcave(concaveElementIdxs(i), edges, faces, vertices, newEdges, convexElements, &
                                   newFaces, newVertices)

    end do

  end subroutine splitConcaveElements

  !! Subroutine 'splitElements'
  !!
  !! Basic description:
  !!   Splits all elements in the original mesh into pyramids. If a given element is already a
  !!   tetrahedron it is not split but simply appended to the list of existing tetrahedra.
  !!
  !! Arguments:
  !!   lastEdgeIdx [inout]        -> Index of the last edge in the mesh.
  !!   lastPyramidIdx [inout]     -> Index of the last pyramid in the mesh.
  !!   lastTetrahedronIdx [inout] -> Index of the last tetrahedron in the mesh.
  !!   lastTriangleIdx [inout]    -> Index of the last triangle in the mesh.
  !!   lastVertexIdx [inout]      -> Index of the last vertex in the mesh.
  !!
  subroutine splitElements(self, elements, faces, lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx, &
                           newCentroids, newEdges, newElements, newFaces, newVertices, tetrahedra, triangles)
    class(unstructuredMesh), intent(inout)        :: self
    type(elementShelf), intent(inout)             :: elements, newElements
    type(faceShelf), intent(inout)                :: faces, newFaces
    integer(shortInt), intent(inout)              :: lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx
    type(edgeShelf), intent(inout)                :: newEdges
    type(vertexShelf), intent(inout)              :: newCentroids, newVertices
    type(elementBox), dimension(:), intent(inout) :: tetrahedra
    type(faceBox), dimension(:), intent(inout)    :: triangles
    integer(shortInt)                             :: i, initialElementIdx, j
                                        
    ! Loop through all original elements and split them into tetrahedra.
    do i = 1, self % nElements
      ! Initialise initialElementIdx then split the current element.
      initialElementIdx = lastNewElementIdx + 1
      call elements % splitElement(i, faces, lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx, &
                                   newEdges, newFaces, newVertices, tetrahedra, triangles)

      ! Set all new tetrahedra.
      do j = initialElementIdx, lastNewElementIdx
        call newElements % addElement(j, tetrahedra(j))
        call newCentroids % initVertex(j, newElements % getElementCentroid(j))
        call newCentroids % addElementIdxToVertex(j, j)

      end do
      
    end do

  end subroutine splitElements

  subroutine splitFaces(self, faces, newEdges, newFaces, newVertices, lastNewEdgeIdx, lastNewFaceIdx, triangles)
    class(unstructuredMesh), intent(inout)     :: self
    type(faceShelf), intent(inout)             :: faces, newFaces
    type(edgeShelf), intent(inout)             :: newEdges
    type(vertexShelf), intent(inout)           :: newVertices
    integer(shortInt), intent(inout)           :: lastNewEdgeIdx, lastNewFaceIdx
    type(faceBox), dimension(:), intent(inout) :: triangles
    integer(shortInt)                          :: i, initialFaceIdx, j

    ! Loop through all original faces in the mesh and split them into triangles.
    do i = 1, self % nFaces
      ! Initialise initialFaceIdx then split the current face.
      initialFaceIdx = lastNewFaceIdx + 1
      call faces % splitFace(i, newEdges, newVertices, lastNewEdgeIdx, lastNewFaceIdx, triangles)

      ! Set all new triangles.
      do j = initialFaceIdx, lastNewFaceIdx
        call newFaces % addFace(j, triangles(j))

      end do

    end do

  end subroutine splitFaces

end module unstructuredMesh_inter