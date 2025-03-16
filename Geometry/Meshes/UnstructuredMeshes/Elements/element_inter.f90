module element_inter

  use edgeShelf_class,     only : edgeShelf
  use face_inter,          only : faceBox
  use faceShelf_class,     only : faceShelf
  use genericProcedures,   only : append, areEqual, computePyramidCentre, computePyramidVolume, &
                                  computeTetrahedronCentre, computeTetrahedronVolume, findCommon, &
                                  fatalError, numToChar
  use numPrecision
  use universalVariables,  only : SURF_TOL, INF, ZERO
  use vertexShelf_class,   only : vertexShelf
  
  implicit none
  private

  ! Extendable procedures.
  public :: kill
  
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
  type, public, abstract                         :: element
    private
    integer(shortInt)                            :: idx = 0, parentIdx = 0
    integer(shortInt), dimension(:), allocatable :: edgeIdxs, faceIdxs, vertexIdxs, tetrahedronIdxs, concaveFaceIdxs
    real(defReal)                                :: volume = ZERO
    real(defReal), dimension(3)                  :: centroid = ZERO
    logical(defBool)                             :: isConvex = .false.
    character(:), allocatable                    :: type
    type(notch), dimension(:), allocatable       :: notches
  contains
    ! Build procedures.
    procedure, non_overridable                   :: addEdgeIdx
    procedure, non_overridable                   :: addFaceIdx
    procedure, non_overridable                   :: addVertexIdx
    procedure, non_overridable                   :: build
    procedure(computeComponents), deferred       :: computeComponents
    procedure, non_overridable                   :: computeConvexity
    procedure, non_overridable                   :: init
    procedure, non_overridable                   :: setIdx
    procedure(split), deferred                   :: split
    procedure(splitConcave), deferred            :: splitConcave
    ! Runtime procedures.
    procedure, non_overridable                   :: buildNotches
    procedure, non_overridable                   :: computeIntersectedFace
    procedure, non_overridable                   :: computePotentialFaces
    procedure, non_overridable                   :: getCentroid
    procedure, non_overridable                   :: getEdgeIdxs
    procedure, non_overridable                   :: getFaceIdxs
    procedure, non_overridable                   :: getIdx
    procedure, non_overridable                   :: getIsConvex
    procedure, non_overridable                   :: getNotches
    procedure, non_overridable                   :: getParentIdx
    procedure, non_overridable                   :: getType
    procedure, non_overridable                   :: getVertexIdxs
    procedure, non_overridable                   :: getVolume
    procedure                                    :: kill
    procedure, non_overridable                   :: testForInclusion
  end type element

  !!
  !!
  !!
  type, public :: notch
    integer(shortInt)               :: edgeIdx
    integer(shortInt), dimension(2) :: faceIdxs
  end type notch

  !!
  !! Small, local container to store polymorphic elements in a single array.
  !!
  !! Public members:
  !!   name -> Name of the mesh.
  !!   ptr  -> Pointer to the mesh.
  !!
  type, public                  :: elementBox
    class(element), allocatable :: item
  end type

  abstract interface

    !!
    !!
    !!
    pure subroutine computeComponents(self, faceIdxs, vertexIdxs, faces, vertices, centroid, volume)
      import                                      :: element, shortInt, faceShelf, vertexShelf, defReal
      class(element), intent(inout)               :: self
      integer(shortInt), dimension(:), intent(in) :: faceIdxs, vertexIdxs
      type(faceShelf), intent(in)                 :: faces
      type(vertexShelf), intent(in)               :: vertices
      real(defReal), dimension(3), intent(out)    :: centroid
      real(defReal), intent(out)                  :: volume

    end subroutine computeComponents

    !!
    !!
    !!
    subroutine split(self, faces, lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx, &
                     newEdges, newFaces, newVertices, tetrahedra, triangles)
      import                                        :: element, edgeShelf, faceShelf, vertexShelf, shortInt, elementBox, &
                                                       faceBox
      class(element), intent(inout)                 :: self
      type(faceShelf), intent(inout)                :: faces, newFaces
      integer(shortInt), intent(inout)              :: lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx
      type(edgeShelf), intent(inout)                :: newEdges
      type(vertexShelf), intent(inout)              :: newVertices
      type(elementBox), dimension(:), intent(inout) :: tetrahedra
      type(faceBox), dimension(:), intent(inout)    :: triangles

    end subroutine split

    !!
    !!
    !!
    subroutine splitConcave(self, edges, faces, vertices, newEdges, convexElements, newFaces, newVertices)
      import                                                   :: element, edgeShelf, faceShelf, vertexShelf, elementBox
      class(element), intent(inout)                            :: self
      type(edgeShelf), intent(inout)                           :: edges, newEdges
      type(faceShelf), intent(inout)                           :: faces, newFaces
      type(vertexShelf), intent(inout)                         :: vertices, newVertices
      type(elementBox), dimension(:), allocatable, intent(out) :: convexElements

    end subroutine splitConcave

  end interface

contains

  !! Subroutine 'addEdgeIdx'
  !!
  !! Basic description:
  !!   Adds the index of an edge sharing the element.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the edge.
  !!
  elemental subroutine addEdgeIdx(self, idx)
    class(element), intent(inout)  :: self
    integer(shortInt), intent(in)  :: idx

    call append(self % edgeIdxs, idx, .true.)

  end subroutine addEdgeIdx
  
  !! Subroutine 'addFaceToElement'
  !!
  !! Basic description:
  !!   Adds the index of a face belonging to the element.
  !!
  !! Arguments:
  !!   faceIdx [in] -> Index of the face.
  !!
  elemental subroutine addFaceIdx(self, faceIdx)
    class(element), intent(inout) :: self
    integer(shortInt), intent(in) :: faceIdx
    
    call append(self % faceIdxs, faceIdx)

  end subroutine addFaceIdx
  
  !! Subroutine 'addVertexToElement'
  !!
  !! Basic description:
  !!   Adds the index of a vertex belonging to the element. Only adds it if the index is not already
  !!   present.
  !!
  !! Arguments:
  !!   vertexIdx [in] -> Index of the vertex.
  !!
  elemental subroutine addVertexIdx(self, vertexIdx)
    class(element), intent(inout) :: self
    integer(shortInt), intent(in) :: vertexIdx

    call append(self % vertexIdxs, vertexIdx, .true.)

  end subroutine addVertexIdx

  !!
  !!
  !!
  subroutine build(self, idx, parentIdx, faceIdxs, vertexIdxs, faces, vertices, type)
    class(element), intent(inout)               :: self
    integer(shortInt), intent(in)               :: idx, parentIdx
    integer(shortInt), dimension(:), intent(in) :: faceIdxs, vertexIdxs
    type(faceShelf), intent(in)                 :: faces
    type(vertexShelf), intent(in)               :: vertices
    character(*), intent(in)                    :: type
    logical(defBool)                            :: isConvex
    real(defReal), dimension(3)                 :: centroid
    real(defReal)                               :: volume

    ! Initialise volume = ZERO and centroid = ZERO
    volume = ZERO
    centroid = ZERO

    if (type == 'Tetrahedron') then
      isConvex = .true.

    else
      isConvex = self % computeConvexity(faceIdxs, vertexIdxs, faces, vertices)

    end if

    if (isConvex) call self % computeComponents(faceIdxs, vertexIdxs, faces, vertices, centroid, volume)

    ! Initialise element.
    call self % init(idx, parentIdx, faceIdxs, vertexIdxs, centroid, volume, isConvex, type)

  end subroutine build

  !!
  !!
  subroutine buildNotches(self, edges, faces, vertices)
    class(element), intent(inout)                  :: self
    type(edgeShelf), intent(in)                    :: edges
    type(faceShelf), intent(in)                    :: faces
    type(vertexShelf), intent(in)                  :: vertices
    integer(shortInt)                              :: i, j, k, faceIdx, absFaceIdx, vertexIdx, nNotches
    integer(shortInt), dimension(:), allocatable   :: faceVertexIdxs, edgeIdxs, edgeFaceIdxs, currentElementFaceIdxs, &
                                                      commonFaceIdxs
    real(defReal), dimension(3)                    :: normal, faceVertexCoords
    type(notch), dimension(:), allocatable         :: tempNotches
    integer(shortInt), dimension(2)                :: edgeVertexIdxs
    
    ! Initialise nNotches = 0
    nNotches = 0
    
    ! Get all the edge indices for this element
    edgeIdxs = self % getEdgeIdxs()

    ! Loop through all the edges in the current element
    do i = 1, size(edgeIdxs)
        ! Retrieve all the faces associated with the edge
        edgeFaceIdxs = edges % getEdgeFaceIdxs(edgeIdxs(i))
        ! Renewing for the current edge (loop variant)
        if (allocated(currentElementFaceIdxs)) deallocate(currentElementFaceIdxs)
        ! Loop through all the faces associated with the edge
        do j = 1, size(edgeFaceIdxs)
            ! If the current face is not part of the element, cycle
            if (.not. any(abs(self % getFaceIdxs()) == edgeFaceIdxs(j))) cycle
            ! O.W keep the face
            call append(currentElementFaceIdxs, edgeFaceIdxs(j))

        end do

        ! Find the notch whose two adjacent faces are defined as the problematic faces that fails convexity tests.
        ! This can be determined by testing if the size of the intersection of the two sets (current face and problematic face) is 2
        ! If notch, store the edgeIdx and faceIdxs
        commonFaceIdxs = findCommon(currentElementFaceIdxs, self % concaveFaceIdxs)
        if (size(commonFaceIdxs) == 2) then
            nNotches = nNotches + 1
            if (nNotches == 1) then
              allocate(self % notches(nNotches))

            else
              tempNotches = self % notches
              if (allocated(self % notches)) deallocate(self % notches)
              allocate(self % notches(nNotches))
              self % notches(1:nNotches - 1) = tempNotches

            end if
            self % notches(nNotches) % edgeIdx = edgeIdxs(i)
            self % notches(nNotches) % faceIdxs = commonFaceIdxs

        end if
        if (allocated(tempNotches)) deallocate(tempNotches)

    end do

  end subroutine buildNotches

  !! Function 'isConvex'
  !!
  !! Basic description:
  !!   Checks whether the element is convex.
  !!
  !! Detailed description:
  !!    Convexity is checked by taking each vertex in the a given face and creating a vector 
  !!    connecting said vertex to each vertex in the element not in the current face. If the element 
  !!    is convex then all the vertices not in the current face must lie on the same side of the 
  !!    face, hence the dot product between the current face's normal vector and the test vector 
  !!    must be negative. If at any point the dot product is found to be positive the check is 
  !!    aborted.
  !!
  !! Arguments:
  !!   vertices [in] -> A vertexShelf.
  !!   faces [in]    -> A faceShelf.
  !!
  !! Result:
  !!   isIt          -> .true. if the element is convex.
  !!
  function computeConvexity(self, faceIdxs, vertexIdxs, faces, vertices) result(isConvex)
    class(element), intent(inout)                :: self
    integer(shortInt), dimension(:), intent(in)  :: faceIdxs, vertexIdxs
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    logical(defBool)                             :: isConvex
    integer(shortInt)                            :: i, j, k, faceIdx, absFaceIdx, vertexIdx
    integer(shortInt), dimension(:), allocatable :: faceVertexIdxs
    real(defReal), dimension(3)                  :: normal, faceVertexCoords
    
    ! Now loop through all the faces in the element.
    do i = 1, size(faceIdxs)
      ! Retrieve the current face's vertices and signed normal vector.
      faceIdx = faceIdxs(i)
      absFaceIdx = abs(faceIdx)
      faceVertexIdxs = faces % getFaceVertexIdxs(absFaceIdx)
      normal = faces % getFaceNormal(faceIdx)
      
      ! Loop through all the vertices in the current face.
      do j = 1, size(faceVertexIdxs)
        ! Retrieve the coordinates of the current face vertex and loop through all the vertices in the element.
        faceVertexCoords = vertices % getVertexCoordinates(faceVertexIdxs(j))
        do k = 1, size(vertexIdxs)
          ! Cycle to the next vertex if the current vertex index corresponds to the index of a vertex in the current face.
          vertexIdx = vertexIdxs(k)
          if (any(faceVertexIdxs == vertexIdx)) cycle
          
          ! Assemble the test vector and check if normal .dot. testVector > ZERO. If yes, the element
          ! is concave. Continue this test to find all the problematic faces in the element.
          if (dot_product(normal, vertices % getVertexCoordinates(vertexIdx) - faceVertexCoords) > ZERO) then
            call append(self % concaveFaceIdxs, absFaceIdx, .true.)

          end if

        end do

      end do

    end do
    
    ! If reached this point the element is convex. Update isIt = .true.
    isConvex = .not. allocated(self % concaveFaceIdxs)

  end function computeConvexity

  !! Subroutine 'computeIntersectedFace'
  !!
  !! Basic description:
  !!   Computes the element's face which is intersected by a line segment.
  !!
  !! Detailed description:
  !!   See Macpherson, et al. (2009). DOI: 10.1002/cnm.1128.
  !!
  !! Arguments:
  !!   startPos [in]            -> Beginning of the line segment.
  !!   endPos [in]              -> End of the line segment.
  !!   potentialFaceIdxs [in]   -> Array of potential faces intersected by the line segment.
  !!   intersectedFaceIdx [out] -> Index of the face intersected by the line segment.
  !!   lambda [out]             -> Fraction of the line segment to be traversed before reaching the
  !!                               intersection point.
  !!   faces [in]               -> A faceShelf.
  !!
  pure subroutine computeIntersectedFace(self, r, rEnd, potentialFaceIdxs, intersectedFaceIdx, lambda, faces)
    class(element), intent(in)                  :: self
    real(defReal), dimension(3), intent(in)     :: r, rEnd
    integer(shortInt), dimension(:), intent(in) :: potentialFaceIdxs
    integer(shortInt), intent(out)              :: intersectedFaceIdx
    real(defReal), intent(out)                  :: lambda
    type(faceShelf), intent(in)                 :: faces
    integer(shortInt)                           :: i, faceIdx, absFaceIdx
    real(defReal), dimension(3)                 :: normal
    real(defReal)                               :: faceLambda
    
    ! Initialise lambda = INF.
    lambda = INF
    
    ! Loop over all potentially intersected triangles.
    do i = 1, size(potentialFaceIdxs)
      ! Retrieve the current face's signed normal vector and flip it if necessary.
      faceIdx = potentialFaceIdxs(i)
      absFaceIdx = abs(faceIdx)
      normal = faces % getFaceNormal(faceIdx)
      
      ! Compute lambda.
      faceLambda = dot_product(faces % getFaceCentroid(absFaceIdx) - r, normal) / dot_product(rEnd - r, normal)
      
      ! If triangleLambda < lambda, update lambda and intersectedTriangleIdx.
      if (faceLambda < lambda) then
        lambda = faceLambda
        intersectedFaceIdx = absFaceIdx

      end if

    end do

  end subroutine computeIntersectedFace

  !! Function 'computePotentialFaces'
  !!
  !! Basic description:
  !!   Computes a list of indices of the potentially intersected faces using the element's 
  !!   centroid and the end of a line segment.
  !!
  !! Detailed description:
  !!   See Macpherson, et al. (2009). DOI: 10.1002/cnm.1128.
  !!
  !! Arguments:
  !!   endPos [in]       -> End of the line segment.
  !!   faces [in]        -> A faceShelf.
  !!
  !! Result:
  !!   potentialFaceIdxs -> Array of indices of the potential faces intersected by the line segment.
  !!
  pure function computePotentialFaces(self, rEnd, faces) result(potentialFaceIdxs)
    class(element), intent(in)                   :: self
    real(defReal), dimension(3), intent(in)      :: rEnd
    type(faceShelf), intent(in)                  :: faces
    integer(shortInt), dimension(:), allocatable :: potentialFaceIdxs
    integer(shortInt)                            :: i, faceIdx, absFaceIdx
    real(defReal)                                :: lambda
    real(defReal), dimension(3)                  :: centroid, normal
    
    ! Retrieve element's centroid then loop over all faces in the tetrahedron.
    allocate(potentialFaceIdxs(0))
    centroid = self % centroid
    do i = 1, size(self % faceIdxs)
      ! Retrieve current face index and create its absolute index.
      faceIdx = self % faceIdxs(i)
      absFaceIdx = abs(faceIdx)
      
      ! Retrieve the signed normal vector of the current face.
      normal = faces % getFaceNormal(faceIdx)
      
      ! Retrieve the centre of the current face and compute lambda.
      lambda = dot_product(faces % getFaceCentroid(absFaceIdx) - centroid, normal) / dot_product(rEnd - centroid, normal)
      
      ! If ZERO <= lambda <= ONE, append the current face to the list of potentially intersected faces.
      if (ZERO <= lambda .and. lambda <= ONE) call append(potentialFaceIdxs, faceIdx)

    end do

  end function computePotentialFaces
  
  !! Function 'getCentroid'
  !!
  !! Basic description:
  !!   Returns the centroid of the element.
  !!
  !! Result:
  !!   centroid -> A vector pointing to the centroid of the element.
  !!
  pure function getCentroid(self) result(centroid)
    class(element), intent(in)  :: self
    real(defReal), dimension(3) :: centroid
    
    centroid = self % centroid

  end function getCentroid

  !! Function 'getEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the edges in the element.
  !!
  !! Result:
  !!   edgeIdxs -> Indices of the edges in the element.
  !!
  pure function getEdgeIdxs(self) result(edgeIdxs)
    class(element), intent(in)                          :: self
    integer(shortInt), dimension(size(self % edgeIdxs)) :: edgeIdxs

    edgeIdxs = self % edgeIdxs

  end function getEdgeIdxs
  
  !! Function 'getFaces'
  !!
  !! Basic description:
  !!   Returns the indices of the faces in the element.
  !!
  !! Result:
  !!   faceIdxs -> An array listing the indices of the faces in the element.
  !!
  pure function getFaceIdxs(self) result(faceIdxs)
    class(element), intent(in)                          :: self
    integer(shortInt), dimension(size(self % faceIdxs)) :: faceIdxs
    
    faceIdxs = self % faceIdxs

  end function getFaceIdxs
  
  !! Function 'getIdx'
  !!
  !! Basic description:
  !!   Returns the index of the element.
  !!
  !! Result:
  !!   idx -> Index of the element.
  !!
  elemental function getIdx(self) result(idx)
    class(element), intent(in) :: self
    integer(shortInt)          :: idx
    
    idx = self % idx

  end function getIdx

  elemental function getIsConvex(self) result(isConvex)
    class(element), intent(in) :: self
    logical(defBool)           :: isConvex

    isConvex = self % isConvex

  end function getIsConvex

  !!
  !!
  !!
  pure function getNotches(self) result(notches)
    class(element), intent(in)             :: self
    type(notch), dimension(:), allocatable :: notches

    if (allocated(self % notches)) then
      notches = self % notches

    else
      allocate(notches(0))

    end if

  end function getNotches

  !! Function 'getParentIdx'
  !!
  !! Basic description:
  !!   Returns the index of the parent element of the element.
  !!
  !! Result:
  !!   parentIdx -> Index of the parent element of the element.
  !!
  elemental function getParentIdx(self) result(parentIdx)
    class(element), intent(in) :: self
    integer(shortInt)          :: parentIdx

    parentIdx = self % parentIdx

  end function getParentIdx

  !!
  !!
  !!
  pure function getType(self) result(type)
    class(element), intent(in) :: self
    character(:), allocatable  :: type

    type = self % type

  end function  getType
  
  !! Function 'getVertices'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in the element.
  !!
  !! Result:
  !!   vertexIdxs -> An array listing indices of the vertices in the element.
  !!
  pure function getVertexIdxs(self) result(vertexIdxs)
    class(element), intent(in)                            :: self
    integer(shortInt), dimension(size(self % vertexIdxs)) :: vertexIdxs
    
    vertexIdxs = self % vertexIdxs

  end function getVertexIdxs
  
  !! Function 'getVolume'
  !!
  !! Basic description:
  !!   Returns the volume of the element.
  !!
  !! Result:
  !!   volume -> Volume of the element.
  !!
  elemental function getVolume(self) result(volume)
    class(element), intent(in) :: self
    real(defReal)              :: volume
    
    volume = self % volume

  end function getVolume

  !!
  !!
  !!
  pure subroutine init(self, idx, parentIdx, faceIdxs, vertexIdxs, centroid, volume, isConvex, type, edgeIdxs)
    class(element), intent(inout)                         :: self
    integer(shortInt), intent(in)                         :: idx, parentIdx
    integer(shortInt), dimension(:), intent(in)           :: faceIdxs, vertexIdxs
    real(defReal), dimension(3), intent(in)               :: centroid
    real(defReal), intent(in)                             :: volume
    logical(defBool), intent(in)                          :: isConvex
    character(*), intent(in)                              :: type
    integer(shortInt), dimension(:), intent(in), optional :: edgeIdxs

    ! Set everything.
    self % idx = idx
    self % parentIdx = parentIdx
    self % faceIdxs = faceIdxs
    self % vertexIdxs = vertexIdxs
    self % centroid = centroid
    self % volume = volume
    self % isConvex = isConvex
    self % type = type
    if (present(edgeIdxs)) self % edgeIdxs = edgeIdxs

  end subroutine init
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(element), intent(inout) :: self
    
    self % idx = 0
    self % parentIdx = 0
    self % volume = ZERO
    self % centroid = ZERO
    self % isConvex = .false.
    if (allocated(self % edgeIdxs)) deallocate(self % edgeIdxs)
    if (allocated(self % vertexIdxs)) deallocate(self % vertexIdxs)
    if (allocated(self % faceIdxs)) deallocate(self % faceIdxs)
    if (allocated(self % tetrahedronIdxs)) deallocate(self % tetrahedronIdxs)
    if (allocated(self % type)) deallocate(self % type)
    if (allocated(self % concaveFaceIdxs)) deallocate(self % concaveFaceIdxs)
    if (allocated(self % notches)) deallocate(self % notches)

  end subroutine kill
  
  !! Subroutine 'setIdx'
  !!
  !! Basic description:
  !!   Sets the index of the element.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element.
  !!
  elemental subroutine setIdx(self, idx)
    class(element), intent(inout) :: self
    integer(shortInt), intent(in) :: idx
    
    self % idx = idx

  end subroutine setIdx

  !! Subroutine 'testForInclusion'
  !!
  !! Basic description:
  !!   Tests whether a set of 3-D coordinates is inside the element.
  !!
  !! Detailed description:
  !!   First retrieves the faces making the element up. For each face, the subroutine then checks
  !!   whether the dot product between the face's normal vector and a second vector going from the 
  !!   set of 3-D coordinates to the face's centroid is positive. If it is, then the two vectors 
  !!   point in the same direction. If this test is successful for all faces then the coordinates 
  !!   are inside the element.
  !!
  !! Notes: OpenFOAM always numbers a given face's vertices such that the normal vector to this
  !!        face points from the owner element to the neighbour one. Since neighbour elements
  !!        always have greater indices than owner ones, if a given element neighbours a given face
  !!        then the negative of this face's index is added to the 'faces' component of the
  !!        'element' structure. Therefore, in the function below if a face has a negative index,
  !!        its normal vector is flipped.
  !!
  !! Arguments:
  !!   faces [in]            -> A faceShelf.
  !!   r [in]                -> A set of 3-D coordinates.
  !!   failedFace [out]      -> Index of the last face for which the inclusion test fails.
  !!   surfTolFaceIdxs [out] -> An array listing faces for which the dot product is below
  !!                            SURF_TOL, meaning that the coordinates are on the face. It is
  !!                            used in the main tracking routine to assign an element to the
  !!                            coordinates in case the coordinates are on one or more face(s).
  !!
  pure subroutine testForInclusion(self, faces, r, failedFaceIdx, surfTolFaceIdxs)
    class(element), intent(in)                                :: self
    type(faceShelf), intent(in)                               :: faces
    real(defReal), dimension(3), intent(in)                   :: r
    integer(shortInt), intent(out)                            :: failedFaceIdx
    integer(shortInt), dimension(:), allocatable, intent(out) :: surfTolFaceIdxs
    integer(shortInt)                                         :: i, faceIdx, absFaceIdx
    real(defReal)                                             :: dotProduct
    
    ! Initialise failedFaceIdx = 0 and loop over all element faces.
    failedFaceIdx = 0
    do i = 1, size(self % faceIdxs)
      ! Create an absolute face index retrieve the face's normal and centroid vectors.
      faceIdx = self % faceIdxs(i)
      absFaceIdx = abs(faceIdx)
      
      ! Make a vector going from the coordinates to the face's centroid and perform the dot
      ! product between this vector and the face's normal vector.
      dotProduct = dot_product(faces % getFaceCentroid(absFaceIdx) - r, faces % getFaceNormal(faceIdx))

      ! Check if the coordinates actually lie on the current face, and if so append zeroDotProductFaces.
      if (areEqual(dotProduct, ZERO)) then
        call append(surfTolFaceIdxs, absFaceIdx)
        cycle

      end if

      ! If dotProduct < ZERO, update failedFace and return early.
      if (dotProduct < ZERO) then
        failedFaceIdx = absFaceIdx
        return

      end if

    end do

  end subroutine testForInclusion

end module element_inter