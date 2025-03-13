module elementShelf_class
  
  use edgeShelf_class,         only : edgeShelf
  use element_inter,           only : element, elementBox
  use face_inter,              only : faceBox
  use faceShelf_class,         only : faceShelf
  use numPrecision
  use polyhedron_class,        only : polyhedron
  use tetrahedron_class,       only : tetrahedron
  use vertexShelf_class,       only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Storage space for elements in a given OpenFOAM mesh.
  !!
  !! Private members:
  !!   shelf -> Array to store elements.
  !!
  type, public                                  :: elementShelf
    private
    type(elementBox), dimension(:), allocatable :: shelf
  contains
    procedure                                   :: addEdgeIdxToElement
    procedure                                   :: addElement
    procedure                                   :: addFaceIdxToElement
    procedure                                   :: addVertexIdxToElement
    procedure                                   :: allocateElement
    procedure                                   :: allocateShelf
    procedure                                   :: buildElement
    procedure                                   :: buildElementNotches
    procedure                                   :: computeFaceIntersection
    procedure                                   :: computePotentialFaceIdxs
    procedure                                   :: getElementCentroid
    procedure                                   :: getElementEdgeIdxs
    procedure                                   :: getElementFaceIdxs
    procedure                                   :: getElementIsConvex
    procedure                                   :: getElementParentIdx
    procedure                                   :: getElementType
    procedure                                   :: getElementVertexIdxs
    procedure                                   :: getElementVolume
    procedure                                   :: getSize
    procedure                                   :: initElement
    procedure                                   :: kill
    procedure                                   :: splitConcave
    procedure                                   :: splitElement
    procedure                                   :: testForInclusion
  end type elementShelf

contains

  !! Subroutine 'addEdgeIdxToElement'
  !!
  !! Basic description:
  !!   Adds the index of an edge to an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the element in the shelf.
  !!   edgeIdx [in] -> Index of the edge in the element.
  !!
  elemental subroutine addEdgeIdxToElement(self, idx, edgeIdx)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx, edgeIdx

    call self % shelf(idx) % item % addEdgeIdx(edgeIdx)

  end subroutine addEdgeIdxToElement

  !!
  !!
  !!
  elemental subroutine addElement(self, idx, item)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx
    type(elementBox), intent(in)       :: item

    self % shelf(idx) = item

  end subroutine addElement

  !! Subroutine 'addFaceIdxToElement'
  !!
  !! Basic description:
  !!   Adds the index of a face to an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the element in the shelf.
  !!   faceIdx [in] -> Index of the face in the element.
  !!
  elemental subroutine addFaceIdxToElement(self, idx, faceIdx)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx, faceIdx

    call self % shelf(idx) % item % addFaceIdx(faceIdx)

  end subroutine addFaceIdxToElement

  !! Subroutine 'addVertexIdxToElement'
  !!
  !! Basic description:
  !!   Adds the index of a vertex to an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]       -> Index of the element in the shelf.
  !!   vertexIdx [in] -> Index of the vertex in the element.
  !!
  elemental subroutine addVertexIdxToElement(self, idx, vertexIdx)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx, vertexIdx

    call self % shelf(idx) % item % addVertexIdx(vertexIdx)

  end subroutine addVertexIdxToElement

  !!
  !!
  !!
  elemental subroutine allocateElement(self, idx, type)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: idx
    character(*), intent(in)           :: type

    select case(type)
      case('Polyhedron')
        allocate(polyhedron :: self % shelf(idx) % item)
      case('Tetrahedron')
        allocate(tetrahedron :: self % shelf(idx) % item)
      case default

    end select

  end subroutine allocateElement

  !! Subroutine 'allocateShelf'
  !!
  !! Basic description:
  !!   Allocates memory in the shelf.
  !!
  !! Arguments:
  !!   nElements [in] -> Number of elements to be included in the shelf.
  !!
  elemental subroutine allocateShelf(self, nElements)
    class(elementShelf), intent(inout) :: self
    integer(shortInt), intent(in)      :: nElements

    allocate(self % shelf(nElements))

  end subroutine allocateShelf

  !!
  !!
  !!
  subroutine buildElement(self, idx, parentIdx, faceIdxs, vertexIdxs, faces, vertices, type)
    class(elementShelf), intent(inout)          :: self
    integer(shortInt), intent(in)               :: idx, parentIdx
    integer(shortInt), dimension(:), intent(in) :: faceIdxs, vertexIdxs
    type(faceShelf), intent(in)                 :: faces
    type(vertexShelf), intent(in)               :: vertices
    character(*), intent(in)                    :: type

    ! Allocate element in shelf then build components.
    call self % allocateElement(idx, type)
    call self % shelf(idx) % item % build(idx, parentIdx, faceIdxs, vertexIdxs, faces, vertices, type)

  end subroutine buildElement

  !!
  !!
  !!
  subroutine buildElementNotches(self, idx, edges, faces, vertices)
    class(elementShelf), intent(inout)          :: self
    integer(shortInt), intent(in)               :: idx
    type(edgeShelf), intent(in)                 :: edges
    type(faceShelf), intent(in)                 :: faces
    type(vertexShelf), intent(in)               :: vertices

    call self % shelf(idx) % item % buildNotches(edges, faces, vertices)

  end subroutine buildElementNotches

  !! Subroutine 'computeFaceIntersection'
  !!
  !! Basic description:
  !!   Computes the intersection of a line segment with the faces of an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]                 -> Index of the element in the shelf.
  !!   r [in]                   -> Line segment's origin coordinates.
  !!   rEnd [in]                -> Line segment's end coordinates.
  !!   potentialFaceIdxs [in]   -> Indices of the element faces potentially intersected by the line segment.
  !!   faces [in]               -> A faceShelf.
  !!   intersectedFaceIdx [out] -> Index of the intersected element face.
  !!   lambda [out]             -> Fraction of the line segment to intersection.
  !!
  pure subroutine computeFaceIntersection(self, idx, r, rEnd, potentialFaceIdxs, faces, intersectedFaceIdx, lambda)
    class(elementShelf), intent(in)             :: self
    integer(shortInt), intent(in)               :: idx
    real(defReal), dimension(3), intent(in)     :: r, rEnd
    integer(shortInt), dimension(:), intent(in) :: potentialFaceIdxs
    type(faceShelf), intent(in)                 :: faces
    integer(shortInt), intent(out)              :: intersectedFaceIdx
    real(defReal), intent(out)                  :: lambda

    call self % shelf(idx) % item % computeIntersectedFace(r, rEnd, potentialFaceIdxs, intersectedFaceIdx, lambda, faces)

  end subroutine computeFaceIntersection

  !! Function 'computePotentialFaceIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the potentially intersected faces of an element in the shelf by a line segment.
  !!
  !! Argument:
  !!   idx [in]          -> Index of the element in the shelf.
  !!   rEnd [in]         -> Line segment's end coordinates.
  !!   faces [in]        -> A faceShelf.
  !!
  !! Result:
  !!   potentialFaceIdxs -> Indices of the element's faces potentially intersected by the line segment.
  !!
  pure function computePotentialFaceIdxs(self, idx, rEnd, faces) result(potentialFaceIdxs)
    class(elementShelf), intent(in)              :: self
    integer(shortInt), intent(in)                :: idx
    real(defReal), dimension(3), intent(in)      :: rEnd
    type(faceShelf), intent(in)                  :: faces
    integer(shortInt), dimension(:), allocatable :: potentialFaceIdxs

    potentialFaceIdxs = self % shelf(idx) % item % computePotentialFaces(rEnd, faces)

  end function computePotentialFaceIdxs

  !! Function 'getElementCentroid'
  !!
  !! Basic description:
  !!   Returns the centroid of an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element in the shelf.
  !!
  !! Result:
  !!   centroid -> 3-D coordinates of the element's centroid.
  !!
  pure function getElementCentroid(self, idx) result(centroid)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    real(defReal), dimension(3)     :: centroid

    centroid = self % shelf(idx) % item % getCentroid()

  end function getElementCentroid

  !! Function 'getElementEdgeIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the edges in an element of the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element in the shelf.
  !!
  !! Result:
  !!   edgeIdxs -> Indices of the edges in the element.
  !!
  pure function getElementEdgeIdxs(self, idx) result(edgeIdxs)
    class(elementShelf), intent(in)              :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: edgeIdxs

    edgeIdxs = self % shelf(idx) % item % getEdgeIdxs()

  end function getElementEdgeIdxs

  !! Function 'getElementFaceIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the faces in an element of the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element in the shelf.
  !!
  !! Result:
  !!   faceIdxs -> Indices of the faces in the element.
  !!
  pure function getElementFaceIdxs(self, idx) result(faceIdxs)
    class(elementShelf), intent(in)              :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: faceIdxs

    faceIdxs = self % shelf(idx) % item % getFaceIdxs()

  end function getElementFaceIdxs

  !! Function 'getElementIsConvex'
  !!
  !! Basic description:
  !!   Returns .true. if an element in the shelf is convex.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element in the shelf.
  !!
  !! Result:
  !!   isConvex -> .true. if the element is convex.
  !!
  elemental function getElementIsConvex(self, idx) result(isConvex)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    logical(defBool)                :: isConvex

    isConvex = self % shelf(idx) % item % getIsConvex()

  end function getElementIsConvex

  !! Function 'getElementParentIdx'
  !!
  !! Basic description:
  !!   Returns the index of the parent element of an element of the shelf.
  !!
  !! Arguments:
  !!   idx [in]  -> Index of the element in the shelf.
  !!
  !! Result:
  !!   parentIdx -> Index of the parent element of the element
  !!
  elemental function getElementParentIdx(self, idx) result(parentIdx)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    integer(shortInt)               :: parentIdx

    parentIdx = self % shelf(idx) % item % getParentIdx()

  end function getElementParentIdx

  !!
  !!
  !!
  pure function getElementType(self, idx) result(type)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    character(:), allocatable       :: type

    type = self % shelf(idx) % item % getType()

  end function getElementType

  !! Function 'getElementVertexIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in an element of the shelf.
  !!
  !! Arguments:
  !!   idx [in]   -> Index of the element in the shelf.
  !!
  !! Result:
  !!   vertexIdxs -> Indices of the vertices in the element.
  !!
  pure function getElementVertexIdxs(self, idx) result(vertexIdxs)
    class(elementShelf), intent(in)              :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: vertexIdxs

    vertexIdxs = self % shelf(idx) % item % getVertexIdxs()

  end function getElementVertexIdxs

  !! Function 'getElementVolume'
  !!
  !! Basic description:
  !!   Returns the volume of an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element in the shelf.
  !!
  !! Result:
  !!   volume   -> Volume of the element.
  !!
  elemental function getElementVolume(self, idx) result(volume)
    class(elementShelf), intent(in) :: self
    integer(shortInt), intent(in)   :: idx
    real(defReal)                   :: volume

    volume = self % shelf(idx) % item % getVolume()

  end function getElementVolume

  !! Function 'getSize'
  !!
  !! Basic description:
  !!   Returns the number of elements in the shelf.
  !!
  !! Result:
  !!   nElements -> Number of elements in the shelf.
  !!
  elemental function getSize(self) result(nElements)
    class(elementShelf), intent(in) :: self
    integer(shortInt)               :: nElements

    nElements = size(self % shelf)

  end function getSize

  !! Subroutine 'initElement'
  !!
  !! Basic description:
  !!   Initialises an element in the shelf.
  !!
  !! Arguments:
  !!   idx [in]      -> Index of the element in the shelf.
  !!   faces [in]    -> A faceShelf.
  !!   vertices [in] -> A vertexShelf.
  !!
  subroutine initElement(self, idx, parentIdx, faceIdxs, vertexIdxs, faces, vertices, centroid, volume, isConvex, type)
    class(elementShelf), intent(inout)          :: self
    integer(shortInt), intent(in)               :: idx, parentIdx
    integer(shortInt), dimension(:), intent(in) :: faceIdxs, vertexIdxs
    type(faceShelf), intent(in)                 :: faces
    type(vertexShelf), intent(in)               :: vertices
    real(defReal), dimension(3), intent(in)     :: centroid
    real(defReal), intent(in)                   :: volume
    logical(defBool), intent(in)                :: isConvex
    character(*), intent(in)                    :: type

    call self % allocateElement(idx, type)
    call self % shelf(idx) % item % init(idx, parentIdx, faceIdxs, vertexIdxs, centroid, volume, isConvex, type)

  end subroutine initElement
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(elementShelf), intent(inout) :: self
    integer(shortInt)                  :: i
    
    if (allocated(self % shelf)) then
      do i = 1, size(self % shelf)
        deallocate(self % shelf(i) % item)

      end do
      deallocate(self % shelf)

    end if

  end subroutine kill

  !!
  !!
  !!
  subroutine splitConcave(self, idx, edges, faces, vertices, newEdges, convexElements, newFaces, newVertices)
    class(elementShelf), intent(inout)            :: self
    integer(shortInt), intent(in)                 :: idx
    type(edgeShelf), intent(inout)                :: edges, newEdges
    type(faceShelf), intent(inout)                :: faces, newFaces
    type(vertexShelf), intent(inout)              :: vertices, newVertices
    type(elementBox), dimension(:), intent(inout) :: convexElements

    call self % shelf(idx) % item % splitConcave(edges, faces, vertices, newEdges, convexElements, newFaces, newVertices)

  end subroutine splitConcave

  !! Subroutine 'splitElement'
  !!
  !! Basic description:
  !!   Splits an element in the shelf into pyramids.
  !!
  !! Arguments:
  !!   idx [in]                -> Index of the element in the shelf.
  !!   faces [in]              -> A faceShelf.
  !!   edges [inout]           -> An edgeShelf.
  !!   vertices [inout]        -> A vertexShelf.
  !!   triangles [inout]       -> A triangleShelf.
  !!   pyramids [inout]        -> A pyramidShelf.
  !!   lastEdgeIdx [inout]     -> Index of the last edge in the edgeShelf.
  !!   lastTriangleIdx [inout] -> Index of the last triangle in the triangleShelf.
  !!   lastPyramidIdx [inout]  -> Index of the last pyramid in the pyramidShelf.
  !!   lastVertexIdx [in]      -> Index of the last vertex in the vertexIdx.
  !!   
  subroutine splitElement(self, idx, faces, lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx, &
                          newEdges, newFaces, newVertices, tetrahedra, triangles)
    class(elementShelf), intent(inout)            :: self
    integer(shortInt), intent(in)                 :: idx
    type(faceShelf), intent(inout)                :: faces, newFaces
    integer(shortInt), intent(inout)              :: lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx
    type(edgeShelf), intent(inout)                :: newEdges
    type(vertexShelf), intent(inout)              :: newVertices
    type(elementBox), dimension(:), intent(inout) :: tetrahedra
    type(faceBox), dimension(:), intent(inout)    :: triangles

    call self % shelf(idx) % item % split(faces, lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx, &
                                          newEdges, newFaces, newVertices, tetrahedra, triangles)

  end subroutine splitElement

  !! Subroutine 'testForInclusion'
  !!
  !! Basic description:
  !!   Tests whether a given element in the shelf contains a point.
  !!
  !! Arguments:
  !!   idx [in]              -> Index of the element in the shelf.
  !!   r [in]                -> 3-D coordinates of the point.
  !!   faces [in]            -> A faceShelf.
  !!   failedFaceIdx [out]   -> Index of the first element's face for which the test fails.
  !!   surfTolFaceIdxs [out] -> Indices of the element's faces on which the point lies.
  !!
  pure subroutine testForInclusion(self, idx, r, faces, failedFaceIdx, surfTolFaceIdxs)
    class(elementShelf), intent(in)                           :: self
    integer(shortInt), intent(in)                             :: idx
    real(defReal), dimension(3), intent(in)                   :: r
    type(faceShelf), intent(in)                               :: faces
    integer(shortInt), intent(out)                            :: failedFaceIdx
    integer(shortInt), dimension(:), allocatable, intent(out) :: surfTolFaceIdxs

    call self % shelf(idx) % item % testForInclusion(faces, r, failedFaceIdx, surfTolFaceIdxs)

  end subroutine testForInclusion
  
end module elementShelf_class