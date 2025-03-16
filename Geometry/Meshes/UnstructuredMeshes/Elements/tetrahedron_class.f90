module tetrahedron_class
  
  use edgeShelf_class,    only : edgeShelf
  use element_inter,      only : element, elementBox, kill_super => kill
  use face_inter,         only : faceBox
  use faceShelf_class,    only : faceShelf
  use genericProcedures,  only : append, computeTetrahedronCentre, computeTetrahedronVolume
  use numPrecision
  use universalVariables, only : INF, ONE, SURF_TOL, ZERO, targetNotFound
  use vertexShelf_class,  only : vertexShelf
  
  implicit none
  private
  
  !!
  !! Tetrahedron of a given OpenFOAM mesh. Results from the decomposition of polyhedral elements
  !! during the mesh importation process. Consists of a list of vertices of the triangles and 
  !! vertices in the tetrahedron. Also lists the index of the parent element from which the 
  !! tetrahedron originates.
  !!
  !! Private members:
  !!   idx          -> Index of the tetrahedron.
  !!   elementIdx   -> Index of the parent element from which the tetrahedron originates.
  !!   triangleIdxs -> Array of triangles indices making the tetrahedron up.
  !!   vertexIdxs   -> Array of vertices indices making the tetrahedron up.
  !!   volume       -> Volume of the tetrahedron.
  !!   centroid     -> Vector pointing to the centroid of the tetrahedron.
  !!
  type, public, extends(element) :: tetrahedron
    private
    integer(shortInt)            :: elementIdx = 0
  contains
    procedure                    :: computeComponents
    procedure                    :: getElement
    procedure                    :: kill
    procedure                    :: setElement
    procedure                    :: split
    procedure                    :: splitConcave
  end type tetrahedron

contains

  !!
  !!
  !!
  pure subroutine computeComponents(self, faceIdxs, vertexIdxs, faces, vertices, centroid, volume)
    class(tetrahedron), intent(inout)           :: self
    integer(shortInt), dimension(:), intent(in) :: faceIdxs, vertexIdxs
    type(faceShelf), intent(in)                 :: faces
    type(vertexShelf), intent(in)               :: vertices
    real(defReal), dimension(3), intent(out)    :: centroid
    real(defReal), intent(out)                  :: volume
    integer(shortInt)                           :: i
    real(defReal), dimension(4, 3)              :: array

    ! Create array then compute centroid and volume.
    do i = 1, 4
      array(i, :) = vertices % getVertexCoordinates(vertexIdxs(i))

    end do

    centroid = computeTetrahedronCentre(array)
    volume = computeTetrahedronVolume(array)

  end subroutine computeComponents
  
  !! Function 'getElement'
  !!
  !! Basic description:
  !!   Returns the index of the parent element from which the tetrahedron originates.
  !!
  !! Result:
  !!   elementIdx -> Index of the parent element from which the tetrahedron originates.
  !!
  elemental function getElement(self) result(elementIdx)
    class(tetrahedron), intent(in) :: self
    integer(shortInt)              :: elementIdx
    
    elementIdx = self % elementIdx

  end function getElement
  
  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an uninitialised state.
  !!
  elemental subroutine kill(self)
    class(tetrahedron), intent(inout) :: self
    
    ! Element.
    call kill_super(self)

    ! Local.
    self % elementIdx = 0

  end subroutine kill
  
  !! Subroutine 'setElement'
  !!
  !! Basic description:
  !!   Sets the index of the parent element from which the tetrahedron originates.
  !!
  !! Arguments:
  !!   elementIdx [in] -> Index of the parent element from which the tetrahedron originates.
  !!
  elemental subroutine setElement(self, elementIdx)
    class(tetrahedron), intent(inout) :: self
    integer(shortInt), intent(in)     :: elementIdx
    
    self % elementIdx = elementIdx
  end subroutine setElement

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
    class(tetrahedron), intent(inout)             :: self
    type(faceShelf), intent(inout)                :: faces, newFaces
    integer(shortInt), intent(inout)              :: lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx
    type(edgeShelf), intent(inout)                :: newEdges
    type(vertexShelf), intent(inout)              :: newVertices
    type(elementBox), dimension(:), intent(inout) :: tetrahedra
    type(faceBox), dimension(:), intent(inout)    :: triangles
    integer(shortInt)                             :: i, j, faceIdx, edgeIdx, absTriangleIdx
    integer(shortInt), dimension(4)               :: faceIdxs, triangleIdxs, vertexIdxs
    integer(shortInt), dimension(6)               :: edgeIdxs
    integer(shortInt), dimension(:), allocatable  :: faceTriangleIdxs
    
    ! Increment lastElementIdx and retrieve the indices of the faces and vertices in the tetrahedron.
    lastNewElementIdx = lastNewElementIdx + 1
    faceIdxs = self % getFaceIdxs()
    vertexIdxs = self % getVertexIdxs()

    ! Loop through all the faces of the new tetrahedron.
    do i = 1, 4
      ! Retrieve the indices of the triangles in the current face.
      faceIdx = faceIdxs(i)
      faceTriangleIdxs = faces % getFaceTriangleIdxs(abs(faceIdx))
      triangleIdxs(i) = sign(faceTriangleIdxs(1), faceIdx)
      call newFaces % addElementIdxToFace(faceTriangleIdxs(1), lastNewElementIdx)
      call newVertices % addElementIdxToVertex(vertexIdxs(i), lastNewElementIdx)

    end do

    edgeIdxs = self % getEdgeIdxs()
    do i = 1, 6
      call newEdges % addElementIdxToEdge(edgeIdxs(i), lastNewElementIdx)

    end do

    ! Initialise a new tetrahedron.
    allocate(tetrahedron :: tetrahedra(lastNewElementIdx) % item)
    call tetrahedra(lastNewElementIdx) % item % init(lastNewElementIdx, self % getIdx(), triangleIdxs, vertexIdxs, &
                                                     self % getCentroid(), self % getVolume(), self % getIsConvex(), &
                                                     'Tetrahedron', edgeIdxs)
  
  end subroutine split

  !!
  !!
  !!
  subroutine splitConcave(self, edges, faces, vertices, newEdges, convexElements, newFaces, newVertices)
    class(tetrahedron), intent(inout)                        :: self
    type(edgeShelf), intent(inout)                           :: edges, newEdges
    type(faceShelf), intent(inout)                           :: faces, newFaces
    type(vertexShelf), intent(inout)                         :: vertices, newVertices
    type(elementBox), dimension(:), allocatable, intent(out) :: convexElements
    
    ! Do nothing.

  end subroutine splitConcave
  
end module tetrahedron_class