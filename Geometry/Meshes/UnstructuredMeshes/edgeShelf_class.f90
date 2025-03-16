module edgeShelf_class
  
  use edge_class,        only : edge
  use genericProcedures, only : fatalError, numToChar
  use numPrecision
  
  implicit none
  private
  
  type, public                            :: edgeShelf
    private
    type(edge), dimension(:), allocatable :: shelf
  contains
    procedure                             :: addChildIdxToEdge
    procedure                             :: addElementIdxToEdge
    procedure                             :: addFaceIdxToEdge
    procedure                             :: allocateShelf
    procedure                             :: collapseShelf
    procedure                             :: expandShelf
    procedure                             :: getEdgeChildrenIdxs
    procedure                             :: getEdgeCutVertexIdx
    procedure                             :: getEdgeElementIdxs
    procedure                             :: getEdgeFaceIdxs
    procedure                             :: getEdgeVertexIdxs
    procedure                             :: getSize
    procedure                             :: initEdge
    procedure                             :: kill
    procedure                             :: setEdgeCutVertexIdx
  end type edgeShelf

contains
  !! Subroutine 'addChildIdxToEdge'
  !!
  !! Basic description:
  !!   Adds the index of a child edge to an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]      -> Index of the edge in the shelf.
  !!   childIdx [in] -> Index of the child edge.
  !!
  elemental subroutine addChildIdxToEdge(self, idx, childIdx)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, childIdx

    call self % shelf(idx) % addChildIdx(childIdx)

  end subroutine addChildIdxToEdge

  !! Subroutine 'addElementIdxToEdge'
  !!
  !! Basic description:
  !!   Adds the index of an element to an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the edge in the shelf.
  !!   elementIdx [in] -> Index of the element containing the edge.
  !!
  elemental subroutine addElementIdxToEdge(self, idx, elementIdx)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, elementIdx

    call self % shelf(idx) % addElementIdx(elementIdx)

  end subroutine addElementIdxToEdge

  !! Subroutine 'addFaceIdxToEdge'
  !!
  !! Basic description:
  !!   Adds the index of a face to an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the edge in the shelf.
  !!   faceIdx [in] -> Index of the face containing the edge.
  !!
  elemental subroutine addFaceIdxToEdge(self, idx, faceIdx)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, faceIdx

    call self % shelf(idx) % addFaceIdx(faceIdx)

  end subroutine addFaceIdxToEdge

  !! Subroutine 'allocateShelf'
  !!
  !! Basic description:
  !!   Allocates memory in the shelf.
  !!
  !! Arguments:
  !!   nEdges [in] -> Number of edges to be included in the shelf.
  !!
  elemental subroutine allocateShelf(self, nEdges)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: nEdges

    allocate(self % shelf(nEdges))

  end subroutine allocateShelf

  !! Subroutine 'collapseShelf'
  !!
  !! Basic description:
  !!   Reduces the size of the shelf to lastIdx.
  !!
  !! Arguments:
  !!   lastIdx [in] -> Index of the last edge in the shelf to be collapsed.
  !!
  elemental subroutine collapseShelf(self, lastIdx)
    class(edgeShelf), intent(inout)       :: self
    integer(shortInt), intent(in)         :: lastIdx
    type(edge), dimension(:), allocatable :: shelf
    
    if (allocated(self % shelf)) then
      ! Create a temporary shelf and copy all the elements up to lastIdx from the 
      ! original shelf.
      shelf = self % shelf(1:lastIdx)
      
      ! Deallocate and reallocate shelf then copy elements back.
      deallocate(self % shelf)
      allocate(self % shelf(lastIdx))
      self % shelf = shelf

    else
      allocate(self % shelf(lastIdx))

    end if

  end subroutine collapseShelf

  !! Subroutine 'expandShelf'
  !!
  !! Basic description:
  !!   Expands the shelf by a specified number of additional edges. Copies elements
  !!   already present. Allocates the shelf if it is not allocated yet.
  !!
  !! Arguments:
  !!   nAdditionalEdges [in] -> Number of additional edges to be included in the shelf.
  !!
  elemental subroutine expandShelf(self, nAdditionalEdges)
    class(edgeShelf), intent(inout)       :: self
    integer(shortInt), intent(in)         :: nAdditionalEdges
    integer(shortInt)                     :: nEdges
    type(edge), dimension(:), allocatable :: shelf

    if (allocated(self % shelf)) then
      ! If shelf is already allocated, compute the number of edges in the shelf to be expanded
      ! and copy elements already present.
      nEdges = size(self % shelf)
      shelf = self % shelf
      
      ! Deallocate shelf and reallocate to new size then copy original elements.
      deallocate(self % shelf)
      allocate(self % shelf(nEdges + nAdditionalEdges))
      self % shelf(1:nEdges) = shelf

    else
      allocate(self % shelf(nAdditionalEdges))

    end if

  end subroutine expandShelf

  !! Function 'getEdgeChildrenIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the children edges of an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the edge in the shelf.
  !!
  !! Result:
  !!   childrenIdxs -> Indices of the children edges of the edge.
  !!
  pure function getEdgeChildrenIdxs(self, idx) result(childrenIdxs)
    class(edgeShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: childrenIdxs

    childrenIdxs = self % shelf(idx) % getChildrenIdxs()

  end function getEdgeChildrenIdxs

  !! Function 'getEdgeCutVertexIdx'
  !!
  !! Basic description:
  !!   Returns the index of the vertex used to cut an edge in the shelf into children edges.
  !!
  !! Arguments:
  !!   idx [in]     -> Index of the edge in the shelf.
  !!
  !! Result:
  !!   cutVertexIdx -> Index of the cut vertex.
  !!
  elemental function getEdgeCutVertexIdx(self, idx) result(cutVertexIdx)
    class(edgeShelf), intent(in)  :: self
    integer(shortInt), intent(in) :: idx
    integer(shortInt)             :: cutVertexIdx

    cutVertexIdx = self % shelf(idx) % getCutVertexIdx()

  end function getEdgeCutVertexIdx

  !! Function 'getEdgeElementIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the elements sharing an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]    -> Index of the edge in the shelf.
  !!
  !! Result:
  !!   elementIdxs -> Indices of the elements sharing the edge.
  !!
  pure function getEdgeElementIdxs(self, idx) result(elementIdxs)
    class(edgeShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: elementIdxs

    elementIdxs = self % shelf(idx) % getElementIdxs()

  end function getEdgeElementIdxs

  !! Function 'getEdgeFaceIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the faces sharing an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the edge in the shelf.
  !!
  !! Result:
  !!   faceIdxs -> Indices of the faces sharing the edge.
  !!
  pure function getEdgeFaceIdxs(self, idx) result(faceIdxs)
    class(edgeShelf), intent(in)                 :: self
    integer(shortInt), intent(in)                :: idx
    integer(shortInt), dimension(:), allocatable :: faceIdxs

    faceIdxs = self % shelf(idx) % getFaceIdxs()

  end function getEdgeFaceIdxs

  !! Function 'getEdgeVertexIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in an edge of the shelf.
  !!
  !! Arguments:
  !!   idx [in]   -> Index of the edge in the shelf.
  !!
  !! Result:
  !!   vertexIdxs -> Indices of the vertices in the edge.
  !!
  pure function getEdgeVertexIdxs(self, idx) result(vertexIdxs)
    class(edgeShelf), intent(in)    :: self
    integer(shortInt), intent(in)   :: idx
    integer(shortInt), dimension(2) :: vertexIdxs

    vertexIdxs = self % shelf(idx) % getVertexIdxs()

  end function getEdgeVertexIdxs

  !! Function 'getSize'
  !!
  !! Basic description:
  !!   Returns the number of edges in the shelf.
  !!
  !! Result:
  !!   nEdges -> Number of edges in the shelf.
  !!
  elemental function getSize(self) result(nEdges)
    class(edgeShelf), intent(in) :: self
    integer(shortInt)            :: nEdges

    nEdges = size(self % shelf)

  end function getSize

  !! Subroutine 'initEdge'
  !!
  !! Basic description:
  !!   Initialises an edge in the shelf.
  !!
  !! Arguments:
  !!   idx [in]        -> Index of the edge in the shelf.
  !!   vertexIdxs [in] -> Indices of the vertices in the edge.
  !!
  pure subroutine initEdge(self, idx, vertexIdxs)
    class(edgeShelf), intent(inout)             :: self
    integer(shortInt), intent(in)               :: idx
    integer(shortInt), dimension(2), intent(in) :: vertexIdxs

    call self % shelf(idx) % setIdx(idx)
    call self % shelf(idx) % setVertexIdxs(vertexIdxs)

  end subroutine initEdge

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  elemental subroutine kill(self)
    class(edgeShelf), intent(inout) :: self

    if (allocated(self % shelf)) deallocate(self % shelf)

  end subroutine kill

  !! Subroutine 'setEdgeCutVertexIdx'
  !!
  !! Basic description:
  !!   Sets the index of the vertex used to cut an edge in the shelf into children edges.
  !!
  !! Arguments:
  !!   idx [in]          -> Index of the edge in the shelf.
  !!   cutVertexIdx [in] -> Index of the cut vertex.
  !!
  elemental subroutine setEdgeCutVertexIdx(self, idx, cutVertexIdx)
    class(edgeShelf), intent(inout) :: self
    integer(shortInt), intent(in)   :: idx, cutVertexIdx

    call self % shelf(idx) % setCutVertexIdx(cutVertexIdx)

  end subroutine setEdgeCutVertexIdx
  
end module edgeShelf_class