module edge_class
  
  use numPrecision
  use genericProcedures, only : append
  
  implicit none
  private
  
  !!
  !! Edge of a mesh linking two vertices.
  !!
  !! Private members:
  !!   idx            -> Index of the edge.
  !!   startVertexIdx -> Index of the first vertex in the edge.
  !!   endVertexIdx   -> Index of the end vertex in the edge.
  !!   edgeToFaces    -> Array that stores edge-to-faces connectivity information.
  !!   edgeToElements -> Array that stores edge-to-elements connectivity information.
  !!
  type, public                                   :: edge
    private
    integer(shortInt)                            :: idx = 0, cutVertexIdx = 0
    integer(shortInt), dimension(2)              :: vertexIdxs = 0
    integer(shortInt), dimension(:), allocatable :: childrenIdxs, faceIdxs, elementIdxs
  contains
    ! Build procedures.#
    procedure                                    :: addChildIdx
    procedure                                    :: addElementIdx
    procedure                                    :: addFaceIdx
    procedure                                    :: kill
    procedure                                    :: setCutVertexIdx
    procedure                                    :: setIdx
    procedure                                    :: setVertexIdxs
    ! Runtime procedures.
    procedure                                    :: getChildrenIdxs
    procedure                                    :: getCutVertexIdx
    procedure                                    :: getElementIdxs
    procedure                                    :: getFaceIdxs
    procedure                                    :: getIdx
    procedure                                    :: getVertexIdxs
  end type edge

contains

  !! Subroutine 'addChildIdx'
  !!
  !! Basic description:
  !!   Adds the index of a child edge to the edge.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the child edge.
  !!
  elemental subroutine addChildIdx(self, idx)
    class(edge), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    call append(self % childrenIdxs, idx)

  end subroutine addChildIdx

  !! Subroutine 'addElementIdx'
  !!
  !! Basic description:
  !!   Adds the index of an element sharing the edge.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the element.
  !!
  elemental subroutine addElementIdx(self, idx)
    class(edge), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    call append(self % elementIdxs, idx, .true.)

  end subroutine addElementIdx

  !! Subroutine 'addFaceIdx'
  !!
  !! Basic description:
  !!   Adds the index of a face sharing the edge.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the face.
  !!
  elemental subroutine addFaceIdx(self, idx)
    class(edge), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    call append(self % faceIdxs, idx)

  end subroutine addFaceIdx

  !! Function 'getChildrenIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the children edges of the edge.
  !!
  !! Result:
  !!   childrenIdxs -> Indices of the children edges of the edge.
  !!
  pure function getChildrenIdxs(self) result(childrenIdxs)
    class(edge), intent(in)                      :: self
    integer(shortInt), dimension(:), allocatable :: childrenIdxs

    childrenIdxs = self % childrenIdxs

  end function getChildrenIdxs

  !! Function 'getCutVertexIdx'
  !!
  !! Basic description:
  !!   Returns the index of the index used to cut the edge into children edges.
  !!
  !! Result:
  !!   cutVertexIdx -> Index of the cut vertex.
  !!
  elemental function getCutVertexIdx(self) result(cutVertexIdx)
    class(edge), intent(in) :: self
    integer(shortInt)       :: cutVertexIdx

    cutVertexIdx = self % cutVertexIdx

  end function getCutVertexIdx

  !! Function 'getElementIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the elements sharing the edge.
  !!
  !! Result:
  !!   elementIdxs -> Indices of the elements sharing the edge.
  !!
  pure function getElementIdxs(self) result(faceIdxs)
    class(edge), intent(in)                                :: self
    integer(shortInt), dimension(size(self % elementIdxs)) :: faceIdxs

    faceIdxs = self % elementIdxs

  end function getElementIdxs

  !! Function 'getFaceIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the faces sharing the edge.
  !!
  !! Result:
  !!   faceIdxs -> Indices of the faces sharing the edge.
  !!
  pure function getFaceIdxs(self) result(faceIdxs)
    class(edge), intent(in)                             :: self
    integer(shortInt), dimension(size(self % faceIdxs)) :: faceIdxs

    faceIdxs = self % faceIdxs

  end function getFaceIdxs

  !! Function 'getIdx'
  !!
  !! Basic description:
  !!   Returns the index of the edge.
  !!
  !! Result:
  !!   idx -> Index of the edge.
  !!
  elemental function getIdx(self) result(idx)
    class(edge), intent(in) :: self
    integer(shortInt)       :: idx

    idx = self % idx

  end function getIdx

  !! Function 'getVertexIdxs'
  !!
  !! Basic description:
  !!   Returns the indices of the vertices in the edge.
  !!
  !! Result:
  !!   vertexIdxs -> Indices of the vertices in the edge.
  !!
  pure function getVertexIdxs(self) result(vertexIdxs)
    class(edge), intent(in)         :: self
    integer(shortInt), dimension(2) :: vertexIdxs

    vertexIdxs = self % vertexIdxs

  end function getVertexIdxs

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  elemental subroutine kill(self)
    class(edge), intent(inout) :: self

    self % idx = 0
    self % vertexIdxs = 0
    self % cutVertexIdx = 0
    if (allocated(self % childrenIdxs)) deallocate(self % childrenIdxs)
    if (allocated(self % faceIdxs)) deallocate(self % faceIdxs)
    if (allocated(self % elementIdxs)) deallocate(self % elementIdxs)

  end subroutine kill

  !! Subroutine 'setCutVertexIdx'
  !!
  !! Basic description:
  !!   Sets the index of the vertex used to cut the edge into children edges.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the cut vertex.
  !!
  elemental subroutine setCutVertexIdx(self, idx)
    class(edge), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    self % cutVertexIdx = idx

  end subroutine setCutVertexIdx

  !! Subroutine 'setIdx'
  !!
  !! Basic description:
  !!   Sets the index of the edge.
  !!
  !! Arguments:
  !!   idx [in] -> Index of the edge.
  !!
  elemental subroutine setIdx(self, idx)
    class(edge), intent(inout)    :: self
    integer(shortInt), intent(in) :: idx

    self % idx = idx

  end subroutine setIdx

  !! Subroutine 'setVertexIdxs'
  !!
  !! Basic description:
  !!   Sets the indices of the vertices in the edge.
  !!
  !! Arguments:
  !!   vertexIdxs [in] -> Array containing the indices of the vertices in the edge.
  !!
  pure subroutine setVertexIdxs(self, vertexIdxs)
    class(edge), intent(inout)                  :: self
    integer(shortInt), dimension(2), intent(in) :: vertexIdxs

    self % vertexIdxs = vertexIdxs

  end subroutine setVertexIdxs

end module edge_class