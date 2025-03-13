module concavePolyhedron_class

    use edgeShelf_class,    only : edgeShelf
    use element_inter,      only : element, elementBox, kill_super => kill
    use face_inter,         only : faceBox
    use faceShelf_class,    only : faceShelf
    use genericProcedures,  only : append, findCommon
    use numPrecision
    use universalVariables, only : INF, ONE, SURF_TOL, ZERO, targetNotFound
    use vertexShelf_class,  only : vertexShelf

    implicit none
    private

    type                                :: notch
        integer(shortInt)               :: edgeIdx
        integer(shortInt), dimension(2) :: faceIdxs
    end type notch

    type, public, extends(element)             :: concavePolyhedron
        private
        type(notch), dimension(:), allocatable :: notches
        contains
          procedure                            :: computeComponents
          procedure, non_overridable           :: decomposeConcave
          procedure, non_overridable           :: findProblemFace
          procedure                            :: split
    end type concavePolyhedron

contains

  pure subroutine computeComponents(self, faceIdxs, vertexIdxs, faces, vertices, centroid, volume)
    class(concavePolyhedron), intent(inout)     :: self
    integer(shortInt), dimension(:), intent(in) :: faceIdxs, vertexIdxs
    type(faceShelf), intent(in)                 :: faces
    type(vertexShelf), intent(in)               :: vertices
    real(defReal), dimension(3), intent(out)    :: centroid
    real(defReal), intent(out)                  :: volume

  end subroutine computeComponents

  !!EDIT
  !!Module definition
  !!variable declaration of subroutine "decomposeConcave"
  !!Redundant vertices and finish extended plane def. (perhaps try new subroutine for definng extended plane?)
  !!Then continue with decomposition and assigning new verticies, edge, plane
  !!
  pure subroutine decomposeConcave(self, faceIdxs, vertexIdxs, faces, vertices, centroid, volume)
    class(concavePolyhedron), intent(inout)      :: self
    integer(shortInt), dimension(:), intent(in)  :: faceIdxs, vertexIdxs
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    real(defReal), dimension(3), intent(out)     :: centroid
    real(defReal), intent(out)                   :: volume
    logical(defBool)                             :: isConvex
    real(defReal), dimension(4, 3)               :: array
    real(defReal), dimension(3)                  :: normalProb, centroidProb
    integer(shortInt), dimension(:), allocatable :: problemFace, absFaceIdxProb, edgeIdxsProb, edgeIdxsElement

    !Thorugh the convexity test, find the index of the problematic face.
    problemFace = 0
    call self % findProblemFace(faceIdxs, vertexIdxs, faces, vertices, problemFace)

    !Retrieve details of the problematic face
    absFaceIdxProb = abs(problemFace)
    !normalProb = faces % getFaceNormal(absFaceIdxProb)
    !centroidProb = faces % getFaceCentroid(absFaceIdxProb)
    !edgeIdxsProb = faces % getFaceEdgeIdxs(absFaceIdxProb)

    !Retrieve details of the element
    edgeIdxsElement = self % getEdgeIdxs()

    !The vertices on the original pland and the intersection points form the extended plant 
    !about which the element is divided into two. Before defining the extended plane,
    !remove redundant vertices from the original plane. Redundant vertices can be tested-
    !if it lies within a line connecting the other two vertices.


  end subroutine decomposeConcave





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
  pure subroutine findProblemFace(self, faceIdxs, vertexIdxs, faces, vertices, problemFace)
    class(concavePolyhedron), intent(in)                        :: self
    integer(shortInt), dimension(:), intent(in)                 :: faceIdxs, vertexIdxs
    type(faceShelf), intent(in)                                 :: faces
    type(vertexShelf), intent(in)                               :: vertices
    integer(shortInt), dimension(:), allocatable, intent(inout) :: problemFace
    integer(shortInt)                                           :: i, j, k, faceIdx, absFaceIdx, vertexIdx
    integer(shortInt), dimension(:), allocatable                :: faceVertexIdxs
    real(defReal), dimension(3)                                 :: normal, faceVertexCoords
    
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
          ! is concave and we can return early.
          if (dot_product(normal, vertices % getVertexCoordinates(vertexIdx) - faceVertexCoords) > ZERO) then
            call append(problemFace, faceIdx)
            cycle
          end if

        end do

      end do

    end do
    

  end subroutine findProblemFace





  !!
  !!
  subroutine notchConstruction(self, edges, faceIdxs, vertexIdxs, faces, vertices, problemFace)
    class(concavePolyhedron), intent(inout)                     :: self
    type(edgeShelf), intent(in)                                 :: edges
    integer(shortInt), dimension(:), intent(in)                 :: faceIdxs, vertexIdxs
    type(faceShelf), intent(in)                                 :: faces
    type(vertexShelf), intent(in)                               :: vertices
    integer(shortInt), dimension(:), intent(inout)              :: problemFace
    integer(shortInt)                                           :: i, j, k, faceIdx, absFaceIdx, vertexIdx, nNotches
    integer(shortInt), dimension(:), allocatable                :: faceVertexIdxs, edgeIdxs, edgeFaceIdxs, currentElementFaceIdxs, &
                                                                   commonFaceIdxs
    real(defReal), dimension(3)                                 :: normal, faceVertexCoords
    type(notch), dimension(:), allocatable                      :: tempNotches

    edgeIdxs = self % getEdgeIdxs()

    do i = 1, size(edgeIdxs)
        edgeFaceIdxs = edges % getEdgeFaceIdxs(edgeIdxs(i))
        if (allocated(currentElementFaceIdxs)) deallocate(currentElementFaceIdxs)
        do j = 1, size(edgeFaceIdxs)
            if (.not. any(abs(self % getFaceIdxs()) == edgeFaceIdxs(j))) cycle
            call append(currentElementFaceIdxs, edgeFaceIdxs(j))

        end do

        commonFaceIdxs = findCommon(currentElementFaceIdxs, problemFace)
        if (size(commonFaceIdxs) == 2) then
            nNotches = nNotches + 1
            tempNotches = self % notches
            if (allocated(self % notches)) deallocate(self % notches)
            allocate(self % notches(nNotches))
            self % notches(1:nNotches - 1) = tempNotches
            self % notches(nNotches) % edgeIdx = edgeIdxs(i)
            self % notches(nNotches) % faceIdxs = commonFaceIdxs

        end if
        if (allocated(tempNotches)) deallocate(tempNotches)

    end do

  end subroutine notchConstruction



  subroutine split(self, faces, lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx, &
    newEdges, newFaces, newVertices, tetrahedra, triangles)
    class(concavePolyhedron), intent(inout)       :: self
    type(faceShelf), intent(inout)                :: faces, newFaces
    integer(shortInt), intent(inout)              :: lastNewEdgeIdx, lastNewElementIdx, lastNewFaceIdx, lastNewVertexIdx
    type(edgeShelf), intent(inout)                :: newEdges
    type(vertexShelf), intent(inout)              :: newVertices
    type(elementBox), dimension(:), intent(inout) :: tetrahedra
    type(faceBox), dimension(:), intent(inout)    :: triangles

  end subroutine split

end module concavePolyhedron_class