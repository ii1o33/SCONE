module concavePolyhedron_class

    use edgeShelf_class,    only : edgeShelf
    use element_inter,      only : element, elementBox, kill_super => kill
    use face_inter,         only : faceBox
    use faceShelf_class,    only : faceShelf
    use genericProcedures,  only : append
    use numPrecision
    use universalVariables, only : INF, ONE, SURF_TOL, ZERO, targetNotFound
    use vertexShelf_class,  only : vertexShelf

    implicit none
    private

    type, public, extends(element) :: concavePolyhedron
        private
        integer(shortInt) :: idx
        contains
          procedure, non_overridable   :: decomposeConcave
          procedure, non_overridable   :: findProblemFace
          procedure, non_overridable   :: append_real
    end type concavePolyhedron

contains

  !!EDIT
  !!Module definition
  !!variable declaration of subroutine "decomposeConcave"
  !!Redundant vertices and finish extended plane def. (perhaps try new subroutine for definng extended plane?)
  !!Then continue with decomposition and assigning new verticies, edge, plane
  !!





  pure subroutine decomposeConcave(self, faceIdxs, vertexIdxs, faces, vertices, centroid, volume)
    class(concavePolyhedron), intent(inout)     :: self
    integer(shortInt), dimension(:), intent(in) :: faceIdxs, vertexIdxs
    type(faceShelf), intent(in)                 :: faces
    type(vertexShelf), intent(in)               :: vertices
    real(defReal), dimension(3), intent(out)    :: centroid
    real(defReal), intent(out)                  :: volume
    integer(shortInt)                           :: problemFace
    logical(defBool)                            :: isConvex
    real(defReal), dimension(4, 3)              :: array

    !Thorugh the convexity test, find the index of the problematic face.
    problemFace = 0
    isConvex = self % computeConvexity(faceIdxs, vertexIdxs, faces, vertices, problemFace)

    !Retrieve details of the problematic face
    absFaceIdxProb = abs(problemFace)
    normalProb = faces % getFaceNormal(absFaceIdxProb)
    centroidProb = faces % getFaceCentroid(absFaceIdxProb)
    edgeIdxsProb = faces % getFaceEdgeIdxs(absFaceIdxProb)

    !Retrieve details of the element
    elementIdx = self % elementIdx
    edgeIdxsElement = self % getElementEdgeIdxs(elementIdx)

    !Allocate arrays of coordinates for the intesection points
    allocate(intersectXcoord(0))
    allocate(intersectYcoord(0))
    allocate(intersectZcoord(0))

    !Loop through all the edges in the element
    do i = 1, size(edgeIdxsElement)
      edgeIdx = edgeIdxsElement(i)
      !If the current loop points to an edge which is in the problematic face, skip.
      if (any(edgeIdxsProb == edge)) cycle
      
      !Retrieve coordinates of the two vertices in the edge
      verticesEdge = self % edges % getEdgeElementIdxs(edgeIdx)
      coord1 = vertices % getVertexCoordinates(verticesEdge(1))
      coord2 = vertices % getVertexCoordinates(verticesEdge(2))
      
      !Calculate the instersection point. Start by calculating t
      const = normalProb(1)*centroidProb(1) + normalProb(2)*centroidProb(2) &
          normalProb(3)*centroidProb(3)
      t = -(normalProb(1)*coord1(1) + normalProb(2)*coord1(2) + normalProb(3)*coord1(3) &
            + const)/(normalProb(1)*(coord2(1)-coord1(1)) + &
            normalProb(2)*(coord2(2)-coord1(2)) + normalProb(3)*(coord2(3)-coord1(3)) +)

      !If t calculated is not valide, skip.
      if (t < ZERO .OR. t > ONE) then
        cycle
      !If t calculated is valid, append the corresponding intersection point to the list.
      else
        call append_real(intersectXcoord, (1-t)*coord1(1) + t*coord2(1))
        call append_real(intersectYcoord, (1-t)*coord1(2) + t*coord2(2))
        call append_real(intersectZcoord, (1-t)*coord1(3) + t*coord2(3))
      end if 
  
    end do 

    !The vertices on the original pland and the intersection points form the extended plant 
    !about which the element is divided into two. Before defining the extended plane,
    !remove redundant vertices from the original plane. Redundant vertices can be tested-
    !if it lies within a line connecting the other two vertices.


  end subroutine computeComponents





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
  pure function findProblemFace(self, faceIdxs, vertexIdxs, faces, vertices, problemFace) result(isConvex)
    class(element), intent(in)                   :: self
    integer(shortInt), dimension(:), intent(in)  :: faceIdxs, vertexIdxs
    type(faceShelf), intent(in)                  :: faces
    type(vertexShelf), intent(in)                :: vertices
    logical(defBool)                             :: isConvex
    integer(shortInt), intent(inout)             :: problemFace
    integer(shortInt)                            :: i, j, k, faceIdx, absFaceIdx, vertexIdx
    integer(shortInt), dimension(:), allocatable :: faceVertexIdxs
    real(defReal), dimension(3)                  :: normal, faceVertexCoords

    ! Initialise isIt = .false.
    isConvex = .false.
    
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
            problemFace = faceIdx
            return
          end if

        end do

      end do

    end do
    
    ! If reached this point the element is convex. Update isIt = .true.
    isConvex = .true.

  end function computeConvexity



    ! Subroutine to append a real number to an allocatable array
  pure subroutine append_real(arr, value)
    real, allocatable, intent(inout) :: arr(:)
    real, intent(in) :: value
    real, allocatable :: temp(:)
    integer :: oldSize

    oldSize = size(arr)

    ! Allocate a new array with an additional space
    allocate(temp(oldSize + 1))

    ! Copy old values if the array is not empty
    if (oldSize > 0) temp(1:oldSize) = arr

    ! Add the new value at the last position
    temp(oldSize + 1) = value

    ! Move the new array to arr, freeing old memory
    call move_alloc(temp, arr)
  end subroutine append_real







end module concavePolyhedron_class
