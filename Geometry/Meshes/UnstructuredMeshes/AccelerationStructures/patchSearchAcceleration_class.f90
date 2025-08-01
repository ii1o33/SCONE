module patchSingleAcceleration_class

  use accelerationStructure_inter, only : accelerationStructure, initAccelerationStructurePayload
  use coord_class,                 only : coord
  use elementShelf_class,          only : elementShelf
  use faceShelf_class,             only : faceShelf
  use numPrecision
  use vertexShelf_class,           only : vertexShelf
  use edgeShelf_class,             only : edgeShelf
  use cartesianGridSingle_class,   only : cartesianGridSingle
  use cartesianGenericProcedures,  only : binarySearchAngle
  
  implicit none
  private

  !!
  !!
  !!
  type, public, extends(accelerationStructure) :: patchSingleAcceleration
    private
    type(cartesianGridSingle)                  :: grid
  contains
    procedure                                  :: findHostElement
    procedure                                  :: init
    procedure                                  :: kill
  end type patchSingleAcceleration

contains

  !!
  !!
  !!
  subroutine findHostElement(self, vertices, edges, faces, elements, coords)
    class(patchSingleAcceleration), intent(in)   :: self
    class(vertexShelf), intent(in)               :: vertices
    class(edgeShelf), intent(in)                 :: edges
    type(faceShelf), intent(in)                  :: faces
    type(elementShelf), intent(in)               :: elements
    type(coord), intent(inout)                   :: coords
    integer(shortInt)                            :: potentialElementIdx, i, edgeIdx, vertexIdx, pointerIdx
    integer(shortInt), dimension(3)              :: cellIdxs
    real(defReal), dimension(3)                  :: r, v_eCoord, rLocalCoord, dummyVector, phiCoord, gridBounds_min
    integer(shortInt), dimension(2)              :: currEdgeVertexIdxs
    real(defReal)                                :: thetaHat, xLocalCoord, yLocalCoord, gridSpacingReciprocal
    integer(shortInt), dimension(:), allocatable :: elementIdxsArray
    
    ! retrieve the coordinates of neutron
    ! (needs to be changed) (needs checking) (is it correct to use "getPositionToNudge" or other coordinates?)
    r = coords % getPositionToNudge()

    !!!!!
    !print*, "coord", r
    !!!!!

    ! !!!
    ! if (.NOT. self % grid % getGridIsOutsideBounds(r)) then
    !   print*, r
    ! end if
    ! !!!

    ! check if the position of neutron is inside the catesian grid bounds
    if (self % grid % getGridIsOutsideBounds(r)) return

    !print*, "here"
    ! retrieve grid properties
    gridBounds_min = self % grid % getGridBounds_min()
    gridSpacingReciprocal = self % grid % getSpacingReciprocal()

    ! find cartesian cell indices
    do i = 1, 3
      cellIdxs(i) = ceiling((r(i) - gridBounds_min(i))*(gridSpacingReciprocal))
    end do

    !!!!!
    !print*, "indices", cellIdxs
    !!!!!

    !!!
    ! print*, "here"
    ! print*, r
    ! print*, self % grid % getGridChi(cellIdxs)
    ! print*, self % grid % getGridPhi(cellIdxs)
    ! print*, self % grid % getGridPhiCapital(cellIdxs)
    !!!



    ! retrieve element index from chi mapping.
    potentialElementIdx = self % grid % getGridChi(cellIdxs)
    
    ! if element index is valid (the current cell, characterised by "cellIdxs", is fully contained within that element)
    if (potentialElementIdx > 0) then
      call coords % setElementIdx(potentialElementIdx)
      call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))
      return

    ! in case the current cell lies outside the computational domain for the unstructured mesh, return.
    ! this is tested after testing if (chi > 0) because that is the most likely case in terms of the number of the cells
    elseif (potentialElementIdx == -1) then
      return
      
    ! otherwise, the current cell intersects with either face(s) or edge(s). Start patch searching.
    else
      edgeIdx = self % grid % getGridPhiCapital(cellIdxs)

      if (edgeIdx == 0) then
        ! push the coordinates away from the current vertex (= phi)
        ! (needs to be changed) (possible improvement/acceleration for the rest of the subroutine below?)
        vertexIdx = self % grid % getGridPhi(cellIdxs)
        phiCoord = vertices % getVertexCoordinates(vertexIdx)
        dummyVector = r - phiCoord
        r = phiCoord + (self % grid % getGridWStar())/(norm2(dummyVector))*(dummyVector)

        ! !!!
        ! print*, "-----------------------------------------------------------------------"
        ! print*, "original vertex", vertexIdx
        ! print*, "original Edge", self % grid % getGridPhiCapital(cellIdxs)
        ! print*, "original Element", self % grid % getGridChi(cellIdxs)
        ! !!! 

        ! find updated cartesian cell indices
        do i = 1, 3
          cellIdxs(i) = ceiling((r(i) - gridBounds_min(i))*(gridSpacingReciprocal))
        end do

        ! get updated edge and element index
        edgeIdx = self % grid % getGridPhiCapital(cellIdxs)
        potentialElementIdx = self % grid % getGridChi(cellIdxs)

        ! if pushed coordinate has direct mapping for element index, use that
        ! (needs to be changed) (possible acceleration for this and other parts of the subroutine)
        ! (needs checking) (is pushed position has direct mapping for element idx, is it guaranteed to lie inside. OW, ">=" not "/=")
        if (potentialElementIdx /= 0) then
          call coords % setElementIdx(potentialElementIdx)
          ! (needs to be changed) (temp:there is no internal subdivision)
          !call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))
          call coords % setParentElementIdx(potentialElementIdx)
          return
        end if

        ! !!!
        ! print*, "new vertex", self % grid % getGridPhi(cellIdxs)
        ! print*, "new Edge", self % grid % getGridPhiCapital(cellIdxs)
        ! print*, "new Element", self % grid % getGridChi(cellIdxs)
        ! !!!

      end if
    
      ! calculate pseudo angle
      currEdgeVertexIdxs = edges % getEdgeVertexIdxs(edgeIdx)
      v_eCoord = vertices % getVertexCoordinates(currEdgeVertexIdxs(2))
      rLocalCoord = r - v_eCoord
      xLocalCoord = dot_product(rLocalCoord, edges % getEdgeLocalBasis1(edgeIdx))
      yLocalCoord = dot_product(rLocalCoord, edges % getEdgeLocalBasis2(edgeIdx))
      thetaHat = SIGN(1 - (xLocalCoord / (abs(xLocalCoord) + abs(yLocalCoord))), yLocalCoord)

      ! perform binary search and return index pointer
      ! (needs checking) (index order and mechanics)
      !isBoundary = edges % getEdgeIsBoundary(edgeIdx) !!! not needed anymore
      pointerIdx = binarySearchAngle(edges % getEdgeAnglesArray(edgeIdx), thetaHat)
      elementIdxsArray = edges % getEdgeElementIdxsArray(edgeIdx)
      potentialElementIdx = elementIdxsArray(pointerIdx)

      ! if the neutron turns out to lie outside the mesh domain, return 
      if (potentialElementIdx == 0) return

      call coords % setElementIdx(potentialElementIdx)
      call coords % setParentElementIdx(elements % getElementParentIdx(potentialElementIdx))
      return




    end if

    !outside of this subroutine:
    ! check if functions at cellClass is callable
    ! check spacingReciprocal
    ! getphi, getPhiCapital, getChi




  end subroutine findHostElement

  !!
  !!
  !!
  subroutine init(self, payload)
    class(patchSingleAcceleration), intent(inout)      :: self
    type(initAccelerationStructurePayload), intent(in) :: payload

    ! Simply initialise the cartesian single-layered grid
    call self % grid % init(payload)

  end subroutine init

  !!
  !!
  !!
  elemental subroutine kill(self)
    class(patchSingleAcceleration), intent(inout) :: self

    ! Local.
    !call self % tree % kill()

  end subroutine kill

end module patchSingleAcceleration_class