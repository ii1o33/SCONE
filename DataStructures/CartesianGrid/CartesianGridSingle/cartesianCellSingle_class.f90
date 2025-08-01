module cartesianCellSingle_class
  
  use numPrecision
  use universalVariables,           only : ZERO
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  use cartesianInitProcedures

  implicit none
  private
  
  !!
  !!
  type, public                                          :: cartesianCellSingle
    private
    integer(shortInt)                                   :: phi = 0, phiCapital = 0, chi = 0
    integer(shortInt), dimension(2)                     :: faceIdxs = 0
    ! (needs to be changed) (cell centre should be a property to avoid repetative calc.)
    ! (due to limited memory, this is calculated each time needed)
    ! (try to avoid adding properties tho due to memory)

  contains
    ! Build procedures.
    procedure                                    :: cellTestEdgeIntersection
    procedure                                    :: cellTestPolyhedronInclusion
    procedure                                    :: cellTestFaceIntersection
    procedure                                    :: cellConstructMapSingleFace
    procedure                                    :: setIsOutsideMesh
    procedure                                    :: cellFinitePrecision
    ! Runtime procedures.
    procedure                                    :: getPhiCapital
    procedure                                    :: getPhi
    procedure                                    :: getChi
    procedure                                    :: getFaceIdxs

  end type cartesianCellSingle

contains

  !!
  !!
  !! (needs to be changed) (Inefficiency due to refactoring original code) (mapping matrices in cartesianCell_class)
  !! (can be moved and stored as a attribute in cartesianGrid. Then call testEdgeIntersection directly from cartesianGrid_class)
  !! (For now this is fine becase only 4% of initialisation time is increased by this - in face no increase in init time observed later.)
  !! (Tried this before and after specifying "elemental" and "pure". Keeping the original architecture (with separte class for cells))
  !! (is faster than without for initialisation, and negligible difference for in-cycle)
  pure subroutine cellTestEdgeIntersection(self, vertices, edges, edgeIdx, circumscribedBallRadius, &
                                      targetDistance, centroid, currEdgeVector, currVertexIdxs, a)
    class(cartesianCellSingle), intent(inout)           :: self
    type(topologicalObjectShelf), intent(in)            :: vertices, edges
    integer(shortInt), intent(in)                       :: edgeIdx
    real(defReal), intent(in)                           :: circumscribedBallRadius, targetDistance, a
    real(defReal), dimension(3), intent(in)             :: centroid, currEdgeVector
    integer(shortInt), dimension(2), intent(in)         :: currVertexIdxs

    call testEdgeIntersection(vertices, edges, edgeIdx, circumscribedBallRadius, targetDistance, centroid, &
                              currEdgeVector, currVertexIdxs, a, self % phi, self % phiCapital)

  end subroutine cellTestEdgeIntersection

  !!
  !!
  !!
  subroutine cellTestPolyhedronInclusion(self, faces, currElementFaceIdxs, centroid, &
                                         faceNormalSigns, elementIdx)
    class(cartesianCellSingle), intent(inout)           :: self
    type(topologicalObjectShelf), intent(in)            :: faces
    integer(shortInt), dimension(:), intent(in)         :: currElementFaceIdxs
    real(defReal), dimension(3), intent(in)             :: centroid
    real(defReal), dimension(:,:), intent(in)           :: faceNormalSigns
    integer(shortInt), intent(in)                       :: elementIdx

    call testPolyhedronInclusion(faces, currElementFaceIdxs, centroid, faceNormalSigns, elementIdx, self % chi)

  end subroutine cellTestPolyhedronInclusion

  !!
  !!
  !!
  subroutine cellTestFaceIntersection(self, vertices, edges, faces, currVertexIdxs, extraDistance, currFaceNormal, &
                                      centroid, cellSpacing, faceIdx, currFaceEdgeIdxs, targetDistance)
    class(cartesianCellSingle), intent(inout)           :: self
    type(topologicalObjectShelf), intent(in)            :: vertices, edges, faces
    integer(shortInt), dimension(:), intent(in)         :: currVertexIdxs, currFaceEdgeIdxs
    real(defReal), dimension(3), intent(in)             :: currFaceNormal, centroid
    real(defReal), intent(in)                           :: extraDistance, cellSpacing, targetDistance
    integer(shortInt), intent(in)                       :: faceIdx

    !!!!!
    ! print*, self%chi
    ! print*, self%Phi
    ! print*, self%phiCapital
    ! print*, self%faceIdxs
    !!!!!

    call testFaceIntersection(vertices, edges, faces, currVertexIdxs, extraDistance, currFaceNormal, &
                                  centroid, cellSpacing, faceIdx, currFaceEdgeIdxs, targetDistance, self % chi, &
                                  self % phi, self % phiCapital, self % faceIdxs)

    !!!!!
    !print*, self%faceIdxs
    !!!!!

  end subroutine cellTestFaceIntersection

  !!
  !!
  !!
  pure subroutine cellconstructMapSingleFace(self, faces)
    class(cartesianCellSingle), intent(inout)           :: self
    type(topologicalObjectShelf), intent(in)            :: faces
    integer(shortInt), dimension(:), allocatable        :: currFaceEdgeIdxs!, currFaceVertexIdxs
    integer(shortInt)     :: temp

    call constructMapSingleFace(faces, self % faceIdxs, self % phiCapital)

  end subroutine cellConstructMapSingleFace

  !!
  !!
  !!
  subroutine setIsOutsideMesh(self, centroid)
    class(cartesianCellSingle), intent(inout)           :: self
    real(defReal), dimension(3), intent(in)             :: centroid

    ! if a cell intersects with neither any edge nor face, and it is not contained in a single polyhedron,
    ! then, this cell lies outside the computational domain for the unstructured mesh
    if (self % chi == 0 .AND. self % phi == 0 .AND. self % phiCapital == 0) then
      ! if lies outside, set the chi value of the cell equal to -1. There are subroutines that test 
      ! if (chi != 0), but since this subroutine is called after all of those, they are unafftected.
      self % chi = -1

      !!!!!
      !print*, "testingggggggggg", centroid
      !!!!!

    end if

  end subroutine setIsOutsideMesh

  !!
  !!
  !!
  subroutine cellFinitePrecision(self, faces, elements, centroid, elementIdx)
    class(cartesianCellSingle), intent(inout)           :: self
    type(topologicalObjectShelf), intent(in)            :: faces, elements
    real(defReal), dimension(3), intent(in)             :: centroid
    integer(shortInt), intent(in)                       :: elementIdx

    ! If the current Cartesian cell does not intersect with any faces nor included in a single mesh element,
    ! Test if centroid lies inside any mesh element. If yes, finite precision error messed it up. Hence, 
    ! update chi mapping info to that mesh element. If not, this cell lies within another element or is outside mesh domain.
    if (self%chi == 0) then                   ! The biggest proportion of cells will be included in a single mesh element
      !if (self%faceIdxs(1)==0) then           ! Then we test face intersection
      if (self%phi == 0) then
        if (self%phiCapital == 0) then
          call coverFinitePrecision(faces, elements, centroid, elementIdx, self%chi)  
        end if
      end if
    end if

  end subroutine cellFinitePrecision

  !!
  !!
  !!
  elemental function getPhiCapital(self) result(phiCapital)
    class(cartesianCellSingle), intent(in)              :: self
    integer(shortInt)                                   :: phiCapital

    phiCapital = self % phiCapital

  end function getPhiCapital

  !!
  !!
  !!
  elemental function getPhi(self) result(phi)
    class(cartesianCellSingle), intent(in)              :: self
    integer(shortInt)                                   :: phi

    phi = self % phi

  end function getPhi

  !!
  !!
  !!
  elemental function getChi(self) result(chi)
    class(cartesianCellSingle), intent(in)              :: self
    integer(shortInt)                                   :: chi

    chi = self % chi

  end function getChi

  !!
  !!
  !!
  function getFaceIdxs(self) result(faceIdxs)
    class(cartesianCellSingle), intent(in)              :: self
    integer(shortInt), dimension(2)                     :: faceIdxs

    faceIdxs = self%faceIdxs

  end function



end module CartesianCellSingle_class