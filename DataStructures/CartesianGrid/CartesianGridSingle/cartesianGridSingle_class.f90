module cartesianGridSingle_class
  
  use accelerationStructure_inter,  only : initAccelerationStructurePayload
  use universalVariables,           only : ZERO, INF
  use genericProcedures,            only : crossProduct, fatalError
  use numPrecision                  
  !!!!
  use cartesianCellSingle_class,    only : cartesianCellSingle
  !use cartesianInitProcedures
  !!!!
  use cartesianGenericProcedures
  use topologicalObjectShelf_class, only : topologicalObjectShelf
  
  implicit none
  private
  
  !!
  !!
  type, public                                                :: cartesianGridSingle
    private
    real(defReal)                                             :: wStar = ZERO, alpha = ZERO, l_min = ZERO, &
                                                                 spacingReciprocal = ZERO, spacing = ZERO
    real(defReal), dimension(3)                               :: gridBounds_max = ZERO, gridBounds_min = ZERO, &
                                                                 meshBounds_max = ZERO, meshBounds_min = ZERO
    integer(shortInt), dimension(3)                           :: n_xyz = 0
    !!!!
    type(cartesianCellSingle), dimension(:,:,:), allocatable  :: grid
    !integer(shortInt), dimension(:,:,:), allocatable          :: phi, phiCapital, chi
    !integer(shortInt), dimension(:,:,:,:), allocatable        :: faceIdxs
    !!!!


  contains

    ! Build procedures.
    procedure                                    :: init
    procedure                                    :: constructMapping
    procedure                                    :: sortAngles
    procedure                                    :: setGridIsOutsideMesh
    procedure                                    :: gridFinitePrecision
    ! Runtime procedures.
    procedure                                    :: getGridBounds_min
    procedure                                    :: getSpacingReciprocal
    procedure                                    :: getGridWStar
    procedure                                    :: getGridPhiCapital
    procedure                                    :: getGridPhi
    procedure                                    :: getGridChi
    procedure                                    :: getGridIsOutsideBounds
  end type cartesianGridSingle

contains

  !!
  !!
  !!
  subroutine init(self, payload)
    class(cartesianGridSingle), intent(inout)           :: self
    type(initAccelerationStructurePayload), intent(in)  :: payload
    real(defReal), dimension(6)                         :: extremalCoordinates
    integer(shortInt)                                   :: i, j !!!!!
    integer(shortInt), dimension(:), allocatable        :: currEdgeVertexIdxs, temp1, temp2 !!!!!
    real(defReal), dimension(3)                         :: currEdgeVector, extraRoom, xyz_max, xyz_min, &
                                                           centroid !!!!!
    real(defReal)                                       :: maxCosValue, tempMaxCosValue, currEdgeLength
    character(*), parameter                             :: here = 'init (cartesianGridSingle_class.f90)'
    !integer(shortInt)                                   :: temp, j, k, temp2

    ! Check payload contents
    if (.not. associated(payload % vertices)) call fatalError(here, 'Unassociated vertices shelf.')
    if (.not. associated(payload % edges)) call fatalError(here, 'Unassociated edges shelf.')
    if (.not. associated(payload % faces)) call fatalError(here, 'Unassociated faces shelf.')
    if (.not. associated(payload % elements)) call fatalError(here, 'Unassociated elements shelf.')

    !-----------------------------------------------------------------------------------------
    ! calculate constants for each face and assign them.
    ! (constant = dot(any point on the plane ⊥ the face, face normal))
    ! (needs to be changed) (set this value in other place and intent(in) not intent(inout))
    !-----------------------------------------------------------------------------------------
    do i = 1, payload % faces % getSize()
      call payload % faces % setFaceConst(i, dot_product(payload % faces % getFaceNormal(i), &
           payload % faces % getFaceCentroid(i))*(-1))
    end do

    !-----------------------------------------------------------------------------------------
    ! set l_min. Concurrently, set edgeLength and edgeUnitVector for all edges
    ! (needs to be changed) (set this value in other place and intent(in) not intent(inout))
    !-----------------------------------------------------------------------------------------
    self % l_min = INF

    do i = 1, payload % edges % getSize()
      currEdgeVertexIdxs = payload % edges % getEdgeVertexIdxs(i)
      currEdgeVector = payload % vertices % getVertexCoordinates(currEdgeVertexIdxs(2)) &
                       - payload % vertices % getVertexCoordinates(currEdgeVertexIdxs(1))
      currEdgeLength = norm2(currEdgeVector)

      call payload % edges % setEdgeUnitVector(i, currEdgeVector/currEdgeLength)
      call payload % edges % setEdgeLength(i, currEdgeLength)

      if (currEdgeLength < self % l_min) self % l_min = currEdgeLength

    end do

    ! (needs to be changed) (written for temp operation; can be optimised further; calculate it during l_min calc.?)
    ! (can probably be a separte subroutine on its own, and the return value can be assigned as a grid attribute?)
    !print*, "Average edge length", calculateAvgEdgeLength(edges) 

    !-----------------------------------------------------------------------------------------
    ! set alpha
    !-----------------------------------------------------------------------------------------
    maxCosValue = findMinFaceAngle(payload % edges, payload % faces)

    tempMaxCosValue = findMinDihedralAngle(payload % edges, payload % faces, payload % elements)
    if (maxCosValue < tempMaxCosValue) maxCosValue = tempMaxCosValue

    self % alpha = ACOS(maxCosValue)

    !-----------------------------------------------------------------------------------------
    ! set grid dimensions
    !-----------------------------------------------------------------------------------------
    ! set cartesian cell spacing (needs to be changed) (times by 0.9999 for wStar?)
    self % wStar = (self % l_min)*min(0.5d0, SIN(self % alpha))
    self % spacing = 2*(self % wStar)*SIN(self % alpha)*SIN((self % alpha)/2)&
                     /sqrt(3.0d0)/(1+SIN(self % alpha))/(1+SIN((self % alpha)/2))
    self % spacingReciprocal = 1/(self % spacing)

    ! set n_xyz
    ! set minimum and max xyz-coordinates of the cartesian grid
    extremalCoordinates = payload % vertices % getExtremalCoordinates()
    self % meshBounds_min = extremalCoordinates(1:3)
    self % meshBounds_max = extremalCoordinates(4:6)

    !(needs to be changed)(change the number of spacing for extra room for diff layers)
    do i = 1, 3
        extraRoom(i) = mod(self % meshBounds_max(i) - self % meshBounds_min(i), self % spacing)
        self % gridBounds_min(i) = self % meshBounds_min(i) - (self % spacing - extraRoom(i))/2
        self % gridBounds_max(i) = self % meshBounds_max(i) + (self % spacing - extraRoom(i))/2

        self % n_xyz(i) = NINT((self % gridBounds_max(i) - self % gridBounds_min(i))/(self % spacing))
        
    end do

    ! print cartesian grid parameters and mesh quality
    print*, "----------------------------------------------------"
    print*, "/\/\ Cartesian grid parameters and mesh quality /\/\"
    print*, "Minimum angle            : ", self % alpha
    print*, "Minimum edge length      : ", self % l_min
    print*, "No. of vertices          : ", payload % vertices % getSize()
    print*, "No. of edges             : ", payload % edges % getSize()
    print*, "No. of faces             : ", payload % faces % getSize()
    print*, "No. of elements          : ", payload % elements % getSize()
    print*, "Grid spacing             : ", self % spacing
    print*, "Grid size in x           : ", self % n_xyz(1)
    print*, "Grid size in y           : ", self % n_xyz(2)
    print*, "Grid size in z           : ", self % n_xyz(3)
    print*, "Grid lower bounds in xyz : ", self % gridBounds_min
    print*, "Grid upper bounds in xyz : ", self % gridBounds_max
    print*, "----------------------------------------------------"

    !!!!
    ! allocate grid matrix
    allocate(self % grid(self % n_xyz(1), self % n_xyz(2), self % n_xyz(3)))
    ! allocate(self % phi(self % n_xyz(1), self % n_xyz(2), self % n_xyz(3)))
    ! allocate(self % phiCapital(self % n_xyz(1), self % n_xyz(2), self % n_xyz(3)))
    ! allocate(self % chi(self % n_xyz(1), self % n_xyz(2), self % n_xyz(3)))
    ! allocate(self % faceIdxs(self % n_xyz(1), self % n_xyz(2), self % n_xyz(3), 2))

    ! ! initialise grid marix
    ! self % phi = 0
    ! self % phiCapital = 0
    ! self % chi = 0
    ! self % faceIdxs = 0
    !!!!

    !-----------------------------------------------------------------------------------------
    !initialise for patch search
    !-----------------------------------------------------------------------------------------
    call self % constructMapping(payload % vertices, payload % edges, payload % faces, payload % elements)
    call self % sortAngles(payload % edges, payload % faces)!, vertices) !!! vertices
    !!!!!
    print*, "TEST BEGINS"
    call self % gridFinitePrecision(payload % vertices, payload % faces, payload % elements) 
    print*, "TEST ENDS"
    call self % setGridIsOutsideMesh()
    !!!!!

    !!!
    ! temp2 = 0
    ! do i = 1, self % n_xyz(1)
    !   do j = 1, self % n_xyz(2)
    !     do k = 1, self % n_xyz(3)
    !       temp = self % grid(i,j,k) % getChi()
    !       if (temp == -1) then
    !         temp2 = temp2 + 1
    !       end if
    !     end do
    !   end do 
    ! end do
    ! print*, temp2
    !!!

    !!!!!
    ! print*, self % getGridChi([34,169,178])
    ! print*, self % getGridPhi([34,169,178])
    ! print*, self % getGridPhiCapital([34,169,178])
    ! call fatalError("ehere","here")
    !!!!!

    !!!!!
    ! print*, self % getGridChi([469,480,1])
    ! print*, self % getGridPhi([469,480,1])
    ! print*, self % getGridPhiCapital([469,480,1])
    ! call fatalError("ehere","here")
    !!!!!

    !!!!!
    ! print*, payload % faces % getFaceEdgeIdxs(591)
    ! print*, payload % faces % getFaceEdgeIdxs(656)
    ! print*, payload % faces % getFaceElementIdxs(591)
    ! print*, payload % faces % getFaceElementIdxs(656)
    ! !call fatalError("ehere","here")
    !!!!!

    !!!!!
    ! print*, "FACEIDX", self % grid(34,169,178) % getFaceIdxs()
    !!!!!

    !!!!!
    ! temp1 = abs(payload % elements % getElementFaceIdxs(620))
    ! print*, "£££££££££££££££££££££££££££££"
    ! do i = 1, size(temp1)
    !   print*, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    !   print*, "face index", temp1(i)
    !   temp2 = abs(payload % faces % getFaceVertexIdxs(temp1(i)))
    !   !temp2 = abs(payload % faces % getFaceVertexIdxs(656))
    !   do j = 1, size(temp2)
    !     print*, payload % vertices % getVertexCoordinates(temp2(j))
    !   end do
    !   print*, temp2
    !   print*, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    ! end do
    ! print*, "£££££££££££££££££££££££££££££"
    ! temp1 = abs(payload % elements % getElementFaceIdxs(630))
    ! print*, "£££££££££££££££££££££££££££££"
    ! do i = 1, size(temp1)
    !   print*, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    !   print*, "face index", temp1(i)
    !   temp2 = abs(payload % faces % getFaceVertexIdxs(temp1(i)))
    !   !temp2 = abs(payload % faces % getFaceVertexIdxs(656))
    !   do j = 1, size(temp2)
    !     print*, payload % vertices % getVertexCoordinates(temp2(j))
    !   end do
    !   print*, temp2
    !   print*, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
    ! end do
    ! print*, "£££££££££££££££££££££££££££££"
    ! call fatalError("ehere","here")
    !!!!!

    !!!!!
    ! centroid(1) = (self % gridBounds_min(1)) + (self % spacing) * (469-0.5)
    ! centroid(2) = (self % gridBounds_min(2)) + (self % spacing) * (480-0.5)
    ! centroid(3) = (self % gridBounds_min(3)) + (self % spacing) * (1-0.5)
    ! print*, centroid
    ! call fatalError("ehere","here")
    !!!!!

    !!!!!
    !print*, "££££££££££££££££££££££££££££££££££££££££££££££££££££££££££££"
    !!!!!

  end subroutine init

  !!
  !!
  !!
  subroutine constructMapping(self, vertices, edges, faces, elements)
    class(cartesianGridSingle), intent(inout)             :: self
    type(topologicalObjectShelf), intent(in)              :: vertices, edges, faces, elements
    integer(shortInt)                                     :: i, j, k, l
    integer(shortInt), dimension(:), allocatable          :: currVertexIdxs, currElementFaceIdxs, currFaceEdgeIdxs
    integer(shortInt), dimension(6)                       :: AABBIndices
    real(defReal)                                         :: circumscribedBallRadius, targetDistance, a, &
                                                             extraDistance, cellSpacing
    real(defReal), dimension(3)                           :: centroid, currEdgeVector, currFaceNormal
    real(defReal), dimension(:,:), allocatable            :: faceNormalSigns
    
    cellSpacing = self % spacing

    !----------------------------------------------------------------------------------------------
    ! edge interesection tests
    !----------------------------------------------------------------------------------------------
    !calculate constants for edge interesection tests
    circumscribedBallRadius = sqrt(3.0d0)*(self % spacing)/2
    targetDistance = (self % wStar) / (1 + SIN(self % alpha))

    do i = 1, edges % getSize()

        ! construct box for candidate cells
        currVertexIdxs = edges % getEdgeVertexIdxs(i)
        AABBIndices = constructAABB(vertices, currVertexIdxs, self % gridBounds_min, self % spacing)

        ! calculate edge-only-dependent properties
        currEdgeVector = (edges % getEdgeUnitvector(i))*(edges % getEdgeLength(i))
        a = dot_product(currEdgeVector, currEdgeVector)

        !Loop over all cartesian cells in the box and test if each cell intersects with the edge
        !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
        do j = AABBIndices(1), AABBIndices(4)
            do k = AABBIndices(2), AABBIndices(5)
                do l = AABBIndices(3), AABBIndices(6)

                    ! (needs to be changed) (store centroid info)
                    centroid(1) = (self % gridBounds_min(1)) + (self % spacing) * (j-0.5)
                    centroid(2) = (self % gridBounds_min(2)) + (self % spacing) * (k-0.5)
                    centroid(3) = (self % gridBounds_min(3)) + (self % spacing) * (l-0.5)

                    !!!!
                    call self % grid(j,k,l) % cellTestEdgeIntersection(vertices, edges, i, circumscribedBallRadius, &
                                                    targetDistance, centroid, currEdgeVector, currVertexIdxs, a)
                    !call testEdgeIntersection(vertices, edges, i, circumscribedBallRadius, targetDistance, centroid, &
                    !                         currEdgeVector, currVertexIdxs, a, self % phi(j,k,l), self % phiCapital(j,k,l))
                    !!!!
                    
                end do 
            end do    
        end do

    end do

    !----------------------------------------------------------------------------------------------
    ! polyhedron inclusion tests
    !----------------------------------------------------------------------------------------------
    allocate(faceNormalSigns(3, 2))
    do i = 1, elements % getSize()

        ! construct box for candidate cells
        currVertexIdxs = elements % getElementVertexIdxs(i)
        AABBIndices = constructAABB(vertices, currVertexIdxs, self % gridBounds_min, self % spacing)

        ! calculate element-only-dependent properties
        currElementFaceIdxs = elements % getElementFaceIdxs(i)

        deallocate(faceNormalSigns)
        allocate(faceNormalSigns(3, size(currElementFaceIdxs)))
        do j = 1, size(currElementFaceIdxs)
          currFaceNormal = faces % getFaceNormal(currElementFaceIdxs(j))
          do k = 1, 3
            if (currFaceNormal(k) > 0) then
              faceNormalSigns(k, j) = 1
            else
              faceNormalSigns(k, j) = -1
            end if
          end do
        end do
        faceNormalSigns = faceNormalSigns * (self % spacing)/2


        !Loop over all cartesian cells in the box and test if each cell is entirely included in the polyhedron
        !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
        do j = AABBIndices(1), AABBIndices(4)
            do k = AABBIndices(2), AABBIndices(5)
                do l = AABBIndices(3), AABBIndices(6)

                  !!!!!
                  ! if (j==469 .AND. k==480 .AND. l==1 .AND. i==630) then
                  !   print*, "began!!!!!!!!!!!"
                  ! end if
                  !!!!!

                    ! (needs to be changed) (store centroid info)
                    centroid(1) = (self % gridBounds_min(1)) + (self % spacing) * (j-0.5)
                    centroid(2) = (self % gridBounds_min(2)) + (self % spacing) * (k-0.5)
                    centroid(3) = (self % gridBounds_min(3)) + (self % spacing) * (l-0.5)

                    !!!!
                    call self % grid(j,k,l) % cellTestPolyhedronInclusion(faces, currElementFaceIdxs, centroid, &
                                                                      faceNormalSigns, i)
                    !call testPolyhedronInclusion(faces, currElementFaceIdxs, centroid, faceNormalSigns, i, &
                    !                             self % chi(j,k,l))
                    !!!!

                  !!!!!
                  ! if (j==469 .AND. k==480 .AND. l==1 .AND. i==630) then
                  !   call fatalError("end", "end")
                  ! end if
                  !!!!!

                end do 
            end do    
        end do

    end do

    !----------------------------------------------------------------------------------------------
    ! face intersection tests
    !----------------------------------------------------------------------------------------------
    targetDistance = targetDistance**2
    do i = 1, faces % getSize()

      ! construct box for candidate cells
      currVertexIdxs = faces % getFaceVertexIdxs(i)
      AABBIndices = constructAABB(vertices, currVertexIdxs, self % gridBounds_min, self % spacing)

      ! calculate face-only-dependent properties
      currFaceEdgeIdxs = faces % getFaceEdgeIdxs(i)
      currFaceNormal = faces % getFaceNormal(i)
      extraDistance = (abs(currFaceNormal(1)) + abs(currFaceNormal(2)) + abs(currFaceNormal(3))) &
                      * (self%spacing) / 2               

      !Loop over all cartesian cells in the box and test if each cell intersect with the current face
      !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
      do j = AABBIndices(1), AABBIndices(4)
          do k = AABBIndices(2), AABBIndices(5)
              do l = AABBIndices(3), AABBIndices(6)

                  ! (needs to be changed) (store centroid info)
                  centroid(1) = (self % gridBounds_min(1)) + (self % spacing) * (j-0.5)
                  centroid(2) = (self % gridBounds_min(2)) + (self % spacing) * (k-0.5)
                  centroid(3) = (self % gridBounds_min(3)) + (self % spacing) * (l-0.5)

                  !!!!!
                  ! if (j==469 .AND. k==480 .AND. l==1 .AND. i==1208) then
                  !   print*, "began!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                  ! end if
                  !!!!!

                  !!!!
                  call self % grid(j,k,l) % cellTestFaceIntersection(vertices, edges, faces, &
                                            currVertexIdxs, extraDistance, currFaceNormal, &
                                            centroid, cellSpacing, i, currFaceEdgeIdxs, targetDistance)
                  !call testFaceIntersection(vertices, edges, faces, currVertexIdxs, extraDistance, &
                  !                          currFaceNormal, centroid, cellSpacing, i, currFaceEdgeIdxs, &
                  !                          targetDistance, self % chi(j,k,l), self % phi(j,k,l), &
                  !                          self % phiCapital(j,k,l), self % faceIdxs(j,k,l,:))
                  !!!!

                  !!!!!
                  ! if (j==34 .AND. k==169 .AND. l==178 .AND. i==656) then
                  !   call fatalError("here", "here")
                  ! end if
                  !!!!!

                  !!!!!
                  ! if (j==469 .AND. k==480 .AND. l==1 .AND. i==1208) then
                  !   print*,"elementIdxs", faces%getFaceElementIdxs(1208)
                  !   call fatalError("end", "end")
                  ! end if
                  !!!!!

              end do 
          end do    
      end do

    end do

    ! for cells that intersec more than one face, mappings constructions are performed during face intersection test
    ! for those that intersect exactly one face, mappings constructions are performed here.

    !!!
    !print*, "beginning single Face case"
    !!!

    ! loop over all cartesian cells and call relevant subroutine
    do i = 1, self % n_xyz(1)
      do j = 1, self % n_xyz(2)
        do k = 1, self % n_xyz(3)
          !!!!
          call self % grid(i,j,k) % cellConstructMapSingleFace(faces)
          !call constructMapSingleFace(faces, self % faceIdxs(i,j,k,:), self % phiCapital(i,j,k))
          !!!!
        end do
      end do 
    end do

    !!!
    !print*, "ending single Face case"
    !!!

  end subroutine constructMapping

  !!
  !!
  !!
  pure subroutine sortAngles(self, edges, faces)!, vertices) !!! vertices
    class(cartesianGridSingle), intent(inout)           :: self
    type(topologicalObjectShelf), intent(in)            :: edges, faces
    integer(shortInt)                                   :: i, j, k, l, m, n, currPhiCapital, v_e, pointerIdx, &
                                                           outer2LoopSize, currEdgeVertex1Idx
    real(defReal), dimension(3)                         :: currEdgeUnitVector, localBasis1, localBasis2, currUnitVector
    integer(shortInt), dimension(:), allocatable        :: currEdgeFaceIdxs, currFaceEdgeIdxs, faceIdxsArray, &
                                                           elementIdxsArray
    integer(shortInt), dimension(2)                     :: currEdgeVertexIdxs, currVertexIdxs, face1ElementIdxs, &
                                                           face2ElementIdxs
    real(defReal)                                       :: x, y, thetaHat
    real(defReal), dimension(:), allocatable            :: anglesArray

    !!!
    !logical                                       :: flagBoundaryFace
    !!!

    !!!
    ! integer(shortInt)                             :: temp
    ! temp = 0
    ! print*, "beginning sorting angles"
    !!!

    ! loop through all cartesian cells so that only [edge index]s, where there exists at least one [cell index] s.t. 
    ! phiCapital([cell index]) = [edge index], are used.
    ! (needs to be changed) (since all edges (most likely) are assigned for phiCapital mapping anyways, just loop through all edges?)
    do i = 1, self % n_xyz(1)
      do j = 1, self % n_xyz(2)
        do k = 1, self % n_xyz(3)

          !!!!
          currPhiCapital = self % grid(i,j,k) % getPhiCapital()
          !currPhiCapital = self % phiCapital(i,j,k)
          !!!!

          ! continue only if the current cell contains valid edge index mapping of phiCapital
          if (currPhiCapital /= 0) then
            ! continue only if the current edge index (= phiCapital) has not been used for sorting angles yet
            if (.NOT. edges % isAllocatedEdgeAnglesArray(currPhiCapital)) then

              !!!
              ! temp = temp + 1
              !!!

              !---------------------------------------------------------------------------------------------------------------
              ! construct 2D local cooridnate system (localBasis1,localBasis2) on the plane whose normal is given as the current edge's unit vector
              ! and contains the second vertex of the edge.
              !---------------------------------------------------------------------------------------------------------------
              currEdgeUnitVector = edges % getEdgeUnitVector(currPhiCapital)
              
              ! construct localBasis1
              if (abs(currEdgeUnitVector(1)) <= abs(currEdgeUnitVector(2)) .AND. &
                  abs(currEdgeUnitVector(1)) <= abs(currEdgeUnitVector(3))) then
                    localBasis1 = [0.0d0, currEdgeUnitVector(3), -currEdgeUnitVector(2)]
              elseif (abs(currEdgeUnitVector(2)) <= abs(currEdgeUnitVector(3))) then 
                    localBasis1 = [-currEdgeUnitVector(3), 0.0d0, currEdgeUnitVector(1)]
              else
                    localBasis1 = [currEdgeUnitVector(2), -currEdgeUnitVector(1), 0.0d0]
              end if
              localBasis1 = localBasis1 / norm2(localBasis1)

              ! construct localBasis2 (cross product gives the normalised vector)
              localBasis2 = crossProduct(currEdgeUnitVector, localBasis1)

              ! store localBasis1 and localBasis2 to each associated edge
              call edges % setEdgeLocalBasis1(currPhiCapital, localBasis1)
              call edges % setEdgeLocalBasis2(currPhiCapital, localBasis2)

              !---------------------------------------------------------------------------------------------------------------
              ! construct (unsorted) arrays for angles and associated elementIdxs
              !---------------------------------------------------------------------------------------------------------------
              ! retrieve relevant information
              currEdgeFaceIdxs = edges % getEdgeFaceIdxs(currPhiCapital)
              currEdgeVertexIdxs = edges % getEdgeVertexIdxs(currPhiCapital)
              v_e = currEdgeVertexIdxs(2)
              currEdgeVertex1Idx = currEdgeVertexIdxs(1)

              ! initialise arrays for angle and face index
              if (allocated(anglesArray)) deallocate(anglesArray)
              if (allocated(faceIdxsArray)) deallocate(faceIdxsArray)
              if (allocated(elementIdxsArray)) deallocate(elementIdxsArray)
              allocate(anglesArray(size(currEdgeFaceIdxs)))
              allocate(faceIdxsArray(size(currEdgeFaceIdxs)))
              allocate(elementIdxsArray(size(currEdgeFaceIdxs)))
              !elementIdxsArray(:) = 0

              ! loop through all faces attached to the current edge (= currPhiCapital)
              outer1: do l = 1, size(currEdgeFaceIdxs)
                currFaceEdgeIdxs = faces % getFaceEdgeIdxs(currEdgeFaceIdxs(l))
                
                ! loop through all edges of the current face (of the current edge = currPhiCapital)
                inner1: do m = 1, size(currFaceEdgeIdxs)
                  currVertexIdxs = edges % getEdgeVertexIdxs(currFaceEdgeIdxs(m))

                  ! test if this current edge of the current face is the one that is connected to edge = currPhiCapital
                  ! and that this current edge is not currPhiCapital itself
                  if (ANY(currVertexIdxs == v_e)) then
                    if (.NOT. ANY(currVertexIdxs == currEdgeVertex1Idx)) then

                      ! retrieve and correct the orientation of the unit vector of the edge
                      if (currVertexIdxs(1) == v_e) then
                        currUnitVector = edges % getEdgeUnitVector(currFaceEdgeIdxs(m))
                      else
                        currUnitVector = edges % getEdgeUnitVector(currFaceEdgeIdxs(m))*(-1)
                      end if

                      ! calculate the 2D local coordinates
                      x = dot_product(currUnitVector, localBasis1)
                      y = dot_product(currUnitVector, localBasis2)

                      ! calculate pseudo angle
                      thetaHat = SIGN(1 - (x / (abs(x) + abs(y))), y)

                      ! add the calculated angle and face index to each corresponding arrays
                      anglesArray(l) = thetaHat
                      faceIdxsArray(l) = currEdgeFaceIdxs(l)

                      ! once the calculation is performed for current face, move on to the next face of the current edge (= currPhiCapital)
                      ! (l = l + 1)
                      exit inner1

                    end if
                  end if

                end do inner1

              end do outer1

              ! perform index sorting on anglesArray and faceIdxsArray
              call sortPairs(anglesArray, faceIdxsArray)

              ! construct a sorted array for element index
              ! pointer of index for sorted array assignment
              pointerIdx = 1

              !!!
              ! initialise flag whether any face is a boundary face
              !flagBoundaryFace = .FALSE.
              !!!

              ! loop through all neighbouring faces
              outer2LoopSize = size(faceIdxsArray)
              outer2: do l = 1, outer2LoopSize
                face1ElementIdxs = faces % getFaceElementIdxs(faceIdxsArray(l))
                face2ElementIdxs = faces % getFaceElementIdxs(faceIdxsArray(mod(l, outer2LoopSize) + 1))

                !!!
                ! checking if any face is a boundary face
                !if (size(faces % getFaceElementIdxs(faceIdxsArray(l))) == 1) flagBoundaryFace = .TRUE.
                !!!

                ! if boundary face, set the second element index = 0, which can be used when constructing elementIdxsArray. 
                if (size(faces % getFaceElementIdxs(faceIdxsArray(l))) == 1) face1ElementIdxs(2) = 0
                if (size(faces % getFaceElementIdxs(faceIdxsArray(mod(l, outer2LoopSize) + 1))) == 1) face2ElementIdxs(2) = 0

                ! find the common element indicies from 2 X arrays of size 2
                ! (needs checking) (check and throw fatal error if there are two common element idxs-for boundary edge?)
                ! (There can be only one element attached to the face for boundary faces?)
                ! (If it is a boundary edge with only two faces attached, no need to calculate this)
                ! (During face intersection, if a cell intersects one or two boundary faces (of the same element))
                ! (then we can directly assign "chi" rather than phi and phiCapital)
                ! (and this removes the issue disscussed above because only [edge idx] s.t. there exists [cell idx] s.t.)
                ! (phiCapital([cell idx]) = [edgeIdx] are used for this subroutine)
                ! (this applies for boundary edges with any number of faces attached to it)
                ! (needs to be changed) (due to the possible acceleration above)
                middle2: do m = 1, 2
                  inner2: do n = 1, 2
                    if (face1ElementIdxs(m) == face2ElementIdxs(n)) then
                      elementIdxsArray(pointerIdx) = face1ElementIdxs(m)
                      pointerIdx = pointerIdx + 1
                      exit middle2
                    end if
                  end do inner2
                end do middle2

              end do outer2

              !!!
              ! if any of the faces is a boundary face, set the last element index = 0
              ! (needs to be changed) (it might not be the last element index in some cases?)
              ! (initially set all indices = 0 so that if face1ElementIdxs(m) != face2ElementIdxs(n))
              ! (for all m and n, space outside the unstructured mesh has element index = 0)
              ! (for this move "pointerIdx = pointerIdx + 1" out of loops) (have tried but not work?)
              ! if (flagBoundaryFace) elementIdxsArray(size(elementIdxsArray)) = 0
              ! solved problem by: if (size(faces % getFaceElementIdxs(faceIdxsArray(l))) == 1) face1ElementIdxs(2) = 0
              !!!

              ! pass and set elementIdxsArray and anglesArray to each corresponding edge
              call edges % setEdgeAnglesArray(currPhiCapital, anglesArray)
              call edges % setEdgeElementIdxsArray(currPhiCapital, elementIdxsArray)

              ! set if the current edge is a boundary edge (containing boundary face with a single element index)
              ! (needs to be changed) (can be accelerated further?)
              ! if (ANY(elementIdxsArray == 0)) then
              !   call edges % setEdgeIsBoundary(currPhiCapital)
              ! end if
              !!! not needed anymore

              ! outside of this subroutine
              !! construct: self % grid(i,j,k) % getCellPhiCapital()
              !! construct: edges % setEdge2dBasis(currPhiCapital, u, v)
              !setEdgeAnglesArray / setEdgeElementIdxsArray


              !!!
              ! print*, "-------------------------------------------------------------------"
              ! do l = 1, size(currEdgeFaceIdxs)
              !   currFaceEdgeIdxs = faces % getFaceVertexIdxs(currEdgeFaceIdxs(l))
              !   print*, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
              !   print*, "Face index", currEdgeFaceIdxs(l)
              !   do m = 1, size(currFaceEdgeIdxs)
              !     print*, vertices % getVertexCoordinates(currFaceEdgeIdxs(m))
              !   end do
              !   print*, faces % getFaceElementIdxs(currEdgeFaceIdxs(l))
              ! end do
              ! print*, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
              ! print*, localBasis1
              ! print*, localBasis2
              ! print*, anglesArray
              ! print*, faceIdxsArray
              ! print*, elementIdxsArray
              ! !print*, edges % getEdgeIsBoundary(currPhiCapital)
              ! print*, "-------------------------------------------------------------------"

              ! print*, edges % getEdgeAnglesArray(currPhiCapital)
              ! print*, edges % getEdgeElementIdxsArray(currPhiCapital)
              ! print*, currEdgeUnitVector
              ! print*, localBasis1
              ! print*, localBasis2
              ! print*, currEdgeFaceIdxs
              ! print*, currEdgeVertexIdxs
              !!!







            end if
          end if




        end do
      end do 
    end do

    !!!
    !print*, "ending sorting angles"
    ! print*, temp
    ! print*, SIGN(1 - (COS(0.0) / (abs(COS(0.0)) + abs(SIN(0.0)))), SIN(0.0))
    ! print*, SIGN(1 - (COS(3.14/4) / (abs(COS(3.14/4)) + abs(SIN(3.14/4)))), SIN(3.14/4))
    ! print*, SIGN(1 - (COS(3.14/2) / (abs(COS(3.14/2)) + abs(SIN(3.14/2)))), SIN(3.14/4))
    ! print*, SIGN(1 - (COS(3*3.14/4) / (abs(COS(3*3.14/4)) + abs(SIN(3*3.14/4)))), SIN(3*3.14/4))
    ! print*, SIGN(1 - (COS(3.14) / (abs(COS(3.14)) + abs(SIN(3.14)))), SIN(3.14))
    ! print*, SIGN(1 - (COS(5*3.14/4) / (abs(COS(5*3.14/4)) + abs(SIN(5*3.14/4)))), SIN(5*3.14/4))
    ! print*, SIGN(1 - (COS(6*3.14/4) / (abs(COS(6*3.14/4)) + abs(SIN(6*3.14/4)))), SIN(6*3.14/4))
    ! print*, SIGN(1 - (COS(7*3.14/4) / (abs(COS(7*3.14/4)) + abs(SIN(7*3.14/4)))), SIN(7*3.14/4))
    !!!



  end subroutine sortAngles

  !!
  !!
  !! 
  subroutine setGridIsOutsideMesh(self)
    class(cartesianGridSingle), intent(inout)              :: self
    integer(shortInt)                                      :: i, j, k
    real(defReal), dimension(3)                            :: centroid

    ! loop over all cartesian cells and call relevant subroutine
    do i = 1, self % n_xyz(1)
      do j = 1, self % n_xyz(2)
        do k = 1, self % n_xyz(3)
          
          !!!!!
          centroid(1) = (self % gridBounds_min(1)) + (self % spacing) * (i-0.5)
          centroid(2) = (self % gridBounds_min(2)) + (self % spacing) * (j-0.5)
          centroid(3) = (self % gridBounds_min(3)) + (self % spacing) * (k-0.5)
          !!!!!

          !!!!
          call self % grid(i,j,k) % setIsOutsideMesh(centroid)

          ! ! if a cell intersects with neither any edge nor face, and it is not contained in a single polyhedron,
          ! ! then, this cell lies outside the computational domain for the unstructured mesh
          ! if (self % chi(i,j,k) == 0 .AND. self % phi(i,j,k) == 0 .AND. self % phiCapital(i,j,k) == 0) then
          !   ! if lies outside, set the chi value of the cell equal to -1. There are subroutines that test 
          !   ! if (chi != 0), but since this subroutine is called after all of those, they are unafftected.
          !   self % chi(i,j,k) = -1
          ! end if

          !!!!

        end do
      end do 
    end do

  end subroutine setGridIsOutsideMesh

  !!
  !!
  !!
  subroutine gridFinitePrecision(self, vertices, faces, elements)
    class(cartesianGridSingle), intent(inout)             :: self
    type(topologicalObjectShelf), intent(in)              :: vertices, faces, elements
    integer(shortInt)                                     :: i, j, k, l
    integer(shortInt), dimension(:), allocatable          :: currVertexIdxs
    integer(shortInt), dimension(6)                       :: AABBIndices
    real(defReal), dimension(3)                           :: centroid


    do i = 1, elements % getSize()

        ! construct box for candidate cells
        currVertexIdxs = elements % getElementVertexIdxs(i)
        AABBIndices = constructAABB(vertices, currVertexIdxs, self % gridBounds_min, self % spacing)

        ! calculate element-only-dependent properties (It is extremly rare that the code has to test this finitePrecision.
        ! Hence, we do not pre-calculate these unlike testPolyhedronInclusion.)

        !Loop over all cartesian cells in the box and test if each cell is entirely included in the polyhedron
        !(needs to be changed) (k and l can be a function of j e.g. k = datum + slope*j so that box is narrowed down)
        do j = AABBIndices(1), AABBIndices(4)
            do k = AABBIndices(2), AABBIndices(5)
                do l = AABBIndices(3), AABBIndices(6)

                    ! (needs to be changed) (store centroid info)
                    centroid(1) = (self % gridBounds_min(1)) + (self % spacing) * (j-0.5)
                    centroid(2) = (self % gridBounds_min(2)) + (self % spacing) * (k-0.5)
                    centroid(3) = (self % gridBounds_min(3)) + (self % spacing) * (l-0.5)

                    call self % grid(j,k,l) % cellFinitePrecision(faces, elements, centroid, i)

                end do 
            end do    
        end do

    end do

  end subroutine gridFinitePrecision

  !!
  !!
  !!
  pure function getGridBounds_min(self) result(gridBounds_min)
    class(cartesianGridSingle), intent(in)              :: self
    real(defReal), dimension(3)                         :: gridBounds_min

    gridBounds_min = self % gridBounds_min

  end function getGridBounds_min

  !!
  !!
  !!
  elemental function getSpacingReciprocal(self) result(spacingReciprocal)
    class(cartesianGridSingle), intent(in)              :: self
    real(defReal)                                       :: spacingReciprocal

    spacingReciprocal = self % spacingReciprocal

  end function getSpacingReciprocal

  !!
  !!
  !!
  elemental function getGridWStar(self) result(wStar)
    class(cartesianGridSingle), intent(in)              :: self
    real(defReal)                                       :: wStar

    wStar = self % wStar

  end function getGridWStar

  !!
  !!
  !! 
  pure function getGridPhiCapital(self, cellIdxs) result(phiCapital)
    class(cartesianGridSingle), intent(in)              :: self
    integer(shortInt), dimension(3), intent(in)         :: cellIdxs
    integer(shortInt)                                   :: phiCapital

    !!!!
    phiCapital = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getPhiCapital()
    !phiCapital = self % phiCapital(cellIdxs(1), cellIdxs(2), cellIdxs(3))
    !!!!

  end function getGridPhiCapital

  !!
  !!
  !! 
  pure function getGridPhi(self, cellIdxs) result(phi)
    class(cartesianGridSingle), intent(in)              :: self
    integer(shortInt), dimension(3), intent(in)         :: cellIdxs
    integer(shortInt)                                   :: phi

    !!!!
    phi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getPhi()
    !phi = self % phi(cellIdxs(1), cellIdxs(2), cellIdxs(3))
    !!!!

  end function getGridPhi

  !!
  !!
  !! 
  pure function getGridChi(self, cellIdxs) result(chi)
    class(cartesianGridSingle), intent(in)              :: self
    integer(shortInt), dimension(3), intent(in)         :: cellIdxs
    integer(shortInt)                                   :: chi

    !!!!
    chi = self % grid(cellIdxs(1), cellIdxs(2), cellIdxs(3)) % getChi()
    !chi = self % chi(cellIdxs(1), cellIdxs(2), cellIdxs(3))
    !!!!

  end function getGridChi

  !!
  !!
  !! (needs to be changed) (possible acceleration?)
  function getGridIsOutsideBounds(self, r) result(isOutside)
    class(cartesianGridSingle), intent(in)              :: self
    real(defReal), dimension(3), intent(in)             :: r
    logical                                             :: isOutside
    integer(shortInt)                                   :: i

    isOutside = .FALSE.

    do i = 1, 3
      !!!!!
      ! if (r(3)  + 0.25098474675295712 < 0.001) then
      !   print*,"///////////////////////////////////////"
      !   print*, r
      !   print*, r(3) < self % gridBounds_min(3)
      !   print*,  -0.25098474675295712 < -0.25031493962876034
      !   print*,"///////////////////////////////////////"
      ! end if
      !!!!!
      if (r(i) > self % meshBounds_max(i) .OR. r(i) < self % meshBounds_min(i)) then
        isOutside = .TRUE.
        return
      end if
    end do

    !!!!!
    ! print*, "boundMax", self%gridBounds_max
    ! print*, "boundMin", self%gridBounds_min
    !!!!!

  end function getGridIsOutsideBounds

end module CartesianGridSingle_class