module OpenFOAMMesh_class

  use coord_class,            only : coord
  use edgeShelf_class,        only : edgeShelf
  use elementShelf_class,     only : elementShelf
  use cellZoneShelf_class,    only : cellZoneShelf
  use faceShelf_class,        only : faceShelf
  use genericProcedures,      only : append, fatalError, numToChar, openToRead, quickSort
  use numPrecision
  use unstructuredMesh_inter, only : unstructuredMesh, &
                                     distanceToBoundaryFace_super => distanceToBoundaryFace, &
                                     distanceToNextFace_super => distanceToNextFace, &
                                     findElementAndParentIdxs_super => findElementAndParentIdxs, &
                                     kill_super => kill
  use vertexShelf_class,      only : vertexShelf

  implicit none
  private

  type, public, extends(unstructuredMesh) :: OpenFOAMMesh
    private
  contains
    ! Superclass procedures.
    procedure                             :: distanceToBoundaryFace
    procedure                             :: distanceToNextFace
    procedure                             :: findElementAndParentIdxs
    ! Local procedures.
    procedure                             :: buildEdges
    procedure                             :: checkFiles
    procedure                             :: getMeshInfo
    procedure                             :: importMesh
    procedure                             :: initElementShelf
    procedure                             :: initElementZoneShelf
    procedure                             :: initFaceShelf
    procedure                             :: initVertexShelf
    procedure                             :: kill
  end type OpenFOAMMesh

  type :: elementInfo
    integer(shortInt), dimension(:), allocatable :: faceIdxs, vertexIdxs
  end type

contains

  !! Subroutine 'buildEdges'
  !!
  !! Basic description:
  !!   Constructs edges from mesh elements, faces and vertices.
  !!
  subroutine buildEdges(self, edges, elements, faces, vertices)
    class(OpenFOAMMesh), intent(inout)           :: self
    type(edgeShelf), intent(inout)               :: edges
    type(elementShelf), intent(inout)            :: elements
    type(faceShelf), intent(inout)               :: faces
    type(vertexShelf), intent(inout)             :: vertices
    integer(shortInt)                            :: i, j, k, lastIdx, nVertices, vertexIdx, nextVertexIdx, &
                                                    elementIdx, edgeIdx
    integer(shortInt), dimension(:), allocatable :: vertexIdxs, elementIdxs
    integer(shortInt), dimension(2)              :: edgeVertexIdxs
    logical(defBool)                             :: createNew

    ! Do a first pass over all elements and allocate memory. This will overshoot the actual number of
    ! edges to be created.
    lastIdx = 0
    do i = 1, self % nElements
      lastIdx = lastIdx + size(elements % getElementFaceIdxs(i)) + size(elements % getElementVertexIdxs(i)) - 2

    end do
    call edges % allocateShelf(lastIdx)

    ! Reset lastIdx = 0 then loop through all the faces in the mesh and construct edges.
    lastIdx = 0
    do i = 1, self % nFaces
      ! Retrieve the vertices in the current face then loop through all vertices in the current face.
      vertexIdxs = faces % getFaceVertexIdxs(i)
      nVertices = size(vertexIdxs)
      do j = 1, nVertices
        ! Retrieve the indices of the two vertices in the edge and sort them.
        vertexIdx = vertexIdxs(j)
        nextVertexIdx = vertexIdxs(mod(j, nVertices) + 1)
        edgeVertexIdxs = [vertexIdx, nextVertexIdx]
        call quickSort(edgeVertexIdxs)

        ! Check if there already exists an edge containing the two vertices.
        edgeIdx = vertices % findCommonEdgeIdx(vertexIdx, nextVertexIdx)
        createNew = edgeIdx == 0

        if (createNew) then
          ! If we need to create a new edge, increment lastIdx, set the indices of the vertices in 
          ! the new edge and add the new edge to the pair of vertices.
          lastIdx = lastIdx + 1
          call edges % initEdge(lastIdx, edgeVertexIdxs)
          call vertices % addEdgeIdxToVertex(vertexIdx, lastIdx)
          call vertices % addEdgeIdxToVertex(nextVertexIdx, lastIdx)
          
          ! Set edgeIdx = lastIdx.
          edgeIdx = lastIdx

        end if

        ! Add the index of the current face to the edge and add the edge to the current face.
        call edges % addFaceIdxToEdge(edgeIdx, i)
        call faces % addEdgeIdxToFace(i, edgeIdx)
        
        ! Retrieve all elements sharing the current face and update connectivity information.
        elementIdxs = faces % getFaceElementIdxs(i)
        do k = 1, size(elementIdxs)
          elementIdx = elementIdxs(k)
          call edges % addElementIdxToEdge(edgeIdx, elementIdx)
          call elements % addEdgeIdxToElement(elementIdx, edgeIdx)

        end do

      end do

    end do

    ! Now collapse the edgeShelf to the number of edges actually created and set number of edges
    ! in the mesh.
    call edges % collapseShelf(lastIdx)
    self % nEdges = lastIdx

  end subroutine buildEdges

  !! Subroutine 'checkFiles'
  !!
  !! Basic description:
  !!   The subroutine 'checkFiles' checks the existence of vital mesh data files. These files
  !!   are the 'points', 'faces' and 'neighbour' files. If any of them are missing, the subroutine
  !!   calls the procedure 'fatalError'.
  !!
  !! Arguments:
  !!   folderPath [in]     -> Path of the folder containing the files for the mesh.
  !!   nInternalFaces [in] -> Number of internal faces in the mesh. If this number is
  !!                          zero then the 'neighbour' file is not checked.
  !!
  !! Errors:
  !!   - fatalError if the 'points' file is missing.
  !!   - fatalError if the 'faces' file is missing.
  !!   - fatalError if nInternalFaces > 0 and the 'neighbour' file is missing.
  !!
  !! Notes:
  !!   The existence of the 'owner' file is checked in 'getMeshInfo'.
  !!
  subroutine checkFiles(self, folderPath, nInternalFaces)
    class(OpenFOAMMesh), intent(inout) :: self
    character(*), intent(in)           :: folderPath
    integer(shortInt), intent(in)      :: nInternalFaces
    logical(defBool)                   :: pointsFile, facesFile, neighbourFile, cellZonesFile
    character(*), parameter            :: Here = 'checkFiles (OpenFOAMMesh_class.f90)'

    ! Check the existence of the 'points' and 'faces' files and report errors if they are not found.
    inquire(file = folderPath//'points', exist = pointsFile)
    if (.not. pointsFile) call fatalError(Here, "Missing 'points' file for OpenFOAM mesh with Id: "//numToChar(self % getId())//'.')
    
    inquire(file = folderPath//'faces', exist = facesFile)
    if (.not. facesFile) call fatalError(Here, "Missing 'faces' file for OpenFOAM mesh with Id: "//numToChar(self % getId())//'.')
    
    ! If nInternal = 0, return early.
    if (nInternalFaces == 0) return

    ! If reached here check that the 'neighbour' file exists and report error if not.
    inquire(file = folderPath//'neighbour', exist = neighbourFile)
    if (.not. neighbourFile) call fatalError(Here, &
    "Missing 'neighbour' file for OpenFOAM mesh with Id: "//numToChar(self % getId())//'.')

    ! Check that the 'cellZones' file exists and register information into mesh.
    inquire(file = folderPath//'cellZones', exist = cellZonesFile)
    call self % setElementZonesFile(cellZonesFile)

  end subroutine checkFiles

  !! Subroutine 'distanceToBoundaryFace'
  !!
  !! Basic description:
  !!   Returns the distance to the mesh boundary face intersected by a particle's path. Also returns the index
  !!   of the parent element containing the intersected boundary face.
  !!
  !! See unstructuredMesh_inter for details.
  !!
  elemental subroutine distanceToBoundaryFace(self, d, coords, parentIdx)
    class(OpenFOAMMesh), intent(in) :: self
    real(defReal), intent(out)      :: d
    type(coord), intent(inout)      :: coords
    integer(shortInt), intent(out)  :: parentIdx

    call distanceToBoundaryFace_super(self, d, coords, parentIdx)

  end subroutine distanceToBoundaryFace

  !! Subroutine 'distanceToNextFace'
  !!
  !! Basic description:
  !!   Returns the distance to the next face intersected by the particle's path.
  !!
  !! See unstructuredMesh_inter for details.
  !!
  elemental subroutine distanceToNextFace(self, d, coords)
    class(OpenFOAMMesh), intent(in) :: self
    real(defReal), intent(out)      :: d
    type(coord), intent(inout)      :: coords

    call distanceToNextFace_super(self, d, coords)

  end subroutine distanceToNextFace

  !! Subroutine 'findElementAndParentIdxs'
  !!
  !! Basic description:
  !!   Returns the index of the mesh element occupied by a particle. Also returns the index of the parent mesh
  !!   element containing the occupied element.
  !!
  !! See unstructuredMesh_inter for details.
  !!
  pure subroutine findElementAndParentIdxs(self, r, u, elementIdx, parentIdx)
    class(OpenFOAMMesh), intent(in)         :: self
    real(defReal), dimension(3), intent(in) :: r, u
    integer(shortInt), intent(out)          :: elementIdx, parentIdx

    call findelementandparentidxs_super(self, r, u, elementIdx, parentIdx)

  end subroutine findElementAndParentIdxs

  !! Subroutine 'getMeshInfo'
  !!
  !! Basic description:
  !!   Retrieves some preliminary information about the mesh (number of vertices, faces, etc.) to be 
  !!   imported using the 'owner' data file.
  !!
  !! Detailed description:
  !!   Opens the 'owner' file and reads it until a line with the keyword 'note' is encountered. From
  !!   this line of text, the number of vertices (nVertices), the number of faces (nFaces), the 
  !!   number of elements (nElements) and the number of internal faces (nInternalFaces) are read.
  !!
  !! Arguments:
  !!   folderPath [in] -> Path of the folder containing the files of the mesh.
  !!
  subroutine getMeshInfo(self, folderPath)
    class(OpenFOAMMesh), intent(inout) :: self
    character(*), intent(in)           :: folderPath
    logical(defBool)                   :: ownerFile
    integer(shortInt)                  :: colIdx, nIdx, endQuoteIdx
    integer(shortInt), parameter       :: unit = 10
    character(100)                     :: string, partialString
    character(:), allocatable          :: ownerPath
    character(*), parameter            :: Here = 'getMeshInfo (OpenFOAMMesh_class.f90)'

    ! Initialise ownerPath and perform a first check for the existence of the 'owner' file. Call fatalError if not found.
    ownerPath = folderPath//'owner'
    inquire(file = ownerPath, exist = ownerFile)
    if (.not. ownerFile) call fatalError(Here, "Missing 'owner' file for OpenFOAM mesh with Id: "//numToChar(self % getId())//'.')
    
    ! Open the 'owner' file and read it. Skip lines until the line with the keyword 'note' is encountered.
    call openToRead(unit, ownerPath)
    read(unit, "(a)") string
    do while (index(string, 'note') == 0)
      read(unit, "(a)") string

    end do
    
    ! Locate the indices corresponding to the first ':' and ' ' symbols, and retrieve the variable
    ! 'nVertices', which is simply the piece of the text line located between the two aforementioned
    ! symbols.
    colIdx = index(string, ':')
    nIdx = index(string(colIdx:len_trim(string)), 'n')
    partialString = trim(adjustl(string(colIdx + 1:colIdx + nIdx)))
    read(partialString, *) self % nVertices
    
    ! Repeat for nElements.
    string(colIdx:colIdx) = ' '
    colIdx = index(string, ':')
    nIdx = index(string(colIdx:len_trim(string)), 'n')
    partialString = trim(adjustl(string(colIdx + 1:colIdx + nIdx)))
    read(partialString, *) self % nElements
    
    ! Repeat for nFaces.
    string(colIdx:colIdx) = ' '
    colIdx = index(string, ':')
    nIdx = index(string(colIdx:len_trim(string)), 'n')
    partialString = trim(adjustl(string(colIdx + 1:colIdx + nIdx)))
    read(partialString, *) self % nFaces
    
    ! Repeat for nInternalFaces. Close 'owner' file once done.
    string(colIdx:colIdx) = ' '
    colIdx = index(string, ':')
    endQuoteIdx = index(string(colIdx:len_trim(string)), '"')
    partialString = trim(adjustl(string(colIdx + 1:colIdx + endQuoteIdx - 2)))
    read(partialString, *) self % nInternalFaces
    close(unit)
    
    ! Now check existence of remaining mesh files.
    call self % checkFiles(folderPath, self % nInternalFaces)

  end subroutine getMeshInfo

  !! Subroutine 'importMesh'
  !!
  !! Basic description:
  !!   Imports an OpenFOAM mesh from the path of the folder containing the mesh files and sets the 
  !!   mesh Id from the dictionary.
  !!
  !! Detailed description:
  !!   'init' first sets the mesh Id from the supplied dictionary. It then checks the existence and 
  !!   consistency of files located in the appropriate mesh folder and allocates memory to the 
  !!   'vertices', 'faces' and 'elements' structures of the 'mesh' structure. It then imports data 
  !!   from the 'points' file into the 'vertices' structures and builds a kd-tree for the mesh from 
  !!   the various vertices' coordinates. The subroutine then proceeds to import data from the 
  !!   'faces', 'owner' and 'neighbour' files and stores it into the appropriate 'faces' and 
  !!   'elements' structure. 'init' also computes the area, centroid and normal vector for each 
  !!   face, as well as the centroid and volume of each element. From the 'faces', 'owner' and 
  !!   'neighbour' files, mesh connectivity information is also assigned to the various structures. 
  !!   The subroutine then imports data from the 'cellZones' file (if it exists) into the 
  !!   'cellZones' structures.
  !!
  !! Arguments:
  !!   folderPath [in] -> Path of the folder containing the files for the mesh geometry.
  !!   dict [in]       -> Input dictionary.
  !!
  !! Errors:
  !!   - fatalError if the mesh contains concave elements.
  !!
  subroutine importMesh(self, folderPath, centroids, edges, elements, elementZones, faces, vertices, &
    concaveElementIdxs)
    class(OpenFOAMMesh), intent(inout) :: self
    character(*), intent(in)           :: folderPath
    type(edgeShelf), intent(out)       :: edges
    type(elementShelf), intent(out)    :: elements
    type(cellZoneShelf), intent(out)   :: elementZones
    type(faceShelf), intent(out)       :: faces
    type(vertexShelf), intent(out)     :: centroids, vertices
    integer(shortInt), dimension(:), allocatable, intent(out) :: concaveElementIdxs
    character(*), parameter            :: Here = 'importMesh (OpenFOAMMesh_class.f90)'
    
    ! Retrieve preliminary information about the mesh and allocate memory.
    call self % getMeshInfo(folderPath)
    call centroids % allocateShelf(self % nElements)
    call elements % allocateShelf(self % nElements)
    call faces % allocateShelf(self % nFaces)
    call vertices % allocateShelf(self % nVertices)
    
    ! Import vertices and set mesh bounding box.
    call self % initVertexShelf(vertices, folderPath)
    call self % setBoundingBox(vertices % getExtremalCoordinates())

    ! Import faces and edges.
    call self % initFaceShelf(faces, vertices, folderPath)
    
    ! Import elements.
    call self % initElementShelf(centroids, elements, faces, vertices, folderPath, concaveElementIdxs)

    ! Build edges from elements, faces and vertices.
    call self % buildEdges(edges, elements, faces, vertices)

    ! Print original mesh composition.
    !    call self % printComposition(nTetrahedra)

    ! If the 'cellZones' file is missing, allocate only one 'cellZones' structure and return.
    if (.not. self % getElementZonesFile()) then
      call self % setElementZonesNumber(1)
      call elementZones % allocateShelf(1)
      call elementZones % initElementZone('Default', 1, 1, self % nElements)
    
    else
      ! If reached here import cell zones.
      call self % initElementZoneShelf(elementZones, folderPath)

    end if

  end subroutine importMesh

  !! Subroutine 'initElementShelf'
  !!
  !! Basic description:
  !!   Initialises the shelf from the 'owner' and 'neighbour' files.
  !!
  !! Arguments:
  !!   folderPath [in]     -> Path of the folder containing the OpenFOAM mesh files.
  !!   nFaces [in]         -> Number of faces in the mesh.
  !!   nInternalFaces [in] -> Number of internal faces in the mesh.
  !!
  subroutine initElementShelf(self, centroids, elements, faces, vertices, folderPath, concaveElementIdxs)
    class(OpenFOAMMesh), intent(inout)             :: self
    class(vertexShelf), intent(inout)              :: centroids, vertices
    class(elementShelf), intent(inout)             :: elements
    class(faceShelf), intent(inout)                :: faces
    character(*), intent(in)                       :: folderPath
    integer(shortInt), dimension(:), allocatable, intent(out) :: concaveElementIdxs
    integer(shortInt)                              :: i, j, elementIdx, vertexIdx
    integer(shortInt), parameter                   :: unit = 10
    integer(shortInt), dimension(:), allocatable   :: elementIdxs, vertexIdxs
    character(100)                                 :: string
    type(elementInfo), dimension(self % nElements) :: elementInfos
    character(:), allocatable                      :: type

    ! Allocate concaveElementIdxs to zero-size.
    allocate(concaveElementIdxs(0))

    ! If there is only one element in the mesh simply add all the faces and vertices to this element, check convexity and return.
    if (self % nElements == 1) then
      allocate(elementInfos(1) % faceIdxs(self % nFaces))
      do i = 1, self % nFaces
        elementInfos(1) % faceIdxs(i) = i
        call faces % addElementIdxToFace(i, 1)

      end do

      allocate(elementInfos(1) % vertexIdxs(self % nVertices))
      do i = 1, self % nVertices
        elementInfos(1) % vertexIdxs(i) = i
        call vertices % addElementIdxToVertex(i, 1)

      end do

      type = 'Polyhedron'
      if (size(elementInfos(1) % vertexIdxs) == 4) type = 'Tetrahedron'
      call elements % buildElement(1, 1, elementInfos(1) % faceIdxs, elementInfos(1) % vertexIdxs, faces, vertices, type)
      call centroids % initVertex(1, elements % getElementCentroid(1))
      call centroids % addElementIdxToVertex(1, 1)
      if (.not. elements % getElementIsConvex(1)) call append(concaveElementIdxs, 1)
      return

    end if

    ! Open the 'owner' file and read it until a line containing the symbol '(' is encountered.
    call openToRead(unit, folderPath//'owner')
    read(unit, "(a)") string
    do while (index(string(1:len_trim(string)), "(") == 0)
      read(unit, "(a)") string

    end do

    ! Skip one more line and loop over all faces.
    read(unit, "(a)") string
    do i = 1, self % nFaces
      ! Read the current element index and add the current face to this element. Note
      ! that we add one since Fortran starts indexing at one and not zero.
      read(string, *) elementIdx

      elementIdx = elementIdx + 1
      call append(elementInfos(elementIdx) % faceIdxs, i)
      call faces % addElementIdxToFace(i, elementIdx)
      vertexIdxs = faces % getFaceVertexIdxs(i)
      do j = 1, size(vertexIdxs)
        vertexIdx = vertexIdxs(j)
        call append(elementInfos(elementIdx) % vertexIdxs, vertexIdx, .true.)
        call vertices % addElementIdxToVertex(vertexIdx, elementIdx)

      end do

      ! Move onto the next line.
      read(unit, "(a)") string

    end do

    ! Close the 'owner' file.
    close(unit)

    ! If nInternalFaces = 0 we can return early here. Else we need to repeat the above procedure
    ! for the 'neighbour file.
    if (self % nInternalFaces > 0) then
      ! Open the 'neighbour' file and read it until a line containing the symbol '(' is encountered.
      call openToRead(unit, folderPath//'neighbour')
      read(unit, "(a)") string
      do while (index(string(1:len_trim(string)), "(") == 0)
        read(unit, "(a)") string

      end do

      ! Check if the current line contains the symbol ')'. If it does, then all element indices
      ! are written on a single line.
      if (index(string(1:len_trim(string)), ")") > 0) then
        ! Allocate the number of entries in the 'elementIndices' array to the number of internal
        ! faces and copy element indices into this array.
        allocate(elementIdxs(self % nInternalFaces))
        read(string(index(string(1:len_trim(string)), "(") + 1:&
        index(string(1:len_trim(string)), ")") - 1), *) elementIdxs

      else
        ! Skip one more line.
        read(unit, "(a)") string

      end if

      ! Loop over all internal faces.
      do i = 1, self % nInternalFaces
        if (allocated(elementIdxs)) then
          elementIdx = elementIdxs(i)

        else
          ! Read the current element index and move onto the next line.
          read(string, *) elementIdx
          read(unit, "(a)") string

        end if
        ! Update connectivity information.
        elementIdx = elementIdx + 1
        call append(elementInfos(elementIdx) % faceIdxs, -i)
        call faces % addElementIdxToFace(i, elementIdx)
        vertexIdxs = faces % getFaceVertexIdxs(i)
        do j = 1, size(vertexIdxs)
          vertexIdx = vertexIdxs(j)
          call append(elementInfos(elementIdx) % vertexIdxs, vertexIdx, .true.)
          call vertices % addElementIdxToVertex(vertexIdx, elementIdx)

        end do

      end do

      ! Close the 'neighbour' file.
      close(unit)

    end if

    ! Build all elements.
    do i = 1, self % nElements
      type = 'Polyhedron'
      if (size(elementInfos(i) % vertexIdxs) == 4) type = 'Tetrahedron'
      call elements % buildElement(i, i, elementInfos(i) % faceIdxs, elementInfos(i) % vertexIdxs, faces, vertices, type)
      call centroids % initVertex(i, elements % getElementCentroid(i))
      call centroids % addElementIdxToVertex(i, i)
      if (.not. elements % getElementIsConvex(i)) call append(concaveElementIdxs, i)

    end do

  end subroutine initElementShelf

  !! Subroutine 'initCellZoneShelf'
  !!
  !! Basic description:
  !!   Initialises the shelf from the cellZones file.
  !!
  !! Notes:
  !!   Due to the structure of the cellZones file the subroutine first scans the file a first time
  !!   to obtain the name of each cell zone and then rewinds it to read the indices of the elements
  !!   contained in each zone.
  !!
  !! Arguments:
  !!   folderPath [in] -> Path of the folder containing the OpenFOAM mesh files.
  !!
  subroutine initElementZoneShelf(self, elementZones, folderPath)
    class(OpenFOAMMesh), intent(inout)           :: self
    class(cellZoneShelf), intent(inout)          :: elementZones
    character(*), intent(in)                     :: folderPath
    integer(shortInt)                            :: i, j, unit = 10, nElements, nElementZones
    integer(shortInt), allocatable, dimension(:) :: elementIndices
    character(100)                               :: string
    character(:), allocatable                    :: name
    logical(defBool)                             :: singleLine

    ! Open the 'cellZones' data file and read it. Skip lines until a blank line is encountered.
    call openToRead(unit, folderPath//'cellZones')
    read(unit, "(a)") string
    do while (len_trim(string) > 0)
      read(unit, "(a)") string
    end do
    
    ! Read the current line and copy the number of cell zones into the variable 'nCellZones'.
    read(unit, "(a)") string
    read(string, *) nElementZones
    
    ! Allocate memory and loop through all cell zones.
    call self % setElementZonesNumber(nElementZones)
    call elementZones % allocateShelf(nElementZones)
    do i = 1, nElementZones
      ! Initialise singleLine = .false. and skip lines until a '{' is encountered.
      singleLine = .false.
      do while (index(string, "{") == 0)
        read(unit, "(a)") string

      end do

      ! Go back to the previous line and read the name of the current cell zone.
      backspace(unit)
      read(unit, "(a)") string
      name = trim(string)
      
      ! Skip lines until the word 'cellLabels' is encoutered.
      do while (index(string, "cellLabels") == 0)
        read(unit, "(a)") string

      end do

      ! If the line contains a ')' then all the elements contained in the current cell zones are
      ! written on the same line.
      if (index(string, ")") > 0) singleLine = .true.

      ! Read the number of elements in the current cell zone depending on whether they are all
      ! written a single line or not.
      if (singleLine) then
        read(string(index(string, ">") + 2:index(string, "(") - 1), *) nElements

      else
        ! Skip one line and retrieve the number of elements in the cell zone..
        read(unit, "(a)") string
        read(string, *) nElements

      end if

      ! Allocate memory.
      allocate(elementIndices(nElements))

      ! Now read the indices of the elements in the current cell zone. Again this depends on
      ! whether they are all written on a single line or not.
      if (singleLine) then
        ! Retrieve the element indices.
        read(string(index(string, "(") + 1:index(string, ")") - 1), *) elementIndices

      else
        ! Skip two lines and retrieve the index of each element in the cell zone line by line..
        do j = 1, 2
          read(unit, "(a)") string

        end do
        do j = 1, nElements
          read(string, *) elementIndices(j)
          read(unit, "(a)") string

        end do

      end if

      ! Set the start and end elements in the cell zone (note that we encrement by one since Fortran
      ! starts indexing at one instead of zero) and reset elementIndices array.
      call elementZones % initElementZone(name, i, elementIndices(1) + 1, elementIndices(nElements) + 1)
      deallocate(elementIndices)

    end do
    
    ! Close the 'cellZones' file.
    close(unit)

  end subroutine initElementZoneShelf

  !! Subroutine 'initFaceShelf'
  !!
  !! Basic description:
  !!   Initialises the shelf from the 'faces' file.
  !!
  !! Detailed description:
  !!   Opens the 'faces' file and reads it until a line with the symbol ')' is encountered, which 
  !!   marks the beginning of the faces' data listing. In the 'faces' file, each line corresponds to
  !!   a single face. The subroutine then retrieves the number of vertices in each face and their
  !!   indices sequentially line-by-line.
  !!
  !! Arguments:
  !!   folderPath [in]     -> Path of the folder containing the mesh files.
  !!   nInternalFaces [in] -> Number of internal faces in the mesh.
  !!
  subroutine initFaceShelf(self, faces, vertices, folderPath)
    class(OpenFOAMMesh), intent(inout)           :: self
    class(faceShelf), intent(inout)              :: faces
    class(vertexShelf), intent(inout)            :: vertices
    character(*), intent(in)                     :: folderPath
    integer(shortInt), parameter                 :: unit = 10
    integer(shortInt)                            :: i, j, nVertices
    integer(shortInt), dimension(:), allocatable :: vertexIdxs
    character(100)                               :: string
    character(:), allocatable                    :: type

    ! Open the 'faces' data file and read it until a line containing the symbol ')' is encountered.
    call openToRead(unit, folderPath//'faces')
    read(unit, "(a)") string
    do while (index(string(1:len_trim(string)), ")") == 0)
      read(unit, "(a)") string

    end do
    
    ! Loop through all faces in the shelf.
    do i = 1, self % nFaces
      ! Skip lines if they are blank (can happen for faces with a large number of vertices).
      if (len_trim(string) == 0) read(unit, "(a)") string
      
      ! Check if all the vertices are listed on a single line and retrieve nVertices accordingly.
      if (index(string(1:len_trim(string)), "(") > 0) then
        read(string(1:index(string, "(") - 1), *) nVertices

      else
        read(string, *) nVertices
        read(unit, "(a)") string

      end if

      ! Allocate memory.
      allocate(vertexIdxs(nVertices))
      
      ! If all vertices are listed on a single line copy them directly.
      if (index(string(1:len_trim(string)), ")") > 0) then
        read(string(index(string, "(") + 1:len_trim(string) - 1), *) vertexIdxs

      else
        ! Skip one line and read each vertex index line-by-line.
        read(unit, "(a)") string
        do j = 1, nVertices
          read(string, *) vertexIdxs(j)
          read(unit, "(a)") string

        end do
        
        ! Skip one line.
        read(unit, "(a)") string

      end if
      
      ! Add one to the vertices indices since Fortran starts indexing at one rather than zero and
      ! update mesh connectivity information.
      vertexIdxs = vertexIdxs + 1
      do j = 1, size(vertexIdxs)
        call vertices % addFaceIdxToVertex(vertexIdxs(j), i)

      end do

      ! Initialise face in the shelf.
      type = 'Polygon'
      if (nVertices == 3) type = 'Triangle'
      call faces % buildFace(i, i, i > self % nInternalFaces, vertexIdxs, vertices, type)
      
      ! Free memory and move onto the next line.
      deallocate(vertexIdxs)
      read(unit, "(a)") string

    end do
    
    ! Close the 'faces' file.
    close(unit)

  end subroutine initFaceShelf

  !! Subroutine 'initVertexShelf'
  !!
  !! Basic description:
  !!   Initialises the shelf from the 'points' file.
  !!
  !! Arguments:
  !!   folderPath [in] -> Path to the folder containing the mesh files.
  !!
  subroutine initVertexShelf(self, vertices, folderPath)
    class(OpenFOAMMesh), intent(in)   :: self
    class(vertexShelf), intent(inout) :: vertices
    character(*), intent(in)          :: folderPath
    integer(shortInt)                 :: i, j, leftBracketIdx, rightBracketIdx
    integer(shortInt), parameter      :: unit = 10
    real(defReal), dimension(3)       :: coordinates, offset
    real(defReal), dimension(6)       :: extremalCoordinates
    logical(defBool)                  :: singleLine
    character(150)                    :: string ! Note: here the string is longer than usual to deal
                                                ! with cases when all vertices are written on a
                                                ! single line.

    ! Open the 'points' file and read it until a line containing the symbol ')' is encountered.
    call openToRead(unit, folderPath//'points')
    read(unit, "(a)") string
    do while (index(string(1:len_trim(string)), ")") == 0)
      read(unit, "(a)") string

    end do

    ! Initialise singleLine = .false. and check if all vertices are written on a single line.
    ! If yes, update singleLine and erase the leftmost bracket from the string.
    singleLine = .false.
    if (string(2:2) == '(') then
      singleLine = .true.
      string(1:1) = ''

    end if

    ! Retrieve shelf offset then loop over all vertices.
    offset = vertices % getOffset()
    extremalCoordinates = vertices % getExtremalCoordinates()
    do i = 1, self % nVertices
      if (singleLine) then
        ! Locate the leftmost and rightmost brackets. Read the coordinates between the
        ! brackets and remove them from the string.
        leftBracketIdx = index(string, '(')
        rightBracketIdx = index(string, ')')
        read(string(leftBracketIdx + 1:rightBracketIdx - 1), *) coordinates
        string(leftBracketIdx:leftBracketIdx) = ''
        string(rightBracketIdx:rightBracketIdx) = ''

      else
        ! Read the coordinates of the current vertex and move onto the next line.
        read(string(2:len_trim(string) - 1), *) coordinates
        read(unit, "(a)") string

      end if

      ! Set the index and coordinates of the current vertex. Apply offset in the process.
      call vertices % initVertex(i, coordinates + offset)

      ! Update extremal coordinates if necessary.
      do j = 1, 3
        if (coordinates(j) < extremalCoordinates(j)) extremalCoordinates(j) = coordinates(j)
        if (coordinates(j) > extremalCoordinates(j + 3)) extremalCoordinates(j + 3) = coordinates(j)

      end do

    end do
    ! Update extremal coordinates in the shelf.
    call vertices % setExtremalCoordinates(extremalCoordinates)

    ! Close the 'points' file.
    close(unit)

  end subroutine initVertexShelf

  !! Subroutine 'kill'
  !!
  !! Basic description:
  !!   Returns to an unitialised state.
  !!
  elemental subroutine kill(self)
    class(OpenFOAMMesh), intent(inout) :: self

    ! Call unstructuredMesh procedure.
    call kill_super(self)

  end subroutine kill

end module OpenFOAMMesh_class