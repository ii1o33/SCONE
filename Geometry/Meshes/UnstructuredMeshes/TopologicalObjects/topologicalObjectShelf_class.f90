module topologicalObjectShelf_class

  use edge_class,                    only : edge, edgeBox
  use element_class,                 only : buildElementPayload, element, elementBox
  use extentTopologicalObject_inter, only : buildExtentTopologicalObjectPayload
  use face_class,                    only : buildFacePayload, face, faceBox
  use genericProcedures,             only : fatalError, numToChar, quickSort
  use iso_fortran_env,               only : int64
  use longIntMap_class,              only : longIntMap
  use numPrecision
  use topologicalObject_inter,       only : buildTopologicalObjectPayload, topologicalObject, topologicalObjectBox
  use topologicalObjectFactory_func, only : newTopologicalObjectBox
  use vertex_class,                  only : buildVertexPayload, vertex, vertexBox

  implicit none
  private

  !!
  !!
  !!
  type, public                                            :: topologicalObjectShelf
    private
    type(topologicalObjectBox), dimension(:), allocatable :: shelf
    type(longIntMap)                                      :: idxMap
    integer(shortInt)                                     :: nObjects = 0
  contains
    procedure          :: addObject
    generic, private   :: generateKey => generateKey_idx, generateKey_payload, generateKey_vertices
    procedure, private :: generateKey_idx
    procedure, private :: generateKey_payload
    procedure, private :: generateKey_vertices
    procedure          :: getActiveObjectIdxs
    generic            :: getEdgeBox => getEdgeBox_shortInt, getEdgeBox_shortIntArray
    procedure, private :: getEdgeBox_shortInt
    procedure, private :: getEdgeBox_shortIntArray
    generic            :: getElementBox => getElementBox_shortInt, getElementBox_shortIntArray
    procedure, private :: getElementBox_shortInt
    procedure, private :: getElementBox_shortIntArray
    generic            :: getFaceBox => getFaceBox_shortInt, getFaceBox_shortIntArray
    procedure, private :: getFaceBox_shortInt
    procedure, private :: getFaceBox_shortIntArray
    generic            :: getObjectBox => getObjectBox_shortInt, getObjectBox_shortIntArray
    procedure, private :: getObjectBox_shortInt
    procedure, private :: getObjectBox_shortIntArray
    generic            :: getObjectBoundingBoxBounds => getObjectBoundingBoxBounds_shortInt, &
                                                        getObjectBoundingBoxBounds_shortIntArray
    procedure, private :: getObjectBoundingBoxBounds_shortInt
    procedure, private :: getObjectBoundingBoxBounds_shortIntArray
    generic            :: getObjectCentroid => getObjectCentroid_shortInt, getObjectCentroid_shortIntArray
    procedure, private :: getObjectCentroid_shortInt
    procedure, private :: getObjectCentroid_shortIntArray
    procedure          :: getObjectIdxOrDefault
    procedure          :: getObjectsNumber
    procedure          :: getShelf
    procedure          :: getSize
    generic            :: getVertexBox => getVertexBox_shortInt, getVertexBox_shortIntArray
    procedure, private :: getVertexBox_shortInt
    procedure, private :: getVertexBox_shortIntArray
    procedure          :: init
    generic            :: initObject => initMultipleObjects, initSingleObject
    procedure, private :: initMultipleObjects
    procedure, private :: initSingleObject
    procedure          :: kill
    procedure          :: shrink
  end type topologicalObjectShelf

contains
  !!
  !!
  !!
  subroutine addObject(self, box)
    class(topologicalObjectShelf), intent(inout)          :: self
    class(topologicalObjectBox), intent(in)               :: box
    integer(shortInt)                                     :: currentSize
    type(topologicalObjectBox), dimension(:), allocatable :: tempShelf

    ! First check if the shelf is allocated.
    if (allocated(self % shelf)) then
      ! Check if shelf is full and double its size if so.
      currentSize = size(self % shelf)
      if (self % nObjects == currentSize) then
        allocate(tempShelf(2 * currentSize))
        tempShelf(1:currentSize) = self % shelf
        call move_alloc(tempShelf, self % shelf)

      end if

    else
      ! Allocate shelf with a reasonable initial size.
      allocate(self % shelf(8))

    end if
    ! Update self % nObjects and add the new item at the next available position.
    self % nObjects = self % nObjects + 1
    self % shelf(self % nObjects) = box

  end subroutine addObject

  !!
  !!
  !!
  function generateKey_idx(self, idx) result(key)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt), intent(in)             :: idx
    integer(longInt)                          :: key

    key = int(idx, int64)

  end function generateKey_idx

  !!
  !!
  !!
  function generateKey_payload(self, payload) result(key)
    class(topologicalObjectShelf), intent(in)           :: self
    class(buildTopologicalObjectPayload), intent(in)    :: payload
    integer(longInt)                                    :: key
    class(buildExtentTopologicalObjectPayload), pointer :: extentPayloadPtr
    type(buildVertexPayload), pointer                   :: vertexPayloadPtr
    character(*), parameter                             :: here = 'generateKey_payload (topologicalObjectShelf_class.f90)'

    select type(ptr => payload)
      class is(buildExtentTopologicalObjectPayload)
        extentPayloadPtr => ptr
        key = self % generateKey_vertices(extentPayloadPtr % vertices)

      type is(buildVertexPayload)
        vertexPayloadPtr => ptr
        key = self % generateKey_idx(vertexPayloadPtr % idx)

      class default
        call fatalError(here, 'Invalid payload type.')

    end select

  end function generateKey_payload

  !!
  !!
  !!
  function generateKey_vertices(self, vertices) result(key)
    class(topologicalObjectShelf), intent(in)    :: self
    type(vertexBox), dimension(:), intent(in)    :: vertices
    integer(longInt)                             :: key
    integer(shortInt)                            :: i, nVertices
    integer(shortInt), dimension(size(vertices)) :: vertexIdxs
    integer(longInt), parameter                  :: prime = 31_longInt
    character(*), parameter                      :: here = 'generateKey_vertices (topologicalObjectShelf_class.f90)'

    ! Get vertex indices and sort them.
    nVertices = size(vertices)
    if (nVertices == 0) call fatalError(here, 'Attempting to generate a key for a zero-sized array.')
    do i = 1, nVertices
      vertexIdxs(i) = vertices(i) % ptr % getIdx()

    end do
    call quickSort(vertexIdxs)

    ! Generate key.
    select case(nVertices)
      case(1)
        key = self % generateKey_idx(vertexIdxs(1))

      case(2)
        key = ishft(int(vertexIdxs(1), int64), 32) + int(vertexIdxs(2), int64)

      case default
        key = 0
        do i = 1, nVertices
          key = key * prime + vertexIdxs(i)

        end do

    end select

  end function generateKey_vertices

  !!
  !!
  !!
  function getActiveObjectIdxs(self) result(activeObjectIdxs)
    class(topologicalObjectShelf), intent(in)    :: self
    integer(shortInt), dimension(:), allocatable :: activeObjectIdxs, tempActiveObjectIdxs
    integer(shortInt)                            :: i, nActiveObjects

    ! Initialise nActiveObjects = 0
    nActiveObjects = 0
    allocate(activeObjectIdxs(self % nObjects))
    do i = 1, self % nObjects
      if (.not. self % shelf(i) % ptr % getIsActive()) cycle
      nActiveObjects = nActiveObjects + 1
      activeObjectIdxs(nActiveObjects) = self % shelf(i) % ptr % getIdx()

    end do

    ! Now resize activeObjectIdxs to correct size.
    allocate(tempActiveObjectIdxs(nActiveObjects))
    tempActiveObjectIdxs = activeObjectIdxs(1:nActiveObjects)
    call move_alloc(tempActiveObjectIdxs, activeObjectIdxs)

  end function getActiveObjectIdxs

  !!
  !!
  !!
  function getEdgeBox_shortInt(self, idx) result(box)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt), intent(in)             :: idx
    type(edgeBox)                             :: box
    type(topologicalObjectBox)                :: objectBox
    character(*), parameter                   :: here = 'getEdgeBox_shortInt (topologicalObjectShelf_class.f90)'

    objectBox = self % shelf(idx)
    select type(ptr => objectBox % ptr)
      type is(edge)
        box % ptr => ptr

      class default
        call fatalError(here, 'Topological object with index: '//numToChar(ptr % getIdx())//' is not an edge.')

    end select

  end function getEdgeBox_shortInt

  !!
  !!
  !!
  function getEdgeBox_shortIntArray(self, idxs) result(boxes)
    class(topologicalObjectShelf), intent(in)   :: self
    integer(shortInt), dimension(:), intent(in) :: idxs
    type(edgeBox), dimension(size(idxs))        :: boxes
    integer(shortInt)                           :: i
    type(topologicalObjectBox)                  :: objectBox
    character(*), parameter                     :: here = 'getEdgeBox_shortIntArray (topologicalObjectShelf_class.f90)'

    do i = 1, size(idxs)
      objectBox = self % shelf(idxs(i))
      select type(ptr => objectBox % ptr)
        type is(edge)
          boxes(i) % ptr => ptr

        class default
          call fatalError(here, 'Topological object with index: '//numToChar(ptr % getIdx())//' is not an edge.')

      end select

    end do

  end function getEdgeBox_shortIntArray

  !!
  !!
  !!
  function getElementBox_shortInt(self, idx) result(box)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt), intent(in)             :: idx
    type(elementBox)                          :: box
    type(topologicalObjectBox)                :: objectBox
    character(*), parameter                   :: here = 'getElementBox_shortInt (topologicalObjectShelf_class.f90)'

    objectBox = self % shelf(idx)
    select type(ptr => objectBox % ptr)
      type is(element)
        box % ptr => ptr

      class default
        call fatalError(here, 'Topological object with index: '//numToChar(ptr % getIdx())//' is not an element.')

    end select

  end function getElementBox_shortInt

  !!
  !!
  !!
  function getElementBox_shortIntArray(self, idxs) result(boxes)
    class(topologicalObjectShelf), intent(in)   :: self
    integer(shortInt), dimension(:), intent(in) :: idxs
    type(elementBox), dimension(size(idxs))     :: boxes
    integer(shortInt)                           :: i
    type(topologicalObjectBox)                  :: objectBox
    character(*), parameter                     :: here = 'getElementBox_shortIntArray (topologicalObjectShelf_class.f90)'

    do i = 1, size(idxs)
      objectBox = self % shelf(idxs(i))
      select type(ptr => objectBox % ptr)
        type is(element)
          boxes(i) % ptr => ptr

        class default
          call fatalError(here, 'Topological object with index: '//numToChar(ptr % getIdx())//' is not an element.')

      end select

    end do

  end function getElementBox_shortIntArray

  !!
  !!
  !!
  function getFaceBox_shortInt(self, idx) result(box)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt), intent(in)             :: idx
    type(faceBox)                             :: box
    type(topologicalObjectBox)                :: objectBox
    character(*), parameter                   :: here = 'getFaceBox_shortInt (topologicalObjectShelf_class.f90)'

    objectBox = self % shelf(idx)
    select type(ptr => objectBox % ptr)
      type is(face)
        box % ptr => ptr

      class default
        call fatalError(here, 'Topological object with index: '//numToChar(ptr % getIdx())//' is not a face.')

    end select

  end function getFaceBox_shortInt

  !!
  !!
  !!
  function getFaceBox_shortIntArray(self, idxs) result(boxes)
    class(topologicalObjectShelf), intent(in)   :: self
    integer(shortInt), dimension(:), intent(in) :: idxs
    type(faceBox), dimension(size(idxs))        :: boxes
    integer(shortInt)                           :: i
    type(topologicalObjectBox)                  :: objectBox
    character(*), parameter                     :: here = 'getFaceBox_shortIntArray (topologicalObjectShelf_class.f90)'

    do i = 1, size(idxs)
      objectBox = self % shelf(idxs(i))
      select type(ptr => objectBox % ptr)
        type is(face)
          boxes(i) % ptr => ptr

        class default
          call fatalError(here, 'Topological object with index: '//numToChar(ptr % getIdx())//' is not a face.')

      end select

    end do

  end function getFaceBox_shortIntArray

  !!
  !!
  !!
  function getObjectBox_shortInt(self, idx) result(box)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt), intent(in)             :: idx
    type(topologicalObjectBox)                :: box

    box = self % shelf(idx)

  end function getObjectBox_shortInt

  !!
  !!
  !!
  function getObjectBox_shortIntArray(self, idxs) result(boxes)
    class(topologicalObjectShelf), intent(in)         :: self
    integer(shortInt), dimension(:), intent(in)       :: idxs
    type(topologicalObjectBox), dimension(size(idxs)) :: boxes
    integer(shortInt)                                 :: i

    do i = 1, size(idxs)
      boxes(i) = self % shelf(idxs(i))

    end do

  end function getObjectBox_shortIntArray

  !!
  !!
  !!
  function getVertexBox_shortInt(self, idx) result(box)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt), intent(in)             :: idx
    type(vertexBox)                           :: box
    type(topologicalObjectBox)                :: objectBox
    character(*), parameter                   :: here = 'getVertexBox_shortInt (topologicalObjectShelf_class.f90)'

    objectBox = self % shelf(idx)
    select type(ptr => objectBox % ptr)
      type is(vertex)
        box % ptr => ptr

      class default
        call fatalError(here, 'Topological object with index: '//numToChar(ptr % getIdx())//' is not a vertex.')

    end select

  end function getVertexBox_shortInt

  !!
  !!
  !!
  function getVertexBox_shortIntArray(self, idxs) result(boxes)
    class(topologicalObjectShelf), intent(in)   :: self
    integer(shortInt), dimension(:), intent(in) :: idxs
    type(vertexBox), dimension(size(idxs))      :: boxes
    integer(shortInt)                           :: i
    type(topologicalObjectBox)                  :: objectBox
    character(*), parameter                     :: here = 'getVertexBox_shortIntArray (topologicalObjectShelf_class.f90)'

    do i = 1, size(idxs)
      objectBox = self % shelf(idxs(i))
      select type(ptr => objectBox % ptr)
        type is(vertex)
          boxes(i) % ptr => ptr

        class default
          call fatalError(here, 'Topological object with index: '//numToChar(ptr % getIdx())//' is not a vertex.')

      end select

    end do

  end function getVertexBox_shortIntArray

  !!
  !!
  !!
  function getObjectBoundingBoxBounds_shortInt(self, idx) result(bounds)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt), intent(in)             :: idx
    real(defReal), dimension(3, 2)            :: bounds

    bounds = self % shelf(idx) % ptr % getBoundingBoxBounds()

  end function getObjectBoundingBoxBounds_shortInt

  !!
  !!
  !!
  function getObjectBoundingBoxBounds_shortIntArray(self, idxs) result(bounds)
    class(topologicalObjectShelf), intent(in)   :: self
    integer(shortInt), dimension(:), intent(in) :: idxs
    real(defReal), dimension(3, 2 * size(idxs)) :: bounds
    integer(shortInt)                           :: i

    do i = 1, size(idxs)
      bounds(:, 2 * (i - 1) + 1:2 * i) = self % shelf(idxs(i)) % ptr % getBoundingBoxBounds()

    end do

  end function getObjectBoundingBoxBounds_shortIntArray

  !!
  !!
  !!
  function getObjectCentroid_shortInt(self, idx) result(centroid)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt), intent(in)             :: idx
    real(defReal), dimension(3)               :: centroid

    centroid = self % shelf(idx) % ptr % getCentroid()

  end function getObjectCentroid_shortInt

  !!
  !!
  !!
  function getObjectCentroid_shortIntArray(self, idxs) result(centroids)
    class(topologicalObjectShelf), intent(in)   :: self
    integer(shortInt), dimension(:), intent(in) :: idxs
    real(defReal), dimension(3, size(idxs))     :: centroids
    integer(shortInt)                           :: i

    do i = 1, size(idxs)
      centroids(:, i) = self % shelf(idxs(i)) % ptr % getCentroid()

    end do

  end function getObjectCentroid_shortIntArray

  !!
  !!
  !!
  function getShelf(self) result(shelf)
    class(topologicalObjectShelf), intent(in)             :: self
    type(topologicalObjectBox), dimension(:), allocatable :: shelf

    if (allocated(self % shelf)) then
      shelf = self % shelf

    else
      allocate(shelf(0))

    end if

  end function getShelf

  !!
  !!
  !!
  function getObjectIdxOrDefault(self, vertices, default) result(idx)
    class(topologicalObjectShelf), intent(in) :: self
    type(vertexBox), dimension(:), intent(in) :: vertices
    integer(shortInt), intent(in)             :: default
    integer(shortInt)                         :: idx

    idx = self % idxMap % getOrDefault(self % generateKey(vertices), default)

  end function getObjectIdxOrDefault

  !!
  !!
  !!
  elemental function getObjectsNumber(self) result(nObjects)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt)                         :: nObjects

    nObjects = self % nObjects

  end function getObjectsNumber

  !!
  !!
  !!
  elemental function getSize(self) result(shelfSize)
    class(topologicalObjectShelf), intent(in) :: self
    integer(shortInt)                         :: shelfSize

    shelfSize = size(self % shelf)

  end function getSize

  !!
  !!
  !!
  subroutine init(self, payloads)
    class(topologicalObjectShelf), intent(inout)                      :: self
    class(buildTopologicalObjectPayload), dimension(:), intent(inout) :: payloads
    integer(shortInt)                                                 :: i, nPayloads

    ! Allocate shelf then add all payloads.
    nPayloads = size(payloads)
    allocate(self % shelf(nPayloads))
    do i = 1, nPayloads
      call newTopologicalObjectBox(payloads(i), self % shelf(i))
      self % nObjects = self % nObjects + 1

      ! Add key to idxMap.
      call self % idxMap % add(self % generateKey(payloads(i)), payloads(i) % idx)

    end do

  end subroutine init

  !!
  !!
  !!
  subroutine initMultipleObjects(self, payloads)
    class(topologicalObjectShelf), intent(inout)                      :: self
    class(buildTopologicalObjectPayload), dimension(:), intent(inout) :: payloads
    integer(shortInt)                                                 :: i
    type(topologicalObjectBox)                                        :: box

    ! Loop through all payloads supplied.
    do i = 1, size(payloads)
      ! Initialise a new object.
      call newTopologicalObjectBox(payloads(i), box)

      ! Add new object to shelf.
      call self % addObject(box)

      ! Add new key to idxMap.
      call self % idxMap % add(self % generateKey(payloads(i)), payloads(i) % idx)

    end do

  end subroutine initMultipleObjects

  !!
  !!
  !!
  subroutine initSingleObject(self, payload)
    class(topologicalObjectShelf), intent(inout)          :: self
    class(buildTopologicalObjectPayload), intent(inout)   :: payload
    type(topologicalObjectBox)                            :: box

    ! Initialise a new object.
    call newTopologicalObjectBox(payload, box)

    ! Add new object to shelf. First check if the shelf is allocated.
    call self % addObject(box)

    ! Add new key to idxMap.
    call self % idxMap % add(self % generateKey(payload), payload % idx)

  end subroutine initSingleObject

  !!
  !!
  !!
  subroutine kill(self)
    class(topologicalObjectShelf), intent(inout) :: self
    integer(shortInt)                            :: i

    ! Local.
    if (allocated(self % shelf)) then
      do i = 1, size(self % shelf)
        call self % shelf(i) % ptr % kill()
        deallocate(self % shelf(i) % ptr)

      end do
      deallocate(self % shelf)

    end if
    call self % idxMap % kill()
    self % nObjects = 0

  end subroutine kill

  !!
  !!
  !!
  subroutine shrink(self)
    class(topologicalObjectShelf), intent(inout)          :: self
    type(topologicalObjectBox), dimension(:), allocatable :: tempShelf
    character(*), parameter                               :: here = 'shrinkShelf (topologicalObjectShelf_inter.f90)'

    ! First check if shelf is unallocated.
    if (.not. allocated(self % shelf)) call fatalError(here, 'Attempting to shrink an unallocated shelf.')
    allocate(tempShelf(self % nObjects))
    tempShelf = self % shelf(1:self % nObjects)
    call move_alloc(tempShelf, self % shelf)

  end subroutine shrink

end module topologicalObjectShelf_class