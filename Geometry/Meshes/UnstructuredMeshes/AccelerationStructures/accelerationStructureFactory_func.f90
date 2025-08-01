module accelerationStructureFactory_func

  use accelerationStructure_inter,   only : accelerationStructure, initAccelerationStructurePayload
  use dictionary_class,              only : dictionary
  use genericProcedures,             only : fatalError
  use noAcceleration_class,          only : noAcceleration
  use numPrecision
  use octreeAcceleration_class,      only : octreeAcceleration
  use patchSingleAcceleration_class, only : patchSingleAcceleration
  use topologicalObjectShelf_class,  only : topologicalObjectShelf

  implicit none
  private

  ! List which contains acceptable types of acceleration structures.
  ! NOTE: It is necessary to adjust trailing blanks so all entries have the same length.
  character(nameLen), dimension(*), parameter :: AVAILABLE_ACCELERATIONS = ['none       ', &
                                                                            'octree     ', &
                                                                            'patchSearch']

  ! Public interface.
  public :: newAccelerationStructurePtr

contains
  !!
  !!
  !!
  subroutine newAccelerationStructurePtr(dict, edges, elements, faces, vertices, ptr)
    class(dictionary), intent(in)                      :: dict
    type(topologicalObjectShelf), target, intent(in)   :: edges, elements, faces, vertices
    class(accelerationStructure), pointer, intent(out) :: ptr
    type(initAccelerationStructurePayload)             :: payload
    character(nameLen)                                 :: type
    character(*), parameter :: here = 'newAccelerationStructurePtr (accelerationStructureFactory_func.f90)'

    ! Get type of acceleration structure from dictionary. Default to 'none' if user has not specified one.
    type = 'none'
    if (dict % isPresent('accelerationMethod')) then
      payload % dict => dict % getDictPtr('accelerationMethod')
      call payload % dict % getOrDefault(type, 'type', 'none')
      payload % edges => edges
      payload % elements => elements
      payload % faces => faces
      payload % vertices => vertices

    end if
    select case(type)
      case('none')
          allocate(noAcceleration :: ptr)

      case('octree')
          allocate(octreeAcceleration :: ptr)

      case('patchSearch')
          allocate(patchSingleAcceleration :: ptr)

      case default
          print '(A)', 'AVAILABLE ACCELERATION STRUCTURES: '
          print '(A)', AVAILABLE_ACCELERATIONS
          call fatalError(here, 'Unrecognised type of acceleration structure: '//trim(type)//'.')

    end select

    ! Initialise acceleration structure.
    call ptr % init(payload)

  end subroutine newAccelerationStructurePtr

end module accelerationStructureFactory_func