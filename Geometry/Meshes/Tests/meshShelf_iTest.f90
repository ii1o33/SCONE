module meshShelf_iTest
  
  use numPrecision
  use dictionary_class, only : dictionary
  use dictParser_func,  only : charToDict
  use mesh_inter,       only : mesh
  use meshShelf_class,  only : meshShelf
  use funit
  
  implicit none
  
  ! Parameters.
  character(*), parameter :: MESHES_DEF = &
  " testMesh1 { id 11; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh1/;} &
   &testMesh2 { id 21; type OpenFOAMMesh; path ./IntegrationTestFiles/Geometry/Meshes/OpenFOAM/testMesh2/;}"
  
  ! Variables.
  type(meshShelf) :: meshes
  
contains
  
  !!
  !! Build shelf.
  !!
@Before
  subroutine setUp()
    type(dictionary) :: dict
    call charToDict(dict, MESHES_DEF)
    call meshes % init(dict)
  end subroutine setUp
  
  !!
  !! Clean after tests.
  !!
@After
  subroutine cleanUp()
    call meshes % kill()
  end subroutine cleanUp
  
  !!
  !! Test shelf.
  !!
@Test
  subroutine test_get()
    class(mesh), pointer :: ptr
    integer(shortInt)    :: idx
    ! Mesh ID 11.
    idx = meshes % getIdx(11)
    ptr => meshes % getPtr(idx)
    @assertEqual(11, ptr % getId())
    @assertEqual(11, meshes % getID(idx))
    ! Mesh ID 21.
    idx = meshes % getIdx(21)
    ptr => meshes % getPtr(idx)
    @assertEqual(21, ptr % getId())
    @assertEqual(21, meshes % getID(idx))    
  end subroutine test_get
end module meshShelf_iTest