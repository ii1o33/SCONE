module dictionaryType

  use numPrecision

  implicit none
  private
  public                            :: dictEntryShortInteger, dictionaryShortInteger, addEntry, getValue, hasKey

  ! Add other dictEntry and dictionary types if needed.

  type :: dictEntryShortInteger
    integer(shortInt)                           :: key
    integer(shortInt)                           :: value
  end type dictEntryShortInteger

  type :: dictionaryShortInteger
    type(dictEntryShortInteger), allocatable    :: entries(:)
  contains
    procedure :: addEntry
    procedure :: getValue
    procedure :: hasKey
  end type dictionaryShortInteger

contains









  ! Add or update an entry
  subroutine addEntry(self, key, value)
    class(dictionaryShortInteger), intent(inout) :: self
    integer(shortInt), intent(in)                :: key
    integer(shortInt), intent(in)                :: value
    integer(shortInt)                            :: i, n
    type(dictEntryShortInteger), allocatable     :: temp(:)

    ! Update value if key already exists
    if (allocated(self%entries)) then
      do i = 1, size(self%entries)
        if (self%entries(i)%key == key) then
          self%entries(i)%value = value
          return
        end if
      end do
    end if

    ! Append new entry !{might be able to make it more efficient?}
    n = 0
    if (allocated(self%entries)) then
      n = size(self%entries)
      allocate(temp(n))
      temp = self%entries
      deallocate(self%entries)
      allocate(self%entries(n + 1))
      self%entries(1:n) = temp
      deallocate(temp)
    else
      allocate(self%entries(1))
    end if

    self%entries(n + 1)%key = key
    self%entries(n + 1)%value = value
  end subroutine addEntry










  
  function getValue(self, key) result(val)
    class(dictionaryShortInteger), intent(in) :: self
    integer(shortInt), intent(in)             :: key
    integer(shortInt)                         :: val
    integer(shortInt)                         :: i

    do i = 1, size(self%entries)
      if (self%entries(i)%key == key) then
        val = self%entries(i)%value
        return
      end if
    end do

    error stop 'Key not found on the dictionary'
  end function getValue






  function hasKey(self, key) result(found)
    class(dictionaryShortInteger), intent(in) :: self
    integer(shortInt), intent(in)             :: key
    logical                                   :: found
    integer(shortInt)                         :: i

    found = .false.

    if (allocated(self%entries)) then
        do i = 1, size(self%entries)
        if (self%entries(i)%key == key) then
            found = .true.
            exit
        end if
        end do
    end if

    return
  end function hasKey

end module dictionaryType


