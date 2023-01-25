      subroutine StripSpaces(string)
      character(len=*) :: string
      integer :: stringLen 
      integer :: last, actual
  
       stringLen = len (string)
        last = 1
      actual = 1
    
      do while (actual.LT.stringLen)
        if (string(last:last) == ' ') then
            actual = actual + 1
            string(last:last) = string(actual:actual)
            string(actual:actual) = ' '
        else
            last = last + 1
            if (actual.lt. last) then
                actual = last
            end if
        endif
        end do

        end subroutine
