!
! Copyright (C) 2017
! Cyrille Barreteau <mailto:cyrille.barreteau@cea.fr>,
! Mathieu Cesar <mailto:mathieu.cesar@cea.fr>,
! Pascal Thibaudeau <mailto:pascal.thibaudeau@cea.fr>.
!
! This software is a computer program whose purpose is DyNaMol.
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and inRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
!
!  string.f90
!  DyNaMol
module string_mod
   use iso_fortran_env, only: error_unit
   use precision_mod, only: il, ip, rl, rp
   implicit none

   !> Working length of character string variables (arbitrary)
   !> 132 is the maximum length of a line free-form Fortran source line
   integer, parameter :: sl = 132

   interface fmt
      procedure :: fmt_txt
      procedure :: fmt_integer
      procedure :: fmt_real4
      procedure :: fmt_real8
   end interface fmt

contains
   !> Eliminate the tab characters and adjust to the left side
   function clean_str(str) result(strout)
      character(len=*) :: str
      character(len=len(str)) :: strin
      character(len=len(str)) :: strout
      integer :: i, stringlen
      integer :: last, actual

      strin = str
      stringlen = len(str)
      last = 1
      actual = 1

      do while (actual < stringLen)
         if (strin(last:last) == ' ') then
            actual = actual + 1
            strin(last:last) = strin(actual:actual)
            strin(actual:actual) = ' '
         else
            last = last + 1
            if (actual < last) &
               actual = last
         end if
      end do

      strout = trim(adjustl(strin))
!    do i=1, len(str)
!      if(ichar(str(i:i)) == 10) then
!        clean_str(i:i) = ' '
!      else
!        clean_str(i:i) = str(i:i)
!      end if
!    end do

!    clean_str = adjustl(clean_str)
   end function clean_str

   function cmplx2str(c) result(str)
      complex(rp), intent(in) :: c
      character(len=:), allocatable :: str

      str = '('//real2str(real(c))//', '//real2str(aimag(c))//')'
   end function cmplx2str

   !> Custom flush subroutine to activate with a flag
   subroutine dynamol_flush(unit)
      integer, intent(in)  :: unit
#if defined(Enable_flush)
      flush (unit=unit)
#endif
   end subroutine dynamol_flush

   !> Convert integer to character string
   function int2str(i, fmt) result(str)
      integer(ip), intent(in) :: i
      character(len=*), intent(in), optional :: fmt
      character(len=il) :: str_temp
      character(len=:), allocatable :: str

      if (present(fmt)) then
         write (str_temp, fmt) i
      else
         write (str_temp, *) i
      end if
      str = trim(adjustl(str_temp))
   end function int2str

   !> Test if a character string is an integer
   function is_str_int(str) result(l)
      character(len=*), intent(in) :: str
      logical :: l
      character, parameter :: char_minus = '-'
      character(len=10), parameter :: char_digit = '0123456789'
      integer :: is, id
      integer :: ls

      ls = len(str)
      if (ls == 0) then
         l = .false.
      elseif (ls == 1) then
         l = .false.
         do id = 1, 10
            l = l .or. (str(1:1) == char_digit(id:id))
            if (l) exit
         end do
      else
         ! Test first character
         l = .false.
         l = l .or. (str(1:1) == char_minus)
         if (.not. l) then
            do id = 1, 10
               l = l .or. (str(1:1) == char_digit(id:id))
               if (l) exit
            end do
         end if
         if (.not. l) return
         ! Test remaining characters
         do is = 2, ls
            l = .false.
            do id = 1, 10
               l = l .or. (str(is:is) == char_digit(id:id))
               if (l) exit
            end do
            if (.not. l) return
         end do
      end if
   end function is_str_int

   !> Test if a character string is an integer
   function is_str_log(str) result(l)
      character(len=*), intent(in) :: str
      logical :: l

      if (str == 't' .or. str == '.true.' .or. str == 'f' .or. str == '.false.') then
         l = .true.
      else
         l = .false.
      end if
   end function is_str_log

   ! Check if a character string is a real
   !function is_str_real(string)

   !> Convert logical to character string
   function log2str(l) result(str)
      logical, intent(in) :: l
      character(len=:), allocatable :: str

      if (l) then
         str = '.true.'
      else
         str = '.false.'
      end if
   end function log2str

   !> Convert lowercase characters to uppercase, copy other characters
   function lower(str1) result(str2)
      character(len=*), intent(in) :: str1
      character(len=len(str1))    :: str2
      character(len=26), parameter :: char_lower = 'abcdefghijklmnopqrstuvwxyz'
      character(len=26), parameter :: char_upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      integer :: is, ic

      str2 = str1
      do is = 1, len(str1)
         do ic = 1, 26
            if (str1(is:is) == char_upper(ic:ic)) then
               str2(is:is) = char_lower(ic:ic)
               exit
            end if
         end do
      end do
   end function lower

   ! Convert real to character string
   function real2str(r, fmt) result(str)
      real(rp), intent(in) :: r
      character(len=*), intent(in), optional :: fmt
      character(len=rl) :: str_temp
      character(len=:), allocatable :: str

      if (present(fmt)) then
         write (str_temp, fmt) r
      else
         write (str_temp, *) r
      end if
      str = trim(adjustl(str_temp))
   end function real2str

   ! Convert real to character string in a fixed form
   function fixedreal2str(r) result(str)
      real(rp), intent(in) :: r
      character(len=rl) :: str_temp
      character(len=:), allocatable :: str

      write (str_temp, "(F18.14)") r
      str = trim(adjustl(str_temp))
   end function fixedreal2str

   ! Convert character string to complex
   !function str2cmplx(str) result(r)

   !> Convert character string to integer
   function str2int(str) result(i)
      character(len=*), intent(in) :: str
      integer :: i

      if (.not. is_str_int(trim(adjustl(str)))) then
         write (error_unit, *) 'string%str2int(): character string ', str, &
            ' is not integer'
         error stop
      end if
      read (str, *) i
   end function str2int

   !> Convert lcharacter string to logical
   function str2log(str) result(l)
      character(len=*), intent(in) :: str
      character(len=len(str)) :: str_temp
      logical :: l

      str_temp = lower(str)
      if (.not. is_str_log(str_temp)) then
         write (error_unit, *) 'string%str2log(): character string ', str, &
            ' is not logical'
         error stop
      end if

      if (str_temp == 't' .or. str_temp == '.true.') then
         l = .true.
      else
         l = .false.
      end if
   end function str2log

   ! Convert character string to real
   !function str2real(str) result(r)

   !> Unique values c in a string array a, index ia is such that c = a(ia) and
   !> index ic is such that a = c(ic)
   subroutine unique_str(a, c, ia, ic)
      character(len=*), dimension(:), intent(in) :: a
      !character(len=len(a)), dimension(:), allocatable, intent(out) :: c
      character(len=*), dimension(:), allocatable, intent(out) :: c
      integer, dimension(:), allocatable, intent(out) :: ia
      integer, dimension(size(a)), intent(out) :: ic
      integer :: sa, sc, ia1, ia2, ic1
      integer, dimension(size(a)) :: index_a
      logical, dimension(size(a)) :: mask

      sa = size(a)
      index_a = (/(ia1, ia1=1, sa)/)
      mask = .false.
      do ia1 = 1, sa - 1
         if (.not. mask(ia1)) then
            do ia2 = ia1 + 1, sa
               mask(ia2) = a(ia1) == a(ia2)
               if (mask(ia2)) index_a(ia2) = index_a(ia1)
            end do
         end if
      end do
      sc = sa - count(mask)
      allocate (c(sc), ia(sc))
      c = pack(a,.not. mask)
      ia = pack(index_a,.not. mask)
      do ia1 = 1, sa
         do ic1 = 1, sc
            if (index_a(ia1) == ia(ic1)) then
               ic(ia1) = ic1
               exit
            end if
         end do
      end do
   end subroutine unique_str

   !> Convert uppercase characters to lowercase, copy other characters
   function upper(str1) result(str2)
      character(len=*), intent(in) :: str1
      character(len=len(str1))    :: str2
      character(len=26), parameter :: char_lower = 'abcdefghijklmnopqrstuvwxyz'
      character(len=26), parameter :: char_upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      integer :: is, ic

      str2 = str1
      do is = 1, len(str1)
         do ic = 1, 26
            if (str1(is:is) == char_lower(ic:ic)) then
               str2(is:is) = char_upper(ic:ic)
               exit
            end if
         end do
      end do
   end function upper

   !> Split 'string' into array of strings separated by 'sep' argument
   subroutine split(string, sep, lstring)
      implicit none
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: sep
      character(len=sl), dimension(:), allocatable, intent(inout) :: lstring
      integer :: length, i, idx0, idx1
      length = count([(string(i:i) .eq. sep, i=1, len(string))]) + 1
      if (allocated(lstring)) deallocate (lstring)
      allocate (lstring(length))
      idx0 = 1
      do i = 1, length
         if (i .eq. length) then
            idx1 = len(string)
         else
            idx1 = max(len(string(:idx0)) - 1 + scan(string(idx0:), sep), idx0)
         end if
         if (string(idx1:idx1) .eq. sep) then
            lstring(i) = string(idx0:idx1 - 1)
         else
            lstring(i) = string(idx0:idx1)
         end if
         idx0 = idx1 + 1
      end do
   end subroutine split

   !> Same as split but function instead
   function fsplit(string, sep)
      implicit none
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: sep
      character(len=sl), dimension(:), allocatable :: fsplit
      call split(string, sep, fsplit)
   end function fsplit

   !> Join elements of an array of strings by sep
   subroutine join(lst_string, sep, string)
      character(len=*), dimension(:), intent(in) :: lst_string
      character(len=*), intent(in), optional :: sep
      character(len=sl), intent(out) :: string
      character(len=sl) :: sep_used
      integer :: i
      sep_used = ''
      if (present(sep)) sep_used = sep
      if (size(lst_string) > 0) then
         string = trim(lst_string(1))
         do i = 2, size(lst_string)
            string = trim(string)//trim(sep_used)//trim(lst_string(i))
         end do
      else
         string = ''
      end if
   end subroutine join

   !> Same as join but function
   function fjoin(lst_string, sep)
      character(len=*), dimension(:), intent(in) :: lst_string
      character(len=*), intent(in), optional :: sep
      character(len=sl) :: fjoin
      ! character(len=sl) :: sep_used
      ! integer :: i
      ! sep_used = ''
      ! if (present(sep)) sep_used = sep
      call join(lst_string, sep, fjoin)
   end function fjoin

   !> Remove all occurrences of a character in string
   subroutine clean_string(string, c)
      character(len=*), intent(inout) :: string
      character(len=1) :: c
      integer :: i
      i = index(string, c)
      do while (i .ne. 0)
         string = string(:i - 1)//string(i + 1:)
         i = index(trim(string), c)
      end do
   end subroutine clean_string

   !> replace all occurrences of 'src' in 'string' to 'dst'
   subroutine replace(string, src, dst)
      character(len=*), intent(inout) :: string
      character(len=*), intent(in) :: src, dst
      integer :: src_len, dst_len, i, i0
      src_len = len(trim(src))
      dst_len = len(trim(dst))
      i = index(string, src)
      i0 = 0
      do while (i .ne. 0)
         string = string(:i0 + i - 1)//trim(dst)//string(i0 + i + src_len:)
         i0 = i0 + i + dst_len - 1
         i = index(string(i0 + 1:), src)
      end do
   end subroutine replace

   !> same as replace but function instead
   function freplace(string, src, dst)
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: src, dst
      character(len=sl) :: freplace
      freplace = string
      call replace(freplace, src, dst)
   end function freplace

   !> checks if string ends with substring
   function endswith(string, substring)
      logical :: endswith
      character(len=*), intent(in) :: string, substring
      endswith = len(trim(string)) >= len(trim(substring))
      if (endswith) endswith = trim(string(len(trim(string)) - len(trim(substring)) + 1:)) == trim(substring)
   end function endswith

   !> checks if string starts with substring
   function startswith(string, substring)
      logical :: startswith
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: substring
      startswith = index(string, trim(substring)) .eq. 1
   end function startswith

   !> retuns the number of times substring occurs in string
   function str_count(string, substring)
      integer :: str_count
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: substring
      integer :: i, len_string, len_substring
      len_string = len(string)
      len_substring = len(substring)
      str_count = 0
      do i = 1, len_string - len_substring
         if (string(i:i + len_substring - 1) == substring) str_count = str_count + 1
      end do
   end function str_count

   !> checks if string contains substring
   function str_contains(string, substring)
      logical :: str_contains
      character(len=*), intent(in) :: string
      character(len=*), intent(in) :: substring
      str_contains = index(string, substring) .ne. 0
   end function str_contains

   function path_join(paths)
      character(len=*), dimension(:), intent(in) :: paths
      character(len=sl) :: path_join
      path_join = fjoin(paths, '/')
      do while (str_contains(path_join, '//'))
         call replace(path_join, '//', '/')
      end do
      if (endswith(path_join, '/')) path_join = path_join(:len(trim(path_join)) - 1)
   end function path_join

   !> Count number or substring in string
   function fcount_str(string, substring)
      implicit none
      character(len=*), intent(in) :: string, substring
      integer :: i, i0, fcount_str
      fcount_str = 0
      i0 = 1
      i = index(string(i0:), substring)
      do while (i .ne. 0)
         fcount_str = fcount_str + 1
         i0 = i0 + i
         i = index(string(i0:), substring)
      end do
   end function fcount_str

   subroutine indent_str(string, amount)
      character(len=*), intent(inout) :: string
      integer, intent(in) :: amount
      string = repeat(' ', amount)//trim(string(:len(string) - amount - 1))
   end subroutine indent_str

   function findent_str(string, amount)
      character(len=*), intent(in) :: string
      integer, intent(in) :: amount
      character(len=sl) :: findent_str
      findent_str = string
      call indent_str(findent_str, amount)
   end function findent_str

   function fmt_txt(format, value) result(txt)
      character(len=*), intent(in) :: format
      character(len=*), intent(in) :: value
      character(len=:), allocatable :: txt
      character(len=sl) :: tmp_txt
      write (tmp_txt, '('//format//')') trim(value)
      txt = trim(tmp_txt)
   end function fmt_txt

   function fmt_integer(format, value) result(txt)
      character(len=*), intent(in) :: format
      integer, intent(in) :: value
      character(len=:), allocatable :: txt
      character(len=sl) :: tmp_txt
      write (tmp_txt, '('//format//')') value
      txt = trim(tmp_txt)
   end function fmt_integer

   function fmt_real4(format, value) result(txt)
      character(len=*), intent(in) :: format
      real(4), intent(in) :: value
      character(len=:), allocatable :: txt
      character(len=sl) :: tmp_txt
      write (tmp_txt, '('//format//')') value
      txt = trim(tmp_txt)
   end function fmt_real4

   function fmt_real8(format, value) result(txt)
      character(len=*), intent(in) :: format
      real(8), intent(in) :: value
      character(len=:), allocatable :: txt
      character(len=sl) :: tmp_txt
      write (tmp_txt, '('//format//')') value
      txt = trim(tmp_txt)
   end function fmt_real8

   ! Function to count lines in a file
   integer function countLines(fileName) result(nLines)
      implicit none
      character(len=*), intent(in) :: fileName
      integer :: iostat2
      character(200) :: line

      nLines = 0
      open (93, file=fileName, status='old', action='read')

      do
         read (93, '(A)', iostat=iostat2) line
         if (iostat2 /= 0) exit
         nLines = nLines + 1
      end do

      close (93)
   end function countLines
end module string_mod
