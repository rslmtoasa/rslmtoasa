
!------------------------------------------------------------------------------
! RS-LMTO-ASA
!------------------------------------------------------------------------------
!
! MODULE: Parameter Parser
!
!> @author
!> Angela Klautau
!> Ramon Cardias
!> Lucas P. Campagna
!> S. Frota-Pessôa
!> Pascoal R. Peduto
!> Anders Bergman
!> S. B. Legoas
!> H. M. Petrilli
!> Ivan P. Miranda
!
! DESCRIPTION:
!> Module to handle with the input call parameters
!------------------------------------------------------------------------------

module os_mod

   use string_mod, only: sl, split, replace, fmt
   ! use math_mod, only: cross_product
   use precision_mod, only: rp
   use logger_mod, only: g_logger
   implicit none

   private

   !> Module's main structure
   type :: argument_parser
      !> Name of input file with all needed namelists
      character(len=sl) :: input_origin
      !> Name of input file with namelists replaced by those in namelists
      character(len=sl) :: input
      !> Output subfolder
      character(len=sl) :: output
      !> A list with namelists and files for this namelists to overwrite the values readed from input
      character(len=sl), dimension(:, :), allocatable :: namelists
   end type argument_parser

   interface argument_parser
      procedure :: argument_parser_constructor
   end interface argument_parser

   public :: argument_parser, exists

contains
   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Constructor
   !
   !> @return type(argument_parser)
   !---------------------------------------------------------------------------
   function argument_parser_constructor() result(obj)
      type(argument_parser) :: obj
      character(len=sl) :: argument
      character(len=sl), dimension(:), allocatable :: keyvalue
      character(len=sl) :: output_input, dummy_temp_file1, dummy_temp_file2
      character(len=sl), dimension(:, :), allocatable :: aux_namelists
      logical :: has_provided_output = .false.
      character(len=sl) :: nml, file
      character(len=sl) :: code
      logical :: checkpoint_exists
      ! iterable variables
      integer :: i, j

      code = generate_code(5)

      ! (1) Read arguments
      obj%input_origin = "input.nml"
      obj%input = ""
      obj%output = "./"
      allocate (obj%namelists(max(command_argument_count() - 1, 0), 2))
      if (command_argument_count() .gt. 0) then
         call get_command_argument(1, obj%input_origin)
      end if
      j = 1
      do i = 2, command_argument_count()
         call get_command_argument(i, argument)
         call split(argument, '=', keyvalue)
         if (keyvalue(1) .eq. 'output') then
            if (has_provided_output) then
               call g_logger%fatal('More than one output provided on argument list.', __FILE__, __LINE__)
            end if
            obj%output = keyvalue(2)
            has_provided_output = .true.
            ! reduce by 1 the amount of namelists on argument list
            if (size(obj%namelists) .gt. 1) then
               call move_alloc(obj%namelists, aux_namelists)
               allocate (obj%namelists(size(aux_namelists(:, 1)) - 1, 2))
               obj%namelists(1:j, :) = aux_namelists(1:j, :)
            end if
         else
            obj%namelists(j, 1) = keyvalue(1) ! namelist
            obj%namelists(j, 2) = keyvalue(2) ! file
            j = j + 1
         end if
      end do
      if (obj%output == './') then
         obj%input = obj%input_origin
         ! if has no namelist to merge then return
         if (size(obj%namelists(:, 1)) .eq. 0) then
            return
         end if
         ! add suffix '_modified' to keep original input safer
         call replace(obj%input, '.nml', '_modified.nml')
         output_input = trim(obj%input)
      else
         obj%input = trim(obj%input_origin)
         output_input = trim(obj%output)//"/"//trim(obj%input)
      end if

      ! (2) create output directory
      call system("mkdir -p "//trim(obj%output)); 
      call system("cp "//trim(obj%input_origin)//" "//trim(output_input)); 
      ! (3) create processed input and save to output directory
      do i = 1, size(obj%namelists(:, 1))
         nml = obj%namelists(i, 1)
         file = obj%namelists(i, 2)
         if (.not. is_string_in_file(file, '&'//trim(nml))) then
            call g_logger%fatal('Namelist '//fmt('A', nml)//' not present in file '//fmt('A', file), __FILE__, __LINE__)
         end if
         dummy_temp_file1 = '.'//trim(nml)//"_namelist_cach_"//trim(code)//".tmp"
         ! (3.a) Copy all namelist from input file but 'nml'
         call system("awk 'BEGIN{flag=1} /&"//trim(nml)//"/{flag=0} /^\s*\//{if(flag==0){flag=1;next;}} flag' "//trim(output_input)//" > "//trim(dummy_temp_file1)//"; mv "//trim(dummy_temp_file1)//" "//trim(output_input))
         ! (3.b) Open namelist field in file
         call system("echo '&"//trim(nml)//"' >> "//trim(output_input))
         ! (3.c) Cache namelist content from obj%input
         call system("awk 'BEGIN{flag=0} /&"//trim(nml)//"/{flag=1;next;} /^\s*\//{if(flag==1){flag=0;next;}} flag' "//trim(obj%input)//" > "//trim(dummy_temp_file1))
         ! (3.d) Cache namelist content from file
         call system("awk 'BEGIN{flag=0} /&"//trim(nml)//"/{flag=1;next;} /^\s*\//{if(flag==1){flag=0;next;}} flag' "//trim(file)//" >> "//trim(dummy_temp_file1))
         ! (3.e) Merge namelists content (keeping the last)
         call system("awk -F= '"//'BEGIN{line=1}{gsub("!.*", "")}{array[$1]=(line++)"="$0}END{for(i in array) print(array[i])}'//"' "//trim(dummy_temp_file1)//" | sort | cut -f2, 3 -d= >> "//trim(output_input))
         ! (3.f) Close namelist field in file
         call system("echo / >> "//trim(output_input))
         ! (3.g) Remove cache files
         call system("rm -f "//trim(dummy_temp_file1))
      end do

      if (obj%output .ne. './') then
         ! (4) copy all needed files to output directory
         ! (4.a) files from atoms namelist
         dummy_temp_file1 = '.RSLMTO_DATABASE_'//trim(code)//'.tmp'
         dummy_temp_file2 = '.RSLMTO_LABELS_'//trim(code)//'.tmp'
         call system('grep -E "^\s*database\s*=" '//trim(output_input)//' | cut -f2 -d= | sed "s/'//"'\|"//'\"'//'//g" > '//trim(dummy_temp_file1))
         call system('test $(cat '//trim(dummy_temp_file1)//') && grep -E "^\s*label\s*.*=" '//trim(output_input)//' | cut -f2 -d= | sed "s/'//"'\|"//'\"'//'//g" > '//trim(dummy_temp_file2))
         call system('test $(cat '//trim(dummy_temp_file1)//') && for i in $(cat '//trim(dummy_temp_file2)//'); do cp $(cat '//trim(dummy_temp_file1)//')${i}.nml '//trim(obj%output)//'; done')
         call system('rm -f '//trim(dummy_temp_file1)//' '//trim(dummy_temp_file2))
      end if

      ! (5) change current working directory
      call chdir(obj%output)
   end function argument_parser_constructor

   function is_string_in_file(fname, string)
      logical :: is_string_in_file
      character(len=*), intent(in) :: fname
      character(len=*), intent(in) :: string
      character(len=sl) :: line
      integer :: iostatus, unit
      is_string_in_file = .false.
      open (newunit=unit, file=fname, iostat=iostatus, action='read')
      if (iostatus /= 0) then
         call g_logger%fatal('file '//fmt('A', fname)//' not found', __FILE__, __LINE__)
      end if
      do while (.true.)
         read (unit, '(a)', iostat=iostatus) line
         is_string_in_file = index(line, string) .ne. 0
         if (iostatus .eq. 0 .and. .not. is_string_in_file) then
            cycle
         end if
         exit
      end do
      close (unit)
   end function is_string_in_file

   function generate_code(code_size) result(code)
      integer :: code_size
      character(len=sl) :: code
      real :: r_code
      call random_number(r_code)
      write (code, '(I0)') int(r_code * 10**code_size)
   end function generate_code

   ! function path_join(path1, path2)
   !   character(len=*) :: path1, path2
   !   character(len=2*sl) :: path_join
   !   if(path1(len(path1)) .eq. '/') path1(len(path1)) = ''
   ! end function path_join

   !---------------------------------------------------------------------------
   ! DESCRIPTION:
   !> @brief
   !> Check whether 'fname' exists
   !
   !> @param[in] fname File to check
   !---------------------------------------------------------------------------
   function exists(fname)
      character(len=*), intent(in) :: fname
      logical :: exists
      exists = .false.
      inquire (FILE=fname, EXIST=exists)
   end function exists
end module os_mod
