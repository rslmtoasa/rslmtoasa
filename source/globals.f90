module globals_mod

   ! File/Folder character size
   integer, parameter :: GLOBAL_CHAR_SIZE = 300

   ! Project folder root (filled by Makefile DO NOT EDIT)
   character(len=*), parameter :: GLOBAL_ROOT_FOLDER = '/home/lucas/projetos-novos/pesquisa/rslmto/rslmto'

   ! Project database folder
   character(len=*), parameter :: GLOBAL_DATABASE_FOLDER = GLOBAL_ROOT_FOLDER//'/database/elements'

end module globals_mod
