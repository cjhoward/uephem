project(uephem)

# Set compiler flags
if(CMAKE_COMPILER_IS_GNUCC)
	set(CMAKE_C_FLAGS "-Wall -Wextra -Os")
elseif(MSVC)
	set(CMAKE_C_FLAGS "/MP /W3 /MT /GR- /GS- /Os /D_WIN32_WINNT=0x0501 /DWINVER=0x0501")
endif()
set(CMAKE_C_FLAGS_DEBUG "")
set(CMAKE_C_FLAGS_RELEASE "")

# Set linker flags
if(MSVC)
	set(CMAKE_EXE_LINKER_FLAGS "/nologo /OPT:REF /OPT:ICF /INCREMENTAL:NO")
endif()
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "")

add_executable(uephem uephem.c)
