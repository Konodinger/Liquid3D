# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.24

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\CMake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\CMake\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation

# Include any dependencies generated for this target.
include dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/compiler_depend.make

# Include the progress variables for this target.
include dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/progress.make

# Include the compile flags for this target's objects.
include dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/flags.make

dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/triangle-vulkan.c.obj: dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/flags.make
dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/triangle-vulkan.c.obj: dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/includes_C.rsp
dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/triangle-vulkan.c.obj: dep/glfw/tests/triangle-vulkan.c
dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/triangle-vulkan.c.obj: dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/triangle-vulkan.c.obj"
	cd /d C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests && C:\Strawberry\c\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/triangle-vulkan.c.obj -MF CMakeFiles\triangle-vulkan.dir\triangle-vulkan.c.obj.d -o CMakeFiles\triangle-vulkan.dir\triangle-vulkan.c.obj -c C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests\triangle-vulkan.c

dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/triangle-vulkan.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/triangle-vulkan.dir/triangle-vulkan.c.i"
	cd /d C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests && C:\Strawberry\c\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests\triangle-vulkan.c > CMakeFiles\triangle-vulkan.dir\triangle-vulkan.c.i

dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/triangle-vulkan.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/triangle-vulkan.dir/triangle-vulkan.c.s"
	cd /d C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests && C:\Strawberry\c\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests\triangle-vulkan.c -o CMakeFiles\triangle-vulkan.dir\triangle-vulkan.c.s

dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/__/deps/glad_vulkan.c.obj: dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/flags.make
dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/__/deps/glad_vulkan.c.obj: dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/includes_C.rsp
dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/__/deps/glad_vulkan.c.obj: dep/glfw/deps/glad_vulkan.c
dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/__/deps/glad_vulkan.c.obj: dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/__/deps/glad_vulkan.c.obj"
	cd /d C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests && C:\Strawberry\c\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/__/deps/glad_vulkan.c.obj -MF CMakeFiles\triangle-vulkan.dir\__\deps\glad_vulkan.c.obj.d -o CMakeFiles\triangle-vulkan.dir\__\deps\glad_vulkan.c.obj -c C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\deps\glad_vulkan.c

dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/__/deps/glad_vulkan.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/triangle-vulkan.dir/__/deps/glad_vulkan.c.i"
	cd /d C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests && C:\Strawberry\c\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\deps\glad_vulkan.c > CMakeFiles\triangle-vulkan.dir\__\deps\glad_vulkan.c.i

dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/__/deps/glad_vulkan.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/triangle-vulkan.dir/__/deps/glad_vulkan.c.s"
	cd /d C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests && C:\Strawberry\c\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\deps\glad_vulkan.c -o CMakeFiles\triangle-vulkan.dir\__\deps\glad_vulkan.c.s

# Object files for target triangle-vulkan
triangle__vulkan_OBJECTS = \
"CMakeFiles/triangle-vulkan.dir/triangle-vulkan.c.obj" \
"CMakeFiles/triangle-vulkan.dir/__/deps/glad_vulkan.c.obj"

# External object files for target triangle-vulkan
triangle__vulkan_EXTERNAL_OBJECTS =

dep/glfw/tests/triangle-vulkan.exe: dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/triangle-vulkan.c.obj
dep/glfw/tests/triangle-vulkan.exe: dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/__/deps/glad_vulkan.c.obj
dep/glfw/tests/triangle-vulkan.exe: dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/build.make
dep/glfw/tests/triangle-vulkan.exe: dep/glfw/src/libglfw3.a
dep/glfw/tests/triangle-vulkan.exe: dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/linklibs.rsp
dep/glfw/tests/triangle-vulkan.exe: dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/objects1.rsp
dep/glfw/tests/triangle-vulkan.exe: dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable triangle-vulkan.exe"
	cd /d C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\triangle-vulkan.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/build: dep/glfw/tests/triangle-vulkan.exe
.PHONY : dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/build

dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/clean:
	cd /d C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests && $(CMAKE_COMMAND) -P CMakeFiles\triangle-vulkan.dir\cmake_clean.cmake
.PHONY : dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/clean

dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests C:\Users\berge\Documents\Cours\IGR\IGR202\Projet_simulation\dep\glfw\tests\CMakeFiles\triangle-vulkan.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : dep/glfw/tests/CMakeFiles/triangle-vulkan.dir/depend

