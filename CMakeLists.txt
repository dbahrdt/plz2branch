cmake_minimum_required(VERSION 3.18)
project(plz2branch)

set(CMAKE_MODULE_PATH
	${CMAKE_SOURCE_DIR}/cmake
)

include(cmake/add_buildtype_ultra_sanitize.cmake)

add_subdirectory(vendor/sserialize sserialize)
add_subdirectory(vendor/libosmtools libosmtools)
add_subdirectory(vendor/memgraph memgraph)

find_package(Qt5Widgets REQUIRED)
find_package(Qt5Gui REQUIRED)
find_package(Marble REQUIRED)
find_package(OpenMP REQUIRED)

#TODO: Find package coin-or-cbc for optimizer which can be found using pkg-config

set(SOURCES_CPP
	src/main.cpp
	src/MainWindow.cpp
	src/State.cpp
	src/BranchTableModel.cpp
	src/MarbleMap.cpp
	src/optimizer.cpp
	src/FileFormat.cpp
)

set(SOURCES_MOC_H
	src/MainWindow.h
	src/State.h
	src/BranchTableModel.h
	src/MarbleMap.h
)

set(MY_LINK_LIBS
	sserialize
	libosmtools
	memgraph
	${MARBLE_LIBRARIES}
	Qt5::Widgets
	Qt5::Gui
	OpenMP::OpenMP_CXX
)

set(MY_INCLUDE_DIRS
	${QT_INCLUDES}
	${MARBLE_INCLUDE_DIR}
)

qt5_wrap_cpp(SOURCES_MOC_CPP ${SOURCES_MOC_H})

add_executable(${PROJECT_NAME} ${SOURCES_CPP} ${SOURCES_MOC_CPP})
target_link_libraries(${PROJECT_NAME} ${MY_LINK_LIBS})
target_include_directories(${PROJECT_NAME} PRIVATE ${MY_INCLUDE_DIRS})
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
