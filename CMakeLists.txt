# CMakeList.txt: vasco-12 的 CMake 项目
cmake_minimum_required (VERSION 3.8)

# 强制设置 vcpkg 工具链路径,CMAKE_TOOLCHAIN_FILE 变量设置成vcpkg的根目录下的脚本路径
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "E:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "Vcpkg toolchain file")
endif()

# 如果支持，请为 MSVC 编译器启用热重载。
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  # 修改：统一使用 ProgramDatabase (/Zi)。
  # 1. 兼容 /MP 多核编译（EditAndContinue 不兼容 /MP）。
  # 2. 确保 Release 模式也生成调试信息。
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "ProgramDatabase")
endif()

project ("vasco-12")

# 调试信息：确认工具链是否生效
message(STATUS "Current Toolchain File: ${CMAKE_TOOLCHAIN_FILE}")

# 要求 C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 开启多处理器编译（MSVC）
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")

# 新增：为 Release 配置启用 PDB 生成
if (MSVC)
  # /DEBUG: 生成 PDB 文件
  # /OPT:REF /OPT:ICF: 保持 Release 模式的优化（移除未引用函数/数据折叠）
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /DEBUG /OPT:REF /OPT:ICF")
endif()

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)

# Polyscope
include(FetchContent)

FetchContent_Declare(
   polyscope
   GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
   GIT_TAG        master
)

FetchContent_MakeAvailable(polyscope)



# 查找第三方依赖
# 注意：顺序很重要。底层依赖（如 ZLIB, TIFF）应先查找。
find_package(Eigen3 CONFIG REQUIRED)
find_package(CGAL CONFIG REQUIRED)
find_package(OpenCV CONFIG REQUIRED)

# 工程源文件
set(SRC_FILES
  data.cpp
  edge.cpp
  GeneralMesh.cpp
  HybridManufacturing.cpp
  layer_graph.cpp
  Mygraph.cpp
  polygon.cpp
  sample_on_ball.cpp
  visual.cpp
  katana/config.cpp
  katana/gcode.cpp
  #katana/infill.cpp
  #katana/katana.cpp
  katana/slicer.cpp
  katana/stl.cpp
  stl2obj/exportobj.cpp
  stl2obj/importstl.cpp
  #stl2obj/stl2obj.cpp
  vasco/generateClosedTube.cpp
  vasco/core/VascoUtils.cpp
  vasco/core/GeometryUtils.cpp
  vasco/core/GraphUtils.cpp
  vasco/core/Voronoi.cpp
  vasco/Visualization.cpp
  gco-v3.0/GCoptimization.cpp
  gco-v3.0/LinkedBlockList.cpp
  main.cpp
)

set(HDR_FILES
  HybridManufacturing.h
  MeshPlaneIntersect.hpp
  katana/config.h
  katana/datastructures.h
  katana/debug.h
  katana/gcode.h
  #katana/infill.h
  #katana/katana.h
  katana/slicer.h
  katana/stl.h
  stl2obj/exportobj.h
  stl2obj/geombase.h
  stl2obj/geometry.h
  stl2obj/importstl.h
  stl2obj/kdtree.h
  stl2obj/vectornd.h
  stl2obj/visitor.h

  vasco/generateClosedTube.h
  vasco/Slicer.h
  vasco/Visualization.h
  vasco/core/Types.h
  vasco/core/VascoUtils.h
  vasco/core/CgalUtils.h
  vasco/core/GraphUtils.h
  vasco/core/Voronoi.h
  vasco/core/ObjIO.h
  vasco/core/Constants.h
  vasco/core/GeometryUtils.h
  Mygraph.h
  helpers.h
  edge.h
  sample_on_ball.h
  gco-v3.0/GCoptimization.h      # 若源码内置
  # 其他头文件可继续补充
)

# 保持 VS 的过滤器视图（按目录树显示）
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SRC_FILES} ${HDR_FILES})

add_executable(HybridManufacturing
  ${SRC_FILES}
  ${HDR_FILES}
)

# 头文件包含路径（相对路径为工程根目录）
target_include_directories(HybridManufacturing
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/vasco
    ${CMAKE_CURRENT_SOURCE_DIR}/vasco/core
    ${CMAKE_CURRENT_SOURCE_DIR}/gco-v3.0
)

# 链接依赖
target_link_libraries(HybridManufacturing
  PUBLIC
    igl::glfw
    polyscope

  PRIVATE
    Eigen3::Eigen
    CGAL::CGAL
    ${OpenCV_LIBS}
)

# 定义源代码根目录宏，用于在代码中构建绝对路径
# target_compile_definitions(HybridManufacturing PRIVATE "VASCO_ROOT_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\"")

# 设置调试工作目录为项目根目录，以便直接读取根目录下的资源文件
set_target_properties(HybridManufacturing PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

if (MSVC)
  #target_compile_definitions(HybridManufacturing PRIVATE _SILENCE_CXX17_OLD_ALLOCATOR_MEMBERS_DEPRECATION_WARNING)
endif()

# 这里的 vasco-12 目标似乎是多余的或者是测试用的？保留原样。
add_executable (vasco-12 "vasco-12.cpp" "vasco-12.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET vasco-12 PROPERTY CXX_STANDARD 20)
endif()

set_target_properties(vasco-12 PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

# TODO: 如有需要，请添加测试并安装目标。
