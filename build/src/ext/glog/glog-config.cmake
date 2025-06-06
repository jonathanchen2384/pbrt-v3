
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was glog-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

####################################################################################

include ("${CMAKE_CURRENT_LIST_DIR}/glog-targets.cmake")
set_and_check (glog_INCLUDE_DIR "/Users/jonathanchen/Desktop/pbrt-v3/src/ext/glog/src")



set (glog_LIBRARY glog)

set (glog_LIBRARIES ${glog_LIBRARY})
set (glog_INCLUDE_DIRS ${glog_INCLUDE_DIR})
