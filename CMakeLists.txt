cmake_minimum_required(VERSION 3.10)

# Определение проекта
project(tfdnoise VERSION 1.0 LANGUAGES CXX)

# Установка стандарта C++ и флагов компиляции
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math")

# Поиск MPI (обязательно для Boost.MPI)
find_package(MPI REQUIRED)
if (MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER}) # Используем mpicxx
else()
    message(FATAL_ERROR "MPI not found. Please install MPI (e.g., OpenMPI or MPICH).")
endif()

# Поиск Boost с компонентами mpi и serialization
find_package(Boost 1.65 REQUIRED COMPONENTS mpi serialization)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
else()
    message(FATAL_ERROR "Boost.MPI or Boost.Serialization not found.")
endif()

# Поиск остальных библиотек
find_library(FFTW3F_LIBRARY fftw3f REQUIRED)
find_library(SEGYIO_LIBRARY segyio REQUIRED)

# Указание директорий с заголовочными файлами
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Добавление исполняемого файла
add_executable(tfdnoise
    tfdnoise.cpp
    include/tfd.cpp
)

# Линковка библиотек
target_link_libraries(tfdnoise
    ${MPI_LIBRARIES}
    ${Boost_LIBRARIES}
    ${FFTW3F_LIBRARY}
    ${SEGYIO_LIBRARY}
)

# Установка правил инсталляции (опционально)
install(TARGETS tfdnoise DESTINATION bin)
