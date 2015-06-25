set(BAMREAD_SRCS bam-read.cpp pileup.cpp segmenter.cpp)
set(BAMREAD_HEADERS pileup.h segmenter.h readerwriterqueue.h atomicops.h)


include_directories(${GAT_SOURCE_DIR}/bamtools/include)
link_directories(${GAT_SOURCE_DIR}/bamtools/lib)

set(LIBRARIES ${GAT_SOURCE_DIR}/bamtools/lib/libbamtools.a ${ZLIB_LIBRARIES})

add_executable(bam-read ${BAMREAD_SRCS} ${BAMREAD_HEADERS})

set_target_properties(bam-read PROPERTIES LINK_SEARCH_START_STATIC 1)
set_target_properties(bam-read PROPERTIES LINK_SEARCH_END_STATIC 1)
target_link_libraries(bam-read ${LIBRARIES})
