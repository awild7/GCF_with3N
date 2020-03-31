add_executable(genPhoto genPhoto.cc ${LIBRARIES}/cross_sections/photoCrossSection.cc ${LIBRARIES}/generator/gcfGenerator.cc ${LIBRARIES}/nucleus/gcfNucleus.cc photoGenerator.cc)
target_link_libraries(genPhoto ${ROOT_LIBRARIES})
