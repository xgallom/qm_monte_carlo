set(files
        main.cpp
        random.cpp
        energyForParameters.cpp
        runSingle.cpp
        runParallel.cpp
        )

add_executable(qm_monte_carlo ${files})
target_link_libraries(qm_monte_carlo pthread)
