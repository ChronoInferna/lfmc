#include "lfmc/timer.hpp"

#include <catch2/catch_test_macros.hpp>
#include <thread>

TEST_CASE("Timer measures elapsed time correctly", "[Timer]") {
    lfmc::Timer timer;
    // Simulate some work with a sleep
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    auto elapsed = timer.elapsedMilliseconds();
    REQUIRE(elapsed >= 100);
}

TEST_CASE("Timer can be restarted", "[Timer]") {
    lfmc::Timer timer;
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    auto firstElapsed = timer.elapsedMilliseconds();
    REQUIRE(firstElapsed >= 50);

    timer.reset(); // Restart the timer
    std::this_thread::sleep_for(std::chrono::milliseconds(70));
    auto secondElapsed = timer.elapsedMilliseconds();
    REQUIRE(secondElapsed >= 70);
}

TEST_CASE("Timer without start returns zero elapsed time", "[Timer]") {
    lfmc::Timer timer;
    auto elapsed = timer.elapsedMilliseconds();
    REQUIRE(elapsed == 0);
}

TEST_CASE("ScopedTimer measures elapsed time correctly", "[ScopedTimer]") {
    long long elapsed = 0;
    {
        lfmc::ScopedTimer scopedTimer(elapsed);
        // Simulate some work with a sleep
        std::this_thread::sleep_for(std::chrono::milliseconds(150));
    }
    REQUIRE(elapsed >= 150);
}
