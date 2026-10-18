#pragma once
#include <algorithm>
#include <atomic>
#include <future>
#include <vector>

namespace slice::detail {
// Each row has one writer. Futures join outstanding workers on exceptions.
// Reductions that determine centers and statistics remain sequential.
template<class Function>
void parallelFor(size_t rows, unsigned threads, Function function) {
    std::atomic<size_t> next{0};
    auto worker = [&] {
        for (size_t i; (i = next.fetch_add(1)) < rows;) function(i);
    };
    std::vector<std::future<void>> workers;
    for (size_t t = 1; t < std::min<size_t>(threads, rows); ++t)
        workers.push_back(std::async(std::launch::async, worker));
    worker();
    for (auto& task : workers) task.get();
}
}
