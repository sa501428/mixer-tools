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
    const size_t chunk=std::max<size_t>(1,std::min<size_t>(16,rows/std::max(1u,threads)/4));
    auto worker = [&] {
        for (size_t begin; (begin = next.fetch_add(chunk)) < rows;)
            for(size_t i=begin;i<std::min(rows,begin+chunk);++i) function(i);
    };
    std::vector<std::future<void>> workers;
    for (size_t t = 1; t < std::min<size_t>(threads, rows); ++t)
        workers.push_back(std::async(std::launch::async, worker));
    worker();
    for (auto& task : workers) task.get();
}
}
