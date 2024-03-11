#pragma once

#include <chrono>
#include <iostream>
#include <string>

class ScopeTimer
{
public:
    ScopeTimer(const std::string& scope_name)
        : m_name(scope_name)
    {
        m_start = std::chrono::system_clock::now();
    }
    ~ScopeTimer()
    {
        auto                          end  = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = end - m_start;
        std::cout << m_name << " consume " << diff.count() << "s" << std::endl;
    }

private:
    std::string                                        m_name;
    std::chrono::time_point<std::chrono::system_clock> m_start;
};

#define BeginScopeTimer(scope_name) ScopeTimer timer##__LINE__(scope_name)